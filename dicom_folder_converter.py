#!/usr/bin/env python3
"""
DICOM PET/CT folder → NIfTI database converter.

Input layout
------------
    <DICOM folder>/            e.g. DICOM_P26235346
        <any>/  *.dcm          series folders, any names, any depth

Series are identified by the DICOM Modality tag, not by folder name:
Modality CT → *_CT.nii.gz, Modality PT → *_PET.nii.gz (+ *_SUV.nii.gz).
If several series share a modality, CT skips localizers and takes the one
with the most files; PET prefers attenuation-corrected, then most files.
The CT must have the same PatientID as the PET. Other series are ignored.

Output layout ("database")
--------------------------
    <output>/
        metadata.csv
        0001/                          integer patient ID (PID)
            202611104821/              {PET acquisition date YYYYMMDD}{4 random digits}
                0001_202611104821_CT.nii.gz
                0001_202611104821_PET.nii.gz
                0001_202611104821_SUV.nii.gz

The exam ID is "{PID}_{date}{rand}", e.g. 0001_202611104821.

metadata.csv has one row per exam. Key columns come first (ID, origin
folder, PET acquisition date/time, injected dose, patient weight, ...),
followed by every non-private, non-binary header element of the first CT and
first PET slice, prefixed with CT_ / PET_ (sequences are flattened).

Modes
-----
  batch  : input is a directory containing many DICOM folders.
  single : input is one DICOM folder.
Both write into the chosen output folder: if it already has metadata.csv,
new exams are appended to it (a patient already present, matched on DICOM
PatientID, keeps the same PID); otherwise a new database is created.

The input DICOM files are only read, never modified.

Usage
-----
  python dicom_folder_converter.py                       # GUI
  python dicom_folder_converter.py batch  IN_DIR OUT_DIR [--start-id 0001]
  python dicom_folder_converter.py single IN_DIR OUT_DIR [--start-id 0001]
"""

import argparse
import csv
import datetime
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import threading
import traceback

try:
    import numpy as np
    NUMPY_OK = True
except ImportError:
    NUMPY_OK = False

try:
    import pydicom
    PYDICOM_OK = True
except ImportError:
    PYDICOM_OK = False

try:
    import nibabel as nib
    NIBABEL_OK = True
except ImportError:
    NIBABEL_OK = False


METADATA_NAME = "metadata.csv"
RAND_DIGITS = 4

# Columns placed first in metadata.csv, in this order.
KEY_COLUMNS = [
    "ID",
    "PID",
    "ExamFolder",
    "Origin",
    "PatientName",
    "PET_ScanDateTime",
    "PET_AcquisitionDate",
    "PET_AcquisitionTime",
    "PET_InjectedDose_Bq",
    "PET_PatientWeight",
    "PET_PatientSize",
    "PET_InjectionDateTime",
    "PET_HalfLife_s",
    "PET_Radiopharmaceutical",
    "PET_Units",
    "PET_DecayCorrection",
    "PET_SeriesDate",
    "PET_SeriesTime",
    "SUV_Status",
    "SUV_ScanDateTimeUsed",
    "SUV_DecayTime_s",
    "SUV_Factor",
    "CT_File",
    "PET_File",
    "SUV_File",
    "CT_NumFiles",
    "PET_NumFiles",
    "ConvertedAt",
]

# Value representations we do not put into the CSV (binary blobs).
_BINARY_VRS = {"OB", "OW", "OF", "OD", "OL", "OV", "UN"}
_PIXEL_DATA_TAGS = {0x7FE00010, 0x7FE00008, 0x7FE00009}


# ---------------------------------------------------------------------------
# DICOM helpers (header only — pixel data is handled by dcm2niix)
# ---------------------------------------------------------------------------

def _series_files(series_dir):
    """All regular files in a series folder (DICOM files may lack .dcm)."""
    out = []
    for f in sorted(os.listdir(series_dir)):
        p = os.path.join(series_dir, f)
        if f.startswith(".") or f.upper() == "DICOMDIR" or not os.path.isfile(p):
            continue
        out.append(p)
    return out


def _read_header(path, tags=None):
    try:
        return pydicom.dcmread(
            path, stop_before_pixels=True, force=True, specific_tags=tags,
        )
    except Exception:
        return None


_MODALITY_ALIASES = {"CT": "CT", "PT": "PT", "PET": "PT"}
_SAMPLE_FILES = 15


def _sample(files, n=_SAMPLE_FILES):
    if len(files) <= n:
        return files
    step = (len(files) - 1) / (n - 1)
    return [files[round(i * step)] for i in range(n)]


def _probe_series(series_dir):
    """Classify a folder from the DICOM tags of a sample of its files.

    Returns a dict (modality, n_files, mixed, localizer, attn_corrected)
    or None if the folder holds no readable DICOM.
    """
    files = _series_files(series_dir)
    if not files:
        return None
    votes = {}
    localizer = attn = False
    patient_id = study_uid = ""
    for f in _sample(files):
        ds = _read_header(f, ["Modality", "ImageType", "CorrectedImage",
                                "PatientID", "StudyInstanceUID"])
        if ds is None or "Modality" not in ds:
            continue
        m = str(ds.Modality).strip().upper()
        votes[m] = votes.get(m, 0) + 1
        itype = [str(v).upper() for v in (getattr(ds, "ImageType", None) or [])]
        if "LOCALIZER" in itype:
            localizer = True
        corr = getattr(ds, "CorrectedImage", None) or []
        if isinstance(corr, str):
            corr = [corr]
        if any("ATTN" in str(v).upper() for v in corr):
            attn = True
        if not patient_id:
            patient_id = str(getattr(ds, "PatientID", "") or "").strip()
        if not study_uid:
            study_uid = str(getattr(ds, "StudyInstanceUID", "") or "").strip()
    if not votes:
        return None
    modality = max(votes, key=votes.get)
    return {
        "modality": _MODALITY_ALIASES.get(modality, modality),
        "n_files": len(files),
        "mixed": len(votes) > 1,
        "localizer": localizer,
        "attn": attn,
        "patient_id": patient_id,
        "study_uid": study_uid,
    }


def _scan_series(series_dir):
    """Return (first_file_by_InstanceNumber, n_files, earliest (AcqDate, AcqTime))."""
    files = _series_files(series_dir)
    entries = []
    earliest = None
    for f in files:
        ds = _read_header(f, ["InstanceNumber", "AcquisitionDate", "AcquisitionTime", "SOPClassUID"])
        if ds is None:
            continue
        try:
            n = int(ds.InstanceNumber)
        except Exception:
            n = 10 ** 9
        entries.append((n, f))
        acq = _parse_da_tm(getattr(ds, "AcquisitionDate", ""), getattr(ds, "AcquisitionTime", ""))
        if acq is not None and (earliest is None or acq < earliest):
            earliest = acq
    if not entries:
        return None, 0, None
    entries.sort()
    return entries[0][1], len(entries), earliest


def find_ct_pet(exam_dir, log_fn=print):
    """Locate the CT and PET series folders inside one DICOM folder.

    Folder names are not trusted: every folder under *exam_dir* (any depth,
    including *exam_dir* itself) that holds DICOM files is classified by its
    Modality tag (CT / PT). When several series share a modality:
      CT  : localizer/scout series are skipped, then the most files wins;
      PET : attenuation-corrected (CorrectedImage contains ATTN) first,
            then the most files.
    Returns (ct_dir, pet_dir, ignored) where *ignored* lists the other
    series folders (relative paths) with their modality.
    """
    series = []
    for root, dirs, _files in os.walk(exam_dir):
        dirs[:] = sorted(d for d in dirs if not d.startswith("."))
        info = _probe_series(root)
        if info is None:
            continue
        info["path"] = root
        info["rel"] = os.path.relpath(root, exam_dir)
        series.append(info)
        if info["mixed"]:
            log_fn(f"  [WARN] {info['rel']}: files with different Modality in one folder; "
                   f"treated as {info['modality']}")

    def pick(modality, rank):
        cands = [s for s in series if s["modality"] == modality]
        return max(cands, key=rank) if cands else None

    pet = pick("PT", lambda s: (s["attn"], s["n_files"]))
    ct_rank = lambda s: (not s["localizer"], s["n_files"])
    ct = None
    if pet is not None:
        # CT must belong to the same patient as the PET; prefer the same study.
        same_study = [s for s in series if s["modality"] == "CT" and pet["study_uid"]
                      and s["study_uid"] == pet["study_uid"]]
        same_patient = [s for s in series if s["modality"] == "CT"
                        and s["patient_id"] == pet["patient_id"]]
        cands = same_study or same_patient
        if cands:
            ct = max(cands, key=ct_rank)
        elif pick("CT", ct_rank) is not None:
            log_fn("  [WARN] CT series found, but none has the same PatientID as the PET")
    else:
        ct = pick("CT", ct_rank)

    for s, label in ((ct, "CT"), (pet, "PET")):
        if s is not None:
            log_fn(f"  {label} series: {s['rel']} (Modality={s['modality']}, {s['n_files']} files)")

    chosen = {id(s) for s in (ct, pet) if s is not None}
    ignored = [f"{s['rel']} ({s['modality']})" for s in series if id(s) not in chosen]
    return (ct["path"] if ct else None), (pet["path"] if pet else None), ignored


def _value_to_str(value):
    if isinstance(value, (bytes, bytearray)):
        return None
    if isinstance(value, (list, tuple)) or type(value).__name__ == "MultiValue":
        parts = []
        for v in value:
            s = _value_to_str(v)
            if s is None:
                return None
            parts.append(s)
        return "\\".join(parts)
    s = str(value)
    return s.replace("\x00", "").replace("\r", " ").replace("\n", " ").strip()


def flatten_header(ds, prefix, out=None):
    """Flatten a pydicom Dataset to {prefix_Keyword: str}. Sequences recurse;
    item k>1 of a sequence gets a _k suffix. Private and binary tags skipped."""
    if out is None:
        out = {}
    for elem in ds:
        tag = int(elem.tag)
        if elem.tag.is_private or tag in _PIXEL_DATA_TAGS or elem.VR in _BINARY_VRS:
            continue
        name = elem.keyword or f"{elem.tag.group:04X}{elem.tag.element:04X}"
        key = f"{prefix}_{name}"
        if elem.VR == "SQ":
            for i, item in enumerate(elem.value or []):
                flatten_header(item, key if i == 0 else f"{key}_{i + 1}", out)
            continue
        try:
            s = _value_to_str(elem.value)
        except Exception:
            s = None
        if s is not None:
            out[key] = s
    return out


# ---------------------------------------------------------------------------
# Date / time + SUV
# ---------------------------------------------------------------------------

def _parse_da_tm(da, tm):
    """DICOM DA ('YYYYMMDD') + TM ('HHMMSS.ffffff') → datetime or None.
    A missing date is allowed and gives 1900-01-01."""
    tm = str(tm or "").strip().replace(":", "")
    if not tm:
        return None
    main, _, frac = tm.partition(".")
    if not main.isdigit() or len(main) < 2:
        return None
    main = main.ljust(6, "0")[:6]
    frac = (frac + "000000")[:6] if frac.isdigit() else "000000"
    da = str(da or "").strip().replace("-", "").replace(".", "")
    try:
        base = datetime.datetime.strptime(da[:8], "%Y%m%d") if len(da) >= 8 else datetime.datetime(1900, 1, 1)
        return base.replace(
            hour=int(main[0:2]), minute=int(main[2:4]), second=int(main[4:6]),
            microsecond=int(frac),
        )
    except ValueError:
        return None


def _parse_dt(dt):
    """DICOM DT ('YYYYMMDDHHMMSS.ffffff&ZZXX') → datetime or None."""
    dt = str(dt or "").strip()
    if len(dt) < 10:
        return None
    for sep in ("+", "-"):
        if sep in dt[8:]:
            dt = dt[:8] + dt[8:].split(sep)[0]
    return _parse_da_tm(dt[:8], dt[8:])


def _fmt_dt(d):
    return d.strftime("%Y%m%d %H%M%S") if d else ""


def compute_suv_params(pet_ds, earliest_acq):
    """Return dict with SUV factor (body-weight SUV) and details.

    SUV = activity[Bq/ml] * weight[g] / (dose[Bq] * 2^(-dt/T½)).
    The scan reference time is SeriesDate/Time (decay correction START),
    or the earliest AcquisitionTime if that is earlier (QIBA rule).
    Without the required tags the factor is None — no SUV guessed.
    """
    info = {
        "PET_InjectedDose_Bq": "", "PET_InjectionDateTime": "", "PET_HalfLife_s": "",
        "PET_Radiopharmaceutical": "", "SUV_Status": "", "SUV_ScanDateTimeUsed": "",
        "SUV_DecayTime_s": "", "SUV_Factor": "",
    }
    factor = None
    units = str(getattr(pet_ds, "Units", "")).strip().upper()
    decay_corr = str(getattr(pet_ds, "DecayCorrection", "")).strip().upper()

    rp = None
    seq = getattr(pet_ds, "RadiopharmaceuticalInformationSequence", None)
    if seq:
        rp = seq[0]

    dose = half_life = None
    inj_dt = None
    if rp is not None:
        info["PET_Radiopharmaceutical"] = str(getattr(rp, "Radiopharmaceutical", "")).strip()
        try:
            dose = float(rp.RadionuclideTotalDose)
            info["PET_InjectedDose_Bq"] = f"{dose:g}"
        except Exception:
            pass
        try:
            half_life = float(rp.RadionuclideHalfLife)
            info["PET_HalfLife_s"] = f"{half_life:g}"
        except Exception:
            pass
        inj_dt = _parse_dt(getattr(rp, "RadiopharmaceuticalStartDateTime", ""))
        if inj_dt is None:
            ref_date = (getattr(pet_ds, "SeriesDate", "") or getattr(pet_ds, "AcquisitionDate", "")
                        or getattr(pet_ds, "StudyDate", ""))
            inj_dt = _parse_da_tm(ref_date, getattr(rp, "RadiopharmaceuticalStartTime", ""))
            inj_time_only = True
        else:
            inj_time_only = False
        info["PET_InjectionDateTime"] = _fmt_dt(inj_dt)

    try:
        weight_kg = float(pet_ds.PatientWeight)
        if weight_kg <= 0:
            weight_kg = None
    except Exception:
        weight_kg = None

    series_dt = _parse_da_tm(getattr(pet_ds, "SeriesDate", ""), getattr(pet_ds, "SeriesTime", ""))
    scan_dt = series_dt
    if earliest_acq is not None and (scan_dt is None or earliest_acq < scan_dt):
        scan_dt = earliest_acq

    missing = []
    if units == "GML":
        factor = 1.0
        info["SUV_Status"] = "ok (PET already in g/ml)"
    elif units not in ("BQML", ""):
        info["SUV_Status"] = f"skipped: unsupported Units={units}"
    else:
        if weight_kg is None:
            missing.append("PatientWeight")
        if dose is None or dose <= 0:
            missing.append("RadionuclideTotalDose")
        if decay_corr != "ADMIN":
            if half_life is None or half_life <= 0:
                missing.append("RadionuclideHalfLife")
            if inj_dt is None:
                missing.append("RadiopharmaceuticalStartTime")
            if scan_dt is None:
                missing.append("SeriesTime/AcquisitionTime")
        if missing:
            info["SUV_Status"] = "skipped: missing " + ", ".join(missing)
        else:
            if decay_corr == "ADMIN":
                decay = 1.0
                dt_s = 0.0
            else:
                if inj_time_only and inj_dt > scan_dt:
                    inj_dt -= datetime.timedelta(days=1)  # injected before midnight
                    info["PET_InjectionDateTime"] = _fmt_dt(inj_dt)
                dt_s = (scan_dt - inj_dt).total_seconds()
                decay = math.exp(-math.log(2) * dt_s / half_life)
            factor = (weight_kg * 1000.0) / (dose * decay)
            info["SUV_Status"] = "ok"
            info["SUV_DecayTime_s"] = f"{dt_s:.0f}"
            if decay_corr != "ADMIN" and not (0 <= dt_s <= 6 * 3600):
                info["SUV_Status"] = f"ok (check: injection→scan = {dt_s / 60:.0f} min)"

    info["SUV_ScanDateTimeUsed"] = _fmt_dt(scan_dt)
    if factor is not None:
        info["SUV_Factor"] = f"{factor:.8g}"
    return factor, info


# ---------------------------------------------------------------------------
# dcm2niix
# ---------------------------------------------------------------------------

def find_dcm2niix():
    exe = "dcm2niix.exe" if sys.platform == "win32" else "dcm2niix"
    cands = []
    if hasattr(sys, "_MEIPASS"):
        cands.append(os.path.join(sys._MEIPASS, exe))
    cands.append(os.path.join(os.path.dirname(os.path.abspath(sys.executable)), exe))
    cands.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), exe))
    for c in cands:
        if os.path.isfile(c):
            if sys.platform != "win32" and not os.access(c, os.X_OK):
                try:
                    os.chmod(c, 0o755)
                except OSError:
                    pass
            return c
    return shutil.which("dcm2niix")


def run_dcm2niix(series_dir, out_nii, log_fn=print):
    """Convert one series folder to *out_nii*. Returns True on success."""
    binary = find_dcm2niix()
    if not binary:
        log_fn("  [ERROR] dcm2niix not found")
        return False
    with tempfile.TemporaryDirectory(prefix="d2n_") as tmp:
        cmd = [binary, "-z", "y", "-b", "n", "-d", "0", "-f", "series", "-o", tmp, series_dir]
        kwargs = {}
        if sys.platform == "win32":
            kwargs["creationflags"] = getattr(subprocess, "CREATE_NO_WINDOW", 0)
        res = subprocess.run(cmd, capture_output=True, text=True, errors="replace", **kwargs)
        niis = [os.path.join(tmp, f) for f in os.listdir(tmp) if f.endswith(".nii.gz")]
        if not niis:
            tail = (res.stdout + res.stderr).strip().splitlines()[-5:]
            log_fn("  [ERROR] dcm2niix produced no output: " + " | ".join(tail))
            return False
        niis.sort(key=os.path.getsize, reverse=True)
        if len(niis) > 1:
            log_fn(f"  [INFO] dcm2niix produced {len(niis)} volumes; keeping the largest")
        shutil.move(niis[0], out_nii)
    return True


def write_suv_nifti(pet_nii, suv_nii, factor):
    img = nib.load(pet_nii)
    data = np.asarray(img.dataobj, dtype=np.float32) * np.float32(factor)
    hdr = img.header.copy()
    hdr.set_data_dtype(np.float32)
    out = nib.Nifti1Image(data, img.affine, hdr)
    out.set_sform(img.get_sform(), int(img.header["sform_code"]))
    out.set_qform(img.get_qform(), int(img.header["qform_code"]))
    out.header.set_slope_inter(1.0, 0.0)
    nib.save(out, suv_nii)


# ---------------------------------------------------------------------------
# metadata.csv database
# ---------------------------------------------------------------------------

def load_metadata(out_dir):
    path = os.path.join(out_dir, METADATA_NAME)
    if not os.path.isfile(path):
        return [], []
    csv.field_size_limit(min(sys.maxsize, 2 ** 31 - 1))
    with open(path, "r", newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        header = list(reader.fieldnames or [])
    return rows, header


def save_metadata(out_dir, rows, prev_header):
    """Write all rows. Column order: KEY_COLUMNS, then previous header order,
    then any new columns in first-seen order. Written atomically."""
    cols = list(KEY_COLUMNS)
    seen = set(cols)
    for c in prev_header:
        if c not in seen:
            cols.append(c)
            seen.add(c)
    for r in rows:
        for c in r:
            if c not in seen:
                cols.append(c)
                seen.add(c)
    path = os.path.join(out_dir, METADATA_NAME)
    tmp = path + ".tmp"
    with open(tmp, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore", restval="")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)
    return cols


def _patient_key(row):
    pid = (row.get("PET_PatientID") or row.get("CT_PatientID") or "").strip()
    if pid:
        return "id:" + pid
    name = (row.get("PET_PatientName") or row.get("CT_PatientName") or "").strip()
    bd = (row.get("PET_PatientBirthDate") or row.get("CT_PatientBirthDate") or "").strip()
    if name:
        return f"name:{name}|{bd}"
    return None


def _exam_key(row):
    uid = (row.get("PET_SeriesInstanceUID") or "").strip()
    return uid or None


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------

def list_exam_dirs(input_dir, mode):
    input_dir = os.path.abspath(input_dir)
    if mode == "single":
        return [input_dir]
    return [
        os.path.join(input_dir, d) for d in sorted(os.listdir(input_dir))
        if os.path.isdir(os.path.join(input_dir, d)) and not d.startswith(".")
    ]


def convert(exam_dirs, out_dir, start_id="0001", log_fn=print, progress_fn=None,
            cancel_event=None):
    """Convert DICOM folders into the database at *out_dir* (created or appended)."""
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    rows, header = load_metadata(out_dir)
    if rows:
        log_fn(f"Appending to existing database ({len(rows)} exam(s)) → {out_dir}")
    else:
        log_fn(f"Creating new database → {out_dir}")

    # Existing patients / exams
    patients = {}
    used_ids = set()
    max_pid = 0
    width = len(str(start_id).strip()) or 4
    for r in rows:
        pid = (r.get("PID") or "").strip()
        used_ids.add((r.get("ID") or "").strip())
        if pid.isdigit():
            max_pid = max(max_pid, int(pid))
            width = max(width, len(pid))
        k = _patient_key(r)
        if k and pid and k not in patients:
            patients[k] = pid
    known_exams = {k for k in (_exam_key(r) for r in rows) if k}
    if rows:
        next_pid = max_pid + 1
    else:
        try:
            next_pid = int(str(start_id).strip())
        except ValueError:
            next_pid = 1

    rng = random.SystemRandom()
    stats = {"ok": 0, "skipped": 0, "failed": 0}
    total = len(exam_dirs)

    for idx, exam_dir in enumerate(exam_dirs, 1):
        if cancel_event is not None and cancel_event.is_set():
            log_fn("\nCancelled.")
            break
        origin = os.path.basename(os.path.normpath(exam_dir))
        log_fn(f"\n[{idx}/{total}] {origin}")
        staging = None
        try:
            ct_dir, pet_dir, ignored = find_ct_pet(exam_dir, log_fn)
            if ignored:
                log_fn(f"  ignoring series folder(s): {', '.join(ignored)}")
            if not ct_dir or not pet_dir:
                log_fn(f"  [SKIP] CT found: {bool(ct_dir)}, PET found: {bool(pet_dir)}")
                stats["skipped"] += 1
                continue

            ct_first, ct_n, _ = _scan_series(ct_dir)
            pet_first, pet_n, pet_earliest = _scan_series(pet_dir)
            if not ct_first or not pet_first:
                log_fn("  [SKIP] no readable DICOM in CT or PET folder")
                stats["skipped"] += 1
                continue

            ct_ds = _read_header(ct_first)
            pet_ds = _read_header(pet_first)
            row = {}
            flatten_header(ct_ds, "CT", row)
            flatten_header(pet_ds, "PET", row)

            ek = _exam_key(row)
            if ek and ek in known_exams:
                log_fn("  [SKIP] already in database (same PET SeriesInstanceUID)")
                stats["skipped"] += 1
                continue

            factor, suv_info = compute_suv_params(pet_ds, pet_earliest)
            row.update(suv_info)

            # Patient name (DICOM "^" separators → spaces), PET first, CT fallback
            pname = row.get("PET_PatientName") or row.get("CT_PatientName") or ""
            row["PatientName"] = " ".join(pname.replace("^", " ").split())

            # Scan datetime: earliest PET slice acquisition, else PET series datetime
            scan_dt = pet_earliest
            if scan_dt is None or scan_dt.year <= 1900:
                scan_dt = _parse_da_tm(row.get("PET_SeriesDate", ""), row.get("PET_SeriesTime", "")) or scan_dt
            row["PET_ScanDateTime"] = scan_dt.strftime("%Y-%m-%d %H:%M:%S") if scan_dt else ""

            # ---- IDs ----
            pkey = _patient_key(row)
            if pkey and pkey in patients:
                pid = patients[pkey]
                log_fn(f"  patient already in database → PID {pid}")
            else:
                pid = str(next_pid).zfill(width)
                next_pid += 1
                if pkey:
                    patients[pkey] = pid

            date = ""
            for k in ("PET_AcquisitionDate", "PET_SeriesDate", "PET_StudyDate"):
                v = (row.get(k) or "").strip()
                if len(v) >= 8 and v[:8].isdigit():
                    date = v[:8]
                    break
            if pet_earliest is not None and pet_earliest.year > 1900:
                date = pet_earliest.strftime("%Y%m%d")
            date = date or "00000000"

            pid_dir = os.path.join(out_dir, pid)
            while True:
                exam_folder = date + str(rng.randrange(10 ** RAND_DIGITS)).zfill(RAND_DIGITS)
                exam_id = f"{pid}_{exam_folder}"
                if exam_id not in used_ids and not os.path.exists(os.path.join(pid_dir, exam_folder)):
                    break

            # ---- Convert into a staging folder, then move into place ----
            staging = tempfile.mkdtemp(prefix=".staging_", dir=out_dir)
            ct_name = f"{exam_id}_CT.nii.gz"
            pet_name = f"{exam_id}_PET.nii.gz"
            suv_name = f"{exam_id}_SUV.nii.gz"

            log_fn(f"  CT  ({ct_n} files) → {ct_name}")
            if not run_dcm2niix(ct_dir, os.path.join(staging, ct_name), log_fn):
                raise RuntimeError("CT conversion failed")
            log_fn(f"  PET ({pet_n} files) → {pet_name}")
            if not run_dcm2niix(pet_dir, os.path.join(staging, pet_name), log_fn):
                raise RuntimeError("PET conversion failed")

            rel = lambda n: "/".join((pid, exam_folder, n))
            if factor is not None:
                write_suv_nifti(os.path.join(staging, pet_name), os.path.join(staging, suv_name), factor)
                log_fn(f"  SUV → {suv_name}  [{suv_info['SUV_Status']}]")
                row["SUV_File"] = rel(suv_name)
            else:
                log_fn(f"  [WARN] SUV not written: {suv_info['SUV_Status']}")

            final_dir = os.path.join(pid_dir, exam_folder)
            os.makedirs(pid_dir, exist_ok=True)
            os.replace(staging, final_dir)
            staging = None

            row.update({
                "ID": exam_id,
                "PID": pid,
                "ExamFolder": exam_folder,
                "Origin": origin,
                "CT_File": rel(ct_name),
                "PET_File": rel(pet_name),
                "CT_NumFiles": str(ct_n),
                "PET_NumFiles": str(pet_n),
                "ConvertedAt": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            })
            rows.append(row)
            used_ids.add(exam_id)
            if ek:
                known_exams.add(ek)
            header = save_metadata(out_dir, rows, header)
            log_fn(f"  ✓ {exam_id}")
            stats["ok"] += 1
        except Exception as e:
            stats["failed"] += 1
            log_fn(f"  [ERROR] {e}")
            log_fn(traceback.format_exc(limit=3))
        finally:
            if staging and os.path.isdir(staging):
                shutil.rmtree(staging, ignore_errors=True)
            if progress_fn:
                progress_fn(idx / total * 100 if total else 100)

    log_fn(
        f"\nDone: {stats['ok']} converted, {stats['skipped']} skipped, "
        f"{stats['failed']} failed. Database: {out_dir}"
    )
    return stats


def missing_dependencies():
    missing = []
    if not NUMPY_OK:
        missing.append("numpy")
    if not PYDICOM_OK:
        missing.append("pydicom")
    if not NIBABEL_OK:
        missing.append("nibabel")
    if not find_dcm2niix():
        missing.append("dcm2niix")
    return missing


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

def _native_browse_dir(title, initial=None):
    """zenity/kdialog folder picker on Linux (Tk's own picker is awkward)."""
    if not sys.platform.startswith("linux"):
        return None
    initial = initial or os.path.expanduser("~")
    if shutil.which("zenity"):
        cmd = ["zenity", "--file-selection", "--directory",
               "--title", title, "--filename", initial.rstrip("/") + "/"]
    elif shutil.which("kdialog"):
        cmd = ["kdialog", "--getexistingdirectory", initial, "--title", title]
    else:
        return None
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except Exception:
        return None
    return (r.stdout.strip() or None) if r.returncode == 0 else ""


def run_gui():
    import tkinter as tk
    from tkinter import ttk, filedialog, scrolledtext, messagebox
    from tkinter import font as tkfont

    if sys.platform == "win32":
        try:
            from ctypes import windll
            try:
                windll.shcore.SetProcessDpiAwareness(2)
            except Exception:
                windll.user32.SetProcessDPIAware()
        except Exception:
            pass

    class App(tk.Tk):
        def __init__(self):
            super().__init__()
            self.title("DICOM PET/CT Folder → NIfTI Database")
            self._scale = self._detect_scale()
            self._apply_scale()
            self.minsize(int(620 * self._scale), int(560 * self._scale))
            self._cancel_event = threading.Event()
            self._running = False
            self._build()
            self._on_mode()
            miss = missing_dependencies()
            if miss:
                self.lbl_deps.config(text="Missing: " + ", ".join(miss))
                self._log(f"[WARN] Missing dependencies: {', '.join(miss)}")

        # ---- DPI ----
        def _detect_scale(self):
            for var in ("DICOM_TOOL_SCALE", "TK_SCALE", "GDK_SCALE"):
                v = os.environ.get(var)
                if v:
                    try:
                        return max(0.5, float(v))
                    except ValueError:
                        pass
            raw = 1.0
            if sys.platform.startswith("linux") and shutil.which("xrdb"):
                try:
                    out = subprocess.check_output(["xrdb", "-query"], text=True, timeout=2)
                    for line in out.splitlines():
                        if line.startswith("Xft.dpi:"):
                            raw = max(raw, float(line.split(":", 1)[1]) / 96.0)
                except Exception:
                    pass
            if raw <= 1.0:
                try:
                    dpi = self.winfo_fpixels("1i")
                    if dpi > 110:
                        raw = dpi / 96.0
                except Exception:
                    pass
            return min(1.4, max(1.0, raw))

        def _apply_scale(self):
            s = self._scale
            self.tk.call("tk", "scaling", (96.0 / 72.0) * s)
            for name in ("TkDefaultFont", "TkTextFont", "TkMenuFont", "TkHeadingFont",
                         "TkCaptionFont", "TkSmallCaptionFont", "TkIconFont",
                         "TkTooltipFont", "TkFixedFont"):
                try:
                    f = tkfont.nametofont(name)
                    f.configure(size=int(round((abs(f.cget("size")) or 10) * s)))
                except tk.TclError:
                    pass

        # ---- layout ----
        def _build(self):
            pad = {"padx": 8, "pady": 4}
            frm_mode = ttk.LabelFrame(self, text="Input")
            frm_mode.pack(fill="x", **pad)
            self.var_mode = tk.StringVar(value="batch")
            ttk.Radiobutton(
                frm_mode, variable=self.var_mode, value="batch", command=self._on_mode,
                text="Directory of DICOM folders → new output database",
            ).pack(anchor="w", padx=10, pady=(6, 2))
            ttk.Radiobutton(
                frm_mode, variable=self.var_mode, value="single", command=self._on_mode,
                text="Single DICOM folder → append to an existing database (or start a new one)",
            ).pack(anchor="w", padx=10, pady=(2, 4))

            row = ttk.Frame(frm_mode)
            row.pack(fill="x", padx=6, pady=(2, 6))
            self.lbl_in = ttk.Label(row, text="Input:", width=8)
            self.lbl_in.pack(side="left")
            self.var_input = tk.StringVar()
            ttk.Entry(row, textvariable=self.var_input).pack(side="left", fill="x", expand=True)
            ttk.Button(row, text="Browse…", command=self._browse_input).pack(side="left", padx=4)

            frm_out = ttk.LabelFrame(self, text="Output database")
            frm_out.pack(fill="x", **pad)
            row = ttk.Frame(frm_out)
            row.pack(fill="x", padx=6, pady=(6, 2))
            ttk.Label(row, text="Output:", width=8).pack(side="left")
            self.var_output = tk.StringVar()
            self.var_output.trace_add("write", lambda *a: self._update_db_info())
            ttk.Entry(row, textvariable=self.var_output).pack(side="left", fill="x", expand=True)
            ttk.Button(row, text="Browse…", command=self._browse_output).pack(side="left", padx=4)

            row = ttk.Frame(frm_out)
            row.pack(fill="x", padx=6, pady=(2, 2))
            ttk.Label(row, text="Start ID:", width=8).pack(side="left")
            self.var_start = tk.StringVar(value="0001")
            self.ent_start = ttk.Entry(row, textvariable=self.var_start, width=10)
            self.ent_start.pack(side="left")
            ttk.Label(row, text="(used only for a new database)", foreground="gray").pack(side="left", padx=6)

            self.lbl_db = ttk.Label(frm_out, text="", foreground="#2060a0", justify="left")
            self.lbl_db.pack(anchor="w", padx=10, pady=(2, 2))
            ttk.Label(
                frm_out, foreground="gray", justify="left",
                text="CT / PET detected from the DICOM Modality tag (folder names are ignored).\n"
                     "Output: {ID}/{date}{rand}/{ID}_{date}{rand}_CT|PET|SUV.nii.gz + metadata.csv",
            ).pack(anchor="w", padx=10, pady=(0, 6))

            frm_btn = ttk.Frame(self)
            frm_btn.pack(fill="x", **pad)
            self.btn_run = ttk.Button(frm_btn, text="▶ Convert", command=self._run, width=14)
            self.btn_run.pack(side="left")
            self.btn_cancel = ttk.Button(frm_btn, text="■ Cancel", command=self._cancel,
                                         width=12, state="disabled")
            self.btn_cancel.pack(side="left", padx=4)
            ttk.Button(frm_btn, text="Clear log", command=self._clear_log).pack(side="left")
            self.lbl_deps = ttk.Label(frm_btn, text="", foreground="red")
            self.lbl_deps.pack(side="right")

            self.progress = ttk.Progressbar(self, mode="determinate")
            self.progress.pack(fill="x", padx=8, pady=(2, 0))
            self.lbl_status = ttk.Label(self, text="Ready", anchor="w")
            self.lbl_status.pack(fill="x", padx=10)

            frm_log = ttk.LabelFrame(self, text="Log")
            frm_log.pack(fill="both", expand=True, padx=8, pady=(2, 8))
            self.log_text = scrolledtext.ScrolledText(
                frm_log, height=10, state="disabled", font="TkFixedFont",
                background="#1e1e1e", foreground="#d4d4d4",
            )
            self.log_text.pack(fill="both", expand=True, padx=4, pady=4)

        # ---- events ----
        def _on_mode(self):
            if self.var_mode.get() == "batch":
                self.lbl_in.config(text="Input dir:")
            else:
                self.lbl_in.config(text="DICOM:")
            self._update_db_info()

        def _update_db_info(self):
            out = self.var_output.get().strip()
            if not out:
                self.lbl_db.config(text="")
                self.ent_start.config(state="normal")
                return
            if os.path.isfile(os.path.join(out, METADATA_NAME)):
                try:
                    n = len(load_metadata(out)[0])
                except Exception:
                    n = "?"
                self.lbl_db.config(text=f"metadata.csv found ({n} exams) → new data will be APPENDED")
                self.ent_start.config(state="disabled")
            else:
                self.lbl_db.config(text="No metadata.csv → a NEW database will be created")
                self.ent_start.config(state="normal")

        def _ask_dir(self, title, initial):
            d = _native_browse_dir(title, initial)
            if d is None:
                d = filedialog.askdirectory(title=title, initialdir=initial or "")
            return d or ""

        def _browse_input(self):
            title = ("Select directory containing DICOM folders" if self.var_mode.get() == "batch"
                     else "Select one DICOM folder (with CT and PET series)")
            d = self._ask_dir(title, self.var_input.get().strip())
            if d:
                self.var_input.set(d)
                if self.var_mode.get() == "batch" and not self.var_output.get().strip():
                    a = os.path.abspath(d)
                    self.var_output.set(os.path.join(os.path.dirname(a), os.path.basename(a) + "_nii"))

        def _browse_output(self):
            d = self._ask_dir("Select output database folder", self.var_output.get().strip())
            if d:
                self.var_output.set(d)

        def _run(self):
            if self._running:
                return
            inp = self.var_input.get().strip()
            out = self.var_output.get().strip()
            mode = self.var_mode.get()
            if not inp or not os.path.isdir(inp):
                messagebox.showerror("Error", "Please select a valid input folder.")
                return
            if not out:
                messagebox.showerror("Error", "Please select an output folder.")
                return
            if os.path.abspath(out).startswith(os.path.abspath(inp) + os.sep) or \
                    os.path.abspath(out) == os.path.abspath(inp):
                messagebox.showerror("Error", "The output folder must not be inside the input folder.")
                return
            miss = missing_dependencies()
            if miss:
                messagebox.showerror("Missing dependencies", "Cannot start. Missing: " + ", ".join(miss))
                return
            has_db = os.path.isfile(os.path.join(out, METADATA_NAME))
            if mode == "batch" and has_db:
                if not messagebox.askyesno(
                        "Existing database",
                        "The output folder already contains metadata.csv.\n\n"
                        "Append these DICOM folders to that database?"):
                    return
            elif not has_db and os.path.isdir(out) and os.listdir(out):
                if not messagebox.askyesno(
                        "Folder not empty",
                        "The output folder is not empty and has no metadata.csv.\n\n"
                        "Create a new database inside it anyway?"):
                    return
            if mode == "single":
                ct, pet, _ = find_ct_pet(inp, log_fn=lambda m: None)
                if not ct or not pet:
                    messagebox.showerror(
                        "Not a DICOM folder",
                        "Could not find both a CT and a PET series (Modality CT / PT) in the selected folder.\n"
                        "To convert a directory of DICOM folders, choose the first option.")
                    return

            exam_dirs = list_exam_dirs(inp, mode)
            start = self.var_start.get().strip() or "0001"
            if not start.isdigit():
                messagebox.showerror("Error", "Start ID must be an integer, e.g. 0001.")
                return

            self._cancel_event.clear()
            self._running = True
            self.btn_run.config(state="disabled")
            self.btn_cancel.config(state="normal")
            self.progress["value"] = 0
            self.lbl_status.config(text="Running…")

            def worker():
                try:
                    convert(exam_dirs, out, start_id=start,
                            log_fn=lambda m: self.after(0, self._log, m),
                            progress_fn=lambda v: self.after(0, self._set_progress, v),
                            cancel_event=self._cancel_event)
                except Exception as e:
                    msg = f"[FATAL] {e}\n{traceback.format_exc()}"
                    self.after(0, self._log, msg)
                finally:
                    self.after(0, self._done)

            threading.Thread(target=worker, daemon=True).start()

        def _done(self):
            self._running = False
            self.btn_run.config(state="normal")
            self.btn_cancel.config(state="disabled")
            self.lbl_status.config(text="Done")
            self._update_db_info()

        def _cancel(self):
            self._cancel_event.set()
            self.lbl_status.config(text="Cancelling after the current exam…")

        def _set_progress(self, v):
            self.progress["value"] = v

        def _log(self, msg):
            self.log_text.config(state="normal")
            self.log_text.insert("end", msg + "\n")
            self.log_text.see("end")
            self.log_text.config(state="disabled")

        def _clear_log(self):
            self.log_text.config(state="normal")
            self.log_text.delete("1.0", "end")
            self.log_text.config(state="disabled")

    App().mainloop()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        run_gui()
        return 0
    if argv == ["--check-deps"]:
        miss = missing_dependencies()
        print("dcm2niix:", find_dcm2niix())
        print("missing:", ", ".join(miss) if miss else "none")
        return 1 if miss else 0
    ap = argparse.ArgumentParser(description="DICOM PET/CT folder → NIfTI database")
    ap.add_argument("mode", choices=["batch", "single"],
                    help="batch: input holds many DICOM folders; single: input is one DICOM folder")
    ap.add_argument("input")
    ap.add_argument("output", help="database folder (appended to if it has metadata.csv)")
    ap.add_argument("--start-id", default="0001", help="first PID for a new database")
    args = ap.parse_args(argv)
    miss = missing_dependencies()
    if miss:
        print("Missing dependencies: " + ", ".join(miss), file=sys.stderr)
        return 2
    stats = convert(list_exam_dirs(args.input, args.mode), args.output, start_id=args.start_id)
    return 0 if stats["failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
