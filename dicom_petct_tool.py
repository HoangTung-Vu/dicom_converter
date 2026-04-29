#!/usr/bin/env python3
"""
DICOM PET/CT → NIfTI Converter

Three-step pipeline:
  Step 1  Group DICOM files by filename prefix (in-place on input dir)
  Step 2  Anonymize (optional) + write metadata.json per patient-exam folder
  Step 3  Convert grouped DICOMs to NIfTI (one file per Step-1 subdir prefix)
          + patient_mapping.csv

Step 3 output (flat under {output_dir}/{ID}/{date}/):
  {ID}_{date}_{PREFIX}.nii.gz       - CT or PET, one per Step-1 subdir
  {ID}_{date}_{PREFIX}_SUV.nii.gz   - SUV (only for PET subdirs)
"""

import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
import sys
import threading
import os
import shutil
import csv
import json
import subprocess
import datetime
import traceback
import re
import tempfile

from typing import Any

import numpy as np

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


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sort_by_instance_number(file_list):
    data = []
    for fpath in file_list:
        try:
            f = pydicom.dcmread(fpath, stop_before_pixels=True)
            n = int(f.InstanceNumber)
        except Exception:
            n = 0
        data.append({"f": fpath, "n": n})
    data.sort(key=lambda x: x["n"])
    return [x["f"] for x in data]


def compute_suv(pet_raw, first_dcm_filepath):
    """Return (suv_array, estimated_bool)."""
    estimated = False
    f = pydicom.dcmread(first_dcm_filepath)

    try:
        weight_grams = float(f.PatientWeight) * 1000
    except Exception:
        weight_grams = 75000
        estimated = True

    try:
        scantime = datetime.datetime.strptime(f.AcquisitionTime, "%H%M%S.%f")
        inj_seq = f.RadiopharmaceuticalInformationSequence[0]
        injection_time = datetime.datetime.strptime(
            inj_seq.RadiopharmaceuticalStartTime, "%H%M%S.%f"
        )
        half_life = float(inj_seq.RadionuclideHalfLife)
        injected_dose = float(inj_seq.RadionuclideTotalDose)
        decay = np.exp(-np.log(2) * (scantime - injection_time).seconds / half_life)
        injected_dose_decay = injected_dose * decay
    except Exception:
        decay = np.exp(-np.log(2) * (1.75 * 3600) / 6588)
        injected_dose_decay = 420_000_000 * decay
        estimated = True

    suv = pet_raw * weight_grams / injected_dose_decay
    return suv, estimated


def run_dcm2niix(src_dir, out_dir, filename_prefix):
    """Returns (nii_path_or_None, success_bool, log_string)."""
    os.makedirs(out_dir, exist_ok=True)
    for f in os.listdir(out_dir):
        if f.endswith(".nii.gz") or f.endswith(".json"):
            os.remove(os.path.join(out_dir, f))

    dcm2niix_bin = (
        shutil.which("dcm2niix")
        or os.path.join(os.path.dirname(sys.executable), "dcm2niix")
    )
    cmd = [
        dcm2niix_bin, "-z", "y",
        "-f", filename_prefix,
        "-o", out_dir,
        src_dir,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    log_str = result.stdout + result.stderr

    nii_files = sorted(
        os.path.join(out_dir, f)
        for f in os.listdir(out_dir)
        if f.endswith(".nii.gz")
    )
    if nii_files:
        nii_files.sort(key=os.path.getsize, reverse=True)
        return nii_files[0], True, log_str
    return None, False, log_str


# Match: <prefix><optional sep><MARKER><optional sep><digits>
# MARKER ∈ {CT, PET, PT}. Longer markers (PET) are listed first so they win.
_SUFFIX_RE = re.compile(
    r'^(?P<prefix>.*?)[_\-]?(?P<marker>PET|CT|PT)[_\-]?(?P<digits>\d+)$',
    re.IGNORECASE,
)
_MODALITY_IN_PREFIX = re.compile(r'PET|CT|PT', re.IGNORECASE)


def extract_prefix(stem):
    """Strip trailing CT###/PT###/PET### (with optional separators).

    Returns the prefix string only if it contains CT, PET, or PT — indicating
    the file belongs to a recognisable imaging series. Returns None otherwise
    so the file lands in others/.
    Empty prefix (e.g. 'pet001.dcm') falls back to the marker name itself.
    """
    m = _SUFFIX_RE.match(stem)
    if not m:
        return None
    prefix = m.group("prefix") or m.group("marker").upper()
    # Discard prefixes that carry no modality hint (localiser, scout, AC map…)
    if not _MODALITY_IN_PREFIX.search(prefix):
        return None
    return prefix

def extract_metadata(dcm_path) -> dict[str, Any]:
    """Read patient/study metadata from a DICOM file."""
    ds = pydicom.dcmread(dcm_path, stop_before_pixels=True)

    def _g(tag, default=""):
        return str(getattr(ds, tag, default)).replace("^^^^", "").strip()

    return {
        "patient_name": _g("PatientName") or "Unknown",
        "patient_id": _g("PatientID"),
        "patient_birthdate": _g("PatientBirthDate"),
        "patient_sex": _g("PatientSex"),
        "patient_age": _g("PatientAge"),
        "patient_weight": _g("PatientWeight"),
        "study_date": _g("StudyDate"),
        "acquisition_date": _g("AcquisitionDate"),
        "study_description": _g("StudyDescription"),
        "referring_physician": _g("ReferringPhysicianName"),
        "performing_physician": _g("PerformingPhysicianName"),
        "reading_physician": _g("NameOfPhysiciansReadingStudy"),
        "operators_name": _g("OperatorsName"),
        "institution_name": _g("InstitutionName"),
    }

_METADATA_FILENAME = "metadata.json"

def detect_modality_from_dir(dir_path):
    """Return 'CT', 'PT', or 'UN' by reading first DICOM's Modality tag.
    Falls back to folder-name heuristic.
    """
    for f in sorted(os.listdir(dir_path)):
        if f.lower().endswith(".dcm"):
            try:
                ds = pydicom.dcmread(
                    os.path.join(dir_path, f), stop_before_pixels=True
                )
                mod = str(getattr(ds, "Modality", "UN")).strip().upper()
                if mod in ("CT", "PT"):
                    return mod
            except Exception:
                pass
            break
    name = os.path.basename(dir_path).upper()
    if name.startswith("PET") or name.startswith("PT"):
        return "PT"
    if name.startswith("CT"):
        return "CT"
    return "UN"

# ---------------------------------------------------------------------------
# Step 1 & 2: Preprocess & Anonymize
# ---------------------------------------------------------------------------

def step1_2_preprocess_and_anonymize(input_dir, mode, anonymize, log_fn=print):
    log_fn(f"\n=== Step 1 & 2: Preprocess & Anonymize (Anonymize={anonymize}) ===")
    input_dir_path = os.path.abspath(input_dir)
    
    if mode == "single":
        mapping_csv = os.path.join(os.path.dirname(input_dir_path), "mapping.csv")
        exam_dirs = [input_dir_path]
    else:
        mapping_csv = os.path.join(input_dir_path, "mapping.csv")
        exam_dirs = sorted(
            os.path.join(input_dir_path, d) for d in os.listdir(input_dir_path) 
            if os.path.isdir(os.path.join(input_dir_path, d))
        )

    if not exam_dirs:
        log_fn("  [ERROR] No exam folders found")
        return None

    patient_mapping = {}
    max_id = 0
    if os.path.exists(mapping_csv):
        try:
            with open(mapping_csv, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    pname = row.get("PatientName", "")
                    pid = row.get("PID", "")
                    if pname and pid:
                        patient_mapping[pname] = pid
                        try:
                            max_id = max(max_id, int(pid))
                        except ValueError:
                            pass
        except Exception as e:
            log_fn(f"  [WARN] Could not read {mapping_csv}: {e}")

    num_digits = 4
    if mode == "batch" and not os.path.exists(mapping_csv) and anonymize:
        num_digits = max(4, len(str(len(exam_dirs) * 100)))
    elif max_id > 0:
        num_digits = max(4, len(str(max_id)))

    movement_plan = []
    new_patients = {}
    new_metadata_rows = []
    target_bases_set = set()
    
    for folder in exam_dirs:
        if os.path.exists(os.path.join(folder, _METADATA_FILENAME)):
            log_fn(f"\n[SKIP] {os.path.basename(folder)} (Already processed)")
            continue
            
        log_fn(f"\n[Scanning {os.path.basename(folder)}]")
        dicom_files = []
        for root, _, files in os.walk(folder):
            for f in files:
                if f.lower().endswith(".dcm"):
                    dicom_files.append(os.path.join(root, f))
                    
        if not dicom_files:
            log_fn("  [SKIP] No .dcm files found")
            continue
            
        prefix_groups = {}
        for f in dicom_files:
            stem = os.path.splitext(os.path.basename(f))[0]
            prefix = extract_prefix(stem)
            if prefix:
                prefix_groups.setdefault(prefix, []).append(f)
            else:
                prefix_groups.setdefault("others", []).append(f)
                
        first_dcm = None
        for pref, files in prefix_groups.items():
            if pref != "others":
                first_dcm = files[0]
                break
        if not first_dcm:
            first_dcm = dicom_files[0]
            
        meta = extract_metadata(first_dcm)
        pname = meta.get("patient_name") or "Unknown"
        sdate = meta.get("acquisition_date") or meta.get("study_date") or "00000000"
        
        if not anonymize:
            target_pid = re.sub(r'[<>:"/\\|?*]', "_", pname).strip() or "Unknown"
        else:
            if pname in patient_mapping:
                target_pid = patient_mapping[pname]
            elif pname in patient_mapping.values():
                target_pid = pname
            elif pname in new_patients:
                target_pid = new_patients[pname]
            else:
                max_id += 1
                target_pid = str(max_id).zfill(num_digits)
                new_patients[pname] = target_pid
                
                new_metadata_rows.append({
                    "PID": target_pid,
                    "PatientName": pname,
                    "PatientBirthDate": meta.get("patient_birthdate", ""),
                    "PatientSex": meta.get("patient_sex", ""),
                    "AcquisitionDate": meta.get("acquisition_date", ""),
                    "StudyDate": meta.get("study_date", ""),
                    "ReferringPhysician": meta.get("referring_physician", ""),
                    "PerformingPhysician": meta.get("performing_physician", ""),
                    "InstitutionName": meta.get("institution_name", ""),
                })
                
        target_base = os.path.join(os.path.dirname(input_dir_path) if mode == "single" else input_dir_path, target_pid)
        target_bases_set.add(target_base)
        
        for pref, files in prefix_groups.items():
            target_dir = os.path.join(target_base, sdate, pref)
            for f in files:
                dest = os.path.join(target_dir, os.path.basename(f))
                if os.path.abspath(f) != os.path.abspath(dest):
                    movement_plan.append((f, dest))
                    
    if movement_plan:
        log_fn(f"\nMoving {len(movement_plan)} files into target structures...")
        for src, dst in movement_plan:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.move(src, dst)
            
    for folder in exam_dirs:
        if os.path.exists(folder):
            for root, dirs, files in os.walk(folder, topdown=False):
                for name in dirs:
                    d = os.path.join(root, name)
                    try:
                        if not os.listdir(d):
                            os.rmdir(d)
                    except: pass
            try:
                if not os.listdir(folder):
                    os.rmdir(folder)
            except: pass

    for tbase in sorted(list(target_bases_set)):
        log_fn(f"\n[Processing metadata & anonymization for {os.path.basename(tbase)}]")
        first_dcm = None
        for root, _, files in os.walk(tbase):
            for f in files:
                if f.lower().endswith(".dcm"):
                    first_dcm = os.path.join(root, f)
                    break
            if first_dcm: break
            
        if not first_dcm:
            continue
            
        meta = extract_metadata(first_dcm)
        meta["anonymized"] = bool(anonymize)
        if anonymize:
            meta["assigned_pid"] = os.path.basename(tbase)
            
        meta_path = os.path.join(tbase, _METADATA_FILENAME)
        with open(meta_path, "w", encoding="utf-8") as fh:
            json.dump(meta, fh, indent=2, ensure_ascii=False, default=str)
        log_fn(f"  metadata.json → {meta_path}")
        
        if anonymize:
            pid_str = meta["assigned_pid"]
            count = 0
            for root, _, files in os.walk(tbase):
                for f in files:
                    if f.lower().endswith(".dcm"):
                        fpath = os.path.join(root, f)
                        try:
                            ds = pydicom.dcmread(fpath)
                            clear_tags = [
                                "InstitutionName", "ReferringPhysicianName", "PerformingPhysicianName",
                                "NameOfPhysiciansReadingStudy", "OperatorsName", "PatientBirthDate",
                                "PatientAddress", "PatientTelephoneNumbers",
                            ]
                            for tag in clear_tags:
                                if hasattr(ds, tag):
                                    setattr(ds, tag, "")
                            if hasattr(ds, "PatientName"): ds.PatientName = pid_str
                            if hasattr(ds, "PatientID"): ds.PatientID = pid_str
                            ds.save_as(fpath)
                            count += 1
                        except Exception as e:
                            log_fn(f"    [WARN] Cannot anonymize {os.path.basename(f)}: {e}")
            log_fn(f"  Anonymized {count} DICOM files in {os.path.basename(tbase)}")

    if anonymize and new_metadata_rows:
        file_exists = os.path.exists(mapping_csv)
        fieldnames = ["PID", "PatientName", "PatientBirthDate", "PatientSex", "AcquisitionDate", "StudyDate", "ReferringPhysician", "PerformingPhysician", "InstitutionName"]
        with open(mapping_csv, "a", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            if not file_exists:
                writer.writeheader()
            writer.writerows(new_metadata_rows)
        log_fn(f"\n  CSV mapping updated → {mapping_csv}")

    log_fn("\n✓ Step 1 & 2 done.")
    
    if mode == "single" and len(target_bases_set) == 1:
        return list(target_bases_set)[0]
    return None

# ---------------------------------------------------------------------------
# Step 3: Convert to NIfTI (CT, PET, SUV)
# ---------------------------------------------------------------------------

def step3_convert_to_nifti(
    input_dir, output_dir, mode, log_fn=print, progress_fn=None, cancel_event=None,
):
    """Step 3 entry point. Converts NIfTI to a separate output directory."""
    log_fn(f"\n=== Step 3: Convert to NIfTI → {output_dir} ===")
    input_dir_path = os.path.abspath(input_dir)
    
    if mode == "single":
        patient_dirs = [input_dir_path]
        mapping_csv = os.path.join(os.path.dirname(input_dir_path), "mapping.csv")
    else:
        patient_dirs = [os.path.join(input_dir_path, d) for d in sorted(os.listdir(input_dir_path)) 
                        if os.path.isdir(os.path.join(input_dir_path, d))]
        mapping_csv = os.path.join(input_dir_path, "mapping.csv")

    if not patient_dirs:
        log_fn("  [ERROR] No patient folders found")
        return

    os.makedirs(output_dir, exist_ok=True)
    if os.path.exists(mapping_csv):
        shutil.copy2(mapping_csv, os.path.join(output_dir, "mapping.csv"))
        log_fn(f"  Copied mapping.csv to output")

    total = len(patient_dirs)
    log_fn(f"  Found {total} patient folder(s)")

    for idx, pdir in enumerate(patient_dirs, 1):
        if cancel_event and cancel_event.is_set():
            log_fn("\nCancelled.")
            break
        log_fn(f"\n[{idx}/{total}] {os.path.basename(pdir)}")
        try:
            pid = os.path.basename(pdir)
            
            meta_src = os.path.join(pdir, "metadata.json")
            if os.path.exists(meta_src):
                out_pid_path = os.path.join(output_dir, pid)
                os.makedirs(out_pid_path, exist_ok=True)
                shutil.copy2(meta_src, os.path.join(out_pid_path, "metadata.json"))
                
            for date_sub in sorted(os.listdir(pdir)):
                date_path = os.path.join(pdir, date_sub)
                if not os.path.isdir(date_path): continue
                
                out_date_path = os.path.join(output_dir, pid, date_sub)
                os.makedirs(out_date_path, exist_ok=True)
                
                for prefix_sub in sorted(os.listdir(date_path)):
                    if prefix_sub == "others": continue
                    prefix_path = os.path.join(date_path, prefix_sub)
                    if not os.path.isdir(prefix_path): continue
                    
                    modality = detect_modality_from_dir(prefix_path)
                    if modality not in ("CT", "PT"):
                        log_fn(f"  [SKIP] {prefix_sub}: unknown modality")
                        continue

                    safe_prefix = re.sub(r'[<>:"/\\|?*]', "_", prefix_sub).strip() or prefix_sub
                    name_prefix = f"{pid}_{date_sub}"
                    nii_out = os.path.join(out_date_path, f"{name_prefix}_{safe_prefix}.nii.gz")

                    if os.path.exists(nii_out):
                        log_fn(f"  Skip (exists): {os.path.basename(nii_out)}")
                    else:
                        with tempfile.TemporaryDirectory() as tmp:
                            tmp_out = os.path.join(tmp, safe_prefix)
                            nii_path, ok, msg = run_dcm2niix(prefix_path, tmp_out, safe_prefix)
                            if nii_path and ok:
                                shutil.move(nii_path, nii_out)
                                log_fn(f"  {modality}  → {os.path.basename(nii_out)}")
                            else:
                                log_fn(f"  [WARN] {prefix_sub} conversion failed: {msg[:200]}")
                                continue

                    if modality == "PT":
                        suv_out = os.path.join(out_date_path, f"{name_prefix}_{safe_prefix}_SUV.nii.gz")
                        if not os.path.exists(suv_out):
                            dcms = sort_by_instance_number([
                                os.path.join(prefix_path, f) for f in os.listdir(prefix_path)
                                if f.lower().endswith(".dcm")
                            ])
                            if dcms:
                                try:
                                    pet_img = nib.load(nii_out)
                                    pet_raw = pet_img.get_fdata()
                                    suv, estimated = compute_suv(pet_raw, dcms[0])
                                    nib.save(
                                        nib.Nifti1Image(suv.astype(np.float32), pet_img.affine, pet_img.header),
                                        suv_out,
                                    )
                                    note = " [estimated params]" if estimated else ""
                                    log_fn(f"  SUV → {os.path.basename(suv_out)}{note}")
                                except Exception as e:
                                    log_fn(f"  [WARN] SUV computation failed for {prefix_sub}: {e}")

        except Exception as e:
            log_fn(f"  [ERROR] {e}")
            traceback.print_exc()
        if progress_fn:
            progress_fn(idx / total * 100)

    log_fn("\n✓ Step 3 done.")


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

def _native_browse_dir(title, initial=None):
    """Use zenity / kdialog on Linux for a sane folder picker.

    Tk's askdirectory on Linux uses an internal dialog where folders open on
    double-click only and many themes ignore that — meaning users often can't
    descend into subdirectories at all. The native GTK/Qt pickers behave
    correctly. Returns selected path or None (None = fall back / cancelled).
    """
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
    if r.returncode == 0:
        out = r.stdout.strip()
        return out or None
    return None  # user cancelled (treat as no selection)


def _enable_dpi_awareness():
    """Tell Windows we handle DPI ourselves so it doesn't bitmap-stretch us."""
    if sys.platform != "win32":
        return
    try:
        from ctypes import windll
        # Per-monitor v2 if available, else system-DPI aware.
        try:
            windll.shcore.SetProcessDpiAwareness(2)
        except Exception:
            windll.user32.SetProcessDPIAware()
    except Exception:
        pass


class App(tk.Tk):
    def __init__(self):
        _enable_dpi_awareness()
        super().__init__()
        self.title("DICOM PET/CT → NIfTI Converter")
        self._apply_dpi_scaling()
        # Base UI size; multiplied by detected scale so the floor grows on
        # HiDPI but the window also stays small on standard DPI.
        s = getattr(self, "_ui_scale", 1.0)
        self.minsize(int(560 * s), int(560 * s))
        self._cancel_event = threading.Event()
        self._running = False
        self._build_ui()
        self._check_deps()

    # Hard cap: even on 4K Xft.dpi=192 setups, raw 2× makes the window
    # taller than the screen. Past 1.4× the layout doesn't fit comfortably.
    _MAX_SCALE = 1.4

    def _detect_scale_factor(self):
        """Best-effort UI scale factor (1.0 = no boost) for the host display.

        Order of preference: explicit env var → Linux Xft.dpi via xrdb →
        physical screen mm → Tk's winfo_fpixels.
        Capped at _MAX_SCALE except when explicitly overridden.
        """
        # Manual override (always wins, no cap)
        for var in ("DICOM_TOOL_SCALE", "TK_SCALE", "GDK_SCALE"):
            v = os.environ.get(var)
            if v:
                try:
                    return max(0.5, float(v))
                except ValueError:
                    pass

        raw = 1.0
        # Linux: Xft.dpi (set by most desktop environments on HiDPI)
        if sys.platform.startswith("linux") and shutil.which("xrdb"):
            try:
                out = subprocess.check_output(
                    ["xrdb", "-query"], text=True, timeout=2,
                )
                for line in out.splitlines():
                    if line.startswith("Xft.dpi:"):
                        dpi = float(line.split(":", 1)[1].strip())
                        if dpi > 0:
                            raw = max(raw, dpi / 96.0)
                        break
            except Exception:
                pass

        if raw <= 1.0:
            try:
                mm = self.winfo_screenmmwidth()
                px = self.winfo_screenwidth()
                if mm and mm > 0:
                    dpi = px / (mm / 25.4)
                    if dpi > 110:
                        raw = max(raw, dpi / 96.0)
            except Exception:
                pass

        if raw <= 1.0:
            try:
                dpi = self.winfo_fpixels("1i")
                if dpi > 110:
                    raw = dpi / 96.0
            except Exception:
                pass

        return min(self._MAX_SCALE, max(1.0, raw))

    def _apply_dpi_scaling(self):
        """Scale Tk + named fonts to match the actual screen DPI."""
        boost = self._detect_scale_factor()
        # Tk's "scaling" sets points→pixels (default 96/72 ≈ 1.333).
        self.tk.call("tk", "scaling", (96.0 / 72.0) * boost)

        from tkinter import font as tkfont
        for name in (
            "TkDefaultFont", "TkTextFont", "TkMenuFont",
            "TkHeadingFont", "TkCaptionFont", "TkSmallCaptionFont",
            "TkIconFont", "TkTooltipFont",
        ):
            try:
                f = tkfont.nametofont(name)
                base = abs(f.cget("size")) or 10
                f.configure(size=int(round(base * boost)))
            except tk.TclError:
                pass
        try:
            fixed = tkfont.nametofont("TkFixedFont")
            base = abs(fixed.cget("size")) or 10
            fixed.configure(size=int(round(base * boost)))
            self._log_font = fixed
        except tk.TclError:
            self._log_font = None

        self._ui_scale = boost

    def _build_ui(self):
        # Mode
        frm_mode = ttk.LabelFrame(self, text="Input Mode")
        frm_mode.pack(fill="x", padx=8, pady=4)
        self.var_mode = tk.StringVar(value="batch")
        ttk.Radiobutton(
            frm_mode,
            text="Batch — input là folder cha; mỗi subfolder = 1 patient-exam",
            variable=self.var_mode, value="batch",
        ).pack(anchor="w", padx=10, pady=(6, 2))
        ttk.Radiobutton(
            frm_mode,
            text="Single — input là 1 folder = 1 patient-exam (chứa .dcm hoặc subdirs đã group)",
            variable=self.var_mode, value="single",
        ).pack(anchor="w", padx=10, pady=(2, 6))

        # Input dir
        frm_in = ttk.LabelFrame(self, text="Input Directory")
        frm_in.pack(fill="x", padx=8, pady=4)
        self.var_input = tk.StringVar()
        ttk.Entry(frm_in, textvariable=self.var_input, width=55).pack(
            side="left", fill="x", expand=True, padx=6, pady=6,
        )
        ttk.Button(frm_in, text="Browse…", command=self._browse_input).pack(
            side="right", padx=6, pady=6,
        )

        # Step 1 & 2
        frm1 = ttk.LabelFrame(self, text="Step 1 & 2 — Preprocess & Anonymize")
        frm1.pack(fill="x", padx=8, pady=4)
        
        self.var_anon = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            frm1,
            text="Anonymize DICOM (xóa PII tags + ghi PatientName/ID = PID — ghi đè in-place)",
            variable=self.var_anon,
        ).pack(anchor="w", padx=10, pady=(6, 2))
        
        ttk.Label(
            frm1,
            text="Gom các file DICOM, phân loại theo prefix, lưu metadata, và tạo mapping.csv.\n"
                 "Files sẽ được gom vào cấu trúc: {PatientName_or_ID}/{StudyDate}/{Prefix}/",
            foreground="gray", justify="left",
        ).pack(anchor="w", padx=10, pady=(6, 2))
        
        self.btn_step1 = ttk.Button(
            frm1, text="▶ Run Step 1 & 2", command=self._run_step1, width=18,
        )
        self.btn_step1.pack(anchor="w", padx=10, pady=(4, 8))

        # Step 3
        frm3 = ttk.LabelFrame(self, text="Step 3 — Convert to NIfTI")
        frm3.pack(fill="x", padx=8, pady=4)
        ttk.Label(
            frm3,
            text="Output: {ID}_{date}_{PREFIX}.nii.gz trong Output Directory.\n"
                 "Tạo ra một cấu trúc giống với input, chỉ thay đổi folder chứa dicom thành file nifti.",
            foreground="gray", justify="left",
        ).pack(anchor="w", padx=10, pady=(6, 2))

        sub3 = ttk.Frame(frm3)
        sub3.pack(fill="x", padx=6, pady=(2, 6))
        ttk.Label(sub3, text="Output:").pack(side="left", padx=(4, 4))
        self.var_output = tk.StringVar()
        ttk.Entry(sub3, textvariable=self.var_output).pack(
            side="left", fill="x", expand=True,
        )
        ttk.Button(sub3, text="Browse…", command=self._browse_output).pack(
            side="left", padx=(4, 0),
        )

        self.btn_step3 = ttk.Button(
            frm3, text="▶ Run Step 3", command=self._run_step3, width=18,
        )
        self.btn_step3.pack(anchor="w", padx=10, pady=(2, 8))

        # Progress + status
        self.progress = ttk.Progressbar(self, length=100, mode="determinate")
        self.progress.pack(fill="x", padx=8, pady=(4, 0))
        self.lbl_status = ttk.Label(self, text="Ready", anchor="w")
        self.lbl_status.pack(fill="x", padx=10)

        # Cancel / clear / deps
        frm_btn = ttk.Frame(self)
        frm_btn.pack(fill="x", padx=8, pady=4)
        self.btn_cancel = ttk.Button(
            frm_btn, text="■ Cancel", command=self._cancel, width=12, state="disabled",
        )
        self.btn_cancel.pack(side="left")
        ttk.Button(frm_btn, text="Clear log", command=self._clear_log).pack(
            side="left", padx=4,
        )
        self.lbl_deps = ttk.Label(frm_btn, text="", foreground="red")
        self.lbl_deps.pack(side="right", padx=4)

        # Log
        frm_log = ttk.LabelFrame(self, text="Log")
        frm_log.pack(fill="both", expand=True, padx=8, pady=(2, 8))
        log_font = getattr(self, "_log_font", None) or ("Monospace", 10)
        self.log_text = scrolledtext.ScrolledText(
            frm_log, height=8, state="disabled",
            font=log_font, background="#1e1e1e", foreground="#d4d4d4",
        )
        self.log_text.pack(fill="both", expand=True, padx=4, pady=4)

    # ---- Helpers ----

    def _check_deps(self):
        missing = []
        if not PYDICOM_OK:
            missing.append("pydicom")
        if not NIBABEL_OK:
            missing.append("nibabel")
        if not shutil.which("dcm2niix"):
            missing.append("dcm2niix")
        if missing:
            self.lbl_deps.config(text=f"Missing: {', '.join(missing)}")
            self._log(f"[WARN] Missing dependencies: {', '.join(missing)}")
            self._log("       Activate the petct/dicom_converter conda env and try again.\n")

    def _ask_dir(self, title, initial=None):
        d = _native_browse_dir(title, initial)
        if d is not None:
            return d
        return filedialog.askdirectory(title=title, initialdir=initial or "")

    def _browse_input(self):
        d = self._ask_dir("Select Input DICOM Directory", self.var_input.get().strip())
        if d:
            self.var_input.set(d)

    def _browse_output(self):
        d = self._ask_dir("Select Output Directory", self.var_output.get().strip())
        if d:
            self.var_output.set(d)

    def _validate_input(self):
        d = self.var_input.get().strip()
        if not d or not os.path.isdir(d):
            messagebox.showerror("Error", "Please select a valid input directory.")
            return None
        return d

    def _set_buttons(self, state):
        self.btn_step1.config(state=state)
        self.btn_step3.config(state=state)

    def _run_in_thread(self, task_fn):
        if self._running:
            return
        self._cancel_event.clear()
        self._running = True
        self._set_buttons("disabled")
        self.btn_cancel.config(state="normal")
        self.progress["value"] = 0
        self.lbl_status.config(text="Running…")

        def worker():
            try:
                task_fn()
            except Exception as e:
                err_msg = f"[FATAL] {e}\n{traceback.format_exc()}"
                self.after(0, lambda m=err_msg: self._log(m))
            finally:
                self._running = False
                self.after(0, lambda: self._set_buttons("normal"))
                self.after(0, lambda: self.btn_cancel.config(state="disabled"))
                self.after(0, lambda: self.lbl_status.config(text="Done"))

        threading.Thread(target=worker, daemon=True).start()

    def _log_async(self, msg):
        self.after(0, lambda m=msg: self._log(m))

    def _progress_async(self, v):
        self.after(0, lambda _v=v: self._set_progress(_v))

    # ---- Actions ----

    def _run_step1(self):
        d = self._validate_input()
        if not d:
            return
        if self.var_anon.get() and not PYDICOM_OK:
            messagebox.showerror("Missing dependency", "pydicom is not installed.")
            return

        def task():
            res = step1_2_preprocess_and_anonymize(
                d, self.var_mode.get(), self.var_anon.get(), log_fn=self._log_async
            )
            if res and self.var_mode.get() == "single":
                self.after(0, lambda r=res: self.var_input.set(r))

        self._run_in_thread(task)

    def _run_step3(self):
        d = self._validate_input()
        if not d:
            return
        out = self.var_output.get().strip()
        if not out:
            messagebox.showerror("Error", "Please select an output directory.")
            return

        missing = []
        if not PYDICOM_OK:
            missing.append("pydicom")
        if not NIBABEL_OK:
            missing.append("nibabel")
        if not shutil.which("dcm2niix"):
            missing.append("dcm2niix")
        if missing:
            messagebox.showerror(
                "Missing dependencies",
                f"Cannot start. Missing: {', '.join(missing)}",
            )
            return

        def task():
            step3_convert_to_nifti(
                d, out, self.var_mode.get(),
                log_fn=self._log_async,
                progress_fn=self._progress_async,
                cancel_event=self._cancel_event,
            )

        self._run_in_thread(task)

    def _cancel(self):
        self._cancel_event.set()
        self.lbl_status.config(text="Cancelling…")

    def _log(self, msg):
        self.log_text.config(state="normal")
        self.log_text.insert("end", msg + "\n")
        self.log_text.see("end")
        self.log_text.config(state="disabled")

    def _clear_log(self):
        self.log_text.config(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.config(state="disabled")

    def _set_progress(self, val):
        self.progress["value"] = val
        self.update_idletasks()


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    app = App()
    app.mainloop()
