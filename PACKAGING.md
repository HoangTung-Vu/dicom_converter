# Packaging Instructions

## Prerequisites

```bash
conda env create -f environment.yml
conda activate dicom_converter
pip install pyinstaller
```

---

## Linux — Build a `.deb` package

The `.deb` approach bundles the app as a standalone script + desktop entry.

### Step 1 — Freeze with PyInstaller

```bash
conda activate dicom_converter

pyinstaller \
  --onefile \
  --windowed \
  --name dicom-petct-converter \
  --add-data "dicom_petct_tool.py:." \
  dicom_petct_tool.py
```

Output: `dist/dicom-petct-converter`

### Step 2 — Create .deb package structure

```bash
APP=dicom-petct-converter
VER=1.0.0
ARCH=amd64
PKG="${APP}_${VER}_${ARCH}"

mkdir -p "${PKG}/usr/local/bin"
mkdir -p "${PKG}/usr/share/applications"
mkdir -p "${PKG}/DEBIAN"

# Copy binary
cp dist/${APP} "${PKG}/usr/local/bin/${APP}"
chmod +x "${PKG}/usr/local/bin/${APP}"

# Desktop entry
cat > "${PKG}/usr/share/applications/${APP}.desktop" <<EOF
[Desktop Entry]
Name=DICOM PET/CT Converter
Comment=Convert DICOM PET/CT to NIfTI
Exec=/usr/local/bin/${APP}
Icon=utilities-system-monitor
Terminal=false
Type=Application
Categories=Science;MedicalSoftware;
EOF

# Control file
cat > "${PKG}/DEBIAN/control" <<EOF
Package: ${APP}
Version: ${VER}
Section: science
Priority: optional
Architecture: ${ARCH}
Depends: dcm2niix
Maintainer: Your Name <you@example.com>
Description: DICOM PET/CT to NIfTI converter with anonymization
 GUI tool to convert raw DICOM PET/CT series to organized NIfTI files,
 compute SUV, resample CT/PET to each other's space, and optionally
 anonymize DICOM headers.
EOF

dpkg-deb --build "${PKG}"
echo "Built: ${PKG}.deb"
```

### Step 3 — Install

```bash
sudo dpkg -i dicom-petct-converter_1.0.0_amd64.deb
# Install dcm2niix separately if needed:
sudo apt install dcm2niix
```

---

## Windows — Build a `.exe` with PyInstaller

Run in Anaconda Prompt (Windows):

```bat
# Option 1: Use the automated script
build_windows.bat

# Option 2: Run manually
conda activate dicom_converter
pip install pyinstaller

# Optional: Place dcm2niix.exe in this folder to bundle it
pyinstaller ^
  --onefile ^
  --windowed ^
  --name DicomPetCtConverter ^
  --add-binary "dcm2niix.exe;." ^
  dicom_petct_tool_end2end.py
```

Output: `dist\DicomPetCtConverter.exe`

### Optional: wrap in Inno Setup installer (.exe installer)

1. Download and install [Inno Setup](https://jrsoftware.org/isinfo.php)
2. Create `installer.iss`:

```ini
[Setup]
AppName=DICOM PET/CT Converter
AppVersion=1.0.0
DefaultDirName={autopf}\DicomPetCtConverter
OutputBaseFilename=DicomPetCtConverter_Setup
Compression=lzma
SolidCompression=yes

[Files]
Source: "dist\DicomPetCtConverter.exe"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{autoprograms}\DICOM PET/CT Converter"; Filename: "{app}\DicomPetCtConverter.exe"
Name: "{autodesktop}\DICOM PET/CT Converter"; Filename: "{app}\DicomPetCtConverter.exe"

[Run]
Filename: "{app}\DicomPetCtConverter.exe"; Description: "Launch app"; Flags: nowait postinstall skipifsilent
```

3. Build: `iscc installer.iss`  
   Output: `Output\DicomPetCtConverter_Setup.exe`

> **Note on dcm2niix (Windows):** The frozen `.exe` calls `dcm2niix` as an external process.  
> Download the Windows binary from https://github.com/rordenlab/dcm2niix/releases  
> and place `dcm2niix.exe` alongside the app, or add it to `PATH`.  
> Alternatively, bundle it in the PyInstaller step with `--add-binary "dcm2niix.exe;."` and update  
> the tool to look for it in `sys._MEIPASS` when frozen.

---

## Quick reference

| Command | Purpose |
|---|---|
| `./run_app.sh` | Run on Linux/macOS (auto-creates env) |
| `run_app.bat` | Run on Windows (auto-creates env) |
| `conda env create -f environment.yml` | Create conda env manually |
| `conda activate dicom_converter` | Activate env |
| `python dicom_petct_tool.py` | Run app directly |
