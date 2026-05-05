@echo off
setlocal

echo ======================================================
echo Building DICOM PET/CT Converter EXE (Windows)
echo ======================================================

:: 1. Check for conda
where conda >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo [ERROR] Conda not found. Please run this from 'Anaconda Prompt'.
    pause
    exit /b 1
)

:: 2. Activate environment
echo [1/4] Activating environment 'dicom_converter'...
call conda activate dicom_converter
if %ERRORLEVEL% neq 0 (
    echo [INFO] Environment 'dicom_converter' not found. Creating it...
    call conda env create -f environment.yml
    call conda activate dicom_converter
)

:: 3. Install PyInstaller if missing
echo [2/4] Ensuring PyInstaller is installed...
pip show pyinstaller >nul 2>nul
if %ERRORLEVEL% neq 0 (
    pip install pyinstaller
)

:: 4. Build EXE
echo [3/4] Building EXE with PyInstaller...
:: We bundle dcm2niix.exe if it exists in the current directory.
:: You should download dcm2niix.exe from:
:: https://github.com/rordenlab/dcm2niix/releases
if exist dcm2niix.exe (
    echo [INFO] Found dcm2niix.exe, bundling it into the EXE.
    set ADD_DATA=--add-binary "dcm2niix.exe;."
) else (
    echo [WARN] dcm2niix.exe not found in current directory. 
    echo [WARN] It will not be bundled. The app will search for it in PATH.
    set ADD_DATA=
)

pyinstaller ^
  --onefile ^
  --windowed ^
  --name DicomPetCtConverter ^
  %ADD_DATA% ^
  dicom_petct_tool_end2end.py

if %ERRORLEVEL% equ 0 (
    echo ======================================================
    echo [4/4] SUCCESS! 
    echo Output: dist\DicomPetCtConverter.exe
    echo ======================================================
) else (
    echo [ERROR] Build failed.
)

pause
