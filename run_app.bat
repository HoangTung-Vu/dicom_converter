@echo off
REM Launch DICOM PET/CT → NIfTI Converter (Windows)
REM Usage: double-click or run from cmd

setlocal
set ENV_NAME=dicom_converter
set SCRIPT_DIR=%~dp0

REM ---- Locate conda ----
where conda >nul 2>&1
if %errorlevel% neq 0 (
    if exist "%USERPROFILE%\miniconda3\Scripts\conda.exe" (
        set CONDA_BASE=%USERPROFILE%\miniconda3
    ) else if exist "%USERPROFILE%\anaconda3\Scripts\conda.exe" (
        set CONDA_BASE=%USERPROFILE%\anaconda3
    ) else (
        echo [ERROR] conda not found. Install Miniconda first.
        pause
        exit /b 1
    )
    call "%CONDA_BASE%\Scripts\activate.bat"
)

REM ---- Create env if missing ----
conda env list | findstr /C:"%ENV_NAME%" >nul 2>&1
if %errorlevel% neq 0 (
    echo [INFO] Creating conda environment "%ENV_NAME%"...
    conda env create -f "%SCRIPT_DIR%environment.yml"
)

REM ---- Activate and run ----
call conda activate %ENV_NAME%
python "%SCRIPT_DIR%dicom_petct_tool.py"

endlocal
