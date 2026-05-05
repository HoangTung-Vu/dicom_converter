#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="dicom_converter"

# ---- Locate conda ----
if command -v conda &>/dev/null; then
    CONDA_BASE=$(conda info --base)
elif [ -f "$HOME/miniconda3/bin/conda" ]; then
    CONDA_BASE="$HOME/miniconda3"
elif [ -f "$HOME/anaconda3/bin/conda" ]; then
    CONDA_BASE="$HOME/anaconda3"
else
    echo "[ERROR] conda not found."
    exit 1
fi

source "$CONDA_BASE/etc/profile.d/conda.sh"

# ---- Create env if missing ----
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "[INFO] Creating conda environment '${ENV_NAME}'..."
    conda env create -f "$SCRIPT_DIR/environment.yml"
fi

conda activate "$ENV_NAME"
python "$SCRIPT_DIR/dicom_petct_tool_end2end.py"
