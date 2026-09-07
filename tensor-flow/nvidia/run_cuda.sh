#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/pinn_tf_cuda_env"
PYTHON_BIN="$VENV_DIR/bin/python"
SITE_PACKAGES="$VENV_DIR/lib/python3.12/site-packages"

if [ ! -x "$PYTHON_BIN" ]; then
    echo "Erro: ambiente virtual não encontrado em $VENV_DIR" >&2
    exit 1
fi

CUDA_LIBRARY_PATH="$(find "$SITE_PACKAGES/nvidia" -type d -name lib -printf '%p:' | sed 's/:$//')"
export LD_LIBRARY_PATH="$CUDA_LIBRARY_PATH${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export CUDA_VISIBLE_DEVICES=0

exec "$PYTHON_BIN" "$@"