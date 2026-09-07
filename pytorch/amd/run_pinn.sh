#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/pinn_amd_env"
PYTHON_BIN="$VENV_DIR/bin/python"
ROCM_PATH="${ROCM_PATH:-/opt/rocm}"

if [[ ! -x "$PYTHON_BIN" ]]; then
    echo "Erro: ambiente virtual nao encontrado em $VENV_DIR" >&2
    exit 1
fi

if [[ ! -d "$ROCM_PATH" ]]; then
    echo "Erro: ROCm nao encontrado em $ROCM_PATH" >&2
    echo "Defina ROCM_PATH para o diretorio correto do ROCm." >&2
    exit 1
fi

export PATH="$ROCM_PATH/bin:$PATH"
export LD_LIBRARY_PATH="$ROCM_PATH/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# A RX 6600 (gfx1032) usa o override gfx1030 para compatibilidade com ROCm.
export HSA_OVERRIDE_GFX_VERSION="${HSA_OVERRIDE_GFX_VERSION:-10.3.0}"
export HIP_VISIBLE_DEVICES="0"
unset CUDA_VISIBLE_DEVICES

if [[ $# -eq 0 ]]; then
    echo "Uso: ./run_pinn.sh <script.py> [argumentos...]" >&2
    exit 2
fi

exec "$PYTHON_BIN" "$@"