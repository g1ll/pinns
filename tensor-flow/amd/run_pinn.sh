#!/bin/bash
# Script para inicializar o ambiente de treinamento de PINNs em GPU AMD (RX 6600)
# Requisitos: ROCm 7.0.2 + Ubuntu 24.04 (Python 3.12)

# Definir o nome do diretório do ambiente virtual (padrão: pinn_env)
VENV_DIR="pinn_env"

if [ -d "$VENV_DIR" ]; then
    echo "=========================================================="
    echo " Ativando o ambiente virtual Python: $VENV_DIR"
    echo "=========================================================="
    source "$VENV_DIR"/bin/activate
else
    echo "Erro: O diretório do ambiente virtual '$VENV_DIR' não foi encontrado."
    echo "Certifique-se de que o script está na mesma pasta que a pasta do ambiente virtual."
    exit 1
fi

# Configuração de Override para a GPU AMD RX 6600 (Navi 23 -> gfx1032)
# Faz com que as bibliotecas do ROCm simulem uma arquitetura gfx1030 (RX 6800) compatível.
echo "Aplicando HSA_OVERRIDE_GFX_VERSION=10.3.0 para compatibilidade..."
export HSA_OVERRIDE_GFX_VERSION=10.3.0

# Se você possuir duas GPUs (NVIDIA + AMD) e quiser isolar a execução apenas na GPU AMD,
# você pode descomentar a linha abaixo:
# export CUDA_VISIBLE_DEVICES=""

# Verifica se o usuário passou algum arquivo python como argumento
if [ -z "$1" ]; then
    echo "Ambiente pronto para uso! Você pode rodar seus scripts usando a GPU AMD."
    echo "Uso recomendado do script: ./run_pinn.sh <seu_script_pinn.py>"
    echo "----------------------------------------------------------"
    # Mantém o terminal no ambiente virtual ativado
    $SHELL
else
    echo "Executando script: $1..."
    echo "----------------------------------------------------------"
    python3 "$1"
fi
