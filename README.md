# PINNs com TensorFlow e PyTorch

Exemplos de redes neurais informadas pela física (PINNs) usando TensorFlow e PyTorch em Linux Ubuntu, com suporte para:

- AMD Radeon RX 6600 via ROCm.
- NVIDIA GeForce RTX 3070 via CUDA.

Os ambientes virtuais são separados. Use apenas o ambiente correspondente à GPU escolhida.

## Requisitos gerais

- Ubuntu 24.04 (ou versão compatível).
- Python 3.12.
- `python3.12-venv`, `pip` e ferramentas básicas de compilação.
- Driver da GPU instalado.
- Aproximadamente 10 GB livres para os pacotes de aprendizado de máquina.

Instalação dos requisitos básicos:

```bash
sudo apt update
sudo apt install -y python3.12 python3.12-venv python3-pip build-essential
```

Verifique a GPU antes de configurar o Python:

```bash
# AMD
rocminfo

# NVIDIA
nvidia-smi
```

## Estrutura principal

```text
pytorch/
├── amd/
│   ├── pinn_amd_env/
│   ├── requirements.txt
│   ├── run_pinn.sh
│   ├── odes/ode_heat_decay.py
│   └── tests/find_device.py
└── nvidia/
        ├── pinn_cuda_env/
        └── odes/ode_heat_decay.py
tensor-flow/
├── amd/
│   ├── pinn_tf_env/
│   ├── requirements.txt
│   ├── run_pinn.sh
│   └── tests/find_device.py
└── nvidia/
    ├── pinn_tf_cuda_env/
    ├── requirements.txt
    ├── run_cuda.sh
    └── tests/find_device.py
```

Os diretórios `pinn_tf_env`, `pinn_tf_cuda_env`, `pinn_amd_env` e `pinn_cuda_env` são ambientes virtuais locais e não precisam ser versionados.

## PyTorch

Os arquivos `requirements.txt` instalam pacotes Python, mas não instalam o driver da GPU, o ROCm, o CUDA Toolkit nem recriam o ambiente virtual. Para reinstalar um ambiente removido, crie um novo `venv` e instale as dependências da GPU correspondente.

Não misture os ambientes AMD e NVIDIA. O PyTorch usa o nome `cuda:0` para dispositivos CUDA e também para dispositivos HIP/ROCm; no ambiente AMD, confirme `torch.version.hip` para garantir que o wheel é ROCm.

### PyTorch com AMD Radeon RX 6600 e ROCm 7.2

Antes de criar o ambiente, instale o ROCm 7.2 compatível com o Ubuntu e confirme a placa:

```bash
rocminfo | grep -E "gfx|Marketing Name"
```

A RX 6600 normalmente aparece como `gfx1032`. Neste projeto, o launcher aplica `HSA_OVERRIDE_GFX_VERSION=10.3.0` para compatibilidade.

Recrie o ambiente e instale as dependências:

```bash
cd pytorch/amd
python3.12 -m venv pinn_amd_env
source pinn_amd_env/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt \
    --extra-index-url https://download.pytorch.org/whl/rocm7.2
```

Teste o carregamento do PyTorch e uma operação na GPU:

```bash
./run_pinn.sh tests/find_device.py
```

Execute a PINN usando o launcher, que configura o ROCm e a RX 6600:

```bash
./run_pinn.sh odes/ode_heat_decay.py
```

### PyTorch com NVIDIA GeForce RTX 3070 e CUDA

Verifique primeiro o driver NVIDIA:

```bash
nvidia-smi
```

O ambiente atual usa PyTorch `2.14.0+cu130` com CUDA `13.0`. Como o diretório NVIDIA não possui um `requirements.txt`, instale explicitamente as dependências Python:

```bash
cd pytorch/nvidia
python3.12 -m venv pinn_cuda_env
source pinn_cuda_env/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy matplotlib
python -m pip install torch torchvision torchaudio \
    --index-url https://download.pytorch.org/whl/cu130
```

Confirme que o PyTorch reconhece a RTX 3070:

```bash
python -c 'import torch; print(torch.__version__); print(torch.version.cuda); print(torch.cuda.is_available()); print(torch.cuda.get_device_name(0))'
```

Execute diretamente, sem launcher adicional:

```bash
python odes/ode_heat_decay.py
```

O script NVIDIA exige CUDA e não possui fallback para CPU. A inicialização explícita do contexto CUDA evita o aviso do cuBLAS na primeira operação de `backward`.

## TensorFlow

### AMD Radeon RX 6600

#### 1. Instalar o ROCm

Instale uma versão do ROCm compatível com a sua distribuição. O script usa `/opt/rocm` como caminho padrão.

Depois confirme que o ROCm enxerga a placa:

```bash
rocminfo | grep -E "gfx|Marketing Name"
```

A RX 6600 deve aparecer como `gfx1032`.

Se o pacote `rocprofiler-register` estiver desalinhado com a versão ativa do ROCm, corrija-o com:

```bash
sudo apt update
sudo apt install --reinstall rocprofiler-register
```

#### 2. Criar o ambiente Python

```bash
cd tensor-flow/amd
python3.12 -m venv pinn_tf_env
source pinn_tf_env/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

O arquivo de dependências usa `tensorflow_rocm` e `tf_keras`.

#### 3. Testar a GPU

Execute a partir do diretório `tensor-flow/amd`:

```bash
./run_pinn.sh tests/find_device.py
```

O resultado esperado contém:

```text
GPUs encontradas: [PhysicalDevice(name='/physical_device:GPU:0', device_type='GPU')]
```

O launcher configura:

- `HSA_OVERRIDE_GFX_VERSION=10.3.0` para a RX 6600.
- `HIP_VISIBLE_DEVICES=0`, que é o índice da GPU AMD no ROCm.
- As bibliotecas em `/opt/rocm/lib`.

#### 4. Executar uma PINN

```bash
./run_pinn.sh odes/ode_heat_decay.py
```

O gráfico gerado será salvo no diretório atual como `resultado_ode_heat_decay_<data>.png`.

### NVIDIA GeForce RTX 3070

####  1. Verificar o driver

```bash
nvidia-smi
```

A RTX 3070 deve aparecer na tabela. O driver NVIDIA fornece `libcuda.so`; o TensorFlow também precisa das bibliotecas CUDA e cuDNN de usuário.

#### 2. Criar o ambiente Python

```bash
cd tensor-flow/nvidia
python3.12 -m venv pinn_tf_cuda_env
source pinn_tf_cuda_env/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install --upgrade "tensorflow[and-cuda]"
```

O segundo comando instala as bibliotecas CUDA/cuDNN no próprio ambiente virtual. Ele é necessário mesmo quando o driver NVIDIA já está instalado.

#### 3. Testar a GPU

Use sempre o launcher, pois ele adiciona ao `LD_LIBRARY_PATH` as bibliotecas CUDA instaladas pelo pip antes de iniciar o Python:

```bash
./run_cuda.sh tests/find_device.py
```

O resultado esperado contém:

```text
GPUs encontradas: [PhysicalDevice(name='/physical_device:GPU:0', device_type='GPU')]
```

#### 4. Executar uma PINN

```bash
./run_cuda.sh odes/ode_heat_decay.py
```

O launcher seleciona a NVIDIA como `CUDA_VISIBLE_DEVICES=0` e configura automaticamente os diretórios `nvidia/*/lib` do ambiente virtual.

### Monitorar o uso da GPU

Em outro terminal, durante o treinamento:

```bash
# NVIDIA
watch -n 0.5 nvidia-smi

# AMD
watch -n 0.5 rocm-smi
```

Uma PINN pequena pode apresentar uso baixo ou intermitente da GPU. Isso ocorre quando as operações são pequenas e o custo do loop Python domina. A memória alocada na GPU, sozinha, não prova que o treinamento está usando a GPU; confirme também o log do TensorFlow, que deve indicar `Created device ... GPU:0`.

### Diagnóstico rápido

#### TensorFlow não encontra a AMD

```bash
rocminfo
ls -l /opt/rocm/lib/librocprofiler-register.so.0
```

Se a biblioteca estiver ausente, alinhe o pacote ROCm:

```bash
sudo apt update
sudo apt install --reinstall rocprofiler-register
```

Verifique também se `HIP_VISIBLE_DEVICES=0` está sendo usado.

#### TensorFlow não encontra a NVIDIA

```bash
nvidia-smi
./run_cuda.sh tests/find_device.py
```

Se aparecer `Error loading CUDA libraries`, reinstale as dependências CUDA no ambiente:

```bash
source pinn_tf_cuda_env/bin/activate
python -m pip install --upgrade "tensorflow[and-cuda]"
```

Não execute o teste diretamente com `pinn_tf_cuda_env/bin/python` quando as bibliotecas CUDA estiverem apenas no ambiente pip; use `./run_cuda.sh`.

#### Conferir a versão e o build do TensorFlow

```bash
./run_cuda.sh -c 'import tensorflow as tf; print(tf.__version__); print(tf.sysconfig.get_build_info()); print(tf.config.list_physical_devices("GPU"))'
```

Para AMD, use o Python do ambiente AMD:

```bash
source tensor-flow/amd/pinn_tf_env/bin/activate
python -c 'import tensorflow as tf; print(tf.__version__); print(tf.config.list_physical_devices("GPU"))'
```

### Observações

- Não misture `tensorflow_rocm` e `tensorflow` CUDA no mesmo ambiente virtual.
- Não use o ambiente AMD para executar os scripts NVIDIA, nem o contrário.
- Sempre execute os launchers a partir dos diretórios `tensor-flow/amd` ou `tensor-flow/nvidia`.
- Para resultados reprodutíveis, registre as versões de Python, TensorFlow, ROCm/CUDA e driver da GPU.
