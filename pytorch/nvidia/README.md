# PINNs NVIDIA

Este diretório usa `pyproject.toml` como fonte das dependências Python e `uv`
para criar ambientes reproduzíveis. O solver SciPy permanece em
`pdes/scipy` como referência independente. Os solvers modernos estão em
`pdes/fenicx` e usam `dolfinx`, NumPy e Gmsh.

## Ambientes

`uv` administra as dependências Python do projeto, incluindo PyTorch/CUDA,
NumPy, SciPy, Matplotlib e Gmsh. O `dolfinx` é uma dependência nativa com
MPI/PETSc: neste projeto ele é instalado dentro do ambiente Conda-forge
isolado, nunca pelo Python do sistema.

Também é possível instalar tudo diretamente pelo Conda. Esse é o fluxo
principal recomendado quando as versões disponíveis nos canais atendem ao
projeto; o `uv` permanece opcional para ambientes Python alternativos.

Crie o ambiente nativo FEniCSx dentro do fluxo do projeto:

```bash
cd pytorch/nvidia
# O VS Code pode injetar PYTHONPATH/PYTHONSTARTUP do Python do sistema.
env -u PYTHONPATH -u PYTHONSTARTUP -u CONDA_NO_PLUGINS -u VIRTUAL_ENV \
  conda env create -f environment-fenicx.yml
conda activate pinns-fenicx
```

Se `conda` não estiver no `PATH`, use o caminho do Miniforge:

```bash
env -u PYTHONPATH -u PYTHONSTARTUP -u CONDA_NO_PLUGINS -u VIRTUAL_ENV \
  /home/g1ll/miniforge3/bin/conda env create -f environment-fenicx.yml
```

Não use `CONDA_NO_PLUGINS=true` neste caso: o Conda 26.7 está configurado
para `libmamba`, que é carregado como plugin. Se precisar desativar plugins,
use também `CONDA_SOLVER=classic`.

O pacote Conda chama-se `fenics-dolfinx`, mas o import Python é `dolfinx`.
Valide que ele está no ambiente isolado:

```bash
python -c "import dolfinx, mpi4py, gmsh; print(dolfinx.__version__)"
conda list | grep -E 'fenics-dolfinx|petsc|slepc|mpi4py|pytorch|cuda'
```

## Execução

SciPy, sem FEniCSx:

```bash
PYTHONPATH=. "$CONDA_PREFIX/bin/python" \
  pdes/scipy/pde_heat_difusion_fem_cavity_2y.py
```

FEniCSx com malha estruturada NumPy, caso 1 (`alpha=beta=0`):

```bash
PYTHONPATH=. "$CONDA_PREFIX/bin/python" \
  -m pdes.fenicx.pde_heat_difusion_fem_cavity_2y_numpy \
  --case 1 --mesh-ref 30
```

FEniCSx com malha Gmsh, casos 1 ou 2:

```bash
PYTHONPATH=. "$CONDA_PREFIX/bin/python" \
  -m pdes.fenicx.pde_heat_difusion_fem_cavity_2y_gmsh \
  --case 2 --mesh-ref 30
```

Os PNGs, XDMF, CSVs e logs são gravados em `pdes/fenicx/results/cavity/`.
`pdes/fenics_old` contém somente a implementação DOLFIN legada para consulta
histórica e não participa do ambiente novo.