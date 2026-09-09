# FEM FEniCS para a cavidade duplo-Y

Os scripts reproduzem o problema de `matlab/pde_cav_2y.m`:

```text
-div(grad(T)) = 1
T = 0 na cavidade
n.grad(T) = 0 no contorno externo
```

`fenics_common.py` concentra a geometria de 18 pontos e a montagem do
problema variacional. Os executáveis são:

- `pde_heat_difusion_fem_cavity_2y_numpy.py`: malha estruturada feita apenas
  com NumPy. Mantém todos os segmentos da cavidade para `alpha=beta=0`.
- `pde_heat_difusion_fem_cavity_2y_gmsh.py`: malha triangular com a cavidade
  definida como um furo no modelo CAD do Gmsh. É a opção indicada para
  `alpha` ou `beta` diferentes de zero.

## Execução

É necessário ter o FEniCS legado (`dolfin`) e, para o segundo script, o
módulo Python `gmsh` disponíveis no mesmo interpretador:

No Ubuntu 24.04, `dolfin` deve ser instalado pelo APT, pois não existe wheel
pip para Python 3.12:

```bash
sudo apt update
sudo apt install python3-dolfin
python3 -m venv --system-site-packages pytorch/nvidia/pinn_cuda_env
pytorch/nvidia/pinn_cuda_env/bin/python -m pip install gmsh
```

O venv atual foi criado sem `--system-site-packages`; portanto, depois da
instalação APT, ele precisa ser recriado ou substituído por um venv criado
com essa opção para enxergar o módulo `dolfin`.

```bash
PYTHONPATH=pytorch/nvidia/pdes/fenics python \
  pytorch/nvidia/pdes/fenics/pde_heat_difusion_fem_cavity_2y_numpy.py

PYTHONPATH=pytorch/nvidia/pdes/fenics python \
  pytorch/nvidia/pdes/fenics/pde_heat_difusion_fem_cavity_2y_gmsh.py
```

Os resultados são gravados em `fenics/results/cavity/`, em XDMF para abrir no
ParaView e CSV com o número de nós, elementos e temperatura máxima.

O refinamento padrão é `mesh_ref=100`, equivalente a `Hmax=1/100` no MATLAB.
Para testar mais rapidamente, importe `run` e use, por exemplo, `run(30)`.