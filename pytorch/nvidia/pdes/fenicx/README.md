# FEniCSx moderno

`common.py` contém a geometria de 18 pontos copiada das equações de
`matlab/pde_cav_2y.m`, a formulação variacional e a exportação dos resultados.

- `pde_heat_difusion_fem_cavity_2y_numpy.py`: malha estruturada NumPy para o
  caso ortogonal.
- `pde_heat_difusion_fem_cavity_2y_gmsh.py`: malha CAD-conforme Gmsh para os
  dois casos, especialmente o caso inclinado.

Os scripts usam `dolfinx.fem.petsc.LinearProblem`, `ufl` e
`dolfinx.io.gmsh.model_to_mesh`. Cada execução salva `resultado_temperatura.png`,
`temperatura.xdmf`, `temperatura_maxima.csv` e `temperatura_maxima.log`.