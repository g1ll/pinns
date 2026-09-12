"""Shared FEniCSx implementation for the MATLAB double-Y cavity."""

import csv
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "cavity"


def cavity_polygon(h0, l0, s1, h1, l1, h2, l2, alpha, beta):
    alpha = np.deg2rad(alpha)
    beta = np.deg2rad(beta)
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta), np.sin(beta)
    x1, y1 = -l0 / 2.0, 0.0
    x2, y2 = x1, s1 - h1 / (2.0 * ca)
    x3, y3 = -(ca * l1 + l0 / 2.0), l1 * sa + y2
    x4, y4 = x3, y3 + h1 / ca
    x5, y5 = x1, y2 + h1 / ca
    x6, y6 = x1, h0 - h2 / cb
    x7, y7 = -(cb * l2 + l0 / 2.0), h0 + l2 * sb - h2 / cb
    x8, y8 = x7, h0 + l2 * sb
    x9, y9 = x1, h0
    left = np.array([
        (x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5),
        (x6, y6), (x7, y7), (x8, y8), (x9, y9),
    ], dtype=np.float64)
    right = left[::-1].copy()
    right[:, 0] *= -1.0
    return np.vstack((left, right))


def case_dimensions(case):
    if case == 1:
        dofs = (1.0, 20.0, 0.14, 14.0, 0.45, 0.0, 0.0)
        fractions = (0.02, 0.02, 0.02)
    elif case == 2:
        dofs = (1.0, 10.0, 0.07, 0.07, 0.5, 13.6875, 13.6875)
        fractions = (0.015, 0.015, 0.04)
    else:
        raise ValueError("case must be 1 or 2")
    hd_l, h0d_l0, h2d_l2, h1d_l1, s1d_h0, alpha, beta = dofs
    fi1, fi2, fi3 = fractions
    length = np.sqrt(1.0 / hd_l)
    height = hd_l * length
    l0 = np.sqrt(fi3 / h0d_l0)
    h0 = h0d_l0 * l0
    l1 = np.sqrt(fi1 / h1d_l1)
    h1 = h1d_l1 * l1
    l2 = np.sqrt(fi2 / h2d_l2)
    h2 = h2d_l2 * l2
    polygon = cavity_polygon(h0, l0, s1d_h0 * h0, h1, l1, h2, l2, alpha, beta)
    return length, height, polygon


def _point_in_polygon(point, polygon):
    inside = False
    for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
        if (start[1] > point[1]) != (end[1] > point[1]):
            crossing = ((end[0] - start[0]) * (point[1] - start[1]) /
                        (end[1] - start[1]) + start[0])
            if point[0] < crossing:
                inside = not inside
    return inside


def numpy_mesh(length, height, mesh_ref, polygon):
    """Build a conforming structured mesh for orthogonal cavity walls."""
    axis_aligned = all(
        np.isclose(start[0], end[0]) or np.isclose(start[1], end[1])
        for start, end in zip(polygon, np.roll(polygon, -1, axis=0))
    )
    if not axis_aligned:
        raise ValueError("The NumPy mesh requires alpha=beta=0; use Gmsh for case 2.")
    x_values = np.unique(np.r_[np.linspace(-length / 2, length / 2, mesh_ref + 1), polygon[:, 0]])
    y_values = np.unique(np.r_[np.linspace(0.0, height, mesh_ref + 1), polygon[:, 1]])
    nodes = np.asarray([(x, y) for y in y_values for x in x_values], dtype=np.float64)
    row_width = len(x_values)
    cells = []
    for row in range(len(y_values) - 1):
        for column in range(row_width - 1):
            center = np.array([(x_values[column] + x_values[column + 1]) / 2,
                               (y_values[row] + y_values[row + 1]) / 2])
            if _point_in_polygon(center, polygon):
                continue
            lower_left = row * row_width + column
            lower_right = lower_left + 1
            upper_left = lower_left + row_width
            upper_right = upper_left + 1
            cells.extend(((lower_left, lower_right, upper_right),
                          (lower_left, upper_right, upper_left)))
    return nodes, np.asarray(cells, dtype=np.int64)


def mesh_from_numpy(nodes, cells):
    from mpi4py import MPI
    import ufl
    from basix.ufl import element
    from dolfinx import mesh

    coordinate_element = element("Lagrange", "triangle", 1, shape=(2,))
    domain = ufl.Mesh(coordinate_element)
    return mesh.create_mesh(MPI.COMM_SELF, cells, domain, nodes)


def mesh_from_gmsh(length, height, polygon, mesh_ref):
    from mpi4py import MPI
    import gmsh
    from dolfinx.io import gmsh as gmsh_io

    gmsh.initialize()
    try:
        gmsh.model.add("double_y_cavity")
        geometry = gmsh.model.geo
        characteristic_length = min(length, height) / mesh_ref
        outer = [
            geometry.addPoint(-length / 2, 0, 0, characteristic_length),
            geometry.addPoint(length / 2, 0, 0, characteristic_length),
            geometry.addPoint(length / 2, height, 0, characteristic_length),
            geometry.addPoint(-length / 2, height, 0, characteristic_length),
        ]
        cavity = [geometry.addPoint(float(x), float(y), 0, characteristic_length)
                  for x, y in polygon]
        boundary = [outer[0], *cavity, outer[1], outer[2], outer[3]]
        lines = [geometry.addLine(boundary[i], boundary[i + 1])
                 for i in range(len(boundary) - 1)]
        lines.append(geometry.addLine(boundary[-1], boundary[0]))
        loop = geometry.addCurveLoop(lines)
        surface = geometry.addPlaneSurface([loop])
        geometry.synchronize()
        gmsh.model.addPhysicalGroup(1, lines, 1)
        gmsh.model.setPhysicalName(1, 1, "boundary")
        gmsh.model.addPhysicalGroup(2, [surface], 2)
        gmsh.model.setPhysicalName(2, 2, "domain")
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(2)
        data = gmsh_io.model_to_mesh(gmsh.model, MPI.COMM_SELF, 0, gdim=2)
        return data.mesh
    finally:
        gmsh.finalize()


def _vertex_values(msh, function):
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim, 0)
    cells = msh.topology.connectivity(tdim, 0).array.reshape(-1, 3)
    coordinates = msh.geometry.x[:, :2]
    points = np.zeros((len(cells) * 3, 3), dtype=coordinates.dtype)
    points[:, :2] = coordinates[cells].reshape(-1, 2)
    cell_indices = np.repeat(np.arange(len(cells), dtype=np.int32), 3)
    sampled = np.asarray(function.eval(points, cell_indices)).reshape(-1)
    values = np.zeros(len(coordinates), dtype=np.float64)
    counts = np.zeros(len(coordinates), dtype=np.int32)
    np.add.at(values, cells.reshape(-1), sampled)
    np.add.at(counts, cells.reshape(-1), 1)
    return coordinates, cells, values / np.maximum(counts, 1)


def save_temperature_plot(msh, function, polygon, output, title):
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    nodes, cells, values = _vertex_values(msh, function)
    maximum = float(np.max(values))
    figure, axis = plt.subplots(figsize=(7, 7))
    triangulation = mtri.Triangulation(nodes[:, 0], nodes[:, 1], cells)
    contour = axis.tricontourf(
        triangulation, values, levels=np.linspace(0, maximum, 31),
        cmap="jet", extend="max",
    )
    closed = np.vstack((polygon, polygon[0]))
    axis.plot(closed[:, 0], closed[:, 1], color="cyan", linewidth=1.0,
              label="Cavidade: T = 0")
    axis.legend(loc="upper right")
    figure.colorbar(contour, ax=axis, label="Temperatura")
    axis.set_title(title)
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_aspect("equal", adjustable="box")
    axis.set_box_aspect(1)
    figure.tight_layout()
    path = output / "resultado_temperatura.png"
    figure.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(figure)
    return path, maximum, len(nodes), len(cells)


def solve(msh, polygon, name):
    from dolfinx import default_scalar_type, fem, io
    from dolfinx.fem.petsc import LinearProblem
    import ufl

    space = fem.functionspace(msh, ("Lagrange", 1))
    facets = mesh_boundary_facets(msh, polygon)
    dofs = fem.locate_dofs_topological(space, msh.topology.dim - 1, facets)
    boundary = fem.dirichletbc(default_scalar_type(0), dofs, space)
    trial = ufl.TrialFunction(space)
    test = ufl.TestFunction(space)
    source = fem.Constant(msh, default_scalar_type(1))
    bilinear = ufl.inner(ufl.grad(trial), ufl.grad(test)) * ufl.dx
    linear = source * test * ufl.dx
    problem = LinearProblem(
        bilinear, linear, bcs=[boundary],
        petsc_options_prefix=f"fenicx_{name}_",
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    )
    temperature = problem.solve()
    output = RESULTS / f"{name}_{datetime.now():%y%m%d_%H%M%S}"
    output.mkdir(parents=True, exist_ok=True)
    with io.XDMFFile(msh.comm, str(output / "temperatura.xdmf"), "w") as file:
        file.write_mesh(msh)
        file.write_function(temperature)
    image, maximum, nodes, elements = save_temperature_plot(
        msh, temperature, polygon, output, f"Distribuicao de temperatura - FEniCSx {name}"
    )
    with (output / "temperatura_maxima.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["nodes", "elements", "temperature_max", "png"])
        writer.writerow([nodes, elements, f"{maximum:.15f}", image.name])
    with (output / "temperatura_maxima.log").open("w", encoding="utf-8") as file:
        file.write(f"nodes={nodes}\n")
        file.write(f"elements={elements}\n")
        file.write(f"temperature_max={maximum:.15f}\n")
        file.write(f"png={image}\n")
    return maximum, output


def mesh_boundary_facets(msh, polygon):
    from dolfinx import mesh

    def on_cavity(x):
        result = np.zeros(x.shape[1], dtype=bool)
        for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
            segment = end - start
            length_squared = np.dot(segment, segment)
            parameter = np.clip(((x[0] - start[0]) * segment[0] +
                                 (x[1] - start[1]) * segment[1]) / length_squared, 0, 1)
            distance = np.sqrt((x[0] - start[0] - parameter * segment[0]) ** 2 +
                               (x[1] - start[1] - parameter * segment[1]) ** 2)
            result |= distance < 1e-8
        return result

    return mesh.locate_entities_boundary(
        msh, msh.topology.dim - 1, on_cavity
    )
