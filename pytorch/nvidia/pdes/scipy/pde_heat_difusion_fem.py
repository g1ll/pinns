from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


# Problema estacionario: -k * div(grad(T)) = q no retangulo.
# T = 0 nas quatro laterais (condicao de Dirichlet).
Lx = 1.0
Ly = 1.0
k = 10.0
q = 1000.0
nx = 100
ny = 100


def build_mesh(length_x, length_y, elements_x, elements_y):
	x_values = np.linspace(0.0, length_x, elements_x + 1)
	y_values = np.linspace(0.0, length_y, elements_y + 1)
	nodes = np.array(
		[(x, y) for y in y_values for x in x_values], dtype=np.float64
	)

	elements = []
	for row in range(elements_y):
		for column in range(elements_x):
			lower_left = row * (elements_x + 1) + column
			lower_right = lower_left + 1
			upper_left = lower_left + elements_x + 1
			upper_right = upper_left + 1
			elements.extend(
				[
					(lower_left, lower_right, upper_right),
					(lower_left, upper_right, upper_left),
				]
			)

	return nodes, np.asarray(elements, dtype=np.int64)


def triangle_matrices(coordinates, conductivity, source):
	x1, y1 = coordinates[0]
	x2, y2 = coordinates[1]
	x3, y3 = coordinates[2]
	signed_twice_area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
	area = abs(signed_twice_area) / 2.0
	if area <= 0.0:
		raise ValueError("A malha contem um elemento triangular degenerado.")

	gradients = np.array(
		[
			[y2 - y3, x3 - x2],
			[y3 - y1, x1 - x3],
			[y1 - y2, x2 - x1],
		],
		dtype=np.float64,
	) / signed_twice_area
	stiffness = conductivity * area * (gradients @ gradients.T)
	load = np.full(3, source * area / 3.0, dtype=np.float64)
	return stiffness, load


def solve_heat_equation():
	nodes, elements = build_mesh(Lx, Ly, nx, ny)
	system = lil_matrix((len(nodes), len(nodes)), dtype=np.float64)
	load = np.zeros(len(nodes), dtype=np.float64)

	for element in elements:
		stiffness, element_load = triangle_matrices(nodes[element], k, q)
		system[np.ix_(element, element)] += stiffness
		load[element] += element_load

	boundary_nodes = np.flatnonzero(
		np.isclose(nodes[:, 0], 0.0)
		| np.isclose(nodes[:, 0], Lx)
		| np.isclose(nodes[:, 1], 0.0)
		| np.isclose(nodes[:, 1], Ly)
	)
	interior_nodes = np.setdiff1d(np.arange(len(nodes)), boundary_nodes)

	temperatures = np.zeros(len(nodes), dtype=np.float64)
	temperatures[interior_nodes] = spsolve(
		system.tocsr()[np.ix_(interior_nodes, interior_nodes)],
		load[interior_nodes],
	)
	return nodes, elements, temperatures, boundary_nodes


def save_temperature_plot(nodes, elements, temperatures):
	figure, axis = plt.subplots(figsize=(7, 6))
	triangulation = mtri.Triangulation(
		nodes[:, 0], nodes[:, 1], elements
	)
	contour = axis.tricontourf(triangulation, temperatures, levels=30, cmap="jet")
	figure.colorbar(contour, ax=axis, label="Temperatura")
	axis.set_title("Distribuicao de temperatura - FEM")
	axis.set_xlabel("x")
	axis.set_ylabel("y")
	axis.set_aspect("equal")
	figure.tight_layout()

	output_directory = Path(__file__).resolve().parent / "results/difusion"
	output_directory.mkdir(parents=True, exist_ok=True)
	timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
	output_path = output_directory / f"resultado_heat_fem_{timestamp}.png"
	figure.savefig(output_path, dpi=150, bbox_inches="tight")
	plt.close(figure)
	return output_path


if __name__ == "__main__":
	mesh_nodes, mesh_elements, temperature, boundary = solve_heat_equation()
	maximum_temperature_node = int(np.argmax(temperature))
	output_file = save_temperature_plot(
		mesh_nodes, mesh_elements, temperature
	)
	print(f"Nos: {len(mesh_nodes)} | Elementos: {len(mesh_elements)}")
	print(f"Nos de contorno: {len(boundary)}")
	print(
		"Temperatura maxima: "
		f"{temperature[maximum_temperature_node]:.6f} "
		f"em ({mesh_nodes[maximum_temperature_node, 0]:.6f}, "
		f"{mesh_nodes[maximum_temperature_node, 1]:.6f})"
	)
	print(f"Resultado salvo em: {output_file}")
