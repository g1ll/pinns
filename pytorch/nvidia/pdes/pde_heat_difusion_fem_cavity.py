import csv
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


# Problema estacionario: -k * div(grad(T)) = q no dominio com cavidade.
# As laterais externas sao adiabaticas (Neumann natural n.grad(T) = 0).
# A cavidade e mantida em T = 0 (Dirichlet), sendo a unica superficie de troca.
Lx = 1.0
Ly = 1.0
k = 1
q = 1
As = Lx * Ly
Ac = 0.1 * As

DomainXMin = -Lx / 2.0
DomainXMax = Lx / 2.0
DomainYMin = 0.0
DomainYMax = Ly
CavityXCenter = 0.0
CavityYStart = 0.0

def build_mesh(length_x, length_y, elements_x, elements_y,h0,l0):
	tx_min = CavityXCenter - l0 / 2.0
	tx_max = CavityXCenter + l0 / 2.0
	ty_min = CavityYStart
	ty_max = CavityYStart + h0
	if tx_min < DomainXMin or tx_max > DomainXMax:
		raise ValueError("A cavidade ultrapassa os limites em x do solido.")
	if ty_min < DomainYMin or ty_max > DomainYMax:
		raise ValueError(
			"A cavidade ultrapassa os limites em y do solido. "
			"Reduza H0L0 ou aumente Ly."
		)

	x_values = np.unique(
		 np.r_[np.linspace(DomainXMin, DomainXMax, elements_x + 1), tx_min, tx_max]
	)
	y_values = np.unique(
		 np.r_[np.linspace(DomainYMin, DomainYMax, elements_y + 1), ty_min, ty_max]
	)
	nodes = np.array(
		[(x, y) for y in y_values for x in x_values], dtype=np.float64
	)

	elements = []
	row_width = len(x_values)
	for row in range(len(y_values) - 1):
		for column in range(len(x_values) - 1):
			x_center = (x_values[column] + x_values[column + 1]) / 2.0
			y_center = (y_values[row] + y_values[row + 1]) / 2.0
			inside_cavity = (
				tx_min < x_center < tx_max and ty_min < y_center < ty_max
			)
			if inside_cavity:
				continue

			lower_left = row * row_width + column
			lower_right = lower_left + 1
			upper_left = lower_left + row_width
			upper_right = upper_left + 1
			elements.extend(
				[
					(lower_left, lower_right, upper_right),
					(lower_left, upper_right, upper_left),
				]
			)

	return nodes, np.asarray(elements, dtype=np.int64), (tx_min, tx_max, ty_min, ty_max)


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


def solve_heat_equation(h0,l0, mesh_ref):
	nx = mesh_ref
	ny = mesh_ref
	nodes, elements, cavity_bounds = build_mesh(Lx, Ly, nx, ny,h0,l0)
	system = lil_matrix((len(nodes), len(nodes)), dtype=np.float64)
	load = np.zeros(len(nodes), dtype=np.float64)

	for element in elements:
		stiffness, element_load = triangle_matrices(nodes[element], k, q)
		system[np.ix_(element, element)] += stiffness
		load[element] += element_load

	tx_min, tx_max, ty_min, ty_max = cavity_bounds
	cavity_nodes = np.flatnonzero(
		(
			np.isclose(nodes[:, 0], tx_min)
			| np.isclose(nodes[:, 0], tx_max)
		)
		& (nodes[:, 1] >= ty_min - 1e-12)
		& (nodes[:, 1] <= ty_max + 1e-12)
		| (
			(
				np.isclose(nodes[:, 1], ty_min)
				| np.isclose(nodes[:, 1], ty_max)
			)
			& (nodes[:, 0] >= tx_min - 1e-12)
			& (nodes[:, 0] <= tx_max + 1e-12)
		)
	)
	used_nodes = np.unique(elements.ravel())
	unknown_nodes = np.setdiff1d(used_nodes, cavity_nodes)

	temperatures = np.zeros(len(nodes), dtype=np.float64)
	temperatures[unknown_nodes] = spsolve(
		system.tocsr()[np.ix_(unknown_nodes, unknown_nodes)],
		load[unknown_nodes],
	)
	return nodes, elements, temperatures, cavity_nodes, cavity_bounds


def save_temperature_plot(
	nodes, elements, temperatures, cavity_bounds, name, dir, temperature_max=0.0
):
	if temperature_max < 0.0:
		raise ValueError("temperature_max deve ser zero ou um valor positivo.")

	plot_temperature_max = (
		float(np.max(temperatures)) if temperature_max == 0.0 else temperature_max
	)
	if plot_temperature_max <= 0.0:
		raise ValueError("A temperatura maxima do dominio deve ser positiva.")
	levels = np.linspace(0.0, plot_temperature_max, 31)
	figure, axis = plt.subplots(figsize=(7, 7))
	triangulation = mtri.Triangulation(
		nodes[:, 0], nodes[:, 1], elements
	)
	contour = axis.tricontourf(
		triangulation,
		temperatures,
		levels=levels,
		cmap="jet",
		extend="max",
	)
	tx_min, tx_max, ty_min, ty_max = cavity_bounds
	cavity_polygon = np.array(
		[
			[tx_min, ty_min],
			[tx_max, ty_min],
			[tx_max, ty_max],
			[tx_min, ty_max],
			[tx_min, ty_min],
		]
	)
	axis.plot(
		cavity_polygon[:, 0], cavity_polygon[:, 1], color="cyan", linewidth=1.0,
		label="Cavidade: T = 0",
	)
	axis.legend(loc="upper right")
	figure.colorbar(contour, ax=axis, label="Temperatura")
	axis.set_title("Distribuicao de temperatura - FEM")
	axis.set_xlabel("x")
	axis.set_ylabel("y")
	axis.set_aspect("equal", adjustable="box")
	axis.set_box_aspect(1)
	figure.tight_layout()


	timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
	output_directory = Path(__file__).resolve().parent / "results/cavity" / dir 
	output_directory.mkdir(parents=True, exist_ok=True)
	output_path = output_directory / f"resultado_cavity_h0l0_{name}_{timestamp}.png"
	figure.savefig(output_path, dpi=150, bbox_inches="tight")
	plt.close(figure)
	return output_path


def solve_dofs_h0l0():
	mesh_ref = 50
	H0L0s = [.125, .25, .5, 1, 2, 4, 6, 8, 9.25]

	timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
	dirname = f"h0l0_{timestamp}"

	results_directory = Path(__file__).resolve().parent / "results/cavity" / dirname 
	results_directory.mkdir(parents=True, exist_ok=True)
	csv_path = results_directory / "temperatura_maxima_h0l0.csv"
	results = []

	for H0L0 in H0L0s:
		H0 = np.sqrt(Ac * H0L0)
		L0 = np.sqrt(Ac / H0L0)
		mesh_nodes, mesh_elements, temperature, cavity, cavity_bounds = solve_heat_equation(H0,L0, mesh_ref)
		maximum_temperature_node = int(np.argmax(temperature))
		temperature_max = float(temperature[maximum_temperature_node])
		results.append((H0L0, temperature_max))

		output_file = save_temperature_plot(
			mesh_nodes, mesh_elements, temperature, cavity_bounds, f"dofs_{H0L0}", dirname, .5
		)
		print(f"H0L0: {H0L0} ")
		print(f"Nos: {len(mesh_nodes)} | Elementos: {len(mesh_elements)}")
		print(f"Area do solido: {As:.6f}")
		print(f"Area da cavidade: {Ac:.6f} ({Ac / As:.1%} de As)")
		print(f"Dimensoes da cavidade (L0 x H0): {L0:.6f} x {H0:.6f}")
		print(f"Nos da cavidade (Dirichlet T=0): {len(cavity)}")
		print(
			"Temperatura maxima: "
			f"{temperature_max:.6f} "
			f"em ({mesh_nodes[maximum_temperature_node, 0]:.6f}, "
			f"{mesh_nodes[maximum_temperature_node, 1]:.6f})"
		)
		print(f"Resultado salvo em: {output_file}")

	with csv_path.open("w", newline="", encoding="utf-8") as csv_file:
		writer = csv.writer(csv_file)
		writer.writerow(["h0l0", "temp_max"])
		writer.writerows(results)
	print(f"Registro CSV salvo em: {csv_path}")
	

if __name__ == "__main__":
	solve_dofs_h0l0()