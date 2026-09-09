import csv
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from scipy.spatial import Delaunay
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

def build_cavity_polygon(h0, l0, s1, h1, l1, h2, l2, alpha, beta):
	"""Return the 18-point double-Y cavity from pde_cav_2y.m."""
	alpha = np.deg2rad(alpha)
	beta = np.deg2rad(beta)
	cos_alpha = np.cos(alpha)
	sin_alpha = np.sin(alpha)
	cos_beta = np.cos(beta)
	sin_beta = np.sin(beta)

	x1 = -l0 / 2.0
	y1 = 0.0
	x2 = x1
	y2 = s1 - h1 / (2.0 * cos_alpha)
	x3 = -(cos_alpha * l1 + l0 / 2.0)
	y3 = l1 * sin_alpha + y2
	x4 = x3
	y4 = y3 + h1 / cos_alpha
	x5 = x1
	y5 = y2 + h1 / cos_alpha
	x6 = x1
	y6 = h0 - h2 / cos_beta
	x7 = -(cos_beta * l2 + l0 / 2.0)
	y7 = h0 + l2 * sin_beta - h2 / cos_beta
	x8 = x7
	y8 = h0 + l2 * sin_beta
	x9 = x1
	y9 = h0

	left_side = np.array(
		[(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5),
		 (x6, y6), (x7, y7), (x8, y8), (x9, y9)],
		dtype=np.float64,
	)
	right_side = left_side[::-1].copy()
	right_side[:, 0] *= -1.0
	return np.vstack((left_side, right_side))


def _points_on_segments(polygon, spacing):
	points = []
	for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
		segment_length = np.linalg.norm(end - start)
		count = max(1, int(np.ceil(segment_length / spacing)))
		points.extend(start + fraction * (end - start) for fraction in np.linspace(0.0, 1.0, count, endpoint=False))
	return np.asarray(points, dtype=np.float64)


def _orientation(first, second, third):
	first_to_second = second - first
	first_to_third = third - first
	return (
		first_to_second[0] * first_to_third[1]
		- first_to_second[1] * first_to_third[0]
	)


def _segments_cross(first, second, third, fourth, tolerance=1e-12):
	values = (
		_orientation(first, second, third),
		_orientation(first, second, fourth),
		_orientation(third, fourth, first),
		_orientation(third, fourth, second),
	)
	return (values[0] * values[1] < -tolerance) and (values[2] * values[3] < -tolerance)


def build_mesh(length_x, length_y, mesh_ref, cavity_polygon):
	spacing = min(length_x, length_y) / mesh_ref
	x_values = np.linspace(-length_x / 2.0, length_x / 2.0, mesh_ref + 1)
	y_values = np.linspace(0.0, length_y, mesh_ref + 1)
	if _is_axis_aligned(cavity_polygon):
		x_values = np.unique(np.r_[x_values, cavity_polygon[:, 0]])
		y_values = np.unique(np.r_[y_values, cavity_polygon[:, 1]])
		nodes = np.asarray(
		[(x, y) for y in y_values for x in x_values], dtype=np.float64
	)
		row_width = len(x_values)
		elements = []
		for row in range(len(y_values) - 1):
			for column in range(len(x_values) - 1):
				cell_center = np.array(
					[
						(x_values[column] + x_values[column + 1]) / 2.0,
						(y_values[row] + y_values[row + 1]) / 2.0,
					]
				)
				if _point_in_polygon(cell_center, cavity_polygon):
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
		return nodes, np.asarray(elements, dtype=np.int64)

	grid = np.asarray([(x, y) for y in y_values for x in x_values], dtype=np.float64)
	boundary_points = _points_on_segments(cavity_polygon, spacing)
	nodes = np.unique(np.vstack((grid, boundary_points)), axis=0)

	triangulation = Delaunay(nodes)
	cavity_edges = list(zip(cavity_polygon, np.roll(cavity_polygon, -1, axis=0)))
	elements = []
	for element in triangulation.simplices:
		coordinates = nodes[element]
		centroid = coordinates.mean(axis=0)
		if _point_in_polygon(centroid, cavity_polygon):
			continue
		edge_midpoints = [
			(coordinates[first] + coordinates[second]) / 2.0
			for first, second in ((0, 1), (1, 2), (2, 0))
		]
		if any(_point_in_polygon(midpoint, cavity_polygon) for midpoint in edge_midpoints):
			continue
		crosses_cavity = any(
			_segments_cross(coordinates[first], coordinates[second], start, end)
			for first, second in ((0, 1), (1, 2), (2, 0))
			for start, end in cavity_edges
		)
		if not crosses_cavity:
			elements.append(element)

	return nodes, np.asarray(elements, dtype=np.int64)


def _is_axis_aligned(polygon, tolerance=1e-12):
	for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
		if (
			abs(start[0] - end[0]) > tolerance
			and abs(start[1] - end[1]) > tolerance
		):
			return False
	return True


def _point_in_polygon(point, polygon):
	inside = False
	for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
		if (start[1] > point[1]) != (end[1] > point[1]):
			x_intersection = (end[0] - start[0]) * (point[1] - start[1]) / (end[1] - start[1]) + start[0]
			if point[0] < x_intersection:
				inside = not inside
	return inside


def cavity_boundary_nodes(nodes, cavity_polygon, tolerance=1e-10):
	boundary_nodes = []
	for node_index, point in enumerate(nodes):
		for start, end in zip(cavity_polygon, np.roll(cavity_polygon, -1, axis=0)):
			if abs(_orientation(start, end, point)) <= tolerance:
				if np.dot(point - start, point - end) <= tolerance:
					boundary_nodes.append(node_index)
					break
	return np.asarray(boundary_nodes, dtype=np.int64)


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


def solve_heat_equation(h0, l0, s1, h1, l1, h2, l2, alpha, beta, mesh_ref):
	cavity_polygon = build_cavity_polygon(
		h0, l0, s1, h1, l1, h2, l2, alpha, beta
	)
	if (
		cavity_polygon[:, 0].min() < DomainXMin
		or cavity_polygon[:, 0].max() > DomainXMax
		or cavity_polygon[:, 1].min() < DomainYMin
		or cavity_polygon[:, 1].max() > DomainYMax
	):
		raise ValueError("A cavidade ultrapassa os limites do solido.")
	nodes, elements = build_mesh(Lx, Ly, mesh_ref, cavity_polygon)
	system = lil_matrix((len(nodes), len(nodes)), dtype=np.float64)
	load = np.zeros(len(nodes), dtype=np.float64)

	for element in elements:
		stiffness, element_load = triangle_matrices(nodes[element], k, q)
		system[np.ix_(element, element)] += stiffness
		load[element] += element_load

	cavity_nodes = cavity_boundary_nodes(nodes, cavity_polygon)
	used_nodes = np.unique(elements.ravel())
	unknown_nodes = np.setdiff1d(used_nodes, cavity_nodes)

	temperatures = np.zeros(len(nodes), dtype=np.float64)
	temperatures[unknown_nodes] = spsolve(
		system.tocsr()[np.ix_(unknown_nodes, unknown_nodes)],
		load[unknown_nodes],
	)
	return nodes, elements, temperatures, cavity_nodes, cavity_polygon


def save_temperature_plot(
	nodes, elements, temperatures, cavity_polygon, name, dir, temperature_max=0.0
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
	closed_cavity_polygon = np.vstack((cavity_polygon, cavity_polygon[0]))
	axis.plot(
		closed_cavity_polygon[:, 0], closed_cavity_polygon[:, 1],
		color="cyan", linewidth=1.0,
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


def solve_dofs_2y(mesh, dofs,fis):
	mesh_ref = mesh
	H0L0 = dofs[0]
	alpha = dofs[4]
	beta = dofs[5]
	s1_h0 = dofs[3]
	fi1 = fis[0]
	fi2 = fis[1]
	fi3 = fis[2]
	h1_l1 = dofs[2]
	h2_l2 = dofs[1]
	psi = 1.0

	timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
	dirname = f"cav2y_{timestamp}"

	results_directory = Path(__file__).resolve().parent / "results/cavity" / dirname 
	results_directory.mkdir(parents=True, exist_ok=True)
	csv_path = results_directory / "temperatura_maxima_h0l0.csv"
	results = []

	# for H0L0 in H0L0s:
	H0 = np.sqrt(fi3 * H0L0)
	L0 = np.sqrt(fi3 / H0L0)
	H1 = np.sqrt(fi1 * h1_l1)
	L1 = np.sqrt(fi1 / h1_l1)
	H2 = np.sqrt(fi2 * h2_l2)
	L2 = np.sqrt(fi2 / h2_l2)
	S1 = s1_h0 * H0
	mesh_nodes, mesh_elements, temperature, cavity, cavity_polygon = solve_heat_equation(
		H0, L0, S1, H1, L1, H2, L2, alpha, beta, mesh_ref
	)
	maximum_temperature_node = int(np.argmax(temperature))
	temperature_max = float(temperature[maximum_temperature_node])
	results.append((H0L0, temperature_max))

	output_file = save_temperature_plot(
		mesh_nodes, mesh_elements, temperature, cavity_polygon,
		f"dofs_{H0L0}", dirname
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
	mesh = 250
	# Case 1 - Forma em T para verificacao de malha
	H0L0 = 20
	alpha = 0.0
	beta = 0.0
	s1_h0 = 0.45
	fi1 = 0.02
	fi2 = 0.02
	fi3 = 0.02
	h1_l1 = 14.0
	h2_l2 = 0.14
 
	# Case 2 - Forma otima para 1dof a=b
	H0L0 = 10
	alpha = 13.687500
	beta = alpha
	s1_h0 = 0.5
	fi1 = 0.015
	fi2 = 0.015
	fi3 = 0.04
	h1_l1 = 0.07
	h2_l2 = h1_l1
 
	dofs = [H0L0, h2_l2, h1_l1, s1_h0, alpha, beta]
	fis = [fi1, fi2, fi3]
	solve_dofs_2y(mesh, dofs, fis)