"""Shared geometry, mesh import and FEniCS solve for the MATLAB double-Y case."""

import csv
import os
import sys
from datetime import datetime
from pathlib import Path


def _use_system_fenics_runtime():
	"""Avoid loading Ubuntu's binary DOLFIN through the Python venv."""
	if os.environ.get("FENICS_SYSTEM_RUNTIME") == "1":
		return
	if sys.executable != "/usr/bin/python3" and Path("/usr/bin/python3").exists():
		environment = os.environ.copy()
		environment.pop("PYTHONPATH", None)
		environment.pop("LD_LIBRARY_PATH", None)
		environment.pop("VIRTUAL_ENV", None)
		environment["FENICS_SYSTEM_RUNTIME"] = "1"
		environment["FENICS_VENV_SITE"] = next(
			(path for path in sys.path if path.endswith("/site-packages")), ""
		)
		if sys.argv[0] == "-c":
			return
		os.execve("/usr/bin/python3", ["/usr/bin/python3", *sys.argv], environment)


_use_system_fenics_runtime()

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
	left = np.array([(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5),
		(x6, y6), (x7, y7), (x8, y8), (x9, y9)])
	right = left[::-1].copy()
	right[:, 0] *= -1.0
	return np.vstack((left, right))


def dimensions(dofs, fis):
	"""Return L, H and cavity dimensions using the equations in dof2y.m."""
	hd_l, h0d_l0, h2d_l2, h1d_l1, s1d_h0, alpha, beta = dofs
	fi1, fi2, fi3 = fis
	length = np.sqrt(1.0 / hd_l)
	height = hd_l * length
	l0 = np.sqrt(fi3 / h0d_l0)
	h0 = h0d_l0 * l0
	l1 = np.sqrt(fi1 / h1d_l1)
	h1 = h1d_l1 * l1
	l2 = np.sqrt(fi2 / h2d_l2)
	h2 = h2d_l2 * l2
	return (length, height, h0, l0, s1d_h0 * h0, h1, l1, h2, l2,
		alpha, beta)


def _point_in_polygon(point, polygon):
	inside = False
	for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
		if (start[1] > point[1]) != (end[1] > point[1]):
			x = ((end[0] - start[0]) * (point[1] - start[1]) /
				(end[1] - start[1]) + start[0])
			if point[0] < x:
				inside = not inside
	return inside


def _orientation(first, second, third):
	vector_a = second - first
	vector_b = third - first
	return vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0]


def _segments_cross(first, second, third, fourth):
	values = (_orientation(first, second, third), _orientation(first, second, fourth),
		_orientation(third, fourth, first), _orientation(third, fourth, second))
	return values[0] * values[1] < -1e-12 and values[2] * values[3] < -1e-12


def _points_on_segments(polygon, spacing):
	points = []
	for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
		count = max(1, int(np.ceil(np.linalg.norm(end - start) / spacing)))
		points.extend(start + fraction * (end - start)
			for fraction in np.linspace(0.0, 1.0, count, endpoint=False))
	return np.asarray(points, dtype=np.float64)


def numpy_mesh(length, height, mesh_ref, polygon):
	"""Create a conforming structured triangular mesh for a=b=0."""
	x = np.unique(np.r_[np.linspace(-length / 2, length / 2, mesh_ref + 1), polygon[:, 0]])
	y = np.unique(np.r_[np.linspace(0.0, height, mesh_ref + 1), polygon[:, 1]])
	axis_aligned = all(np.isclose(start[0], end[0]) or np.isclose(start[1], end[1])
		for start, end in zip(polygon, np.roll(polygon, -1, axis=0)))
	if not axis_aligned:
		from scipy.spatial import Delaunay
		grid = np.asarray([(xi, yi) for yi in y for xi in x])
		nodes = np.unique(np.vstack((grid, _points_on_segments(polygon, min(length, height) / mesh_ref))), axis=0)
		triangulation = Delaunay(nodes)
		cavity_edges = list(zip(polygon, np.roll(polygon, -1, axis=0)))
		elements = []
		for element in triangulation.simplices:
			coordinates = nodes[element]
			if _point_in_polygon(coordinates.mean(axis=0), polygon):
				continue
			if any(_point_in_polygon((coordinates[first] + coordinates[second]) / 2.0, polygon)
				for first, second in ((0, 1), (1, 2), (2, 0))):
				continue
			if any(_segments_cross(coordinates[first], coordinates[second], start, end)
				for first, second in ((0, 1), (1, 2), (2, 0))
				for start, end in cavity_edges):
				continue
			elements.append(element)
		return nodes, np.asarray(elements, dtype=np.int64)
	nodes = np.asarray([(xi, yi) for yi in y for xi in x], dtype=np.float64)
	width = len(x)
	elements = []
	for row in range(len(y) - 1):
		for column in range(width - 1):
			center = np.array([(x[column] + x[column + 1]) / 2,
				(y[row] + y[row + 1]) / 2])
			if _point_in_polygon(center, polygon):
				continue
			ll = row * width + column
			lr, ul = ll + 1, ll + width
			ur = ul + 1
			elements.extend(((ll, lr, ur), (ll, ur, ul)))
	return nodes, np.asarray(elements, dtype=np.int64)


def dolfin_mesh(nodes, elements):
	try:
		from dolfin import Mesh, MeshEditor
	except ImportError as error:
		raise RuntimeError("Instale FEniCS (dolfin) para executar estes scripts.") from error
	mesh = Mesh()
	editor = MeshEditor()
	editor.open(mesh, "triangle", 2, 2)
	editor.init_vertices(len(nodes))
	editor.init_cells(len(elements))
	for index, point in enumerate(nodes):
		editor.add_vertex(index, [float(point[0]), float(point[1])])
	for index, element in enumerate(elements):
		editor.add_cell(index, [int(value) for value in element])
	editor.close()
	return mesh


def save_temperature_plot(nodes, elements, temperatures, polygon, output_directory, title):
	venv_site = os.environ.get("FENICS_VENV_SITE")
	if venv_site and venv_site not in sys.path:
		sys.path.append(venv_site)
	import matplotlib.pyplot as plt
	import matplotlib.tri as mtri

	plot_temperature_max = float(np.max(temperatures))
	if plot_temperature_max <= 0.0:
		raise ValueError("A temperatura maxima do dominio deve ser positiva.")
	levels = np.linspace(0.0, plot_temperature_max, 31)
	figure, axis = plt.subplots(figsize=(7, 7))
	triangulation = mtri.Triangulation(nodes[:, 0], nodes[:, 1], elements)
	contour = axis.tricontourf(
		triangulation,
		temperatures,
		levels=levels,
		cmap="jet",
		extend="max",
	)
	closed_polygon = np.vstack((polygon, polygon[0]))
	axis.plot(
		closed_polygon[:, 0], closed_polygon[:, 1],
		color="cyan", linewidth=1.0, label="Cavidade: T = 0",
	)
	axis.legend(loc="upper right")
	figure.colorbar(contour, ax=axis, label="Temperatura")
	axis.set_title(title)
	axis.set_xlabel("x")
	axis.set_ylabel("y")
	axis.set_aspect("equal", adjustable="box")
	axis.set_box_aspect(1)
	figure.tight_layout()
	output_path = output_directory / "resultado_temperatura.png"
	figure.savefig(output_path, dpi=150, bbox_inches="tight")
	plt.close(figure)
	return output_path


def solve_fenics(nodes, elements, polygon, mesh_name):
	try:
		from dolfin import (Constant, DirichletBC, Function, FunctionSpace,
			SubDomain, TestFunction, TrialFunction, XDMFFile, dot, dx, grad, solve)
	except ImportError as error:
		raise RuntimeError("Instale FEniCS (dolfin) para executar estes scripts.") from error

	used_nodes = np.unique(elements.reshape(-1))
	old_to_new = np.full(len(nodes), -1, dtype=np.int64)
	old_to_new[used_nodes] = np.arange(len(used_nodes))
	nodes = nodes[used_nodes]
	elements = old_to_new[elements]
	mesh = dolfin_mesh(nodes, elements)
	space = FunctionSpace(mesh, "P", 1)

	class CavityBoundary(SubDomain):
		def inside(self, point, on_boundary):
			if not on_boundary:
				return False
			for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
				segment = end - start
				length_squared = float(np.dot(segment, segment))
				parameter = np.clip(np.dot(point - start, segment) / length_squared, 0.0, 1.0)
				if np.linalg.norm(point - (start + parameter * segment)) < 1e-9:
					return True
			return False

	bc = DirichletBC(space, Constant(0.0), CavityBoundary())
	trial, test = TrialFunction(space), TestFunction(space)
	form = dot(grad(trial), grad(test)) * dx
	rhs = Constant(1.0) * test * dx
	temperature = Function(space)
	solve(form == rhs, temperature, bc)
	values = temperature.vector().get_local()
	output = RESULTS / f"{mesh_name}_{datetime.now():%y%m%d_%H%M%S}"
	output.mkdir(parents=True, exist_ok=True)
	with XDMFFile(mesh.mpi_comm(), str(output / "temperatura.xdmf")) as xdmf:
		xdmf.write(temperature)
	with (output / "temperatura_maxima.csv").open("w", newline="", encoding="utf-8") as file:
		writer = csv.writer(file)
		writer.writerow(["nodes", "elements", "temperature_max"])
		writer.writerow([len(nodes), len(elements), float(values.max())])
	with (output / "temperatura_maxima.log").open("w", encoding="utf-8") as file:
		file.write(f"nodes={len(nodes)}\n")
		file.write(f"elements={len(elements)}\n")
		file.write(f"temperature_max={float(values.max()):.15f}\n")
	save_temperature_plot(
		nodes,
		elements,
		values,
		polygon,
		output,
		f"Distribuicao de temperatura - FEniCS {mesh_name}",
	)
	return mesh, temperature, float(values.max()), output


def default_case(case=1):
	if case == 1:
		dofs = [1.0, 20.0, 0.14, 14.0, 0.45, 0.0, 0.0]
		fis = [0.02, 0.02, 0.02]
	elif case == 2:
		dofs = [1.0, 10.0, 0.07, 0.07, 0.5, 13.6875, 13.6875]
		fis = [0.015, 0.015, 0.04]
	else:
		raise ValueError("case deve ser 1 (forma em T) ou 2 (forma otima).")
	return dimensions(dofs, fis)