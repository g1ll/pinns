"""FEniCS heat solver using a CAD-conforming Gmsh mesh."""

import argparse
import numpy as np
import os
import sys

from fenics_common import cavity_polygon, default_case, solve_fenics


def gmsh_mesh(length, height, polygon, mesh_ref):
	venv_site = os.environ.get("FENICS_VENV_SITE")
	if venv_site and venv_site not in sys.path:
		sys.path.append(venv_site)
	try:
		import gmsh
	except ImportError as error:
		raise RuntimeError("Instale o modulo Python gmsh para executar este script.") from error
	gmsh.initialize()
	try:
		gmsh.model.add("cavity_2y")
		geo = gmsh.model.geo
		lc = min(length, height) / mesh_ref
		outer = [geo.addPoint(-length / 2, 0, 0, lc),
			geo.addPoint(length / 2, 0, 0, lc), geo.addPoint(length / 2, height, 0, lc),
			geo.addPoint(-length / 2, height, 0, lc)]
		cavity_points = [geo.addPoint(float(x), float(y), 0, lc) for x, y in polygon]
		boundary_points = [outer[0], *cavity_points, outer[1], outer[2], outer[3]]
		boundary_lines = [geo.addLine(boundary_points[index], boundary_points[index + 1])
			for index in range(len(boundary_points) - 1)]
		boundary_lines.append(geo.addLine(boundary_points[-1], boundary_points[0]))
		boundary_loop = geo.addCurveLoop(boundary_lines)
		geo.addPlaneSurface([boundary_loop])
		geo.synchronize()
		gmsh.option.setNumber("Mesh.Algorithm", 6)
		gmsh.model.mesh.generate(2)
		node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
		nodes = np.asarray(coordinates, dtype=np.float64).reshape(-1, 3)[:, :2]
		order = np.argsort(node_tags)
		nodes = nodes[order]
		types, _, node_blocks = gmsh.model.mesh.getElements(2)
		triangles = []
		for element_type, block in zip(types, node_blocks):
			if element_type != 2:
				continue
			triangles.extend(np.asarray(block, dtype=np.int64).reshape(-1, 3))
		tag_to_index = {int(tag): index for index, tag in enumerate(np.asarray(node_tags)[order])}
		elements = np.asarray([[tag_to_index[int(tag)] for tag in triangle] for triangle in triangles])
		return nodes, elements
	finally:
		gmsh.finalize()


def run(mesh_ref=100, case=1):
	length, height, h0, l0, s1, h1, l1, h2, l2, alpha, beta = default_case(case)
	polygon = cavity_polygon(h0, l0, s1, h1, l1, h2, l2, alpha, beta)
	nodes, elements = gmsh_mesh(length, height, polygon, mesh_ref)
	_, _, maximum, output = solve_fenics(nodes, elements, polygon, f"gmsh_case{case}")
	print(f"Malha Gmsh: {len(nodes)} nos, {len(elements)} elementos")
	print(f"Temperatura maxima: {maximum:.12f}")
	print(f"Resultado salvo em: {output}")
	return maximum


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--mesh-ref", type=int, default=20)
	args = parser.parse_args()
	for case in (1, 2):
		run(args.mesh_ref, case)