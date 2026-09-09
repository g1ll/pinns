"""FEniCS heat solver using a conforming NumPy mesh."""

import argparse

from fenics_common import cavity_polygon, default_case, numpy_mesh, solve_fenics


def run(mesh_ref=100, case=1):
	length, height, h0, l0, s1, h1, l1, h2, l2, alpha, beta = default_case(case)
	polygon = cavity_polygon(h0, l0, s1, h1, l1, h2, l2, alpha, beta)
	nodes, elements = numpy_mesh(length, height, mesh_ref, polygon)
	_, _, maximum, output = solve_fenics(nodes, elements, polygon, f"numpy_case{case}")
	print(f"Malha NumPy: {len(nodes)} nos, {len(elements)} elementos")
	print(f"Temperatura maxima: {maximum:.12f}")
	print(f"Resultado salvo em: {output}")
	return maximum


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--mesh-ref", type=int, default=20)
	args = parser.parse_args()
	for case in (1, 2):
		run(args.mesh_ref, case)