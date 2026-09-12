"""Solve the double-Y heat PDE with a NumPy-generated mesh and FEniCSx."""

import argparse

try:
    from .common import case_dimensions, mesh_from_numpy, numpy_mesh, solve
except ImportError:
    from common import case_dimensions, mesh_from_numpy, numpy_mesh, solve


def run(case=1, mesh_ref=30):
    length, height, polygon = case_dimensions(case)
    nodes, cells = numpy_mesh(length, height, mesh_ref, polygon)
    mesh = mesh_from_numpy(nodes, cells)
    maximum, output = solve(mesh, polygon, f"numpy_case{case}")
    print(f"NumPy case {case}: {len(nodes)} nodes, {len(cells)} elements")
    print(f"temperature_max={maximum:.15f}")
    print(f"output={output}")
    return maximum


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=int, choices=(1, 2), default=None)
    parser.add_argument("--mesh-ref", type=int, default=30)
    args = parser.parse_args()
    cases = (args.case,) if args.case else (1,)
    for case in cases:
        run(case, args.mesh_ref)


if __name__ == "__main__":
    main()
