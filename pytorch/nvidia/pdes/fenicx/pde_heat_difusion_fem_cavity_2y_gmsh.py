"""Solve the double-Y heat PDE with Gmsh and modern FEniCSx."""

import argparse

try:
    from .common import case_dimensions, mesh_from_gmsh, solve
except ImportError:
    from common import case_dimensions, mesh_from_gmsh, solve


def run(case=1, mesh_ref=30):
    length, height, polygon = case_dimensions(case)
    mesh = mesh_from_gmsh(length, height, polygon, mesh_ref)
    maximum, output = solve(mesh, polygon, f"gmsh_case{case}")
    print(f"Gmsh case {case}: temperature_max={maximum:.15f}")
    print(f"output={output}")
    return maximum


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=int, choices=(1, 2), default=None)
    parser.add_argument("--mesh-ref", type=int, default=30)
    args = parser.parse_args()
    cases = (args.case,) if args.case else (1, 2)
    for case in cases:
        run(case, args.mesh_ref)


if __name__ == "__main__":
    main()
