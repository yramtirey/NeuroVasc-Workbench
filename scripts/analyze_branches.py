import argparse

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.graph.branch_graph import build_branch_graph
from neurovasc.graph.vascular_graph import build_vascular_graph
from neurovasc.io.volume import load_nifti


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze branch-level vessel structure."
    )
    parser.add_argument(
        "path",
        help="Path to vascular NIfTI mask.",
    )
    args = parser.parse_args()

    print()
    print("NeuroVasc Workbench")
    print("=" * 60)

    volume = load_nifti(args.path)
    centerline = extract_centerline(volume)
    vascular_graph = build_vascular_graph(centerline)
    branch_graph = build_branch_graph(vascular_graph)

    print("Branch-Level Analysis")
    print("-" * 60)
    print(f"Raw centerline voxels:     {centerline.n_points}")
    print(f"Raw graph nodes:           {vascular_graph.n_nodes}")
    print(f"Raw graph edges:           {vascular_graph.n_edges}")
    print(f"Logical nodes:             {branch_graph.n_logical_nodes}")
    print(f"Logical edges:             {branch_graph.n_logical_edges}")
    print(f"Logical endpoints:         {len(branch_graph.endpoint_nodes)}")
    print(f"Logical junctions:         {len(branch_graph.junction_nodes)}")
    print(f"Branches:                  {branch_graph.n_branches}")
    print("-" * 60)

    for branch in branch_graph.branches:
        print(
            f"{branch.branch_id}: "
            f"{branch.start_node} → {branch.end_node} | "
            f"Length = {branch.length_mm:.2f} mm | "
            f"Chord = {branch.chord_length_mm:.2f} mm | "
            f"Tortuosity = {branch.tortuosity:.3f}"
        )

    print("=" * 60)


if __name__ == "__main__":
    main()