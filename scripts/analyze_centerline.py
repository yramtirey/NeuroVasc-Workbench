import argparse

from neurovasc.geometry.centerline import (
    extract_centerline,
)
from neurovasc.graph.vascular_graph import (
    build_vascular_graph,
)
from neurovasc.io.volume import load_nifti


def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Extract and analyze a vascular "
            "centerline with NeuroVasc Workbench."
        )
    )

    parser.add_argument(
        "path",
        help="Path to vascular NIfTI mask.",
    )

    args = parser.parse_args()

    print()
    print("NeuroVasc Workbench")
    print("=" * 55)

    volume = load_nifti(
        args.path
    )

    print(
        f"Volume: {volume.path.name}"
    )

    print(
        "Extracting 3D centerline..."
    )

    centerline = extract_centerline(
        volume
    )

    print(
        "Building vascular graph..."
    )

    vascular_graph = build_vascular_graph(
        centerline
    )

    print()
    print("Centerline Analysis")
    print("-" * 55)

    print(
        f"Centerline voxels:   "
        f"{centerline.n_points:,}"
    )

    print(
        f"Graph nodes:         "
        f"{vascular_graph.n_nodes:,}"
    )

    print(
        f"Graph edges:         "
        f"{vascular_graph.n_edges:,}"
    )

    print(
        f"Connected components:"
        f" {vascular_graph.n_components:,}"
    )

    print(
        f"Endpoints:           "
        f"{len(vascular_graph.endpoints):,}"
    )

    print(
        f"Junction voxels:     "
        f"{len(vascular_graph.branchpoints):,}"
    )

    print(
        f"Centerline length:   "
        f"{vascular_graph.total_length_mm:,.2f} mm"
    )

    print("=" * 55)


if __name__ == "__main__":
    main()