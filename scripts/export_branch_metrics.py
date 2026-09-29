import argparse
from pathlib import Path

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.graph.branch_graph import build_branch_graph
from neurovasc.graph.metrics import branch_metrics_table
from neurovasc.graph.vascular_graph import build_vascular_graph
from neurovasc.io.volume import load_nifti


def main() -> None:

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "path",
        help="Path to vascular NIfTI mask.",
    )

    parser.add_argument(
        "--output",
        default="outputs/measurements/branch_metrics.csv",
    )

    args = parser.parse_args()

    volume = load_nifti(
        args.path
    )

    centerline = extract_centerline(
        volume
    )

    vascular_graph = build_vascular_graph(
        centerline
    )

    branch_graph = build_branch_graph(
        vascular_graph
    )

    table = branch_metrics_table(
        branch_graph
    )

    output = Path(
        args.output
    )

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    table.to_csv(
        output,
        index=False,
    )

    print()
    print("NeuroVasc Workbench")
    print("=" * 60)

    print(table.to_string(index=False))

    print()
    print(
        f"Saved branch measurements → {output}"
    )

    print("=" * 60)


if __name__ == "__main__":
    main()