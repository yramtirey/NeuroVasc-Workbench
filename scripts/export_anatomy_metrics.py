import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from neurovasc.io.volume import load_nifti
from neurovasc.segmentation.anatomy_metrics import (
    measure_anatomical_label,
)


def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Measure TopCoW anatomical vessel labels "
            "with NeuroVasc Workbench."
        )
    )

    parser.add_argument(
        "path",
        help="Path to TopCoW multiclass segmentation.",
    )

    parser.add_argument(
        "--output",
        default=(
            "outputs/measurements/"
            "case_001_anatomy_metrics.csv"
        ),
    )

    args = parser.parse_args()

    volume = load_nifti(
        args.path
    )

    labels = sorted(
        int(label)
        for label in np.unique(volume.data)
        if int(label) != 0
    )

    records = []

    print()
    print("NeuroVasc Workbench")
    print("=" * 75)
    print("Measuring anatomical vessels...")
    print()

    for label in labels:

        print(
            f"Processing label {label}..."
        )

        metrics = measure_anatomical_label(
            volume,
            label,
        )

        records.append(
            {
                "label": metrics.label,
                "abbreviation": metrics.abbreviation,
                "name": metrics.name,
                "voxel_count": metrics.voxel_count,
                "voxel_volume_mm3": metrics.voxel_volume_mm3,
                "mesh_volume_mm3": metrics.mesh_volume_mm3,
                "surface_area_mm2": metrics.surface_area_mm2,
                "centerline_length_mm": metrics.centerline_length_mm,
                "centerline_points": metrics.centerline_points,
                "connected_components": metrics.connected_components,
                "mean_diameter_mm": metrics.mean_diameter_mm,
                "median_diameter_mm": metrics.median_diameter_mm,
                "min_diameter_mm": metrics.min_diameter_mm,
                "max_diameter_mm": metrics.max_diameter_mm,
                "p10_diameter_mm": metrics.p10_diameter_mm,
                "p90_diameter_mm": metrics.p90_diameter_mm,
            }
        )

    table = pd.DataFrame(
        records
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

    display_columns = [
        "abbreviation",
        "centerline_length_mm",
        "mean_diameter_mm",
        "median_diameter_mm",
        "min_diameter_mm",
        "max_diameter_mm",
        "connected_components",
    ]

    print()
    print(
        table[
            display_columns
        ].to_string(
            index=False,
            float_format=lambda x: f"{x:.2f}",
        )
    )

    print()
    print(
        f"Saved → {output}"
    )
    print("=" * 75)


if __name__ == "__main__":
    main()