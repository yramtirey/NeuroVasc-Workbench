import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.diameter_profile import build_diameter_profile
from neurovasc.geometry.main_path import extract_main_path
from neurovasc.graph.vascular_graph import build_vascular_graph
from neurovasc.io.volume import Volume, load_nifti
from neurovasc.segmentation.labels import (
    TOPCOW_FULL_NAMES,
    TOPCOW_LABELS,
)


def make_label_volume(
    volume: Volume,
    label: int,
) -> Volume:

    binary = (
        np.asarray(volume.data)
        == label
    ).astype(np.uint8)

    return Volume(
        data=binary,
        affine=volume.affine,
        spacing=volume.spacing,
        path=volume.path,
    )


def main() -> None:

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "path",
    )

    parser.add_argument(
        "--label",
        type=int,
        required=True,
    )

    parser.add_argument(
        "--trim-mm",
        type=float,
        default=2.0,
    )

    args = parser.parse_args()

    volume = load_nifti(
        args.path
    )

    vessel = make_label_volume(
        volume,
        args.label,
    )

    centerline = extract_centerline(
        vessel
    )

    graph = build_vascular_graph(
        centerline
    )

    main_path = extract_main_path(
        graph
    )

    profile = build_diameter_profile(
        vessel,
        main_path,
        trim_mm=args.trim_mm,
    )

    abbreviation = TOPCOW_LABELS.get(
        args.label,
        str(args.label),
    )

    full_name = TOPCOW_FULL_NAMES.get(
        args.label,
        "Unknown vessel",
    )

    print()
    print("NeuroVasc Workbench")
    print("=" * 60)
    print(f"Vessel:          {abbreviation}")
    print(f"Anatomy:         {full_name}")
    print(f"Main path:       {main_path.length_mm:.2f} mm")
    print(f"Trimmed profile: {profile.length_mm:.2f} mm")
    print(f"Samples:         {len(profile.diameter_mm)}")
    print()
    print(
        f"Mean diameter:   "
        f"{profile.mean_diameter_mm:.2f} mm"
    )
    print(
        f"Median diameter: "
        f"{profile.median_diameter_mm:.2f} mm"
    )
    print(
        f"Minimum diameter:"
        f" {profile.min_diameter_mm:.2f} mm"
    )
    print(
        f"Maximum diameter:"
        f" {profile.max_diameter_mm:.2f} mm"
    )
    print("=" * 60)

    #
    # Save measurements
    #
    output_dir = Path(
        "outputs/measurements/profiles"
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    table = pd.DataFrame(
        {
            "distance_mm": profile.distance_mm,
            "diameter_mm": profile.diameter_mm,
            "x_mm": profile.world_points[:, 0],
            "y_mm": profile.world_points[:, 1],
            "z_mm": profile.world_points[:, 2],
        }
    )

    csv_path = (
        output_dir
        / f"{abbreviation}_diameter_profile.csv"
    )

    table.to_csv(
        csv_path,
        index=False,
    )

    #
    # Plot diameter vs vessel distance.
    #
    plt.figure(
        figsize=(9, 5)
    )

    plt.plot(
        profile.distance_mm,
        profile.diameter_mm,
        marker="o",
        markersize=3,
    )

    plt.xlabel(
        "Distance along centerline (mm)"
    )

    plt.ylabel(
        "Local diameter (mm)"
    )

    plt.title(
        f"{full_name} Diameter Profile"
    )

    plt.grid(
        alpha=0.25
    )

    plt.tight_layout()

    figure_path = (
        output_dir
        / f"{abbreviation}_diameter_profile.png"
    )

    plt.savefig(
        figure_path,
        dpi=200,
    )

    plt.show()

    print()
    print(
        f"Saved profile → {csv_path}"
    )

    print(
        f"Saved figure  → {figure_path}"
    )


if __name__ == "__main__":
    main()