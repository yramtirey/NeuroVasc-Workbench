import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.diameter_profile import build_diameter_profile
from neurovasc.geometry.main_path import extract_main_path
from neurovasc.geometry.profile_analysis import analyze_profile
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
        np.asarray(volume.data) == label
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

    parser.add_argument(
        "--sigma",
        type=float,
        default=1.25,
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

    path = extract_main_path(
        graph
    )

    profile = build_diameter_profile(
        vessel,
        path,
        trim_mm=args.trim_mm,
    )

    analysis = analyze_profile(
        profile,
        sigma=args.sigma,
    )

    abbreviation = TOPCOW_LABELS.get(
        args.label,
        str(args.label),
    )

    name = TOPCOW_FULL_NAMES.get(
        args.label,
        "Unknown vessel",
    )

    print()
    print("NeuroVasc Workbench")
    print("=" * 65)

    print(f"Vessel:                 {abbreviation}")
    print(f"Anatomy:                {name}")

    print(
        f"Reference diameter:     "
        f"{analysis.reference_diameter_mm:.2f} mm"
    )

    print(
        f"Minimum smoothed caliber:"
        f" {analysis.minimum_diameter_mm:.2f} mm"
    )

    print(
        f"Minimum location:       "
        f"{analysis.minimum_distance_mm:.2f} mm"
    )

    print(
        f"Relative reduction:     "
        f"{analysis.maximum_reduction_percent:.1f}%"
    )

    print()
    print(
        "NOTE: Reduction is a segmentation-derived geometric "
        "metric, not a clinical stenosis diagnosis."
    )

    print("=" * 65)

    plt.figure(
        figsize=(10, 5)
    )

    plt.plot(
        analysis.distance_mm,
        analysis.raw_diameter_mm,
        marker="o",
        markersize=3,
        alpha=0.45,
        label="Raw diameter",
    )

    plt.plot(
        analysis.distance_mm,
        analysis.smooth_diameter_mm,
        linewidth=2.5,
        label="Smoothed diameter",
    )

    plt.axhline(
        analysis.reference_diameter_mm,
        linestyle="--",
        label="Reference caliber",
    )

    plt.scatter(
        [analysis.minimum_distance_mm],
        [analysis.minimum_diameter_mm],
        s=80,
        zorder=5,
        label="Minimum caliber",
    )

    plt.xlabel(
        "Distance along centerline (mm)"
    )

    plt.ylabel(
        "Diameter (mm)"
    )

    plt.title(
        f"{name} — Caliber Analysis"
    )

    plt.legend()

    plt.grid(
        alpha=0.25
    )

    plt.tight_layout()

    output_dir = Path(
        "outputs/measurements/profiles"
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output = (
        output_dir
        / f"{abbreviation}_caliber_analysis.png"
    )

    plt.savefig(
        output,
        dpi=200,
    )

    plt.show()

    print(
        f"\nSaved → {output}"
    )


if __name__ == "__main__":
    main()