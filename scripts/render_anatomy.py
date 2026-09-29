import argparse

import numpy as np
import pyvista as pv

from neurovasc.geometry.mesh import mask_to_mesh
from neurovasc.io.volume import Volume, load_nifti
from neurovasc.segmentation.labels import (
    TOPCOW_LABELS,
    TOPCOW_FULL_NAMES,
)


def make_label_volume(
    volume: Volume,
    label: int,
) -> Volume:
    """
    Create a binary Volume containing only one anatomical label.
    """

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

    parser = argparse.ArgumentParser(
        description=(
            "Render anatomical TopCoW vessel labels "
            "with NeuroVasc Workbench."
        )
    )

    parser.add_argument(
        "path",
        help="Path to TopCoW multiclass segmentation.",
    )

    args = parser.parse_args()

    volume = load_nifti(args.path)

    labels = [
        int(label)
        for label in np.unique(volume.data)
        if int(label) != 0
    ]

    print()
    print("NeuroVasc Workbench")
    print("=" * 65)
    print(f"Case: {volume.path.name}")
    print()
    print("Anatomical vessels detected")
    print("-" * 65)

    plotter = pv.Plotter()

    plotter.set_color_cycler("default")

    legend_entries = []

    for label in labels:

        name = TOPCOW_LABELS.get(
            label,
            f"Unknown-{label}",
        )

        full_name = TOPCOW_FULL_NAMES.get(
            label,
            "Unknown vessel",
        )

        label_volume = make_label_volume(
            volume,
            label,
        )

        voxel_count = int(
            np.count_nonzero(
                label_volume.data
            )
        )

        print(
            f"{label:2d}  "
            f"{name:<8} "
            f"{voxel_count:>6,} voxels  "
            f"{full_name}"
        )

        vessel = mask_to_mesh(
            label_volume
        )

        actor = plotter.add_mesh(
            vessel.mesh,
            smooth_shading=True,
            opacity=0.90,
            label=name,
        )

    print("-" * 65)

    plotter.add_legend(
        bcolor="white",
        face="circle",
        size=(0.18, 0.35),
    )

    plotter.add_axes()

    plotter.camera_position = "iso"

    plotter.show()


if __name__ == "__main__":
    main()