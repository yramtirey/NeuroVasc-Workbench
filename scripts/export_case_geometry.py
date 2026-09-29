import argparse
from pathlib import Path

import numpy as np

from neurovasc.geometry.mesh import mask_to_mesh
from neurovasc.io.volume import Volume, load_nifti
from neurovasc.segmentation.labels import TOPCOW_LABELS


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
        help="TopCoW segmentation",
    )

    parser.add_argument(
        "--case-id",
        default="case_001",
    )

    args = parser.parse_args()

    volume = load_nifti(args.path)

    labels = sorted(
        int(value)
        for value in np.unique(volume.data)
        if int(value) != 0
    )

    output_dir = (
        Path("outputs")
        / "cases"
        / args.case_id
        / "geometry"
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print()
    print("NeuroVasc Geometry Export")
    print("=" * 60)

    for label in labels:

        abbreviation = TOPCOW_LABELS.get(
            label,
            f"label_{label}",
        )

        print(
            f"Exporting {abbreviation}..."
        )

        label_volume = make_label_volume(
            volume,
            label,
        )

        vessel = mask_to_mesh(
            label_volume
        )

        output = (
            output_dir
            / f"label_{label}.vtp"
        )

        vessel.mesh.save(
            output,
            binary=False,
        )

    print()
    print(
        f"Geometry saved → {output_dir}"
    )

    print("=" * 60)


if __name__ == "__main__":
    main()