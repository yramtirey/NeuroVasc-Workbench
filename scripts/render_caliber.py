import argparse

import numpy as np
import pyvista as pv

from neurovasc.geometry.caliber import measure_caliber
from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.mesh import mask_to_mesh
from neurovasc.io.volume import Volume, load_nifti
from neurovasc.segmentation.labels import TOPCOW_LABELS


def label_volume(
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

    parser = argparse.ArgumentParser(
        description=(
            "Render vessel centerlines colored "
            "by local diameter."
        )
    )

    parser.add_argument(
        "path",
        help="TopCoW segmentation",
    )

    parser.add_argument(
        "--label",
        type=int,
        default=4,
        help="TopCoW vessel label to inspect",
    )

    args = parser.parse_args()

    volume = load_nifti(
        args.path
    )

    vessel_volume = label_volume(
        volume,
        args.label,
    )

    centerline = extract_centerline(
        vessel_volume
    )

    caliber = measure_caliber(
        vessel_volume,
        centerline=centerline,
    )

    mesh = mask_to_mesh(
        vessel_volume
    )

    print()
    print("NeuroVasc Workbench")
    print("=" * 60)

    print(
        f"Vessel: "
        f"{TOPCOW_LABELS.get(args.label, args.label)}"
    )

    print(
        f"Centerline samples: "
        f"{caliber.n_samples}"
    )

    print(
        f"Mean diameter:   "
        f"{caliber.mean_diameter_mm:.2f} mm"
    )

    print(
        f"Median diameter: "
        f"{caliber.median_diameter_mm:.2f} mm"
    )

    print(
        f"Minimum diameter:"
        f" {caliber.min_diameter_mm:.2f} mm"
    )

    print(
        f"Maximum diameter:"
        f" {caliber.max_diameter_mm:.2f} mm"
    )

    print("=" * 60)

    points = pv.PolyData(
        caliber.world_points
    )

    points["diameter_mm"] = (
        caliber.diameter_mm
    )

    plotter = pv.Plotter()

    plotter.add_mesh(
        mesh.mesh,
        opacity=0.18,
        smooth_shading=True,
    )

    plotter.add_mesh(
        points,
        scalars="diameter_mm",
        point_size=12,
        render_points_as_spheres=True,
        scalar_bar_args={
            "title": "Diameter (mm)"
        },
    )

    plotter.add_axes()

    plotter.camera_position = "iso"

    plotter.show()


if __name__ == "__main__":
    main()