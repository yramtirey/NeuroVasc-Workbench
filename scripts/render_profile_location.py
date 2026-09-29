import argparse

import numpy as np
import pyvista as pv

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.diameter_profile import build_diameter_profile
from neurovasc.geometry.main_path import extract_main_path
from neurovasc.geometry.mesh import mask_to_mesh
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


def make_polyline(
    points: np.ndarray,
) -> pv.PolyData:

    if len(points) < 2:
        return pv.PolyData()

    return pv.lines_from_points(
        points,
        close=False,
    )


def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Map a caliber-profile minimum "
            "back into 3D vessel anatomy."
        )
    )

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

    #
    # Load case
    #
    volume = load_nifti(
        args.path
    )

    vessel = make_label_volume(
        volume,
        args.label,
    )

    #
    # Geometry pipeline
    #
    mesh = mask_to_mesh(
        vessel
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

    analysis = analyze_profile(
        profile,
        sigma=args.sigma,
    )

    #
    # Locate minimum in physical coordinates
    #
    min_index = analysis.minimum_index

    minimum_world = (
        profile.world_points[min_index]
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
    print("=" * 65)

    print(
        f"Vessel:                "
        f"{abbreviation}"
    )

    print(
        f"Anatomy:               "
        f"{full_name}"
    )

    print(
        f"Minimum caliber:       "
        f"{analysis.minimum_diameter_mm:.2f} mm"
    )

    print(
        f"Distance along vessel: "
        f"{analysis.minimum_distance_mm:.2f} mm"
    )

    print(
        "World coordinate:      "
        f"({minimum_world[0]:.2f}, "
        f"{minimum_world[1]:.2f}, "
        f"{minimum_world[2]:.2f}) mm"
    )

    print(
        f"Reference caliber:     "
        f"{analysis.reference_diameter_mm:.2f} mm"
    )

    print(
        f"Relative reduction:    "
        f"{analysis.maximum_reduction_percent:.1f}%"
    )

    print()
    print(
        "Candidate geometric caliber reduction only — "
        "not a clinical diagnosis."
    )

    print("=" * 65)

    #
    # Ordered centerline
    #
    centerline_poly = make_polyline(
        profile.world_points
    )

    #
    # Attach diameter measurements to points
    #
    centerline_poly[
        "diameter_mm"
    ] = profile.diameter_mm

    #
    # Minimum-caliber marker
    #
    minimum_point = pv.PolyData(
        np.asarray(
            [minimum_world],
            dtype=float,
        )
    )

    #
    # Viewer
    #
    plotter = pv.Plotter()

    # Transparent artery
    plotter.add_mesh(
        mesh.mesh,
        opacity=0.18,
        smooth_shading=True,
    )

    # Centerline colored by caliber
    plotter.add_mesh(
        centerline_poly,
        scalars="diameter_mm",
        line_width=8,
        scalar_bar_args={
            "title": "Diameter (mm)"
        },
    )

    # Candidate minimum
    plotter.add_mesh(
        minimum_point,
        point_size=28,
        render_points_as_spheres=True,
    )

    # Text label at candidate point
    plotter.add_point_labels(
        np.asarray(
            [minimum_world]
        ),
        [
            (
                f"Minimum caliber\n"
                f"{analysis.minimum_diameter_mm:.2f} mm\n"
                f"{analysis.minimum_distance_mm:.1f} mm along path"
            )
        ],
        font_size=12,
        point_size=0,
        always_visible=True,
    )

    plotter.add_axes()

    plotter.camera_position = "iso"

    plotter.show()


if __name__ == "__main__":
    main()