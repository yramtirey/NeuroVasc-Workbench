from pathlib import Path

import numpy as np
import pandas as pd

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.diameter_profile import build_diameter_profile
from neurovasc.geometry.main_path import extract_main_path
from neurovasc.geometry.mesh import mask_to_mesh
from neurovasc.geometry.profile_analysis import analyze_profile
from neurovasc.graph.vascular_graph import build_vascular_graph
from neurovasc.io.volume import Volume
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


def analyze_vessel(
    volume: Volume,
    label: int,
    trim_mm: float = 2.0,
) -> dict:

    vessel = make_label_volume(
        volume,
        label,
    )

    voxel_count = int(
        np.count_nonzero(vessel.data)
    )

    voxel_volume_mm3 = (
        voxel_count
        * volume.voxel_volume_mm3
    )

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
        trim_mm=trim_mm,
    )

    analysis = analyze_profile(
        profile
    )

    min_world = profile.world_points[
        analysis.minimum_index
    ]

    return {
        "label": label,

        "abbreviation":
            TOPCOW_LABELS.get(
                label,
                f"Label-{label}",
            ),

        "name":
            TOPCOW_FULL_NAMES.get(
                label,
                "Unknown",
            ),

        "voxel_count":
            voxel_count,

        "voxel_volume_mm3":
            float(voxel_volume_mm3),

        "mesh_volume_mm3":
            float(mesh.volume_mm3),

        "surface_area_mm2":
            float(mesh.surface_area_mm2),

        "centerline_components":
            int(graph.n_components),

        "main_path_length_mm":
            float(main_path.length_mm),

        "mean_diameter_mm":
            float(profile.mean_diameter_mm),

        "median_diameter_mm":
            float(profile.median_diameter_mm),

        "minimum_diameter_mm":
            float(analysis.minimum_diameter_mm),

        "maximum_diameter_mm":
            float(profile.max_diameter_mm),

        "reference_diameter_mm":
            float(analysis.reference_diameter_mm),

        "candidate_reduction_percent":
            float(
                analysis.maximum_reduction_percent
            ),

        "minimum_location_mm":
            float(
                analysis.minimum_distance_mm
            ),

        "minimum_x_mm":
            float(min_world[0]),

        "minimum_y_mm":
            float(min_world[1]),

        "minimum_z_mm":
            float(min_world[2]),
    }


def analyze_case(
    volume: Volume,
    trim_mm: float = 2.0,
) -> tuple[pd.DataFrame, list[str]]:

    labels = sorted(
        int(value)
        for value in np.unique(volume.data)
        if int(value) != 0
    )

    records = []
    warnings = []

    for label in labels:

        try:

            record = analyze_vessel(
                volume,
                label,
                trim_mm=trim_mm,
            )

            records.append(
                record
            )

        except Exception as error:

            name = TOPCOW_LABELS.get(
                label,
                str(label),
            )

            warnings.append(
                f"{name}: {error}"
            )

    table = pd.DataFrame(
        records
    )

    return table, warnings