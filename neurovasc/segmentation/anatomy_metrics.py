from dataclasses import dataclass

import numpy as np
from neurovasc.geometry.caliber import measure_caliber
from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.mesh import mask_to_mesh
from neurovasc.graph.vascular_graph import build_vascular_graph
from neurovasc.io.volume import Volume
from neurovasc.segmentation.labels import (
    TOPCOW_FULL_NAMES,
    TOPCOW_LABELS,
)


@dataclass
class AnatomyMetrics:
    label: int
    abbreviation: str
    name: str
    voxel_count: int
    voxel_volume_mm3: float
    mesh_volume_mm3: float
    surface_area_mm2: float
    centerline_length_mm: float
    centerline_points: int
    connected_components: int
    mean_diameter_mm: float
    median_diameter_mm: float
    min_diameter_mm: float
    max_diameter_mm: float
    p10_diameter_mm: float
    p90_diameter_mm: float


def binary_label_volume(
    volume: Volume,
    label: int,
) -> Volume:
    """Create a binary NeuroVasc Volume for one anatomical label."""

    binary = (
        np.asarray(volume.data) == label
    ).astype(np.uint8)

    return Volume(
        data=binary,
        affine=volume.affine,
        spacing=volume.spacing,
        path=volume.path,
    )


def measure_anatomical_label(
    volume: Volume,
    label: int,
) -> AnatomyMetrics:
    """Measure geometry for one TopCoW anatomical vessel label."""

    label_volume = binary_label_volume(
        volume,
        label,
    )

    voxel_count = int(
        np.count_nonzero(
            label_volume.data
        )
    )

    if voxel_count == 0:
        raise ValueError(
            f"Label {label} is not present."
        )

    voxel_volume_mm3 = (
        voxel_count
        * volume.voxel_volume_mm3
    )

    mesh = mask_to_mesh(
        label_volume
    )

    centerline = extract_centerline(
        label_volume
    )
    caliber = measure_caliber(
        label_volume,
        centerline=centerline,
    )

    graph = build_vascular_graph(
        centerline
    )

    return AnatomyMetrics(
        label=label,
        abbreviation=TOPCOW_LABELS.get(
            label,
            f"Label-{label}",
        ),
        name=TOPCOW_FULL_NAMES.get(
            label,
            "Unknown",
        ),
        voxel_count=voxel_count,
        voxel_volume_mm3=float(
            voxel_volume_mm3
        ),
        mesh_volume_mm3=float(
            mesh.volume_mm3
        ),
        surface_area_mm2=float(
            mesh.surface_area_mm2
        ),
        centerline_length_mm=float(
            graph.total_length_mm
        ),
        centerline_points=int(
            centerline.n_points
        ),
        connected_components=int(
            graph.n_components
        ),
        mean_diameter_mm=caliber.mean_diameter_mm,
        median_diameter_mm=caliber.median_diameter_mm,
        min_diameter_mm=caliber.min_diameter_mm,
        max_diameter_mm=caliber.max_diameter_mm,
        p10_diameter_mm=caliber.p10_diameter_mm,
        p90_diameter_mm=caliber.p90_diameter_mm,
    )