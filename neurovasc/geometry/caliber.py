from dataclasses import dataclass

import numpy as np
from scipy.ndimage import distance_transform_edt

from neurovasc.geometry.centerline import (
    Centerline,
    extract_centerline,
)
from neurovasc.io.volume import Volume


@dataclass
class CaliberProfile:
    """Diameter measurements sampled along a vessel centerline."""

    world_points: np.ndarray
    voxel_points: np.ndarray
    radius_mm: np.ndarray
    diameter_mm: np.ndarray

    @property
    def n_samples(self) -> int:
        return len(self.diameter_mm)

    @property
    def mean_diameter_mm(self) -> float:
        return float(np.mean(self.diameter_mm))

    @property
    def median_diameter_mm(self) -> float:
        return float(np.median(self.diameter_mm))

    @property
    def min_diameter_mm(self) -> float:
        return float(np.min(self.diameter_mm))

    @property
    def max_diameter_mm(self) -> float:
        return float(np.max(self.diameter_mm))

    @property
    def p10_diameter_mm(self) -> float:
        return float(
            np.percentile(
                self.diameter_mm,
                10,
            )
        )

    @property
    def p90_diameter_mm(self) -> float:
        return float(
            np.percentile(
                self.diameter_mm,
                90,
            )
        )


def compute_radius_map(
    volume: Volume,
    threshold: float = 0.5,
) -> np.ndarray:
    """
    Compute distance from each vessel voxel to the nearest
    background voxel in physical millimeters.
    """

    mask = (
        np.asarray(volume.data)
        > threshold
    )

    if not np.any(mask):
        raise ValueError(
            "Cannot compute caliber: vessel mask is empty."
        )

    radius_map = distance_transform_edt(
        mask,
        sampling=volume.spacing,
    )

    return radius_map


def measure_caliber(
    volume: Volume,
    centerline: Centerline | None = None,
) -> CaliberProfile:
    """
    Measure local vessel diameter along the skeleton centerline.
    """

    if centerline is None:
        centerline = extract_centerline(
            volume
        )

    radius_map = compute_radius_map(
        volume
    )

    voxels = centerline.voxel_points

    radii = radius_map[
        voxels[:, 0],
        voxels[:, 1],
        voxels[:, 2],
    ]

    diameters = 2.0 * radii

    return CaliberProfile(
        world_points=centerline.world_points,
        voxel_points=centerline.voxel_points,
        radius_mm=np.asarray(
            radii,
            dtype=float,
        ),
        diameter_mm=np.asarray(
            diameters,
            dtype=float,
        ),
    )