from dataclasses import dataclass

import numpy as np

from neurovasc.geometry.caliber import compute_radius_map
from neurovasc.geometry.main_path import MainCenterlinePath
from neurovasc.io.volume import Volume


@dataclass
class DiameterProfile:
    distance_mm: np.ndarray
    diameter_mm: np.ndarray
    world_points: np.ndarray

    @property
    def length_mm(self) -> float:
        if len(self.distance_mm) == 0:
            return 0.0

        return float(self.distance_mm[-1])

    @property
    def mean_diameter_mm(self) -> float:
        return float(
            np.mean(self.diameter_mm)
        )

    @property
    def median_diameter_mm(self) -> float:
        return float(
            np.median(self.diameter_mm)
        )

    @property
    def min_diameter_mm(self) -> float:
        return float(
            np.min(self.diameter_mm)
        )

    @property
    def max_diameter_mm(self) -> float:
        return float(
            np.max(self.diameter_mm)
        )


def build_diameter_profile(
    volume: Volume,
    path: MainCenterlinePath,
    trim_mm: float = 2.0,
    min_samples: int = 5,
) -> DiameterProfile:
    """
    Sample local vessel diameter along an ordered centerline.

    Endpoint trimming is applied when enough samples remain.
    For short vessels, NeuroVasc falls back to the full path
    rather than returning an empty profile.
    """

    radius_map = compute_radius_map(
        volume
    )

    voxels = path.voxel_points

    radii = radius_map[
        voxels[:, 0],
        voxels[:, 1],
        voxels[:, 2],
    ]

    diameters = (
        2.0
        * np.asarray(
            radii,
            dtype=float,
        )
    )

    distance = np.asarray(
        path.distance_mm,
        dtype=float,
    )

    world_points = np.asarray(
        path.world_points,
        dtype=float,
    )

    if len(distance) == 0:
        raise ValueError(
            "Main centerline contains no samples."
        )

    total_length = path.length_mm

    if trim_mm > 0 and total_length > 2 * trim_mm:

        keep = (
            (distance >= trim_mm)
            &
            (distance <= total_length - trim_mm)
        )

        #
        # Only trim if enough samples survive.
        #
        if np.count_nonzero(keep) >= min_samples:

            distance = distance[keep]
            diameters = diameters[keep]
            world_points = world_points[keep]

            distance = (
                distance
                - distance[0]
            )

    return DiameterProfile(
        distance_mm=distance,
        diameter_mm=diameters,
        world_points=world_points,
    )