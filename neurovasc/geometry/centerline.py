from dataclasses import dataclass

import nibabel as nib
import numpy as np
from skimage.morphology import skeletonize

from neurovasc.io.volume import Volume


@dataclass
class Centerline:
    """One-voxel-wide centerline extracted from a vessel mask."""

    mask: np.ndarray
    voxel_points: np.ndarray
    world_points: np.ndarray

    @property
    def n_points(self) -> int:
        return int(len(self.voxel_points))


def extract_centerline(
    volume: Volume,
    threshold: float = 0.5,
) -> Centerline:
    """
    Skeletonize a binary vascular volume.

    The resulting skeleton is one voxel wide and represents
    the approximate centerline of the vascular network.
    """

    vessel_mask = np.asarray(volume.data) > threshold

    if not np.any(vessel_mask):
        raise ValueError(
            "Cannot extract centerline: vessel mask is empty."
        )

    skeleton = skeletonize(
        vessel_mask,
        method="lee",
    ).astype(bool)

    voxel_points = np.argwhere(skeleton)

    if len(voxel_points) == 0:
        raise ValueError(
            "Skeletonization produced no centerline points."
        )

    world_points = nib.affines.apply_affine(
        volume.affine,
        voxel_points,
    )

    return Centerline(
        mask=skeleton,
        voxel_points=voxel_points,
        world_points=world_points,
    )