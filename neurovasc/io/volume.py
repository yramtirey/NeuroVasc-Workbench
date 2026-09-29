from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np


@dataclass
class Volume:
    """A 3D medical imaging volume and its spatial metadata."""

    data: np.ndarray
    affine: np.ndarray
    spacing: tuple[float, float, float]
    path: Path

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.data.shape

    @property
    def voxel_volume_mm3(self) -> float:
        """Physical volume represented by one voxel."""
        return float(np.prod(self.spacing))

    @property
    def physical_size_mm(self) -> tuple[float, float, float]:
        """Physical field of view along each image axis."""
        return tuple(
            float(size * spacing)
            for size, spacing in zip(
                self.shape,
                self.spacing,
            )
        )


def load_nifti(path: str | Path) -> Volume:
    """
    Load a 3D NIfTI medical image.

    Parameters
    ----------
    path:
        Path to a .nii or .nii.gz file.

    Returns
    -------
    Volume
        Image data and spatial metadata.
    """

    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(
            f"Volume not found: {path}"
        )

    if not (
        path.name.endswith(".nii")
        or path.name.endswith(".nii.gz")
    ):
        raise ValueError(
            "Expected a .nii or .nii.gz file."
        )

    image = nib.load(path)

    data = np.asarray(image.dataobj)

    if data.ndim != 3:
        raise ValueError(
            f"Expected a 3D image, "
            f"but received shape {data.shape}."
        )

    spacing = tuple(
        float(value)
        for value in image.header.get_zooms()[:3]
    )

    affine = np.asarray(
        image.affine,
        dtype=float,
    )

    return Volume(
        data=data,
        affine=affine,
        spacing=spacing,
        path=path,
    )