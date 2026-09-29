from dataclasses import dataclass

import nibabel as nib
import numpy as np
import pyvista as pv
from skimage.measure import marching_cubes

from neurovasc.io.volume import Volume


@dataclass
class VesselMesh:
    """Triangular surface representation of a vessel mask."""

    mesh: pv.PolyData

    @property
    def surface_area_mm2(self) -> float:
        return float(self.mesh.area)

    @property
    def volume_mm3(self) -> float:
        return float(self.mesh.volume)

    @property
    def n_vertices(self) -> int:
        return int(self.mesh.n_points)

    @property
    def n_faces(self) -> int:
        return int(self.mesh.n_cells)


def mask_to_mesh(
    volume: Volume,
    threshold: float = 0.5,
) -> VesselMesh:
    """
    Convert a 3D binary vessel mask into a surface mesh.
    """

    mask = np.asarray(volume.data) > threshold

    if not np.any(mask):
        raise ValueError(
            "Cannot create mesh: vessel mask is empty."
        )

    vertices, faces, _, _ = marching_cubes(
        mask.astype(np.float32),
        level=0.5,
    )

    # Convert voxel coordinates -> physical/world coordinates.
    world_vertices = nib.affines.apply_affine(
        volume.affine,
        vertices,
    )

    # PyVista face format:
    # [3, v0, v1, v2, 3, v0, v1, v2, ...]
    pyvista_faces = np.column_stack(
        (
            np.full(len(faces), 3),
            faces,
        )
    ).ravel()

    mesh = pv.PolyData(
        world_vertices,
        pyvista_faces,
    )

    mesh.clean(inplace=True)

    return VesselMesh(mesh=mesh)