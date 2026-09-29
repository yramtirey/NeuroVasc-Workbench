"""V2 robustness additions without modifying the v1 phantom generator."""
from dataclasses import dataclass

import nibabel as nib
import numpy as np

from scripts.run_caliber_validation import experiment_matrix, SPACINGS
from .geometries import SyntheticVessel, straight_cylinder, oblique_cylinder


@dataclass
class Experiment:
    vessel: SyntheticVessel
    cohort: str = "core"
    phase_id: str = "centered"
    phase_voxels: tuple = (0.0, 0.0, 0.0)
    orientation_id: str = "z"


def phase_shifted_cylinder(diameter_mm, spacing, phase_voxels):
    """Move the voxel grid, not the analytic vessel; affine includes the shift."""
    vessel = straight_cylinder(diameter_mm=diameter_mm, spacing=spacing)
    phase = np.asarray(phase_voxels, dtype=float)
    if phase.shape != (3,) or not np.isfinite(phase).all() or np.any(np.abs(phase) > 0.5):
        raise ValueError("Phase must contain three finite offsets within +/- half a voxel")
    vessel.volume.affine[:3, 3] += phase * spacing
    indices = np.moveaxis(np.indices(vessel.volume.shape), 0, -1)
    world = nib.affines.apply_affine(vessel.volume.affine, indices)
    vessel.volume.data = ((np.abs(world[..., 2]) <= vessel.length_mm / 2 + 1e-10) &
                          (np.hypot(world[..., 0], world[..., 1]) <= diameter_mm / 2 + 1e-10)).astype(np.uint8)
    return vessel


def experiment_matrix_v2():
    for vessel in experiment_matrix():
        yield Experiment(vessel, orientation_id="111" if vessel.geometry == "oblique" else "arc_xz" if vessel.geometry == "curved" else "z")
    for spacing in SPACINGS:
        for diameter in (1, 2, 3, 4, 5):
            for phase in (0.25, 0.5):
                offset = (phase, phase, phase)
                yield Experiment(phase_shifted_cylinder(diameter, spacing, offset), "phase", f"shift_{phase:g}", offset)
        for diameter in (2, 3, 4):
            for direction in ((1, 2, 1), (2, 1, 3)):
                yield Experiment(oblique_cylinder(diameter_mm=diameter, spacing=spacing, direction=direction), "orientation", orientation_id="".join(map(str, direction)))
