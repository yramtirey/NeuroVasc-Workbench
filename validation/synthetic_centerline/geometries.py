"""Reuse the exact caliber-v2 matrix, plus an explicitly branching phantom."""
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from neurovasc.io.volume import Volume
from validation.synthetic_caliber.experiments_v2 import Experiment, experiment_matrix_v2


@dataclass
class Bifurcation:
    volume: Volume
    branches: dict[str, np.ndarray]
    diameters: dict[str, float]
    geometry: str = 'bifurcation'
    proximal_diameter_mm: float = 4.0
    distal_diameter_mm: float = 3.0

    @property
    def length_mm(self):
        return sum(sorted(np.linalg.norm(segment[1] - segment[0]) for segment in self.branches.values())[-2:])

    def metadata(self):
        return {
            'geometry': self.geometry, 'branches': {k: v.tolist() for k, v in self.branches.items()},
            'branch_diameters_mm': self.diameters, 'junction_mm': [0, 0, 0],
            'spacing_mm': list(self.volume.spacing), 'true_principal_length_mm': self.length_mm,
            'voxelization': 'union of three flat-ended analytic tubes; binary voxel-center inclusion',
        }


def y_bifurcation(spacing=(0.5, 0.5, 0.5), parent_diameter_mm=4.0,
                  daughter_diameters_mm=(3.0, 3.0), parent_length_mm=22.0,
                  daughter_length_mm=18.0, branch_angle_deg=35.0):
    spacing = np.asarray(spacing, dtype=float)
    values = [parent_diameter_mm, *daughter_diameters_mm, parent_length_mm, daughter_length_mm]
    if spacing.shape != (3,) or not np.isfinite(spacing).all() or np.any(spacing <= 0) or len(daughter_diameters_mm) != 2 or not all(np.isfinite(x) and x > 0 for x in values) or not 0 < branch_angle_deg < 90:
        raise ValueError('invalid_bifurcation_parameters')
    angle = np.deg2rad(branch_angle_deg)
    branches = {
        'parent': np.array([[0, 0, 0], [0, 0, -parent_length_mm]], dtype=float),
        'left': np.array([[0, 0, 0], [-daughter_length_mm * np.sin(angle), 0, daughter_length_mm * np.cos(angle)]]),
        'right': np.array([[0, 0, 0], [daughter_length_mm * np.sin(angle), 0, daughter_length_mm * np.cos(angle)]]),
    }
    diameters = dict(zip(branches, [parent_diameter_mm, *daughter_diameters_mm]))
    extent = np.max(np.abs(np.concatenate(list(branches.values()))), axis=0)
    half = np.ceil((extent + max(diameters.values()) / 2 + 3 * spacing.max()) / spacing).astype(int)
    axes = [np.arange(-n, n + 1) * step for n, step in zip(half, spacing)]
    points = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1)
    mask = np.zeros(points.shape[:-1], dtype=bool)
    for name, segment in branches.items():
        length = np.linalg.norm(segment[1])
        direction = segment[1] / length
        along = points @ direction
        radial = np.linalg.norm(points - along[..., None] * direction, axis=-1)
        mask |= (along >= -1e-10) & (along <= length + 1e-10) & (radial <= diameters[name] / 2 + 1e-10)
    affine = np.diag([*spacing, 1.0])
    affine[:3, 3] = -half * spacing
    volume = Volume(mask.astype(np.uint8), affine, tuple(spacing), Path('synthetic_y.nii.gz'))
    return Bifurcation(volume, branches, diameters, proximal_diameter_mm=parent_diameter_mm, distal_diameter_mm=daughter_diameters_mm[0])


def matrix():
    for index, experiment in enumerate(experiment_matrix_v2(), 1):
        vessel = experiment.vessel
        # Identical IDs allow a direct join with caliber-v2 failures.
        geometry_id = f'{index:02d}_{vessel.geometry}_{vessel.proximal_diameter_mm:g}_{experiment.phase_id}_{experiment.orientation_id}_z{vessel.volume.spacing[2]:g}'
        yield geometry_id, experiment
    for index, spacing in enumerate(((0.5, 0.5, 0.5), (0.46875, 0.46875, 0.8)), 59):
        yield f'{index}_bifurcation_z{spacing[2]}', Experiment(y_bifurcation(spacing), cohort='branching', orientation_id='y_xz')
