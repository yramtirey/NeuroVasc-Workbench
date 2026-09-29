"""Analytic, voxel-center sampled, flat-ended tubes in physical millimeters."""
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from neurovasc.io.volume import Volume


@dataclass
class SyntheticVessel:
    volume: Volume
    geometry: str
    length_mm: float
    proximal_diameter_mm: float
    distal_diameter_mm: float
    direction: np.ndarray
    bend_radius_mm: float | None = None

    def normalized_position(self, points: np.ndarray) -> np.ndarray:
        """Project physical points onto the analytic axis, independent of path order."""
        if self.bend_radius_mm is not None:
            r = self.bend_radius_mm
            along = r * np.arctan2(points[..., 2], points[..., 0] + r)
        else:
            along = points @ self.direction
        return np.clip(along / self.length_mm + 0.5, 0, 1)

    def true_diameter(self, points: np.ndarray) -> np.ndarray:
        u = self.normalized_position(points)
        return self.proximal_diameter_mm + u * (
            self.distal_diameter_mm - self.proximal_diameter_mm
        )

    def metadata(self) -> dict:
        return {
            "geometry": self.geometry,
            "length_mm": self.length_mm,
            "proximal_diameter_mm": self.proximal_diameter_mm,
            "distal_diameter_mm": self.distal_diameter_mm,
            "spacing_mm": list(self.volume.spacing),
            "direction": self.direction.tolist(),
            "bend_radius_mm": self.bend_radius_mm,
            "voxelization": "binary inclusion at voxel centers; boundary included; no partial volume",
            "grid_phase": "odd-sized grid with origin on a voxel center",
            "end_caps": "flat, perpendicular to analytic centerline",
        }


def _generate(geometry, proximal, distal, length, spacing, direction, bend_radius=None):
    spacing = np.asarray(spacing, dtype=float)
    direction = np.asarray(direction, dtype=float)
    if spacing.shape != (3,) or not np.all(np.isfinite(spacing)) or np.any(spacing <= 0):
        raise ValueError("Spacing must contain three positive finite values")
    if not all(np.isfinite(v) and v > 0 for v in (proximal, distal, length)):
        raise ValueError("Diameters and length must be positive and finite")
    if direction.shape != (3,) or not np.all(np.isfinite(direction)) or np.linalg.norm(direction) == 0:
        raise ValueError("Direction must be a nonzero finite 3-vector")
    direction = direction / np.linalg.norm(direction)
    radius = max(proximal, distal) / 2
    if bend_radius is not None:
        if not np.isfinite(bend_radius) or bend_radius <= radius or length / bend_radius >= np.pi:
            raise ValueError("Arc requires bend radius > tube radius and angle < pi")
        angle = length / (2 * bend_radius)
        extent = np.array([bend_radius * (1 - np.cos(angle)), 0, bend_radius * np.sin(angle)])
    else:
        extent = np.abs(direction) * length / 2
    # Background padding prevents the EDT from interacting with the array border.
    half_shape = np.ceil((extent + radius + 3 * spacing.max()) / spacing).astype(int)
    axes = [np.arange(-n, n + 1) * step for n, step in zip(half_shape, spacing)]
    points = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)
    affine = np.diag([*spacing, 1.0])
    affine[:3, 3] = -half_shape * spacing
    volume = Volume(np.zeros(points.shape[:-1], dtype=np.uint8), affine, tuple(spacing), Path(f"{geometry}.nii.gz"))
    vessel = SyntheticVessel(volume, geometry, length, proximal, distal, direction, bend_radius)
    if bend_radius is None:
        along = points @ direction
        radial = np.linalg.norm(points - along[..., None] * direction, axis=-1)
    else:
        along = bend_radius * np.arctan2(points[..., 2], points[..., 0] + bend_radius)
        radial = np.hypot(np.hypot(points[..., 0] + bend_radius, points[..., 2]) - bend_radius, points[..., 1])
    volume.data = ((np.abs(along) <= length / 2 + 1e-10) &
                   (radial <= vessel.true_diameter(points) / 2 + 1e-10)).astype(np.uint8)
    return vessel


def straight_cylinder(diameter_mm=4.0, length_mm=30.0, spacing=(0.5, 0.5, 0.5)):
    return _generate("straight", diameter_mm, diameter_mm, length_mm, spacing, (0, 0, 1))


def oblique_cylinder(diameter_mm=3.0, length_mm=30.0, spacing=(0.5, 0.5, 0.5), direction=(1, 1, 1)):
    return _generate("oblique", diameter_mm, diameter_mm, length_mm, spacing, direction)


def tapered_vessel(proximal_diameter_mm=4.0, distal_diameter_mm=2.0, length_mm=30.0, spacing=(0.5, 0.5, 0.5)):
    return _generate("tapered", proximal_diameter_mm, distal_diameter_mm, length_mm, spacing, (0, 0, 1))


def curved_vessel(diameter_mm=3.0, length_mm=30.0, spacing=(0.5, 0.5, 0.5), bend_radius_mm=24.0):
    return _generate("curved", diameter_mm, diameter_mm, length_mm, spacing, (0, 0, 1), bend_radius_mm)
