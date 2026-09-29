"""Experimental mesh-plane equivalent caliber. Not used by case analysis.

A single closed contour enclosing the sampling point is required. Multiple or
open contours are flagged, not summed or silently replaced with EDT estimates.
"""
from dataclasses import dataclass

import numpy as np
import pyvista as pv

from neurovasc.geometry.main_path import MainCenterlinePath
from neurovasc.geometry.mesh import mask_to_mesh
from neurovasc.io.volume import Volume


@dataclass
class CrossSectionalProfile:
    distance_mm: np.ndarray
    world_points: np.ndarray
    tangent_vectors: np.ndarray
    area_mm2: np.ndarray
    diameter_mm: np.ndarray
    valid: np.ndarray
    failure_reasons: list[str]


def equivalent_diameter(area_mm2):
    """Convert positive finite area to equivalent circular diameter; else NaN."""
    area = np.asarray(area_mm2, dtype=float)
    return 2 * np.sqrt(np.where(np.isfinite(area) & (area > 0), area, np.nan) / np.pi)


def estimate_tangents(world_points, window_mm=6.0):
    """Local linear fit of xyz against physical arc length over +/- window/2.

    End windows are clipped (one-sided fits). If fewer than two distinct points
    fall in the window, use the two nearest arc-length samples. Coincident,
    nonfinite or directionless neighborhoods return NaN vectors.
    """
    points = np.asarray(world_points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("world_points must have shape (N, 3)")
    if not np.isfinite(window_mm) or window_mm <= 0:
        raise ValueError("window_mm must be positive and finite")
    tangents = np.full_like(points, np.nan)
    if len(points) < 2 or not np.isfinite(points).all():
        return tangents
    distance = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    for i in range(len(points)):
        chosen = np.flatnonzero(np.abs(distance - distance[i]) <= window_mm / 2 + 1e-10)
        if len(chosen) < 2:
            chosen = np.argsort(np.abs(distance - distance[i]), kind="stable")[:2]
        s = distance[chosen] - distance[chosen].mean()
        denominator = s @ s
        if denominator <= 1e-20:
            continue
        slope = s @ (points[chosen] - points[chosen].mean(axis=0)) / denominator
        norm = np.linalg.norm(slope)
        if norm > 1e-12:
            tangents[i] = slope / norm
    return tangents


def plane_basis(tangent):
    """Return orthonormal (u,v) with u cross v parallel to the tangent."""
    normal = np.asarray(tangent, dtype=float)
    if normal.shape != (3,) or not np.isfinite(normal).all() or np.linalg.norm(normal) < 1e-12:
        raise ValueError("Plane requires a finite nonzero tangent")
    normal = normal / np.linalg.norm(normal)
    axis = np.eye(3)[np.argmin(np.abs(normal))]
    u = np.cross(normal, axis)
    u /= np.linalg.norm(u)
    return u, np.cross(normal, u)


def section_area(mesh: pv.PolyData, origin, tangent, tolerance_mm=1e-6):
    """Return (area, reason); invalid area is NaN. No rendering is required."""
    if not np.isfinite(tolerance_mm) or tolerance_mm <= 0:
        raise ValueError("Contour tolerance must be positive and finite")
    origin = np.asarray(origin, dtype=float)
    if origin.shape != (3,) or not np.isfinite(origin).all():
        return np.nan, "invalid_origin"
    try:
        u, v = plane_basis(tangent)
    except ValueError:
        return np.nan, "invalid_tangent"
    try:
        section = mesh.slice(normal=np.cross(u, v), origin=origin).clean(
            tolerance=tolerance_mm, absolute=True,
        )
    except (ValueError, RuntimeError) as error:
        return np.nan, f"slice_failed:{type(error).__name__}"
    lines = section.lines
    edges = set()
    cursor = 0
    while cursor < len(lines):
        count = int(lines[cursor])
        ids = lines[cursor + 1:cursor + 1 + count]
        for a, b in zip(ids[:-1], ids[1:]):
            if a != b:
                edges.add(tuple(sorted((int(a), int(b)))))
        cursor += count + 1
    if not edges:
        return np.nan, "empty_intersection"
    neighbors = {}
    for a, b in sorted(edges):
        neighbors.setdefault(a, set()).add(b)
        neighbors.setdefault(b, set()).add(a)
    if any(len(connected) != 2 for connected in neighbors.values()):
        return np.nan, "open_or_branched_contour"
    unseen, contours = set(neighbors), []
    while unseen:
        start = min(unseen)
        loop, previous, current = [], None, start
        while current not in loop:
            loop.append(current)
            unseen.discard(current)
            nxt = min(neighbors[current] - ({previous} if previous is not None else set()))
            previous, current = current, nxt
        if current != start or len(loop) < 3:
            return np.nan, "degenerate_contour"
        contours.append(loop)
    if len(contours) != 1:
        return np.nan, "multiple_contours"
    relative = np.asarray(section.points)[contours[0]] - origin
    polygon = np.column_stack((relative @ u, relative @ v))
    x, y = polygon.T
    area = abs(float(x @ np.roll(y, -1) - y @ np.roll(x, -1))) / 2
    if not np.isfinite(area) or area <= tolerance_mm ** 2:
        return np.nan, "degenerate_area"
    # Ray crossing at the sample origin. Reject unrelated/off-vessel contours.
    inside = False
    for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
        if (a[1] > 0) != (b[1] > 0):
            crossing_x = a[0] + (b[0] - a[0]) * (-a[1]) / (b[1] - a[1])
            if crossing_x > 0:
                inside = not inside
    if not inside:
        return np.nan, "origin_outside_contour"
    return area, "ok"


def measure_cross_sectional_caliber(volume: Volume, path: MainCenterlinePath,
                                     tangent_window_mm=6.0) -> CrossSectionalProfile:
    """Build the production 0.5 isosurface once, then cut at each path sample.

    Tangents are fit on the full path before validation endpoint trimming.
    All sample slots are retained, including failures. Mesh-construction failures
    invalidate the profile without aborting the caller's experiment matrix.
    """
    points = np.asarray(path.world_points, dtype=float)
    tangents = estimate_tangents(points, tangent_window_mm)
    areas = np.full(len(points), np.nan)
    reasons = ["invalid_tangent"] * len(points)
    try:
        mesh = mask_to_mesh(volume).mesh
    except (ValueError, RuntimeError) as error:
        reasons = [f"mesh_failed:{type(error).__name__}"] * len(points)
    else:
        for i, (point, tangent) in enumerate(zip(points, tangents, strict=True)):
            areas[i], reasons[i] = section_area(mesh, point, tangent)
    return CrossSectionalProfile(
        np.asarray(path.distance_mm, dtype=float).copy(), points.copy(), tangents,
        areas, equivalent_diameter(areas), np.isfinite(areas), reasons,
    )
