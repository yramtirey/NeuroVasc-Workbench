"""Coordinate-based finite analytic centerlines; used only for evaluation."""
import numpy as np


def segment_projection(points, segment):
    start, end = np.asarray(segment, dtype=float)
    direction = end - start
    length = np.linalg.norm(direction)
    direction /= length
    position = np.clip((points - start) @ direction / length, 0, 1)
    return start + position[:, None] * (end - start), np.tile(direction, (len(points), 1)), position


def project(vessel, points, branch_name=None):
    points = np.asarray(points, dtype=float)
    if vessel.geometry == 'bifurcation':
        names = list(branch_name) if isinstance(branch_name, (list, tuple)) else [branch_name] if branch_name else list(vessel.branches)
        candidates = [segment_projection(points, vessel.branches[name]) for name in names]
        errors = np.stack([np.linalg.norm(points - c[0], axis=1) for c in candidates])
        chosen = np.argmin(errors, axis=0)
        mapped = np.array([candidates[k][0][i] for i, k in enumerate(chosen)])
        tangents = np.array([candidates[k][1][i] for i, k in enumerate(chosen)])
        position = np.array([candidates[k][2][i] for i, k in enumerate(chosen)])
        labels = [names[k] for k in chosen]
    elif vessel.bend_radius_mm is not None:
        radius = vessel.bend_radius_mm
        theta = np.clip(np.arctan2(points[:, 2], points[:, 0] + radius), -vessel.length_mm / (2 * radius), vessel.length_mm / (2 * radius))
        mapped = np.column_stack((radius * (np.cos(theta) - 1), np.zeros(len(theta)), radius * np.sin(theta)))
        tangents = np.column_stack((-np.sin(theta), np.zeros(len(theta)), np.cos(theta)))
        position = theta * radius / vessel.length_mm + 0.5
        labels = ['main'] * len(points)
    else:
        mapped, tangents, position = segment_projection(points, np.array([-0.5, 0.5])[:, None] * vessel.length_mm * vessel.direction)
        labels = ['main'] * len(points)
    return mapped, tangents, position, labels


def sample_truth(vessel, branch_name=None, samples=200):
    if vessel.geometry == 'bifurcation':
        segment = vessel.branches[branch_name]
        return segment[0] + np.linspace(0, 1, samples)[:, None] * (segment[1] - segment[0])
    if vessel.bend_radius_mm is not None:
        r = vessel.bend_radius_mm
        theta = np.linspace(-vessel.length_mm / (2 * r), vessel.length_mm / (2 * r), samples)
        return np.column_stack((r * (np.cos(theta) - 1), np.zeros(samples), r * np.sin(theta)))
    return np.linspace(-0.5, 0.5, samples)[:, None] * vessel.length_mm * vessel.direction
