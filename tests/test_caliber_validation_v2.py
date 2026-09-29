"""V2 geometry, tangent, failure-accounting and frozen v1 regression checks."""
import csv
from pathlib import Path

import numpy as np
import pyvista as pv
import pytest

from neurovasc.geometry.cross_sectional_caliber import (
    equivalent_diameter, estimate_tangents, plane_basis, section_area,
    measure_cross_sectional_caliber,
)
from neurovasc.geometry.main_path import MainCenterlinePath
from neurovasc.io.volume import Volume
from scripts.run_caliber_validation import experiment_matrix
from validation.synthetic_caliber.comparison_v2 import prepare, score, matched_metrics, trimmed_indices
from validation.synthetic_caliber.evaluate import evaluate
from validation.synthetic_caliber.experiments_v2 import experiment_matrix_v2, phase_shifted_cylinder
from validation.synthetic_caliber.geometries import straight_cylinder, oblique_cylinder, curved_vessel

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('direction', [(0, 0, 1), (1, 2, 1), (2, 1, 3)])
def test_physical_straight_tangent_and_endpoints(direction):
    direction = np.asarray(direction, dtype=float)
    direction /= np.linalg.norm(direction)
    # Nonuniform physical distances, independent of voxel indices.
    points = np.array([0, 0.4, 1.1, 2.6, 5.1, 9])[:, None] * direction + [12, -8, 30]
    tangent = estimate_tangents(points)
    np.testing.assert_allclose(tangent, np.tile(direction, (len(points), 1)), atol=1e-12)
    np.testing.assert_allclose(estimate_tangents(points[::-1]), -tangent, atol=1e-12)


def test_curved_tangents_follow_analytic_arc():
    theta = np.linspace(-0.625, 0.625, 61)
    points = np.column_stack((24 * (np.cos(theta) - 1), np.zeros(61), 24 * np.sin(theta)))
    truth = np.column_stack((-np.sin(theta), np.zeros(61), np.cos(theta)))
    fitted = estimate_tangents(points, window_mm=6)
    np.testing.assert_allclose(np.linalg.norm(fitted, axis=1), 1)
    angles = np.rad2deg(np.arccos(np.clip(np.sum(fitted * truth, axis=1), -1, 1)))
    assert angles.max() < 5  # One-sided endpoint windows have a small predictable angular bias.


@pytest.mark.parametrize('normal', [(1, 0, 0), (0, 0, 1), (1, 2, 3), (-1, 1, 2)])
def test_plane_basis(normal):
    u, v = plane_basis(normal)
    n = np.asarray(normal) / np.linalg.norm(normal)
    np.testing.assert_allclose(np.array([u, v, n]) @ np.array([u, v, n]).T, np.eye(3), atol=1e-12)
    np.testing.assert_allclose(np.cross(u, v), n, atol=1e-12)


def test_equivalent_diameter_and_invalid_areas():
    np.testing.assert_allclose(equivalent_diameter([np.pi, 4 * np.pi]), [2, 4])
    assert np.isnan(equivalent_diameter([0, -1, np.nan, np.inf])).all()


def test_moderate_binary_cylinder_area():
    context = prepare(straight_cylinder(4))
    profile = measure_cross_sectional_caliber(context.vessel.volume, context.path)
    assert profile.valid.all()
    # Observed area=12.125 vs pi*2^2=12.5664 (3.51% low);
    # limits are 0.6 mm² (~4.8%) and 0.12 mm (~3%), not fit-to-pass broad bounds.
    assert np.max(np.abs(profile.area_mm2 - 4 * np.pi)) < 0.6
    assert np.max(np.abs(profile.diameter_mm - 4)) < 0.12


def test_analytic_surface_circle_area():
    # Isolates the cutter/contour calculation from binary voxelization error.
    mesh = pv.Cylinder(direction=(0, 0, 1), radius=2, height=20, resolution=256).triangulate()
    area, reason = section_area(mesh, (0, 0, 0), (0, 0, 1))
    assert reason == 'ok'
    assert area == pytest.approx(4 * np.pi, rel=0.001)


@pytest.mark.parametrize('spacing', [(0.5, 0.5, 0.5), (0.46875, 0.46875, 0.8)])
def test_oblique_and_anisotropic_transforms(spacing):
    context = prepare(oblique_cylinder(4, spacing=spacing, direction=(1, 2, 1)))
    profile = measure_cross_sectional_caliber(context.vessel.volume, context.path)
    assert profile.valid.all()
    assert np.mean(np.abs(profile.diameter_mm - 4)) < 0.2  # 5% diameter bound.
    np.testing.assert_allclose(np.linalg.norm(profile.tangent_vectors, axis=1), 1)


def test_rotated_affine_equivariance():
    context = prepare(straight_cylinder(4, spacing=(0.46875, 0.46875, 0.8)))
    original = measure_cross_sectional_caliber(context.vessel.volume, context.path)
    angle = 0.7
    rotation = np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
    transform = np.eye(4)
    transform[:3, :3] = rotation
    transform[:3, 3] = [15, -32, 8]
    volume = Volume(context.vessel.volume.data, transform @ context.vessel.volume.affine, context.vessel.volume.spacing, Path('rotated.nii.gz'))
    points = context.path.world_points @ rotation.T + transform[:3, 3]
    path = MainCenterlinePath(context.path.voxel_points, points, context.path.distance_mm)
    rotated = measure_cross_sectional_caliber(volume, path)
    assert rotated.valid.all()
    np.testing.assert_allclose(rotated.diameter_mm, original.diameter_mm, atol=1e-6)


def test_curved_profile_and_trim_samples():
    context = prepare(curved_vessel(3))
    profile = measure_cross_sectional_caliber(context.vessel.volume, context.path)
    assert profile.valid.all()
    indices, distance = trimmed_indices(context, 2)
    assert 0 < len(indices) < len(profile.diameter_mm)
    assert distance[0] == 0
    assert np.isfinite(profile.diameter_mm[indices]).all()


@pytest.mark.parametrize('origin,normal,reason', [
    ((0, 0, 10), (0, 0, 1), 'empty_intersection'),
    ((0, 0, 0), (0, 0, 0), 'invalid_tangent'),
    ((3, 0, 0), (0, 0, 1), 'origin_outside_contour'),
    ((np.nan, 0, 0), (0, 0, 1), 'invalid_origin'),
])
def test_failed_sections(origin, normal, reason):
    area, status = section_area(pv.Sphere(radius=1), origin, normal)
    assert np.isnan(area)
    assert status == reason


def test_open_and_disconnected_contours():
    area, reason = section_area(pv.Plane(direction=(1, 0, 0)), (0, 0, 0), (0, 0, 1))
    assert np.isnan(area) and reason == 'open_or_branched_contour'
    mesh = pv.Sphere(radius=1) + pv.Sphere(radius=1, center=(4, 0, 0))
    area, reason = section_area(mesh, (0, 0, 0), (0, 0, 1))
    assert np.isnan(area) and reason == 'multiple_contours'


def test_degenerate_path_and_empty_volume_remain_explicit():
    vessel = straight_cylinder()
    path = MainCenterlinePath(np.zeros((3, 3), dtype=int), np.zeros((3, 3)), np.zeros(3))
    result = measure_cross_sectional_caliber(vessel.volume, path)
    assert not result.valid.any()
    assert result.failure_reasons == ['invalid_tangent'] * 3
    vessel.volume.data[:] = 0
    result = measure_cross_sectional_caliber(vessel.volume, path)
    assert np.isnan(result.diameter_mm).all()
    assert all(reason.startswith('mesh_failed:') for reason in result.failure_reasons)


def test_matched_comparison_uses_common_support_and_counts_missingness():
    result = matched_metrics([4.2, 9], [4.05, np.nan], [4, 4])
    assert result['winner'] == 'cross_sectional'
    assert result['common_valid_samples'] == 1
    assert result['common_valid_fraction'] == 0.5
    assert result['comparison_status'] == 'partial'
    assert result['cross_sectional_mae_minus_edt_mae'] == pytest.approx(-0.15)
    assert matched_metrics([4.2], [4.205], [4])['winner'] == 'tie'
    assert matched_metrics([4.0], [4.2], [4])['winner'] == 'edt'
    assert matched_metrics([4], [np.nan], [4])['winner'] == 'unscorable'
    stats = score(np.array([4.05, np.nan]), np.array([4, 4]))
    assert stats['sample_count'] == 2 and stats['invalid_sample_count'] == 1
    assert stats['valid_fraction'] == 0.5
    assert score([], [])['valid_fraction'] is None


def test_grid_phase_and_shared_skeleton_failure_are_not_hidden():
    centered = straight_cylinder(4)
    shifted = phase_shifted_cylinder(4, (0.5, 0.5, 0.5), (0.5, 0.5, 0.5))
    np.testing.assert_allclose(shifted.volume.affine[:3, 3] - centered.volume.affine[:3, 3], [0.25] * 3)
    assert shifted.volume.data.any()
    with pytest.raises(ValueError, match='no centerline points'):
        prepare(shifted)
    matrix = list(experiment_matrix_v2())
    assert len(matrix) == 58
    assert sum(e.cohort == 'core' for e in matrix) == 26


def test_all_frozen_v1_edt_results_unchanged():
    with (ROOT / 'validation/public_results/caliber_v1/results.csv').open() as file:
        frozen = list(csv.DictReader(file))
    generated = [evaluate(v, trim)[0] for v in experiment_matrix() for trim in (0.0, 2.0)]
    assert len(generated) == len(frozen) == 52
    for actual, expected in zip(generated, frozen, strict=True):
        for key in ('estimated_mean_mm', 'mae_mm', 'rmse_mm', 'signed_bias_mm', 'percent_error', 'sample_count'):
            assert actual[key] == pytest.approx(float(expected[key]), abs=1e-12, rel=1e-12)
