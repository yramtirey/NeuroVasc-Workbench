"""Independent truth/metric checks plus experimental centerline regression guards."""
from copy import deepcopy
from pathlib import Path

import nibabel as nib
import networkx as nx
import numpy as np
import pytest

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.main_path import extract_main_path
from neurovasc.geometry.geodesic_centerline import extract_geodesic_tree, GeodesicParameters
from neurovasc.graph.vascular_graph import build_vascular_graph
from validation.synthetic_caliber.geometries import straight_cylinder, oblique_cylinder, curved_vessel, tapered_vessel
from validation.synthetic_caliber.experiments_v2 import phase_shifted_cylinder
from validation.synthetic_centerline.geometries import y_bifurcation, matrix
from validation.synthetic_centerline.truth import project, sample_truth
from validation.synthetic_centerline.evaluate import lee, evaluate, score_points, angular_errors, validate_path, topology
from validation.synthetic_centerline.report import compare, aggregate
from scripts.run_centerline_validation import output_directory


@pytest.mark.parametrize('vessel', [straight_cylinder(3), oblique_cylinder(3, direction=(1, 2, 1)), curved_vessel(3), tapered_vessel(4, 2)])
def test_analytic_truth_and_reversal(vessel):
    points = sample_truth(vessel, samples=101)
    mapped, tangent, u, _ = project(vessel, points)
    np.testing.assert_allclose(mapped, points, atol=1e-12)
    np.testing.assert_allclose(u, np.linspace(0, 1, 101), atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(tangent, axis=1), 1)
    metrics, _ = score_points(vessel, points[::-1])
    assert metrics['coverage_fraction'] == pytest.approx(1)
    assert metrics['mean_position_error_mm'] < 1e-12
    assert metrics['monotonic_order']
    assert metrics['normalized_start'] == pytest.approx(1)
    if vessel.geometry == 'curved':
        # One-sided 3 mm half-window fits of a radius-24 mm arc can differ
        # from the endpoint tangent by at most ~3.6 degrees.
        assert metrics['max_tangent_error_deg'] < 4
    else:
        assert metrics['mean_tangent_error_deg'] < 1e-5


def test_position_coverage_and_short_path():
    vessel = straight_cylinder(3)
    points = np.array([[.3, .4, -3], [.3, .4, 0], [.3, .4, 3]])
    metrics, _ = score_points(vessel, points)
    assert metrics['mean_position_error_mm'] == pytest.approx(.5)
    assert metrics['position_rmse_mm'] == pytest.approx(.5)
    assert metrics['coverage_fraction'] == pytest.approx(.2)
    assert metrics['length_error_mm'] == pytest.approx(24)
    # A finite centerline includes end-cap distance, not an infinite-axis shortcut.
    mapped, _, u, _ = project(vessel, np.array([[0, 0, 17]]))
    np.testing.assert_allclose(mapped, [[0, 0, 15]])
    assert u[0] == 1


def test_angular_metric_sign_and_invalid_vectors():
    np.testing.assert_allclose(angular_errors([[0, 0, 1], [0, 0, -1], [1, 0, 0]], [[0, 0, 1]] * 3), [0, 0, 90])
    assert np.isnan(angular_errors([[0, 0, 0]], [[0, 0, 1]])[0])


@pytest.mark.parametrize('spacing', [(.5, .5, .5), (.46875, .46875, .8)])
def test_y_truth_and_recovered_topology(spacing):
    vessel = y_bifurcation(spacing)
    assert vessel.length_mm == pytest.approx(40)
    assert sum(np.linalg.norm(s[1]) for s in vessel.branches.values()) == pytest.approx(58)
    for name in vessel.branches:
        points = sample_truth(vessel, name)
        mapped, _, position, _ = project(vessel, points, name)
        np.testing.assert_allclose(mapped, points, atol=1e-12)
        assert np.ptp(position) == pytest.approx(1)
    for method in ('lee', 'alternative'):
        result, _, branches, _, _ = evaluate(vessel, method)
        assert result['extraction_success'] and result['topology_correct']
        assert 0 <= result['coverage_fraction'] <= 1
        assert result['endpoint_count'] == 3 and result['junction_count'] == 1 and result['cycle_count'] == 0
        # A voxel-snapped junction may lie one spatial voxel diagonal away.
        # Observed worst error is sqrt(.5²+.5²)=0.7071 mm.
        assert result['junction_error_mm'] <= np.linalg.norm(spacing) + 1e-12
        assert all(b['branch_path_found'] for b in branches)


def test_baseline_is_unchanged_adapter():
    vessel = straight_cylinder(4)
    direct = extract_main_path(build_vascular_graph(extract_centerline(vessel.volume)))
    wrapped = extract_main_path(lee(vessel.volume))
    np.testing.assert_array_equal(wrapped.voxel_points, direct.voxel_points)
    np.testing.assert_array_equal(wrapped.world_points, direct.world_points)


@pytest.mark.parametrize('spacing', [(.5, .5, .5), (.46875, .46875, .8)])
@pytest.mark.parametrize('diameter', [1, 2, 3, 4, 5])
def test_all_half_phase_failures_recovered(spacing, diameter):
    vessel = phase_shifted_cylinder(diameter, spacing, (.5, .5, .5))
    baseline = evaluate(vessel, 'lee')[0]
    alternative = evaluate(vessel, 'alternative')[0]
    assert not baseline['extraction_success']
    assert alternative['extraction_success'] and alternative['correct_principal_path']
    assert alternative['coverage_fraction'] >= .95


def test_thin_half_phase_regression():
    vessel = phase_shifted_cylinder(1, (.5, .5, .5), (.5, .5, .5))
    row = evaluate(vessel, 'alternative')[0]
    # The four nearest axes are sqrt(.25²+.25²) mm from analytic truth.
    assert row['mean_position_error_mm'] == pytest.approx(np.sqrt(.125))
    assert row['coverage_fraction'] == pytest.approx(29.5 / 30)
    assert row['length_error_mm'] < .5  # Less than one voxel, observed 0.134 mm.


def test_anisotropic_coordinates_and_rigid_affine_rotation():
    vessel = oblique_cylinder(3, spacing=(.46875, .46875, .8), direction=(1, 2, 1))
    original = extract_geodesic_tree(vessel.volume)
    angle = .71
    transform = np.eye(4)
    transform[:3, :3] = [[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]]
    transform[:3, 3] = [12, -8, 3]
    rotated = deepcopy(vessel.volume)
    rotated.affine = transform @ rotated.affine
    changed = extract_geodesic_tree(rotated)
    assert set(original.graph.nodes) == set(changed.graph.nodes)
    assert {frozenset(e) for e in original.graph.edges} == {frozenset(e) for e in changed.graph.edges}
    for node in original.graph:
        np.testing.assert_allclose(original.graph.nodes[node]['world'], nib.affines.apply_affine(vessel.volume.affine, node), atol=1e-12)
        np.testing.assert_allclose(changed.graph.nodes[node]['world'], nib.affines.apply_affine(transform, original.graph.nodes[node]['world']), atol=1e-12)


def test_empty_disconnected_shear_and_limits_fail_explicitly():
    volume = straight_cylinder(3).volume
    empty = deepcopy(volume)
    empty.data[:] = 0
    with pytest.raises(ValueError, match='empty'):
        extract_geodesic_tree(empty)
    disconnected = deepcopy(empty)
    disconnected.data[1, 1, 1] = disconnected.data[-2, -2, -2] = 1
    with pytest.raises(ValueError, match='disconnected'):
        extract_geodesic_tree(disconnected)
    sheared = deepcopy(volume)
    sheared.affine[0, 1] = .1
    with pytest.raises(ValueError, match='affine'):
        extract_geodesic_tree(sheared)
    with pytest.raises(ValueError, match='size_limit'):
        extract_geodesic_tree(volume, GeodesicParameters(maximum_foreground_voxels=10))
    with pytest.raises(ValueError, match='invalid_geodesic_limits'):
        extract_geodesic_tree(volume, GeodesicParameters(maximum_branches=0))


@pytest.mark.parametrize('corruption', ['distance', 'nonfinite', 'edge', 'nodes'])
def test_malformed_graph_path_rejected(corruption):
    graph = lee(straight_cylinder(3).volume)
    path = extract_main_path(graph)
    if corruption == 'distance':
        path.distance_mm[1] = -1
    elif corruption == 'nonfinite':
        path.world_points[1, 0] = np.nan
    elif corruption == 'edge':
        a, b = next(iter(graph.graph.edges))
        graph.graph[a][b]['length_mm'] = -1
    else:
        path.voxel_points[1] = path.voxel_points[0]
    with pytest.raises(ValueError):
        validate_path(graph.graph, path)


def test_cycle_and_logical_junction_accounting():
    graph = nx.cycle_graph(4)
    stats, groups = topology(graph)
    assert stats['cycle_count'] == 1 and stats['endpoint_count'] == 0
    graph.add_edge(0, 4)
    graph.add_edge(1, 5)
    stats, groups = topology(graph)
    assert stats['junction_voxel_count'] == 2 and stats['junction_count'] == 1


def test_matrix_ids_mask_identity_and_failure_reporting():
    import csv
    from validation.synthetic_caliber.experiments_v2 import experiment_matrix_v2
    root = Path(__file__).resolve().parents[1]
    experiments = list(matrix())
    assert len(experiments) == 60
    for (_, current), old in zip(experiments, experiment_matrix_v2()):
        np.testing.assert_array_equal(current.vessel.volume.data, old.vessel.volume.data)
        np.testing.assert_array_equal(current.vessel.volume.affine, old.vessel.volume.affine)
    with (root / 'validation/public_results/caliber_v2/results.csv').open() as stream:
        known = {r['experiment_id'] for r in csv.DictReader(stream) if r['status'] == 'upstream_failed'}
    assert len(known) == 12
    assert known <= {eid for eid, _ in experiments}
    for eid, experiment in experiments:
        if eid in known:
            assert not evaluate(experiment.vessel, 'lee')[0]['extraction_success']
            assert evaluate(experiment.vessel, 'alternative')[0]['extraction_success']


def test_outputs_cannot_overwrite_archives():
    root = Path(__file__).resolve().parents[1]
    for path in ['caliber_v1', 'caliber_v1/nested', 'caliber_v2', '.', 'centerline_v1/../../outside']:
        with pytest.raises(ValueError):
            output_directory(root / 'validation/outputs' / path)
    assert output_directory(root / 'validation/outputs/centerline_v1')


def test_deterministic_mask_inferred_endpoints_and_comparison():
    vessel = straight_cylinder(3)
    first = extract_geodesic_tree(vessel.volume)
    second = extract_geodesic_tree(vessel.volume)
    assert list(first.graph.edges) == list(second.graph.edges)
    endpoints = [n for n, degree in first.graph.degree if degree == 1]
    assert len(endpoints) == 2
    coordinates = np.array([first.graph.nodes[n]['world'] for n in endpoints])
    assert np.ptp(coordinates[:, 2]) >= vessel.length_mm - 1
    # A failed method has no metric win and no invented zero positional error.
    row = evaluate(vessel, 'lee')[0]
    row.update(geometry_id='test', geometry='straight', phase_id='centered', orientation_id='z', spacing_z_mm=.5, known_v2_failure=False)
    failed = dict(row, extraction_success=False, correct_principal_path=False)
    pair = compare(failed, row)
    assert pair['outcome'] == 'lee_failed_alternative_succeeded'
    assert pair['position_error_mm_winner'] is None
    summary = aggregate([row, failed])
    assert summary['success_rate'] == .5
    assert summary['mean_position_error_mm'] == row['mean_position_error_mm']
    assert summary['failure_inclusive_mean_coverage'] == row['coverage_fraction'] / 2
