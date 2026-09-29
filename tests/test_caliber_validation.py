"""Regression bounds characterize this phantom, not clinical acceptability."""
import nibabel as nib
import numpy as np
import pytest

from neurovasc.geometry.diameter_profile import DiameterProfile
from scripts.run_caliber_validation import experiment_matrix
from validation.synthetic_caliber.evaluate import evaluate
from validation.synthetic_caliber.geometries import (
    straight_cylinder, oblique_cylinder, tapered_vessel, curved_vessel,
)
from validation.synthetic_caliber.metrics import comparison_metrics


def test_straight_physical_dimensions_and_metadata(tmp_path):
    vessel = straight_cylinder(4, 30)
    volume = vessel.volume
    points = nib.affines.apply_affine(volume.affine, np.argwhere(volume.data))
    np.testing.assert_allclose(points.min(axis=0), [-2, -2, -15])
    np.testing.assert_allclose(points.max(axis=0), [2, 2, 15])
    np.testing.assert_allclose(points.mean(axis=0), 0, atol=1e-12)
    assert set(np.unique(volume.data)) == {0, 1}
    for axis in range(3):
        assert not np.take(volume.data, [0, -1], axis=axis).any()
    meta = vessel.metadata()
    assert meta['length_mm'] == 30
    assert meta['proximal_diameter_mm'] == meta['distal_diameter_mm'] == 4
    assert meta['spacing_mm'] == [0.5, 0.5, 0.5]
    path = tmp_path / 'phantom.nii.gz'
    nib.save(nib.Nifti1Image(volume.data, volume.affine), path)
    restored = nib.load(path)
    np.testing.assert_array_equal(np.asarray(restored.dataobj), volume.data)
    np.testing.assert_allclose(restored.affine, volume.affine)


def test_analytic_taper_truth_is_independent_of_order_and_trimming():
    vessel = tapered_vessel(4, 2)
    points = np.array([[0, 0, 15], [0, 0, -15], [0, 0, 0], [0, 0, 7.5]])
    np.testing.assert_allclose(vessel.normalized_position(points), [1, 0, 0.5, 0.75])
    np.testing.assert_allclose(vessel.true_diameter(points), [2, 4, 3, 2.5])
    # A reverse, already-trimmed future estimator must still receive correct truth.
    def known_estimator(vessel, trim):
        return DiameterProfile(np.arange(4.0), vessel.true_diameter(points), points)
    row, _ = evaluate(vessel, estimator=known_estimator, estimator_name='known_truth')
    assert row['mae_mm'] == 0
    assert row['correlation'] == pytest.approx(1)


def test_oblique_and_curved_analytic_coordinates():
    oblique = oblique_cylinder()
    points = np.array([-15, 0, 15])[:, None] * np.array([1, 1, 1]) / np.sqrt(3)
    np.testing.assert_allclose(oblique.normalized_position(points), [0, 0.5, 1], atol=1e-12)
    curve = curved_vessel()
    theta = np.array([-15 / 24, 0, 15 / 24])
    points = np.column_stack((24 * (np.cos(theta) - 1), np.zeros(3), 24 * np.sin(theta)))
    np.testing.assert_allclose(curve.normalized_position(points), [0, 0.5, 1], atol=1e-12)
    np.testing.assert_allclose(curve.true_diameter(points), 3)


@pytest.mark.parametrize('vessel', list(experiment_matrix()), ids=lambda v: f'{v.geometry}-{v.proximal_diameter_mm}-{v.volume.spacing[2]}')
def test_full_matrix_pipeline_and_trim(vessel):
    untrimmed, full = evaluate(vessel, 0)
    trimmed, kept = evaluate(vessel, 2)
    assert untrimmed['sample_count'] >= trimmed['sample_count'] >= 5
    assert trimmed['sample_count'] < untrimmed['sample_count']  # All chosen phantoms permit the default trim.
    assert np.isfinite(trimmed['rmse_mm'])
    full_points = {(s['x_mm'], s['y_mm'], s['z_mm']): s for s in full}
    for sample in kept:
        original = full_points[sample['x_mm'], sample['y_mm'], sample['z_mm']]
        assert sample['true_diameter_mm'] == original['true_diameter_mm']
        assert sample['raw_diameter_mm'] == original['raw_diameter_mm']
    assert kept[0]['distance_mm'] == 0


def test_moderate_cylinder_regression():
    row, _ = evaluate(straight_cylinder(4), 2)
    # Observed MAE 0.1231 mm. Bound 0.25 mm is half one voxel, gives
    # numerical/version headroom, and catches a spacing/radius-vs-diameter bug.
    assert row['mae_mm'] < 0.25
    assert abs(row['estimated_mean_mm'] - 4) < 0.25
    assert row['sample_count'] >= 40
    assert row['normalized_position_max'] - row['normalized_position_min'] > 0.7


def test_metrics_distinguish_bias_and_absolute_errors():
    metrics = comparison_metrics([1, 3], [2, 2])
    assert metrics['absolute_error_mm'] == 0
    assert metrics['mae_mm'] == metrics['rmse_mm'] == 1
    assert metrics['percent_error'] == 50
    assert metrics['correlation'] is None
    with pytest.raises(ValueError):
        comparison_metrics([], [])


@pytest.mark.parametrize('kwargs', [{'diameter_mm': -1}, {'spacing': (0, 0.5, 0.5)}, {'spacing': (0.5, 0.5)}, {'length_mm': float('nan')}])
def test_invalid_geometries(kwargs):
    with pytest.raises(ValueError):
        straight_cylinder(**kwargs)
