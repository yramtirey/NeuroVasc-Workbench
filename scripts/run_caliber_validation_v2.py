"""Reproducible matched mesh-section vs EDT validation; v1 outputs are read-only."""
import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'validation/outputs/.matplotlib'))

import matplotlib
matplotlib.use('Agg')
import nibabel as nib
import numpy as np

from scripts.run_caliber_validation import write_csv
from validation.synthetic_caliber.experiments_v2 import experiment_matrix_v2
from validation.synthetic_caliber.comparison_v2 import ESTIMATORS, Estimates, prepare, trimmed_indices, score, matched_metrics
from validation.synthetic_caliber.figures_v2 import create_figures


def aggregate(rows):
    valid = [r for r in rows if r.get('valid_sample_count', 0) > 0]
    def mean(key, selected=valid):
        return float(np.mean([r[key] for r in selected])) if selected else None
    total = sum(r['sample_count'] for r in rows)
    return {
        'overall_mae': mean('mae_mm'),
        'overall_rmse': float(np.sqrt(np.mean([r['rmse_mm'] ** 2 for r in valid]))) if valid else None,
        'bias': mean('signed_bias_mm'),
        'isotropic_mae': mean('mae_mm', [r for r in valid if r['spacing_z_mm'] == 0.5]),
        'anisotropic_mae': mean('mae_mm', [r for r in valid if r['spacing_z_mm'] == 0.8]),
        'valid_fraction': sum(r['valid_sample_count'] for r in rows) / total if total else None,
        'attempted_samples': total, 'valid_samples': sum(r['valid_sample_count'] for r in rows),
        'invalid_samples': sum(r['invalid_sample_count'] for r in rows),
        'runs': len(rows), 'scorable_runs': len(valid),
        'upstream_failed_runs': sum(r['status'] == 'upstream_failed' for r in rows),
    }


def clean_json(value):
    if isinstance(value, dict):
        return {k: clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(v) for v in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'validation/outputs/caliber_v2')
    args = parser.parse_args()
    output = args.output.resolve()
    v1 = ROOT / 'validation/outputs/caliber_v1'
    # In particular, prevent using an ancestor of the frozen v1 output directory.
    if not output.is_relative_to(ROOT / 'validation/outputs') or output == ROOT / 'validation/outputs' or output.is_relative_to(v1):
        parser.error('Choose an output subdirectory under validation/outputs other than caliber_v1')
    output.mkdir(parents=True, exist_ok=True)
    (output / 'masks').mkdir(exist_ok=True)
    rows, pairs, samples, geometries = [], [], [], []
    for index, experiment in enumerate(experiment_matrix_v2(), 1):
        vessel = experiment.vessel
        eid = f'{index:02d}_{vessel.geometry}_{vessel.proximal_diameter_mm:g}_{experiment.phase_id}_{experiment.orientation_id}_z{vessel.volume.spacing[2]:g}'
        metadata = {
            'experiment_id': eid, 'cohort': experiment.cohort, 'geometry': vessel.geometry,
            'spacing_x_mm': vessel.volume.spacing[0], 'spacing_y_mm': vessel.volume.spacing[1], 'spacing_z_mm': vessel.volume.spacing[2],
            'phase_id': experiment.phase_id, 'orientation_id': experiment.orientation_id,
            'true_diameter_mm': vessel.proximal_diameter_mm if vessel.geometry != 'tapered' else None,
            'proximal_diameter_mm': vessel.proximal_diameter_mm, 'distal_diameter_mm': vessel.distal_diameter_mm,
        }
        geometries.append({**vessel.metadata(), **metadata, 'grid_phase': experiment.phase_id,
                           'phase_voxels': experiment.phase_voxels, 'affine': vessel.volume.affine.tolist(), 'shape': vessel.volume.shape})
        nib.save(nib.Nifti1Image(vessel.volume.data, vessel.volume.affine), output / 'masks' / f'{eid}.nii.gz')
        try:
            context = prepare(vessel)
        except (ValueError, RuntimeError) as error:
            for trim in (0.0, 2.0):
                base = {**metadata, 'trim_mm': trim}
                for name in ESTIMATORS:
                    rows.append({**base, 'estimator_name': name, **score(np.array([]), np.array([])),
                                 'status': 'upstream_failed', 'failure_reason': str(error), 'trim_applied': None, 'samples_removed': None})
                pairs.append({**base, **matched_metrics([], [], []), 'failure_reason': str(error)})
            print(f'{eid}: shared centerline failure: {error}', flush=True)
            continue
        estimates = {}
        for name, adapter in ESTIMATORS.items():
            try:
                estimates[name] = adapter(context)
            except (ValueError, RuntimeError) as error:
                n = len(context.path.world_points)
                estimates[name] = Estimates(np.full(n, np.nan), np.full(n, np.nan), np.full((n, 3), np.nan), [f'adapter_failed:{error}'] * n)
        for trim in (0.0, 2.0):
            indices, distances = trimmed_indices(context, trim)
            points = context.path.world_points[indices]
            truth = vessel.true_diameter(points)
            positions = vessel.normalized_position(points)
            base = {**metadata, 'trim_mm': trim}
            for name, estimates_for_method in estimates.items():
                values = estimates_for_method.diameter_mm[indices]
                reasons = [estimates_for_method.reasons[i] for i in indices]
                row = {
                    **base, 'estimator_name': name, **score(values, truth),
                    'failure_reason': ';'.join(sorted(set(r for r in reasons if r != 'ok'))),
                    'trim_applied': len(indices) < len(context.path.world_points),
                    'samples_removed': len(context.path.world_points) - len(indices),
                }
                rows.append(row)
                for j, i in enumerate(indices):
                    valid = bool(np.isfinite(values[j]) and values[j] > 0)
                    samples.append({
                        **base, 'estimator_name': name, 'sample_index': int(i), 'distance_mm': float(distances[j]),
                        'normalized_position': float(positions[j]), 'true_diameter_mm': float(truth[j]),
                        'diameter_mm': float(values[j]) if valid else None,
                        'area_mm2': float(estimates_for_method.area_mm2[i]) if np.isfinite(estimates_for_method.area_mm2[i]) else None,
                        'error_mm': float(values[j] - truth[j]) if valid else None,
                        **dict(zip(('x_mm', 'y_mm', 'z_mm'), map(float, points[j]))),
                        **{axis: float(t) if np.isfinite(t) else None for axis, t in zip(('tangent_x', 'tangent_y', 'tangent_z'), estimates_for_method.tangents[i])},
                        'valid': valid, 'status': reasons[j],
                    })
            pairs.append({**base, **matched_metrics(estimates['edt'].diameter_mm[indices], estimates['cross_sectional'].diameter_mm[indices], truth)})
        print(f'{eid}: evaluated both estimators', flush=True)

    scored = [p for p in pairs if p['winner'] != 'unscorable']
    cross_rows = [r for r in rows if r['estimator_name'] == 'cross_sectional']
    cross_valid = [r for r in cross_rows if r['valid_sample_count'] > 0]
    winners = Counter(p['winner'] for p in pairs)
    by_geometry = {g: float(np.mean([p['cross_sectional_mae_minus_edt_mae'] for p in scored if p['geometry'] == g])) for g in sorted(set(p['geometry'] for p in scored))}
    def compare_subset(subset):
        return {
            'edt': aggregate([r for r in subset if r['estimator_name'] == 'edt']),
            'cross_sectional': aggregate([r for r in subset if r['estimator_name'] == 'cross_sectional']),
        }
    summary = {
        'version': 'caliber_v2', 'geometry_count': len(geometries), 'experiment_count': len(pairs),
        'estimator_run_count': len(rows), 'matched_experiment_count': len(scored),
        'shared_centerline_failure_geometry_count': len({r['experiment_id'] for r in rows if r['status'] == 'upstream_failed'}),
        **compare_subset(rows),
        'core_v1_phantoms': compare_subset([r for r in rows if r['cohort'] == 'core']),
        'by_trim': {str(t): compare_subset([r for r in rows if r['trim_mm'] == t]) for t in (0.0, 2.0)},
        'by_geometry': {g: compare_subset([r for r in rows if r['geometry'] == g]) for g in ('straight', 'oblique', 'curved', 'tapered')},
        'by_geometry_core': {g: compare_subset([r for r in rows if r['geometry'] == g and r['cohort'] == 'core']) for g in ('straight', 'oblique', 'curved', 'tapered')},
        'comparison': {
            'mae_delta': float(np.mean([p['cross_sectional_mae_minus_edt_mae'] for p in scored])) if scored else None,
            'rmse_delta': float(np.sqrt(np.mean([p['cross_sectional_rmse'] ** 2 for p in scored])) - np.sqrt(np.mean([p['edt_rmse'] ** 2 for p in scored]))) if scored else None,
            'edt_wins': winners['edt'], 'cross_sectional_wins': winners['cross_sectional'], 'ties': winners['tie'],
            'unscorable_pairs': winners['unscorable'], 'tie_tolerance_mm': 0.01,
            'best_geometry_for_cross_sectional': min(by_geometry, key=by_geometry.get) if by_geometry else None,
            'worst_geometry_for_cross_sectional': max(by_geometry, key=by_geometry.get) if by_geometry else None,
            'mean_mae_delta_by_geometry': by_geometry,
            'worst_cross_sectional_percent_error': max((r['percent_error'] for r in cross_valid), default=None),
            'worst_cross_sectional_run': max(cross_valid, key=lambda r: r['percent_error']) if cross_valid else None,
        },
        'cross_section_failure_reasons': dict(Counter(s['status'] for s in samples if s['estimator_name'] == 'cross_sectional' and not s['valid'])),
        'parameters': {'tangent_window_mm': 6.0, 'contour_merge_tolerance_mm': 1e-6, 'trim_mm': [0, 2], 'partial_volume_experiment': False},
        'runtime': {'python': platform.python_version(), 'packages': {p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'scikit-image', 'networkx', 'nibabel', 'matplotlib', 'pyvista', 'vtk')}},
        'source_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for pattern in ('neurovasc/geometry/*.py', 'neurovasc/graph/*.py', 'neurovasc/io/*.py', 'validation/synthetic_caliber/*.py', 'scripts/run_caliber_validation*.py') for path in sorted(ROOT.glob(pattern))},
        'notes': [
            'Geometric/synthetic validation only. Experimental method is not connected to real-case analysis.',
            'Overall metrics are equal-run MAE and sqrt(mean run MSE); both trim variants included as paired, nonindependent runs.',
            'Per-estimator metrics use its valid samples; matched deltas/wins use common valid sample support and exclude unscorable pairs.',
            'Valid fraction is sample-weighted over recovered centerlines, including both trim variants. Shared skeleton failures have zero attempted samples and are counted separately as upstream_failed_runs.',
            'Raw estimators are compared without post-hoc diameter smoothing. Tangents are fitted on the full path before applying production trimming.',
            'Phase offsets shift grid voxel centers in physical space while analytic vessels stay fixed; fractional shifts apply along all three axes.',
            'No optional partial-volume experiment; no calibration to truth and no silent EDT fallback for failed sections.',
        ],
    }
    write_csv(output / 'results.csv', rows)
    write_csv(output / 'matched_comparison.csv', pairs)
    write_csv(output / 'profiles.csv', samples)
    for filename, value in (('summary.json', summary), ('geometries.json', geometries)):
        (output / filename).write_text(json.dumps(clean_json(value), indent=2, allow_nan=False) + '\n')
    create_figures(output, rows, pairs, samples)
    (output / 'manifest.txt').write_text('\n'.join(sorted({str(p.relative_to(output)) for p in output.rglob('*') if p.is_file()} | {'manifest.txt'})) + '\n')
    print(json.dumps({key: summary[key] for key in ('geometry_count', 'experiment_count', 'matched_experiment_count', 'edt', 'cross_sectional', 'comparison')}, indent=2))


if __name__ == '__main__':
    main()
