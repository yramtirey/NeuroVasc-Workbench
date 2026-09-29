"""Deterministic synthetic centerline validation, with failures recorded as data."""
import argparse
import csv
from dataclasses import asdict
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

from neurovasc.geometry.geodesic_centerline import GeodesicParameters
from scripts.run_caliber_validation import write_csv
from validation.synthetic_centerline.geometries import matrix
from validation.synthetic_centerline.evaluate import evaluate, ESTIMATORS
from validation.synthetic_centerline.report import summarize, compare, TIES
from validation.synthetic_centerline.figures import create_figures


def output_directory(path):
    path = Path(path).resolve()
    allowed = (ROOT / 'validation/outputs/centerline_v1').resolve()
    if not path.is_relative_to(allowed):
        raise ValueError('Output must be centerline_v1 or one of its subdirectories; caliber archives are protected')
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'validation/outputs/centerline_v1')
    args = parser.parse_args()
    try:
        output = output_directory(args.output)
    except ValueError as error:
        parser.error(str(error))
    output.mkdir(parents=True, exist_ok=True)
    (output / 'masks').mkdir(exist_ok=True)
    with (ROOT / 'validation/outputs/caliber_v2/results.csv').open() as stream:
        known = {r['experiment_id'] for r in csv.DictReader(stream) if r['status'] == 'upstream_failed'}
    rows, profiles, branches, pairs, geometries, examples = [], [], [], [], [], {}
    for eid, experiment in matrix():
        vessel = experiment.vessel
        meta = {'geometry_id': eid, 'geometry': vessel.geometry, 'cohort': experiment.cohort,
                'spacing_x_mm': vessel.volume.spacing[0], 'spacing_y_mm': vessel.volume.spacing[1], 'spacing_z_mm': vessel.volume.spacing[2],
                'phase_id': experiment.phase_id, 'orientation_id': experiment.orientation_id,
                'diameter_mm': vessel.proximal_diameter_mm if vessel.geometry not in ('tapered', 'bifurcation') else None,
                'proximal_diameter_mm': vessel.proximal_diameter_mm, 'distal_diameter_mm': vessel.distal_diameter_mm,
                'mask_voxel_count': int(np.count_nonzero(vessel.volume.data)), 'known_v2_failure': eid in known}
        geometries.append({**vessel.metadata(), **meta, 'affine': vessel.volume.affine.tolist(), 'shape': list(vessel.volume.shape), 'phase_voxels': experiment.phase_voxels})
        nib.save(nib.Nifti1Image(vessel.volume.data, vessel.volume.affine), output / 'masks' / f'{eid}.nii.gz')
        example = {'vessel': vessel, 'known_failure': eid in known, 'cohort': experiment.cohort}
        pair = []
        for method in ESTIMATORS:
            result, samples, branch_rows, graph, path = evaluate(vessel, method)
            row = {**meta, 'estimator_name': method, **result}
            rows.append(row)
            pair.append(row)
            profiles.extend({**meta, 'estimator_name': method, **s} for s in samples)
            branches.extend({**meta, 'estimator_name': method, **b} for b in branch_rows)
            example[method] = (graph, path)
        pairs.append(compare(*pair))
        examples[eid] = example
        print(f"{eid}: Lee={pair[0]['status']}, alternative={pair[1]['status']}", flush=True)
    if not known.issubset(examples):
        raise RuntimeError('Known caliber-v2 failures missing from validation matrix')
    summary = {'version': 'centerline_v1', 'geometry_count': len(geometries), 'run_count': len(rows),
               'nonbranching_geometry_count': 58, 'branching_geometry_count': 2,
               **summarize(rows, pairs), 'parameters': {**asdict(GeodesicParameters()), 'tangent_window_mm': 6.0,
               'junction_tangent_exclusion_mm': 3.0, 'principal_endpoint_normalized_tolerance': .1,
               'branch_minimum_coverage': .8, 'branch_endpoint_tolerance_max_voxels': 2,
               'comparison_tie_tolerances': TIES, 'trim_mm': 0},
               'counting': 'One untrimmed run per estimator and geometry. Successful-only errors, equal geometry weight; common-geometry support also reported. Failures are never zero error.',
               'dependencies': {p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'scikit-image', 'networkx', 'nibabel', 'matplotlib')},
               'python': platform.python_version()}
    sources = [*ROOT.glob('neurovasc/**/*.py'), *ROOT.glob('validation/synthetic_centerline/*.py'),
               *ROOT.glob('validation/synthetic_caliber/*.py'), ROOT / 'scripts/run_centerline_validation.py',
               ROOT / 'scripts/run_caliber_validation.py', ROOT / 'scripts/run_caliber_validation_v2.py']
    summary['source_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)}
    for filename, content in [('results.csv', rows), ('profiles.csv', profiles), ('branch_results.csv', branches), ('matched_comparison.csv', pairs)]:
        write_csv(output / filename, content)
    for filename, content in [('summary.json', summary), ('geometries.json', geometries)]:
        (output / filename).write_text(json.dumps(content, indent=2, allow_nan=False) + '\n')
    create_figures(output, rows, profiles, examples)
    # This run's own artifacts only; isolated regression runs are not attributed
    # to this runner, even if they exist in a nested directory on a later rerun.
    artifacts = [p for p in output.iterdir() if p.is_file() and p.name != 'manifest.txt']
    artifacts += list((output / 'figures').glob('*.png')) + list((output / 'masks').glob('*.nii.gz'))
    (output / 'manifest.txt').write_text('\n'.join(sorted([str(p.relative_to(output)) for p in artifacts] + ['manifest.txt'])) + '\n')
    print(json.dumps({key: summary[key] for key in ('geometry_count', 'run_count', 'lee', 'alternative', 'comparison')}, indent=2))


if __name__ == '__main__':
    main()
