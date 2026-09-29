"""Equal-geometry summaries: failed runs never become zero-error runs."""
import numpy as np

METRICS = {'coverage': 'coverage_fraction', 'position_error_mm': 'mean_position_error_mm',
           'length_error_mm': 'length_error_mm', 'tangent_error_deg': 'mean_tangent_error_deg'}
TIES = {'coverage': .01, 'position_error_mm': .01, 'length_error_mm': .1, 'tangent_error_deg': 1.0}


def aggregate(rows):
    valid = [r for r in rows if r['extraction_success']]
    result = {'run_count': len(rows), 'success_count': len(valid), 'success_rate': len(valid) / len(rows) if rows else None,
              'failure_inclusive_mean_coverage': sum(r.get('coverage_fraction', 0) for r in valid) / len(rows) if rows else None,
              'correct_principal_count': sum(r['correct_principal_path'] for r in rows),
              'topology_correct_count': sum(r.get('topology_correct', False) for r in rows)}
    for key, column in METRICS.items():
        values = [r[column] for r in valid if r.get(column) is not None]
        result[f'mean_{key}'] = float(np.mean(values)) if values else None
    return result


def compare(a, b):
    out = {'geometry_id': a['geometry_id'], 'geometry': a['geometry'], 'phase_id': a['phase_id'],
           'orientation_id': a['orientation_id'], 'spacing_z_mm': a['spacing_z_mm'],
           'lee_success': a['extraction_success'], 'alternative_success': b['extraction_success'],
           'known_v2_failure': a['known_v2_failure']}
    out['outcome'] = ('both_succeeded' if a['extraction_success'] and b['extraction_success'] else
                      'lee_failed_alternative_succeeded' if b['extraction_success'] else
                      'alternative_failed_lee_succeeded' if a['extraction_success'] else 'both_failed')
    for metric, column in METRICS.items():
        delta = b[column] - a[column] if out['outcome'] == 'both_succeeded' and a.get(column) is not None and b.get(column) is not None else None
        out[f'{metric}_alternative_minus_lee'] = delta
        signed = -delta if metric == 'coverage' and delta is not None else delta
        out[f'{metric}_winner'] = None if delta is None else 'tie' if abs(delta) <= TIES[metric] else 'alternative' if signed < 0 else 'lee'
    return out


def summarize(rows, pairs):
    def methods(selected):
        return {name: aggregate([r for r in selected if r['estimator_name'] == name]) for name in ('lee', 'alternative')}
    result = methods(rows)
    result['comparison'] = {key: sum(p['outcome'] == key for p in pairs) for key in ('both_succeeded', 'both_failed', 'lee_failed_alternative_succeeded', 'alternative_failed_lee_succeeded')}
    for metric in METRICS:
        result['comparison'][f'{metric}_wins'] = {name: sum(p[f'{metric}_winner'] == name for p in pairs) for name in ('lee', 'alternative', 'tie')}
    shared = {p['geometry_id'] for p in pairs if p['outcome'] == 'both_succeeded'}
    result['common_geometry_support'] = methods([r for r in rows if r['geometry_id'] in shared])
    for field in ('geometry', 'phase_id', 'orientation_id', 'spacing_z_mm', 'cohort'):
        result[f'by_{field}'] = {str(value): methods([r for r in rows if r[field] == value]) for value in sorted(set(r[field] for r in rows))}
    result['known_v2_failures'] = [p for p in pairs if p['known_v2_failure']]
    return result
