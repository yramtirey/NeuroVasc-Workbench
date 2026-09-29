"""Experimental selection from segmentation only. Never imported by production."""
from .centerline import extract_centerline
from .main_path import extract_main_path
from .geodesic_centerline import extract_geodesic_tree
from .centerline_qc import Recovery, quality_control
from .centerline_refinement import Refiner, RefinementParameters
from neurovasc.graph.vascular_graph import build_vascular_graph

CANDIDATE = RefinementParameters()


def recover(volume, method):
    try:
        graph = build_vascular_graph(extract_centerline(volume)) if method == 'lee' else extract_geodesic_tree(volume)
        path = extract_main_path(graph).world_points
        return Recovery(graph.graph, path, supports_cycles=method == 'lee')
    except ValueError as error:
        return Recovery(error=str(error), supports_cycles=method == 'lee')


def hybrid(volume, parameters=CANDIDATE, mode='network', lee_result=None, raw_result=None, refined_result=None):
    baseline = lee_result if lee_result is not None else recover(volume, 'lee')
    lee_qc = quality_control(volume, baseline, mode)
    diagnostics = {'selected_method': 'lee', 'lee_extraction_success': baseline.path is not None,
                   'lee_qc_pass': lee_qc['qc_pass'], 'lee_qc_fail_reasons': lee_qc['qc_fail_reasons'],
                   'production_eligible': False, 'fallback_supports_cycles': False, 'fallback_used': False, 'fallback_reason': '', 'geodesic_success': None, 'refinement_success': None}
    if lee_qc['qc_pass']:
        return baseline, {**diagnostics, 'final_qc_pass': True}, lee_qc
    raw = raw_result if raw_result is not None else recover(volume, 'geodesic_raw')
    diagnostics.update(selected_method='geodesic_refined', fallback_used=True,
                       fallback_reason=baseline.error or lee_qc['qc_fail_reasons'], geodesic_success=raw.path is not None)
    refined = refined_result
    if refined is None:
        try:
            refined = Refiner(volume).refine(raw, parameters) if raw.path is not None else Recovery(error=raw.error, supports_cycles=False)
        except ValueError as error:
            refined = Recovery(error=str(error), supports_cycles=False)
    qc = quality_control(volume, refined, mode)
    diagnostics.update(refinement_success=refined.path is not None, final_qc_pass=qc['qc_pass'])
    # A failed final QC is returned explicitly, never silently replaced by raw.
    return refined, diagnostics, qc
