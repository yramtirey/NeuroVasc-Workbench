"""Mask-only hybrid QC, geometric refinement and preserved regression contracts."""
from copy import deepcopy
from dataclasses import replace
import ast
from pathlib import Path
import networkx as nx
import numpy as np
import pytest
from neurovasc.geometry.centerline_qc import Recovery, quality_control, occupancy, smoothness
from neurovasc.geometry.centerline_refinement import Refiner, RefinementParameters, resample, distance
from neurovasc.geometry.hybrid_centerline import recover, hybrid, CANDIDATE
from validation.synthetic_caliber.geometries import straight_cylinder, curved_vessel, oblique_cylinder
from validation.synthetic_caliber.experiments_v2 import phase_shifted_cylinder
from validation.synthetic_centerline.geometries import y_bifurcation, matrix
from validation.synthetic_centerline.experiments_v2 import matrix_v2, loop_vessel
from validation.synthetic_centerline.evaluate import score_points, topology
from validation.synthetic_centerline.comparison_v2 import measure
from scripts.run_centerline_validation_v2 import configurations, output_directory


def test_good_lee_qc_and_deterministic_selection():
    volume=straight_cylinder(3).volume
    baseline=recover(volume,'lee')
    qc=quality_control(volume,baseline)
    assert qc['qc_pass']
    first,d,q=hybrid(volume)
    second,e,_=hybrid(volume)
    assert d==e and d['selected_method']=='lee' and not d['fallback_used']
    np.testing.assert_array_equal(first.path,baseline.path)
    np.testing.assert_array_equal(first.path,second.path)
    assert not d['production_eligible'] and not d['fallback_supports_cycles']


def test_short_failed_and_malformed_path_qc():
    volume=straight_cylinder(3).volume
    assert not quality_control(volume,Recovery(error='empty'))['qc_pass']
    baseline=recover(volume,'lee')
    shortened=baseline.path[20:25]
    selected={tuple(p) for p in shortened}
    nodes=[n for n,d in baseline.graph.nodes(data=True) if tuple(d['world']) in selected]
    short=Recovery(baseline.graph.subgraph(nodes).copy(),shortened)
    qc=quality_control(volume,short)
    assert not qc['qc_pass']
    assert 'mask_extent_coverage' in qc['qc_fail_reasons']
    recovered, decision, final_qc = hybrid(volume, lee_result=short)
    assert decision['lee_extraction_success'] and decision['fallback_used']
    assert 'mask_extent_coverage' in decision['fallback_reason']
    assert final_qc['qc_pass'] and recovered.path is not None
    broken=deepcopy(baseline)
    broken.path[2]=broken.path[1]
    assert 'repeated_path_points' in quality_control(volume,broken)['qc_fail_reasons']
    broken=deepcopy(baseline)
    broken.path=broken.path[::2]
    assert 'ordered_graph_path' in quality_control(volume,broken)['qc_fail_reasons']


@pytest.mark.parametrize('spacing',[(.5,.5,.5),(.46875,.46875,.8)])
def test_known_failure_triggers_fallback_and_improves(spacing):
    vessel=phase_shifted_cylinder(1,spacing,(.5,.5,.5))
    raw=recover(vessel.volume,'geodesic_raw')
    result,decision,qc=hybrid(vessel.volume)
    assert not decision['lee_extraction_success'] and decision['fallback_used']
    assert decision['refinement_success'] and qc['qc_pass']
    assert not result.supports_cycles
    before=score_points(vessel,raw.path)[0]
    after=score_points(vessel,result.path)[0]
    assert after['mean_position_error_mm']<.05  # 10% of smallest voxel, observed <0.005 mm.
    assert after['mean_position_error_mm']<before['mean_position_error_mm']
    assert after['mean_tangent_error_deg']<1  # Observed <0.1 degree, vs voxel-step raw path.
    assert after['coverage_fraction']>.95
    assert np.isfinite(result.path).all()
    assert len(np.unique(result.path,axis=0))==len(result.path)
    assert np.all(occupancy(vessel.volume,resample(result.path,.1))>=.5-1e-8)


def test_selection_does_not_import_or_call_truth(monkeypatch):
    import validation.synthetic_centerline.truth as truth
    monkeypatch.setattr(truth,'project',lambda *args: (_ for _ in ()).throw(AssertionError('truth called')))
    volume=phase_shifted_cylinder(1,(.5,.5,.5),(.5,.5,.5)).volume
    assert hybrid(volume)[1]['fallback_used']
    root=Path(__file__).resolve().parents[1]
    for name in ('centerline_qc','centerline_refinement','hybrid_centerline'):
        tree=ast.parse((root/f'neurovasc/geometry/{name}.py').read_text())
        for node in ast.walk(tree):
            if isinstance(node,ast.ImportFrom):
                assert not (node.module or '').startswith('validation')


@pytest.mark.parametrize('step',[.5,1.])
def test_uniform_physical_resampling(step):
    points=np.array([[0.,0,0],[0,0,.3],[0,0,1.7],[0,0,4.3]])
    result=resample(points,step)
    np.testing.assert_allclose(result[[0,-1]],points[[0,-1]])
    gaps=np.diff(distance(result))
    np.testing.assert_allclose(gaps,gaps[0],atol=1e-12)
    assert gaps.max()<=step+1e-12
    np.testing.assert_allclose(result[:,:2],0)


def test_straight_refinement_and_endpoint_control():
    vessel=straight_cylinder(1)
    raw=recover(vessel.volume,'geodesic_raw')
    refiner=Refiner(vessel.volume)
    without=refiner.refine(raw,replace(CANDIDATE,refine_endpoints=False))
    with_end=refiner.refine(raw,CANDIDATE)
    np.testing.assert_allclose(with_end.path[:,:2],0,atol=1e-12)
    assert smoothness(with_end.path)['mean_turn_deg']<1e-5
    assert distance(with_end.path)[-1]>=distance(without.path)[-1]
    assert np.all(occupancy(vessel.volume,with_end.path)>=.5-1e-8)


def test_curved_shape_is_not_flattened():
    vessel=curved_vessel(3)
    result=Refiner(vessel.volume).refine(recover(vessel.volume,'geodesic_raw'))
    metrics=score_points(vessel,result.path)[0]
    assert metrics['mean_position_error_mm']<.125  # Quarter of isotropic voxel, observed 0.028 mm.
    assert metrics['mean_tangent_error_deg']<2  # Observed 0.79 degrees.
    assert np.ptp(result.path[:,0])>3  # True arc sagittal excursion is 4.54 mm.
    assert smoothness(result.path)['cumulative_turn_deg']>40


@pytest.mark.parametrize('spacing',[(.5,.5,.5),(.46875,.46875,.8)])
def test_y_junction_is_fixed_and_topology_preserved(spacing):
    vessel=y_bifurcation(spacing)
    raw=recover(vessel.volume,'geodesic_raw')
    refined=Refiner(vessel.volume).refine(raw)
    assert topology(raw.graph)[0]['endpoint_count']==topology(refined.graph)[0]['endpoint_count']==3
    assert topology(refined.graph)[0]['cycle_count']==0
    before=[d['world'] for n,d in raw.graph.nodes(data=True) if raw.graph.degree[n]>=3]
    after=[d['world'] for n,d in refined.graph.nodes(data=True) if refined.graph.degree[n]>=3]
    np.testing.assert_array_equal(before,after)
    assert measure(vessel,refined)[0]['topology_preserved']


@pytest.mark.parametrize('attached',[False,True])
def test_loop_is_not_reduced_to_tree_truth(attached):
    vessel=loop_vessel(attached=attached)
    lee=recover(vessel.volume,'lee')
    raw=recover(vessel.volume,'geodesic_raw')
    assert measure(vessel,lee)[3]['topology_preserved']
    loop=measure(vessel,raw)[3]
    assert loop['true_cycle_count']==1 and loop['cycle_count']==0
    assert not loop['topology_preserved'] and not raw.supports_cycles
    assert (loop['mean_junction_error_mm'] is not None) == attached
    assert loop['unmatched_junction_count'] == (0 if attached else 2)
    assert loop['network_coverage_fraction']<1
    with pytest.raises(ValueError,match='supports_cycles=false'):
        Refiner(vessel.volume).refine(lee)
    selected,decision,_=hybrid(vessel.volume)
    assert decision['selected_method']=='lee'
    assert measure(vessel,selected)[3]['topology_preserved']
    assert not decision['production_eligible']


def test_matrix_reuses_all_v1_ids_and_masks():
    old=list(matrix());new=list(matrix_v2())
    assert len(new)==64
    for (oldid,a),(newid,b) in zip(old,new):
        assert oldid==newid
        np.testing.assert_array_equal(a.vessel.volume.data,b.vessel.volume.data)
        np.testing.assert_array_equal(a.vessel.volume.affine,b.vessel.volume.affine)


def test_declared_sweep_and_archive_protection():
    configs=configurations()
    assert len(configs)==15 and sum(p==CANDIDATE for p in configs.values())==1
    root=Path(__file__).resolve().parents[1]/'validation/outputs'
    for forbidden in ('centerline_v1','caliber_v1','caliber_v2','centerline_v2/../centerline_v1'):
        with pytest.raises(ValueError):output_directory(root/forbidden)
    assert output_directory(root/'centerline_v2/check')


def test_qc_simple_mode_rejects_cycle():
    vessel=loop_vessel()
    qc=quality_control(vessel.volume,recover(vessel.volume,'lee'),mode='simple')
    assert 'simple_mode_cycles' in qc['qc_fail_reasons']


def test_invalid_resampling_and_qc_modes():
    with pytest.raises(ValueError):resample(np.zeros((3,3)),.5)
    with pytest.raises(ValueError):resample(np.array([[0,0,0],[0,0,1]]),0)
    with pytest.raises(ValueError):quality_control(straight_cylinder(1).volume,Recovery(),mode='unknown')


def test_qc_accepts_turn_threshold_with_roundoff():
    volume=oblique_cylinder(3,direction=(1,2,1)).volume
    baseline=recover(volume,'lee')
    qc=quality_control(volume,baseline)
    turning=next(c for c in qc['checks'] if c['criterion']=='mean_turn_deg')
    assert turning['value']==pytest.approx(60,abs=1e-10)
    assert turning['passed']
    assert hybrid(volume)[1]['selected_method']=='lee'
