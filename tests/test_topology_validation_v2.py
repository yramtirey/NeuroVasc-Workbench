"""Fixed topology contracts, geometric containment and archived network regressions."""
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import ast
import networkx as nx
import nibabel as nib
import numpy as np
import pytest
from neurovasc.geometry.topology_centerline import extract_topology_network
from neurovasc.geometry.topology_refinement import refine_network, CANDIDATE, geometric_checks
from neurovasc.graph.topology_network import reduce_network, cycle_rank
from neurovasc.graph.topology_invariants import check_invariants, materialize
from neurovasc.graph.junction_refinement import classify_spurs
from neurovasc.geometry.centerline_refinement import Refiner, resample
from neurovasc.geometry.centerline_qc import occupancy
from validation.synthetic_topology.geometries import phantom
from validation.synthetic_topology.experiments_v2 import junction_phantom, matrix_v2
from validation.synthetic_topology.evaluate_v2 import measure, prf
from scripts.run_topology_validation_v2 import configurations, output_directory

ROOT=Path(__file__).resolve().parents[1]


def recover(p): return reduce_network(extract_topology_network(p.volume))


@pytest.mark.parametrize('kind', ['simple_ring','figure_eight','incomplete_loop','ring_branches'])
@pytest.mark.parametrize('spacing',[(.5,.5,.5),(.46875,.46875,.8)])
def test_fixed_topology_and_inside_lumen(kind,spacing):
    p=phantom(kind,spacing);original=recover(p);snapshot=deepcopy(original)
    result=refine_network(p.volume,original)
    assert result.accepted,result.failure_reason
    assert result.invariants['topology_invariant_pass']
    assert set(result.network.branches)==set(original.branches)
    assert result.network.cycles==original.cycles
    assert result.network.counts()=={**original.counts(),'nodes':len(result.network.raw),'edges':result.network.raw.number_of_edges()}
    for key,b in result.network.branches.items():
        np.testing.assert_array_equal(original.branches[key]['points'],snapshot.branches[key]['points'])
        np.testing.assert_array_equal(b['points'][0],result.network.logical.nodes[b['start']]['world'])
        np.testing.assert_array_equal(b['points'][-1],result.network.logical.nodes[b['end']]['world'])
        assert np.all(occupancy(p.volume,resample(b['points'],.05))>=.5-1e-8)


@pytest.mark.parametrize('phase',[0.,.25,.5])
@pytest.mark.parametrize('spacing',[(.5,.5,.5),(.46875,.46875,.8)])
@pytest.mark.parametrize('gap',[.4,1.2])
def test_entire_near_touch_cohort_no_false_bridge(phase,spacing,gap):
    p=phantom('near_touch',spacing,phase,separation=gap);n=recover(p);r=refine_network(p.volume,n)
    assert r.network.counts()['components']==2
    assert measure(p,r.network)[0]['false_connection_count']==0
    if r.accepted: assert r.invariants['topology_invariant_pass']


@pytest.mark.parametrize('diameter',[1.,1.5,2.])
@pytest.mark.parametrize('phase',[0.,.25,.5])
@pytest.mark.parametrize('spacing',[(.5,.5,.5),(.46875,.46875,.8)])
def test_all_small_bridges_remain_cyclic(diameter,phase,spacing):
    p=phantom('small_bridge',spacing,phase,bridge_diameter=diameter);n=recover(p);r=refine_network(p.volume,n)
    row,*_=measure(p,r.network)
    assert r.network.counts()['cycle_rank']==2
    assert row['communicating_bridge_recovered']
    assert not any(s['spur_removed'] for s in r.spur_diagnostics)


@pytest.mark.parametrize('kind,count,rank',[('four_way',1,0),('double_junction',2,0),('asymmetric_figure_eight',1,2),('short_communicator',2,2),('near_junction_branch',2,0)])
def test_harder_truth_and_junction_identity(kind,count,rank):
    p=junction_phantom(kind);assert p.metadata()['true_junctions']==count
    assert cycle_rank(p.graph)==rank
    n=recover(p);r=refine_network(p.volume,n)
    assert r.network.counts()['junctions']==n.counts()['junctions']
    if kind in ('four_way','double_junction'):assert n.counts()['junctions']==count
    assert not any(d['merge_performed'] for d in r.junction_diagnostics)
    assert set(n.branches)==set(r.network.branches)
    assert not any(s['spur_removed'] for s in r.spur_diagnostics)


def test_reject_topology_corruption_and_keep_original():
    p=phantom('ring_branches');n=recover(p);r=refine_network(p.volume,n)
    assert r.accepted
    corrupt=deepcopy(r.network);key=next(iter(corrupt.branches));b=corrupt.branches.pop(key)
    corrupt.logical.remove_edge(b['start'],b['end'],key)
    assert not check_invariants(n,corrupt)['topology_invariant_pass']
    corrupt=deepcopy(r.network);corrupt.branches[key]['points'][0]+=10
    assert not check_invariants(n,corrupt)['topology_invariant_pass']
    corrupt=deepcopy(r.network);corrupt.cycles=[]
    assert not check_invariants(n,corrupt)['topology_invariant_pass']


def test_rejected_geometry_is_explicit_and_original_unmodified(monkeypatch):
    import neurovasc.geometry.topology_refinement as module
    p=phantom('simple_ring');n=recover(p)
    monkeypatch.setattr(module,'geometric_checks',lambda *args: (_ for _ in ()).throw(ValueError('controlled_outside_lumen')))
    result=module.refine_network(p.volume,n)
    assert not result.accepted and result.network is n and result.original is n
    assert result.failure_reason=='controlled_outside_lumen'


def test_obvious_short_spur_marked_but_identity_protected():
    p=phantom('ring_branches');n=recover(p)
    b=next(b for b in n.branches.values() if n.logical.nodes[b['end']]['node_kind']=='endpoint')
    origin=np.array(n.logical.nodes[b['start']]['world']);n.logical.nodes[b['end']]['world']=(origin+np.array([.1,0,0])).tolist()
    b['points']=np.array([origin,origin+[.1,0,0]]);b['length_mm']=.1
    records=classify_spurs(n,Refiner(p.volume),1.)
    item=next(s for s in records if s['branch_id']==b['branch_id'])
    assert item['spur_candidate'] and not item['spur_removed']
    assert 'identity' in item['decision']


def test_branch_adjacency_and_junction_degree_metrics():
    assert prf(2,1,2)==pytest.approx((2/3,1/2,4/7))
    p=junction_phantom('four_way');n=recover(p);row,branches,cycles,nodes,_=measure(p,n)
    assert row['branch_f1']==row['adjacency_f1']==1
    assert row['adjacency_true_positive']==4 and row['adjacency_false_positive']==0
    j=next(x for x in nodes if x['node_kind']=='junction')
    assert j['degree_true']==j['degree_recovered']==4 and j['degree_correct']
    assert all(b['assignment_correct'] for b in branches)


def test_physical_resampling_determinism_and_affine_rotation():
    p=phantom('ring_branches');n=recover(p)
    a=refine_network(p.volume,n);b=refine_network(p.volume,n)
    assert a.accepted and b.accepted
    angle=.43;rotation=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    volume=deepcopy(p.volume);transform=np.eye(4);transform[:3,:3]=rotation;transform[:3,3]=[5,-4,9]
    volume.affine=transform@volume.affine
    rotated=recover(type('P',(),{'volume':volume})());c=refine_network(volume,rotated)
    assert c.accepted
    for key in a.network.branches:
        np.testing.assert_array_equal(a.network.branches[key]['points'],b.network.branches[key]['points'])
        expected=nib.affines.apply_affine(transform,a.network.branches[key]['points'])
        np.testing.assert_allclose(c.network.branches[key]['points'],expected,atol=1e-8)
        sample=resample(a.network.branches[key]['points'],.5)
        assert np.linalg.norm(np.diff(sample,axis=0),axis=1).max()<=.5+1e-10


def test_exact_76_archived_masks_and_output_protection():
    import hashlib, json
    import numpy as np
    archive=ROOT/'validation/outputs/topology_v1'
    cases=list(matrix_v2(archive if (archive/'masks').is_dir() else None))
    fingerprints=json.loads((ROOT/'validation/public_results/topology_v1/mask_fingerprints.json').read_text())
    for eid,phantom,cohort in cases:
        if cohort!='core_v1':continue
        expected=fingerprints[eid]
        data=np.asarray(phantom.volume.data,dtype=np.uint8)
        assert list(data.shape)==expected['shape']
        assert hashlib.sha256(data.tobytes(order='C')).hexdigest()==expected['array_sha256']
        np.testing.assert_allclose(phantom.volume.affine,expected['affine'],atol=1e-6)
    assert len(cases)==106 and sum(c=='core_v1' for _,_,c in cases)==76
    for archive in ('topology_v1','centerline_v1','centerline_v2','caliber_v1','caliber_v2'):
        with pytest.raises(ValueError):output_directory(ROOT/'validation/outputs'/archive)
    assert len(configurations())==7


def test_no_truth_access_and_invalid_parameters():
    for file in ('neurovasc/geometry/topology_refinement.py','neurovasc/graph/junction_refinement.py','neurovasc/graph/topology_invariants.py'):
        tree=ast.parse((ROOT/file).read_text())
        assert all(not (n.module or '').startswith('validation') for n in ast.walk(tree) if isinstance(n,ast.ImportFrom))
    p=phantom('simple_ring');n=recover(p)
    with pytest.raises(ValueError):refine_network(p.volume,n,replace(CANDIDATE,resample_spacing_mm=0))
    with pytest.raises(ValueError):refine_network(p.volume,n,replace(CANDIDATE,endpoint_mode='unknown'))


def test_all_76_frozen_raw_counts_and_refined_cycle_ranks():
    import csv
    with (ROOT/'validation/public_results/topology_v1/results.csv').open() as stream:
        archived={r['geometry_id']:r for r in csv.DictReader(stream) if r['method']=='topology'}
    for eid,p,cohort in matrix_v2():
        if cohort!='core_v1':continue
        n=recover(p);counts=n.counts();old=archived[eid]
        for metric,field in [('components','components'),('endpoints','endpoint_count'),('junctions','junction_count'),('branches','branch_count'),('cycle_rank','cycle_rank')]:
            assert counts[metric]==int(old['recovered_'+field]),eid
        result=refine_network(p.volume,n)
        assert result.accepted,(eid,result.failure_reason)
        assert result.invariants['topology_invariant_pass'],eid
        assert result.network.counts()['cycle_rank']==int(old['true_cycle_rank']),eid


def test_geometric_contact_rejected_even_with_distinct_node_ids():
    p=phantom('near_touch',separation=.4);n=recover(p);r=refine_network(p.volume,n)
    corrupt=deepcopy(r.network);branches=list(corrupt.branches.values())
    a,b=branches[:2];shift=b['points'][0]-a['points'][0]
    a['points']+=shift
    for node in (a['start'],a['end']):corrupt.logical.nodes[node]['world']=(np.array(corrupt.logical.nodes[node]['world'])+shift).tolist()
    corrupt.raw=materialize(corrupt)
    assert nx.number_connected_components(corrupt.raw)==2
    with pytest.raises(ValueError,match='cross_component_geometric_contact'):
        geometric_checks(p.volume,corrupt,n)


def test_all_endpoint_modes_keep_topology_and_report_failures():
    p=phantom('incomplete_loop');n=recover(p)
    for mode in ('raw','recentered','tangent','boundary'):
        result=refine_network(p.volume,n,replace(CANDIDATE,endpoint_mode=mode))
        assert result.network.counts()['endpoints']==2
        assert result.network.counts()['cycle_rank']==0
        if result.accepted:assert result.invariants['topology_invariant_pass']
        else:assert result.failure_reason and result.network is n
