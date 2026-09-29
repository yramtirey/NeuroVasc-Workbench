"""Analytic graph truth, digital connectivity, multigraph reduction and matching."""
from copy import deepcopy
from pathlib import Path
import ast
import networkx as nx
import nibabel as nib
import numpy as np
import pytest
from scipy.ndimage import label,generate_binary_structure
from skimage.measure import euler_number
from neurovasc.geometry.topology_centerline import extract_topology_network,simple_point
from neurovasc.graph.topology_network import reduce_network,cycle_rank,fundamental_cycles
from neurovasc.graph.topology_qc import topology_qc
from validation.synthetic_topology.geometries import phantom,matrix
from validation.synthetic_topology.evaluate import score,recover_all,match_structures
from scripts.run_topology_validation import output_directory


@pytest.mark.parametrize('kind,expected',[
 ('simple_ring',(1,0,0,1,1)),('ring_branches',(1,2,2,4,1)),('figure_eight',(1,0,1,2,2)),
 ('communicating_bridge',(1,0,2,3,2)),('incomplete_loop',(1,2,0,1,0)),('y',(1,3,1,3,0)),
 ('near_touch',(2,4,0,2,0)),('small_bridge',(1,0,2,3,2)),('ring_multiple',(1,4,4,8,1))])
def test_analytic_truth(kind,expected):
    v=phantom(kind);m=v.metadata()
    assert tuple(m[k] for k in ('true_components','true_endpoints','true_junctions','true_branch_count','true_cycle_rank'))==expected
    assert len(m['cycle_basis'])==expected[-1]
    for cycle in m['cycle_basis']:
        g=nx.MultiGraph()
        for key in cycle:
            b=v.branches[key];g.add_edge(b['start'],b['end'])
        assert nx.is_connected(g) and all(degree==2 for _,degree in g.degree)
    if kind=='simple_ring':assert m['true_length_mm']==pytest.approx(16*np.pi)
    if kind=='communicating_bridge':assert m['true_length_mm']==50
    if kind=='y':assert m['true_length_mm']==pytest.approx(58)


@pytest.mark.parametrize('diameter',[1.,1.5,2.])
def test_small_bridge_truth_and_presence(diameter):
    v=phantom('small_bridge',bridge_diameter=diameter)
    bridges=[b for b in v.branches.values() if b['communicating_bridge']]
    assert len(bridges)==1 and bridges[0]['diameter_mm']==diameter
    assert bridges[0]['length_mm']==6
    assert v.volume.data.any()
    row,*_=score(v,extract_topology_network(v.volume))
    assert row['cycle_rank_correct'] and row['communicating_bridge_recovered']


def test_near_touch_digital_adjacency_threshold():
    near=phantom('near_touch',separation=.4)
    assert label(near.volume.data,generate_binary_structure(3,1))[1]==2
    assert label(near.volume.data,np.ones((3,3,3)))[1]==1
    g=extract_topology_network(near.volume)
    assert nx.number_connected_components(g)==2
    assert score(near,g)[0]['branch_connectivity_correct']
    far=phantom('near_touch',separation=1.2)
    assert label(far.volume.data,np.ones((3,3,3)))[1]==2


def make_graph(points,edges):
    g=nx.Graph()
    for i,p in enumerate(points):g.add_node(i,world=tuple(p))
    for a,b in edges:g.add_edge(a,b,length_mm=float(np.linalg.norm(points[a]-points[b])))
    return g


def test_closed_component_parallel_branches_and_cycle_basis():
    theta=np.linspace(0,2*np.pi,64,endpoint=False)
    points=np.column_stack((8*np.cos(theta),np.zeros(64),8*np.sin(theta)))
    g=make_graph(points,[(i,(i+1)%64) for i in range(64)])
    net=reduce_network(g)
    assert net.counts()['branches']==1 and len(net.cycles)==1
    assert nx.number_of_selfloops(net.logical)==1
    assert net.counts()['endpoints']==net.counts()['junctions']==0
    row,branches,cycles,*_=score(phantom('simple_ring'),g)
    assert row['branch_connectivity_correct'] and row['topology_correct']
    assert cycles[0]['membership_correct']
    # A theta graph must retain all three parallel branches, not one nx.Graph edge.
    points=np.array([[-2,0,0],[2,0,0],[0,0,2],[0,0,0],[0,0,-2]],float)
    g=make_graph(points,[(0,2),(2,1),(0,3),(3,1),(0,4),(4,1)])
    net=reduce_network(g)
    assert net.counts()['branches']==3 and cycle_rank(net.logical)==2
    assert net.logical.number_of_edges('J0','J1')==3


def test_junction_cluster_internal_cycle_is_not_pruned():
    points=np.array([[0,0,0],[1,0,0],[.5,1,0],[-1,0,0],[2,0,0],[.5,2,0]],float)
    g=make_graph(points,[(0,1),(1,2),(2,0),(0,3),(1,4),(2,5)])
    net=reduce_network(g)
    assert net.counts()['junctions']==1 and net.counts()['cycle_rank']==1
    assert len(net.metadata['cycles_pruned'])==0
    assert any(b['kind']=='junction_internal_cycle' for b in net.branches.values())
    assert all(c['passed'] for c in topology_qc(net) if c['critical'])


@pytest.mark.parametrize('spacing',[(.5,.5,.5),(.46875,.46875,.8)])
@pytest.mark.parametrize('kind',['simple_ring','ring_branches','figure_eight','communicating_bridge','incomplete_loop'])
def test_cycles_and_euler_preserved_on_declared_phantoms(kind,spacing):
    v=phantom(kind,spacing)
    g=extract_topology_network(v.volume);net=reduce_network(g)
    assert net.counts()['cycle_rank']==v.metadata()['true_cycle_rank']
    skeleton=np.zeros_like(v.volume.data)
    for p in g:skeleton[p]=1
    assert euler_number(skeleton,connectivity=1)==euler_number(v.volume.data,connectivity=1)
    assert net.counts()['components']==label(v.volume.data,generate_binary_structure(3,1))[1]


def test_deterministic_rotated_affine():
    v=phantom('ring_branches',spacing=(.46875,.46875,.8))
    a=extract_topology_network(v.volume);b=extract_topology_network(v.volume)
    assert list(a.edges)==list(b.edges)
    angle=.43
    transform=np.eye(4);transform[:3,:3]=[[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]]
    transform[:3,3]=[5,-4,9]
    volume=deepcopy(v.volume);volume.affine=transform@volume.affine
    rotated=extract_topology_network(volume)
    assert list(a.nodes)==list(rotated.nodes) and list(a.edges)==list(rotated.edges)
    for n in a:np.testing.assert_allclose(rotated.nodes[n]['world'],nib.affines.apply_affine(transform,a.nodes[n]['world']),atol=1e-12)


def test_no_truth_or_lee_dependency_in_new_extractor(monkeypatch):
    import neurovasc.geometry.centerline as baseline
    monkeypatch.setattr(baseline,'extract_centerline',lambda *a:(_ for _ in ()).throw(AssertionError('Lee called')))
    v=phantom('communicating_bridge')
    assert cycle_rank(extract_topology_network(v.volume))==2
    tree=ast.parse((Path(__file__).resolve().parents[1]/'neurovasc/geometry/topology_centerline.py').read_text())
    for node in ast.walk(tree):
        if isinstance(node,ast.ImportFrom):
            assert not any(s in (node.module or '') for s in ('validation','geodesic','centerline'))


def test_assignment_is_independent_of_node_and_branch_order():
    v=phantom('ring_branches')
    g=recover_all(v.volume)['lee'][0];net=reduce_network(g)
    first=score(v,g)[0]
    reordered=nx.Graph()
    for n,d in reversed(list(g.nodes(data=True))):reordered.add_node(('renamed',n),**d)
    for a,b,d in reversed(list(g.edges(data=True))):reordered.add_edge(('renamed',a),('renamed',b),**d)
    second=score(v,reordered)[0]
    assert first['topology_correct'] and second['topology_correct']
    assert first['endpoint_position_error_mm']==pytest.approx(second['endpoint_position_error_mm'])
    assert first['junction_position_error_mm']==pytest.approx(second['junction_position_error_mm'])


def test_empty_and_size_limit():
    volume=phantom('simple_ring').volume
    with pytest.raises(ValueError,match='size_limit'):extract_topology_network(volume,maximum_voxels=1)
    volume.data[:]=0
    with pytest.raises(ValueError,match='empty'):extract_topology_network(volume)


def test_matrix_counts_and_output_guard():
    experiments=list(matrix())
    assert len(experiments)==76 and len({e for e,_ in experiments})==76
    assert {v.phase for _,v in experiments}=={0.,.25,.5}
    assert sum(v.orientation=='rotated' for _,v in experiments)==4
    root=Path(__file__).resolve().parents[1]/'validation/outputs'
    for folder in ('caliber_v1','caliber_v2','centerline_v1','centerline_v2','topology_v1/../centerline_v1'):
        with pytest.raises(ValueError):output_directory(root/folder)
    assert output_directory(root/'topology_v1')


def test_missing_vessel_is_not_a_false_connection():
    v=phantom('near_touch',separation=1.2)
    graph=extract_topology_network(v.volume)
    one=graph.subgraph(next(iter(nx.connected_components(graph)))).copy()
    row,*_=score(v,one)
    assert row['recovered_components']==1 and row['true_components']==2
    assert row['false_connection_count']==0 and row['missed_true_components']==1
    assert row['spurious_cross_component_edges']==0


def test_failed_extractions_retain_missed_structures_and_qc():
    v=phantom('small_bridge',bridge_diameter=1.)
    row,branches,cycles,qc,network,nodes=score(v,None,'controlled_failure',False)
    assert not row['extraction_success'] and not row['communicating_bridge_recovered']
    assert len(branches)==3 and len(cycles)==2 and network is None
    assert all(not c['cycle_preserved'] for c in cycles)
    assert qc[0]['critical'] and not qc[0]['passed']
