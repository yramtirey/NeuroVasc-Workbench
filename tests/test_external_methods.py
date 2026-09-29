"""Comparator contracts tested without requiring an installed VMTK runtime."""
import json
from pathlib import Path
import networkx as nx
import numpy as np
import pytest
import vtk
from vtk.util.numpy_support import numpy_to_vtk
from external.vmtk.surface_preparation import export_surface,read_vtp,write_vtp,snap_seeds
from external.vmtk.centerlines import parse
from external.vmtk.cross_sections import status
from validation.external_methods.matching import errors,matched,sample_on_truth,endpoints
from validation.external_methods.evaluate import scope,DIRECT,PARTIAL,MISMATCH,graph_from_lines,centerline
from validation.external_methods.report import output_directory,ALLOWED,aggregates
from validation.synthetic_caliber.geometries import straight_cylinder
from validation.synthetic_centerline.truth import sample_truth
from scripts.run_external_methods_benchmark import automatic_seeds,oracle_seeds,probe,execute
from neurovasc.geometry.centerline_qc import Recovery


def line_data(radius=True):
    data=vtk.vtkPolyData();points=vtk.vtkPoints()
    for p in ((0,0,-1),(0,0,0),(0,0,1)):points.InsertNextPoint(p)
    data.SetPoints(points);cells=vtk.vtkCellArray();cells.InsertNextCell(3)
    for i in range(3):cells.InsertCellPoint(i)
    data.SetLines(cells)
    if radius:
        array=numpy_to_vtk(np.array([1.,2.,1.]));array.SetName('MaximumInscribedSphereRadius');data.GetPointData().AddArray(array)
    return data


def test_vtp_roundtrip(tmp_path):
    path=tmp_path/'line.vtp';write_vtp(line_data(),path)
    result=parse(read_vtp(path))[0]
    assert result['radius_mm']==[1.,2.,1.]
    np.testing.assert_allclose(result['points'],[[0,0,-1],[0,0,0],[0,0,1]])


def test_missing_radius_rejected():
    with pytest.raises(ValueError,match='missing_Maximum'):parse(line_data(False))


def test_invalid_radius_rejected():
    data=line_data();data.GetPointData().GetArray('MaximumInscribedSphereRadius').SetValue(0,float('nan'))
    with pytest.raises(ValueError,match='invalid_centerline'):parse(data)


def test_degenerate_external_line_is_not_success():
    data=line_data()
    for i in range(3):data.GetPoints().SetPoint(i,0,0,0)
    with pytest.raises(ValueError,match='degenerate_centerline'):parse(data)


def test_surface_world_coordinates(tmp_path):
    vessel=straight_cylinder(diameter_mm=4,spacing=(.46875,.46875,.8))
    rotation=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
    vessel.volume.affine[:3,:3]=rotation@vessel.volume.affine[:3,:3]
    vessel.volume.affine[:3,3]=rotation@vessel.volume.affine[:3,3]+[11,22,33]
    path=tmp_path/'surface.vtp';parameters=export_surface(vessel.volume,path);surface=read_vtp(path)
    bounds=np.array(surface.GetBounds()).reshape(3,2)
    np.testing.assert_allclose(bounds.mean(axis=1),[11,22,33],atol=1e-6)
    assert parameters['isovalue']==.5 and not parameters['smoothing'] and not parameters['decimation']
    assert surface.GetNumberOfPolys()>0


def test_seed_conversion():
    data=line_data();ids,points=snap_seeds(data,[[.01,0,-1],[0,0,1]])
    assert ids==[0,2] and points==[[0.,0.,-1.],[0.,0.,1.]]


@pytest.mark.parametrize('points',[[[0,0,0]],[[0,0,0],[0,0,0]],[[0,0,np.nan],[0,0,1]],[[0,0],[1,1]]])
def test_invalid_seeds(points):
    with pytest.raises(ValueError):snap_seeds(line_data(),points)


def test_automatic_endpoints_need_no_truth():
    graph=nx.path_graph(3)
    for n in graph:graph.nodes[n]['world']=[0,0,2-n]
    seeds,error=automatic_seeds(Recovery(graph=graph))
    assert seeds==[[0.,0.,0.],[0.,0.,2.]] and not error


def test_automatic_failure_explicit():
    seeds,error=automatic_seeds(Recovery(error='empty_mask'))
    assert seeds==[] and error=='empty_mask'


def test_oracle_seeds_are_analytic():
    vessel=straight_cylinder(length_mm=30)
    assert oracle_seeds(vessel)==[[0.,0.,-15.],[0.,0.,15.]]


@pytest.mark.parametrize('kind,expected',[('straight',DIRECT),('oblique',DIRECT),('curved',DIRECT),('tapered',DIRECT),('bifurcation',PARTIAL),('simple_ring',MISMATCH),('ring_branches',MISMATCH),('figure_eight',MISMATCH),('communicating_bridge',MISMATCH)])
def test_scope(kind,expected):assert scope(kind)==expected


def test_common_support_excludes_missing_samples():
    r=matched([2.,np.nan,4.],[3.,3.,np.nan],[2.,2.,2.])
    assert r['common_sample_count']==1 and r['left_mae_mm']==0 and r['right_mae_mm']==1
    assert r['common_support_fraction']==pytest.approx(1/3)


def test_empty_match_has_no_accuracy():
    r=matched([np.nan],[2.],[2.]);assert r['status']=='unscorable' and r['left_mae_mm'] is None


def test_radius_error_formula():
    r=errors([1.,3.],[2.,2.])
    assert r['mae_mm']==1 and r['rmse_mm']==1 and r['bias_mm']==0 and r['percent_error']==50


def test_coordinate_sampling_reversal_and_missingness():
    vessel=straight_cylinder();p=sample_truth(vessel,samples=121);values=np.full(121,4.)
    a=sample_on_truth(vessel,p,values,p,.5);b=sample_on_truth(vessel,p[::-1],values[::-1],p,.5)
    np.testing.assert_allclose(a,b)
    short=sample_on_truth(vessel,p[40:80],values[40:80],p,.5)
    assert np.isnan(short[:40]).all() and np.isnan(short[80:]).all()
    values[60]=np.nan
    assert np.isnan(sample_on_truth(vessel,p,values,p,.5)[60])


def test_endpoint_counts_expose_missing():
    r=endpoints([[0,0,0]],[[0,0,0],[0,0,1]])
    assert r['matched_endpoint_count']==1 and r['missing_endpoint_count']==1


def test_shared_path_length_not_double_counted():
    lines=[dict(points=[[0,0,0],[0,0,1],[0,0,2]]),dict(points=[[0,0,0],[0,0,1],[1,0,2]])]
    g=graph_from_lines(lines)
    assert g.number_of_edges()==3
    assert sum(d['length_mm'] for *_,d in g.edges(data=True))==pytest.approx(2+np.sqrt(2))


def test_exact_centerline_metrics():
    vessel=straight_cylinder();p=sample_truth(vessel,samples=121)
    row,_=centerline(vessel,graph_from_lines([dict(points=p)]),[p])
    assert row['coverage_fraction']==1 and row['length_error_mm']<1e-10 and row['mean_position_error_mm']<1e-10


def test_failure_inclusive_and_unavailable_separate():
    base=dict(method='vmtk',oracle_seed=True,comparison_scope=DIRECT)
    rows=[dict(base,status='success',extraction_success=True,coverage_fraction=.8),dict(base,status='failed',extraction_success=False)]
    r=next(iter(aggregates(rows).values()))
    assert r['failure_inclusive_coverage']==.4 and r['successful_support_coverage']==.8
    r=next(iter(aggregates([dict(base,status='unavailable',extraction_success=None)]).values()))
    assert r['success_rate'] is None and r['failure_inclusive_coverage'] is None


def test_missing_external_environment(tmp_path):
    assert not probe(tmp_path/'missing-python')['available']


def test_no_fabricated_execution(tmp_path,monkeypatch):
    import scripts.run_external_methods_benchmark as runner
    monkeypatch.setattr(runner,'probe',lambda _:dict(available=False,reason='test_missing_vmtk'))
    request=dict(geometry_id='test',oracle_seed=True,comparison_scope=DIRECT,surface='unused')
    (tmp_path/'request.json').write_text(json.dumps(request));(tmp_path/'requests.json').write_text('["request.json"]')
    execute(tmp_path,tmp_path/'missing')
    r=json.loads((tmp_path/'vmtk_outputs/test/oracle/result.json').read_text())
    assert r['status']=='unavailable' and r['extraction_success'] is None
    assert 'lines' not in r and 'mae_mm' not in r and r['oracle_seed'] is True


def test_rerun_removes_only_known_worker_outputs(tmp_path):
    from scripts.run_external_methods_benchmark import clear_job_outputs
    for name in ('surface.vtp','centerlines.vtp','result.json','notes.txt'):
        (tmp_path/name).write_text('existing')
    clear_job_outputs(tmp_path)
    assert not (tmp_path/'centerlines.vtp').exists()
    assert (tmp_path/'result.json').exists()
    clear_job_outputs(tmp_path,include_result=True)
    assert not (tmp_path/'result.json').exists()
    assert (tmp_path/'surface.vtp').read_text()=='existing' and (tmp_path/'notes.txt').read_text()=='existing'


def test_slicer_ce_explicitly_unavailable():
    r=status();assert not r['executed'] and r['version'] is None and 'ce_diameter_mm' in r['expected_output_columns']


def test_absent_external_branch_groups_have_no_fabricated_connectivity():
    from validation.external_methods.evaluate import vmtk_y_groups
    from validation.synthetic_centerline.geometries import y_bifurcation
    details,metrics=vmtk_y_groups(y_bifurcation(),{})
    assert not details and metrics['branch_status']=='failed'
    assert metrics['branch_f1'] is None and metrics['recovered_branch_count'] is None


def test_repeatability_comparison_exposes_missing_and_changed_values():
    from scripts.verify_external_methods_repeatability import compare
    assert not compare({'a':[1.,2.]},{'a':[1.,2.+1e-12]})
    assert compare({'a':[1.,2.]},{'a':[1.]})
    assert compare({'a':[1.,2.]},{'a':[1.,2.1]})
    assert compare({'a':None},{'a':0.})


@pytest.mark.parametrize('path',['validation/outputs/caliber_v1','validation/outputs/topology_v2','outputs/cases/case_001','validation/outputs/external_methods_v1/../topology_v1'])
def test_output_protection(path):
    with pytest.raises(ValueError):output_directory(Path(__file__).resolve().parents[1]/path)


def test_output_allowed():assert output_directory(ALLOWED/'new_run')==ALLOWED/'new_run'


def test_symlink_escape(tmp_path):
    import tempfile
    ALLOWED.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=ALLOWED) as d:
        link=Path(d)/'escape';link.symlink_to(tmp_path,target_is_directory=True)
        with pytest.raises(ValueError):output_directory(link/'run')
