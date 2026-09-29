"""Actual subprocess integration; explicit skip when optional runtime unavailable."""
import json
import subprocess
import tempfile
from pathlib import Path
import pytest
from scripts.run_external_methods_benchmark import probe,ROOT,oracle_seeds
from validation.external_methods.report import ALLOWED
from validation.external_methods.evaluate import scope
from validation.synthetic_caliber.geometries import straight_cylinder
from validation.synthetic_centerline.geometries import y_bifurcation
from external.vmtk.surface_preparation import export_surface


@pytest.fixture(scope='module')
def external_python():
    path=ROOT/'external/vmtk/.venv/bin/python'
    info=probe(path)
    if not info.get('available'):pytest.skip('Actual VMTK unavailable: '+info.get('reason','probe failed'))
    return path


@pytest.mark.parametrize('kind',['straight','y'])
def test_actual_compiled_centerlines(external_python,kind):
    vessel=straight_cylinder(diameter_mm=4) if kind=='straight' else y_bifurcation()
    ALLOWED.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='integration_',dir=ALLOWED) as d:
        out=Path(d);surface=out/'surface.vtp';export_surface(vessel.volume,surface)
        request=dict(geometry_id='integration_'+kind,surface=str(surface),oracle_seed=True,
                     comparison_scope=scope(vessel.geometry),seeds_mm=oracle_seeds(vessel))
        path=out/'request.json';path.write_text(json.dumps(request))
        p=subprocess.run([str(external_python),'-m','external.vmtk.adapter','--request',str(path),'--output',str(out/'result')],cwd=ROOT,capture_output=True,text=True,timeout=120)
        assert p.returncode==0,p.stderr
        r=json.loads((out/'result/result.json').read_text())
        assert r['extraction_success'],r.get('failure_reason')
        assert r['oracle_seed'] is True and len(r['lines'])==(1 if kind=='straight' else 2)
        assert all(v>0 for line in r['lines'] for v in line['radius_mm'])
        assert not r.get('branch_failure_reason')
        geometry=r['branches']['geometry']
        assert len(geometry['GroupIds'])==(1 if kind=='straight' else 3)
        if kind=='straight':assert abs(geometry['Tortuosity'][0])<1e-6 # L/chord - 1, not L/chord
        else:
            assert len(r['branches']['junctions_mm'])==1
            from validation.external_methods.evaluate import vmtk_y_groups
            details,metrics=vmtk_y_groups(vessel,r['branches'])
            assert metrics['recovered_branch_count']==3
            assert metrics['branch_connectivity_correct'] and metrics['branch_f1']==1
            assert all(len(d['attached_junction_group_ids'])==1 for d in details)
