import json
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import pyvista as pv
import pytest
from fastapi.testclient import TestClient

from backend.app import main as api
from scripts import prepare_synthetic_demo as demo


@pytest.fixture(scope='module')
def prepared(tmp_path_factory):
    return demo.prepare(tmp_path_factory.mktemp('demo')/'cases')


def test_deterministic():
    a, affine, parameters = demo.synthetic_volume()
    b, affine_b, parameters_b = demo.synthetic_volume()
    assert np.array_equal(a,b) and np.array_equal(affine,affine_b)
    assert parameters == parameters_b
    assert set(np.unique(a)) == set(range(13))


def test_complete_consistent_physical_artifacts(prepared):
    report = json.loads((prepared/'case_report.json').read_text())
    csv = pd.read_csv(prepared/'vessel_metrics.csv')
    image = nib.load(prepared/'synthetic_demo.nii.gz')
    assert report['vessel_count'] == len(csv) == 12
    assert not report['warnings']
    assert np.allclose(image.header.get_zooms(),demo.SPACING)
    assert np.allclose(image.affine,report['provenance']['affine'])
    digest = hashlib.sha256(np.asarray(image.dataobj).tobytes()+image.affine.astype('<f8').tobytes()).hexdigest()
    assert digest == report['provenance']['segmentation_fingerprint']
    for vessel in report['vessels']:
        assert all(np.isfinite(v) for v in vessel.values() if isinstance(v,(int,float)))
        label = vessel['label']
        profile = json.loads((prepared/'profiles'/f'label_{label}.json').read_text())
        mesh = pv.read(prepared/'geometry'/f'label_{label}.vtp')
        assert mesh.n_points and np.isfinite(mesh.points).all()
        assert len(profile['samples']) > 10
        points = np.array([[s['x_mm'],s['y_mm'],s['z_mm']] for s in profile['samples']])
        indices = np.rint(nib.affines.apply_affine(np.linalg.inv(image.affine),points)).astype(int)
        assert np.all(np.asarray(image.dataobj)[tuple(indices.T)] == label)
        bounds = np.array(mesh.bounds).reshape(3,2)
        assert np.all(points >= bounds[:,0]) and np.all(points <= bounds[:,1])
        for key in ('minimum_diameter_mm','minimum_location_mm','reference_diameter_mm','candidate_reduction_percent'):
            assert profile[key] == pytest.approx(vessel[key])
            assert csv.loc[csv.label==label,key].iloc[0] == pytest.approx(vessel[key])
        minimum = min(profile['samples'], key=lambda s: abs(s['distance_mm']-profile['minimum_location_mm']))
        assert minimum['smooth_diameter_mm'] == pytest.approx(vessel['minimum_diameter_mm'])
        assert [minimum[k] for k in ('x_mm','y_mm','z_mm')] == pytest.approx([vessel[k] for k in ('minimum_x_mm','minimum_y_mm','minimum_z_mm')])
        diameters = [s['raw_diameter_mm'] for s in profile['samples']]
        assert max(diameters)-min(diameters) > 1
        assert all(np.isfinite(list(sample.values())).all() for sample in profile['samples'])
        assert 'Synthetic' in profile['note'] and 'Synthetic tube' in profile['name']


def test_routes_and_legacy(prepared, monkeypatch):
    monkeypatch.setattr(api,'CASES_DIR',prepared.parent)
    legacy = prepared.parent/'legacy_case'
    legacy.mkdir(exist_ok=True)
    report = json.loads((prepared/'case_report.json').read_text())
    report.pop('provenance')
    report['case_id']='legacy_case'
    (legacy/'case_report.json').write_text(json.dumps(report))
    client = TestClient(api.app)
    assert {c['case_id'] for c in client.get('/api/cases').json()} == {'synthetic_demo','legacy_case'}
    assert client.get('/api/cases/legacy_case').status_code == 200
    assert client.get('/api/cases/synthetic_demo').json()['vessel_count'] == 12
    assert client.get('/api/cases/synthetic_demo/geometry/1').status_code == 200
    assert client.get('/api/cases/synthetic_demo/profiles/1').json()['label'] == 1
    for route in ('missing','missing/geometry/1','missing/profiles/1','synthetic_demo/profiles/99'):
        assert client.get('/api/cases/'+route).status_code == 404
    assert client.get('/api/cases').status_code == 200


def test_safe_regeneration_and_failure(prepared, monkeypatch):
    before = (prepared/'case_report.json').read_bytes()
    other = prepared.parent/'other'
    other.mkdir(exist_ok=True)
    (other/'sentinel').write_text('preserve')
    with pytest.raises(FileExistsError):
        demo.prepare(prepared.parent)
    def fail(stage):
        (stage/'case_report.json').write_text('{}')
        raise RuntimeError('simulated failure')
    with monkeypatch.context() as m:
        m.setattr(demo,'build_case',fail)
        with pytest.raises(RuntimeError): demo.prepare(prepared.parent,force=True)
    assert (prepared/'case_report.json').read_bytes() == before
    assert not list(prepared.parent.parent.glob('.synthetic-build-*'))
    demo.prepare(prepared.parent,force=True)
    assert (prepared/'case_report.json').read_bytes() == before
    assert (other/'sentinel').read_text() == 'preserve'
    provenance = prepared/'provenance.json'
    original = provenance.read_text()
    provenance.write_text('{}')
    with pytest.raises(ValueError): demo.prepare(prepared.parent,force=True)
    provenance.write_text(original)


def test_failed_first_build(tmp_path,monkeypatch):
    def fail(stage): raise RuntimeError('failed')
    monkeypatch.setattr(demo,'build_case',fail)
    with pytest.raises(RuntimeError): demo.prepare(tmp_path/'cases')
    assert not (tmp_path/'cases'/'synthetic_demo').exists()
    assert not list(tmp_path.glob('.synthetic-build-*'))


def test_concurrent_generator_is_rejected(tmp_path):
    lock = tmp_path/'.synthetic-demo.lock'
    lock.write_text('another build owns this lock')
    with pytest.raises(FileExistsError): demo.prepare(tmp_path/'cases')
    assert lock.read_text() == 'another build owns this lock'


def test_symlink_target_is_protected(tmp_path):
    cases = tmp_path/'cases'
    cases.mkdir()
    elsewhere = tmp_path/'elsewhere'
    elsewhere.mkdir()
    (cases/'synthetic_demo').symlink_to(elsewhere, target_is_directory=True)
    with pytest.raises(FileExistsError): demo.prepare(cases,force=True)
    assert (cases/'synthetic_demo').is_symlink()
    assert not (tmp_path/'.synthetic-demo.lock').exists()
