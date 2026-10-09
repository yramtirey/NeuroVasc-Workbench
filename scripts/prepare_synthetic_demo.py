"""Prepare a deterministic stylized tube demo using production Lee/EDT analysis."""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import tempfile
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import nibabel as nib
import numpy as np
from neurovasc.io.volume import load_nifti
from neurovasc.report.case_report import analyze_case
from scripts.export_case_geometry import export_geometry
from scripts.export_case_profiles import export_profiles

CASE_ID = 'synthetic_demo'
OWNER = 'neurovasc.synthetic_demo.v1'
NOTE = ('Synthetic demonstration: stylized generated tubes, not real anatomical imaging. '
        'Production Lee skeletonization and EDT caliber; no anatomical or clinical accuracy is claimed.')
SPACING = (0.6, 0.7, 0.8)
SHAPE = (96, 96, 64)


def synthetic_volume():
    affine = np.diag([*SPACING, 1.0])
    affine[:3, 3] = (-28.8, -33.6, -25.6)
    x, y, z = np.indices(SHAPE) * np.array(SPACING)[:, None, None, None]
    data = np.zeros(SHAPE, dtype=np.uint8)
    parameters = []
    # Twelve separate gently curved tubes, in a 4 by 3 arrangement with depth.
    for label in range(1, 13):
        col, row = (label - 1) % 4, (label - 1) // 4
        cx = 8 + col * 13 + 1.5 * np.sin((z - 8) / 7)
        cy = 9 + row * 22 + 1.2 * np.cos((z - 8) / 8 + col)
        radius = 2.6 - 1.15 * np.exp(-((z - (24 + row * 2)) / 4.5)**2)
        mask = ((x-cx)**2 + (y-cy)**2 <= radius**2) & (z >= 8 + col) & (z <= 42 - row)
        if np.any(data[mask]):
            raise ValueError('Synthetic labels overlap')
        data[mask] = label
        parameters.append(dict(label=label, base_x_mm=8+col*13, base_y_mm=9+row*22,
                               z_start_mm=8+col, z_end_mm=42-row, base_radius_mm=2.6,
                               radius_drop_mm=1.15, radius_center_mm=24+row*2, radius_width_mm=4.5,
                               x_wave_amplitude_mm=1.5, x_wave_scale_mm=7,
                               y_wave_amplitude_mm=1.2, y_wave_scale_mm=8, y_phase=col))
    return data, affine, parameters


def build_case(stage):
    data, affine, parameters = synthetic_volume()
    source = stage / 'synthetic_demo.nii.gz'
    nib.save(nib.Nifti1Image(data, affine), source)
    volume = load_nifti(source)
    table, warnings = analyze_case(volume, trim_mm=2.0)
    if warnings or len(table) != 12:
        raise ValueError(f'Incomplete synthetic analysis: {warnings}')
    table['abbreviation'] = [f'SYN-{label:02d}' for label in table.label]
    table['name'] = [f'Synthetic tube {label:02d}' for label in table.label]
    table.to_csv(stage / 'vessel_metrics.csv', index=False)
    provenance = dict(generator=OWNER, synthetic=True, description=NOTE,
                      geometry_parameters=parameters, shape=list(SHAPE), spacing_mm=list(volume.spacing),
                      affine=volume.affine.tolist(), trim_mm=2.0,
                      fingerprint_algorithm='sha256: uint8 C-order voxels + float64 little-endian affine',
                      segmentation_fingerprint=hashlib.sha256(volume.data.tobytes()+volume.affine.astype('<f8').tobytes()).hexdigest(),
                      dependencies={name:version(name) for name in
                                    ('numpy','scipy','nibabel','scikit-image','networkx','pyvista','vtk','pandas')},
                      python=platform.python_version())
    report = dict(case_id=CASE_ID, source_file=source.name, shape=list(volume.shape),
                  spacing_mm=list(volume.spacing), labels_present=list(range(1,13)), vessel_count=12,
                  vessels=table.to_dict(orient='records'), warnings=warnings,
                  measurement_note=NOTE, provenance=provenance)
    (stage/'case_report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    (stage/'provenance.json').write_text(json.dumps(provenance, indent=2, allow_nan=False)+'\n')
    (stage/'geometry').mkdir()
    export_geometry(volume, stage/'geometry', range(1,13), {label:f'Synthetic tube {label:02d}' for label in range(1,13)})
    export_profiles(source, stage, report)


def _prepare_locked(cases_dir, force=False):
    cases_dir = Path(cases_dir)
    cases_dir.mkdir(parents=True, exist_ok=True)
    target = cases_dir/CASE_ID
    if target.exists() or target.is_symlink():
        if target.is_symlink() or not force:
            raise FileExistsError(f'{target} exists; --force only regenerates a generator-owned demo')
        provenance = json.loads((target/'provenance.json').read_text())
        if provenance.get('generator') != OWNER or provenance.get('synthetic') is not True:
            raise ValueError('Refusing to replace a case not owned by this generator')
    # Staging lives outside cases so the API cannot discover a half-built report.
    with tempfile.TemporaryDirectory(prefix='.synthetic-build-', dir=cases_dir.parent) as temporary:
        stage = Path(temporary)/CASE_ID
        stage.mkdir()
        build_case(stage)
        backup = Path(temporary)/'previous'
        if target.exists():
            target.rename(backup)
        try:
            stage.rename(target)
        except BaseException:
            if backup.exists():
                backup.rename(target)
            raise
    return target


def prepare(cases_dir, force=False):
    cases_dir = Path(cases_dir)
    cases_dir.mkdir(parents=True, exist_ok=True)
    lock = cases_dir.parent / '.synthetic-demo.lock'
    # Exclusive creation prevents concurrent generators from replacing each other.
    with lock.open('x'):
        try:
            return _prepare_locked(cases_dir, force)
        finally:
            lock.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()
    print(prepare(Path(__file__).resolve().parents[1]/'outputs'/'cases', args.force))


if __name__ == '__main__':
    main()
