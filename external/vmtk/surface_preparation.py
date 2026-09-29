"""VTP interchange in unmodified NIfTI world millimeters (no RAS/LPS flip)."""
from pathlib import Path

PARAMETERS = dict(isovalue=0.5, source='unchanged neurovasc.geometry.mesh.mask_to_mesh',
                  cleaning='production PyVista clean, then VTK clean tolerance 0',
                  triangulation=True, smoothing=False, decimation=False,
                  capping='existing binary phantom end caps retained; no clipping or new caps',
                  open_profiles=False, coordinates='NIfTI affine world mm; no axis flip')


def read_vtp(path):
    import vtk
    path = Path(path)
    if not path.is_file():
        raise ValueError('missing_vtp')
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path)); reader.Update()
    data = reader.GetOutput()
    if reader.GetErrorCode() or data.GetNumberOfPoints() == 0:
        raise ValueError('empty_or_invalid_vtp')
    return data


def write_vtp(data, path):
    import vtk
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path)); writer.SetInputData(data)
    if not writer.Write():
        raise OSError(f'VTP write failed: {path}')


def export_surface(volume, path):
    import vtk
    from neurovasc.geometry.mesh import mask_to_mesh
    mesh = mask_to_mesh(volume).mesh
    clean = vtk.vtkCleanPolyData(); clean.SetInputData(mesh); clean.SetTolerance(0); clean.Update()
    triangles = vtk.vtkTriangleFilter(); triangles.SetInputData(clean.GetOutput()); triangles.Update()
    surface = triangles.GetOutput()
    write_vtp(surface, path)
    return {**PARAMETERS, 'points': surface.GetNumberOfPoints(), 'triangles': surface.GetNumberOfCells(),
            'bounds_mm': list(surface.GetBounds())}


def snap_seeds(surface, points):
    import numpy as np
    import vtk
    xyz = np.asarray(points, float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) < 2 or not np.isfinite(xyz).all():
        raise ValueError('need_two_finite_xyz_seeds')
    locator = vtk.vtkPointLocator(); locator.SetDataSet(surface); locator.BuildLocator()
    ids = [int(locator.FindClosestPoint(p)) for p in xyz]
    if len(set(ids)) != len(ids):
        raise ValueError('seeds_snap_to_same_vertex')
    return ids, [list(surface.GetPoint(i)) for i in ids]
