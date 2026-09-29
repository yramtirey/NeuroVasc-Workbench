"""Calls actual VMTK; no substitute implementation or truth access."""
from .surface_preparation import snap_seeds

SETTINGS = dict(cost_function='1/R', append_endpoints=False, resampling=False,
                resampling_step_mm=0.5, seed_selector='idlist', surface_smoothing=False,
                vtk_random_seed=0, vtk_smp_threads=1)


def compute(surface, seeds):
    import vtk
    from vmtk import vtkvmtkComputationalGeometryPython as cg
    vtk.vtkMath.RandomSeed(SETTINGS['vtk_random_seed'])
    vtk.vtkSMPTools.Initialize(SETTINGS['vtk_smp_threads'])
    ids, snapped = snap_seeds(surface, seeds)
    source = vtk.vtkIdList(); source.InsertNextId(ids[0])
    targets = vtk.vtkIdList()
    for index in ids[1:]: targets.InsertNextId(index)
    algorithm = cg.vtkvmtkPolyDataCenterlines()
    algorithm.SetInputData(surface)
    algorithm.SetSourceSeedIds(source); algorithm.SetTargetSeedIds(targets)
    algorithm.SetRadiusArrayName('MaximumInscribedSphereRadius')
    algorithm.SetCostFunction(SETTINGS['cost_function'])
    algorithm.SetFlipNormals(0); algorithm.SetSimplifyVoronoi(0)
    algorithm.SetDelaunayTolerance(0.001)
    algorithm.SetAppendEndPointsToCenterlines(0)
    algorithm.SetCenterlineResampling(0)
    algorithm.SetResamplingStepLength(SETTINGS['resampling_step_mm'])
    algorithm.Update()
    data = algorithm.GetOutput()
    if data is None or data.GetNumberOfLines() == 0:
        raise ValueError('vmtk_returned_no_centerlines')
    return data, algorithm.GetVoronoiDiagram(), dict(
        **SETTINGS, requested_seeds_mm=seeds, surface_seed_ids=ids, snapped_seeds_mm=snapped,
        invocation='actual vtkvmtkPolyDataCenterlines compiled filter',
        flip_normals=False, simplify_voronoi=False, delaunay_tolerance=0.001)


def parse(data):
    import numpy as np
    from vtk.util.numpy_support import vtk_to_numpy
    radius = data.GetPointData().GetArray('MaximumInscribedSphereRadius')
    if radius is None:
        raise ValueError('missing_MaximumInscribedSphereRadius')
    points = vtk_to_numpy(data.GetPoints().GetData())
    radii = vtk_to_numpy(radius)
    if len(radii) != len(points) or not np.isfinite(points).all() or not np.isfinite(radii).all() or np.any(radii < 0):
        raise ValueError('invalid_centerline_arrays')
    lines = []
    for i in range(data.GetNumberOfCells()):
        cell = data.GetCell(i)
        if cell.GetCellType() not in (3, 4):
            continue
        ids = [cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())]
        if len(ids) >= 2:
            if np.linalg.norm(np.diff(points[ids],axis=0),axis=1).sum() <= 1e-7:
                raise ValueError('degenerate_centerline_zero_physical_length')
            lines.append(dict(points=points[ids].tolist(), radius_mm=radii[ids].tolist()))
    if not lines:
        raise ValueError('no_polyline_cells')
    return lines
