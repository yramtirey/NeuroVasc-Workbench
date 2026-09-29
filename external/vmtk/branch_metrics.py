"""Actual branch extractor, geometry and bifurcation reference-system outputs."""
from .surface_preparation import write_vtp


def arrays(data, association='point'):
    from vtk.util.numpy_support import vtk_to_numpy
    attributes = data.GetPointData() if association == 'point' else data.GetCellData()
    return {attributes.GetArrayName(i): vtk_to_numpy(attributes.GetArray(i)).tolist()
            for i in range(attributes.GetNumberOfArrays()) if attributes.GetArray(i) is not None}


def compute(centerlines, output):
    from vmtk import vtkvmtkComputationalGeometryPython as cg
    from .centerlines import parse
    extractor = cg.vtkvmtkCenterlineBranchExtractor()
    extractor.SetInputData(centerlines)
    for key in ('Radius', 'Blanking', 'GroupIds', 'CenterlineIds', 'TractIds'):
        getattr(extractor, 'Set'+key+'ArrayName')('MaximumInscribedSphereRadius' if key == 'Radius' else key)
    extractor.Update(); split = extractor.GetOutput()
    write_vtp(split, output / 'branches.vtp')
    geometry = cg.vtkvmtkCenterlineBranchGeometry(); geometry.SetInputData(split)
    for key in ('Radius','GroupIds','Blanking','Length','Curvature','Torsion','Tortuosity'):
        getattr(geometry, 'Set'+key+'ArrayName')('MaximumInscribedSphereRadius' if key == 'Radius' else key)
    geometry.SetLineSmoothing(0); geometry.Update()
    write_vtp(geometry.GetOutput(), output / 'branch_geometry.vtp')
    refs = cg.vtkvmtkCenterlineBifurcationReferenceSystems(); refs.SetInputData(split)
    for key in ('Radius','Blanking','GroupIds','Normal','UpNormal'):
        getattr(refs, 'Set'+key+'ArrayName')('MaximumInscribedSphereRadius' if key == 'Radius' else key)
    refs.Update()
    write_vtp(refs.GetOutput(), output / 'bifurcation_references.vtp')
    cells = arrays(split, 'cell')
    # Export all tracts, including blanked junction tracts. Do not silently
    # merge these into NeuroVasc branches; VMTK GroupIds remain explicit.
    tracts = parse(split)
    for i, tract in enumerate(tracts):
        tract.update({key: value[i] for key, value in cells.items()})
    return dict(tracts=tracts, geometry=arrays(geometry.GetOutput()),
                geometry_cell_arrays=arrays(geometry.GetOutput(), 'cell'),
                junctions_mm=[list(refs.GetOutput().GetPoint(i)) for i in range(refs.GetOutput().GetNumberOfPoints())],
                junction_arrays=arrays(refs.GetOutput()),
                tortuosity_definition='VMTK output L/chord - 1; add 1 only for explicitly labelled conversion')
