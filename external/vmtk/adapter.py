"""Optional VMTK subprocess worker. JSON/VTP boundary isolates incompatible VTK."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def environment():
    result = dict(python=sys.version, executable=sys.executable, platform=platform.platform(),
                  architecture=platform.machine(), adapter_sha256={p.name:sha256(p) for p in Path(__file__).parent.glob('*.py')}, installation_source='PyPI; isolated pip environment',
                  available=False, versions={}, module_sha256={})
    try:
        import vtk
        import vmtk
        from vmtk import vtkvmtkComputationalGeometryPython as cg
        for package in ('vmtk','vtk','numpy'):
            result['versions'][package] = importlib.metadata.version(package)
        package_dir = Path(vmtk.__file__).parent
        for name in ('vmtkcenterlines.py','vmtkbranchextractor.py','vmtkbranchgeometry.py',
                     'vmtkbifurcationreferencesystems.py',Path(cg.__file__).name,'libvtkvmtkComputationalGeometry.dylib'):
            path = package_dir/name
            if path.exists(): result['module_sha256'][name] = sha256(path)
        cg.vtkvmtkPolyDataCenterlines()
        result.update(available=True, vtk_runtime=vtk.vtkVersion.GetVTKVersion(),
                      invocation='direct compiled VMTK computational geometry API')
        try:
            from vmtk import vtkvmtk
            result['umbrella_import_available'] = True
        except (ImportError,OSError) as error:
            result['umbrella_import_available'] = False
            result['umbrella_import_error'] = str(error)
    except (ImportError, OSError, importlib.metadata.PackageNotFoundError) as error:
        result['reason'] = f'{type(error).__name__}: {error}'
    return result


def run(request, output):
    from .surface_preparation import read_vtp, write_vtp
    from .centerlines import compute, parse
    from .branch_metrics import compute as branches
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    env = environment()
    result = dict(schema_version=1, geometry_id=request['geometry_id'], oracle_seed=request['oracle_seed'],
                  comparison_scope=request['comparison_scope'], extraction_success=False,
                  status='unavailable', failure_reason=env.get('reason',''), environment=env,
                  surface_sha256=sha256(request['surface']), request_sha256=sha256(request['request_file']))
    if env['available']:
        try:
            if request.get('seed_failure'):
                raise ValueError(request['seed_failure'])
            surface = read_vtp(request['surface'])
            centerline, voronoi, settings = compute(surface, request['seeds_mm'])
            lines = parse(centerline)
            write_vtp(centerline, output/'centerlines.vtp')
            if voronoi is not None:
                write_vtp(voronoi, output/'voronoi.vtp')
            result.update(extraction_success=True, status='success', failure_reason='', lines=lines, settings=settings)
            if request['comparison_scope'] != 'METHOD_ASSUMPTION_MISMATCH':
                try:
                    result['branches'] = branches(centerline, output)
                except Exception as error:
                    result['branch_failure_reason'] = f'{type(error).__name__}: {error}'
        except Exception as error:
            result.update(status='failed', failure_reason=f'{type(error).__name__}: {error}')
    (output/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', action='store_true')
    parser.add_argument('--request', type=Path)
    parser.add_argument('--output', type=Path)
    args=parser.parse_args()
    if args.probe:
        print(json.dumps(environment(), allow_nan=False)); return
    if args.request is None or args.output is None:
        parser.error('--request and --output required')
    # Worker is intentionally restricted to this milestone as well as the runner.
    root=Path(__file__).resolve().parents[2]
    allowed=(root/'validation/outputs/external_methods_v1').resolve()
    if not args.output.resolve().is_relative_to(allowed):
        parser.error('output outside external_methods_v1')
    request=json.loads(args.request.read_text()); request['request_file']=str(args.request.resolve())
    run(request,args.output)


if __name__ == '__main__':
    main()
