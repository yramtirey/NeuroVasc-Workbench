"""Isolated actual-VMTK benchmark: export, run, report; production is read-only."""
import argparse
import itertools
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'validation/outputs/external_methods_v1/.matplotlib'))
import nibabel as nib
import numpy as np
from external.vmtk.surface_preparation import export_surface,PARAMETERS
from external.vmtk.cross_sections import status as ce_status
from external.vmtk.centerlines import SETTINGS
from validation.external_methods import evaluate as ev
from validation.external_methods.matching import errors,matched,sample_on_truth
from validation.external_methods.report import ALLOWED,output_directory,json_write,csv_write,aggregates,inventory,sha,clean
from validation.synthetic_centerline.truth import sample_truth


def matrix():
    from validation.synthetic_centerline.geometries import matrix as centerlines
    from validation.synthetic_topology.geometries import matrix as topology
    for gid,experiment in centerlines():
        yield gid,experiment.vessel,dict(cohort=experiment.cohort,phase=experiment.phase_id,
              orientation=experiment.orientation_id,archive='centerline_v1',archive_id=gid)
    for gid,vessel in topology():
        if vessel.geometry in ('simple_ring','ring_branches','figure_eight','communicating_bridge') and vessel.phase==0 and vessel.orientation=='canonical':
            yield 'topology_'+gid,vessel,dict(cohort='loop_exploratory',phase='centered',orientation='canonical',archive='topology_v1',archive_id=gid)


def meta(gid,vessel,extra):
    return dict(geometry_id=gid,geometry=vessel.geometry,spacing=list(map(float,vessel.volume.spacing)),
                spacing_x_mm=float(vessel.volume.spacing[0]),spacing_y_mm=float(vessel.volume.spacing[1]),spacing_z_mm=float(vessel.volume.spacing[2]),
                comparison_scope=ev.scope(vessel.geometry),true_diameter_mm=(vessel.proximal_diameter_mm if hasattr(vessel,'proximal_diameter_mm') and vessel.proximal_diameter_mm==vessel.distal_diameter_mm else None),**extra)


def automatic_seeds(recovery):
    if recovery.graph is None: return [],recovery.error or 'automatic_endpoint_recovery_failed'
    points=[list(map(float,d['world'])) for n,d in recovery.graph.nodes(data=True) if recovery.graph.degree[n]==1]
    points.sort(key=lambda p:(p[2],p[0],p[1]))
    return (points,'') if len(points)>=2 else (points,'fewer_than_two_segmentation_derived_endpoints')


def oracle_seeds(vessel):
    if vessel.geometry=='bifurcation':
        points=[s[-1].tolist() for s in vessel.branches.values()]
    else: points=sample_truth(vessel,samples=2).tolist()
    return sorted(points,key=lambda p:(p[2],p[0],p[1]))


def probe(python):
    if not Path(python).is_file(): return dict(available=False,reason='external_python_not_found',executable=str(python))
    try:
        p=subprocess.run([str(python),'-m','external.vmtk.adapter','--probe'],cwd=ROOT,capture_output=True,text=True,timeout=120)
        if p.returncode: return dict(available=False,reason='probe_failed',stderr=p.stderr,returncode=p.returncode)
        return json.loads(p.stdout)
    except (OSError,subprocess.TimeoutExpired,json.JSONDecodeError) as e:
        return dict(available=False,reason=f'{type(e).__name__}: {e}')


def clear_job_outputs(destination,include_result=False):
    """Remove only this worker's known generated files, never its input surface.

    A rejected or interrupted rerun must not leave old VTPs looking like current
    successful output, even though the evaluator already follows result.json.
    """
    names=['centerlines.vtp','voronoi.vtp','branches.vtp','branch_geometry.vtp','bifurcation_references.vtp']
    if include_result:names.append('result.json')
    for name in names:(destination/name).unlink(missing_ok=True)


def export(output,limit=None):
    from neurovasc.geometry.main_path import extract_main_path
    from neurovasc.graph.vascular_graph import VascularGraph
    from neurovasc.geometry.diameter_profile import build_diameter_profile
    from neurovasc.geometry.cross_sectional_caliber import measure_cross_sectional_caliber
    from neurovasc.geometry.topology_centerline import extract_topology_network
    from neurovasc.geometry.centerline_qc import Recovery
    from neurovasc.geometry.centerline_refinement import principal_points
    rows=[];caliber={};supports={};branches=[];geometries=[];requests=[];examples={};surface_parameters={}
    chosen=list(matrix())[:limit]
    for gid,vessel,extra in chosen:
        info=meta(gid,vessel,extra);geometries.append({**info,**clean(vessel.metadata()),'affine':vessel.volume.affine.tolist()})
        archive=ROOT/'validation/outputs'/extra['archive']/'masks'/f"{extra['archive_id']}.nii.gz"
        archived=nib.load(archive)
        if not np.array_equal(np.asarray(archived.dataobj),vessel.volume.data) or not np.allclose(archived.affine,vessel.volume.affine,atol=1e-6):
            raise ValueError(f'archive_mask_affine_mismatch:{gid}')
        geometries[-1]['archived_mask_sha256']=sha(archive)
        surface=output/'inputs'/gid/'surface.vtp';surface_parameters[gid]=export_surface(vessel.volume,surface)
        if info['comparison_scope']==ev.MISMATCH:
            recoveries={};auto=[];auto_error='cyclic_network_outside_declared_tree_comparison'
        else:
            recoveries=ev.neurovasc_recoveries(vessel.volume)
            auto,auto_error=automatic_seeds(recoveries['geodesic_raw'])
            if vessel.geometry=='bifurcation':
                try:
                    network=extract_topology_network(vessel.volume)
                    recoveries['topology']=Recovery(network,principal_points(network))
                except ValueError as e: recoveries['topology']=Recovery(error=str(e))
            for method,recovery in recoveries.items():
                row={**info,'method':method,'oracle_seed':False,'extraction_success':False,'status':'failed','failure_reason':recovery.error}
                if recovery.graph is not None and recovery.path is not None:
                    try:
                        paths=ev.network_paths(recovery.graph) if vessel.geometry=='bifurcation' else [recovery.path]
                        metrics,support=ev.centerline(vessel,recovery.graph,paths)
                        row.update(metrics,extraction_success=True,status='success',failure_reason='')
                        supports[f'{gid}|{method}|False']=clean(support)
                        if vessel.geometry=='bifurcation':
                            bm,details=ev.y_branches(vessel,recovery.graph)
                            branches.append({**row,**bm,'branch_definition':'NeuroVasc logical reduction','details':details})
                            examples[f'{gid}|{method}']={'paths':[p.tolist() for p in paths]}
                    except ValueError as e: row['failure_reason']=str(e)
                rows.append(row)
            if vessel.geometry!='bifurcation':
                stations=sample_truth(vessel,samples=121)[8:-8] # fixed analytic 2 mm cap exclusion; .25 mm stations
                truth=vessel.true_diameter(stations)
                item=dict(stations=stations.tolist(),truth_diameter=truth.tolist(),methods={})
                try:
                    if recoveries['lee'].graph is None: raise ValueError(recoveries['lee'].error)
                    path=extract_main_path(VascularGraph(recoveries['lee'].graph))
                    profile=build_diameter_profile(vessel.volume,path,trim_mm=0)
                    section=measure_cross_sectional_caliber(vessel.volume,path)
                    for method,value in [('edt',profile.diameter_mm),('neurovasc_ce',section.diameter_mm)]:
                        estimate=sample_on_truth(vessel,path.world_points,value,stations,max(vessel.volume.spacing))
                        item['methods'][method]={'values':clean(estimate),'status':'success' if np.isfinite(estimate).any() else 'failed','failure_reason':''}
                    item['section_reasons']=section.failure_reasons
                except ValueError as e:
                    for method in ('edt','neurovasc_ce'):item['methods'][method]={'values':[None]*len(stations),'status':'failed','failure_reason':str(e)}
                caliber[gid]=item
        for oracle in (True,False):
            mismatch=info['comparison_scope']==ev.MISMATCH
            request={**info,'oracle_seed':oracle,'surface':str(surface.resolve()),
                     'seeds_mm':[] if mismatch else oracle_seeds(vessel) if oracle else auto,
                     'seed_strategy':'analytic_endpoints_diagnostic' if oracle else 'unchanged_raw_geodesic_mask_endpoints_lexicographic_zxy',
                     'seed_failure':'cyclic_network_outside_declared_tree_comparison' if mismatch else '' if oracle else auto_error}
            path=surface.parent/f'request_{"oracle" if oracle else "automatic"}.json'
            json_write(path,request); requests.append(str(path.relative_to(output)))
        print('export',gid,flush=True)
    json_write(output/'matrix.json',geometries);json_write(output/'requests.json',requests)
    json_write(output/'neurovasc_results.json',rows);json_write(output/'neurovasc_branches.json',branches)
    json_write(output/'caliber_samples.json',caliber);json_write(output/'centerline_support.json',supports)
    json_write(output/'examples.json',examples);json_write(output/'surface_parameters.json',dict(common=PARAMETERS,by_geometry=surface_parameters,centerline_settings=SETTINGS))


def execute(output,python,reuse=False):
    env=probe(python);json_write(output/'environment.json',dict(external=env,slicervmtk_ce=ce_status()))
    def job(relative):
        path=output/relative;request=json.loads(path.read_text());mode='oracle' if request['oracle_seed'] else 'automatic'
        dest=output/'vmtk_outputs'/request['geometry_id']/mode;dest.mkdir(parents=True,exist_ok=True)
        result_path=dest/'result.json'
        if request['comparison_scope']==ev.MISMATCH:
            result=dict(status='assumption_mismatch',extraction_success=None,failure_reason=request['seed_failure'])
        elif not env['available']:
            result=dict(status='unavailable',extraction_success=None,failure_reason=env.get('reason','external_runtime_unavailable'))
        else:
            if reuse and result_path.exists():
                prior=json.loads(result_path.read_text())
                if prior.get('request_sha256')==sha(path) and prior.get('surface_sha256')==sha(request['surface']) and prior.get('environment',{}).get('module_sha256')==env.get('module_sha256') and prior.get('environment',{}).get('adapter_sha256')==env.get('adapter_sha256'):
                    if not prior.get('extraction_success'):clear_job_outputs(dest)
                    return
            try:
                clear_job_outputs(dest,include_result=True)
                proc=subprocess.run([str(python),'-m','external.vmtk.adapter','--request',str(path),'--output',str(dest)],cwd=ROOT,capture_output=True,text=True,timeout=120)
                (dest/'stdout.log').write_text(proc.stdout);(dest/'stderr.log').write_text(proc.stderr)
                if proc.returncode or not result_path.exists():
                    result=dict(status='failed',extraction_success=False,failure_reason=f'external_process_exit_{proc.returncode}')
                else:
                    result=json.loads(result_path.read_text())
                    if result.get('request_sha256')!=sha(path) or result.get('surface_sha256')!=sha(request['surface']):
                        raise ValueError('stale_or_mismatched_external_result')
            except (OSError,ValueError,subprocess.TimeoutExpired) as e:
                result=dict(status='failed',extraction_success=False,failure_reason=f'{type(e).__name__}: {e}')
        if not result.get('extraction_success'):clear_job_outputs(dest)
        result.update(geometry_id=request['geometry_id'],oracle_seed=request['oracle_seed'],comparison_scope=request['comparison_scope'])
        json_write(result_path,result)
        print('external',request['geometry_id'],mode,result['status'],flush=True)

    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(job,json.loads((output/'requests.json').read_text())))


def report(output):
    rows=json.loads((output/'neurovasc_results.json').read_text());branch_rows=json.loads((output/'neurovasc_branches.json').read_text())
    caliber=json.loads((output/'caliber_samples.json').read_text());supports=json.loads((output/'centerline_support.json').read_text())
    examples=json.loads((output/'examples.json').read_text());geometries=json.loads((output/'matrix.json').read_text())
    lookup={gid:(vessel,meta(gid,vessel,extra)) for gid,vessel,extra in matrix() if gid in {g['geometry_id'] for g in geometries}}
    for relative in json.loads((output/'requests.json').read_text()):
        request=json.loads((output/relative).read_text());gid=request['geometry_id'];oracle=request['oracle_seed'];mode='oracle' if oracle else 'automatic'
        path=output/'vmtk_outputs'/gid/mode/'result.json'
        result=json.loads(path.read_text()) if path.exists() else dict(status='unavailable',extraction_success=None,failure_reason='external_result_missing')
        vessel,info=lookup[gid];row={**info,'method':'vmtk','oracle_seed':oracle,**{k:result.get(k) for k in ('status','extraction_success','failure_reason')}}
        if result.get('extraction_success'):
            # Only import results bound to current input bytes, never arbitrary stale VTPs.
            if result.get('surface_sha256')!=sha(request['surface']) or result.get('request_sha256')!=sha(output/relative):
                raise ValueError(f'external_result_provenance_mismatch:{gid}')
            graph=ev.graph_from_lines(result['lines']);paths=[np.array(line['points']) for line in result['lines']]
            metrics,support=ev.centerline(vessel,graph,paths);row.update(metrics)
            supports[f'{gid}|vmtk|{oracle}']=clean(support)
            if vessel.geometry=='bifurcation':
                bm,details=ev.y_branches(vessel,graph)
                b=result.get('branches',{});geometry=b.get('geometry',{})
                refs=np.asarray(b.get('junctions_mm',[]),float)
                group_details,group_metrics=ev.vmtk_y_groups(vessel,b)
                branch_rows.append({**row,**bm,'duplicate_path_embedding_branch_count':bm['recovered_branch_count'],'duplicate_path_embedding_f1':bm['branch_f1'],**group_metrics,'vmtk_group_details':group_details,'branch_definition':'actual VMTK nonblanked GroupIds attached through blanked junction groups; partial length comparison',
                    'vmtk_branch_group_count':len(geometry.get('GroupIds',[])),
                    'vmtk_reference_junction_count':len(refs),
                    'vmtk_reference_junction_error_mm':float(np.linalg.norm(refs,axis=1).mean()) if len(refs) else None,
                    'vmtk_group_geometry':geometry,'details':[{k:v for k,v in d.items() if k not in ('radius_definition','mean_radius_mm','true_radius_mm')} for d in details],'branch_failure_reason':result.get('branch_failure_reason',group_metrics.get('branch_failure_reason',''))})
                examples[f'{gid}|vmtk_{mode}']={'paths':[p.tolist() for p in paths],'junctions_mm':b.get('junctions_mm',[])}
            elif gid in caliber:
                line=max(result['lines'],key=lambda l:sum(np.linalg.norm(np.diff(l['points'],axis=0),axis=1)))
                item=caliber[gid]
                values=sample_on_truth(vessel,line['points'],2*np.asarray(line['radius_mm']),item['stations'],max(vessel.volume.spacing))
                item['methods'][f'vmtk_{mode}']={'values':clean(values),'status':'success' if np.isfinite(values).any() else 'failed',
                    'failure_reason':'' if np.isfinite(values).any() else 'no_valid_gated_analytic_station_support'}
        if gid in caliber and f'vmtk_{mode}' not in caliber[gid]['methods']:
            caliber[gid]['methods'][f'vmtk_{mode}']={'values':[None]*len(caliber[gid]['stations']),'status':row['status'],'failure_reason':row['failure_reason']}
        rows.append(row)
    caliber_rows=[];comparisons=[];sample_rows=[]
    for gid,item in caliber.items():
        vessel,info=lookup[gid];truth=np.array(item['truth_diameter'])
        for method,entry in item['methods'].items():
            oracle=method=='vmtk_oracle';values=np.array(entry['values'],dtype=float)
            for quantity,scale in [('diameter',1.),('radius',.5)]:
                if quantity=='radius' and method=='neurovasc_ce': continue
                caliber_rows.append({**info,'method':method,'oracle_seed':oracle,'quantity':quantity,**errors(values*scale,truth*scale),
                                     'status':entry['status'],'failure_reason':entry['failure_reason'],'requested_sample_count':len(values),
                                     'valid_fraction':float(np.isfinite(values).mean()),'support':'analytic .25 mm stations, 2 mm cap exclusion'})
            sample_rows.extend({**info,'method':method,'oracle_seed':oracle,'sample_index':i,'x_mm':p[0],'y_mm':p[1],'z_mm':p[2],
                               'true_diameter_mm':t,'estimated_diameter_mm':v} for i,(p,t,v) in enumerate(zip(item['stations'],truth,values)))
        for left,right in itertools.combinations(item['methods'],2):
            a=np.array(item['methods'][left]['values'],float);b=np.array(item['methods'][right]['values'],float)
            comparisons.append({**info,'metric_family':'caliber','left_method':left,'right_method':right,
                                'left_oracle_seed':left=='vmtk_oracle','right_oracle_seed':right=='vmtk_oracle',**matched(a,b,truth)})
        caliber_rows.append({**info,'method':'slicervmtk_ce','oracle_seed':None,'quantity':'diameter','status':'unavailable',
                            'failure_reason':ce_status()['reason'],'sample_count':None,'mae_mm':None,'estimated_mean_mm':None})
    for gid,(_,info) in lookup.items():
        if info['comparison_scope']==ev.MISMATCH: continue
        candidates=[r for r in rows if r['geometry_id']==gid]
        for left_row,right_row in itertools.combinations(candidates,2):
            left=f"{left_row['method']}|{left_row['oracle_seed']}";right=f"{right_row['method']}|{right_row['oracle_seed']}"
            a=supports.get(gid+'|'+left);b=supports.get(gid+'|'+right)
            common=dict(**info,metric_family='centerline_truth_support',left_method=left,right_method=right,
                        left_status=left_row['status'],right_status=right_row['status'])
            if a is None or b is None:
                comparisons.append({**common,'status':'unscorable','common_sample_count':0,
                                    'requested_sample_count':len((a or b or {}).get('covered',[])) or None,
                                    'common_support_fraction':None,'left_mean_truth_to_path_mm':None,'right_mean_truth_to_path_mm':None})
                continue
            cover=np.array(a['covered'])&np.array(b['covered'])
            comparisons.append({**common,'status':'matched' if cover.any() else 'unscorable',
                'common_sample_count':int(cover.sum()),'requested_sample_count':len(cover),'common_support_fraction':float(cover.mean()),
                'left_coverage':float(np.mean(a['covered'])),'right_coverage':float(np.mean(b['covered'])),
                'left_mean_truth_to_path_mm':float(np.array(a['nearest_error'])[cover].mean()) if cover.any() else None,
                'right_mean_truth_to_path_mm':float(np.array(b['nearest_error'])[cover].mean()) if cover.any() else None})
    failures=[r for r in rows+caliber_rows if r.get('status')!='success']
    failures.extend({**r,'status':'branch_failed','failure_reason':r['branch_failure_reason']} for r in branch_rows if r.get('branch_failure_reason'))
    csv_write(output/'centerline_results.csv',rows,('geometry_id','method','oracle_seed','comparison_scope'))
    csv_write(output/'results.csv',[dict(metric_family='centerline',**r) for r in rows]+[dict(metric_family='caliber',**r) for r in caliber_rows]+[dict(metric_family='branch',**r) for r in branch_rows])
    csv_write(output/'caliber_results.csv',caliber_rows);csv_write(output/'branch_results.csv',branch_rows,('geometry_id','method','comparison_scope'))
    csv_write(output/'matched_comparison.csv',comparisons,('geometry_id','metric_family','left_method','right_method','common_sample_count'))
    csv_write(output/'failures.csv',failures,('geometry_id','method','status','failure_reason'))
    csv_write(output/'caliber_profiles.csv',sample_rows)
    env=json.loads((output/'environment.json').read_text()) if (output/'environment.json').exists() else {'external':{'available':False}}
    env['neurovasc']=dict(python=sys.version,executable=sys.executable,
        versions={name:importlib.metadata.version(name) for name in
                  ('numpy','scipy','scikit-image','nibabel','networkx','pyvista','vtk','matplotlib')})
    env['installation_records']={str(p.relative_to(ROOT)):sha(p) for p in
        (ROOT/'external/vmtk/install-report.json',ROOT/'external/vmtk/requirements.lock') if p.exists()}
    json_write(output/'environment.json',env)
    summary=dict(status='executed' if any(r['method']=='vmtk' and r.get('extraction_success') for r in rows) else 'External benchmark infrastructure complete; VMTK execution pending',
                 geometry_count=len(geometries),scope_counts={s:sum(g['comparison_scope']==s for g in geometries) for s in (ev.DIRECT,ev.PARTIAL,ev.MISMATCH)},
                 centerline=aggregates(rows),slicervmtk_ce=ce_status(),production_recommendation='No production default change',
                 external_environment=env['external'],source_sha256={str(p.relative_to(ROOT)):sha(p) for folder in ('external/vmtk','validation/external_methods') for p in (ROOT/folder).glob('*.py')},
                 conventions=dict(ground_truth='existing analytic phantoms, not VMTK',oracle_and_automatic_separate=True,
                                  mismatch_excluded=True,caliber_support='analytic 0.25 mm stations with 2 mm cap exclusion; common finite gated support',
                                  radius_and_area_estimands_distinct=True,automatic_seed_dependency='unchanged NeuroVasc raw geodesic endpoints'))
    def caliber_aggregate(subset):
        out={}
        for row in subset:
            key=row['method']+'|'+row['quantity'];out.setdefault(key,[]).append(row)
        result={}
        for key,group in out.items():
            good=[r for r in group if r.get('sample_count') and r.get('mae_mm') is not None]
            result[key]=dict(requested_runs=len(group),scorable_runs=len(good),
                unavailable_runs=sum(r['status']=='unavailable' for r in group),
                mean_valid_fraction=float(np.mean([r['valid_fraction'] for r in group if 'valid_fraction' in r])) if any('valid_fraction' in r for r in group) else None,
                **{k:float(np.mean([r[k] for r in good])) if good else None for k in ('estimated_mean_mm','mae_mm','bias_mm','percent_error')},
                rmse_mm=float(np.sqrt(np.mean([r['rmse_mm']**2 for r in good]))) if good else None)
        return result
    summary['caliber_own_support']=caliber_aggregate(caliber_rows)
    summary['by_spacing']={str(z):dict(centerline=aggregates([r for r in rows if r['spacing_z_mm']==z]),
                                      caliber=caliber_aggregate([r for r in caliber_rows if r['spacing_z_mm']==z])) for z in (.5,.8)}
    summary['by_orientation']={o:dict(centerline=aggregates([r for r in rows if r['geometry']=='oblique' and r['orientation']==o]),
        caliber=caliber_aggregate([r for r in caliber_rows if r['geometry']=='oblique' and r['orientation']==o])) for o in ('111','121','213')}
    summary['one_mm_cases']=[r for r in caliber_rows if r.get('true_diameter_mm')==1. and r['quantity']=='diameter']
    summary['y_bifurcation']=branch_rows
    paired={}
    for r in comparisons:
        if r['metric_family']!='caliber':continue
        key=r['left_method']+' vs '+r['right_method'];paired.setdefault(key,[]).append(r)
    summary['caliber_common_support']={k:dict(requested_pairs=len(rs),scorable_pairs=sum(r['status']=='matched' for r in rs),
        mean_common_fraction=float(np.mean([r['common_support_fraction'] for r in rs])),
        **{field:float(np.mean([r[field] for r in rs if r[field] is not None])) if any(r[field] is not None for r in rs) else None
           for field in ('left_mae_mm','right_mae_mm','right_minus_left_mae_mm')}) for k,rs in paired.items()}
    summary['source_sha256']['scripts/run_external_methods_benchmark.py']=sha(Path(__file__))
    summary['neurovasc_environment']=env['neurovasc']
    summary['installation_records']=env['installation_records']
    summary['production_source_sha256']={str(p.relative_to(ROOT)):sha(p) for folder in
        ('neurovasc/geometry','neurovasc/graph') for p in (ROOT/folder).glob('*.py')}
    parameters=json.loads((output/'surface_parameters.json').read_text())
    parameters['centerline_settings']=SETTINGS
    json_write(output/'surface_parameters.json',parameters)
    json_write(output/'summary.json',summary);json_write(output/'examples.json',examples)
    from validation.external_methods.figures import figures
    figures(output,rows,caliber_rows,branch_rows,examples,lookup)
    inventory(output)
    print(json.dumps({'status':summary['status'],'geometry_count':summary['geometry_count'],'scope_counts':summary['scope_counts']}),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ALLOWED)
    p.add_argument('--external-python',type=Path,default=ROOT/'external/vmtk/.venv/bin/python')
    p.add_argument('--stage',choices=('export','run','report','all'),default='all')
    p.add_argument('--reuse',action='store_true');p.add_argument('--limit',type=int)
    a=p.parse_args();output=output_directory(a.output)
    if a.limit is not None and (a.limit<1 or output==ALLOWED): p.error('--limit requires a positive count and a separate descendant output')
    output.mkdir(parents=True,exist_ok=True)
    if a.stage in ('export','all'): export(output,a.limit)
    if a.stage in ('run','all'): execute(output,a.external_python.absolute(),a.reuse)
    if a.stage in ('report','all'): report(output)


if __name__=='__main__':main()
