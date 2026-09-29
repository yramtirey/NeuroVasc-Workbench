"""Isolated topology-preserving refinement validation; no production mutation."""
import argparse
from dataclasses import asdict, replace
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'validation/outputs/.matplotlib'))
import matplotlib
matplotlib.use('Agg')
import networkx as nx
import nibabel as nib
import numpy as np
from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.centerline_refinement import Refiner, resample, distance
from neurovasc.geometry.centerline_qc import smoothness
from neurovasc.geometry.topology_centerline import extract_topology_network
from neurovasc.geometry.topology_refinement import CANDIDATE, refine_network, clean_path
from neurovasc.graph.topology_network import reduce_network
from neurovasc.graph.vascular_graph import build_vascular_graph
from validation.synthetic_topology.experiments_v2 import matrix_v2
from validation.synthetic_topology.evaluate_v2 import METHODS, measure, summarize
from validation.synthetic_topology.figures_v2 import figures, overlay
from scripts.run_caliber_validation import write_csv


def output_directory(path):
    path=Path(path).resolve();allowed=(ROOT/'validation/outputs/topology_v2').resolve()
    if not path.is_relative_to(allowed):raise ValueError('Output must be topology_v2 or a descendant; archives protected')
    return path


def serial(value):
    if isinstance(value,dict):return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple,set)):return [serial(v) for v in value]
    if isinstance(value,np.ndarray):return serial(value.tolist())
    if isinstance(value,np.generic):return serial(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def configurations():
    # Controlled one-factor-at-a-time sweep, not a per-geometry optimizer.
    return [('candidate',CANDIDATE)]+[(name,replace(CANDIDATE,**changes)) for name,changes in (
        ('radius_0.5',dict(junction_cluster_radius_mm=.5)),('radius_1.5',dict(junction_cluster_radius_mm=1.5)),
        ('resample_1.0',dict(resample_spacing_mm=1.)),('smooth_0',dict(smoothing_scale_mm=0.)),
        ('smooth_2',dict(smoothing_scale_mm=2.)),('spur_off',dict(spur_handling='off')))]


def in_sweep(p,cohort):
    return p.orientation=='rotated' or (cohort=='core_v1' and (p.geometry in ('figure_eight','ring_branches') or p.geometry=='small_bridge' and p.bridge_diameter_mm==1. or p.geometry=='near_touch' and p.minimum_separation_mm==.4)) or (p.geometry=='double_junction' and p.phase==0)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=ROOT/'validation/outputs/topology_v2')
    args=parser.parse_args()
    try:output=output_directory(args.output)
    except ValueError as error:parser.error(str(error))
    output.mkdir(parents=True,exist_ok=True)
    for folder in ('masks','networks'): (output/folder).mkdir(exist_ok=True)
    tables={name:[] for name in ('results','junction_results','branch_results','cycle_results','connectivity_results','spur_results','length_results','matched_comparison','refinement_sweep','endpoint_sweep','qc_results','invariant_results','junction_decisions')}
    metadata=[]
    for eid,p,cohort in matrix_v2(ROOT/'validation/outputs/topology_v1'):
        key=dict(geometry_id=eid,geometry=p.geometry,cohort=cohort,phase=p.phase,orientation=p.orientation,spacing_z_mm=float(p.volume.spacing[2]),bridge_diameter_mm=p.bridge_diameter_mm,minimum_separation_mm=p.minimum_separation_mm)
        metadata.append({**key,**p.metadata()});nib.save(nib.Nifti1Image(p.volume.data,p.volume.affine),output/'masks'/f'{eid}.nii.gz')
        networks={};errors={}
        for method,extractor in [('lee',lambda:build_vascular_graph(extract_centerline(p.volume)).graph),('topology_v1',lambda:extract_topology_network(p.volume))]:
            try:networks[method]=reduce_network(extractor())
            except ValueError as error:networks[method]=None;errors[method]=str(error)
        original=networks['topology_v1'];refiner=Refiner(p.volume);refinements={}
        if original is not None:
            for method,stage in [('topology_v2_junction','junction'),('topology_v2_refined','full')]:
                result=refine_network(p.volume,original,stage=stage,refiner=refiner);refinements[method]=result;networks[method]=result.network
                tables['invariant_results'].append({**key,'method':method,**result.invariants,'accepted':result.accepted,'failure_reason':result.failure_reason})
                for table,data in [('spur_results',result.spur_diagnostics),('junction_decisions',result.junction_diagnostics)]:tables[table].extend({**key,'method':method,**r} for r in data)
        else:
            for m in METHODS[2:]:networks[m]=None;errors[m]='upstream_extraction_failed'
        group={}
        for method in METHODS:
            result=refinements.get(method);network=networks[method]
            metrics,branches,cycles,junctions,qc=measure(p,network,errors.get(method,''))
            if result is not None:
                qc.extend(dict(criterion=k,value=bool(v),threshold='true',passed=bool(v),critical=True)
                          for k,v in result.invariants.items() if isinstance(v,(bool,np.bool_)))
                qc.append(dict(criterion='dense_lumen_and_component_contact',value=result.accepted,
                               threshold='occupancy >=0.5 at min_spacing/8; component separation > sampling step',
                               passed=result.accepted,critical=True))
            row={**key,'method':method,**metrics,'refinement_accepted':result.accepted if result else method not in METHODS[2:],
                 'refinement_failure_reason':result.failure_reason if result else '',
                 'output_state':'refined' if result and result.accepted else 'original_retained_after_rejection' if result else 'baseline'}
            group[method]=row;tables['results'].append(row)
            for table,data in [('branch_results',branches),('cycle_results',cycles),('junction_results',junctions),('qc_results',qc)]:tables[table].extend({**key,'method':method,**r} for r in data)
            tables['connectivity_results'].append({**key,'method':method,**{k:v for k,v in row.items() if k.startswith(('adjacency_','branch_')) or k in ('topology_correct','branch_connectivity_correct','output_state')}})
            if network is not None:
                raw=sum(d['length_mm'] for _,_,d in network.raw.edges(data=True))
                source=original if method in METHODS[2:] else network
                raw_source=sum(d['length_mm'] for _,_,d in source.raw.edges(data=True))
                sampled=sum(distance(resample(clean_path(b['points']),CANDIDATE.resample_spacing_mm))[-1] for b in source.branches.values())
                for estimate,length in [('raw_voxel_graph',raw_source),('resampled_branch_polylines',sampled),('returned_network',raw)]:
                    truth=row['true_network_length_mm'];tables['length_results'].append({**key,'method':method,'estimate':estimate,'length_mm':length,'true_length_mm':truth,'signed_bias_mm':length-truth,'absolute_error_mm':abs(length-truth),'percent_error':100*abs(length-truth)/truth,'output_state':row['output_state']})
                payload=dict(nodes=dict(network.logical.nodes(data=True)),branches=network.branches,cycles=network.cycles,counts=network.counts(),metadata=network.metadata,output_state=row['output_state'])
                (output/'networks'/f'{eid}_{method}.json').write_text(json.dumps(serial(payload),indent=2,allow_nan=False)+'\n')
        if original is not None:
            new=refinements['topology_v2_refined']
            if not group['topology_v2_refined']['topology_correct'] or not new.accepted or p.phase==0:
                overlay(output,eid,p,original,new.network,new.accepted)
            if in_sweep(p,cohort):
                for name,params in configurations():
                    result=new if name=='candidate' else refine_network(p.volume,original,params,refiner=refiner)
                    metrics,*_=measure(p,result.network)
                    tables['refinement_sweep'].append({**key,'configuration':name,**asdict(params),**metrics,'refinement_accepted':result.accepted,'failure_reason':result.failure_reason,'topology_invariant_pass':result.invariants['topology_invariant_pass']})
            if p.geometry in ('y','incomplete_loop','near_junction_branch'):
                for mode in ('raw','recentered','tangent','boundary'):
                    result=refine_network(p.volume,original,replace(CANDIDATE,endpoint_mode=mode),refiner=refiner)
                    metrics,*_=measure(p,result.network)
                    turns=[smoothness(b['points'])['tangent_step_max_deg'] for b in result.network.branches.values()]
                    tables['endpoint_sweep'].append({**key,'endpoint_mode':mode,**metrics,'refinement_accepted':result.accepted,'failure_reason':result.failure_reason,'max_tangent_step_deg':max((t for t in turns if t is not None),default=None)})
        for baseline in METHODS[:-1]:
            old,new=group[baseline],group['topology_v2_refined'];pair={**key,'baseline':baseline}
            for metric in ('cycle_rank_correct','branch_connectivity_correct','branch_f1','adjacency_f1','network_length_error_mm','junction_position_error_mm','endpoint_position_error_mm','network_coverage_fraction','mean_network_position_error_mm'):
                pair[metric+'_delta']=float(new[metric])-float(old[metric]) if new.get(metric) is not None and old.get(metric) is not None else None
            tables['matched_comparison'].append(pair)
        print(eid,group['topology_v2_refined']['refinement_accepted'],group['topology_v2_refined']['network_length_error_mm'],flush=True)
    rows=tables['results']
    def by_method(selected):return {m:summarize([r for r in selected if r['method']==m]) for m in METHODS}
    summary=dict(version='topology_v2',geometry_count=len(metadata),run_count=len(rows),methods=by_method(rows),core_v1=by_method([r for r in rows if r['cohort']=='core_v1']),parameters=asdict(CANDIDATE),production_eligible=False)
    for field in ('geometry','phase','spacing_z_mm','orientation','cohort'):
        summary['by_'+field]={str(v):by_method([r for r in rows if r[field]==v]) for v in sorted({r[field] for r in rows})}
    summary['phase_core_canonical']={str(v):by_method([r for r in rows if r['phase']==v and r['cohort']=='core_v1' and r['orientation']=='canonical']) for v in (0.,.25,.5)}
    summary['figure_eight']={m:{k:sum(bool(r.get(k)) for r in rows if r['method']==m and r['geometry']=='figure_eight') for k in ('topology_correct','junction_count_correct','cycle_rank_correct','central_junction_structure_correct')} for m in METHODS}
    summary['small_bridge']={str(d):{m:sum(bool(r.get('communicating_bridge_recovered')) for r in rows if r['geometry']=='small_bridge' and r['method']==m and r['bridge_diameter_mm']==d) for m in METHODS} for d in (1.,1.5,2.)}
    summary['near_touch']={m:{'runs':len(g),'extracted':sum(r['extraction_success'] for r in g),'false_connection_runs':sum(bool(r.get('false_connection_count')) for r in g),'false_connection_rate':np.mean([bool(r.get('false_connection_count')) for r in g if r['extraction_success']])} for m in METHODS for g in [[r for r in rows if r['geometry']=='near_touch' and r['method']==m]]}
    summary['sweep']={name:summarize([r for r in tables['refinement_sweep'] if r['configuration']==name]) for name,_ in configurations()}
    summary['endpoint_sweep']={mode:summarize([r for r in tables['endpoint_sweep'] if r['endpoint_mode']==mode]) for mode in ('raw','recentered','tangent','boundary')}
    summary['refinement_sweep_run_count']=len(tables['refinement_sweep'])
    summary['endpoint_sweep_run_count']=len(tables['endpoint_sweep'])
    summary['candidate_selection']='Declared central settings retained as an experimental compromise; no uniform sweep winner, no per-geometry tuning. See global sweep and coverage regression.'
    summary['invariants']={'proposals':len(tables['invariant_results']),'passed':sum(r['topology_invariant_pass'] for r in tables['invariant_results']),'accepted':sum(r['accepted'] for r in tables['invariant_results'])}
    summary['dependencies']={p:importlib.metadata.version(p) for p in ('numpy','scipy','scikit-image','networkx','nibabel','matplotlib')}
    source=[*ROOT.glob('neurovasc/**/*.py'),*ROOT.glob('validation/synthetic_topology/*.py'),ROOT/'scripts/verify_topology_v2_regressions.py',ROOT/'tests/test_topology_validation_v2.py',Path(__file__).resolve()]
    summary['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source)}
    for name,table in tables.items():write_csv(output/f'{name}.csv',table)
    for name,value in [('summary',summary),('geometries',metadata)]: (output/f'{name}.json').write_text(json.dumps(serial(value),indent=2,allow_nan=False)+'\n')
    figures(output,rows,tables['invariant_results'],tables['spur_results'])
    files=[str(p.relative_to(output)) for p in output.rglob('*') if p.is_file() and p.name!='manifest.txt' and 'regressions' not in p.relative_to(output).parts]
    (output/'manifest.txt').write_text('\n'.join(sorted(files+['manifest.txt']))+'\n')
    print(json.dumps(serial(summary['core_v1']),indent=2))

if __name__=='__main__':main()
