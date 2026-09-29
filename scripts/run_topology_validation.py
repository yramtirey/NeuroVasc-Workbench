"""Deterministic topology validation; no production or archived output mutations."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'validation/outputs/.matplotlib'))
import matplotlib
matplotlib.use('Agg')
import nibabel as nib
import numpy as np
from scipy.ndimage import label,generate_binary_structure
from scripts.run_caliber_validation import write_csv
from scripts.run_caliber_validation_v2 import clean_json
from validation.synthetic_topology.geometries import matrix
from validation.synthetic_topology.evaluate import METHODS,MATCH_TOLERANCE_MM,recover_all,score,aggregate
from validation.synthetic_topology.figures import summaries,diagnostic


def output_directory(path):
    p=Path(path).resolve();allowed=(ROOT/'validation/outputs/topology_v1').resolve()
    if not p.is_relative_to(allowed):raise ValueError('Output must be topology_v1 or a descendant; archives protected')
    return p


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=ROOT/'validation/outputs/topology_v1')
    args=parser.parse_args()
    try:output=output_directory(args.output)
    except ValueError as e:parser.error(str(e))
    output.mkdir(parents=True,exist_ok=True)
    for name in ('masks','networks'): (output/name).mkdir(exist_ok=True)
    rows,branches,cycles,qcs,nodes,geometries=[],[],[],[],[],[]
    for eid,phantom in matrix():
        volume=phantom.volume
        meta={'geometry_id':eid,'geometry':phantom.geometry,'phase':phantom.phase,'orientation':phantom.orientation,
              'spacing_x_mm':float(volume.spacing[0]),'spacing_y_mm':float(volume.spacing[1]),'spacing_z_mm':float(volume.spacing[2]),
              'bridge_diameter_mm':phantom.bridge_diameter_mm,'minimum_separation_mm':phantom.minimum_separation_mm,
              'mask_components_6':label(volume.data>0,generate_binary_structure(3,1))[1],
              'mask_components_26':label(volume.data>0,np.ones((3,3,3)))[1]}
        geometries.append({**meta,**phantom.metadata()})
        nib.save(nib.Nifti1Image(volume.data,volume.affine),output/'masks'/f'{eid}.nii.gz')
        for method,(graph,error,capable,metadata) in recover_all(volume).items():
            result,branch_rows,cycle_rows,checks,network,node_rows=score(phantom,graph,error,capable)
            key={**meta,'method':method};row={**key,**result,'selected_method':metadata.get('selected_method')}
            rows.append(row);branches.extend({**key,**b} for b in branch_rows);cycles.extend({**key,**c} for c in cycle_rows)
            qcs.extend({**key,**c} for c in checks);nodes.extend({**key,**n} for n in node_rows)
            if network is not None:
                raw_ids={n:i for i,n in enumerate(graph)}
                memberships={}
                for branch_id,b in network.branches.items():
                    for point in b['points']:memberships.setdefault(tuple(point),set()).add(branch_id)
                raw_nodes=[{'node_id':raw_ids[n],'node_kind':'endpoint' if graph.degree[n]==1 else 'junction' if graph.degree[n]>=3 else 'regular',
                            'world':list(d['world']),'voxel':list(d['voxel']) if d.get('voxel') is not None else None,
                            'clearance_mm':d.get('clearance_mm'),'branch_ids':sorted(memberships.get(tuple(d['world']),set()))} for n,d in graph.nodes(data=True)]
                payload={'raw_nodes':raw_nodes,'raw_edges':[{'start':raw_ids[a],'end':raw_ids[b],'length_mm':d['length_mm']} for a,b,d in graph.edges(data=True)],'counts':network.counts(),'metadata':{**network.metadata,**graph.graph},
                         'nodes':{n:{k:v for k,v in d.items() if k!='raw_nodes'} for n,d in network.logical.nodes(data=True)},
                         'branches':{k:{kk:vv.tolist() if isinstance(vv,np.ndarray) else vv for kk,vv in b.items()} for k,b in network.branches.items()},
                         'cycles':network.cycles,'supports_cycles':capable,'production_eligible':False}
                (output/'networks'/f'{eid}_{method}.json').write_text(json.dumps(clean_json(payload),indent=2,allow_nan=False)+'\n')
            if not row['topology_correct'] or phantom.phase==0 and phantom.orientation=='canonical' and volume.spacing[2]==.5:
                diagnostic(output,eid,phantom,method,row,network,branch_rows)
        print(eid,':',', '.join(f"{r['method']}={r.get('recovered_cycle_rank','failed')}/{r['true_cycle_rank']}" for r in rows[-5:]),flush=True)
    comparisons=[]
    for eid in {r['geometry_id'] for r in rows}:
        group={r['method']:r for r in rows if r['geometry_id']==eid};new=group['topology']
        for method in METHODS[:-1]:
            old=group[method];pair={'geometry_id':eid,'baseline':method,'baseline_success':old['extraction_success'],'topology_success':new['extraction_success']}
            for k in ('cycle_rank_correct','branch_connectivity_correct','topology_correct'):pair[k+'_delta']=int(new[k])-int(old[k])
            for k in ('network_coverage_fraction','mean_network_position_error_mm','network_length_error_mm'):
                pair[k+'_delta']=new[k]-old[k] if new.get(k) is not None and old.get(k) is not None else None
            comparisons.append(pair)
    comparisons.sort(key=lambda r:(r['geometry_id'],r['baseline']))
    def by_method(selected):return {m:aggregate([r for r in selected if r['method']==m]) for m in METHODS}
    summary={'version':'topology_v1','geometry_count':len(geometries),'run_count':len(rows),'per_method':by_method(rows),
             'parameters':{'foreground_connectivity':6,'background_connectivity':26,'matching_tolerance_mm':MATCH_TOLERANCE_MM,
                           'maximum_foreground_voxels':150000,'truth_sampling_max_mm':.08,'recovered_sampling_max_mm':.15,
                           'cycle_basis':'deterministic fundamental; not minimum; parallel/self branches retained','cycles_pruned':0,'production_eligible':False},
             'dependencies':{p:importlib.metadata.version(p) for p in ('numpy','scipy','scikit-image','networkx','nibabel','matplotlib')},'python':platform.python_version()}
    for field in ('geometry','phase','spacing_z_mm','orientation'):
        summary['by_'+field]={str(v):by_method([r for r in rows if r[field]==v]) for v in sorted(set(r[field] for r in rows))}
    summary['small_bridge_results']={str(d):{m:{'runs':len(g),'bridge_recovered':sum(bool(r['communicating_bridge_recovered']) for r in g),'cycle_rank_correct':sum(r['cycle_rank_correct'] for r in g)} for m in METHODS for g in [[r for r in rows if r['geometry']=='small_bridge' and r['bridge_diameter_mm']==d and r['method']==m]]} for d in (1.,1.5,2.)}
    summary['near_touch_results']={str(d):{m:{'runs':len(g),'extracted':sum(r['extraction_success'] for r in g),'false_connection_runs':sum(r['extraction_success'] and r['false_connection_count']>0 for r in g),'components_correct':sum(r['component_count_correct'] for r in g),'missing_component_runs':sum(r['extraction_success'] and r['missed_true_components']>0 for r in g)} for m in METHODS for g in [[r for r in rows if r['geometry']=='near_touch' and r['minimum_separation_mm']==d and r['method']==m]]} for d in (.4,1.2)}
    summary['known_loop_results']={k:summary['by_geometry'][k] for k in ('simple_ring','ring_branches','figure_eight','communicating_bridge','incomplete_loop')}
    for group in summary['near_touch_results'].values():
        for result in group.values():
            result['false_connection_rate']=result['false_connection_runs']/result['extracted'] if result['extracted'] else None
            result['false_connection_rate_denominator']='successful extractions; extraction failures reported separately'
    source=[*ROOT.glob('neurovasc/**/*.py'),*ROOT.glob('validation/synthetic_topology/*.py'),*ROOT.glob('validation/synthetic_centerline/*.py'),ROOT/'scripts/run_topology_validation.py']
    summary['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source)}
    for name,data in [('results',rows),('network_metrics',rows),('branch_results',branches),('cycle_results',cycles),('qc_results',qcs),('matched_comparison',comparisons),('node_matches',nodes)]:write_csv(output/f'{name}.csv',data)
    for name,data in [('summary',summary),('geometries',geometries)]: (output/f'{name}.json').write_text(json.dumps(clean_json(data),indent=2,allow_nan=False)+'\n')
    summaries(output,rows)
    files=[str(p.relative_to(output)) for p in output.rglob('*') if p.is_file() and p.name!='manifest.txt' and 'regressions' not in p.relative_to(output).parts]
    (output/'manifest.txt').write_text('\n'.join(sorted(files+['manifest.txt']))+'\n')
    print(json.dumps(summary['per_method'],indent=2))

if __name__=='__main__':main()
