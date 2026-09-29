"""Evaluation-only geometry truth, kept separate from QC and hybrid decisions."""
import networkx as nx
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from neurovasc.geometry.centerline_qc import smoothness, quality_control, occupancy
from neurovasc.geometry.centerline_refinement import resample
from neurovasc.geometry.cross_sectional_caliber import estimate_tangents
from .evaluate import score_points, topology, branch_scores, angular_errors
from .truth import sample_truth
from .experiments_v2 import loop_truth, loop_samples

METHODS=('lee','geodesic_raw','geodesic_refined','hybrid')


def endpoint_metrics(vessel, recovery):
    g=recovery.graph
    found=np.array([d['world'] for n,d in g.nodes(data=True) if g.degree[n]==1])
    if vessel.geometry.startswith('loop'):
        truth=np.array([s[-1] for s in vessel.branches.values()])
    elif vessel.geometry=='bifurcation':
        truth=np.array([s[-1] for s in vessel.branches.values()])
    else:
        truth=sample_truth(vessel,samples=2)
    if len(truth)==0 or len(found)==0:
        return {'mean_endpoint_error_mm':None,'max_endpoint_error_mm':None,'matched_endpoint_count':0}
    distances=np.linalg.norm(found[:,None]-truth[None],axis=2)
    a,b=linear_sum_assignment(distances)
    errors=distances[a,b]
    return {'mean_endpoint_error_mm':float(errors.mean()),'max_endpoint_error_mm':float(errors.max()),'matched_endpoint_count':len(a)}


def measure(vessel, recovery):
    qc=quality_control(vessel.volume,recovery)
    row={'extraction_success':recovery.path is not None,'failure_reason':recovery.error,
         'qc_pass':qc['qc_pass'],'qc_fail_count':qc['qc_fail_count'],'qc_fail_reasons':qc['qc_fail_reasons'],
         'supports_cycles':recovery.supports_cycles,'recovered_sample_count':0}
    if recovery.path is None:
        return row,[],[],None,qc
    graph,p=recovery.graph,recovery.path
    stats,groups=topology(graph)
    row.update(stats,**smoothness(p),**endpoint_metrics(vessel,recovery),recovered_sample_count=len(p),
               points_outside_lumen_count=int(np.sum(occupancy(vessel.volume,p)<.5-1e-8)))
    branches,loop=[],None
    if vessel.geometry.startswith('loop'):
        points=np.array([d['world'] for _,d in graph.nodes(data=True)])
        truth,_=loop_truth(vessel,points)
        errors=np.linalg.norm(points-truth,axis=1)
        dense=np.concatenate([resample(np.array([graph.nodes[a]['world'],graph.nodes[b]['world']]),.2) for a,b in graph.edges])
        reference=np.concatenate(loop_samples(vessel))
        distance=cKDTree(dense).query(reference)[0]
        coverage=float(np.mean(distance<=max(vessel.volume.spacing)))
        length=sum(d['length_mm'] for _,_,d in graph.edges(data=True))
        _,tangent=loop_truth(vessel,p)
        angles=angular_errors(estimate_tangents(p,6.),tangent)
        if vessel.branches:
            junctions=np.array([s[0] for s in vessel.branches.values()])
            angles[np.min(np.linalg.norm(p[:,None]-junctions[None],axis=2),axis=1)<=3]=np.nan
        valid=angles[np.isfinite(angles)]
        row.update(true_length_mm=vessel.length_mm,recovered_length_mm=float(length),length_error_mm=abs(length-vessel.length_mm),
                   length_error_percent=100*abs(length-vessel.length_mm)/vessel.length_mm,coverage_fraction=coverage,
                   normalized_start=None,normalized_end=None,mean_position_error_mm=float(errors.mean()),
                   median_position_error_mm=float(np.median(errors)),max_position_error_mm=float(errors.max()),position_rmse_mm=float(np.sqrt(np.mean(errors**2))),
                   mean_tangent_error_deg=float(valid.mean()) if len(valid) else None,
                   median_tangent_error_deg=float(np.median(valid)) if len(valid) else None,max_tangent_error_deg=float(valid.max()) if len(valid) else None,
                   interior_position_error_mm=None,metric_scope='network; tangent/smoothness on principal path')
        correct=stats['cycle_count']==1 and stats['endpoint_count']==len(vessel.branches) and stats['junction_count']==len(vessel.branches)
        expected_junctions=np.array([segment[0] for segment in vessel.branches.values()])
        recovered_junctions=np.array([np.mean([graph.nodes[n]['world'] for n in group],axis=0) for group in groups])
        junction_errors=[]
        if len(expected_junctions) and len(recovered_junctions):
            costs=np.linalg.norm(recovered_junctions[:,None]-expected_junctions[None],axis=2)
            a,b=linear_sum_assignment(costs)
            junction_errors=costs[a,b]
        loop={**stats,'true_cycle_count':1,
              'mean_junction_error_mm':float(np.mean(junction_errors)) if len(junction_errors) else None,
              'max_junction_error_mm':float(np.max(junction_errors)) if len(junction_errors) else None,
              'unmatched_junction_count':abs(len(expected_junctions)-len(recovered_junctions)),'topology_preserved':correct,'network_coverage_fraction':coverage,
              'coverage_tolerance_mm':max(vessel.volume.spacing),'true_network_length_mm':vessel.length_mm,'network_length_mm':float(length),
              'mean_position_error_mm':float(errors.mean()),'supports_cycles':recovery.supports_cycles}
        samples=[{'sample_index':i,'x_mm':q[0],'y_mm':q[1],'z_mm':q[2],'position_error_mm':e,'path_kind':'network'} for i,(q,e) in enumerate(zip(points,errors))]
    else:
        metrics,samples=score_points(vessel,p)
        row.update(metrics,metric_scope='principal path')
        arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]
        interior=(arc>=3)&(arc[-1]-arc>=3)
        errors=np.array([s['position_error_mm'] for s in samples])
        row['interior_position_error_mm']=float(errors[interior].mean()) if interior.any() else None
        samples=[{**s,'path_kind':'main'} for s in samples]
        if vessel.geometry=='bifurcation':
            branches,branch_samples=branch_scores(vessel,graph,groups)
            samples+=branch_samples
            row['junction_error_mm']=float(np.linalg.norm(np.mean([graph.nodes[n]['world'] for n in groups[0]],axis=0))) if len(groups)==1 else None
            row['topology_preserved']=stats['endpoint_count']==3 and stats['junction_count']==1 and stats['cycle_count']==0
        else:
            row['topology_preserved']=stats['endpoint_count']==2 and stats['junction_count']==0 and stats['cycle_count']==0
    row['absolute_length_error_mm']=row['length_error_mm']
    row.update({k:v for k,v in recovery.metadata.items() if k in ('mean_displacement_mm','max_displacement_mm','points_outside_lumen_count')})
    return row,samples,branches,loop,qc


def aggregate(rows):
    valid=[r for r in rows if r['extraction_success']]
    keys=['mean_position_error_mm','coverage_fraction','length_error_mm','mean_tangent_error_deg','mean_endpoint_error_mm','interior_position_error_mm','mean_turn_deg','mean_curvature_per_mm']
    result={'runs':len(rows),'success_count':len(valid),'success_rate':len(valid)/len(rows) if rows else None,
            'qc_pass_rate':sum(r['qc_pass'] for r in rows)/len(rows) if rows else None,
            'failure_inclusive_coverage':sum(r.get('coverage_fraction',0) for r in valid)/len(rows) if rows else None}
    for key in keys:
        values=[r[key] for r in valid if r.get(key) is not None]
        result[key]=float(np.mean(values)) if values else None
    return result
