"""Independent evaluation of physical line/graph outputs using existing truth."""
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
from neurovasc.geometry.centerline_refinement import resample
from neurovasc.geometry.cross_sectional_caliber import estimate_tangents
from validation.synthetic_centerline.truth import project, sample_truth
from validation.synthetic_centerline.evaluate import angular_errors
from .matching import endpoints

DIRECT='SUPPORTED_DIRECT_COMPARISON'
PARTIAL='PARTIAL_COMPARISON'
MISMATCH='METHOD_ASSUMPTION_MISMATCH'


def scope(geometry):
    if geometry in ('simple_ring','ring_branches','figure_eight','communicating_bridge'):
        return MISMATCH
    return PARTIAL if geometry=='bifurcation' else DIRECT


def graph_from_lines(lines):
    """Unify duplicate path segments, not an anatomical branch inference.

    1e-7 mm coordinate keys only identify duplicate VMTK source-target paths.
    Raw VMTK GroupIds/blanked tracts are retained separately.
    """
    graph=nx.Graph()
    for line in lines:
        p=np.asarray(line['points'],float)
        if p.ndim!=2 or p.shape[1]!=3 or len(p)<2 or not np.isfinite(p).all():
            raise ValueError('invalid_external_polyline')
        keys=[tuple(np.round(q,7)) for q in p]
        for k,q in zip(keys,p): graph.add_node(k,world=q.tolist())
        for a,b,p0,p1 in zip(keys[:-1],keys[1:],p[:-1],p[1:]):
            length=float(np.linalg.norm(p1-p0))
            if a!=b and length>1e-9: graph.add_edge(a,b,length_mm=length)
    if graph.number_of_edges()==0: raise ValueError('empty_external_graph')
    return graph


def clean_line(points):
    p=np.asarray(points,float)
    return p[np.r_[True,np.linalg.norm(np.diff(p,axis=0),axis=1)>1e-9]]


def centerline(vessel, graph, paths):
    dense=[]; angles=[]
    for path in paths:
        p=clean_line(path)
        if len(p)<2: continue
        q=resample(p,.1);dense.append(q)
        mapped,t,_,_=project(vessel,q)
        angle=angular_errors(estimate_tangents(q,6.),t)
        if vessel.geometry=='bifurcation': angle[np.linalg.norm(mapped,axis=1)<=3]=np.nan
        angles.extend(angle[np.isfinite(angle)])
    if not dense: raise ValueError('no_finite_paths')
    points=np.unique(np.round(np.concatenate(dense),9),axis=0)
    truthpoints=project(vessel,points)[0]
    error=np.linalg.norm(points-truthpoints,axis=1)
    if vessel.geometry=='bifurcation':
        truth=np.concatenate([sample_truth(vessel,name,samples=401) for name in vessel.branches])
        terminals=np.array([s[-1] for s in vessel.branches.values()])
        true_length=sum(np.linalg.norm(s[-1]-s[0]) for s in vessel.branches.values())
    else:
        truth=sample_truth(vessel,samples=601);terminals=truth[[0,-1]];true_length=vessel.length_mm
    tolerance=max(vessel.volume.spacing)
    support_dist=cKDTree(points).query(truth)[0]
    actual_length=sum(d['length_mm'] for *_,d in graph.edges(data=True))
    found=[d['world'] for n,d in graph.nodes(data=True) if graph.degree[n]==1]
    groups=list(nx.connected_components(graph.subgraph([n for n in graph if graph.degree[n]>=3])))
    row=dict(mean_position_error_mm=float(error.mean()),position_rmse_mm=float(np.sqrt(np.mean(error**2))),
             coverage_fraction=float(np.mean(support_dist<=tolerance)),coverage_tolerance_mm=tolerance,
             recovered_length_mm=float(actual_length),true_length_mm=float(true_length),
             length_error_mm=float(abs(actual_length-true_length)),mean_tangent_error_deg=float(np.mean(angles)) if angles else None,
             **endpoints(found,terminals), recovered_endpoint_count=len(found),recovered_junction_count=len(groups),
             recovered_cycle_rank=graph.number_of_edges()-len(graph)+nx.number_connected_components(graph),
             metric_scope='unique network edges; physical sampled coverage')
    return row, dict(truth_points=truth,nearest_error=support_dist,covered=support_dist<=tolerance)


def neurovasc_recoveries(volume):
    from neurovasc.geometry.hybrid_centerline import recover,hybrid
    from neurovasc.geometry.centerline_refinement import Refiner
    from neurovasc.geometry.centerline_qc import Recovery
    lee=recover(volume,'lee');raw=recover(volume,'geodesic_raw')
    try: refined=Refiner(volume).refine(raw) if raw.path is not None else Recovery(error=raw.error)
    except ValueError as e: refined=Recovery(error=str(e))
    chosen,diag,_=hybrid(volume,lee_result=lee,raw_result=raw,refined_result=refined)
    chosen.metadata={**chosen.metadata,'hybrid_selection':diag}
    return dict(lee=lee,geodesic_raw=raw,geodesic_refined=refined,hybrid=chosen)


def network_paths(graph):
    from neurovasc.graph.topology_network import reduce_network
    return [b['points'] for b in reduce_network(graph).branches.values() if len(b['points'])>=2 and b['length_mm']>0]


def y_branches(vessel,graph):
    from validation.synthetic_topology.geometries import phantom
    from validation.synthetic_topology.evaluate import score
    truth=phantom('y',vessel.volume.spacing)
    if not np.array_equal(truth.volume.data,vessel.volume.data): raise ValueError('Y_truth_mask_mismatch')
    row,branches,_,_,network,nodes=score(truth,graph)
    from scipy.ndimage import distance_transform_edt, map_coordinates
    import nibabel as nib
    radii=distance_transform_edt(vessel.volume.data>0,sampling=vessel.volume.spacing)
    radius_errors=[]
    for detail in branches:
        if detail['status']!='matched': continue
        branch=network.branches[detail['recovered_branch']]
        reference=truth.branches[detail['true_branch']]
        p=branch['points']; chord=np.linalg.norm(p[-1]-p[0])
        detail['distance_metric_tortuosity']=branch['length_mm']/chord if chord>1e-8 else None
        # Exclude 3 mm around the junction and 2 mm from analytic caps.
        axis=reference['points'][-1]/reference['length_mm'];along=p@axis
        keep=(along>=3)&(along<=reference['length_mm']-2)
        xyz=nib.affines.apply_affine(np.linalg.inv(vessel.volume.affine),p[keep])
        sampled=map_coordinates(radii,xyz.T,order=1,mode='constant',cval=np.nan)
        detail['radius_definition']='EDT interpolated at recovered branch points; descriptive on own support'
        detail['mean_radius_mm']=float(np.nanmean(sampled)) if len(sampled) else None
        detail['true_radius_mm']=reference['diameter_mm']/2
        if detail['mean_radius_mm'] is not None:
            radius_errors.append(abs(detail['mean_radius_mm']-detail['true_radius_mm']))
    tp=sum(b['status']=='matched' for b in branches)
    total=len(network.branches); precision=tp/total if total else 0.;recall=tp/3
    out=dict(true_branch_count=3,recovered_branch_count=total,branch_precision=precision,branch_recall=recall,
             branch_f1=2*precision*recall/(precision+recall) if precision+recall else 0.,
             branch_connectivity_correct=row['branch_connectivity_correct'],
             mean_radius_error_mm=float(np.mean(radius_errors)) if radius_errors else None,
             mean_branch_length_error_mm=float(np.mean([b['length_error_mm'] for b in branches if b['status']=='matched'])) if tp else None,
             junction_error_mm=float(np.mean([n['error_mm'] for n in nodes if n['node_kind']=='junction' and n['status']=='matched'])) if any(n['node_kind']=='junction' and n['status']=='matched' for n in nodes) else None)
    return out,branches


def vmtk_y_groups(vessel,external):
    """Match actual nonblanked GroupIds; retain their distinct cut boundaries."""
    from scipy.optimize import linear_sum_assignment
    geometry=external.get('geometry',{})
    group_ids=geometry.get('GroupIds',[]);targets=list(vessel.branches)
    unavailable=dict(recovered_branch_count=None,branch_precision=None,branch_recall=None,branch_f1=None,
        branch_connectivity_correct=None,mean_branch_length_error_mm=None,mean_radius_error_mm=None,
        junction_error_mm=None,branch_status='failed',branch_failure_reason='no_actual_vmtk_branch_groups')
    if not group_ids: return [],unavailable
    groups=[]
    for i,gid in enumerate(group_ids):
        tracts=[t for t in external['tracts'] if t.get('GroupIds')==gid and not t.get('Blanking')]
        if not tracts: continue
        points=np.concatenate([t['points'] for t in tracts]);radii=np.concatenate([t['radius_mm'] for t in tracts])
        terminal=points[np.argmax(np.linalg.norm(points,axis=1))]
        groups.append((i,gid,points,radii,terminal,tracts))
    if not groups:return [],unavailable
    costs=np.array([[np.linalg.norm(g[4]-vessel.branches[n][-1]) for n in targets] for g in groups])
    extended=np.full((len(groups),len(targets)+len(groups)),100.)
    extended[:,:len(targets)]=np.where(costs<=2.5,costs,1e6)
    ii,jj=linear_sum_assignment(extended);mapping={i:j for i,j in zip(ii,jj) if j<len(targets)}
    details=[]
    for index,(i,gid,points,radii,terminal,tracts) in enumerate(groups):
        raw=geometry['Tortuosity'][i]
        entry=dict(group_id=gid,length_mm=geometry['Length'][i],vmtk_tortuosity_raw=raw,
                   vmtk_tortuosity_plus_one=raw+1,curvature=geometry.get('Curvature',[None]*len(group_ids))[i],
                   torsion=geometry.get('Torsion',[None]*len(group_ids))[i],status='unmatched')
        if index in mapping:
            name=targets[mapping[index]];segment=vessel.branches[name];length=np.linalg.norm(segment[-1]-segment[0])
            along=points@(segment[-1]/length);keep=(along>=3)&(along<=length-2)
            mean=float(np.mean(radii[keep])) if keep.any() else None
            p=np.array(tracts[0]['points']);chord=np.linalg.norm(p[-1]-p[0])
            entry.update(true_branch=name,status='matched',true_length_mm=float(length),
                         length_error_mm=abs(entry['length_mm']-length),true_radius_mm=vessel.diameters[name]/2,
                         mean_radius_mm=mean,mean_radius_error_mm=abs(mean-vessel.diameters[name]/2) if mean is not None else None,
                         first_tract_distance_metric=float(np.linalg.norm(np.diff(p,axis=0),axis=1).sum()/chord) if chord>0 else None)
        details.append(entry)
    valid=[d for d in details if d['status']=='matched']
    radius=[d['mean_radius_error_mm'] for d in valid if d['mean_radius_error_mm'] is not None]
    # Actual VMTK connectivity is defined by GroupIds and blanked junction
    # tracts, not by reducing coincident/nearby points in duplicate paths.
    attachments={gid:set() for gid in group_ids}
    centerline_ids={t['CenterlineIds'] for t in external['tracts']}
    for cid in centerline_ids:
        ordered=sorted([t for t in external['tracts'] if t['CenterlineIds']==cid],key=lambda t:t['TractIds'])
        for i,tract in enumerate(ordered):
            if tract['Blanking']:continue
            for j in (i-1,i+1):
                if 0<=j<len(ordered) and ordered[j]['Blanking']:
                    attachments[tract['GroupIds']].add(ordered[j]['GroupIds'])
    refs=external.get('junctions_mm',[]);ref_ids=external.get('junction_arrays',{}).get('GroupIds',[])
    junction_valid=len(refs)==1 and len(ref_ids)==1 and np.linalg.norm(refs[0])<=2.5
    tp=sum(d['status']=='matched' and junction_valid and attachments[d['group_id']]=={ref_ids[0]} for d in details)
    precision=tp/len(group_ids) if group_ids else 0.;recall=tp/3
    for d in details:d['attached_junction_group_ids']=sorted(attachments[d['group_id']])
    return details,dict(branch_status='success',recovered_branch_count=len(group_ids),branch_precision=precision,branch_recall=recall,
        branch_f1=2*precision*recall/(precision+recall) if precision+recall else 0.,
        branch_connectivity_correct=tp==3 and len(group_ids)==3,
        mean_branch_length_error_mm=float(np.mean([d['length_error_mm'] for d in valid])) if valid else None,
        junction_error_mm=float(np.linalg.norm(refs[0])) if len(refs)==1 else None,
        vmtk_matched_group_count=len(valid),vmtk_missed_true_branches=3-len(valid),
        vmtk_group_length_error_mm=float(np.mean([d['length_error_mm'] for d in valid])) if valid else None,
        mean_radius_error_mm=float(np.mean(radius)) if radius else None)
