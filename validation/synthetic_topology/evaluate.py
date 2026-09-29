"""Physical one-to-one branch matching and full-network metrics; truth only here."""
from collections import Counter
import networkx as nx
import nibabel as nib
import numpy as np
from scipy.ndimage import distance_transform_edt,map_coordinates,label,generate_binary_structure
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from neurovasc.geometry.centerline import extract_centerline
from neurovasc.graph.vascular_graph import build_vascular_graph
from neurovasc.geometry.hybrid_centerline import recover,hybrid
from neurovasc.geometry.centerline_refinement import Refiner,resample
from neurovasc.geometry.centerline_qc import Recovery
from neurovasc.geometry.topology_centerline import extract_topology_network
from neurovasc.graph.topology_network import reduce_network,fundamental_cycles
from neurovasc.graph.topology_qc import topology_qc

METHODS=('lee','geodesic_raw','geodesic_refined','hybrid','topology')
MATCH_TOLERANCE_MM=2.5


def recover_all(volume):
    results={}
    # Lee full graph is retained even for disconnected masks and zero-endpoint
    # rings. Its current hybrid still applies the original single-component QC.
    try:results['lee']=(build_vascular_graph(extract_centerline(volume)).graph,'',True,{})
    except ValueError as e:results['lee']=(None,str(e),True,{})
    raw=recover(volume,'geodesic_raw')
    results['geodesic_raw']=(raw.graph,raw.error,False,{})
    try:refined=Refiner(volume).refine(raw) if raw.path is not None else Recovery(error=raw.error,supports_cycles=False)
    except ValueError as e:refined=Recovery(error=str(e),supports_cycles=False)
    results['geodesic_refined']=(refined.graph,refined.error,False,refined.metadata)
    chosen,diagnostics,qc=hybrid(volume,raw_result=raw,refined_result=refined)
    results['hybrid']=(chosen.graph,chosen.error,chosen.supports_cycles,diagnostics)
    try:results['topology']=(extract_topology_network(volume),'',True,{})
    except ValueError as e:results['topology']=(None,str(e),True,{})
    return results


def dense_edges(graph,step=.15):
    arrays=[resample(np.array([graph.nodes[a]['world'],graph.nodes[b]['world']]),step) for a,b in graph.edges if np.linalg.norm(np.asarray(graph.nodes[a]['world'])-graph.nodes[b]['world'])>1e-10]
    return np.concatenate(arrays) if arrays else np.array([d['world'] for _,d in graph.nodes(data=True)])


def branch_points(branch):
    return resample(branch['points'],.15) if len(branch['points'])>1 and branch['length_mm']>0 else branch['points']


def match_structures(network,phantom):
    actual,truth=network.logical,phantom.graph
    mapping={};node_rows=[]
    def location(graph,branches,node):
        if graph.nodes[node]['node_kind']!='regular':return np.array(graph.nodes[node]['world'])
        return np.mean(np.concatenate([branch['points'] for branch in branches.values() if node in (branch['start'],branch['end'])]),axis=0)
    for kind in ('endpoint','junction','regular'):
        a=[n for n,d in actual.nodes(data=True) if d['node_kind']==kind]
        t=[n for n,d in truth.nodes(data=True) if d['node_kind']==kind]
        if not a or not t:continue
        costs=np.array([[np.linalg.norm(location(actual,network.branches,n)-location(truth,phantom.branches,m)) for m in t] for n in a])
        # Dummy assignments enforce the gate during optimization, not after an
        # unconstrained assignment that could sacrifice an available valid match.
        extended=np.full((len(a),len(t)+len(a)),(len(a)+len(t)+1)*MATCH_TOLERANCE_MM)
        extended[:,:len(t)]=np.where(costs<=MATCH_TOLERANCE_MM,costs,1e6)
        ia,it=linear_sum_assignment(extended)
        for i,j in zip(ia,it):
            if j<len(t):
                mapping[a[i]]=t[j];node_rows.append(dict(node_kind=kind,recovered_node=a[i],true_node=t[j],error_mm=float(costs[i,j]),status='matched'))
    node_rows.extend(dict(node_kind=d['node_kind'], recovered_node=n, true_node=None, error_mm=None, status='spurious_or_unmatched') for n,d in actual.nodes(data=True) if n not in mapping)
    node_rows.extend(dict(node_kind=d['node_kind'], recovered_node=None, true_node=n, error_mm=None, status='missed') for n,d in truth.nodes(data=True) if n not in mapping.values())
    candidates=list(network.branches);targets=list(phantom.branches)
    cost=np.full((len(candidates),len(targets)+len(candidates)),(len(candidates)+len(targets)+1)*MATCH_TOLERANCE_MM)
    cost[:,:len(targets)]=1e6
    for i,key in enumerate(candidates):
        b=network.branches[key]
        if b['start'] not in mapping or b['end'] not in mapping:continue
        ends=sorted([mapping[b['start']],mapping[b['end']]])
        points=branch_points(b)
        for j,target in enumerate(targets):
            t=phantom.branches[target]
            if ends!=sorted([t['start'],t['end']]):continue
            tpoints=t['points']
            metric=float((cKDTree(tpoints).query(points)[0].mean()+cKDTree(points).query(tpoints)[0].mean())/2)
            if metric<=MATCH_TOLERANCE_MM:cost[i,j]=metric
    branch_map={};rows=[]
    if candidates:
        ia,it=linear_sum_assignment(cost)
        for i,j in zip(ia,it):
            if j<len(targets):
                branch_map[candidates[i]]=targets[j]
                b,t=network.branches[candidates[i]],phantom.branches[targets[j]]
                rows.append(dict(recovered_branch=candidates[i],true_branch=targets[j],status='matched',position_error_mm=float(cost[i,j]),
                                 length_error_mm=abs(b['length_mm']-t['length_mm']),communicating_bridge=t['communicating_bridge']))
    rows += [dict(recovered_branch=k,true_branch=None,status='spurious_or_unmatched',communicating_bridge=False) for k in candidates if k not in branch_map]
    rows += [dict(recovered_branch=None,true_branch=k,status='missed',communicating_bridge=phantom.branches[k]['communicating_bridge']) for k in targets if k not in branch_map.values()]
    return mapping,branch_map,node_rows,rows


def score(phantom,graph,error='',supports_cycles=True):
    truth=phantom.metadata()
    row={'extraction_success':graph is not None and len(graph)>0,'failure_reason':error,'supports_cycles':supports_cycles,
         'true_components':truth['true_components'],'true_endpoint_count':truth['true_endpoints'],
         'true_junction_count':truth['true_junctions'],'true_branch_count':truth['true_branch_count'],'true_cycle_rank':truth['true_cycle_rank'],
         'true_network_length_mm':truth['true_length_mm']}
    if not row['extraction_success']:
        row.update({k:False for k in ('component_count_correct','endpoint_count_correct','junction_count_correct','cycle_rank_correct','branch_count_correct','branch_connectivity_correct','topology_correct','qc_pass')})
        row['communicating_bridge_recovered']=False if any(b['communicating_bridge'] for b in phantom.branches.values()) else None
        row.update(topology_failure_reason='extraction_failed',number_of_cycles_found=0,
                   number_of_cycles_missed=truth['true_cycle_rank'],number_of_spurious_cycles=None,
                   false_connection_count=None,spurious_cross_component_edges=None,missed_true_components=truth['true_components'])
        branches=[dict(recovered_branch=None,true_branch=k,status='missed_extraction_failed',communicating_bridge=b['communicating_bridge']) for k,b in phantom.branches.items()]
        cycles=[dict(cycle_kind='analytic',cycle_id=i,true_members=';'.join(c),cycle_preserved=False,membership_correct=False,cycle_coverage=None) for i,c in enumerate(truth['cycle_basis'])]
        checks=[dict(criterion='extraction_exists',value=False,threshold='nonempty graph',passed=False,critical=True)]
        return row,branches,cycles,checks,None,[]
    network=reduce_network(graph);counts=network.counts()
    for source,key in [('components','components'),('endpoints','endpoint_count'),('junctions','junction_count'),('branches','branch_count'),('cycle_rank','cycle_rank')]:
        row['recovered_'+key]=counts[source]
    for name,field in [('component','components'),('endpoint','endpoint_count'),('junction','junction_count'),('branch','branch_count'),('cycle_rank','cycle_rank')]:
        row[name+'_count_correct' if name!='cycle_rank' else 'cycle_rank_correct']=row['recovered_'+field]==row['true_'+field]
    mapping,bmap,nodes,branches=match_structures(network,phantom)
    missed=len(phantom.branches)-len(bmap);spurious=len(network.branches)-len(bmap)
    row.update(missed_connections=missed,spurious_connections=spurious,branch_connectivity_correct=missed==0 and spurious==0,
               unmatched_recovered_nodes=len(network.logical)-len(mapping),unmatched_true_nodes=len(phantom.graph)-len(mapping))
    for kind in ('endpoint','junction'):
        errors=[n['error_mm'] for n in nodes if n['node_kind']==kind and n['error_mm'] is not None]
        row[kind+'_position_error_mm']=float(np.mean(errors)) if errors else None
        row[kind+'_matched_count']=len(errors)
    row['communicating_bridge_recovered']=all(k in bmap.values() for k,b in phantom.branches.items() if b['communicating_bridge']) if any(b['communicating_bridge'] for b in phantom.branches.values()) else None
    dense=dense_edges(graph);actual_tree=cKDTree(dense)
    reference=np.concatenate([b['points'] for b in phantom.branches.values()]);reference_tree=cKDTree(reference)
    errors=reference_tree.query(dense)[0]
    tolerance=max(phantom.volume.spacing)
    row.update(total_network_length_mm=float(sum(d['length_mm'] for _,_,d in graph.edges(data=True))),
               network_coverage_fraction=float(np.mean(actual_tree.query(reference)[0]<=tolerance)),mean_network_position_error_mm=float(errors.mean()),
               coverage_tolerance_mm=tolerance,graph_complexity_ratio=counts['branches']/max(1,len(phantom.branches)))
    row['network_length_error_mm']=abs(row['total_network_length_mm']-row['true_network_length_mm'])
    # Fewer recovered components can mean a missed vessel, not a false bridge.
    # Attribute foreground graph nodes to nearest analytic connected components
    # and count actual edges crossing those labels separately from missing support.
    true_components=list(nx.connected_components(phantom.graph))
    component_samples=[np.concatenate([b['points'] for b in phantom.branches.values() if b['start'] in component]) for component in true_components]
    graph_nodes=list(graph)
    graph_points=np.array([graph.nodes[n]['world'] for n in graph_nodes])
    labels=np.argmin(np.stack([cKDTree(p).query(graph_points)[0] for p in component_samples]),axis=0)
    attribution=dict(zip(graph_nodes,map(int,labels)))
    row['spurious_cross_component_edges']=sum(attribution[a]!=attribution[b] for a,b in graph.edges)
    row['false_connection_count']=sum(len({attribution[n] for n in component})>1 for component in nx.connected_components(graph))
    row['missed_true_components']=len(true_components)-len(set(attribution.values()))
    radii=distance_transform_edt(np.pad(phantom.volume.data>.5,1),sampling=phantom.volume.spacing)[1:-1,1:-1,1:-1]
    vox=nib.affines.apply_affine(np.linalg.inv(phantom.volume.affine),dense)
    clearance=map_coordinates(radii,vox.T,order=1,mode='constant',cval=0)
    row.update(mean_clearance_mm=float(clearance.mean()),minimum_clearance_mm=float(clearance.min()))
    cycles=[];inverse={v:k for k,v in bmap.items()}
    for i,cycle in enumerate(truth['cycle_basis']):
        points=np.concatenate([phantom.branches[k]['points'] for k in cycle])
        present=all(k in inverse for k in cycle)
        if present:
            selected=nx.MultiGraph()
            for k in cycle:
                b=network.branches[inverse[k]];selected.add_edge(b['start'],b['end'])
            present=nx.is_connected(selected) and all(d==2 for _,d in selected.degree())
        perimeter=sum(phantom.branches[k]['length_mm'] for k in cycle)
        recovered_perimeter=sum(network.branches[inverse[k]]['length_mm'] for k in cycle) if present else None
        cycles.append(dict(cycle_kind='analytic',cycle_id=i,true_members=';'.join(cycle),recovered_members=';'.join(inverse[k] for k in cycle if k in inverse),
                           membership_correct=present,cycle_preserved=present,cycle_coverage=float(np.mean(actual_tree.query(points)[0]<=tolerance)),
                           true_perimeter_mm=perimeter,recovered_perimeter_mm=recovered_perimeter,perimeter_error_mm=abs(recovered_perimeter-perimeter) if present else None))
    # Any recovered basis cycle involving an unmatched branch cannot be certified
    # as a true cycle. Record it explicitly without calling basis differences alone spurious.
    for i,cycle in enumerate(network.cycles):
        certified=all(k in bmap for k in cycle)
        cycles.append(dict(cycle_kind='recovered_basis',cycle_id=i,recovered_members=';'.join(cycle),membership_correct=certified,
                           true_members=';'.join(bmap[k] for k in cycle if k in bmap),recovered_perimeter_mm=sum(network.branches[k]['length_mm'] for k in cycle)))
    found=sum(c['cycle_preserved'] for c in cycles if c['cycle_kind']=='analytic')
    row.update(number_of_cycles_found=found,number_of_cycles_missed=truth['true_cycle_rank']-found,
               number_of_spurious_cycles=max(0,counts['cycle_rank']-truth['true_cycle_rank']),
               uncertified_recovered_basis_cycles=sum(not c['membership_correct'] for c in cycles if c['cycle_kind']=='recovered_basis'))
    row['topology_correct']=all(row[k] for k in ('component_count_correct','endpoint_count_correct','junction_count_correct','branch_count_correct','cycle_rank_correct','branch_connectivity_correct'))
    reasons=[]
    for k,text in [('component_count_correct','component_mismatch'),('endpoint_count_correct','wrong_endpoint_count'),('junction_count_correct','wrong_junction_count'),('branch_count_correct','wrong_branch_count')]:
        if not row[k]:reasons.append(text)
    if counts['cycle_rank']<truth['true_cycle_rank']:reasons.append('missed_cycle' if supports_cycles else 'topology_collapse_into_tree')
    if counts['cycle_rank']>truth['true_cycle_rank']:reasons.append('spurious_cycle')
    if row['false_connection_count']:reasons.append('false_bridge')
    if row['missed_true_components']:reasons.append('missing_component')
    if missed:reasons.append('missed_connections')
    if spurious:reasons.append('spurious_connections')
    if row['communicating_bridge_recovered'] is False:reasons.append('missed_bridge')
    row['topology_failure_reason']=';'.join(reasons)
    qc=topology_qc(network,supports_cycles,max(phantom.volume.spacing))
    row['qc_pass']=all(c['passed'] for c in qc if c['critical'])
    return row,branches,cycles,qc,network,nodes


def aggregate(rows):
    valid=[r for r in rows if r['extraction_success']]
    result={'runs':len(rows),'extraction_success_rate':len(valid)/len(rows) if rows else None}
    for field in ('component_count_correct','endpoint_count_correct','junction_count_correct','branch_count_correct','cycle_rank_correct','branch_connectivity_correct','topology_correct'):
        result[field]=sum(r[field] for r in rows)/len(rows) if rows else None
    for field in ('network_coverage_fraction','mean_network_position_error_mm','network_length_error_mm'):
        values=[r[field] for r in valid if r.get(field) is not None]
        result[field]=float(np.mean(values)) if values else None
    result['failure_inclusive_coverage']=sum(r.get('network_coverage_fraction',0) for r in valid)/len(rows) if rows else None
    aliases={'component_accuracy':'component_count_correct','endpoint_count_accuracy':'endpoint_count_correct','junction_count_accuracy':'junction_count_correct','cycle_rank_accuracy':'cycle_rank_correct','branch_count_accuracy':'branch_count_correct','connectivity_accuracy':'branch_connectivity_correct','mean_network_coverage':'network_coverage_fraction','mean_position_error_mm':'mean_network_position_error_mm','mean_network_length_error_mm':'network_length_error_mm'}
    result.update({name:result[key] for name,key in aliases.items()})
    return result
