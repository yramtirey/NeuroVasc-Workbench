"""Evaluate fixed logical topology separately from its physical embedding."""
from collections import Counter
import networkx as nx
import numpy as np
from validation.synthetic_topology.evaluate import score, match_structures
from neurovasc.graph.topology_qc import topology_qc

METHODS=('lee','topology_v1','topology_v2_junction','topology_v2_refined')


def prf(tp,fp,fn):
    precision=tp/(tp+fp) if tp+fp else float(not fn)
    recall=tp/(tp+fn) if tp+fn else float(not fp)
    f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 1.
    return precision,recall,f1


def measure(phantom,network,error=''):
    if network is None:
        row,br,cy,qc,_,nodes=score(phantom,None,error)
        row.update(branch_precision=0.,branch_recall=0.,branch_f1=0.,adjacency_precision=0.,adjacency_recall=0.,adjacency_f1=0.)
        return row,br,cy,[],qc
    # Reuse v1 physical metrics on the embedding; structural evaluation always
    # uses the supplied fixed logical topology, not a reduction of smoothed points.
    row,_,_,qc,_,_=score(phantom,network.raw)
    qc=topology_qc(network,True,max(phantom.volume.spacing))
    mapping,bmap,nodes,branches=match_structures(network,phantom)
    truth=phantom.metadata();counts=network.counts()
    for actual,true,metric in (('components','components','component_count_correct'),('endpoints','endpoint_count','endpoint_count_correct'),('junctions','junction_count','junction_count_correct'),('branches','branch_count','branch_count_correct'),('cycle_rank','cycle_rank','cycle_rank_correct')):
        row['recovered_'+true]=counts[actual];row[metric]=counts[actual]==row['true_'+true]
    tp=len(bmap);fp=len(network.branches)-tp;fn=len(phantom.branches)-tp
    row.update(missed_connections=fn,spurious_connections=fp,branch_connectivity_correct=fp==fn==0,branch_true_positive=tp,branch_false_positive=fp,branch_false_negative=fn)
    row['branch_precision'],row['branch_recall'],row['branch_f1']=prf(tp,fp,fn)
    actual=Counter(tuple(sorted((mapping[b['start']],mapping[b['end']]))) for b in network.branches.values() if b['start'] in mapping and b['end'] in mapping)
    expected=Counter(tuple(sorted((b['start'],b['end']))) for b in phantom.branches.values())
    atp=sum((actual & expected).values());afp=len(network.branches)-atp;afn=len(phantom.branches)-atp
    row.update(adjacency_true_positive=atp,adjacency_false_positive=afp,adjacency_false_negative=afn)
    row['adjacency_precision'],row['adjacency_recall'],row['adjacency_f1']=prf(atp,afp,afn)
    junctions=[]
    for n in nodes:
        if n['node_kind'] not in ('junction','endpoint'):continue
        entry=dict(n)
        entry['degree_true']=phantom.graph.degree[n['true_node']] if n['true_node'] is not None else None
        entry['degree_recovered']=network.logical.degree[n['recovered_node']] if n['recovered_node'] is not None else None
        entry['degree_correct']=n['status']=='matched' and entry['degree_true']==entry['degree_recovered']
        if n['node_kind']=='junction':
            entry.update(junction_degree_true=entry['degree_true'],junction_degree_recovered=entry['degree_recovered'],junction_degree_correct=entry['degree_correct'])
        junctions.append(entry)
    for kind in ('junction','endpoint'):
        subset=[n for n in junctions if n['node_kind']==kind];matched=[n for n in subset if n['status']=='matched']
        errors=[n['error_mm'] for n in matched]
        row[kind+'_position_error_mm']=float(np.mean(errors)) if errors else None
        row['max_'+kind+'_position_error_mm']=max(errors) if errors else None
        row['matched_'+kind+'_count']=len(matched)
        row['missed_'+kind+'s']=sum(n['status']=='missed' for n in subset)
        row['spurious_'+kind+'s']=sum(n['status']=='spurious_or_unmatched' for n in subset)
        row[kind+'_degree_correct']=all(n['degree_correct'] for n in subset)
    cycles=[];inverse={v:k for k,v in bmap.items()}
    for i,cycle in enumerate(truth['cycle_basis']):
        present=all(k in inverse for k in cycle)
        if present:
            g=nx.MultiGraph()
            for key in cycle:
                b=network.branches[inverse[key]];g.add_edge(b['start'],b['end'])
            present=nx.is_connected(g) and all(d==2 for _,d in g.degree)
        exact=sum(phantom.branches[k]['length_mm'] for k in cycle)
        length=sum(network.branches[inverse[k]]['length_mm'] for k in cycle) if present else None
        cycles.append(dict(cycle_id=i,true_members=';'.join(cycle),recovered_members=';'.join(inverse[k] for k in cycle if k in inverse),membership_correct=present,cycle_preserved=present,true_perimeter_mm=exact,recovered_perimeter_mm=length,perimeter_error_mm=abs(length-exact) if present else None))
    for entry in branches:
        key,target=entry['recovered_branch'],entry['true_branch']
        entry['recovered_length_mm']=network.branches[key]['length_mm'] if key is not None else None
        entry['true_length_mm']=phantom.branches[target]['length_mm'] if target is not None else None
        entry['assignment_correct']=entry['status']=='matched'
        entry['cycle_membership_correct']=entry['status']=='matched' and all(c['membership_correct'] for c in cycles if target in c['true_members'].split(';'))
        if entry['status']=='matched':
            entry['signed_length_bias_mm']=entry['recovered_length_mm']-entry['true_length_mm']
            entry['percent_length_error']=100*abs(entry['signed_length_bias_mm'])/entry['true_length_mm']
    row['uncertified_recovered_basis_cycles']=sum(not all(k in bmap for k in c) for c in network.cycles)
    row['number_of_cycles_found']=sum(c['cycle_preserved'] for c in cycles)
    row['number_of_cycles_missed']=truth['true_cycle_rank']-row['number_of_cycles_found']
    row['communicating_bridge_recovered']=all(k in inverse for k,b in phantom.branches.items() if b['communicating_bridge']) if any(b['communicating_bridge'] for b in phantom.branches.values()) else None
    row['topology_correct']=all(row[k] for k in ('component_count_correct','endpoint_count_correct','junction_count_correct','branch_count_correct','cycle_rank_correct','branch_connectivity_correct'))
    row['central_junction_structure_correct']=row['junction_count_correct'] and row['junction_degree_correct'] and row['branch_connectivity_correct']
    row['signed_length_bias_mm']=row['total_network_length_mm']-row['true_network_length_mm']
    row['percent_length_error']=100*abs(row['signed_length_bias_mm'])/row['true_network_length_mm']
    row['topology_failure_reason']=';'.join(k for k in ('component_count_correct','endpoint_count_correct','junction_count_correct','branch_count_correct','cycle_rank_correct','branch_connectivity_correct') if not row[k])
    return row,branches,cycles,junctions,qc


def summarize(rows):
    valid=[r for r in rows if r['extraction_success']]
    out={'runs':len(rows),'extraction_success_rate':len(valid)/len(rows) if rows else None,'refinement_acceptance_rate':sum(r.get('refinement_accepted',True) for r in rows)/len(rows) if rows else None}
    accuracy={'component_accuracy':'component_count_correct','endpoint_count_accuracy':'endpoint_count_correct','junction_count_accuracy':'junction_count_correct','branch_count_accuracy':'branch_count_correct','cycle_rank_accuracy':'cycle_rank_correct','strict_connectivity_accuracy':'branch_connectivity_correct','strict_topology_accuracy':'topology_correct'}
    for name,field in accuracy.items():out[name]=sum(bool(r.get(field)) for r in rows)/len(rows) if rows else None
    means={'mean_network_coverage':'network_coverage_fraction','mean_network_position_error_mm':'mean_network_position_error_mm','mean_network_length_error_mm':'network_length_error_mm','mean_junction_error_mm':'junction_position_error_mm','mean_endpoint_error_mm':'endpoint_position_error_mm','branch_precision':'branch_precision','branch_recall':'branch_recall','branch_f1':'branch_f1','adjacency_precision':'adjacency_precision','adjacency_recall':'adjacency_recall','adjacency_f1':'adjacency_f1'}
    for name,field in means.items():
        values=[r[field] for r in valid if r.get(field) is not None];out[name]=float(np.mean(values)) if values else None
    out['failure_inclusive_coverage']=sum(r.get('network_coverage_fraction',0) for r in valid)/len(rows) if rows else None
    return out
