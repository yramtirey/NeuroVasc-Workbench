"""Separate Matplotlib figures for fixed-topology refinement evidence."""
import matplotlib.pyplot as plt
import numpy as np
from .evaluate_v2 import METHODS, summarize


def figures(output, rows, invariants, spurs):
    folder=output/'figures';folder.mkdir(exist_ok=True)
    core=[r for r in rows if r['cohort']=='core_v1']
    for title,field in [('A_strict_topology','strict_topology_accuracy'),('B_cycle_rank','cycle_rank_accuracy'),('C_junction_count','junction_count_accuracy'),('D_branch_F1','branch_f1'),('E_adjacency_F1','adjacency_f1'),('F_length_error_mm','mean_network_length_error_mm'),('G_junction_error_mm','mean_junction_error_mm'),('H_endpoint_error_mm','mean_endpoint_error_mm'),('I_coverage','mean_network_coverage'),('R_network_position_mm','mean_network_position_error_mm')]:
        values=[summarize([r for r in core if r['method']==m])[field] for m in METHODS]
        plt.figure(figsize=(9,4));plt.bar(range(4),values,color=[f'C{i}' for i in range(4)])
        plt.xticks(range(4),METHODS,fontsize=8);plt.ylabel(field);plt.title('Exact 76-geometry v1 cohort; rejected refinements retain labelled originals')
        plt.tight_layout();plt.savefig(folder/f'{title}.png',dpi=130);plt.close()
    for title,field,metric in [('J_phase','phase','topology_correct'),('K_rotation','orientation','network_length_error_mm'),('P_small_bridge','bridge_diameter_mm','communicating_bridge_recovered'),('Q_near_touch','minimum_separation_mm','false_connection_count')]:
        chosen=[r for r in core if r.get(field) is not None and (field!='phase' or r['orientation']=='canonical')]
        categories=sorted({r[field] for r in chosen});plt.figure(figsize=(8,4))
        for m in METHODS:
            values=[]
            for category in categories:
                group=[r for r in chosen if r['method']==m and r[field]==category and r.get(metric) is not None]
                values.append(np.mean([bool(r[metric]) if metric=='false_connection_count' else r[metric] for r in group]) if group else np.nan)
            plt.plot(range(len(categories)),values,'o-',label=m)
        plt.xticks(range(len(categories)),categories);plt.xlabel(field);plt.ylabel(metric);plt.legend(fontsize=8);plt.tight_layout();plt.savefig(folder/f'{title}.png',dpi=130);plt.close()
    plt.figure(figsize=(8,4))
    for i,m in enumerate(METHODS[2:]):
        group=[r for r in invariants if r['method']==m]
        plt.bar(i,np.mean([r['topology_invariant_pass'] for r in group]) if group else 0,color=f'C{i}',label=m)
    plt.xticks(range(2),METHODS[2:],fontsize=8);plt.ylim(0,1.05);plt.ylabel('Proposed invariant pass fraction');plt.tight_layout();plt.savefig(folder/'S_invariance.png',dpi=130);plt.close()
    plt.figure(figsize=(7,4))
    group=[s for s in spurs if s['method']=='topology_v2_refined']
    plt.bar(['terminal candidates','marked','removed'],[len(group),sum(s['spur_candidate'] for s in group),sum(s['spur_removed'] for s in group)])
    plt.ylabel('Branches across full matrix');plt.title('Conservative spur classification: mark only');plt.tight_layout();plt.savefig(folder/'O_spur_decisions.png',dpi=130);plt.close()


def overlay(output,eid,phantom,before,after,accepted):
    folder=output/'figures'/'examples';folder.mkdir(parents=True,exist_ok=True)
    plt.figure(figsize=(7,6));labels=set()
    for network,color,label in ((None,'C2','analytic truth'),(before,'C1','v1 raw'),(after,'C0','v2 accepted' if accepted else 'v2 rejected; original retained')):
        branches=phantom.branches if network is None else network.branches
        for b in branches.values():
            p=b['points'];plt.plot(p[:,0],p[:,2],'--' if network is None else '-',color=color,linewidth=1.4,label=label if label not in labels else None);labels.add(label)
    plt.xlabel('x (mm), projected');plt.ylabel('z (mm), projected');plt.axis('equal');plt.title(eid,fontsize=9);plt.legend(fontsize=8);plt.tight_layout();plt.savefig(folder/f'{eid}.png',dpi=120);plt.close()
