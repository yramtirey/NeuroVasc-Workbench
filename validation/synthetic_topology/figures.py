"""Separate Matplotlib summary figures and explicitly marked topology failures."""
import matplotlib.pyplot as plt
import numpy as np
from .evaluate import METHODS,aggregate


def summaries(output,rows):
    folder=output/'figures';folder.mkdir(exist_ok=True)
    for name,key,label in [('A_cycles','cycle_rank_correct','Cycle-rank accuracy'),('B_endpoints','endpoint_count_correct','Endpoint-count accuracy'),
                           ('C_junctions','junction_count_correct','Junction-count accuracy'),('D_branches','branch_count_correct','Branch-count accuracy'),
                           ('E_connectivity','branch_connectivity_correct','Connectivity accuracy'),('F_coverage','network_coverage_fraction','Coverage on extracted networks'),
                           ('G_position','mean_network_position_error_mm','Network position error (mm)'),('H_length','network_length_error_mm','Network length error (mm)')]:
        plt.figure(figsize=(9,4.5))
        for i,m in enumerate(METHODS):
            a=aggregate([r for r in rows if r['method']==m]);plt.bar(i,a[key] or 0,color=f'C{i}')
            plt.text(i,a[key] or 0,f"{a['extraction_success_rate']:.0%} extracted",ha='center',va='bottom',fontsize=8)
        plt.xticks(range(5),METHODS,rotation=12);plt.ylabel(label);plt.tight_layout();plt.savefig(folder/f'{name}.png',dpi=130);plt.close()
    for name,field,key in [('I_phase','phase','cycle_rank_correct'),('J_spacing','spacing_z_mm','cycle_rank_correct'),
                           ('K_bridge','bridge_diameter_mm','communicating_bridge_recovered'),('L_near_touch','minimum_separation_mm','false_connection')]:
        selected=[r for r in rows if r.get(field) is not None]
        categories=sorted(set(r[field] for r in selected));plt.figure(figsize=(8,4.5))
        for method in METHODS:
            values=[]
            for category in categories:
                group=[r for r in selected if r[field]==category and r['method']==method]
                if key=='false_connection':
                    valid=[r for r in group if r['extraction_success']]
                    values.append(np.mean([r['false_connection_count']>0 for r in valid]) if valid else np.nan)
                else:values.append(np.mean([bool(r.get(key,False)) for r in group]))
            plt.plot(range(len(categories)),values,'o-',label=method)
        plt.xticks(range(len(categories)),categories);plt.xlabel(field);plt.ylabel(key+' fraction')
        if key=='false_connection':plt.title('Extracted networks only; failed extractions are not counted as correct')
        plt.legend();plt.tight_layout();plt.savefig(folder/f'{name}.png',dpi=130);plt.close()


def diagnostic(output,eid,phantom,method,row,network,branch_rows):
    folder=output/'figures'/'diagnostics';folder.mkdir(parents=True,exist_ok=True)
    plt.figure(figsize=(7,6))
    matched_truth={r['true_branch'] for r in branch_rows if r['status']=='matched'}
    matched_actual={r['recovered_branch'] for r in branch_rows if r['status']=='matched'}
    labels=set()
    for key,b in phantom.branches.items():
        label='matched analytic branch' if key in matched_truth else 'missed/unmatched analytic branch'
        p=b['points'];plt.plot(p[:,0],p[:,2],'--',color='C2' if key in matched_truth else 'C1',linewidth=2,label=label if label not in labels else None);labels.add(label)
    if network is not None:
        for key,b in network.branches.items():
            label='matched recovery' if key in matched_actual else 'spurious/unmatched recovery'
            p=b['points'];plt.plot(p[:,0],p[:,2],color='C0' if key in matched_actual else 'C3',linewidth=1,label=label if label not in labels else None);labels.add(label)
    reason=row.get('topology_failure_reason') or row.get('failure_reason') or 'topology/connectivity match'
    plt.title(f'{eid} / {method}\n{reason}',fontsize=8,wrap=True)
    plt.xlabel('x (mm), projected');plt.ylabel('z (mm), projected');plt.axis('equal');plt.legend(fontsize=7)
    plt.tight_layout();plt.savefig(folder/f'{eid}_{method}.png',dpi=110);plt.close()
