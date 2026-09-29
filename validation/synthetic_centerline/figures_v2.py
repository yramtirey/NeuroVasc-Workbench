"""Single-panel Matplotlib diagnostics; raw and failed results remain visible."""
import matplotlib.pyplot as plt
import numpy as np
from .comparison_v2 import METHODS, aggregate
from .experiments_v2 import loop_samples
from .truth import sample_truth


def create_figures(output, rows, sweep, examples):
    folder=output/'figures'
    folder.mkdir(exist_ok=True)
    nonloop=[r for r in rows if not r['is_loop']]
    def save(name):
        plt.tight_layout();plt.savefig(folder/f'{name}.png',dpi=150);plt.close()
    for name,key,label in [('A_success','success_rate','Extraction fraction'),('B_qc','qc_pass_rate','QC pass fraction'),
                           ('C_position','mean_position_error_mm','Position error (mm)'),('D_tangent','mean_tangent_error_deg','Tangent error (degrees)'),
                           ('E_coverage','coverage_fraction','Coverage fraction'),('F_length','length_error_mm','Absolute length error (mm)'),
                           ('M_endpoints','mean_endpoint_error_mm','Endpoint error (mm)'),('N_smoothness','mean_turn_deg','Mean local turn (degrees)')]:
        plt.figure(figsize=(8,4.5))
        for i,m in enumerate(METHODS):
            group=[r for r in nonloop if r['estimator_name']==m]
            a=aggregate(group)
            plt.bar(i,a[key],color=f'C{i}')
            plt.text(i,a[key],f"{a['success_count']}/60",ha='center',va='bottom',fontsize=9)
        plt.xticks(range(4),METHODS,rotation=12)
        plt.ylabel(label);plt.title('V1 cohort: errors on own successful support; labels show extraction counts')
        save(name)
    for name,key in [('G_raw_refined_position','mean_position_error_mm'),('H_raw_refined_tangent','mean_tangent_error_deg')]:
        plt.figure(figsize=(5.5,5.5))
        pairs=[]
        for eid in examples:
            group={r['estimator_name']:r for r in nonloop if r['geometry_id']==eid}
            if group and all(group[m].get(key) is not None for m in ('geodesic_raw','geodesic_refined')):
                pairs.append((group['geodesic_raw'][key],group['geodesic_refined'][key]))
        values=np.array(pairs)
        plt.scatter(*values.T);maximum=values.max()*1.05
        plt.plot([0,maximum],[0,maximum],'--',label='equal error')
        plt.xlabel('Raw '+key);plt.ylabel('Refined '+key);plt.legend();save(name)
    plt.figure(figsize=(7,4))
    displacements=[d for x in examples.values() if not x['is_loop'] for d in x['recoveries']['geodesic_refined'].metadata.get('displacements_mm',[])]
    plt.hist(displacements,bins=35)
    plt.xlabel('Displacement at normalized arc correspondence (mm)');plt.ylabel('Refined branch sample count');save('I_displacements')
    for name,field in [('J_phase','phase_id'),('K_spacing','spacing_z_mm'),('L_orientation','orientation_id')]:
        selected=[r for r in nonloop if r['geometry']==('straight' if field=='phase_id' else 'oblique')] if field!='spacing_z_mm' else nonloop
        categories=sorted(set(r[field] for r in selected))
        plt.figure(figsize=(7,4.5))
        for m in METHODS:
            values=[aggregate([r for r in selected if r['estimator_name']==m and r[field]==v])['mean_position_error_mm'] for v in categories]
            plt.plot(range(len(categories)),[np.nan if x is None else x for x in values],'o-',label=m)
        plt.xticks(range(len(categories)),categories);plt.xlabel(field);plt.ylabel('Mean position error (mm)')
        plt.title('Successful support only; missing Lee half-phase error is not zero');plt.legend();save(name)
    plt.figure(figsize=(7,4))
    hybrid=[r for r in nonloop if r['estimator_name']=='hybrid']
    plt.bar(['Lee retained','Refined fallback'],[sum(not r['fallback_used'] for r in hybrid),sum(r['fallback_used'] for r in hybrid)])
    plt.ylabel('Geometries');save('R_selection')
    for metric in ('mean_position_error_mm','mean_tangent_error_deg'):
        plt.figure(figsize=(10,7))
        variants=list(dict.fromkeys(r['variant'] for r in sweep))
        values=[aggregate([r for r in sweep if r['variant']==v and not r['is_loop']])[metric] for v in variants]
        plt.barh(variants,values);plt.xlabel(metric);plt.title('Declared refinement sweep; every variant retained');save('sweep_'+metric)
    for eid,item in examples.items():
        v=item['vessel']
        if not (item['known_v1_failure'] or item['is_loop'] or v.geometry=='bifurcation'):
            continue
        plt.figure(figsize=(8,6))
        if item['known_v1_failure']:
            # All diagnostic phantoms use diagonal affines; show foreground
            # projected over y, explicitly labeling this as an x-z projection.
            mask=np.max(v.volume.data,axis=1).T
            a=v.volume.affine
            extent=[a[0,3]-.5*v.volume.spacing[0],a[0,3]+(mask.shape[1]-.5)*v.volume.spacing[0],a[2,3]-.5*v.volume.spacing[2],a[2,3]+(mask.shape[0]-.5)*v.volume.spacing[2]]
            plt.imshow(mask,origin='lower',extent=extent,cmap='Greys',alpha=.3)
        truths=loop_samples(v) if item['is_loop'] else [sample_truth(v,n) for n in v.branches] if v.geometry=='bifurcation' else [sample_truth(v)]
        for i,t in enumerate(truths):plt.plot(t[:,0],t[:,2],'k--',linewidth=1,label='analytic truth' if i==0 else None)
        for index,(method,recovery) in enumerate(item['recoveries'].items()):
            if recovery.graph is None:
                plt.plot([],[],color=f'C{index}',label=f'{method}: failed')
                continue
            graph=recovery.graph
            for i,(a,b) in enumerate(graph.edges):
                line=np.array([graph.nodes[a]['world'],graph.nodes[b]['world']])
                plt.plot(line[:,0],line[:,2],color=f'C{index}',alpha=.8,linewidth=1.4 if method!='hybrid' else .8,
                         linestyle=':' if method=='hybrid' else '-',label=method if i==0 else None)
        plt.xlabel('x (mm)');plt.ylabel('z (mm)');plt.axis('equal');plt.title(eid+' — x-z projection',fontsize=10);plt.legend(fontsize=8)
        save(('O_' if item['known_v1_failure'] else 'Q_' if item['is_loop'] else 'P_')+eid)
