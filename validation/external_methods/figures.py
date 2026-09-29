"""Separate Matplotlib figures; unavailable measurements are never zero bars."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .evaluate import DIRECT


def label(row):
    return row['method']+(' oracle' if row.get('oracle_seed') else ' automatic' if row['method']=='vmtk' else '')


def figures(output,centers,calibers,branches,examples,lookup):
    directory=output/'figures';directory.mkdir(exist_ok=True)
    def save(name,title,ylabel,values):
        groups={}
        for key,value in values:
            if value is not None and np.isfinite(value):groups.setdefault(key,[]).append(value)
        horizontal=len(groups)>8
        fig,ax=plt.subplots(figsize=(10,max(5,.4*len(groups))) if horizontal else (9,5))
        if groups:
            if horizontal:
                ax.barh(range(len(groups)),[np.mean(v) for v in groups.values()])
                ax.set_yticks(range(len(groups)),[f'{k} (n={len(v)})' for k,v in groups.items()])
                ax.invert_yaxis()
            else:
                ax.bar(list(groups),[np.mean(v) for v in groups.values()])
                ax.set_xticks(range(len(groups)),[f'{k}\n(n={len(v)})' for k,v in groups.items()],rotation=35,ha='right')
        else:ax.text(.5,.5,'No executed comparable measurements',ha='center',transform=ax.transAxes)
        ax.set_title(title)
        (ax.set_xlabel if horizontal else ax.set_ylabel)(ylabel)
        fig.tight_layout();fig.savefig(directory/name,dpi=140);plt.close(fig)
    direct=[r for r in centers if r['comparison_scope']==DIRECT]
    good=[r for r in direct if r.get('extraction_success')]
    for name,field,title,unit in [
        ('01_position.png','mean_position_error_mm','Centerline position: successful direct runs','Error (mm)'),
        ('02_coverage.png','coverage_fraction','Centerline coverage: successful direct runs','Fraction'),
        ('03_length.png','length_error_mm','Unique physical edge length error','Absolute error (mm)')]:
        save(name,title,unit,[(label(r),r.get(field)) for r in good])
    for number,quantity in [(4,'radius'),(5,'diameter')]:
        save(f'{number:02d}_{quantity}.png',f'{quantity.capitalize()} MAE: own valid station support','MAE (mm)',
             [(label(r),r.get('mae_mm')) for r in calibers if r['quantity']==quantity and r['comparison_scope']==DIRECT])
    fig,ax=plt.subplots(figsize=(6,6))
    for method in ('vmtk_oracle','vmtk_automatic'):
        a={r['geometry_id']:r for r in calibers if r['method']=='edt' and r['quantity']=='radius' and r.get('estimated_mean_mm') is not None}
        b={r['geometry_id']:r for r in calibers if r['method']==method and r['quantity']=='radius' and r.get('estimated_mean_mm') is not None}
        ids=sorted(a.keys()&b.keys());ax.scatter([a[i]['estimated_mean_mm'] for i in ids],[b[i]['estimated_mean_mm'] for i in ids],label=f'{method} (n={len(ids)})',alpha=.65)
    ax.set(xlabel='EDT mean radius (mm)',ylabel='VMTK MIS mean radius (mm)',title='Per-run means; own support, not pointwise agreement');ax.legend();fig.tight_layout();fig.savefig(directory/'06_radius_scatter.png',dpi=140);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4));ax.axis('off');ax.text(.5,.5,'External SlicerVMTK CE not executed\nNo CE comparison values generated',ha='center',va='center');fig.savefig(directory/'07_external_ce_unavailable.png',dpi=140);plt.close(fig)
    save('08_orientation.png','Oblique centerline error by orientation (successful runs)','Error (mm)',[(label(r)+' '+r['orientation'],r.get('mean_position_error_mm')) for r in good if r['geometry']=='oblique'])
    save('09_spacing.png','Centerline error by spacing (successful runs)','Error (mm)',[(label(r)+f" z={r['spacing_z_mm']:g}",r.get('mean_position_error_mm')) for r in good])
    for gid,(vessel,_) in lookup.items():
        if vessel.geometry!='bifurcation':continue
        fig,ax=plt.subplots(figsize=(8,6))
        for segment in vessel.branches.values():ax.plot(segment[:,0],segment[:,2],color='black',linewidth=3,alpha=.3)
        colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
        selected=[(key,example) for key,example in examples.items() if key.split('|')[0]==gid]
        for index,(key,example) in enumerate(selected):
            for i,path in enumerate(example['paths']):
                p=np.array(path);ax.plot(p[:,0],p[:,2],color=colors[index%len(colors)],label=key.split('|')[1] if i==0 else None,alpha=.7)
        ax.set(xlabel='World x (mm)',ylabel='World z (mm)',title=f'Y network x-z projection: {gid}');ax.axis('equal');ax.legend(fontsize=7);fig.tight_layout();fig.savefig(directory/f'10_y_{gid}.png',dpi=140);plt.close(fig)
    save('11_junction.png','Y localization: graph junction / VMTK reference origin','Error (mm)',[(label(r),r.get('vmtk_reference_junction_error_mm') if r['method']=='vmtk' else r.get('junction_error_mm')) for r in branches])
    save('12_success.png','Extraction success: direct attempted runs (unavailable excluded)','Fraction',[(label(r),int(bool(r.get('extraction_success')))) for r in direct if r['status']!='unavailable'])
    save('13_seeds.png','VMTK seed strategy: successful direct position error','Error (mm)',[(label(r),r.get('mean_position_error_mm')) for r in good if r['method']=='vmtk'])
