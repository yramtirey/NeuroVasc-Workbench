"""Synthetic hybrid centerline validation; existing archives and cases are read-only."""
import argparse
from collections import Counter
from dataclasses import asdict, replace
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'validation/outputs/.matplotlib'))
import matplotlib
matplotlib.use('Agg')
import nibabel as nib
import numpy as np
from neurovasc.geometry.centerline_qc import Recovery, QC_THRESHOLDS
from neurovasc.geometry.centerline_refinement import Refiner, RefinementParameters
from neurovasc.geometry.hybrid_centerline import recover, hybrid, CANDIDATE
from validation.synthetic_centerline.experiments_v2 import matrix_v2
from validation.synthetic_centerline.comparison_v2 import METHODS,measure,aggregate
from scripts.run_caliber_validation import write_csv
from scripts.run_caliber_validation_v2 import clean_json


def configurations():
    params=[RefinementParameters(recenter=center,spacing_mm=step,smoothing_scale_mm=scale) for center in (False,True) for step in (.5,1.) for scale in (0.,1.,2.)]
    params += [replace(CANDIDATE,refine_endpoints=False)]
    params += [RefinementParameters(recenter=False,spacing_mm=step,smoothing_scale_mm=0.,refine_endpoints=False) for step in (.5,1.)]
    return {f'center{int(p.recenter)}_step{p.spacing_mm:g}_smooth{p.smoothing_scale_mm:g}_end{int(p.refine_endpoints)}':p for p in params}


def output_directory(path):
    path=Path(path).resolve()
    allowed=(ROOT/'validation/outputs/centerline_v2').resolve()
    if not path.is_relative_to(allowed):
        raise ValueError('Choose centerline_v2 or a descendant; earlier archives are protected')
    return path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'validation/outputs/centerline_v2')
    args=parser.parse_args()
    try:output=output_directory(args.output)
    except ValueError as error:parser.error(str(error))
    output.mkdir(parents=True,exist_ok=True)
    (output/'masks').mkdir(exist_ok=True)
    with (ROOT/'validation/outputs/caliber_v2/results.csv').open() as stream:
        known={r['experiment_id'] for r in csv.DictReader(stream) if r['status']=='upstream_failed'}
    variants=configurations()
    candidate=next(k for k,p in variants.items() if p==CANDIDATE)
    rows,profiles,qcs,sweep,branches,loops,geometries,examples=[],[],[],[],[],[],[],{}
    for eid,experiment in matrix_v2():
        vessel=experiment.vessel
        meta={'geometry_id':eid,'geometry':vessel.geometry,'cohort':experiment.cohort,
              'spacing_z_mm':float(vessel.volume.spacing[2]),'spacing_x_mm':float(vessel.volume.spacing[0]),'spacing_y_mm':float(vessel.volume.spacing[1]),
              'phase_id':experiment.phase_id,'orientation_id':experiment.orientation_id,'known_v1_failure':eid in known,
              'diameter_mm':vessel.proximal_diameter_mm,'is_loop':vessel.geometry.startswith('loop')}
        geometries.append({**meta,**vessel.metadata(),'affine':vessel.volume.affine.tolist(),'shape':list(vessel.volume.shape)})
        nib.save(nib.Nifti1Image(vessel.volume.data,vessel.volume.affine),output/'masks'/f'{eid}.nii.gz')
        lee= recover(vessel.volume,'lee')
        raw= recover(vessel.volume,'geodesic_raw')
        refiner=Refiner(vessel.volume)
        refined_results={}
        for variant,parameters in variants.items():
            try:result=refiner.refine(raw,parameters) if raw.path is not None else Recovery(error=raw.error,supports_cycles=False)
            except ValueError as error:result=Recovery(error=str(error),supports_cycles=False)
            metrics,_,_,_,_=measure(vessel,result)
            sweep.append({**meta,'variant':variant,**asdict(parameters),**metrics})
            if variant==candidate:refined_results[variant]=result
        refined=refined_results[candidate]
        selected,diagnostics,_=hybrid(vessel.volume,lee_result=lee,raw_result=raw,refined_result=refined)
        recoveries=dict(zip(METHODS,(lee,raw,refined,selected)))
        examples[eid]={'vessel':vessel,'recoveries':recoveries,**meta}
        for method,recovery in recoveries.items():
            metrics,samples,branch_rows,loop,qc=measure(vessel,recovery)
            key={**meta,'estimator_name':method}
            rows.append({**key,**metrics,**(diagnostics if method=='hybrid' else {}),
                         'refinement_variant':candidate if method=='geodesic_refined' or method=='hybrid' and diagnostics['fallback_used'] else None})
            profiles.extend({**key,**s} for s in samples)
            qcs.extend({**key,**c} for c in qc['checks'])
            branches.extend({**key,**r} for r in branch_rows)
            if loop is not None:loops.append({**key,**loop})
        print(f"{eid}: {diagnostics['selected_method']}; final QC={diagnostics['final_qc_pass']}",flush=True)
    pairs=[]
    for eid in examples:
        group={r['estimator_name']:r for r in rows if r['geometry_id']==eid}
        for comparator in METHODS[1:]:
            a,b=group['lee'],group[comparator]
            pair={'geometry_id':eid,'comparator':comparator,'lee_success':a['extraction_success'],'comparator_success':b['extraction_success'],'known_v1_failure':a['known_v1_failure']}
            for key in ('mean_position_error_mm','coverage_fraction','length_error_mm','mean_tangent_error_deg','mean_endpoint_error_mm','mean_turn_deg'):
                pair[f'{key}_minus_lee']=b[key]-a[key] if a.get(key) is not None and b.get(key) is not None else None
            pairs.append(pair)
    nonloop=[r for r in rows if not r['is_loop']]
    def grouped(selected):return {m:aggregate([r for r in selected if r['estimator_name']==m]) for m in METHODS}
    common={r['geometry_id'] for r in nonloop if r['estimator_name']=='lee' and r['extraction_success']}
    hybrid_rows=[r for r in nonloop if r['estimator_name']=='hybrid']
    summary={'version':'centerline_v2','geometry_count':len(examples),'run_count':len(rows),'sweep_run_count':len(sweep),
             'headline_scope':'60 unchanged v1 geometries; 4 loop networks reported separately',**grouped(nonloop),
             'all_geometries_including_loops':grouped(rows),
             'shared_lee_success':grouped([r for r in nonloop if r['geometry_id'] in common]),
             'known_failures':{'count':len(known),'hybrid_recovered':sum(r['extraction_success'] for r in hybrid_rows if r['known_v1_failure']),
                               **grouped([r for r in nonloop if r['known_v1_failure']])},
             'selection_counts':dict(Counter(r['selected_method'] for r in hybrid_rows)),
             'fallback_rate':sum(r['fallback_used'] for r in hybrid_rows)/len(hybrid_rows),
             'refinement_sweep':{k:aggregate([r for r in sweep if r['variant']==k and not r['is_loop']]) for k in variants},
             'parameters':{'candidate':candidate,'configurations':{k:asdict(p) for k,p in variants.items()},'qc_thresholds':QC_THRESHOLDS,'qc_mode':'network','production_eligible':False,'fallback_supports_cycles':False},
             'topology':{'loop_results':loops,'y_results':[{k:r.get(k) for k in ('geometry_id','estimator_name','topology_preserved','junction_error_mm','endpoint_count','junction_count','cycle_count')} for r in rows if r['geometry']=='bifurcation']},
             'dependencies':{p:importlib.metadata.version(p) for p in ('numpy','scipy','scikit-image','networkx','nibabel','matplotlib')},'python':platform.python_version()}
    for field in ('phase_id','spacing_z_mm','orientation_id','geometry'):
        summary[f'by_{field}']={str(value):grouped([r for r in nonloop if r[field]==value]) for value in sorted(set(r[field] for r in nonloop))}
    sources=[*ROOT.glob('neurovasc/**/*.py'),*ROOT.glob('validation/synthetic_centerline/*.py'),*ROOT.glob('validation/synthetic_caliber/*.py'),*ROOT.glob('scripts/run_*validation*.py')]
    summary['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)}
    for name,data in [('results',rows),('matched_comparison',pairs),('profiles',profiles),('qc_results',qcs),('refinement_sweep',sweep),('branch_results',branches),('loop_results',loops)]:
        write_csv(output/f'{name}.csv',data)
    for name,data in [('summary',summary),('geometries',geometries)]:
        (output/f'{name}.json').write_text(json.dumps(clean_json(data),indent=2,allow_nan=False)+'\n')
    from validation.synthetic_centerline.figures_v2 import create_figures
    create_figures(output,rows,sweep,examples)
    files=[str(p.relative_to(output)) for p in output.rglob('*') if p.is_file() and 'regressions' not in p.relative_to(output).parts and p.name!='manifest.txt']
    (output/'manifest.txt').write_text('\n'.join(sorted(files+['manifest.txt']))+'\n')
    print(json.dumps({k:summary[k] for k in (*METHODS,'selection_counts','fallback_rate')},indent=2))

if __name__=='__main__':main()
