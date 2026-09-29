"""Repeat actual external jobs; compare numeric outputs, preserving primary runs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
from validation.external_methods.report import ALLOWED,output_directory,json_write,inventory,sha
from scripts.run_external_methods_benchmark import probe


def compare(a,b,path=''):
    """Exact discrete identities; 1e-10 mm/array tolerance, never curve realignment."""
    if isinstance(a,dict):
        if not isinstance(b,dict) or a.keys()!=b.keys():return [path+':keys']
        return [issue for k in a for issue in compare(a[k],b[k],path+'/'+k)]
    if isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):return [path+':length']
        return [issue for i,(x,y) in enumerate(zip(a,b)) for issue in compare(x,y,path+f'/{i}')]
    if isinstance(a,(int,float)) and not isinstance(a,bool):
        return [] if isinstance(b,(int,float)) and np.isclose(a,b,rtol=0,atol=1e-10) else [path+':numeric']
    return [] if a==b else [path+':value']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ALLOWED)
    parser.add_argument('--external-python',type=Path,default=ROOT/'external/vmtk/.venv/bin/python')
    args=parser.parse_args();output=output_directory(args.output);python=args.external_python.absolute()
    env=probe(python)
    if not env.get('available'):raise SystemExit('Actual VMTK unavailable: '+env.get('reason','probe failed'))
    repeat=output/'repeatability';repeat.mkdir(exist_ok=True)
    def job(relative):
        path=output/relative;request=json.loads(path.read_text())
        if request['comparison_scope']=='METHOD_ASSUMPTION_MISMATCH':return None
        gid=request['geometry_id'];mode='oracle' if request['oracle_seed'] else 'automatic'
        original=output/'vmtk_outputs'/gid/mode/'result.json';a=json.loads(original.read_text())
        if a.get('request_sha256')!=sha(path) or a.get('surface_sha256')!=sha(request['surface']):
            raise ValueError('primary result does not match current input')
        for field in ('module_sha256','adapter_sha256'):
            if a.get('environment',{}).get(field)!=env.get(field):
                raise ValueError('repeatability requires identical adapter and external environment')
        dest=repeat/gid/mode;dest.mkdir(parents=True,exist_ok=True)
        process=subprocess.run([str(python),'-m','external.vmtk.adapter','--request',str(path),'--output',str(dest)],
            cwd=ROOT,capture_output=True,text=True,timeout=120)
        (dest/'stdout.log').write_text(process.stdout);(dest/'stderr.log').write_text(process.stderr)
        if process.returncode:raise RuntimeError(f'{gid}: external process failed: {process.stderr}')
        b=json.loads((dest/'result.json').read_text())
        fields=('status','extraction_success','failure_reason','settings','lines','branches','branch_failure_reason')
        issues=compare({k:a.get(k) for k in fields},{k:b.get(k) for k in fields})
        print('repeat',gid,mode,'equal' if not issues else 'DIFFERENT',flush=True)
        return dict(geometry_id=gid,oracle_seed=request['oracle_seed'],original_status=a['status'],
                    repeated_status=b['status'],equal=not issues,differences=issues,
                    primary_result_sha256=sha(original),repeat_result_sha256=sha(dest/'result.json'))
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows=[r for r in pool.map(job,json.loads((output/'requests.json').read_text())) if r is not None]
    record=dict(attempted=len(rows),equal=sum(r['equal'] for r in rows),atol=1e-10,rtol=0,
        comparator='Exact nested identities/shapes; numeric values within absolute tolerance; no registration or reordering',
        environment=env,runs=rows,source_sha256=sha(Path(__file__)))
    json_write(output/'repeatability.json',record);inventory(output)
    print(json.dumps({k:record[k] for k in ('attempted','equal','atol')}))
    if record['equal']!=record['attempted']:raise SystemExit(1)


if __name__=='__main__':main()
