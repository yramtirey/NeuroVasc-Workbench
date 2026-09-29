"""Typed empty results, failure-inclusive summaries and safe artifact paths."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
ALLOWED=ROOT/'validation/outputs/external_methods_v1'


def output_directory(path):
    path=Path(path).absolute()
    # Check against a lexical trusted root; resolving the root itself would
    # accidentally authorize a symlink to an old archive.
    if not path.resolve().is_relative_to(ALLOWED.absolute()):
        raise ValueError('output must remain inside validation/outputs/external_methods_v1')
    for parent in [path,*path.parents]:
        if parent==ROOT: break
        if parent.is_symlink(): raise ValueError('symlink output path rejected')
    return path.resolve()


def clean(value):
    if isinstance(value,dict): return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple,np.ndarray)): return [clean(v) for v in value]
    if isinstance(value,np.generic): value=value.item()
    if isinstance(value,float) and not np.isfinite(value): return None
    return value


def json_write(path,value):
    Path(path).write_text(json.dumps(clean(value),indent=2,allow_nan=False)+'\n')


def csv_write(path,rows,required=()):
    fields=list(dict.fromkeys([*required,*[k for r in rows for k in r]]))
    with Path(path).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(clean(v)) if isinstance(v,(dict,list,tuple,np.ndarray)) else clean(v) for k,v in row.items()})


def aggregates(rows):
    groups={}
    for row in rows:
        if row['comparison_scope']=='METHOD_ASSUMPTION_MISMATCH': continue
        key=f"{row['method']}|oracle={row.get('oracle_seed')}|{row['comparison_scope']}"
        groups.setdefault(key,[]).append(row)
    out={}
    for key,group in groups.items():
        successes=[r for r in group if r.get('extraction_success')]
        available=[r for r in group if r.get('status')!='unavailable']
        result=dict(requested_runs=len(group),attempted_runs=len(available),successful_runs=len(successes),
                    success_rate=len(successes)/len(available) if available else None,
                    successful_support_coverage=float(np.mean([r['coverage_fraction'] for r in successes])) if successes else None,
                    failure_inclusive_coverage=sum(r.get('coverage_fraction',0) for r in successes)/len(available) if available else None)
        for name in ('mean_position_error_mm','length_error_mm','endpoint_error_mm','mean_tangent_error_deg'):
            values=[r[name] for r in successes if r.get(name) is not None]
            result[name]=float(np.mean(values)) if values else None
        out[key]=result
    return out


def inventory(output):
    # Development probes and caches are not evidence from the declared matrix.
    excluded={'.matplotlib','__pycache__','smoke','development'}
    paths=sorted(p.relative_to(output).as_posix() for p in output.rglob('*')
                 if p.is_file() and p.name!='manifest.txt'
                 and not any(part in excluded or part.startswith('.') or part.startswith('integration_')
                             for part in p.relative_to(output).parts))
    (output/'manifest.txt').write_text('\n'.join(paths+['manifest.txt'])+'\n')
    return len(paths)+1


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
