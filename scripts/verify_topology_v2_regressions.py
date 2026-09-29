"""Rerun five archived milestones in a temporary source mirror, never in place."""
import csv
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
SUITES=(('caliber_v1','run_caliber_validation.py'),('caliber_v2','run_caliber_validation_v2.py'),
        ('centerline_v1','run_centerline_validation.py'),('centerline_v2','run_centerline_validation_v2.py'),('topology_v1','run_topology_validation.py'))


def semantic_rows(path):
    # V1 opaque recovered IDs can vary with Python hash iteration. Preserve all
    # anatomical assignments/status/metrics while ignoring those names and order.
    def canonical(row):
        result={}
        for key,value in row.items():
            if key in ('recovered_branch','recovered_node','recovered_members','cycle_id'):continue
            if key=='true_members':value=';'.join(sorted(value.split(';')))
            try:value=round(float(value),10)
            except ValueError:pass
            result[key]=value
        return json.dumps(result,sort_keys=True)
    with path.open() as stream:return Counter(canonical(row) for row in csv.DictReader(stream))


def main():
    destination=ROOT/'validation/outputs/topology_v2/regressions'
    destination.mkdir(parents=True,exist_ok=True)
    archived={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for suite,_ in SUITES for f in (ROOT/'validation/outputs'/suite).rglob('*') if f.is_file()}
    checks={};semantic={}
    with tempfile.TemporaryDirectory(prefix='neurovasc-topology-v2-') as folder:
        mirror=Path(folder)
        for name in ('neurovasc','scripts'):
            shutil.copytree(ROOT/name,mirror/name,ignore=shutil.ignore_patterns('__pycache__'))
        for source in (ROOT/'validation').glob('synthetic_*'):
            shutil.copytree(source,mirror/'validation'/source.name,ignore=shutil.ignore_patterns('__pycache__'))
        for suite,script in SUITES:
            with (destination/f'{suite}.log').open('w') as log:
                subprocess.run([sys.executable,str(mirror/'scripts'/script)],cwd=mirror,stdout=log,stderr=subprocess.STDOUT,check=True)
            generated=mirror/'validation/outputs'/suite
            for original in (ROOT/'validation/outputs'/suite).glob('*.csv'):
                new=generated/original.name
                key=f'{suite}/{original.name}'
                checks[key]=original.read_bytes()==new.read_bytes()
                if suite=='topology_v1':
                    semantic[key]=semantic_rows(original)==semantic_rows(new)
            shutil.copytree(generated,destination/suite,dirs_exist_ok=True)
            print(suite,'completed',flush=True)
    changed=[f for f,h in archived.items() if not (ROOT/f).exists() or hashlib.sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    report=dict(archived_files_checked=len(archived),changed_archive_files=changed,csv_byte_equality=checks,topology_v1_semantic_equality=semantic,semantic_convention='Ignore opaque recovered IDs and row order; retain truth assignments, counts, statuses and metrics rounded to 1e-10')
    (destination/'checks.json').write_text(json.dumps(report,indent=2)+'\n')
    files=sorted(str(f.relative_to(destination)) for f in destination.rglob('*') if f.is_file() and f.name!='manifest.txt')
    (destination/'manifest.txt').write_text('\n'.join(files+['manifest.txt'])+'\n')
    if changed or not all(value or semantic.get(key,False) for key,value in checks.items()):raise SystemExit('Regression verification failed; see checks.json')
    print(f'{sum(checks.values())}/{len(checks)} CSV files byte-identical; opaque-ID tables semantically equal; archives unchanged')

if __name__=='__main__':main()
