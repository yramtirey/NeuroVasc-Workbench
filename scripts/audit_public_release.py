"""Audit staged public files without printing suspected secret values."""
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT=Path(__file__).resolve().parents[1]


def audit():
    names=subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().split('\0')
    violations=[];files=[]
    private_parts={'.venv','node_modules','__pycache__','.local-release','.pytest_cache','.git','dist'}
    secrets=[re.compile(rb'gh[pousr]_[A-Za-z0-9]{20,}'),re.compile(rb'github_pat_[A-Za-z0-9_]{30,}'),re.compile(rb'AKIA[0-9A-Z]{16}')]
    for name in filter(None,names):
        p=Path(name)
        # Inspect index bytes: this audits exactly what would be committed.
        data=subprocess.check_output(['git','show',':'+name],cwd=ROOT)
        reasons=[]
        if private_parts.intersection(p.parts) or name.startswith(('outputs/','validation/outputs/')):reasons.append('private/generated path')
        if name.endswith(('.nii','.nii.gz','.dcm','.zip','.tar','.tar.gz','.pyc','.DS_Store')):reasons.append('excluded file type')
        if len(data)>5_000_000:reasons.append('file exceeds 5 MB public limit')
        if b'\0' not in data:
            if b'/' + b'Users/' in data or b'/' + b'private/tmp/' in data:reasons.append('local absolute path')
            if any(rule.search(data) for rule in secrets):reasons.append('credential-like value')
            if b'-----BEGIN '+b'PRIVATE KEY-----' in data or b'-----BEGIN '+b'RSA PRIVATE KEY-----' in data:reasons.append('private key')
        if reasons:violations.append({'path':name,'reasons':reasons})
        files.append({'path':name,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
    return {'files':files,'total_bytes':sum(f['bytes'] for f in files),'violations':violations,'passed':bool(files) and not violations}


if __name__=='__main__':
    result=audit()
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['passed'] else 1)
