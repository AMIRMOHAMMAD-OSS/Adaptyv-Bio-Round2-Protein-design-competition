import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

def execute(argv):
    p = subprocess.run(argv, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode:
        raise RuntimeError(f'Failed: {argv}\n{p.stderr[-2000:]}')
    return p.stdout.strip()

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(2 ** 20), b''):
            h.update(b)
    return h.hexdigest()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', required=True)
    p.add_argument('--rfd-env', default='binder-rfd-wsl')
    p.add_argument('--mpnn-env', default='binder-mpnn-wsl')
    a = p.parse_args()
    record = {'python': sys.version, 'repos': {}, 'weights': {}, 'environments': {}}
    for name, var in [('RFdiffusion', 'RFD_ROOT'), ('ProteinMPNN', 'MPNN_ROOT')]:
        path = Path(os.environ[var]).resolve()
        record['repos'][name] = {'path': str(path), 'commit': execute(['git', '-C', str(path), 'rev-parse', 'HEAD']), 'working_tree': execute(['git', '-C', str(path), 'status', '--porcelain', '--untracked-files=no'])}
        for w in sorted(path.glob('models/*.pt')) + sorted(path.glob('vanilla_model_weights/*.pt')):
            print('Hashing', w, flush=True)
            record['weights'][str(w)] = sha(w)
    for env in [a.rfd_env, a.mpnn_env]:
        record['environments'][env] = execute(['conda', 'run', '-n', env, 'python', '-m', 'pip', 'freeze'])
    cf = Path(os.environ['COLABFOLD_BATCH']).resolve()
    cf_python = cf.parent / 'python'
    record['environments']['colabfold'] = execute([str(cf_python), '-m', 'pip', 'freeze'])
    paramroot = Path(os.environ.get('COLABFOLD_PARAMS', str(Path.home() / '.cache/colabfold')))
    params = sorted(paramroot.glob('**/params_model_*multimer*.npz'))
    if not params:
        raise RuntimeError(f'No AlphaFold multimer weights found in {paramroot}; predownload weights and/or set COLABFOLD_PARAMS')
    for w in params:
        print('Hashing', w, flush=True)
        record['weights'][str(w)] = sha(w)
    out = Path(a.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=2) + '\n')
    print(out)
if __name__ == '__main__':
    main()
