from __future__ import annotations
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(2 ** 20), b''):
            h.update(block)
    return h.hexdigest()

def save_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    os.replace(tmp, path)

def read_json(path):
    return json.loads(Path(path).read_text())

def csv_write(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else [])
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def fasta_read(path):
    records, header, sequence = ([], None, [])
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                records.append((header, ''.join(sequence)))
            header, sequence = (line[1:], [])
        else:
            if header is None:
                raise ValueError('FASTA sequence precedes its header')
            sequence.append(line)
    if header is not None:
        records.append((header, ''.join(sequence)))
    return records

def fasta_write(path, records):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(''.join((f'>{name}\n{seq}\n' for name, seq in records)))

def run_command(argv, log, cwd=None):
    log = Path(log)
    log.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    save_json(str(log) + '.command.json', {'argv': argv, 'cwd': str(cwd) if cwd else None})
    env = dict(os.environ)
    env.setdefault('XLA_PYTHON_CLIENT_PREALLOCATE', 'false')
    with log.open('w') as out:
        p = subprocess.run(argv, cwd=cwd, stdout=out, stderr=subprocess.STDOUT, env=env)
    save_json(str(log) + '.status.json', {'returncode': p.returncode, 'seconds': time.time() - started})
    if p.returncode:
        tail = '\n'.join(log.read_text(errors='replace').splitlines()[-15:])
        raise RuntimeError(f'Command failed ({p.returncode}); log: {log}\n{tail}')

class Stage:

    def __init__(self, directory, signature, inputs=()):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.marker = self.directory / 'complete.json'
        self.signature = json.loads(json.dumps({'settings': signature, 'inputs': {str(Path(p).resolve()): digest(p) for p in inputs}}))

    def reusable(self):
        if not self.marker.exists():
            return False
        old = read_json(self.marker)
        if old['signature'] != self.signature:
            raise RuntimeError(f'Inputs/settings changed in {self.directory}; use a NEW output directory.')
        for name, sha in old['outputs'].items():
            if not Path(name).is_file() or digest(name) != sha:
                raise RuntimeError(f'Completed output missing/changed: {name}; use a NEW output directory.')
        return True

    def finish(self, outputs):
        save_json(self.marker, {'signature': self.signature, 'outputs': {str(Path(p).resolve()): digest(p) for p in outputs}})
