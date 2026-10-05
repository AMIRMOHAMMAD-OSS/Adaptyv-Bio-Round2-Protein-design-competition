from __future__ import annotations
import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import shlex
import sys
from .design import load_config, prepare, design_spec, rf_command, generate, sequences
from .io import save_json, read_json, digest
from .scoring import predict, score_one, select_and_report

@contextmanager
def exclusive(out):
    path = Path(out)
    path.mkdir(parents=True, exist_ok=True)
    lock = path / 'RUNNING.lock'
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 384)
    except FileExistsError:
        raise RuntimeError(f'Another process may be using {path}. If its job died, verify it is stopped, then remove {lock}.')
    try:
        os.write(fd, f'pid={os.getpid()}\n'.encode())
        os.close(fd)
        yield
    finally:
        lock.unlink(missing_ok=True)

def provenance(c):
    paths = c.get('provenance_files', [])
    records = {}
    for value in paths:
        path = Path(os.path.expanduser(os.path.expandvars(value))).resolve()
        if not path.is_file():
            raise FileNotFoundError(f'Missing tool installation record: {path}; run scripts/record_tools.py')
        records[str(path)] = digest(path)
    return records

def run(c, through='all'):
    out = Path(c['run']['output'])
    with exclusive(out):
        sig = {'config': c, 'structure_sha256': digest(c['structure']['pdb']), 'tool_records': provenance(c)}
        manifest = out / 'run_manifest.json'
        if manifest.exists() and read_json(manifest) != sig:
            raise RuntimeError('Run settings, source, tool record, or pipeline code changed. Choose a NEW run.output directory.')
        save_json(manifest, sig)
        meta = prepare(c)
        if through == 'prepare':
            return
        specs = [design_spec(c, meta, i) for i in range(c['design']['num_backbones'])]
        save_json(out / 'design_specs.json', specs)
        candidates = []
        for i, spec in enumerate(specs, 1):
            print(f"[{i}/{len(specs)}] Backbone {spec['id']}, length={spec['length']}, seed={spec['seed']}", flush=True)
            pdb = generate(c, meta, spec)
            if through != 'generate':
                candidates.extend(sequences(c, meta, spec, pdb))
        if through == 'generate':
            return
        unique = {}
        dups = []
        for candidate in sorted(candidates, key=lambda r: (r['mpnn_score'], r['id'])):
            seq = candidate['sequence']
            if seq in unique:
                dups.append({'duplicate': candidate['id'], 'representative': unique[seq]['id']})
            else:
                unique[seq] = candidate
        save_json(out / 'candidates.json', list(unique.values()))
        save_json(out / 'duplicates.json', dups)
        if through == 'sequences':
            return
        rows = []
        for i, candidate in enumerate(unique.values(), 1):
            print(f"[{i}/{len(unique)}] Predict and score {candidate['id']}", flush=True)
            pairs = predict(c, meta, candidate)
            rows.append(score_one(c, meta, candidate, pairs))
        chosen = select_and_report(c, rows)
        print(f"Finished: {len(rows)} unique candidates; {len(chosen)} selected. Results: {out / '06_results'}", flush=True)

def main(argv=None):
    p = argparse.ArgumentParser(description='Reproducible binder design and BLI fitting')
    sub = p.add_subparsers(dest='cmd', required=True)
    for name in ['prepare', 'plan', 'run']:
        sp = sub.add_parser(name)
        sp.add_argument('config')
        if name == 'run':
            sp.add_argument('--through', choices=['prepare', 'generate', 'sequences', 'all'], default='all')
    sp = sub.add_parser('inspect')
    sp.add_argument('pdb')
    sp = sub.add_parser('bli')
    sp.add_argument('csv')
    sp.add_argument('--out', required=True)
    sp.add_argument('--shared-rmax', action='store_true')
    args = p.parse_args(argv)
    try:
        if args.cmd == 'inspect':
            from .structure import read_pdb, sequence
            rs = read_pdb(args.pdb)
            for chain in dict.fromkeys((r.chain for r in rs)):
                cr = [r for r in rs if r.chain == chain]
                print(f'Chain {chain}: {len(cr)} observed residues, {cr[0].label}–{cr[-1].label}\n{sequence(cr)}')
        elif args.cmd == 'bli':
            from .bli import fit_file
            print(json.dumps(fit_file(args.csv, args.out, args.shared_rmax), indent=2))
        elif args.cmd == 'prepare':
            c = load_config(args.config)
            with exclusive(c['run']['output']):
                print(json.dumps(prepare(c), indent=2))
        elif args.cmd == 'plan':
            c = load_config(args.config)
            with exclusive(c['run']['output']):
                meta = prepare(c)
            print('RFdiffusion commands (MPNN and prediction commands depend on their output):')
            for i in range(c['design']['num_backbones']):
                spec = design_spec(c, meta, i)
                cmd, cwd = rf_command(c, meta, spec, strict=False)
                print(f"{spec['id']} • length {spec['length']} • cwd {cwd}\n{shlex.join(cmd)}")
            print(f"Maximum sequences: {c['design']['num_backbones'] * c['mpnn']['sequences_per_backbone']}")
        else:
            run(load_config(args.config), args.through)
    except (ValueError, RuntimeError, FileNotFoundError, KeyError) as e:
        print(f'ERROR: {e}', file=sys.stderr)
        return 2
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
