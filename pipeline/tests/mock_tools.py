from pathlib import Path
from dataclasses import replace
import json
import random
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from binderflow.structure import read_pdb, write_pdb, AA1
from binderflow.io import fasta_read, fasta_write, save_json
from scripts.demo import helix
name = sys.argv[1]
args = sys.argv[2:]
if name == 'rfdiffusion':
    d = dict((a.split('=', 1) for a in args))
    assert d['inference.deterministic'] == 'True'
    src = read_pdb(d['inference.input_pdb'])
    assert d['contigmap.contigs'] == '[3/A4-4/4/A9-9/3/0 B1-12]'
    p = Path(d['inference.output_prefix'] + '_' + d['inference.design_startnum'] + '.pdb')
    write_pdb(p, [[r for r in src if r.chain == 'A'], [r for r in src if r.chain == 'B']])
elif name == 'proteinmpnn':
    d = dict(zip(args[::2], args[1::2]))
    obj = json.loads(Path(d['--jsonl_path']).read_text())
    chains = json.loads(Path(d['--chain_id_jsonl']).read_text())
    assert chains[obj['name']] == [['A'], ['B']]
    fixed = json.loads(Path(d['--fixed_positions_jsonl']).read_text())[obj['name']]['A']
    native = obj['seq_chain_A']
    rng = random.Random(int(d['--seed']))
    records = [('native', native)]
    for i in range(int(d['--num_seq_per_target'])):
        seq = [rng.choice('ADEKQRSTV') for _ in native]
        for pos in fixed:
            seq[pos - 1] = native[pos - 1]
        records.append((f'T=0.1, sample={i + 1}, score={0.9 + 0.1 * i}, global_score=1.0', ''.join(seq)))
    fasta_write(Path(d['--out_folder']) / 'seqs' / (obj['name'] + '.fa'), records)
elif name == 'colabfold':
    fa, out = (Path(args[0]), Path(args[1]))
    out.mkdir(parents=True, exist_ok=True)
    d = dict(zip(args[2::2], args[3::2]))
    assert d['--model-type'] == 'alphafold2_multimer_v3'
    label, seq = fasta_read(fa)[0]
    b, t = seq.split(':')
    good = not label.endswith('s001')
    for i in range(int(d['--num-models']) * int(d['--num-seeds'])):
        tag = f'rank_{i + 1:03d}_alphafold2_multimer_v3_model_{i + 1}_seed_000'
        write_pdb(out / f'{label}_unrelaxed_{tag}.pdb', [helix(b, 'A'), helix(t, 'B', (8, 0, 0))])
        n = len(b) + len(t)
        pae = np.full((n, n), 4.0 if good else 18.0)
        np.fill_diagonal(pae, 0)
        save_json(out / f'{label}_scores_{tag}.json', {'plddt': [90.0] * n, 'pae': pae.tolist(), 'iptm': 0.85 - 0.03 * i if good else 0.3, 'ptm': 0.85})
else:
    raise ValueError(name)
