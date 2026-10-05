from pathlib import Path
import argparse
import sys
import numpy as np
import yaml
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from binderflow.structure import Residue, AA1, write_pdb
from binderflow.io import csv_write, save_json
from binderflow.bli import response, fit_file
from binderflow.cli import run
from binderflow.design import load_config

def helix(seq, chain, shift=(0, 0, 0)):

    def point(i):
        angle = np.deg2rad(i * 100)
        return np.array([2.3 * np.cos(angle), 2.3 * np.sin(angle), 1.5 * i]) + np.array(shift)
    rows = []
    for i, aa in enumerate(seq):
        x = point(i)
        prev = point(i - 1)
        nxt = point(i + 1)
        v = (nxt - x) / np.linalg.norm(nxt - x)
        w = (x - prev) / np.linalg.norm(x - prev)
        atoms = {'N': x - 1.2 * w, 'CA': x, 'C': x + 1.3 * v, 'O': x + 1.3 * v + np.array([0, 0, 1.2])}
        rows.append(Residue(chain, i + 1, ' ', AA1[aa], atoms, 90.0))
    return rows

def setup(out):
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    c = yaml.safe_load((root / 'configs/smoke.yaml').read_text())
    pdb = out / 'synthetic.pdb'
    write_pdb(pdb, [helix('ASEKQSTVDEQR', 'D'), helix('QERSTKADVEQS', 'B', (8, 0, 0))])
    c['structure'].update(pdb=str(pdb), target_chain='B', target_ranges=[[1, 12]], donor_chain='D', donor_ranges=[[1, 12]], hotspot_residues=[4, 8])
    c['run'].update(output=str(out / 'run'), synthetic_demo=True)
    c['anchors']['source_residues'] = [4, 9]
    c['design'].update(num_backbones=2, n_flank=[3, 3], c_flank=[3, 3], loop_length_delta=0, length=[12, 12])
    c['mpnn']['sequences_per_backbone'] = 2
    c['prediction'].update(num_models=2, num_seeds=1, num_recycles=1)
    c['filters'].update(min_contact_residue_pairs=1, max_interchain_clash_pairs=100, max_hydrophobic_fraction=1, max_hydrophobic_run=20, max_abs_charge=100)
    c['selection'].update(top_k=6, max_identity=0.95)
    c['provenance_files'] = []
    script = str(root / 'tests/mock_tools.py')
    c['tools'] = {name: {'command': [sys.executable, script, name]} for name in ['rfdiffusion', 'proteinmpnn', 'colabfold']}
    (out / 'demo.yaml').write_text(yaml.safe_dump(c, sort_keys=False))
    rng = np.random.default_rng(7)
    rows = []
    for i, conc in enumerate([5, 20, 80, 320]):
        t = np.arange(0, 601, 3, dtype=float)
        y = response(t, conc * 1e-09, 180.0, 200000.0, 0.004, 1.0 + i * 0.05, 0.01) + rng.normal(0, 0.001, len(t))
        rows.extend(({'trace_id': f'synthetic_{i}', 'time_s': float(a), 'response_nm': float(b), 'concentration_nM': conc, 'association_end_s': 180} for a, b in zip(t, y)))
    csv_write(out / 'SYNTHETIC_bli.csv', rows)
    save_json(out / 'SYNTHETIC_truth.json', {'kon': 200000.0, 'koff': 0.004, 'KD_nM': 20.0, 'note': 'Generated test data; no experiments were performed.'})
    return out / 'demo.yaml'

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', default='demo_output')
    args = ap.parse_args()
    cfg = setup(args.out)
    print('SYNTHETIC DEMO: RFdiffusion, ProteinMPNN and ColabFold are replaced by deterministic test programs.', flush=True)
    run(load_config(cfg))
    result = fit_file(Path(args.out) / 'SYNTHETIC_bli.csv', Path(args.out) / 'SYNTHETIC_bli_fit')
    print(f"Synthetic BLI: recovered KD={result['KD_nM']:.3f} nM (truth: 20 nM)")
if __name__ == '__main__':
    main()
