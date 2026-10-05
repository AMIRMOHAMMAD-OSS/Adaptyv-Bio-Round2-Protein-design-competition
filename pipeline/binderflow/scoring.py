from __future__ import annotations
import html
import re
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from .io import Stage, read_json, save_json, run_command, fasta_write, csv_write
from .structure import read_pdb, select, sequence, ca, fit_transform, rmsd, interface_metrics, context_atoms
from .design import tool_command

def sequence_metrics(seq, ph=7.4):
    hydrophobic = set('AVILMFWY')
    charge = 1 / (1 + 10 ** (ph - 8.0)) - 1 / (1 + 10 ** (3.1 - ph))
    for aa, pka in [('K', 10.5), ('R', 12.5), ('H', 6.0)]:
        charge += seq.count(aa) / (1 + 10 ** (ph - pka))
    for aa, pka in [('D', 3.9), ('E', 4.1), ('C', 8.3), ('Y', 10.1)]:
        charge -= seq.count(aa) / (1 + 10 ** (pka - ph))
    stretches = re.findall('[AVILMFWY]+', seq)
    sequons = len(re.findall('(?=N[^P][ST])', seq))
    return {'length': len(seq), 'charge_pH7_4': float(charge), 'hydrophobic_fraction': sum((x in hydrophobic for x in seq)) / len(seq), 'longest_hydrophobic_run': max(map(len, stretches), default=0), 'n_glycosylation_sequons': sequons, 'cysteines': seq.count('C')}

def prediction_command(c, candidate, fasta, directory):
    cmd, cwd = tool_command(c, 'colabfold')
    p = c['prediction']
    cmd += [str(fasta), str(directory), '--model-type', 'alphafold2_multimer_v3', '--msa-mode', p['msa_mode'], '--pair-mode', 'unpaired', '--num-models', str(p['num_models']), '--num-recycle', str(p['num_recycles']), '--num-seeds', str(p['num_seeds']), '--random-seed', str(candidate['spec']['seed']), '--rank', 'multimer']
    return (cmd, cwd)

def predict(c, meta, candidate):
    out = Path(c['run']['output']) / '04_predictions' / candidate['id']
    if candidate['sequence'] == meta['target_sequence']:
        raise ValueError('Identical binder and target sequence would invoke homomer handling; unsupported')
    st = Stage(out, {'candidate': candidate, 'prediction': c['prediction'], 'tool': c['tools']['colabfold'], 'resolved_tool': tool_command(c, 'colabfold'), 'code': c['_implementation_sha256']}, [candidate['backbone']])
    if not st.reusable():
        model_dir = out / 'models'
        if model_dir.exists() and any(model_dir.iterdir()):
            raise RuntimeError(f'Incomplete prediction directory: {model_dir}. Move it aside before retrying this candidate.')
        fa = out / 'complex.fasta'
        fasta_write(fa, [(candidate['id'], candidate['sequence'] + ':' + meta['target_sequence'])])
        cmd, cwd = prediction_command(c, candidate, fa, model_dir)
        run_command(cmd, out / 'colabfold.log', cwd)
        pairs = model_pairs(model_dir, c['prediction']['num_models'] * c['prediction']['num_seeds'])
        for pdb, js in pairs:
            validate_prediction(pdb, js, candidate['sequence'], meta['target_sequence'])
        st.finish([fa] + [x for pair in pairs for x in pair])
    pairs = model_pairs(out / 'models', c['prediction']['num_models'] * c['prediction']['num_seeds'])
    return pairs

def model_pairs(directory, expected):
    files = sorted(Path(directory).glob('*_scores_rank_*.json'))
    if len(files) != expected:
        raise ValueError(f'Expected {expected} prediction score files in {directory}, got {len(files)}')
    pairs = []
    for js in files:
        pdb = js.with_name(js.name.replace('_scores_rank_', '_unrelaxed_rank_').replace('.json', '.pdb'))
        if not pdb.exists():
            raise FileNotFoundError(f'Missing model paired with scores: {pdb}')
        pairs.append((pdb, js))
    if not any(('_rank_001_' in p.name for p, _ in pairs)):
        raise ValueError('Missing rank 001 model')
    return pairs

def validate_prediction(pdb, scores, binder_seq, target_seq):
    rs = read_pdb(pdb)
    if {r.chain for r in rs} != {'A', 'B'}:
        raise ValueError('ColabFold must return exactly binder A and target B')
    b, t = (select(rs, 'A'), select(rs, 'B'))
    if sequence(b) != binder_seq or sequence(t) != target_seq:
        raise ValueError(f'Prediction sequence or chain order mismatch: {pdb}')
    j = read_json(scores)
    n = len(b) + len(t)
    pl, pae = (np.asarray(j['plddt'], float), np.asarray(j['pae'], float))
    if pl.shape != (n,) or pae.shape != (n, n) or (not np.isfinite(pl).all()) or (not np.isfinite(pae).all()):
        raise ValueError('Prediction confidence arrays have invalid size or values')
    if np.any((pl < 0) | (pl > 100)) or np.any(pae < 0) or (not np.isfinite(j['iptm'])) or (not 0 <= j['iptm'] <= 1):
        raise ValueError('Prediction confidence values out of bounds')
    return (b, t, j, pl, pae)

def score_one(c, meta, candidate, pairs):
    pdb, js = next((pair for pair in pairs if '_rank_001_' in pair[0].name))
    b, t, j, pl, pae = validate_prediction(pdb, js, candidate['sequence'], meta['target_sequence'])
    n = len(b)
    metrics, d = interface_metrics(b, t, c['filters']['contact_cutoff_A'], c['filters']['clash_cutoff_A'])
    contacts = d < c['filters']['contact_cutoff_A']
    bi = (pae[:n, n:] + pae[n:, :n].T) / 2
    contact_pae = float(bi[contacts].mean()) if contacts.any() else None
    br = read_pdb(candidate['backbone'])
    bref, tref = (select(br, 'A'), select(br, 'B'))
    rot, trans = fit_transform(ca(t), ca(tref))
    bp = ca(b) @ rot + trans
    fold_rot, fold_trans = fit_transform(ca(b), ca(bref))
    anchors = [int(i) - 1 for i in candidate['spec']['anchors']]
    motif_rmsd = rmsd(bp[anchors], ca(bref)[anchors]) if anchors else None
    source = read_pdb(c['structure']['pdb'])
    ts = select(source, c['structure']['target_chain'], c['structure']['target_ranges'])
    r2, t2 = fit_transform(ca(t), ca(ts))
    heavy = np.array([v @ r2 + t2 for r in b for v in r.atoms.values()])
    outside, glycans = context_atoms(c['structure']['pdb'], c['structure']['target_chain'], meta['target_original_numbers'])
    cutoff = c['filters']['clash_cutoff_A']
    context_n = int((cdist(heavy, outside) < cutoff).sum()) if len(outside) else 0
    glycan_n = int((cdist(heavy, glycans) < cutoff).sum()) if len(glycans) else 0
    all_iptm = [read_json(s)['iptm'] for _, s in pairs]
    row = {'id': candidate['id'], 'sequence': candidate['sequence'], 'mpnn_score': candidate['mpnn_score'], **sequence_metrics(candidate['sequence']), **metrics, 'binder_plddt': float(pl[:n].mean()), 'iptm': float(j['iptm']), 'iptm_min_across_models': float(min(all_iptm)), 'interchain_pae_A': float(bi.mean()), 'contact_pae_A': contact_pae, 'target_aligned_binder_rmsd_A': rmsd(bp, ca(bref)), 'binder_fold_rmsd_A': rmsd(ca(b) @ fold_rot + fold_trans, ca(bref)), 'target_rmsd_A': rmsd(ca(t) @ rot + trans, ca(tref)), 'motif_rmsd_A': motif_rmsd, 'receptor_context_clash_pairs': context_n, 'observed_glycan_clash_pairs': glycan_n, 'observed_glycan_atoms_checked': len(glycans), 'prediction_pdb': str(pdb)}
    row['failure_reasons'] = filter_reasons(row, c['filters'])
    row['passed'] = not bool(row['failure_reasons'])
    out = Path(c['run']['output']) / '05_scores' / candidate['id']
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / 'residue_min_distances_A.npy', d)
    csv_write(out / 'contacts.csv', [{'binder_position': int(i) + 1, 'target_position': int(k) + 1, 'source_target_number': meta['target_original_numbers'][k], 'min_distance_A': float(d[i, k]), 'bidirectional_pae_A': float(bi[i, k])} for i, k in zip(*np.where(contacts))], ['binder_position', 'target_position', 'source_target_number', 'min_distance_A', 'bidirectional_pae_A'])
    save_json(out / 'metrics.json', row)
    return row

def filter_reasons(row, f):
    failures = []
    checks = [('binder_plddt', '>=', f['min_binder_plddt']), ('iptm', '>=', f['min_iptm']), ('contact_pae_A', '<=', f['max_contact_pae_A']), ('target_aligned_binder_rmsd_A', '<=', f['max_pose_rmsd_A']), ('binder_fold_rmsd_A', '<=', f['max_fold_rmsd_A']), ('contact_residue_pairs', '>=', f['min_contact_residue_pairs']), ('interchain_clash_pairs', '<=', f['max_interchain_clash_pairs']), ('receptor_context_clash_pairs', '<=', f['max_context_clash_pairs']), ('observed_glycan_clash_pairs', '<=', f['max_observed_glycan_clash_pairs']), ('hydrophobic_fraction', '<=', f['max_hydrophobic_fraction']), ('longest_hydrophobic_run', '<=', f['max_hydrophobic_run']), ('n_glycosylation_sequons', '<=', f['max_glycosylation_sequons'])]
    if row['motif_rmsd_A'] is not None:
        checks.append(('motif_rmsd_A', '<=', f['max_motif_rmsd_A']))
    for key, op, limit in checks:
        value = row[key]
        if value is None or not np.isfinite(value) or (op == '>=' and value < limit) or (op == '<=' and value > limit):
            failures.append(f'{key} {op} {limit} (observed {value})')
    if abs(row['charge_pH7_4']) > f['max_abs_charge']:
        failures.append('absolute charge exceeds limit')
    return '; '.join(failures)

def identity(a, b):
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(cur[-1] + 1, prev[j] + 1, prev[j - 1] + (x != y)))
        prev = cur
    return 1 - prev[-1] / max(len(a), len(b))

def select_and_report(c, rows):
    out = Path(c['run']['output']) / '06_results'
    out.mkdir(parents=True, exist_ok=True)
    order = sorted(rows, key=lambda r: (not r['passed'], -r['iptm'], r['contact_pae_A'] if r['contact_pae_A'] is not None else float('inf'), r['target_aligned_binder_rmsd_A'], r['mpnn_score'], r['id']))
    chosen = []
    for row in order:
        row['selected'] = False
        row['selection_note'] = 'failed filters' if not row['passed'] else 'top_k reached'
        if row['passed'] and len(chosen) < c['selection']['top_k']:
            if any((identity(row['sequence'], p['sequence']) >= c['selection']['max_identity'] for p in chosen)):
                row['selection_note'] = 'similar to a higher-ranked selected sequence'
            else:
                row['selected'] = True
                row['selection_note'] = 'selected'
                chosen.append(row)
    csv_write(out / 'all_candidates.csv', order)
    fields = list(order[0]) if order else ['id', 'sequence']
    csv_write(out / 'selected_candidates.csv', chosen, fields)
    fasta_write(out / 'selected.fasta', [(r['id'], r['sequence']) for r in chosen])
    save_json(out / 'summary.json', {'synthetic_demo': c['run'].get('synthetic_demo', False), 'unique_evaluated': len(rows), 'passed': sum((r['passed'] for r in rows)), 'selected': len(chosen), 'requested': c['selection']['top_k']})
    title = 'SYNTHETIC SOFTWARE DEMO — NOT BIOLOGICAL RESULTS' if c['run'].get('synthetic_demo') else 'Computational candidate report'
    tr = ''.join(('<tr>' + ''.join((f'<td>{html.escape(str(r[k]))}</td>' for k in ['id', 'passed', 'selected', 'iptm', 'binder_plddt', 'contact_pae_A', 'target_aligned_binder_rmsd_A', 'failure_reasons'])) + '</tr>' for r in order))
    (out / 'report.html').write_text(f"""<!doctype html><meta charset="utf-8"><title>{title}</title>\n<style>body{{font:15px system-ui;max-width:1400px;margin:40px auto;padding:20px}}table{{border-collapse:collapse;width:100%}}td,th{{padding:10px;text-align:left;border-bottom:1px solid #ddd}}h1{{color:#163b58}}.note{{padding:16px;background:#fff1ce}}</style>\n<h1>{title}</h1><p>{len(rows)} unique candidates evaluated; {sum((r['passed'] for r in rows))} passed; {len(chosen)} selected.</p>\n<p class="note">Scores are screening heuristics, not measured affinities. No passing candidates are invented to fill the quota.\nThe observed-glycan check does not model missing glycans, dynamics, or glycoform diversity.</p>\n<table><thead><tr><th>ID</th><th>Pass</th><th>Selected</th><th>ipTM</th><th>Binder pLDDT</th><th>Contact PAE (Å)</th><th>Pose RMSD (Å)</th><th>Failures</th></tr></thead><tbody>{tr}</tbody></table>""")
    return chosen
