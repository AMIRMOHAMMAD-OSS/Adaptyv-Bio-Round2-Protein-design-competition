from __future__ import annotations
from dataclasses import replace
import os
from pathlib import Path
import random
import re
import numpy as np
import yaml
from .io import Stage, digest, save_json, read_json, csv_write, fasta_read, fasta_write, run_command
from .structure import AA1, BACKBONE, read_pdb, select, canonical, sequence, ca, write_pdb, distances, fit_transform, rmsd

def load_config(path):
    path = Path(path).resolve()
    c = yaml.safe_load(path.read_text())
    if c.get('schema_version') != 1:
        raise ValueError('Expected schema_version: 1')
    if c['design']['mode'] not in {'motif', 'denovo'}:
        raise ValueError('design.mode must be motif or denovo')
    c['structure']['pdb'] = str((path.parent / os.path.expandvars(c['structure']['pdb'])).resolve())
    c['run']['output'] = str((path.parent / os.path.expandvars(c['run']['output'])).resolve())
    if not Path(c['structure']['pdb']).is_file():
        raise FileNotFoundError(f"Missing structure {c['structure']['pdb']}; run scripts/fetch_example.py first")
    d = c['design']
    if d['num_backbones'] < 1 or c['mpnn']['sequences_per_backbone'] < 1:
        raise ValueError('Design counts must be positive')
    if not isinstance(d['seed'], int) or d['seed'] < 1:
        raise ValueError('Use a positive integer seed (ProteinMPNN treats zero as random)')
    for bounds in (d['n_flank'], d['c_flank'], d['length']):
        if len(bounds) != 2 or bounds[0] < 0 or bounds[1] < bounds[0]:
            raise ValueError('Invalid design length interval')
    if not 0 < c['mpnn']['temperature'] <= 1:
        raise ValueError('Expected MPNN temperature in (0,1]')
    if c['prediction']['msa_mode'] not in {'single_sequence', 'mmseqs2_uniref_env', 'mmseqs2_uniref'}:
        raise ValueError('Unsupported msa_mode')
    if not 1 <= c['prediction']['num_models'] <= 5 or c['prediction']['num_seeds'] < 1:
        raise ValueError('Invalid prediction model/seed count')
    if c['selection']['top_k'] < 1 or not 0 < c['selection']['max_identity'] <= 1:
        raise ValueError('Invalid selection settings')
    c['_implementation_sha256'] = digest_sources()
    return c

def digest_sources():
    import hashlib
    h = hashlib.sha256()
    for p in sorted(Path(__file__).parent.glob('*.py')):
        h.update(p.name.encode())
        h.update(p.read_bytes())
    return h.hexdigest()

def prepare(c):
    out = Path(c['run']['output']) / '01_prepared'
    st = Stage(out, {'structure': c['structure'], 'anchors': c['anchors'], 'design': c['design'], 'code': c['_implementation_sha256']}, [c['structure']['pdb']])
    if st.reusable():
        return read_json(out / 'metadata.json')
    residues = read_pdb(c['structure']['pdb'])
    s = c['structure']
    target = select(residues, s['target_chain'], s['target_ranges'])
    donor = select(residues, s['donor_chain'], s['donor_ranges']) if s.get('donor_chain') else []
    if donor and s['target_chain'] == s['donor_chain']:
        raise ValueError('Target and motif donor must be different chains')
    ac = c['anchors']
    anchor_rows, fixed = ([], [])
    if donor:
        d, pd = (distances(donor, target), distances(donor, target, polar=True))
        for i, r in enumerate(donor):
            anchor_rows.append({'source_residue': r.label, 'canonical_position': i + 1, 'amino_acid': r.aa, 'nearest_target': target[int(d[i].argmin())].label, 'distance_A': float(d[i].min()), 'polar_distance_A': float(pd[i].min())})
        explicit = ac.get('source_residues') or []
        if explicit:
            by_num = {r.number: i + 1 for i, r in enumerate(donor)}
            if not set(explicit) <= by_num.keys():
                raise ValueError('An explicitly requested anchor is absent from the donor crop')
            fixed = sorted({by_num[n] for n in explicit})
        else:
            pool = [r for r in anchor_rows if r['distance_A'] <= ac['contact_cutoff_A']]
            if ac['require_polar_contact']:
                pool = [r for r in pool if r['polar_distance_A'] <= ac['polar_cutoff_A']]
            pool.sort(key=lambda r: (r['polar_distance_A'] > ac['polar_cutoff_A'], r['distance_A'], r['canonical_position']))
            fixed = sorted((r['canonical_position'] for r in pool[:ac['max_anchors']]))
    if c['design']['mode'] == 'motif' and len(fixed) < 1:
        raise ValueError('Motif mode has no anchors; check chain pairing, crop, and cutoffs')
    for r in anchor_rows:
        r['selected'] = r['canonical_position'] in fixed
    hot_source = s.get('hotspot_residues') or []
    target_index = {r.number: i + 1 for i, r in enumerate(target)}
    if hot_source:
        if not set(hot_source) <= target_index.keys():
            raise ValueError('Hotspot absent from target crop')
        hotspots = [target_index[x] for x in hot_source]
    elif donor:
        td = distances(target, donor).min(1)
        hotspots = [int(i) + 1 for i in np.argsort(td)[:3] if td[i] <= ac['contact_cutoff_A']]
    else:
        hotspots = []
    if not hotspots:
        raise ValueError('No target hotspots; supply structure.hotspot_residues in original PDB numbering')
    map_rows = [{'role': role, 'source': r.label, 'canonical': f'{chain}:{i + 1}', 'amino_acid': r.aa} for role, chain, rs in [('donor', 'A', donor), ('target', 'B', target)] for i, r in enumerate(rs)]
    write_pdb(out / 'input.pdb', [canonical(donor, 'A'), canonical(target, 'B')] if donor else [canonical(target, 'B')])
    fasta_write(out / 'target.fasta', [('target', sequence(target))])
    csv_write(out / 'residue_mapping.csv', map_rows)
    csv_write(out / 'anchor_candidates.csv', anchor_rows, ['source_residue', 'canonical_position', 'amino_acid', 'nearest_target', 'distance_A', 'polar_distance_A', 'selected'])
    meta = {'target_sequence': sequence(target), 'donor_sequence': sequence(donor), 'fixed_donor_positions': fixed, 'hotspots': hotspots, 'target_original_numbers': [r.number for r in target], 'source_sha256': digest(s['pdb']), 'note': 'N/O distance is a polar-contact proxy, not a hydrogen-bond assignment.'}
    save_json(out / 'metadata.json', meta)
    st.finish([out / x for x in ['input.pdb', 'target.fasta', 'residue_mapping.csv', 'anchor_candidates.csv', 'metadata.json']])
    return meta

def design_spec(c, meta, index):
    cfg = c['design']
    seed = cfg['seed'] + index
    rng = random.Random(seed)
    tokens, anchors, position = ([], {}, 0)
    if cfg['mode'] == 'denovo':
        length = rng.randint(*cfg['length'])
        tokens = [str(length)]
    else:

        def gap(n):
            nonlocal position
            if n:
                tokens.append(str(n))
                position += n
        gap(rng.randint(*cfg['n_flank']))
        previous = None
        for p in meta['fixed_donor_positions']:
            if previous is not None:
                n = p - previous - 1
                if n:
                    delta = cfg['loop_length_delta']
                    gap(rng.randint(max(1, n - delta), n + delta))
            tokens.append(f'A{p}-{p}')
            position += 1
            anchors[str(position)] = {'donor_position': p, 'aa': meta['donor_sequence'][p - 1]}
            previous = p
        gap(rng.randint(*cfg['c_flank']))
        length = position
    if not cfg['length'][0] <= length <= cfg['length'][1]:
        raise ValueError(f'Sampled binder length {length} outside design.length; revise flanks/loops or allowed interval')
    return {'id': f'bb{index:04d}', 'seed': seed, 'length': length, 'anchors': anchors, 'contig': '/'.join(tokens) + f"/0 B1-{len(meta['target_sequence'])}"}

def tool_command(c, name, strict=True):
    t = c['tools'][name]
    argv = [os.path.expanduser(os.path.expandvars(str(v))) for v in t['command']]
    cwd = os.path.expanduser(os.path.expandvars(t.get('cwd', ''))) or None
    if strict and any(('${' in v for v in argv + ([cwd] if cwd else []))):
        raise ValueError(f'Unset environment variable in {name} command; see README.md')
    return (argv, cwd)

def rf_command(c, meta, spec, strict=True):
    out = Path(c['run']['output']) / '02_backbones' / spec['id']
    cmd, cwd = tool_command(c, 'rfdiffusion', strict)
    inp = Path(c['run']['output']) / '01_prepared' / 'input.pdb'
    cmd += [f'inference.input_pdb={inp}', f"inference.output_prefix={out / 'rf'}", 'inference.num_designs=1', f"inference.design_startnum={spec['seed']}", 'inference.deterministic=True', 'inference.cautious=False', f"diffuser.T={c['design']['diffusion_steps']}", f"contigmap.contigs=[{spec['contig']}]", 'ppi.hotspot_res=[' + ','.join((f'B{x}' for x in meta['hotspots'])) + ']']
    return (cmd, cwd)

def generate(c, meta, spec):
    out = Path(c['run']['output']) / '02_backbones' / spec['id']
    inp = Path(c['run']['output']) / '01_prepared' / 'input.pdb'
    cmd, cwd = rf_command(c, meta, spec)
    st = Stage(out, {'spec': spec, 'argv': cmd, 'cwd': cwd, 'code': c['_implementation_sha256']}, [inp])
    if st.reusable():
        return out / 'complex.pdb'
    run_command(cmd, out / 'rfdiffusion.log', cwd)
    raw = out / f"rf_{spec['seed']}.pdb"
    residues = read_pdb(raw)
    if {r.chain for r in residues} != {'A', 'B'}:
        raise ValueError('RFdiffusion output must contain binder A and target B')
    binder = [r for r in residues if r.chain == 'A']
    target = [r for r in residues if r.chain == 'B']
    if len(binder) != spec['length'] or len(target) != len(meta['target_sequence']):
        raise ValueError('RFdiffusion chain length disagrees with explicit contig')
    source = read_pdb(inp)
    source_target = [r for r in source if r.chain == 'B']
    source_donor = [r for r in source if r.chain == 'A']
    rot, trans = fit_transform(ca(target), ca(source_target))
    trmsd = rmsd(ca(target) @ rot + trans, ca(source_target))
    if trmsd > c['design']['fixed_geometry_tolerance_A']:
        raise ValueError(f'RFdiffusion changed target geometry: {trmsd:.3f} A')
    for pos, a in spec['anchors'].items():
        p = int(pos) - 1
        ar = source_donor[a['donor_position'] - 1]
        err = rmsd(np.array([binder[p].atoms[k] @ rot + trans for k in sorted(BACKBONE)]), np.array([ar.atoms[k] for k in sorted(BACKBONE)]))
        if err > c['design']['fixed_geometry_tolerance_A']:
            raise ValueError(f'Fixed motif geometry moved at binder {pos}: {err:.3f} A')
        binder[p] = replace(binder[p], name=AA1[a['aa']])
    binder = [replace(r, atoms={k: v for k, v in r.atoms.items() if k in BACKBONE}) for r in binder]
    target = [replace(r, name=AA1[aa], atoms={k: v for k, v in r.atoms.items() if k in BACKBONE}) for r, aa in zip(target, meta['target_sequence'])]
    write_pdb(out / 'complex.pdb', [canonical(binder, 'A'), canonical(target, 'B')])
    check = read_pdb(out / 'complex.pdb')
    select(check, 'A')
    select(check, 'B')
    save_json(out / 'spec.json', spec)
    st.finish([raw, out / 'complex.pdb', out / 'spec.json'])
    return out / 'complex.pdb'

def mpnn_input(pdb, name):
    rs = read_pdb(pdb)
    obj = {'name': name, 'num_of_chains': 2}
    for chain in ('A', 'B'):
        cr = select(rs, chain)
        obj[f'seq_chain_{chain}'] = sequence(cr)
        obj[f'coords_chain_{chain}'] = {f'{atom}_chain_{chain}': [r.atoms[atom].tolist() for r in cr] for atom in ('N', 'CA', 'C', 'O')}
    obj['seq'] = obj['seq_chain_A'] + obj['seq_chain_B']
    return obj

def sequences(c, meta, spec, pdb):
    out = Path(c['run']['output']) / '03_sequences' / spec['id']
    cmd, cwd = tool_command(c, 'proteinmpnn')
    mp = c['mpnn']
    st = Stage(out, {'spec': spec, 'mpnn': mp, 'command': cmd, 'cwd': cwd, 'code': c['_implementation_sha256']}, [pdb])
    if st.reusable():
        return read_json(out / 'candidates.json')
    name = spec['id']
    (out / 'parsed.jsonl').write_text(__import__('json').dumps(mpnn_input(pdb, name)) + '\n')
    save_json(out / 'chains.json', {name: [['A'], ['B']]})
    save_json(out / 'fixed.json', {name: {'A': sorted(map(int, spec['anchors'])), 'B': []}})
    cmd += ['--jsonl_path', str(out / 'parsed.jsonl'), '--chain_id_jsonl', str(out / 'chains.json'), '--fixed_positions_jsonl', str(out / 'fixed.json'), '--out_folder', str(out), '--num_seq_per_target', str(mp['sequences_per_backbone']), '--batch_size', '1', '--sampling_temp', str(mp['temperature']), '--seed', str(spec['seed']), '--model_name', mp['model_name'], '--omit_AAs', mp['omit_amino_acids'], '--save_score', '1']
    run_command(cmd, out / 'proteinmpnn.log', cwd)
    fa = out / 'seqs' / f'{name}.fa'
    samples = [(h, s) for h, s in fasta_read(fa) if re.search('(?:^|,\\s*)sample=', h)]
    if len(samples) != mp['sequences_per_backbone']:
        raise ValueError(f"Expected {mp['sequences_per_backbone']} MPNN samples, got {len(samples)}")
    candidates = []
    for i, (header, seq) in enumerate(samples):
        if len(seq) != spec['length'] or not set(seq) <= AA1.keys():
            raise ValueError('Invalid MPNN sequence; expected ONLY designed chain A in FASTA')
        for pos, a in spec['anchors'].items():
            if seq[int(pos) - 1] != a['aa']:
                raise ValueError(f'MPNN mutated a fixed anchor: {name}, position {pos}')
        score = re.search('(?:^|,\\s*)score=([0-9.eE+-]+)', header)
        if not score or not np.isfinite(float(score[1])):
            raise ValueError('MPNN score missing or nonfinite')
        candidates.append({'id': f'{name}_s{i:03d}', 'sequence': seq, 'backbone': str(pdb), 'mpnn_score': float(score[1]), 'spec': spec})
    save_json(out / 'candidates.json', candidates)
    st.finish([fa, out / 'parsed.jsonl', out / 'chains.json', out / 'fixed.json', out / 'candidates.json'])
    return candidates
