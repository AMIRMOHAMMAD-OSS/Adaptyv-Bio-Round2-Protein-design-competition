from __future__ import annotations
from dataclasses import dataclass, replace
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
AA3 = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(), 'ARNDCQEGHILKMFPSTWYV'))
AA1 = {v: k for k, v in AA3.items()}
BACKBONE = {'N', 'CA', 'C', 'O'}

@dataclass
class Residue:
    chain: str
    number: int
    icode: str
    name: str
    atoms: dict
    bfactor: float = 0.0

    @property
    def aa(self):
        return AA3[self.name]

    @property
    def label(self):
        return f'{self.chain}:{self.number}{self.icode.strip()}'

def read_pdb(path):
    residues = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if not line.startswith('ATOM  ') or line[16] not in ' A':
            continue
        if line[17:20] not in AA3:
            raise ValueError(f'Noncanonical ATOM residue in {path}: {line[17:20]}')
        element = line[76:78].strip() or line[12:16].strip()[0]
        if element in {'H', 'D'}:
            continue
        key = (line[21], int(line[22:26]), line[26])
        if key not in residues:
            residues[key] = Residue(*key, line[17:20], {}, float(line[60:66].strip() or 0))
        atom = line[12:16].strip()
        if atom in residues[key].atoms:
            continue
        xyz = np.array([float(line[i:i + 8]) for i in (30, 38, 46)])
        if not np.isfinite(xyz).all():
            raise ValueError('Nonfinite PDB coordinates')
        residues[key].atoms[atom] = xyz
    if not residues:
        raise ValueError(f'No canonical protein atoms in {path}')
    return list(residues.values())

def select(residues, chain, spans=None):
    out = [r for r in residues if r.chain == chain and (not spans or any((lo <= r.number <= hi for lo, hi in spans)))]
    if not out:
        raise ValueError(f'Empty selection: chain {chain}, ranges {spans}')
    if any((r.icode.strip() for r in out)):
        raise ValueError('Insertion codes are not supported in selected residues; explicitly renumber your PDB first.')
    ids = [r.number for r in out]
    if spans and (not {i for lo, hi in spans for i in range(lo, hi + 1)} <= set(ids)):
        raise ValueError('Requested residue range includes unobserved residues; repair or choose a complete crop')
    if len(ids) != len(set(ids)) or ids != sorted(ids):
        raise ValueError('Selected chain has repeated or unordered residue numbers')
    for r in out:
        if not BACKBONE <= r.atoms.keys():
            raise ValueError(f'Incomplete backbone at {r.label}; repair structure before design')
    for a, b in zip(out, out[1:]):
        if b.number != a.number + 1 or np.linalg.norm(a.atoms['C'] - b.atoms['N']) > 2.2:
            raise ValueError(f'Peptide gap at {a.label} -> {b.label}; choose a continuous crop or repair it')
    return out

def canonical(residues, chain):
    return [replace(r, chain=chain, number=i + 1, icode=' ') for i, r in enumerate(residues)]

def sequence(residues):
    return ''.join((r.aa for r in residues))

def ca(residues):
    return np.array([r.atoms['CA'] for r in residues])

def write_pdb(path, chains):
    lines, serial = ([], 1)
    for chain in chains:
        for r in chain:
            for name, xyz in r.atoms.items():
                if serial > 99999 or r.number > 9999:
                    raise ValueError('PDB field overflow')
                lines.append(f'ATOM  {serial:5d} {name:^4s} {r.name:3s} {r.chain}{r.number:4d}{r.icode:1s}   {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{r.bfactor:6.2f}          {name[0]:>2s}')
                serial += 1
        lines.append('TER')
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text('\n'.join(lines) + '\nEND\n')

def distances(a, b, polar=False):
    out = np.full((len(a), len(b)), np.inf)
    for i, ra in enumerate(a):
        ax = [v for k, v in ra.atoms.items() if not polar or k[0] in 'NO']
        for j, rb in enumerate(b):
            bx = [v for k, v in rb.atoms.items() if not polar or k[0] in 'NO']
            if ax and bx:
                out[i, j] = cdist(ax, bx).min()
    return out

def fit_transform(mobile, reference):
    mobile, reference = (np.asarray(mobile), np.asarray(reference))
    if mobile.shape != reference.shape or len(mobile) < 3:
        raise ValueError('Kabsch requires matching arrays with at least three points')
    if np.linalg.matrix_rank(mobile - mobile.mean(0)) < 2:
        raise ValueError('Cannot determine rotation from collinear target coordinates')
    mc, rc = (mobile.mean(0), reference.mean(0))
    u, _, vt = np.linalg.svd((mobile - mc).T @ (reference - rc))
    d = np.eye(3)
    d[-1, -1] = np.linalg.det(u @ vt)
    rot = u @ d @ vt
    return (rot, rc - mc @ rot)

def rmsd(a, b):
    return float(np.sqrt(np.mean(np.sum((np.asarray(a) - np.asarray(b)) ** 2, axis=1))))

def interface_metrics(binder, target, cutoff=5.0, clash_cutoff=2.0):
    d = distances(binder, target)
    contact = d < cutoff
    heavy_b = np.array([v for r in binder for v in r.atoms.values()])
    heavy_t = np.array([v for r in target for v in r.atoms.values()])
    return ({'contact_residue_pairs': int(contact.sum()), 'binder_interface_residues': int(contact.any(1).sum()), 'target_interface_residues': int(contact.any(0).sum()), 'interchain_clash_pairs': int((cdist(heavy_b, heavy_t) < clash_cutoff).sum())}, d)

def context_atoms(path, target_chain, crop_numbers):
    outside, glycans = ([], [])
    glycan_names = {'NAG', 'BMA', 'MAN', 'FUC', 'GAL', 'GLC', 'SIA', 'NDG', 'BGC'}
    for line in Path(path).read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if not line.startswith(('ATOM  ', 'HETATM')) or line[21] != target_chain or line[16] not in ' A':
            continue
        if (line[76:78].strip() or line[12:16].strip()[0]) in {'H', 'D'}:
            continue
        xyz = [float(line[i:i + 8]) for i in (30, 38, 46)]
        if line.startswith('ATOM') and int(line[22:26]) not in crop_numbers:
            outside.append(xyz)
        elif line.startswith('HETATM') and line[17:20] in glycan_names:
            glycans.append(xyz)
    return (np.array(outside).reshape(-1, 3), np.array(glycans).reshape(-1, 3))
