# Scientific choices and limits

## Target identity and historical reproduction

EGF is a ligand; EGFR is its receptor. The competition organizer says participants
designed EGFR binders. The supplied entry instead labels EGF in several places,
while describing a template-guided EGFR–TGF-α procedure. This package adopts
EGFR explicitly. To bind EGF, supply an actual EGF-containing complex, select EGF
as target, and choose the appropriate donor interface. Renaming EGFR to EGF
without changing the structure is not valid.

Exact historical reproduction additionally requires the original structure and
assembly, chains, numbering map, anchor list and hydrogen-bond definition,
seeds, model/checkpoint revisions, backbone/sequence counts, folding settings,
filter thresholds, submitted sequences, and raw BLI traces with assay metadata.
Those are absent from the provided README. The new default parameters are
research starting points, not recovered competition settings.

## Preparation and anchor constraints

The checked example retains an observed, continuous EGFR fragment B310–501.
It is a structural example crop, not a claim about the precise experimental
construct boundaries. The ligand donor is D3–50 in 1MOX. Waters and heteroatoms
are excluded from model input. Protein atoms and observed glycans on the same
original receptor chain outside the crop are checked later for clashes.
Incomplete backbones, noncanonical ATOM residues, selected insertion codes,
and peptide gaps are rejected instead of silently concatenated. Missing loops
must be explicitly repaired before use.

Automatic donor anchor selection uses minimum heavy-atom distance to the target,
prioritizing residues with N/O close contacts. `require_polar_contact` can require
a 3.5 Å N/O proximity. **This is not a chemically assigned hydrogen bond:** it
lacks reliable protonation, donor/acceptor typing, and donor–H–acceptor angles.
For a vetted hydrogen-bond-based motif, supply explicit donor residue numbers
in `anchors.source_residues` after structural review.

Only donor anchors enter the RFdiffusion contig as fixed motif residues. Other
positions are generated. Native gaps between anchors provide initial linker
lengths, with configurable integer variation; N/C flanks are sampled separately.
This does not copy and freeze an entire receptor fragment by accident. All
lengths are sampled before launching RFdiffusion so the mapping to binder
positions is explicit and does not rely on guesswork or untrusted pickle files.
Sequence and geometry are different constraints: RFdiffusion preserves motif
backbone coordinates; ProteinMPNN fixes the corresponding amino acids. Native
anchor side-chain rotamers are not guaranteed and are not advertised as fixed.
The target is fixed and excluded from MPNN sequence design. Fixed anchor
identities are verified after MPNN, and predicted full-chain sequences are
verified after ColabFold.

## Model stages

RFdiffusion uses deterministic mode and the design index as its documented seed.
Changing hardware/software can still change floating-point results. MPNN uses
full-backbone `v_48_020` with a configured sampling temperature and positive seed.
ColabFold predicts a heterodimer with binder first and target second, using
AlphaFold2-Multimer-v3 without a supplied design template. Single-sequence mode
is the default to avoid MSA-server dependence; its predictive quality may differ
from MSA-based validation. If enabling the remote MSA mode, sequences are sent
to the configured upstream MSA service; private/local MSA support should be
configured in your ColabFold installation.

ProteinMPNN produces sequences, not all-atom repacked structures. Its output is
passed to structure prediction before interface scoring. The score stage uses
the best-ranked Multimer prediction (0.8 ipTM + 0.2 pTM in upstream Multimer
ranking) and reports the minimum ipTM across all requested models. This is not
a requirement that every model pass. The pipeline does not calculate Rosetta
dG, experimental affinity, ESM scores, exact hydrogen-bond energies, or a
secondary-structure-guided scaffold search. Those were not specified sufficiently
in the source README; they are not replaced with fabricated values.

## Metrics and selection

- Contacts: residue pairs with any heavy atoms closer than 5 Å.
- Clashes: cross-chain heavy-atom pairs closer than 2 Å (a simple overlap screen,
  not a van der Waals energy). Unrelaxed predictions can fail this strict screen.
- Contact PAE: mean of both directional PAE blocks over contacting residue pairs.
  Whole cross-chain PAE is reported separately. Zero contacts fails the PAE screen.
- Pose RMSD: binder Cα RMSD after fitting the target to the designed target.
- Fold RMSD: binder Cα RMSD after fitting the binder to its own designed backbone.
- Motif RMSD: fixed-anchor Cα RMSD in the target-aligned frame.
- Context clashes: fit the predicted target to the original source crop and check
  the binder against the rest of that receptor and its observed glycan atoms.
  This cannot detect missing glycan atoms, unobserved loops, conformational
  changes, other oligomer partners, or full glycoform diversity. A final
  full-construct/assembly review remains necessary.
- Sequence screens: hydrophobic fraction and runs, N-X-S/T sequons (X ≠ P),
  cysteine count, and a Henderson–Hasselbalch charge estimate at pH 7.4. These
  are screening descriptors, not validated solubility or expression predictions.

Thresholds are visible in the YAML and intentionally not represented as universal.
Passing candidates sort by decreasing ipTM, then lower contact PAE, pose RMSD,
and MPNN negative-log-likelihood. Diversity uses
`1 - Levenshtein_distance / max(sequence_lengths)`, with an explicit threshold.
This length-aware heuristic is not structural clustering. Duplicate sequences
are evaluated only once, using the best MPNN-scoring backbone representative;
other backbone associations remain in `duplicates.json`. This can discard
alternative poses and is a declared efficiency tradeoff.

## BLI model

For association at concentration C and zero initial bound occupancy:

`R(t) = offset + Rmax * kon*C/(kon*C + koff) * (1 - exp(-(kon*C + koff)*t))`

For dissociation after time ta:

`R(t) = offset + (R(ta) - offset) * exp(-koff*(t - ta))`

The implementation optimizes positive parameters in log space with several
starting points, using ordinary least squares. Concentrations enter in molar
units after conversion from the CSV's explicit nM column. Time is in seconds;
kon is M⁻¹s⁻¹, koff s⁻¹, and KD molar. The response is continuous at the switch.

This 1:1 model assumes no major mass transport limitation, rebinding, avidity,
heterogeneous sites, drift, carryover occupancy, or concentration uncertainty.
Inspect residuals and assay controls; the code cannot establish that these
assumptions hold. Inputs must already be reference-subtracted. Offsets do not
replace drift correction or reference subtraction. Per-trace Rmax is useful
for unequal sensor loading but can weaken identifiability. The local-Jacobian
confidence interval retains the covariance between log-kon and log-koff;
it does not cover serially correlated noise or systematic errors. The numerical
bounds are guardrails, not physical assertions. Boundary or poorly identified
solutions are flagged for review, and no binding verdict is generated.

## Sources checked for this implementation

- Original project: https://github.com/AMIRMOHAMMAD-OSS/Adaptyv-Bio-Round2-Protein-design-competition
- Competition organizer: https://www.adaptyvbio.com/blog/benchbb/
- Competition paper: https://doi.org/10.1101/2025.04.17.648362
- Structure record: https://www.rcsb.org/structure/1MOX
- RFdiffusion documentation: https://github.com/RosettaCommons/RFdiffusion
- ProteinMPNN code/interface: https://github.com/dauparas/ProteinMPNN
- ColabFold code/interface: https://github.com/sokrypton/ColabFold
- LocalColabFold installation: https://github.com/YoshitakaMo/localcolabfold

Source SHAs used by the installer are in `upstream_revisions.json`. Python
package transitive dependencies and system libraries are not fully locked until
you record the actual installation; therefore this is a source-pinned reference
setup, not a tested portable container image.
