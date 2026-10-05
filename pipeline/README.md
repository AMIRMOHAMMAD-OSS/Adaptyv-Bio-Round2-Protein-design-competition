# Ubuntu binder generation

This implementation runs RFdiffusion, ProteinMPNN, and ColabFold with real model weights. It prepares an EGFR interface, generates motif-constrained backbones, designs sequences with fixed anchors, predicts complexes, scores them, and selects passing candidates. Optional BLI fitting accepts experimental trace CSV files.

The supplied example targets **EGFR**, using TGF-alpha as the motif donor in PDB 1MOX. EGF and EGFR are different proteins. The original project report uses both names; the target in this implementation is explicit in the configuration.

This is a new implementation of the described approach. The repository does not contain the original structure choice, design-numbering map, anchor list, seeds, or complete run settings needed to reproduce the competition sequences exactly. The example settings do not recover those missing inputs.

## Requirements

Use Ubuntu x86_64 or Ubuntu under WSL2, Conda, Git, curl, and a working NVIDIA driver. The setup uses separate Conda environments for RFdiffusion, ProteinMPNN, and ColabFold. The environment names retain the `-wsl` suffix to match the existing installation; they also work on native Ubuntu.

A fresh installation downloads several GB of Python and GPU dependencies, about 0.97 GB of RFdiffusion weights, and a 4.10 GB AlphaFold parameter archive. Budget at least 30 GiB of actual free disk space for installation and extraction, plus room for outputs. This is a planning allowance, not a fixed installation size. On WSL, check the Windows drive that stores Ubuntu, as well as `df -h /`. The downloader checks `/mnt/c` by default on WSL; set `BINDER_HOST_DRIVE` to another mounted host drive if the distribution is stored there.

The 6 GB RTX 3060 laptop GPU passed the reported import and arithmetic checks. Full protein inference on that GPU remains unverified. Start with the one-candidate smoke configuration and close other GPU applications. Memory use depends on sequence length and model settings.

## New installation

Run these commands from an Ubuntu terminal with Conda available:

```bash
git clone https://github.com/AMIRMOHAMMAD-OSS/Adaptyv-Bio-Round2-Protein-design-competition.git
cd Adaptyv-Bio-Round2-Protein-design-competition/pipeline
bash scripts/install.sh
bash scripts/download_models.sh
bash scripts/run.sh smoke
```

The scripts resolve their project directory automatically and write local tool paths to `tools.env`. Use `bash` to execute the scripts. If a download is interrupted, rerun `bash scripts/download_models.sh`; completed downloads are reused and partial downloads resume.

After a successful smoke run:

```bash
bash scripts/run.sh full
```

The full configuration requests 5 backbones and 4 sequences per backbone. Duplicate sequences are evaluated once. It selects up to 6 candidates that pass the filters and diversity criterion. A result with zero selected candidates is valid.

## Use the existing installation

Uploading these files does not require running the installer or downloading weights again. Keep an ongoing download in its current project directory. After it finishes, the existing environments and model checkouts can be linked to this copy:

```bash
cd Adaptyv-Bio-Round2-Protein-design-competition/pipeline
bash scripts/configure.sh "$HOME/projects/egfr_binder_pipeline/external_wsl"
bash scripts/check_gpu.sh
bash scripts/download_models.sh
bash scripts/run.sh smoke
```

The setup records paths and checkpoint hashes in `.local/tool_versions.json`. Reusing model files does not reuse old run directories. Each new configuration or code version needs a new `run.output` path. If a tool record changes, the setup preserves both versions for review instead of overwriting the record attached to earlier runs.

## Configuration and outputs

Edit `configs/smoke.yaml` or `configs/full.yaml` before beginning a new run. Both use PDB 1MOX, EGFR chain B residues 310–501, and TGF-alpha chain D residues 3–50. Source residue numbers and generated binder positions are recorded separately. The polar-contact screen uses an N/O distance proxy; it does not assign hydrogen bonds.

| Stage | Output directory |
| --- | --- |
| Prepared structure, residue mapping, anchors | `01_prepared/` |
| Generated backbones and RFdiffusion logs | `02_backbones/` |
| ProteinMPNN sequences and fixed positions | `03_sequences/` |
| Predicted complexes and ColabFold logs | `04_predictions/` |
| Interface measurements | `05_scores/` |
| Ranked candidates and report | `06_results/` |

These directories are under `runs/egfr_smoke_v1/` or `runs/egfr_full_v1/`. The results include `all_candidates.csv`, `selected_candidates.csv`, `selected.fasta`, `summary.json`, and `report.html`. Read per-stage logs while a model is running; detailed model output is not streamed to the top-level terminal.

Rerunning the same configuration verifies and reuses completed stages. If a ColabFold job failed, the pipeline reports the partial `models` directory that must be moved aside before retrying. Do this only after the old process has stopped. Changed settings, code, inputs, or recorded tools require a fresh output directory.

Single-sequence prediction is the default and does not contact an MSA server. Enabling either `mmseqs2` mode sends sequences to the upstream MSA service. Prediction scores and heuristic filters do not establish experimental binding or affinity. See [METHODS.md](docs/METHODS.md) for metric definitions and limitations.

## Tests and optional BLI fitting

```bash
conda activate binderflow
cd Adaptyv-Bio-Round2-Protein-design-competition/pipeline
python -m unittest discover -s tests -v
```

The tests use deterministic synthetic model programs. They check orchestration, fixed residues, chain assignments, output validation, resume behavior, BLI parameter recovery, and interrupted downloads. They do not validate real model inference or biological performance.

For reference-subtracted experimental BLI traces, use the columns in `examples/bli_template.csv`:

```bash
conda activate binderflow
cd Adaptyv-Bio-Round2-Protein-design-competition/pipeline
python -m binderflow bli /absolute/path/to/traces.csv --out runs/bli_fit
```

No wet-lab measurements are generated by this pipeline.

## Versions and licensing

The installer pins upstream source revisions in [upstream_revisions.json](docs/upstream_revisions.json) and the main GPU packages. Transitive dependencies are not fully locked. Installed package lists and checkpoint SHA-256 hashes are recorded locally before inference.

Project code uses the repository's [MIT license](../LICENSE). Downloaded upstream software and model weights retain their own licenses. Model weights, third-party checkouts, local environments, caches, logs, and generated results are excluded from Git.

Sources: [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion), [ProteinMPNN](https://github.com/dauparas/ProteinMPNN), [ColabFold](https://github.com/sokrypton/ColabFold), [PDB 1MOX](https://www.rcsb.org/structure/1MOX).
