set -eo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
init_conda
set -u
tools_root="${1:-$project_root/external}"
tools_root="$(cd "$tools_root" && pwd)"
[[ -f "$tools_root/RFdiffusion/scripts/run_inference.py" ]]
[[ -f "$tools_root/ProteinMPNN/protein_mpnn_run.py" ]]
cf_bin="$(conda run -n binder-colabfold-wsl python -c 'import sys; from pathlib import Path; print(Path(sys.executable).parent / "colabfold_batch")')"
cf_params="$(conda run -n binder-colabfold-wsl python -c 'from colabfold.download import default_data_dir; print(default_data_dir)')"
[[ -x "$cf_bin" ]]
rfd_library_path="$(conda run -n binder-rfd-wsl python -c 'import sys, torch; from pathlib import Path; p=Path(torch.__file__).parent; dirs=[p.parent/"nvidia"/n/"lib" for n in ("cuda_runtime","cusparse","curand")]; dirs += [p/"lib",Path(sys.prefix)/"lib","/usr/lib/wsl/lib"]; print(":".join(str(d) for d in dirs if Path(d).is_dir()))')"
conda env config vars set -n binder-rfd-wsl DGLBACKEND=pytorch "LD_LIBRARY_PATH=$rfd_library_path"
mkdir -p "$project_root/.local"
printf 'export RFD_ROOT=%q\nexport MPNN_ROOT=%q\nexport COLABFOLD_BATCH=%q\nexport COLABFOLD_PARAMS=%q\nexport BINDER_TOOLS_RECORD=%q\n' \
    "$tools_root/RFdiffusion" "$tools_root/ProteinMPNN" "$cf_bin" "$cf_params" "$project_root/.local/tool_versions.json" \
    > "$project_root/tools.env"
printf 'Configured: %s\n' "$project_root/tools.env"
