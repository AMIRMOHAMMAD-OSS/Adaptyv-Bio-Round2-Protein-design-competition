set -eo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
init_conda
set -u
for executable in git curl nvidia-smi; do command -v "$executable" >/dev/null; done
[[ "$(uname -m)" == x86_64 ]] || { printf '%s\n' 'This setup requires Ubuntu x86_64.' >&2; exit 1; }
tools_root="${1:-$project_root/external}"
mkdir -p "$tools_root"
tools_root="$(cd "$tools_root" && pwd)"
[[ "$tools_root" != *' '* ]] || { printf '%s\n' 'Use a tools directory without spaces.' >&2; exit 1; }
export PIP_NO_CACHE_DIR=1
export PIP_DEFAULT_TIMEOUT=120
export PIP_RETRIES=10
export DGLBACKEND=pytorch

clone_at() {
    local url="$1" destination="$2" revision="$3"
    if [[ ! -d "$destination/.git" ]]; then git clone "$url" "$destination"; fi
    if [[ -n "$(git -C "$destination" status --porcelain --untracked-files=no)" ]]; then
        printf 'Modified checkout: %s\n' "$destination" >&2
        return 1
    fi
    if ! git -C "$destination" cat-file -e "$revision^{commit}" 2>/dev/null; then
        git -C "$destination" fetch origin "$revision"
    fi
    git -C "$destination" checkout --detach "$revision"
}

ensure_env() {
    local name="$1" version="$2"
    if ! conda run -n "$name" python -c 'import sys' >/dev/null 2>&1; then
        conda create -y -n "$name" --override-channels -c conda-forge "python=$version" pip
    fi
    conda run -n "$name" python -c 'import sys; assert ".".join(map(str,sys.version_info[:2])) == sys.argv[1], "Wrong Python version in existing environment"' "$version"
}

clone_at https://github.com/RosettaCommons/RFdiffusion.git "$tools_root/RFdiffusion" 86507b6538f51fce57b5a72477165f03999ed7ae
clone_at https://github.com/dauparas/ProteinMPNN.git "$tools_root/ProteinMPNN" 8907e6671bfbfc92303b5f79c4b5e6ce47cdef57
ensure_env binderflow 3.11
conda run --no-capture-output -n binderflow python -m pip install -e "$project_root"
ensure_env binder-rfd-wsl 3.9
conda run --no-capture-output -n binder-rfd-wsl python -m pip install 'pip<26' 'setuptools<81' wheel
conda run --no-capture-output -n binder-rfd-wsl python -m pip install \
    'numpy==1.23.5' 'scipy==1.10.1' 'protobuf==3.20.3' \
    'torch==1.12.1+cu116' --extra-index-url https://download.pytorch.org/whl/cu116
conda run --no-capture-output -n binder-rfd-wsl python -m pip install \
    'dgl==1.0.2+cu116' -f https://data.dgl.ai/wheels/cu116/repo.html \
    'e3nn==0.3.3' 'wandb==0.12.0' 'pynvml==11.0.0' \
    'decorator==5.1.0' 'hydra-core==1.3.2' 'pyrsistent==0.19.3' \
    'nvidia-cuda-runtime-cu11==11.7.99' 'nvidia-cusparse-cu11==11.7.4.91' \
    'nvidia-curand-cu11==10.2.10.91' \
    'git+https://github.com/NVIDIA/dllogger@0478734ff7be75adde8d160e04872664d1c62e5f'
conda run --no-capture-output -n binder-rfd-wsl python -m pip install "$tools_root/RFdiffusion/env/SE3Transformer"
conda run --no-capture-output -n binder-rfd-wsl python -m pip install -e "$tools_root/RFdiffusion" --no-deps
ensure_env binder-mpnn-wsl 3.10
conda run --no-capture-output -n binder-mpnn-wsl python -m pip install 'numpy==1.26.4'
conda run --no-capture-output -n binder-mpnn-wsl python -m pip install 'torch==2.5.1' --index-url https://download.pytorch.org/whl/cu121
ensure_env binder-colabfold-wsl 3.11
conda run --no-capture-output -n binder-colabfold-wsl python -m pip install \
    'colabfold[alphafold] @ git+https://github.com/sokrypton/ColabFold.git@efbf31c37cedb38cd09c69c1b991910a9866480e' \
    'jax[cuda12]==0.6.2' 'numpy==2.1.3'
for environment in binderflow binder-rfd-wsl binder-mpnn-wsl binder-colabfold-wsl; do
    conda run --no-capture-output -n "$environment" python -m pip check
done
bash "$project_root/scripts/configure.sh" "$tools_root"
bash "$project_root/scripts/check_gpu.sh"
printf '%s\n' 'Installation complete. Run bash scripts/download_models.sh next.'
