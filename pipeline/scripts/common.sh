project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

init_conda() {
    local conda_base
    if command -v conda >/dev/null 2>&1; then
        conda_base="$(conda info --base)"
    elif [[ -x "$HOME/miniconda3/bin/conda" ]]; then
        conda_base="$HOME/miniconda3"
    else
        printf '%s\n' 'Conda is required. Open a Conda-enabled Ubuntu terminal.' >&2
        return 1
    fi
    source "$conda_base/etc/profile.d/conda.sh"
}

load_tools() {
    if [[ ! -r "$project_root/tools.env" ]]; then
        printf '%s\n' 'Run bash scripts/install.sh or bash scripts/configure.sh first.' >&2
        return 1
    fi
    source "$project_root/tools.env"
    : "${RFD_ROOT:?}" "${MPNN_ROOT:?}" "${COLABFOLD_BATCH:?}"
    : "${COLABFOLD_PARAMS:?}" "${BINDER_TOOLS_RECORD:?}"
    export DGLBACKEND=pytorch
    export XLA_PYTHON_CLIENT_PREALLOCATE=false
    export XLA_PYTHON_CLIENT_ALLOCATOR=platform
    export XLA_PYTHON_CLIENT_MEM_FRACTION=0.85
    export TF_FORCE_UNIFIED_MEMORY=0
}
