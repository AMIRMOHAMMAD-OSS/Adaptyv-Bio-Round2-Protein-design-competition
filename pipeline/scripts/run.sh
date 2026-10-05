set -eo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
init_conda
conda activate binderflow
load_tools
set -u
cd "$project_root"
mode="${1:-smoke}"
case "$mode" in
    smoke|full) config="configs/$mode.yaml" ;;
    *) printf '%s\n' 'Usage: bash scripts/run.sh smoke|full' >&2; exit 2 ;;
esac
[[ -s "$BINDER_TOOLS_RECORD" ]] || { printf '%s\n' 'Run bash scripts/download_models.sh first.' >&2; exit 1; }
mkdir -p logs
python -m binderflow plan "$config"
python -m binderflow run "$config" 2>&1 | tee "logs/${mode}_$(date +%Y%m%d_%H%M%S).log"
