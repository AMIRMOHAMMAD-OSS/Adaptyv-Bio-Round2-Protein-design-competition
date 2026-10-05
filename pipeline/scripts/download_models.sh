set -eo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
init_conda
conda activate binderflow
load_tools
set -u
cd "$project_root"
command -v curl >/dev/null
command -v flock >/dev/null
mkdir -p logs .local
exec 9>".local/setup.lock"
flock -n 9 || { printf '%s\n' 'Another model setup is running.' >&2; exit 1; }
[[ -s "$MPNN_ROOT/vanilla_model_weights/v_48_020.pt" ]]
python scripts/download_models.py --rfd-root "$RFD_ROOT" --colabfold-params "$COLABFOLD_PARAMS" \
    2>&1 | tee "logs/download_$(date +%Y%m%d_%H%M%S).log"
python scripts/fetch_example.py
record_next="${BINDER_TOOLS_RECORD}.next"
python scripts/record_tools.py --out "$record_next"
if [[ -f "$BINDER_TOOLS_RECORD" ]]; then
    if cmp -s "$record_next" "$BINDER_TOOLS_RECORD"; then
        rm -- "$record_next"
    else
        printf 'Tool installation changed. Review %s and %s before starting a new run.\n' "$BINDER_TOOLS_RECORD" "$record_next" >&2
        exit 1
    fi
else
    mv -- "$record_next" "$BINDER_TOOLS_RECORD"
fi
printf '%s\n' 'Models are ready. Run bash scripts/run.sh smoke next.'
