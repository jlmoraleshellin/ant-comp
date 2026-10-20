#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --job-name=propka

# One batch per task; line: sub-manifest binder_chains input_col centres_source keep_inputs ph_values release_factors reference_ph free_only
# Sub-manifest row: name structure centres expect_col expect_value (the last three may be empty).

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site activation (required off Modal): `python` / $PIPELINE_PYTHON must be able to `import propka`.
sapia_activate SAPIA_ACTIVATE_PROPKA

PY=${PIPELINE_PYTHON:-python}

# Never name a field after a bash special (GROUPS, UID...): under `set -e` the task dies with empty .out.
# Every field is non-empty (the builder always writes all nine), so none can vanish from `cut`.
field() { echo "$SAPIA_LINE" | cut -f"$1"; }
TASK_FILE=$(field 1)
BINDER_CHAINS=$(field 2)
INPUT_COLUMN=$(field 3)
CENTRES_SOURCE=$(field 4)
KEEP_INPUTS=$(field 5)
PH_VALUES=$(field 6)
RELEASE_FACTORS=$(field 7)
REFERENCE_PH=$(field 8)
FREE_ONLY=$(field 9)
: "${FREE_ONLY:?manifest field 9 (free_only) is empty: re-submit with the current run_propka.py}"

FLAGS=()
[ "$KEEP_INPUTS" = "1" ] && FLAGS+=(--keep-propka-input)
[ "$FREE_ONLY" = "1" ] && FLAGS+=(--free-only)

echo "[$(date +%T)] task $SAPIA_TASK_ID: propka on $(wc -l <"$TASK_FILE") designs" \
     "(binder chains $BINDER_CHAINS, centres $CENTRES_SOURCE; free_only=$FREE_ONLY;" \
     "F at pH $PH_VALUES vs reference $REFERENCE_PH, release factors $RELEASE_FACTORS)"

# A worker crash writes a fallback error TSV per design without one; its rc is re-raised so .exit is honest.
WORKER_RC=0
"$PY" "${SAPIA_TOOL_DIR:?}/propka_worker.py" \
        --task-file "$TASK_FILE" \
        --binder-chains "$BINDER_CHAINS" \
        --input-column "$INPUT_COLUMN" \
        --centres-source "$CENTRES_SOURCE" \
        --ph-values "$PH_VALUES" \
        --release-factors "$RELEASE_FACTORS" \
        --reference-ph "$REFERENCE_PH" \
        ${FLAGS[@]+"${FLAGS[@]}"} \
        --out-dir "$OUT_DIR" || WORKER_RC=$?

if [ "$WORKER_RC" -ne 0 ]; then
    while IFS=$'\t' read -r NAME _REST || [ -n "$NAME" ]; do
        [ -n "$NAME" ] || continue
        if [ ! -f "$OUT_DIR/${NAME}.tsv" ]; then
            printf 'name\tstatus\tpath\n%s\terror: worker crashed\t\n' "$NAME" \
                >"$OUT_DIR/${NAME}.tsv"
        fi
    done <"$TASK_FILE"
fi

exit "$WORKER_RC"
