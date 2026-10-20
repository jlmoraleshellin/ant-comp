#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --job-name=hisbond

# One batch per task; line: sub-manifest binder target min_dist max_dist max_angle env_radius input_col centers_arg

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required off Modal) — see docs/configuration.md.
sapia_activate SAPIA_ACTIVATE_HISBOND

PY=${PIPELINE_PYTHON:-python}

# Never name a field after a bash special (GROUPS, UID...): under `set -e` the task dies with empty .out.
TASK_FILE=$(echo "$SAPIA_LINE" | cut -f1)
BINDER_CHAIN=$(echo "$SAPIA_LINE" | cut -f2)
TARGET_CHAINS=$(echo "$SAPIA_LINE" | cut -f3)
MIN_DIST=$(echo "$SAPIA_LINE" | cut -f4)
MAX_DIST=$(echo "$SAPIA_LINE" | cut -f5)
MAX_ANGLE=$(echo "$SAPIA_LINE" | cut -f6)
ENV_RADIUS=$(echo "$SAPIA_LINE" | cut -f7)
INPUT_COLUMN=$(echo "$SAPIA_LINE" | cut -f8)
CENTERS_ARG=$(echo "$SAPIA_LINE" | cut -f9)

echo "[$(date +%T)] task $SAPIA_TASK_ID: hisbond on $(wc -l <"$TASK_FILE") designs" \
     "(binder=$BINDER_CHAIN target=$TARGET_CHAINS ${MIN_DIST}<=dist<=$MAX_DIST A" \
     "angle<=$MAX_ANGLE deg, heavy-atom contacts; hydrogens read for the tautomer" \
     "only; env radius $ENV_RADIUS A)"

# A worker crash writes a fallback error TSV per design without one; its rc is re-raised so .exit is honest.
WORKER_RC=0
"$PY" "${SAPIA_TOOL_DIR:?}/hisbond_worker.py" \
        --task-file "$TASK_FILE" \
        --binder-chain "$BINDER_CHAIN" \
        --target-chains "$TARGET_CHAINS" \
        --min-dist "$MIN_DIST" \
        --max-dist "$MAX_DIST" \
        --max-angle "$MAX_ANGLE" \
        --env-radius "$ENV_RADIUS" \
        --input-column "$INPUT_COLUMN" \
        --centers-arg "$CENTERS_ARG" \
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
