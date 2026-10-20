#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --job-name=chainsel

# Extract one design's chain subset per array task, and write a per-design TSV
# (name, status, path, n_chains, n_res, chains, n_atoms[, chain_n_res,
# chain_mean_b in split mode]) for collect_chainsel.py
# to merge back. CPU-only. The worker needs gemmi, so it runs under
# $PIPELINE_PYTHON (defaulting to `python`).
#
# The manifest carries the selection in its canonical group form
# ('<chain>[+<chain>...]:<out_id>', comma-joined), so a plain --chains run and a
# --merge-groups run are one code path here and in the worker. A --split-groups run
# leaves the group field empty and carries its canonical split spec in field 5
# instead; the worker takes its own split path when that field is non-empty.

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required) — see docs/configuration.md.
sapia_activate SAPIA_ACTIVATE_CHAINSEL

# The interpreter defaults to `python` (made right by the activation hook, or already
# on PATH — it must carry gemmi). Set PIPELINE_PYTHON to point straight at a specific
# interpreter.
PY=${PIPELINE_PYTHON:-python}

NAME=$(echo "$SAPIA_LINE" | cut -f1)
SRC=$(echo "$SAPIA_LINE" | cut -f2)
# NOTE: do NOT name this variable GROUPS -- bash pre-sets GROUPS as a special
# indexed array, the assignment fails with rc=1, and under `set -e` the task dies
# with EMPTY .out and .err. Same hazard: UID, EUID, PPID, PIPESTATUS, SECONDS.
SEL_GROUPS=$(echo "$SAPIA_LINE" | cut -f3)
RENUMBER_FROM=$(echo "$SAPIA_LINE" | cut -f4)
SPLIT_SPEC=$(echo "$SAPIA_LINE" | cut -f5)
OUT_FORMAT=$(echo "$SAPIA_LINE" | cut -f6)
HET_MODE=$(echo "$SAPIA_LINE" | cut -f7)
RENUMBER_MODE=$(echo "$SAPIA_LINE" | cut -f8)

RESULT_TSV="$OUT_DIR/${NAME}.tsv"
OUT_STRUCT="$OUT_DIR/${NAME}.${OUT_FORMAT}"

# store_true flags: pass them only when the manifest asked for them. Written as
# `if` blocks, not `[[ ... ]] && ...`, which would return 1 under `set -e` and kill
# the task whenever the flag is off (i.e. on the defaults).
EXTRA=()
if [[ "$HET_MODE" == "keep" ]]; then
    EXTRA+=(--keep-het)
fi
if [[ "$RENUMBER_MODE" == "renumber" ]]; then
    EXTRA+=(--renumber)
fi

echo "[$(date +%T)] task $SAPIA_TASK_ID: chainsel ${SEL_GROUPS:-split $SPLIT_SPEC} for $NAME"

# The worker lives next to this .sh; SAPIA_TOOL_DIR (exported by the driver) points
# there, independent of the submit cwd. If the worker itself crashes (before it can
# record error-as-data), write a fallback error TSV so collect still sees this design.
WORKER_RC=0
"$PY" "${SAPIA_TOOL_DIR:?}/chainsel_worker.py" \
        --name "$NAME" \
        --src "$SRC" \
        --groups "$SEL_GROUPS" \
        --renumber-from "$RENUMBER_FROM" \
        --split-groups "$SPLIT_SPEC" \
        --out-format "$OUT_FORMAT" \
        ${EXTRA[@]+"${EXTRA[@]}"} \
        --out "$OUT_STRUCT" \
        --result-tsv "$RESULT_TSV" || WORKER_RC=$?

if [ "$WORKER_RC" -ne 0 ]; then
    printf 'name\tstatus\tpath\n' >"$RESULT_TSV"
    printf '%s\terror: worker crashed\t\n' "$NAME" >>"$RESULT_TSV"
fi

exit "$WORKER_RC"
