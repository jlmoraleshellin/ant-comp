#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --job-name=posrmsd

# CA RMSD of structure A onto structure B over explicit chain positions, one
# design per array task. Writes <OUT_DIR>/<prefix>/<name>.tsv for
# collect_posrmsd.py and <name>_sup.cif (A moved into B's frame). CPU-only; the
# worker needs gemmi + numpy and runs under $PIPELINE_PYTHON (default `python`).

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required) — see docs/configuration.md.
sapia_activate SAPIA_ACTIVATE_POSRMSD

PY=${PIPELINE_PYTHON:-python}

# Manifest fields. None of these names shadows a bash special variable (PATH,
# GROUPS, UID, ...), which would fail the assignment under `set -e` and kill the
# task with empty .out and .err.
NAME=$(printf '%s\n' "$SAPIA_LINE" | cut -f1)
STRUCT_A=$(printf '%s\n' "$SAPIA_LINE" | cut -f2)
STRUCT_B=$(printf '%s\n' "$SAPIA_LINE" | cut -f3)
# May be empty (the default): an empty field mid-line still holds its place.
SEL_CHAIN_A=$(printf '%s\n' "$SAPIA_LINE" | cut -f4)
SEL_CHAIN_B=$(printf '%s\n' "$SAPIA_LINE" | cut -f5)
SEG_SPEC=$(printf '%s\n' "$SAPIA_LINE" | cut -f6)
CMP_PREFIX=$(printf '%s\n' "$SAPIA_LINE" | cut -f7)

RESULTS_DIR="$OUT_DIR/$CMP_PREFIX"
mkdir -p "$RESULTS_DIR"
RESULT_TSV="$RESULTS_DIR/${NAME}.tsv"
OUT_STRUCT="$RESULTS_DIR/${NAME}_sup.cif"

echo "[$(date +%T)] task $SAPIA_TASK_ID: posrmsd $CMP_PREFIX for $NAME"

# If the worker itself crashes (before it can record error-as-data), write a
# fallback error TSV so collect still sees this design. Its return code is
# captured explicitly and re-raised, so the .exit file means what it says.
WORKER_RC=0
"$PY" "${SAPIA_TOOL_DIR:?}/posrmsd_worker.py" \
        --name "$NAME" \
        --a "$STRUCT_A" \
        --b "$STRUCT_B" \
        --chain-a "$SEL_CHAIN_A" \
        --chain-b "$SEL_CHAIN_B" \
        --segments "$SEG_SPEC" \
        --out-struct "$OUT_STRUCT" \
        --result-tsv "$RESULT_TSV" || WORKER_RC=$?

if [ "$WORKER_RC" -ne 0 ]; then
    printf 'name\tstatus\tsup_path\n' >"$RESULT_TSV"
    printf '%s\terror: worker crashed\t\n' "$NAME" >>"$RESULT_TSV"
fi

exit "$WORKER_RC"
