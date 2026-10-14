#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --job-name=motifrmsd

# Measure one design's motif CA RMSD against its reference per array task, and
# write a per-design TSV (name, status, aligned_path, metrics) for
# collect_motifrmsd.py to merge back. CPU-only. The worker needs gemmi + numpy,
# so it runs under $PIPELINE_PYTHON (defaulting to `python`).

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required) — see docs/configuration.md.
sapia_activate SAPIA_ACTIVATE_MOTIFRMSD

# The interpreter defaults to `python` (made right by the activation hook, or already
# on PATH — it must carry gemmi + numpy). Set PIPELINE_PYTHON to point straight at a
# specific interpreter.
PY=${PIPELINE_PYTHON:-python}

# Manifest fields. None of these names shadows a bash special variable (GROUPS,
# UID, PIPESTATUS, ...), which would fail the assignment under `set -e` and kill
# the task with empty .out and .err.
NAME=$(echo "$SAPIA_LINE" | cut -f1)
DESIGN=$(echo "$SAPIA_LINE" | cut -f2)
SIDECAR=$(echo "$SAPIA_LINE" | cut -f3)
REF_STRUCT=$(echo "$SAPIA_LINE" | cut -f4)
# May be empty (the default): an empty field in the MIDDLE of the line is still a
# field, so the cut positions below it do not shift.
MOTIF_CHAIN=$(echo "$SAPIA_LINE" | cut -f5)
ALLOW_PARTIAL=$(echo "$SAPIA_LINE" | cut -f6)

RESULT_TSV="$OUT_DIR/${NAME}.tsv"
OUT_PDB="$OUT_DIR/${NAME}_onmotif.pdb"

PARTIAL_FLAG=()
if [ "$ALLOW_PARTIAL" = "yes" ]; then
    PARTIAL_FLAG=(--allow-partial)
fi

echo "[$(date +%T)] task $SAPIA_TASK_ID: motifrmsd for $NAME"

# The worker lives next to this .sh; SAPIA_TOOL_DIR (exported by the driver) points
# there, independent of the submit cwd. If the worker itself crashes (before it can
# record error-as-data), write a fallback error TSV so collect still sees this design.
# The worker's return code is captured EXPLICITLY and re-raised at the end: with a
# bare `if ! worker; then <fallback>; fi` the script's exit status would be the status
# of the FALLBACK's last command, so a task could exit 1 with every row written
# correctly (or, worse, exit 0 after a genuine crash). The .exit file is the only
# completion signal on Modal, so it has to mean what it says.
WORKER_RC=0
"$PY" "${SAPIA_TOOL_DIR:?}/motifrmsd_worker.py" \
        --name "$NAME" \
        --design "$DESIGN" \
        --sidecar "$SIDECAR" \
        --ref "$REF_STRUCT" \
        --motif-chain "$MOTIF_CHAIN" \
        ${PARTIAL_FLAG[@]+"${PARTIAL_FLAG[@]}"} \
        --out-pdb "$OUT_PDB" \
        --result-tsv "$RESULT_TSV" || WORKER_RC=$?

if [ "$WORKER_RC" -ne 0 ]; then
    printf 'name\tstatus\taligned_path\n' >"$RESULT_TSV"
    printf '%s\terror: worker crashed\t\n' "$NAME" >>"$RESULT_TSV"
fi

exit "$WORKER_RC"
