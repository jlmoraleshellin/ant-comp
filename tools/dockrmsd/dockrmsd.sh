#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --job-name=dockrmsd

# Docking RMSD of a predicted complex against its reference design, for one batch
# of designs per array task. The manifest line is:
#   sub-manifest<TAB>binder<TAB>targets<TAB>ref binder<TAB>ref targets<TAB>
#   fit atoms<TAB>target match<TAB>min target res<TAB>min target res frac<TAB>
#   min target identity<TAB>fit tolerance<TAB>max target fit rmsd<TAB>
#   pin mapping<TAB>prefix
# The sub-manifest holds name<TAB>model<TAB>reference per design. The worker writes
# <OUT_DIR>/<prefix>/<name>.tsv for collect_dockrmsd.py plus <name>_sup.cif (the
# model moved into the reference frame). CPU-only; the worker needs gemmi + numpy
# and runs under $PIPELINE_PYTHON (default `python`).

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required off Modal) — see docs/configuration.md.
sapia_activate SAPIA_ACTIVATE_DOCKRMSD

PY=${PIPELINE_PYTHON:-python}

# NOTE: never name a manifest field after a bash special variable (GROUPS, UID,
# PIPESTATUS, ...) — the assignment fails with rc=1 and, under `set -e`, the task
# dies with EMPTY .out and .err.
TASK_FILE=$(printf '%s\n' "$SAPIA_LINE" | cut -f1)
BINDER_CHAIN=$(printf '%s\n' "$SAPIA_LINE" | cut -f2)
TARGET_CHAINS=$(printf '%s\n' "$SAPIA_LINE" | cut -f3)
REF_BINDER_CHAIN=$(printf '%s\n' "$SAPIA_LINE" | cut -f4)
REF_TARGET_CHAINS=$(printf '%s\n' "$SAPIA_LINE" | cut -f5)
FIT_ATOMS=$(printf '%s\n' "$SAPIA_LINE" | cut -f6)
TARGET_MATCH=$(printf '%s\n' "$SAPIA_LINE" | cut -f7)
MIN_TARGET_RES=$(printf '%s\n' "$SAPIA_LINE" | cut -f8)
MIN_TARGET_RES_FRAC=$(printf '%s\n' "$SAPIA_LINE" | cut -f9)
MIN_TARGET_IDENTITY=$(printf '%s\n' "$SAPIA_LINE" | cut -f10)
FIT_TOLERANCE=$(printf '%s\n' "$SAPIA_LINE" | cut -f11)
MAX_TARGET_FIT_RMSD=$(printf '%s\n' "$SAPIA_LINE" | cut -f12)
PIN_MAPPING=$(printf '%s\n' "$SAPIA_LINE" | cut -f13)
CMP_PREFIX=$(printf '%s\n' "$SAPIA_LINE" | cut -f14)

RESULTS_DIR="$OUT_DIR/$CMP_PREFIX"
mkdir -p "$RESULTS_DIR"

echo "[$(date +%T)] task $SAPIA_TASK_ID: dockrmsd $CMP_PREFIX on" \
     "$(wc -l <"$TASK_FILE") designs (model binder=$BINDER_CHAIN targets=$TARGET_CHAINS;" \
     "reference binder=$REF_BINDER_CHAIN targets=$REF_TARGET_CHAINS; atoms=$FIT_ATOMS;" \
     "target match=$TARGET_MATCH; pin mapping=$PIN_MAPPING)"

# If the worker crashes before recording error-as-data, write a fallback error TSV
# for every design of this task that has none, so collect still sees them. The
# worker's return code is re-raised, so the .exit file means what it says.
WORKER_RC=0
"$PY" "${SAPIA_TOOL_DIR:?}/dockrmsd_worker.py" \
        --task-file "$TASK_FILE" \
        --binder-chain "$BINDER_CHAIN" \
        --target-chains "$TARGET_CHAINS" \
        --ref-binder-chain "$REF_BINDER_CHAIN" \
        --ref-target-chains "$REF_TARGET_CHAINS" \
        --fit-atoms "$FIT_ATOMS" \
        --target-match "$TARGET_MATCH" \
        --min-target-res "$MIN_TARGET_RES" \
        --min-target-res-frac "$MIN_TARGET_RES_FRAC" \
        --pin-mapping "$PIN_MAPPING" \
        --min-target-identity "$MIN_TARGET_IDENTITY" \
        --fit-tolerance "$FIT_TOLERANCE" \
        --max-target-fit-rmsd "$MAX_TARGET_FIT_RMSD" \
        --results-dir "$RESULTS_DIR" || WORKER_RC=$?

if [ "$WORKER_RC" -ne 0 ]; then
    while IFS=$'\t' read -r NAME _REST || [ -n "$NAME" ]; do
        [ -n "$NAME" ] || continue
        if [ ! -f "$RESULTS_DIR/${NAME}.tsv" ]; then
            printf 'name\tstatus\tsup_path\n%s\terror: worker crashed\t\n' "$NAME" \
                >"$RESULTS_DIR/${NAME}.tsv"
        fi
    done <"$TASK_FILE"
fi

exit "$WORKER_RC"
