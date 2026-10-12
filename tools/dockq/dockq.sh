#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --job-name=dockq

# DockQ v2 (the published package) of a predicted binder complex against a reference
# complex, for one batch of designs per array task. The manifest line is:
# sub-manifest<TAB>binder<TAB>targets<TAB>ref binder<TAB>ref targets<TAB>
# allowed mismatches<TAB>residue match<TAB>min binder identity<TAB>input column<TAB>
# ref column<TAB>min native contacts. The sub-manifest holds name<TAB>model<TAB>reference per design (either
# path may be empty; the worker turns that into an error row). The worker writes
# <name>.tsv per design for collect_dockq.py. CPU-only; runs under $PIPELINE_PYTHON
# (default `python`), which must import the DockQ package.

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required off Modal) — see docs/configuration.md.
sapia_activate SAPIA_ACTIVATE_DOCKQ

PY=${PIPELINE_PYTHON:-python}

# NOTE: never name a manifest field after a bash special variable (GROUPS, UID,
# PIPESTATUS, ...) — the assignment fails with rc=1 and, under `set -e`, the task
# dies with EMPTY .out and .err.
TASK_FILE=$(echo "$SAPIA_LINE" | cut -f1)
BINDER_CHAIN=$(echo "$SAPIA_LINE" | cut -f2)
TARGET_CHAINS=$(echo "$SAPIA_LINE" | cut -f3)
REF_BINDER_CHAIN=$(echo "$SAPIA_LINE" | cut -f4)
REF_TARGET_CHAINS=$(echo "$SAPIA_LINE" | cut -f5)
ALLOWED_MISMATCHES=$(echo "$SAPIA_LINE" | cut -f6)
RESIDUE_MATCH=$(echo "$SAPIA_LINE" | cut -f7)
MIN_BINDER_IDENTITY=$(echo "$SAPIA_LINE" | cut -f8)
INPUT_COLUMN=$(echo "$SAPIA_LINE" | cut -f9)
REF_COLUMN=$(echo "$SAPIA_LINE" | cut -f10)
# Field 11, appended: defaulted to 1 (DockQ's own rule) for a 10-field manifest.
MIN_NATIVE_CONTACTS=$(echo "$SAPIA_LINE" | cut -f11)
MIN_NATIVE_CONTACTS=${MIN_NATIVE_CONTACTS:-1}

echo "[$(date +%T)] task $SAPIA_TASK_ID: dockq on $(wc -l <"$TASK_FILE") designs" \
     "(model binder=$BINDER_CHAIN targets=$TARGET_CHAINS; reference binder=$REF_BINDER_CHAIN" \
     "targets=$REF_TARGET_CHAINS; residue match=$RESIDUE_MATCH)"

# If the worker crashes before recording error-as-data, write a fallback error TSV for
# every design of this task that has none, so collect still sees them. The worker's
# return code is re-raised at the end, so the task's exit status means what it says.
WORKER_RC=0
"$PY" "${SAPIA_TOOL_DIR:?}/dockq_worker.py" \
        --task-file "$TASK_FILE" \
        --binder-chain "$BINDER_CHAIN" \
        --target-chains "$TARGET_CHAINS" \
        --ref-binder-chain "$REF_BINDER_CHAIN" \
        --ref-target-chains "$REF_TARGET_CHAINS" \
        --allowed-mismatches "$ALLOWED_MISMATCHES" \
        --residue-match "$RESIDUE_MATCH" \
        --min-binder-identity "$MIN_BINDER_IDENTITY" \
        --min-native-contacts "$MIN_NATIVE_CONTACTS" \
        --input-column "$INPUT_COLUMN" \
        --ref-col "$REF_COLUMN" \
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
