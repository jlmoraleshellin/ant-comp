#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --job-name=epitope

# Score one batch of designs per array task: which target residues the binder
# touches, which of ITS OWN residues do the touching, and how much of the requested
# epitope it covers. The manifest line points at a sub-manifest
# (name<TAB>structure<TAB>hotspots<TAB>sequence per design; the sequence field is
# empty unless --sequence-column was given, and the structure field is empty when the
# builder could resolve no path, which the worker turns into an error row). The
# worker writes <name>.tsv and <name>_per_residue.tsv per design for
# collect_epitope.py. CPU-only. The worker needs gemmi + numpy, so it runs under
# $PIPELINE_PYTHON (default `python`).

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required) — see docs/configuration.md.
sapia_activate SAPIA_ACTIVATE_EPITOPE

PY=${PIPELINE_PYTHON:-python}

# NOTE: never name a manifest field after a bash special variable (GROUPS, UID,
# PIPESTATUS, ...) — the assignment fails with rc=1 and, under `set -e`, the task
# dies with EMPTY .out and .err.
TASK_FILE=$(echo "$SAPIA_LINE" | cut -f1)
DESIGN_CHAINS=$(echo "$SAPIA_LINE" | cut -f2)
TARGET_CHAINS=$(echo "$SAPIA_LINE" | cut -f3)
CONTACT_CUTOFF=$(echo "$SAPIA_LINE" | cut -f4)
# Appended 2026-09-30. Defaulted, so a manifest written by an older build (4 fields)
# still runs with the original behaviour instead of passing empty strings.
SEQ_SOURCE=$(echo "$SAPIA_LINE" | cut -f5)
SEQ_SOURCE=${SEQ_SOURCE:-structure}
INPUT_COLUMN=$(echo "$SAPIA_LINE" | cut -f6)
INPUT_COLUMN=${INPUT_COLUMN:-the input structure column}
# Appended 2026-10-08. Defaulted, so a 6-field manifest from an older build still
# runs (with the default overlap cutoff) instead of passing an empty string.
CLASH_CUTOFF=$(echo "$SAPIA_LINE" | cut -f7)
CLASH_CUTOFF=${CLASH_CUTOFF:-2.5}

echo "[$(date +%T)] task $SAPIA_TASK_ID: epitope on $(wc -l <"$TASK_FILE") designs" \
     "(design=$DESIGN_CHAINS target=$TARGET_CHAINS cutoff=$CONTACT_CUTOFF A" \
     "identities=$SEQ_SOURCE)"

# If the worker itself crashes (before it can record error-as-data), write a fallback
# error TSV for every design of this task that has none, so collect still sees them.
# The worker's return code is captured EXPLICITLY and re-raised at the end: with a
# bare `if ! worker; then <fallback>; fi` the script's exit status would be the status
# of the FALLBACK's last command, so a task could exit 1 with every row written
# correctly (or, worse, exit 0 after a genuine crash). The .exit file is the only
# completion signal on Modal, so it has to mean what it says.
WORKER_RC=0
"$PY" "${SAPIA_TOOL_DIR:?}/epitope_worker.py" \
        --task-file "$TASK_FILE" \
        --design-chains "$DESIGN_CHAINS" \
        --target-chains "$TARGET_CHAINS" \
        --contact-cutoff "$CONTACT_CUTOFF" \
        --clash-cutoff "$CLASH_CUTOFF" \
        --seq-source "$SEQ_SOURCE" \
        --input-column "$INPUT_COLUMN" \
        --out-dir "$OUT_DIR" || WORKER_RC=$?

if [ "$WORKER_RC" -ne 0 ]; then
    while IFS=$'\t' read -r NAME _SRC _HOTSPOTS _SEQ || [ -n "$NAME" ]; do
        [ -n "$NAME" ] || continue
        if [ ! -f "$OUT_DIR/${NAME}.tsv" ]; then
            printf 'name\tstatus\tpath\n%s\terror: worker crashed\t\n' "$NAME" \
                >"$OUT_DIR/${NAME}.tsv"
        fi
    done <"$TASK_FILE"
fi

exit "$WORKER_RC"
