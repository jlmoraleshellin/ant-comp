#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --job-name=ssprofile

# Profile one batch of designs per array task: CA radius of gyration and 3-state DSSP
# fractions of the --design-chains selection. The manifest line points at a
# sub-manifest (name<TAB>structure per design; the structure field is empty when the
# builder could resolve no path, which the worker turns into an error row). The worker
# writes <name>.tsv and <name>_residues.tsv per design for collect_ssprofile.py.
# CPU-only. The worker needs gemmi + numpy + mdtraj, so it runs under
# $PIPELINE_PYTHON (default `python`).

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required) — see docs/configuration.md.
sapia_activate SAPIA_ACTIVATE_SSPROFILE

PY=${PIPELINE_PYTHON:-python}

# NOTE: never name a manifest field after a bash special variable (GROUPS, UID,
# PIPESTATUS, ...) — the assignment fails with rc=1 and, under `set -e`, the task
# dies with EMPTY .out and .err.
TASK_FILE=$(echo "$SAPIA_LINE" | cut -f1)
DESIGN_CHAINS=$(echo "$SAPIA_LINE" | cut -f2)
INPUT_COLUMN=$(echo "$SAPIA_LINE" | cut -f3)
INPUT_COLUMN=${INPUT_COLUMN:-the input structure column}

echo "[$(date +%T)] task $SAPIA_TASK_ID: ssprofile on $(wc -l <"$TASK_FILE") designs" \
     "(chains=$DESIGN_CHAINS)"

# Clear every member's previous result BEFORE the worker runs. Without this, a task
# killed mid-run (OOM, timeout, pre-emption) leaves an earlier run's <name>.tsv in
# place, the fallback below sees it exists and writes nothing, and collect reads the
# stale file as this run's answer -- an 'OK' describing a different input.
while IFS=$'\t' read -r NAME _SRC || [ -n "$NAME" ]; do
    [ -n "$NAME" ] || continue
    rm -f "$OUT_DIR/${NAME}.tsv" "$OUT_DIR/${NAME}_residues.tsv"
done <"$TASK_FILE"

# If the worker itself crashes (before it can record error-as-data), write a fallback
# error TSV for every design of this task that has none, so collect still sees them.
# The worker's return code is captured EXPLICITLY and re-raised at the end, so the
# task's exit status means what it says.
WORKER_RC=0
"$PY" "${SAPIA_TOOL_DIR:?}/ssprofile_worker.py" \
        --task-file "$TASK_FILE" \
        --design-chains "$DESIGN_CHAINS" \
        --input-column "$INPUT_COLUMN" \
        --out-dir "$OUT_DIR" || WORKER_RC=$?

if [ "$WORKER_RC" -ne 0 ]; then
    while IFS=$'\t' read -r NAME _SRC || [ -n "$NAME" ]; do
        [ -n "$NAME" ] || continue
        if [ ! -f "$OUT_DIR/${NAME}.tsv" ]; then
            printf 'name\tstatus\tpath\n%s\terror: worker crashed\t\n' "$NAME" \
                >"$OUT_DIR/${NAME}.tsv"
        fi
    done <"$TASK_FILE"
fi

exit "$WORKER_RC"
