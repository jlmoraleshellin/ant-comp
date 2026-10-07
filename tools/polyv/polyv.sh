#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --job-name=polyv

# Triage one batch of backbones per array task: erase the binder to poly-valine,
# repack only its side chains against a frozen target, and measure interface
# geometry with PyRosetta. The manifest line points at a sub-manifest
# (name<TAB>structure per design); the worker writes <name>.pdb and <name>.tsv per
# design for collect_polyv.py. The worker needs pyrosetta, so it runs under
# $PYROSETTA_PYTHON (default `python`), matching the bundled pyrosetta tool.

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required) — see docs/configuration.md.
sapia_activate SAPIA_ACTIVATE_POLYV

PY=${PYROSETTA_PYTHON:-python}

TASK_FILE=$(echo "$SAPIA_LINE" | cut -f1)
BINDER_CHAINS=$(echo "$SAPIA_LINE" | cut -f2)
TARGET_CHAINS=$(echo "$SAPIA_LINE" | cut -f3)
KEEP_GP=$(echo "$SAPIA_LINE" | cut -f4)
EXTRA_ROT=$(echo "$SAPIA_LINE" | cut -f5)
EXCLUDE=$(echo "$SAPIA_LINE" | cut -f6)
DIST_W=$(echo "$SAPIA_LINE" | cut -f7)

# Absolute, because we are about to leave this directory.
OUT_DIR=$(cd "$OUT_DIR" && pwd)

# Hard requirement: nothing is ever written next to the input structure. Rosetta can
# drop files (ROSETTA_CRASH.log, tracer output) into the CURRENT directory, so the
# task runs from a scratch dir of its own and writes results only under $OUT_DIR.
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR"

# store_true flags as `if` blocks: `[[ ... ]] && ...` returns 1 under `set -e` and
# would kill the task whenever the flag is off.
EXTRA=()
if [[ "$KEEP_GP" == "keep" ]]; then
    EXTRA+=(--keep-gly-pro)
fi

echo "[$(date +%T)] task $SAPIA_TASK_ID: polyv on $(wc -l <"$TASK_FILE") designs"

# If the worker itself crashes (before it can record error-as-data), write a fallback
# error TSV for every design of this task that has none, so collect still sees them.
WORKER_RC=0
"$PY" "${SAPIA_TOOL_DIR:?}/polyv_worker.py" \
        --task-file "$TASK_FILE" \
        --binder-chains "$BINDER_CHAINS" \
        --target-chains "$TARGET_CHAINS" \
        --extra-rotamers "$EXTRA_ROT" \
        --distance-weight "$DIST_W" \
        --exclude-resnames "$EXCLUDE" \
        ${EXTRA[@]+"${EXTRA[@]}"} \
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
