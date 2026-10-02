#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --job-name=rpxdock

# Dock one scaffold per array task into a one-component symmetric architecture,
# and write a per-design TSV (one row per kept dock) plus the dock structures for
# collect_rpxdock.py to mint child rows from. CPU-only; the manifest builder sets
# gpus_per_task = 0.
#
# The manifest line is (name, scaffold pdb, per-design config json). Everything
# that varies per run lives in the config, so this script stays a thin wrapper.

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required) — see docs/configuration.md.
# Must put an interpreter carrying `rpxdock` (and its compiled extensions) on PATH.
sapia_activate SAPIA_ACTIVATE_RPXDOCK

# The interpreter defaults to `python` (made right by the activation hook, or
# already on PATH — under Modal the image provides it). Set PIPELINE_PYTHON to
# point straight at a specific interpreter.
PY=${PIPELINE_PYTHON:-python}

# NOTE: manifest-field locals are prefixed. Bash pre-sets names like GROUPS, UID
# and PIPESTATUS; assigning a command substitution to one fails with rc=1 and,
# under `set -e`, kills the task with EMPTY .out and .err.
RPX_NAME=$(echo "$SAPIA_LINE" | cut -f1)
RPX_INPUT=$(echo "$SAPIA_LINE" | cut -f2)
RPX_CONFIG=$(echo "$SAPIA_LINE" | cut -f3)

RESULT_TSV="$OUT_DIR/${RPX_NAME}.tsv"

echo "[$(date +%T)] task $SAPIA_TASK_ID: rpxdock on $RPX_NAME ($RPX_INPUT)"

# cppimport compiles RPXdock's C++ extensions on first import and writes them next
# to the sources. The Modal image pre-builds them, so nothing should compile here;
# if a site install has not, this makes the write target explicit rather than
# depending on the submit cwd.
export CPPIMPORT_RELEASE_MODE=${CPPIMPORT_RELEASE_MODE:-0}

# The worker records its own errors as data (status 'error: ...') and exits 1. If
# it dies before it can do even that, leave a fallback TSV so collect still sees
# this scaffold said something.
if ! "$PY" "${SAPIA_TOOL_DIR:?}/rpxdock_worker.py" \
        --name "$RPX_NAME" \
        --input "$RPX_INPUT" \
        --config "$RPX_CONFIG" \
        --out-dir "$OUT_DIR"; then
    if [ ! -s "$RESULT_TSV" ]; then
        printf 'name\tparent\tstatus\tpath\n' >"$RESULT_TSV"
        printf '\t%s\terror: worker crashed\t\n' "$RPX_NAME" >>"$RESULT_TSV"
    fi
    exit 1
fi
