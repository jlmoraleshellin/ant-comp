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

# --- manifest-line guard -----------------------------------------------------
# The prelude does `SAPIA_LINE="$(sed -n "${SAPIA_TASK_ID}p" "$MANIFEST")"`, and
# `sed -n Np` on a file that does not (yet) show N lines EXITS 0 WITH EMPTY
# OUTPUT. `set -e` does not trip. Three empty fields then reach the worker,
# `--config ''` becomes Path('') == '.', and the task dies with a baffling
# `IsADirectoryError: Is a directory: '.'`.
#
# Measured 2026-10-02, run 20261002_143419_dimer_phase2: 18 of 37 tasks hit this
# at once. The manifest on the Volume was provably intact afterwards — 37 lines,
# 3 fields each, zero NUL bytes — so this is a READ-side visibility failure in
# the task container, not a corrupt writer. prosapia's publish_manifest already
# uses batch_upload against exactly this hazard (its docstring names it) and it
# was not sufficient. One task re-read the file ~2 minutes later and STILL saw
# nothing, so the stale view can outlive a short backoff.
#
# So: retry briefly in case the view settles, then FAIL LOUDLY. A wrong-but-
# plausible run is far worse than a dead one, and the silent version of this bug
# cost a full 37-task batch.
sapia_read_manifest_line() {
    local attempt line nlines
    for attempt in 1 2 3 4 5; do
        line=$(sed -n "${SAPIA_TASK_ID}p" "$MANIFEST" || true)
        if [ -n "$line" ]; then
            printf '%s' "$line"
            return 0
        fi
        nlines=$(wc -l <"$MANIFEST" 2>/dev/null || echo "?")
        echo "[$(date +%T)] WARNING: manifest line $SAPIA_TASK_ID empty on attempt" \
             "$attempt/5 (manifest '$MANIFEST' shows $nlines lines); retrying" >&2
        sleep $((attempt * 5))
    done
    return 1
}

if [ -z "${SAPIA_LINE:-}" ]; then
    echo "[$(date +%T)] manifest line $SAPIA_TASK_ID came back EMPTY; re-reading" >&2
    SAPIA_LINE=$(sapia_read_manifest_line) || {
        echo "FATAL: task $SAPIA_TASK_ID could not read line $SAPIA_TASK_ID of" \
             "'$MANIFEST' after 5 attempts. The manifest is present but this" \
             "container does not see line $SAPIA_TASK_ID (stale Volume view)." \
             "Refusing to run with empty inputs." >&2
        exit 1
    }
fi

# NOTE: manifest-field locals are prefixed. Bash pre-sets names like GROUPS, UID
# and PIPESTATUS; assigning a command substitution to one fails with rc=1 and,
# under `set -e`, kills the task with EMPTY .out and .err.
RPX_NAME=$(echo "$SAPIA_LINE" | cut -f1)
RPX_INPUT=$(echo "$SAPIA_LINE" | cut -f2)
RPX_CONFIG=$(echo "$SAPIA_LINE" | cut -f3)

# Field COUNT first: `cut -f2` on a line with NO TABS returns the whole line, not
# an empty string, so a one-field line would otherwise sail past the empty checks
# below with name == input == config and fail later as a confusing "missing file".
RPX_NFIELDS=$(printf '%s' "$SAPIA_LINE" | awk -F'\t' '{print NF}')
if [ "$RPX_NFIELDS" -ne 3 ]; then
    echo "FATAL: task $SAPIA_TASK_ID read a manifest line with $RPX_NFIELDS" \
         "field(s), expected 3 (name, scaffold, config). line='$SAPIA_LINE'" >&2
    exit 1
fi

# Every field must be present AND the paths must exist. Catches a short line and
# a half-visible Volume as separate, named failures rather than one confusing one.
if [ -z "$RPX_NAME" ] || [ -z "$RPX_INPUT" ] || [ -z "$RPX_CONFIG" ]; then
    echo "FATAL: task $SAPIA_TASK_ID read a malformed manifest line." \
         "name='$RPX_NAME' input='$RPX_INPUT' config='$RPX_CONFIG'" >&2
    exit 1
fi
for _f in "$RPX_INPUT" "$RPX_CONFIG"; do
    if [ ! -s "$_f" ]; then
        echo "FATAL: task $SAPIA_TASK_ID ($RPX_NAME): '$_f' is missing or empty." \
             "The manifest line parsed, but this container cannot see the file." >&2
        exit 1
    fi
done
# --- end manifest-line guard -------------------------------------------------

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
