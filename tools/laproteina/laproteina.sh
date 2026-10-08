#!/bin/bash
#SBATCH --job-name=laproteina
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE. Our manifest
# line is one generation job: its design-group name, the LPC config, the job id, the
# Hydra overrides, and the directory this task must leave its output in.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required) — see docs/configuration.md.
# Must put `complexa` on PATH and set COMPLEXA_HOME to the LPC checkout.
# A no-op under Modal, where modal_image.py provides both.
sapia_activate SAPIA_ACTIVATE_LAPROTEINA

# Manifest columns (tab-separated). Names are prefixed to stay clear of bash's own
# special variables — a plain assignment to one of those fails under `set -e` with
# no output at all. See the trap in the authoring-a-tool skill.
# The manifest can still be invisible to this container when the task starts.
# Every volume mount is created with allow_background_commits=True, so although
# the submitter publishes the manifest server-side before spawning, a task that
# starts immediately (warm image, no build to wait through) can read the file
# through a mount view that has not caught up -- getting a sparse read, which
# shows up as the prelude's "ignored null byte in input" and an empty field 1.
#
# Observed: of two otherwise identical runs submitted a minute apart, the first
# read a null-filled manifest and died, the second read it correctly. So: if the
# prelude's line looks empty, re-read our own line straight from $MANIFEST,
# backing off, before giving up.
_lpc_read_line() {
    local i line
    line="$SAPIA_LINE"
    for i in 1 2 3 4 5 6; do
        # tr -d '\0' so a partially-materialised read cannot smuggle NULs into
        # the variable and silently truncate every field.
        line=$(sed -n "${SAPIA_TASK_ID}p" "${MANIFEST:?}" 2>/dev/null | tr -d '\0')
        [ -n "$(printf '%s' "$line" | cut -f1)" ] && { printf '%s' "$line"; return 0; }
        echo "  manifest line $SAPIA_TASK_ID not readable yet (attempt $i), waiting ${i}0s..." >&2
        sleep "${i}0"
    done
    printf '%s' "$line"
}

LPC_NAME=$(printf '%s' "$SAPIA_LINE" | tr -d '\0' | cut -f1)
if [ -z "$LPC_NAME" ]; then
    SAPIA_LINE=$(_lpc_read_line)
    LPC_NAME=$(printf '%s' "$SAPIA_LINE" | cut -f1)
fi
LPC_CONFIG=$(printf '%s' "$SAPIA_LINE" | cut -f2)
LPC_JOB_ID=$(printf '%s' "$SAPIA_LINE" | cut -f3)
LPC_OVERRIDES=$(printf '%s' "$SAPIA_LINE" | cut -f4)
LPC_DEST=$(printf '%s' "$SAPIA_LINE" | cut -f5)

LPC_HOME=${COMPLEXA_HOME:?set COMPLEXA_HOME to the Proteina-Complexa checkout}
LPC_BIN=${COMPLEXA_BIN:-complexa}

# A short manifest line would leave these empty, and an empty LPC_DEST would reach
# `rm -rf` below. Fail here, where the message says what is actually wrong.
: "${LPC_NAME:?manifest field 1 (name) is empty}"
: "${LPC_DEST:?manifest field 5 (dest) is empty}"

echo "[$(date +%T)] task $SAPIA_TASK_ID: LPC generate '$LPC_NAME' (job $LPC_JOB_ID)"
echo "  config:    $LPC_CONFIG"
echo "  overrides: ${LPC_OVERRIDES:-(none)}"
echo "  dest:      $LPC_DEST"

# LPC resolves its config, checkpoints and target paths relative to its own
# checkout, and writes to a hardcoded ./inference/<...> under the cwd. So run from
# there and relocate the result afterwards, rather than fighting the relative paths.
cd "$LPC_HOME"

# --- the target registry -----------------------------------------------------
# LPC resolves --task-name against configs/targets/targets_dict.yaml *inside this
# checkout*. On Modal the checkout is a fresh clone, so a target added for this
# campaign is simply absent and the run dies at config resolution. Overlay the
# assets volume (targets_dict.yaml plus data/target_data/) if one is mounted.
if [ -n "${LPC_ASSETS:-}" ] && [ -d "${LPC_ASSETS}" ]; then
    echo "  overlaying assets from $LPC_ASSETS"
    cp -r "${LPC_ASSETS}/." "$LPC_HOME/"
fi

# --- checkpoints -------------------------------------------------------------
# The shipped config points at ./ckpts inside the checkout. Off Modal that may be
# right; on Modal they live on a volume. Redirect only when told to, so a local
# run that already has them keeps working.
if [ -n "${COMPLEXA_CKPT_DIR:-}" ]; then
    echo "  checkpoints from $COMPLEXA_CKPT_DIR"
    LPC_OVERRIDES="$LPC_OVERRIDES ++ckpt_path=${COMPLEXA_CKPT_DIR}"
    LPC_OVERRIDES="$LPC_OVERRIDES ++autoencoder_ckpt_path=${COMPLEXA_CKPT_DIR}/complexa_ae.ckpt"
fi

# --- the silent-skip trap ----------------------------------------------------
# generate.py exits 0 immediately if ./inference/results_<config_name>_<job_id>.csv
# already exists ("Results already exist ... Exiting generate.py"). That file is
# named by job id, NOT by run_name, so a second run reusing the same --job-id would
# silently produce nothing and still look successful. Clear this job's marker first;
# it belongs to a previous run, never to a concurrent one, because job ids are
# unique within a run.
find ./inference -maxdepth 1 -name "results_*_${LPC_JOB_ID}.csv" -print -delete 2>/dev/null || true

# Hydra swallows the real cause of an instantiation failure, reporting only
# "Error locating target '<dotted.path>'" with no hint as to which import inside
# that module actually failed. On a remote task that is the difference between a
# diagnosis and another build cycle, so always ask for the chained traceback.
export HYDRA_FULL_ERROR=1

# The pair representation is O(L^2) per sample and a binder-on-trimer complex is
# ~520 residues, so a full batch is the memory bottleneck, not the model. Measured
# on a 40 GiB A100: batch 2 is comfortable, batch 16 dies asking for 4.86 GiB with
# 2.93 free. expandable_segments reclaims the "reserved but unallocated" slack the
# allocator otherwise strands (8.08 GiB of it in that failure) -- it does not
# raise the ceiling, so keep --set generation.dataloader.batch_size sane too.
export PYTORCH_CUDA_ALLOC_CONF=${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}

# LPC_OVERRIDES is intentionally unquoted: it is a space-separated list of Hydra
# tokens (++key=value) that must each become its own argv entry.
"$LPC_BIN" generate "$LPC_CONFIG" --job-id "$LPC_JOB_ID" --verbose $LPC_OVERRIDES

# The output root is ./inference/<config_name>_<task_name>_<run_name>, where
# config_name comes from Hydra's own resolution. Rather than reproduce that logic,
# find the directory by the run token, which the submitter set to the design-group
# name precisely so it would be unique and findable.
LPC_SRC=$(find ./inference -maxdepth 1 -type d -name "*_${LPC_NAME}" -print -quit)

if [ -z "$LPC_SRC" ] || [ ! -d "$LPC_SRC" ]; then
    echo "ERROR: no output directory matching '*_${LPC_NAME}' under $LPC_HOME/inference" >&2
    echo "       LPC exited 0 but produced nothing findable. Check its log above:" >&2
    echo "       a missing target, an unresolvable \${oc.env:...} or a stale results" >&2
    echo "       CSV all end here." >&2
    exit 1
fi

echo "  found:     $LPC_SRC"
mkdir -p "$(dirname "$LPC_DEST")"
rm -rf "$LPC_DEST"
# Move rather than copy: the samples are the only artefact and leaving a second copy
# inside the LPC checkout would grow without bound across runs.
mv "$LPC_SRC" "$LPC_DEST"

LPC_N=$(find "$LPC_DEST" -name '*.pdb' | wc -l)
echo "[$(date +%T)] task $SAPIA_TASK_ID: done — $LPC_N pdb(s) in $LPC_DEST"

if [ "$LPC_N" -eq 0 ]; then
    echo "ERROR: job produced no structures" >&2
    exit 1
fi
