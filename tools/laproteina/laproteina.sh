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
LPC_NAME=$(echo "$SAPIA_LINE" | cut -f1)
LPC_CONFIG=$(echo "$SAPIA_LINE" | cut -f2)
LPC_JOB_ID=$(echo "$SAPIA_LINE" | cut -f3)
LPC_OVERRIDES=$(echo "$SAPIA_LINE" | cut -f4)
LPC_DEST=$(echo "$SAPIA_LINE" | cut -f5)

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
