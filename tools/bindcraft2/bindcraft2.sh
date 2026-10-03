#!/bin/bash
#SBATCH --job-name=bindcraft2
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE. Our manifest
# line is one campaign: its design-group name, the mode, the settings JSON the
# submitter wrote (which carries project_folder), and an args file holding the
# run-wide `bindcraft design` flags, one argv token per line.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required) — see docs/configuration.md.
# Must put `bindcraft` on PATH (and, off a shared cache, either point
# BINDCRAFT_AF2_PARAMS at the AlphaFold parameters, or let the cache under
# BINDCRAFT_WEIGHTS / $XDG_CACHE_HOME/bindcraft fill on first use).
sapia_activate SAPIA_ACTIVATE_BINDCRAFT2

# Manifest columns (tab-separated). Names are prefixed to stay clear of bash's own
# special variables — a plain assignment to one of those fails under `set -e` with
# no output at all.
BC2_NAME=$(echo "$SAPIA_LINE" | cut -f1)
BC2_MODE=$(echo "$SAPIA_LINE" | cut -f2)
BC2_SETTINGS=$(echo "$SAPIA_LINE" | cut -f3)
BC2_ARGS=$(echo "$SAPIA_LINE" | cut -f4)

BC2_BIN=${BINDCRAFT_BIN:-bindcraft}

# Warm-up mode: download and verify the AlphaFold parameters into the cache, run no
# campaign. `fetch-weights` exits non-zero when a checkpoint is missing or
# unfinished, which is the point — a dummy campaign would not say so.
if [ "$BC2_MODE" = "fetch-weights" ]; then
    echo "[$(date +%T)] task $SAPIA_TASK_ID: bindcraft fetch-weights"
    "$BC2_BIN" fetch-weights
    exit 0
fi

# The run-wide design flags, one argv token per line. Read as an array rather than
# word-split from a string, so a value carrying spaces or JSON survives intact
# (`--set 'binder_lengths=[70, 90]'` is one token, not two).
BC2_ARGV=()
if [ -f "$BC2_ARGS" ]; then
    mapfile -t BC2_ARGV < "$BC2_ARGS"
fi

echo "[$(date +%T)] task $SAPIA_TASK_ID: bindcraft campaign '$BC2_NAME'"
echo "  settings: $BC2_SETTINGS"
echo "  flags:    ${BC2_ARGV[*]:-(none)}"

# ${arr[@]+"${arr[@]}"} so an empty array is not an unbound-variable error under
# `set -u`.
"$BC2_BIN" design "$BC2_SETTINGS" ${BC2_ARGV[@]+"${BC2_ARGV[@]}"}
