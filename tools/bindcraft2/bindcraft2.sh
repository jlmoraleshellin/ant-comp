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
# line is one campaign: its design-group name, the settings JSON the submitter wrote
# (which carries project_folder), and the run-wide `bindcraft design` flags.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required) — see docs/configuration.md.
# Must put `bindcraft` on PATH (and, off a shared cache, either point
# BINDCRAFT_AF2_PARAMS at the AlphaFold parameters or let ~/.cache/bindcraft fill).
sapia_activate SAPIA_ACTIVATE_BINDCRAFT2

# Manifest columns (tab-separated). Names are prefixed to stay clear of bash's own
# special variables — a plain assignment to one of those fails under `set -e` with
# no output at all.
BC2_NAME=$(echo "$SAPIA_LINE" | cut -f1)
BC2_SETTINGS=$(echo "$SAPIA_LINE" | cut -f2)
BC2_FLAGS=$(echo "$SAPIA_LINE" | cut -f3)

# Point bindcraft's AlphaFold-parameter cache at the mounted Volume, at RUNTIME.
#
# BC2_CACHE_DIR is set by the Modal image; it is unset elsewhere (e.g. vib),
# where the site activation script decides the cache instead, so this is a
# no-op there. It is deliberately NOT called XDG_CACHE_HOME in the image:
# setting that name at image-build time makes the builder's own `uv` write
# $XDG_CACHE_HOME/uv into the image, which then blocks the Volume from
# mounting and kills the container before it starts -- with no .out/.err/.exit
# written at all. Converting it here keeps the variable out of the build.
if [ -n "${BC2_CACHE_DIR:-}" ]; then
    export XDG_CACHE_HOME="$BC2_CACHE_DIR"
    echo "  cache:    XDG_CACHE_HOME=$XDG_CACHE_HOME (from BC2_CACHE_DIR)"
fi

BC2_BIN=${BINDCRAFT_BIN:-bindcraft}

echo "[$(date +%T)] task $SAPIA_TASK_ID: bindcraft campaign '$BC2_NAME'"
echo "  settings: $BC2_SETTINGS"
echo "  flags:    ${BC2_FLAGS:-(none)}"

# BC2_FLAGS is intentionally unquoted: it is a space-separated list of CLI tokens
# (--modality X, --humanize, --set k=v) that must each become its own argv entry.
"$BC2_BIN" design "$BC2_SETTINGS" $BC2_FLAGS
