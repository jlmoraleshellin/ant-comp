#!/bin/bash
#SBATCH --job-name=atomiumph
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Slurm-only activation: the script must export ATOMIUMPH and, if needed, ATOMIUMPH_PYTHON.
sapia_activate SAPIA_ACTIVATE_ATOMIUMPH

ATOMIUMPH_PYTHON=${ATOMIUMPH_PYTHON:-python}
: "${ATOMIUMPH:?ATOMIUMPH must point at the AtomiUM-pH checkout (the image sets it on Modal; the activation script elsewhere)}"

# PYG_HOME (a per-task PyG cache) is exported by the vib ACTIVATION SCRIPT, not here.

# One task = one sub-manifest; row: grp_dir binder_chains target_chains lam extra
run_one_group() {
    # Lowercase locals: an uppercase bash special (GROUPS, UID...) kills the task under `set -e`.
    local grp_dir=$1 binder_chains=$2 target_chains=$3 lam=$4 extra=$5

    local inputs="$grp_dir/inputs"
    local out="$grp_dir/out"
    mkdir -p "$out"

    echo "Task ${SAPIA_TASK_ID}: $(basename "$grp_dir") ($(ls "$inputs" | wc -l) structure(s))"

    # Provenance: the image's COMMIT becomes a table column (collect walks up to find it).
    if [ -f "${ATOMIUMPH:?}/COMMIT" ]; then
        cp "${ATOMIUMPH:?}/COMMIT" "$grp_dir/image_commit.txt"
        echo "image commit: $(cat "${ATOMIUMPH:?}/COMMIT")"
    fi
    # set -f from here until design_ph returns: unquoted $extra would otherwise glob a * ? [ in --set.
    set -f
    local noise
    # `|| true`: an unmatched grep must not kill the task (set -e) over a log line.
    noise=$(printf '%s\n' $extra | grep -A1 -x -- '--checkpoint_model_noise' \
        | tail -1) || true
    # Log the checkpoint design_ph opens: ft_weights_for_pH/ft_model_epoch10_<noise>.pt (30-class pH head).
    echo "--- pH checkpoints in $ATOMIUMPH/ft_weights_for_pH ---"
    ls -la "${ATOMIUMPH:?}"/ft_weights_for_pH 2>/dev/null \
        || echo "(no ft_weights_for_pH dir — design_ph.py will fail to load weights)"
    if [ -n "$noise" ]; then
        local ckpt="${ATOMIUMPH:?}/ft_weights_for_pH/ft_model_epoch10_${noise}.pt"
        if [ -f "$ckpt" ]; then
            echo "checkpoint for --checkpoint_model_noise ${noise}: ${ckpt} ($(wc -c < "$ckpt") bytes)"
        else
            echo "checkpoint for --checkpoint_model_noise ${noise}: ${ckpt} IS MISSING"
        fi
    else
        echo "(could not read --checkpoint_model_noise out of the argv tail)"
    fi
    echo "-----------------------------------------------------"

    # Chains stay quoted (one argument each); $lam and $extra are unquoted, whitespace-free by construction.
    "$ATOMIUMPH_PYTHON" "${ATOMIUMPH:?}/design_ph.py" \
        --pdb_dir "$inputs" \
        --pattern '*.pdb' \
        --out_folder "$out" \
        --binder_chains "$binder_chains" \
        --target_chains "$target_chains" \
        --lam $lam \
        $extra < /dev/null  # stdin closed: a read would swallow the remaining sub-manifest rows
    set +f
}

# Upstream's --help (dumped at image build) is the only authority on forwarded flags; never fatal.
echo "===== design_ph.py --help (image build, pinned commit) ====="
if [ -s "${ATOMIUMPH:?}/DESIGN_PH_HELP.txt" ]; then
    cat "${ATOMIUMPH:?}/DESIGN_PH_HELP.txt" || true
else
    echo "(no help dump in the image — it could not be rendered at build time)"
fi
echo "===== end design_ph.py --help ====="

SUBMANIFEST=$(printf '%s' "$SAPIA_LINE" | cut -f1)

while IFS= read -r row || [ -n "$row" ]; do
    [ -z "$row" ] && continue
    run_one_group \
        "$(printf '%s' "$row" | cut -f1)" \
        "$(printf '%s' "$row" | cut -f2)" \
        "$(printf '%s' "$row" | cut -f3)" \
        "$(printf '%s' "$row" | cut -f4)" \
        "$(printf '%s' "$row" | cut -f5)"
done < "$SUBMANIFEST"
