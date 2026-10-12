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

# Site-specific activation (required off Modal) — a no-op under Modal, where the
# image supplies $ATOMIUMPH. This tool is Modal-only today: there is no
# SAPIA_ACTIVATE_ATOMIUMPH anywhere, so a slurm task would stop here with a clear
# message rather than running against the wrong environment.
sapia_activate SAPIA_ACTIVATE_ATOMIUMPH

ATOMIUMPH_PYTHON=${ATOMIUMPH_PYTHON:-python}

# One task = one sub-manifest listing the subgroups packed onto it. Each subgroup is a
# grp_<g>/ dir whose inputs/ holds symlinked complexes that share parameters, so it runs
# as ONE batched design_ph.py call in --pdb_dir mode (checkpoint loaded once). Subgroups
# run sequentially on the task's single GPU. Sub-manifest rows are tab-separated:
#   grp_dir  binder_chains  target_chains  lam  extra
#
# Locals are lowercase on purpose: an uppercase name colliding with a bash special
# (GROUPS, UID, PIPESTATUS…) fails its assignment under `set -e` and kills the task
# before it prints anything. See the trap note in authoring-a-tool.
run_one_group() {
    local grp_dir=$1 binder_chains=$2 target_chains=$3 lam=$4 extra=$5

    local inputs="$grp_dir/inputs"
    local out="$grp_dir/out"
    mkdir -p "$out"

    echo "Task ${SAPIA_TASK_ID}: $(basename "$grp_dir") ($(ls "$inputs" | wc -l) structure(s))"

    # Provenance. The image's COMMIT becomes a table column (collect walks up from each
    # output file to find it), so two runs from different pins of a moving private
    # branch are distinguishable in the table rather than only in a build log.
    if [ -f "${ATOMIUMPH:?}/COMMIT" ]; then
        cp "${ATOMIUMPH:?}/COMMIT" "$grp_dir/image_commit.txt"
        echo "image commit: $(cat "${ATOMIUMPH:?}/COMMIT")"
    fi
    # WHICH WEIGHTS. READ FROM SOURCE at 7ffb756 (design_ph.py:2041-2042): the
    # checkpoint is ft_weights_for_pH/ft_model_epoch10_<noise>.pt -- the fine-tuned
    # 30-class (protonation-aware) head. model_weights/model_noised_=_*.pt is the
    # 21-class base that only atomium.py reads; earlier runs listed THAT directory,
    # which said nothing about what design_ph opened.
    #
    # So print the resolved path and its byte size. The noise level is in $extra, as
    # `--checkpoint_model_noise <level>`, which is also where the builder puts it.
    #
    # `set -f` STARTS HERE, not at the design_ph call: the line below expands $extra
    # unquoted to word-split it, and pathname expansion rides along with word
    # splitting — so a `*`, `?` or `[` inside a --set value would be replaced by
    # whatever files happen to match in the cwd, HERE, before the call the guard was
    # written for. It stays off until after design_ph returns.
    #
    # `|| true`: the task runs under `set -euo pipefail`, so an unmatched grep in this
    # pipeline would kill the task before design_ph ever starts — over a log line.
    set -f
    local noise
    noise=$(printf '%s\n' $extra | grep -A1 -x -- '--checkpoint_model_noise' \
        | tail -1) || true
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

    # Quoting, field by field — this is where a wrong choice silently changes the run:
    #   $binder_chains / $target_chains are ONE argument each ("B C D" is a single
    #     space-separated string upstream), so they are quoted.
    #   $lam is several argv entries (--lam takes a list of floats), so it is expanded
    #     UNQUOTED. Values are formatted whitespace-free by the builder.
    #   $extra (every other typed flag + user --set tokens) is expanded unquoted too
    #     and is whitespace-free by construction.
    #
    # `set -f` is already on (see the checkpoint block above) and covers every
    # unquoted expansion in this function, including this one. Restored right after.
    #
    # stdin from /dev/null: this function is called from a `while read` loop fed by the
    # sub-manifest on stdin. If design_ph ever reads stdin it would swallow the
    # remaining rows, those subgroups would never run, and the task would still exit 0.
    "$ATOMIUMPH_PYTHON" "${ATOMIUMPH:?}/design_ph.py" \
        --pdb_dir "$inputs" \
        --pattern '*.pdb' \
        --out_folder "$out" \
        --binder_chains "$binder_chains" \
        --target_chains "$target_chains" \
        --lam $lam \
        $extra < /dev/null
    set +f
}

# Upstream's own --help, dumped into the image at build time. Echoed ONCE per task,
# before any structure, because it is the only authority on what the flags this tool
# forwards actually do — several (notably --protect_titratable) are passed on the
# user's word alone and are documented as unverified. Read it on the first run of a
# new pin and correct the skill. Never fatal: a missing or empty dump is a note, not
# a failure, and `|| true` keeps `set -e` out of it.
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
