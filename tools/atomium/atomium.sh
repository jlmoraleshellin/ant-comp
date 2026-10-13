#!/bin/bash
#SBATCH --job-name=atomium
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required off Modal) — see docs/configuration.md.
# Under Modal this is a no-op and the image supplies $ATOMIUM.
sapia_activate SAPIA_ACTIVATE_ATOMIUM

ATOMIUM_PYTHON=${ATOMIUM_PYTHON:-python}
: "${ATOMIUM:?ATOMIUM must point at the AtomiUM checkout (the image sets it on Modal; the activation script elsewhere)}"

# ProteinMPNN's helper scripts. On Modal they ship inside the AtomiUM checkout. On vib that
# checkout has no proteinMPNN_helper_scripts/, so the vib ACTIVATION SCRIPT exports PROTEINMPNN_HELPERS,
# the cluster's ProteinMPNN `helper_scripts` directory itself. 
HELPERS=${PROTEINMPNN_HELPERS:-"$ATOMIUM/proteinMPNN_helper_scripts"}
# All four are checked up front: a missing one would otherwise kill the task mid-loop
for h in parse_multiple_chains assign_fixed_chains make_fixed_positions_dict make_tied_positions_dict; do
    if [ ! -f "$HELPERS/$h.py" ]; then
        echo "atomium: no $h.py in HELPERS='$HELPERS' (PROTEINMPNN_HELPERS='${PROTEINMPNN_HELPERS:-}', ATOMIUM='$ATOMIUM'). PROTEINMPNN_HELPERS must be ProteinMPNN's helper_scripts directory itself, not its install root." >&2
        exit 1
    fi
done
# Provenance: which helper source this run used (the .out is the only place it is recorded).
echo "atomium: HELPERS=$HELPERS"

# One task = one sub-manifest listing the subgroups packed onto it. Each subgroup is
# a grp_<g>/ dir whose inputs/ holds symlinked PDBs that share params, so it runs as
# ONE batched atomium.py call (model loaded once). Subgroups run sequentially on the
# task's single GPU. Sub-manifest rows are tab-separated:
#   grp_dir  chains  fixed_positions  tie_mode  tied_positions  sampling_temp  extra
# chains, fixed_positions, tie_mode and tied_positions may be empty AND are adjacent,
# so fields are pulled with `cut` (IFS=$'\t' read collapses empty tab fields).
#
# Locals are lowercase on purpose: an uppercase name that collides with a bash
# special variable (GROUPS, UID, PIPESTATUS…) fails its assignment under `set -e` and
# kills the task before it prints anything. See the trap note in authoring-a-tool.
run_one_group() {
    local grp_dir=$1 chains=$2 fixed_positions=$3 tie_mode=$4
    local tied_positions=$5 sampling_temp=$6 extra=$7

    local inputs="$grp_dir/inputs"
    local parsed="$grp_dir/parsed_pdbs.jsonl"
    local assigned="$grp_dir/assigned_pdbs.jsonl"
    local fixed="$grp_dir/fixed_pdbs.jsonl"
    local tied="$grp_dir/tied_pdbs.jsonl"

    echo "Task ${SAPIA_TASK_ID}: $(basename "$grp_dir") ($(ls "$inputs" | wc -l) design(s))"

    "$ATOMIUM_PYTHON" "$HELPERS/parse_multiple_chains.py" \
        --input_path="$inputs" --output_path="$parsed"

    # Optional side-input jsonls, each appended to the run only when produced. Every
    # per-chain helper shares the one $chains list as its --chain_list.
    local side=()

    if [ -n "$chains" ]; then
        "$ATOMIUM_PYTHON" "$HELPERS/assign_fixed_chains.py" \
            --input_path="$parsed" --output_path="$assigned" \
            --chain_list "$chains"
        side+=(--chain_id_jsonl "$assigned")
    fi

    if [ -n "$fixed_positions" ]; then
        "$ATOMIUM_PYTHON" "$HELPERS/make_fixed_positions_dict.py" \
            --input_path="$parsed" --output_path="$fixed" \
            --chain_list "$chains" --position_list "$fixed_positions"
        side+=(--fixed_positions_jsonl "$fixed")
    fi

    # Tying: homo => --homooligomer 1 (all chains, auto-detected); explicit => per-chain
    # --chain_list/--position_list. --homooligomer defaults to 0, so it is never passed.
    if [ "$tie_mode" = "homo" ]; then
        "$ATOMIUM_PYTHON" "$HELPERS/make_tied_positions_dict.py" \
            --input_path="$parsed" --output_path="$tied" --homooligomer 1
        side+=(--tied_positions_jsonl "$tied")
    elif [ "$tie_mode" = "explicit" ]; then
        "$ATOMIUM_PYTHON" "$HELPERS/make_tied_positions_dict.py" \
            --input_path="$parsed" --output_path="$tied" \
            --chain_list "$chains" --position_list "$tied_positions"
        side+=(--tied_positions_jsonl "$tied")
    fi

    # $extra (typed run flags + optional --bias_AA_jsonl + user --set tokens) is
    # unquoted so each space-separated token becomes its own argv entry; it is built
    # whitespace-free for exactly that reason. $sampling_temp is quoted instead --
    # it may hold several temperatures ("0.1 0.2") that must stay ONE argument.
    # Outputs land in $grp_dir/seqs/<design>.fa (one fasta per staged input).
    "$ATOMIUM_PYTHON" "$ATOMIUM/atomium.py" \
        --jsonl_path "$parsed" \
        --out_folder "$grp_dir" \
        --sampling_temp "$sampling_temp" \
        $extra \
        ${side[@]+"${side[@]}"}
}

SUBMANIFEST=$(printf '%s' "$SAPIA_LINE" | cut -f1)

while IFS= read -r row || [ -n "$row" ]; do
    [ -z "$row" ] && continue
    run_one_group \
        "$(printf '%s' "$row" | cut -f1)" \
        "$(printf '%s' "$row" | cut -f2)" \
        "$(printf '%s' "$row" | cut -f3)" \
        "$(printf '%s' "$row" | cut -f4)" \
        "$(printf '%s' "$row" | cut -f5)" \
        "$(printf '%s' "$row" | cut -f6)" \
        "$(printf '%s' "$row" | cut -f7)"
done < "$SUBMANIFEST"
