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

# Site activation (required off Modal; a no-op there, the image supplies $ATOMIUM).
sapia_activate SAPIA_ACTIVATE_ATOMIUM

ATOMIUM_PYTHON=${ATOMIUM_PYTHON:-python}
: "${ATOMIUM:?ATOMIUM must point at the AtomiUM checkout (the image sets it on Modal; the activation script elsewhere)}"

# PYG_HOME (a per-task PyG cache) is exported by the vib ACTIVATION SCRIPT

# ProteinMPNN helpers: in the checkout on Modal; on vib the activation script exports PROTEINMPNN_HELPERS.
HELPERS=${PROTEINMPNN_HELPERS:-"$ATOMIUM/proteinMPNN_helper_scripts"}
# All four are checked up front: a missing one would otherwise kill the task mid-loop
for h in parse_multiple_chains assign_fixed_chains make_fixed_positions_dict make_tied_positions_dict; do
    if [ ! -f "$HELPERS/$h.py" ]; then
        echo "atomium: no $h.py in HELPERS='$HELPERS' (PROTEINMPNN_HELPERS='${PROTEINMPNN_HELPERS:-}', ATOMIUM='$ATOMIUM'). PROTEINMPNN_HELPERS must be ProteinMPNN's helper_scripts directory itself, not its install root." >&2
        exit 1
    fi
done
# Provenance: the .out is the only place the helper source is recorded.
echo "atomium: HELPERS=$HELPERS"

# One task = one sub-manifest; row: grp_dir chains fixed_positions tie_mode tied_positions sampling_temp extra
run_one_group() {
    # Lowercase locals: an uppercase bash special (GROUPS, UID...) kills the task under `set -e`.
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

    # Side-input jsonls, appended only when produced; every per-chain helper shares $chains.
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

    # Tying: homo => --homooligomer 1 (all chains); explicit => per-chain lists.
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

    # $extra is unquoted (one argv entry per token); $sampling_temp is quoted (several temperatures, ONE argument).
    "$ATOMIUM_PYTHON" "$ATOMIUM/atomium.py" \
        --jsonl_path "$parsed" \
        --out_folder "$grp_dir" \
        --sampling_temp "$sampling_temp" \
        $extra \
        ${side[@]+"${side[@]}"}
}

SUBMANIFEST=$(printf '%s' "$SAPIA_LINE" | cut -f1)

# cut, not read: IFS=$'\t' read collapses the empty tab fields (chains, positions, tie_mode).
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
