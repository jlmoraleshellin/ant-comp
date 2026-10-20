#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --job-name=mkcomplex

# One design's complex sequence into <name>.tsv: prepend x REPEAT, design x DESIGN_COPIES, append x REPEAT.

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Activation is OPTIONAL here (bash + awk only): sourced only when SAPIA_ACTIVATE_MKCOMPLEX is set.
if [[ -n "${SAPIA_ACTIVATE_MKCOMPLEX:-}" ]]; then
    sapia_activate SAPIA_ACTIVATE_MKCOMPLEX
fi

NAME=$(echo "$SAPIA_LINE" | cut -f1)
SEQ=$(echo "$SAPIA_LINE" | cut -f2)
PREPEND_SEQS=$(echo "$SAPIA_LINE" | cut -f3)
PREPEND_FASTA=$(echo "$SAPIA_LINE" | cut -f4)
APPEND_SEQS=$(echo "$SAPIA_LINE" | cut -f5)
APPEND_FASTA=$(echo "$SAPIA_LINE" | cut -f6)
SEP=$(echo "$SAPIA_LINE" | cut -f7)
REPEAT=$(echo "$SAPIA_LINE" | cut -f8)
DESIGN_COPIES=$(echo "$SAPIA_LINE" | cut -f9)

RESULT_TSV="$OUT_DIR/${NAME}.tsv"

# The 20 standard uppercase one-letter codes; a stray header or masked residue would corrupt predictions.
AA_RE='^[ACDEFGHIKLMNPQRSTVWY]+$'

write_result() { # status sequence n_chains chain_lens design_chain_index total_len design_copies
    printf 'name\tstatus\tsequence\tn_chains\tchain_lens\tdesign_chain_index\ttotal_len\tdesign_copies\n' \
        >"$RESULT_TSV"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$NAME" "$1" "$2" "$3" "$4" "$5" "$6" "$7" >>"$RESULT_TSV"
}

# Errors are DATA (status 'error: ...', exit 0) so one bad row never fails the array; values stay empty -> NA.
fail() {
    echo "$NAME: $1" >&2
    write_result "error: $1" "" "" "" "" "" ""
    exit 0
}

# Counts are validated at submit; this only catches a malformed manifest line (error row, not a wrong count).
require_count() { # flag value
    if [[ ! "$2" =~ ^[0-9]+$ ]] || ((10#$2 < 1)); then
        fail "$1 must be a positive integer, got '${2:0:20}' (malformed manifest line)"
    fi
}

# One sequence per record, whitespace stripped; a header-only record yields an EMPTY line so add_fixed rejects it.
fasta_seqs() {
    awk '
        /^>/ { if (seen) print seq; seq = ""; seen = 1; next }
        { gsub(/[ \t\r]/, "", $0); seq = seq $0 }
        END { if (seen) print seq }
    ' "$1"
}

CHAINS=()
N_FIXED=0

add_fixed() { # one fixed chain -> REPEAT consecutive copies
    local s=$1
    s="${s//[[:space:]]/}"
    # An empty element of a non-empty list ('T,,T', a header-only record) is bad input, not skipped.
    if [[ -z "$s" ]]; then
        fail "empty chain in the fixed-chain list (an interior ',,' or a FASTA record with no sequence)"
    fi
    if [[ ! "$s" =~ $AA_RE ]]; then
        fail "invalid residue in fixed chain '${s:0:40}' (only the 20 standard uppercase one-letter codes)"
    fi
    local i
    for ((i = 0; i < REPEAT; i++)); do
        CHAINS+=("$s")
    done
    N_FIXED=$((N_FIXED + 1))
}

add_side() { # csv fasta -> the side's fixed chains, in order
    local csv=$1 fasta=$2 chain
    local raw=()
    if [[ -n "$csv" ]]; then
        mapfile -t raw < <(printf '%s' "$csv" | tr ',' '\n')
    elif [[ -n "$fasta" ]]; then
        [[ -f "$fasta" ]] || fail "missing fasta ($fasta)"
        mapfile -t raw < <(fasta_seqs "$fasta")
        if ((${#raw[@]} == 0)); then
            fail "no records in fasta ($fasta)"
        fi
    else
        return 0
    fi
    for chain in ${raw[@]+"${raw[@]}"}; do
        add_fixed "$chain"
    done
}

echo "[$(date +%T)] task $SAPIA_TASK_ID: mkcomplex for $NAME"

require_count "--repeat" "$REPEAT"
require_count "--design-copies" "$DESIGN_COPIES"

add_side "$PREPEND_SEQS" "$PREPEND_FASTA"
# The index of the FIRST design copy, whatever DESIGN_COPIES is.
DESIGN_INDEX=${#CHAINS[@]}

SEQ="${SEQ//[[:space:]]/}"
case "$SEQ" in
    "" | NA | NaN | nan | None | none | null)
        fail "empty input sequence" ;;
esac
if [[ "$SEQ" == *"$SEP"* ]]; then
    fail "input sequence already contains the separator '$SEP' (already a multi-chain string)"
fi
# The design's own chain gets the same alphabet check: a masked/'X' residue would be copied into every prediction.
if [[ ! "$SEQ" =~ $AA_RE ]]; then
    fail "invalid residue in design sequence '${SEQ:0:40}' (only the 20 standard uppercase one-letter codes)"
fi
for ((i = 0; i < DESIGN_COPIES; i++)); do
    CHAINS+=("$SEQ")
done

add_side "$APPEND_SEQS" "$APPEND_FASTA"

if ((N_FIXED == 0)); then
    fail "no fixed chains to add (empty --prepend-*/--append-*)"
fi

# Chains map positionally onto A-Z downstream (select_chains), which truncates silently past 26.
if ((${#CHAINS[@]} > 26)); then
    fail "${#CHAINS[@]} chains exceeds the 26-chain ceiling: chains are mapped positionally onto A-Z and anything past chain 26 is silently dropped downstream"
fi

ASSEMBLED=""
CHAIN_LENS=""
TOTAL_LEN=0
for chain in "${CHAINS[@]}"; do
    ASSEMBLED="${ASSEMBLED:+${ASSEMBLED}${SEP}}${chain}"
    CHAIN_LENS="${CHAIN_LENS:+${CHAIN_LENS},}${#chain}"
    TOTAL_LEN=$((TOTAL_LEN + ${#chain}))
done

write_result "OK" "$ASSEMBLED" "${#CHAINS[@]}" "$CHAIN_LENS" "$DESIGN_INDEX" "$TOTAL_LEN" "$DESIGN_COPIES"
echo "[$(date +%T)] $NAME: ${#CHAINS[@]} chains, $TOTAL_LEN residues, ${DESIGN_COPIES} design copies starting at index $DESIGN_INDEX"
