#!/bin/bash
#SBATCH --job-name=hbdesigner
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# One array task = one input backbone. Runs HBDesigner on it, optionally grafts
# omitted chains back from the tool's own input PDB, then runs the worker, which
# writes the per-design TSV collect_hbdesigner.py reads.
#
# Layout written per design:
#   $OUT_DIR/<name>/<name>.pdb                      staged input (= graft reference)
#   $OUT_DIR/<name>/<name>_HBDes_rank_<i>.pdb       one per kept design
#   $OUT_DIR/<name>/<name>_HBDes_stats.csv          upstream's metrics table
#   $OUT_DIR/<name>/grafted/<name>_HBDes_rank_<i>.pdb   only when grafting ran
#   $OUT_DIR/<name>.tsv                             one row per kept design
#
# The input is COPIED in as <name>.pdb on purpose: upstream names every output
# after the input file's stem, so staging it under the design name is what makes
# the output filenames predictable (and collision-free between designs).

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required off Modal). Under Modal this is a no-op and
# the image supplies run_hbdesigner, PyRosetta, torch and the baked model weights.
sapia_activate SAPIA_ACTIVATE_HBDESIGNER

# wandb is a hard import dependency (inference imports the training module, which
# imports wandb at module level). Keep it from trying to reach the network from a
# task container; the image sets these too, this is the belt to that braces.
export WANDB_MODE=${WANDB_MODE:-disabled}
export WANDB_SILENT=${WANDB_SILENT:-true}

PY=${HBDESIGNER_PYTHON:-python}
# Upstream's console script (pyproject [project.scripts]). Overridable for a site
# that installs it under another name or wraps it.
HBD_BIN=${HBDESIGNER_BIN:-run_hbdesigner}

# Locals cut from the manifest line. Do NOT rename any of these to a bash special
# variable (GROUPS, UID, PIPESTATUS, SECONDS...): the assignment fails under
# `set -e` and the task dies with EMPTY .out and .err.
NAME=$(printf '%s' "$SAPIA_LINE" | cut -f1)
SRC=$(printf '%s' "$SAPIA_LINE" | cut -f2)
GRAFT_CHAINS=$(printf '%s' "$SAPIA_LINE" | cut -f3)
HBD_ARGS=$(printf '%s' "$SAPIA_LINE" | cut -f4)

DESIGN_DIR="$OUT_DIR/$NAME"
STAGED="$DESIGN_DIR/${NAME}.pdb"
GRAFT_DIR="$DESIGN_DIR/grafted"
RESULT_TSV="$OUT_DIR/${NAME}.tsv"

mkdir -p "$DESIGN_DIR"
cp -f "$SRC" "$STAGED"

echo "[$(date +%T)] task $SAPIA_TASK_ID: hbdesigner $NAME ($HBD_ARGS)"

# $HBD_ARGS is unquoted so each space-separated token becomes its own argv entry;
# the builder refuses any token containing whitespace for exactly that reason.
# A non-zero exit is NOT fatal here: the worker still runs, so the failure is
# recorded in the table as a status instead of vanishing into a log, and the task
# then exits with the real code so the orchestrator sees it too.
set +e
"$HBD_BIN" \
    --pdb "$STAGED" \
    --out_dir "$DESIGN_DIR" \
    $HBD_ARGS
RUN_RC=$?
set -e

# One-sided interface design returns the omitted chain as POLY-GLYCINE. Upstream's
# graft_seq.py puts its sequence and sidechains back, copying only onto positions
# that came back as glycine (so designed/anchor residues survive). Without this
# the next tool is handed a target chain that is not there.
if [ "$RUN_RC" -eq 0 ] && [ -n "$GRAFT_CHAINS" ]; then
    mkdir -p "$GRAFT_DIR"
    shopt -s nullglob
    for rank_pdb in "$DESIGN_DIR"/${NAME}_HBDes_rank_*.pdb; do
        echo "[$(date +%T)] grafting chain(s) $GRAFT_CHAINS onto $(basename "$rank_pdb")"
        "$PY" -m hbdesigner.scripts.graft_seq \
            --target_pdb "$rank_pdb" \
            --ref_pdb "$STAGED" \
            --out_pdb "$GRAFT_DIR/$(basename "$rank_pdb")" \
            --graft_chains "$GRAFT_CHAINS"
    done
    shopt -u nullglob
fi

# The worker lives next to this .sh; SAPIA_TOOL_DIR points there regardless of cwd.
# It always runs -- a design that produced nothing must still be recorded.
if ! "$PY" "${SAPIA_TOOL_DIR:?}/hbdesigner_worker.py" \
        --name "$NAME" \
        --design-dir "$DESIGN_DIR" \
        --ref-pdb "$STAGED" \
        --graft-chains "$GRAFT_CHAINS" \
        --graft-dir "$GRAFT_DIR" \
        --run-rc "$RUN_RC" \
        --result-tsv "$RESULT_TSV"; then
    printf 'name\trank\tstatus\tpath\n' >"$RESULT_TSV"
    printf '%s\t0\terror: worker crashed\t\n' "$NAME" >>"$RESULT_TSV"
fi

exit "$RUN_RC"
