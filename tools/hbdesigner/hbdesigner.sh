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

# Locals cut from the manifest line. Do NOT rename any of these to a bash special
# variable (GROUPS, UID, PIPESTATUS, SECONDS...): the assignment fails under
# `set -e` and the task dies with EMPTY .out and .err.
NAME=$(printf '%s' "$SAPIA_LINE" | cut -f1)
SRC=$(printf '%s' "$SAPIA_LINE" | cut -f2)
GRAFT_CHAINS=$(printf '%s' "$SAPIA_LINE" | cut -f3)
OMIT_CHAINS=$(printf '%s' "$SAPIA_LINE" | cut -f4)
# Kept last and always non-empty, so the manifest line never ends on an empty field.
HBD_ARGS=$(printf '%s' "$SAPIA_LINE" | cut -f5)

# GRAFT_CHAINS and OMIT_CHAINS are legitimately empty (monomer, two-sided design),
# so they are not tested; the other three are never empty in a well-formed line.
# A short line means the builder and this script disagree about the manifest
# format -- fail loudly rather than run with fields silently shifted by one.
if [ -z "$NAME" ] || [ -z "$SRC" ] || [ -z "$HBD_ARGS" ]; then
    echo "hbdesigner: manifest line has fewer than 5 fields: $SAPIA_LINE" >&2
    exit 2
fi

# NAME comes off the manifest line and is about to be used in `rm -rf`. The
# builder already refuses an unusable name at submit, but this script must not
# depend on that: validate it here too, before anything deletes. A name is one
# path component of letters, digits, '.', '-' and '_', starting with a letter or
# digit -- so it can hold no '/', no whitespace and cannot be '..'.
case "$NAME" in
    *[!A-Za-z0-9._-]* | [!A-Za-z0-9]* | *..* )
        echo "hbdesigner: refusing unusable design name '$NAME': a name must start with a letter or digit and hold only letters, digits, '.', '-' and '_' (it becomes a directory this task deletes and recreates)" >&2
        exit 2
        ;;
esac

# Clear this design's previous output BEFORE anything else -- before activation,
# before the interpreter is even chosen. HBDesigner writes NO stats CSV and NO
# rank PDB when nothing survives scoring (the zero-output success), so on a
# --force rerun that finds nothing the worker would otherwise read the PREVIOUS
# run's CSV and rank PDBs and emit full OK rows: real-looking metrics from a
# parameter set that is no longer the one in the manifest. Removing them first
# makes "this run found nothing" look like what it is. It sits above
# sapia_activate for the same reason: if activation fails the task dies here, and
# a surviving <design>.tsv from the previous run would be collected as this run's
# result. The only code before this point is the prelude and the two guards.
# Both paths use the ${VAR:?} form: an empty OUT_DIR or NAME aborts the task
# instead of expanding to a path that means something else. Only this design's
# own files are touched.
rm -f "${OUT_DIR:?}/${NAME:?}.tsv"
rm -rf "${OUT_DIR:?}/${NAME:?}"

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
# then exits with the real code so the worker sees it too.
set +e
"$HBD_BIN" \
    --pdb "$STAGED" \
    --out_dir "$DESIGN_DIR" \
    $HBD_ARGS
RUN_RC=$?
set -e

# With --guide_res the installed HBDesigner appends a virtual guide atom to every
# rank PDB as a HETATM (atom V1, residue ORI, element V, in a new trailing chain).
# It is a marker, not structure, and upstream's own graft_seq.py asserts the target
# and reference have the same residue count -- which the extra entry can break. So
# strip it first, keeping the untouched file under <design>/guide_raw/. Non-fatal:
# a failure here leaves `guide_atoms_removed` NA (not 0) and the row still collects.
if [ "$RUN_RC" -eq 0 ]; then
    if ! "$PY" "${SAPIA_TOOL_DIR:?}/hbdesigner_worker.py" \
            --mode strip-guide \
            --name "$NAME" \
            --design-dir "$DESIGN_DIR"; then
        echo "[$(date +%T)] WARNING: guide-atom strip failed for $NAME" >&2
    fi
fi

# One-sided interface design returns the omitted chain as POLY-GLYCINE. Upstream's
# graft_seq.py puts its sequence and sidechains back, copying only onto positions
# that came back as glycine (so designed/anchor residues survive). Without this
# the next tool is handed a target chain that is not there.
#
# Like the run itself, this must NOT be fatal: graft_seq.py asserts
# `p_target.n_res == p_ref.n_res`, so a stray atom (a guide atom the strip above
# failed to remove) aborts it. Under `set -e` that killed the task before the
# worker ran, leaving no TSV at all -- and collect treats a missing TSV as "no
# evidence" and skips the design silently, which is the one outcome this script
# is written to avoid. So the loop runs under `set +e`, the first failing rank
# sets GRAFT_RC, and the worker turns that into an `error: graft failed` row.
GRAFT_RC=0
if [ "$RUN_RC" -eq 0 ] && [ -n "$GRAFT_CHAINS" ]; then
    mkdir -p "$GRAFT_DIR"
    shopt -s nullglob
    set +e
    for rank_pdb in "$DESIGN_DIR"/${NAME}_HBDes_rank_*.pdb; do
        echo "[$(date +%T)] grafting chain(s) $GRAFT_CHAINS onto $(basename "$rank_pdb")"
        "$PY" -m hbdesigner.scripts.graft_seq \
            --target_pdb "$rank_pdb" \
            --ref_pdb "$STAGED" \
            --out_pdb "$GRAFT_DIR/$(basename "$rank_pdb")" \
            --graft_chains "$GRAFT_CHAINS"
        RC=$?
        if [ "$RC" -ne 0 ]; then
            echo "[$(date +%T)] ERROR: graft_seq exited $RC on $(basename "$rank_pdb")" >&2
            GRAFT_RC=$RC
        fi
    done
    set -e
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
        --omit-chains "$OMIT_CHAINS" \
        --run-rc "$RUN_RC" \
        --graft-rc "$GRAFT_RC" \
        --result-tsv "$RESULT_TSV"; then
    printf 'name\trank\tstatus\tpath\n' >"$RESULT_TSV"
    printf '%s\t0\terror: worker crashed\t\n' "$NAME" >>"$RESULT_TSV"
fi

# The TSV is written either way; the exit code is what makes the task visibly
# failed in the ledger, so a graft failure is both a row AND a red task.
if [ "$RUN_RC" -ne 0 ]; then
    exit "$RUN_RC"
fi
exit "$GRAFT_RC"
