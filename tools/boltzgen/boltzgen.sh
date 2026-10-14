#!/bin/bash
#SBATCH --job-name=boltzgen
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE. Our manifest
# line is one BoltzGen run: its design-group name, the spec YAML the submitter wrote,
# the --output dir, the run-wide flags, the PROTOCOL and the STEP LIST.
#
# 64G of host RAM is not a guess: the official slurm example asks for it
# (slurm-example/run_job_array.slurm) and upstream issue #208 reports the `analysis`
# step OOM-ing at 16G with only 100 designs.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation (required off Modal) -- must put `boltzgen` on PATH.
# A no-op under --executor modal, where modal_image.py provides the environment.
sapia_activate SAPIA_ACTIVATE_BOLTZGEN

# Manifest columns (tab-separated). Names are prefixed to stay clear of bash's own
# special variables -- a plain assignment to one of those fails under `set -e` with
# no output at all.
BG_NAME=$(echo "$SAPIA_LINE" | cut -f1)
BG_SPEC=$(echo "$SAPIA_LINE" | cut -f2)
BG_OUT=$(echo "$SAPIA_LINE" | cut -f3)
BG_FLAGS=$(echo "$SAPIA_LINE" | cut -f4)
BG_PROTOCOL=$(echo "$SAPIA_LINE" | cut -f5)
BG_STEPS=$(echo "$SAPIA_LINE" | cut -f6)

# Fields 5 and 6 used to be the literals `nanobody-anything` and
# `design inverse_folding folding analysis` in this script. They moved into the
# manifest when --protocol was added, because the two protocols need DIFFERENT step
# lists (protein-anything also needs design_folding; see run_boltzgen.py's docstring
# for why omitting it empties the metrics CSV rather than crashing). They are still
# the TOOL'S CONTRACT, not a user flag: run_boltzgen.py derives both from --protocol
# and BLOCKED_EXTRA_FLAGS refuses --protocol/--steps inside --extra-flags.
#
# Refuse an empty one rather than running `boltzgen run --protocol --steps ...`, which
# argparse would read as a protocol literally named "--steps". An old manifest written
# by the pre-protocol version of this tool has four fields, so this is exactly the
# case to catch loudly.
: "${BG_PROTOCOL:?manifest field 5 (protocol) is empty -- re-run \`sapia run boltzgen\` to regenerate the manifest}"
: "${BG_STEPS:?manifest field 6 (steps) is empty -- re-run \`sapia run boltzgen\` to regenerate the manifest}"

BG_BIN=${BOLTZGEN_BIN:-boltzgen}

# Where the ~6.4 GB of checkpoints live. Set by modal_image.py to the mounted Volume;
# unset elsewhere, where the site activation script decides (or BoltzGen's ~/.cache
# default applies). Converted to the flag at RUNTIME, never baked into the image as an
# env var the image builder could key off -- see the modal_image.py comment.
BG_CACHE_ARGS=()
if [ -n "${BOLTZGEN_CACHE:-}" ]; then
    mkdir -p "$BOLTZGEN_CACHE"
    BG_CACHE_ARGS=(--cache "$BOLTZGEN_CACHE")
fi

mkdir -p "$BG_OUT"

echo "[$(date +%T)] task $SAPIA_TASK_ID: boltzgen run '$BG_NAME'"
echo "  spec:   $BG_SPEC"
echo "  output: $BG_OUT"
echo "  cache:  ${BOLTZGEN_CACHE:-(tool default)}"
echo "  proto:  $BG_PROTOCOL"
echo "  steps:  $BG_STEPS"
echo "  flags:  ${BG_FLAGS:-(none)}"

# Validate the spec BEFORE spending four GPU steps on it.
#
# This is the only INDEPENDENT test that the YAML run_boltzgen.py writes is the shape
# BoltzGen's schema actually wants, and that the epitope landed where we think it did.
# Without it, a schema parser that ignores an unknown key would turn a nesting mistake
# into "no epitope was requested" -- a run that steers nothing while
# boltzgen_binding_resnames still reports an epitope, which is exactly the class of
# failure that column exists to prevent. `check` also writes a colour-coded mmCIF of
# the binding site to $BG_OUT/spec_check/ for a human to open.
#
# `set -e` makes a non-zero exit kill the task, which is what we want: a spec that does
# not validate must not become a 12-hour run.
echo "[$(date +%T)] validating spec with 'boltzgen check'"
"$BG_BIN" check "$BG_SPEC" --output "$BG_OUT/spec_check" "${BG_CACHE_ARGS[@]}"
echo "[$(date +%T)] spec OK"

# The step list is the TOOL'S CONTRACT, derived from the protocol by run_boltzgen.py
# and carried in the manifest -- NOT a per-run choice:
#   * esmfold2_scoring is omitted because it needs ~33 GB of weights this tool never
#     stages and contributes nothing to generation;
#   * analysis would crash with esmfold2_metrics left at its true default
#     (analyze.py:638-652 reads esmfold2_scores/<id>.input.json unconditionally), so
#     run_boltzgen.py always appends `--config analysis esmfold2_metrics=false`;
#   * filtering is omitted because it CANNOT be made to run: filter.py:218 sets
#     use_esmfold2 = not use_affinity with no config key, and use_affinity is False
#     for BOTH protocols here, so filter.py:421 raises. Running default steps is a
#     guaranteed crash at the last step, under either protocol;
#   * design_folding is INCLUDED for protein-anything and omitted for
#     nanobody-anything, because upstream derives designfolding_metrics from the
#     protocol while --steps only filters the step list. See run_boltzgen.py.
# BG_FLAGS and BG_STEPS are intentionally unquoted: both are space-separated lists of
# CLI tokens that must each become its own argv entry.
#
# `set -f` for the duration, because unquoted expansion is also GLOB expansion: a flag
# value containing *, ? or [ would be matched against the task's cwd and silently
# replaced by filenames (or left alone if nothing matched, which is worse -- it works
# until the day a file happens to match). Word-splitting is what we want; globbing is
# not, and the two are the same unquoting.
set -f
"$BG_BIN" run "$BG_SPEC" \
    --output "$BG_OUT" \
    --protocol "$BG_PROTOCOL" \
    --steps $BG_STEPS \
    "${BG_CACHE_ARGS[@]}" \
    $BG_FLAGS
set +f
