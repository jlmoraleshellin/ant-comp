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
# the --output dir, and the run-wide flags.
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

# --protocol and --steps are the TOOL'S CONTRACT, not a per-run choice:
#   * esmfold2_scoring is omitted because it needs ~33 GB of weights this tool never
#     stages and contributes nothing to generation;
#   * analysis would crash with esmfold2_metrics left at its true default
#     (analyze.py:638-652 reads esmfold2_scores/<id>.input.json unconditionally), so
#     run_boltzgen.py always appends `--config analysis esmfold2_metrics=false`;
#   * filtering is omitted because it CANNOT be made to run: filter.py:218 sets
#     use_esmfold2 from the protocol with no config key, and filter.py:421 then
#     raises. Running default steps is a guaranteed crash at the last step.
# See run_boltzgen.py's docstring. BG_FLAGS is intentionally unquoted: it is a
# space-separated list of CLI tokens that must each become its own argv entry.
#
# `set -f` for the duration, because unquoted expansion is also GLOB expansion: a flag
# value containing *, ? or [ would be matched against the task's cwd and silently
# replaced by filenames (or left alone if nothing matched, which is worse -- it works
# until the day a file happens to match). Word-splitting is what we want; globbing is
# not, and the two are the same unquoting.
set -f
"$BG_BIN" run "$BG_SPEC" \
    --output "$BG_OUT" \
    --protocol nanobody-anything \
    --steps design inverse_folding folding analysis \
    "${BG_CACHE_ARGS[@]}" \
    $BG_FLAGS
set +f
