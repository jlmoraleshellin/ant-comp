# Activation for the pure-Python custom tools on the VIB DataCore:
# `mkcomplex`, `chainsel`, `ringfit`, `ifacegeom`, `dimerfit`, `af3conf`, `boltzconf`,
# `epitope`, `hisbond`, `linkpath`, `trifuse`. Point each one's SAPIA_ACTIVATE_<NAME> here.
# NOT for `propka` (needs the propka package) or `dockq` (needs DockQ, which pins numpy < 2).
#
# This file is the SLURM counterpart of a tool's `modal_image.py` — same job
# (provide the runtime), different mechanism. It is sourced by the tool's task
# script via `sapia_activate`, on the compute node, before the worker runs.
#
# These three tools need nothing beyond prosapia's own dependencies, which already
# include gemmi, numpy and pandas. So there is no separate env to build: the
# workspace venv that runs `sapia` is also the one that runs the workers.
#
#   mkcomplex  bash + awk only (its sapia_activate call is optional)
#   chainsel   gemmi
#   ringfit    gemmi + numpy
#
# SAPIA_VIB_WORKSPACE comes from the workspace .env, which the prelude sources
# before calling sapia_activate. Deriving from it keeps this script free of any
# one user's paths, so the repo stays portable between collaborators.

source "${SAPIA_VIB_WORKSPACE:?set SAPIA_VIB_WORKSPACE in the workspace .env}/.venv/bin/activate"

# Pin the interpreter rather than relying on `PY=${PIPELINE_PYTHON:-python}` in the
# task scripts: `python` on a compute node is whatever the node's PATH happens to
# resolve, which is not necessarily this venv.
export PIPELINE_PYTHON="${SAPIA_VIB_WORKSPACE}/.venv/bin/python"
