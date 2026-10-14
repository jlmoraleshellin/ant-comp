"""User-side fork of the bundled ``proteinmpnn``, adding per-row column flags.

Mechanism: rung 5 of the ``editing-a-tool`` ladder --
``get_builtin("proteinmpnn").with_overrides(...)``. Only the two hooks that had
to change are replaced:

  * ``build_manifest_fn`` -- reads ``--fixed-positions-column`` /
    ``--chains-to-design-column`` per row;
  * ``add_run_args_fn``   -- upstream's flags plus those two.

NOT overridden, and therefore still tracking upstream: ``collect_fn``,
``default_script`` (the bundled ``proteinmpnn.sh``), its sibling
``modal_image.py``, ``action``, and ``default_input_column``. The manifest's
field layout is unchanged, so the bundled task script reads it as before.

``name`` is deliberately kept as ``proteinmpnn``, so this SHADOWS the built-in
rather than registering alongside it. Shadowing keeps the output dir leaf and
every column name identical (``proteinmpnn_sequence`` stays Boltz's default
input column), which is what makes the change additive for existing callers:
with neither new flag passed the manifest rows, the sub-manifest fields and the
staged symlinks are the built-in's.

ONE deliberate departure from upstream's on-disk layout: the staged subgroup
directories are named ``grp_<digest>`` and the task sub-manifests
``task_<digest>.tsv``, instead of upstream's ``grp_0..N`` / ``task_0..N``
counters, which restart at 0 on every submit and therefore COLLIDE between two
submits to the same table+leaf. See ``run_proteinmpnn.py``'s docstring section
"SUBGROUP DIRECTORY NAMES". No downstream code reads those names: the task
script takes the grp dir from field 1 of its sub-manifest and the bundled
collector globs ``*/seqs/*.fa``.

``default_input_column`` is inherited as ``rfdiffusion_path`` -- the upstream
trap is unchanged. Coming from an rfd3 table you still need ``-i
rfdiffusion3_path``.
"""

from prosapia.core import get_builtin

from .run_proteinmpnn import (
    add_run_proteinmpnn_args,
    build_proteinmpnn_manifest,
)

TOOL = get_builtin("proteinmpnn").with_overrides(
    description=(
        "Run ProteinMPNN on a set of input PDBs. Fork: --fixed-positions-column "
        "/ --chains-to-design-column read a per-design spec from a table column."
    ),
    build_manifest_fn=build_proteinmpnn_manifest,
    add_run_args_fn=add_run_proteinmpnn_args,
)
