#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) recording, per design, the SHAPE and SECONDARY-
STRUCTURE MAKE-UP of a chosen chain selection: the CA radius of gyration (mass-
unweighted) and the 3-state DSSP helix / strand / loop fractions.

PREMISE AND SCOPE. One structure file per design, and a chain selection
(``--design-chains``) that is a protein: one or more chains, each ONE CONTINUOUS
POLYPEPTIDE with a complete N/CA/C/O backbone. The numbers describe that selection's
OWN backbone and nothing else:

  * It never looks at the unselected chains. No interface, no clash, no pose check,
    no claim about binding. A binder strand that pairs only with a TARGET strand is
    assigned as loop here, because the target is not part of the DSSP calculation.
  * No energetics, no relaxation, no designability. A compact all-helix chain can
    be a terrible design.
  * The numbers are meaningless if the selection is the WRONG chain (read
    ``ssprofile_sequence`` and ``ssprofile_n_res`` against what you expect; chain
    letters are AUTHOR chain IDs for PDB and mmCIF alike), or if the structure is a
    prediction whose fold you have not yet checked (the fractions then describe a
    fiction).

``action: update`` -- Rg and SS composition are properties of designs that already
exist; this annotates the table in place and mints no child. Run it twice with
different ``-l/--dir-label`` to profile the designed backbone and its prediction onto
one table.

``default_input_column`` is the ``"not applicable"`` sentinel (as in cms, chainsel,
epitope): ``rfdiffusion3_path``, ``boltz_path``, ``bindcraft2_path`` and
``chainsel_path`` are all normal inputs, so there is no honest default and the
builder raises unless ``-i/--input-column`` is given.

THE STRUCTURE COLUMN IS RESOLVED UP THE LINEAGE, exactly as in epitope: an update
tool is handed a plain read of its own table (no join_lineage), so a sequence-design
child table does not carry its parent's ``rfdiffusion3_path``. Each path is resolved
row first, then up ``parent_name``/``parent_table``. A design whose structure
resolves NOWHERE is submitted anyway with an empty path and lands in the table as
``error: could not resolve ...`` -- never dropped from the manifest.

Chain selection uses the prosapia chain mini-language (``,`` separates, ``:`` is an
inclusive letter range): ``D``, ``A,B``, ``A:D``.

Usage:
    sapia run ssprofile outputs/20261009_binders \
        --table table0 \
        --input-column rfdiffusion3_path \
        --design-chains B \
        --dir-label bb
"""

from argparse import ArgumentParser
from pathlib import Path
from typing import cast

import pandas as pd

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.core.naming import status_column
from prosapia.utils import ensure_pdb, expand_chain_spec

# The sentinel default_input_column (see the module docstring).
NO_DEFAULT_COLUMN = "not applicable"


class SsprofileArgs(CommonArgs):
    design_chains: str
    designs_per_task: int


def _clean(value: object) -> str | None:
    """A table cell as a non-empty string, or None for null/blank/absent."""
    if value is None:
        return None
    try:
        if pd.isna(value):  # type: ignore
            return None
    except (TypeError, ValueError):  # arrays and other non-scalars
        return None
    text = str(value).strip()
    return text or None


def _resolve(ctx: ManifestCtx[SsprofileArgs], name: str, column: str) -> str | None:
    """The row's own value for ``column``, else the nearest ancestor's, else None.
    Row first, lineage only as a fallback -- the epitope / usalign --col-b precedent."""
    if column in ctx.df.columns and name in ctx.df.index:
        if (own := _clean(ctx.df.at[name, column])) is not None:
            return own
    return _clean(ctx.lookup(name, column))


def parse_design_chains(spec: str) -> list[str]:
    """Chain mini-language -> ordered, validated chain letters. Raises on an empty
    selection or a repeated chain (which would double-count its residues)."""
    chains = expand_chain_spec(spec)
    if not chains:
        raise ValueError("--design-chains must name at least one chain (e.g. 'B').")
    if len(set(chains)) != len(chains):
        raise ValueError(f"--design-chains lists a chain twice: {chains}.")
    return chains


def add_run_ssprofile_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--design-chains",
        type=str,
        required=True,
        help="The chain(s) to profile, in the chain mini-language (',' separates, "
        "':' is an inclusive letter range): 'D', 'A,B', 'A:D'. Required -- there is "
        "no default chain. Only these chains' CA atoms enter Rg and only their "
        "backbone enters DSSP. A chain named here but ABSENT from a design is an "
        "'error:' for that design with every column NA -- never a fallback to "
        "another chain. With several chains, Rg is of their union and the fractions "
        "are over all their residues together.",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=20,
        help="Designs profiled per task (default 20). The work is well under a "
        "second per design, so a task per design would be nearly all container "
        "start-up. Lower it for more parallelism.",
    )


def build_ssprofile_manifest(ctx: ManifestCtx[SsprofileArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool: callers never need -g 0

    chains = parse_design_chains(ctx.args.design_chains)
    if ctx.args.designs_per_task < 1:
        raise ValueError("--designs-per-task must be >= 1.")

    column = ctx.args.input_column
    if column == NO_DEFAULT_COLUMN:
        available = ", ".join(str(c) for c in ctx.df.columns if str(c).endswith("_path"))
        raise ValueError(
            f"ssprofile has no default input column: pass -i/--input-column with the "
            f"structure column to profile. Structure columns in table "
            f"'{ctx.args.table}': {available or '(none)'} -- a PARENT table's column "
            f"(e.g. rfdiffusion3_path) is also valid and is resolved up the lineage."
        )

    # NOT ctx.ready: filter_ready drops rows whose input cell is null or blank, so a
    # design whose structure lives on a parent (column absent here) would vanish and
    # the run would silently be smaller than the table. Every row is a candidate; the
    # resume skip below is copied from ManifestCtx.ready so --force behaves the same.
    ready = ctx.df
    status_col = status_column(ctx.out_dir.name)
    if not ctx.args.force and status_col in ready.columns:
        ready = ready[ready[status_col] != "OK"]

    members: list[tuple[str, str]] = []
    n_unresolved = 0
    for name in ready.index:
        name = cast(str, name)
        raw = _resolve(ctx, name, column)
        if raw is None:
            # Submitted with an EMPTY structure field on purpose: the worker turns it
            # into 'error: could not resolve ...' IN THE TABLE.
            n_unresolved += 1
            print(
                f"{name}: could not resolve {column} (row or lineage) -- submitting "
                f"anyway so it lands as an error, not as a missing row"
            )
            members.append((name, ""))
            continue
        src = Path(raw)
        # CIF (and .cif.gz) -> PDB up front, cached under run_dir/.cif_to_pdb. A path
        # that does not exist is passed through RAW so the worker records it as
        # 'error: structure missing: ...' rather than the run dropping the design.
        staged = ensure_pdb(src, ctx.args.run_dir) if src.exists() else src
        members.append((name, str(volume_path(staged))))

    if n_unresolved:
        print(
            f"ssprofile: {n_unresolved}/{len(members)} design(s) had no {column} on "
            f"their row or any ancestor; each will collect as 'error: could not "
            f"resolve ...'."
        )
    if members and n_unresolved == len(members):
        raise ValueError(
            f"ssprofile resolved no structure at all for any of the {len(members)} "
            f"ready designs via column {column!r} (searched table "
            f"'{ctx.args.table}' and its ancestors). Check -i/--input-column."
        )

    # One sub-manifest per task, and a top-level row per task pointing at it.
    tasks_dir = ctx.out_dir / "ssprofile_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    per_task = ctx.args.designs_per_task
    manifest_rows: list[tuple[str, ...]] = []
    for t, i in enumerate(range(0, len(members), per_task)):
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            for name, src in members[i : i + per_task]:
                f.write(f"{name}\t{src}\n")
        manifest_rows.append(
            (
                str(volume_path(task_file)),
                ",".join(chains),
                # Always non-empty, so it cannot vanish from `cut`. New fields go on
                # the END: ssprofile.sh cuts by index.
                column,
            )
        )
    return manifest_rows
