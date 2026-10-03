#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) placing each C2 dock back into the binder/target
frame and measuring occlusion of the target binding site.

Premise, which is the tool's whole contract: *a binder that is known to bind the
target was docked against itself with C2 symmetry; place that dock back into the
binder/target frame through protomer A and ask whether the partner protomer
occludes the target binding site.*

Each task runs dimerfit_worker.py on a batch of designs. Per design it

  1. reads the C2 dock (exactly two protomers of the same binder),
  2. Kabsch-superposes dock protomer A onto the reference complex's BINDER chain,
  3. applies that ONE transform to the WHOLE dimer -- so protomer B lands wherever
     the C2 operator put it relative to the target,
  4. measures occlusion of the target binding site by protomer B and the C2
     interface's position relative to the epitope,
  5. writes the transformed dimer (no target) and the dimer plus target.

``action: update`` -- this is a property of docks that already exist, so it annotates
the dock table in place rather than minting a generation. ``dimerfit_path`` (the
transformed dimer) is what the next step consumes; ``dimerfit_complex_path`` exists
to be looked at.

What it does NOT do
-------------------
* It does **not** detect any interface on the target: the epitope lists are read from
  columns (``--epitope-column`` / ``--target-epitope-column``, normally ifacegeom's),
  resolved up the lineage from the dock table to the complex table.
* It does **not** compute an energy. ``n_clash``/``clash_frac`` are heavy-atom
  distance counts -- not ``fa_rep``, and there is no Rosetta in this image.
* It does **not** relax, repack or rescore anything.
* It does **not** say whether the two protomers can be linked. That is ``linkpath``,
  which consumes ``dimerfit_path`` (the transformed dimer, no target) and measures an
  obstruction-aware path through solvent-accessible space rather than a straight
  line.
* It is **not** a symmetric-docking scorer (that is ``rpxdock``) and **not** an
  assembly-fit test (that is ``ringfit``, whose premise -- a binder straddling two
  adjacent protomers of a larger oligomer -- is a different question entirely).

``default_input_column`` is the ``"not applicable"`` sentinel (as in ifacegeom, cms
and chainsel): there is no honest default dock column, and a wrong one fails
silently, so the builder refuses the run unless ``-i`` names a column the table has.

Several designs are packed into one task (``--designs-per-task``): the work is a
Kabsch fit and a few distance scans, well under a second per design, while a
container cold start costs ~10 s.

Usage:
    sapia run dimerfit outputs/20261002_143419_dimer_phase2 \
        --table table1 \
        --input-column rpxdock_path \
        --ref-column input_path \
        --epitope-column ifacegeom_binder_res \
        --target-epitope-column ifacegeom_target_res \
        --binder-chain-in-ref B --target-chains-in-ref A
"""

from argparse import ArgumentParser
from pathlib import Path
from typing import Any, cast

import pandas as pd
from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path

# The sentinel default_input_column (see the module docstring): matches no column, so
# the builder's check below always fires unless -i was given.
NO_DEFAULT_COLUMN = "not applicable"
MATCH_MODES = ("auto", "ordinal", "resnum")


class DimerfitArgs(CommonArgs):
    ref_column: str
    epitope_column: str
    target_epitope_column: str
    binder_chain_in_ref: str
    target_chains_in_ref: str
    dock_chains: str
    resnum_match: str
    clash_cutoff: float
    occlusion_cutoff: float
    dimer_contact_cutoff: float
    designs_per_task: int


def add_run_dimerfit_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--ref-column",
        type=str,
        default="input_path",
        help="Column holding the REFERENCE binder/target complex -- the frame the "
        "dock is placed back into. Resolved up the lineage, so a dock table may "
        "name a column that lives on its parent complex table. Default 'input_path'.",
    )
    parser.add_argument(
        "--epitope-column",
        type=str,
        default="ifacegeom_binder_res",
        help="Column holding the BINDER-side epitope as 'B:12,B:15,...' (chain:"
        "resnum, ifacegeom's format -- NOT 'B12'), in the REFERENCE binder's "
        "numbering. Mapped onto the dock through the residue match, and used for "
        "epitope_com_dist / n_overlap_res / frac_overlap. Resolved up the lineage. "
        "Default 'ifacegeom_binder_res'.",
    )
    parser.add_argument(
        "--target-epitope-column",
        type=str,
        default="ifacegeom_target_res",
        help="Column holding the TARGET-side binding site as 'A:12,A:15,...' in the "
        "reference target's own numbering -- the residues protomer B has to cover. "
        "Drives occluded_frac / n_occluded_res. Read straight off the reference "
        "structure (no mapping involved). Resolved up the lineage. Default "
        "'ifacegeom_target_res'.",
    )
    parser.add_argument(
        "--binder-chain-in-ref",
        type=str,
        default="B",
        help="The BINDER's chain ID in the reference complex. Default 'B' -- the "
        "inherited hEGFR complexes have the binder on chain B and the target on "
        "chain A, the opposite of the usual convention, and getting this wrong "
        "superposes the dock onto the target. A chain named but absent is an error "
        "for that design, never a smaller selection.",
    )
    parser.add_argument(
        "--target-chains-in-ref",
        type=str,
        default="A",
        help="The TARGET's chain IDs in the reference complex, comma-joined. 'auto' "
        "takes every protein chain that is not the binder chain. Default 'A'.",
    )
    parser.add_argument(
        "--dock-chains",
        type=str,
        default="auto",
        help="The dock's two protomer chain IDs, comma-joined and ORDERED: the first "
        "is protomer A (superposed onto the reference binder), the second protomer B "
        "(the partner whose occlusion is measured). 'auto' (default) uses the dock "
        "file's own chain order. The dock must hold exactly two protein chains "
        "either way.",
    )
    parser.add_argument(
        "--resnum-match",
        type=str,
        choices=MATCH_MODES,
        default="auto",
        help="How dock residues are paired with reference-binder residues for the "
        "superposition. 'ordinal' pairs the i-th CA with the i-th (right for a dump "
        "renumbered 1..N, which is what rpxdock writes); 'resnum' pairs equal "
        "residue numbers; 'auto' (default) takes ordinal when the two chains hold "
        "the same number of CA atoms and resnum otherwise. The pairing is always "
        "checked by residue identity and reported as dimerfit_seq_match_frac.",
    )
    parser.add_argument(
        "--clash-cutoff",
        type=float,
        default=2.5,
        help="Heavy-atom distance (A) below which a protomer-B atom counts as "
        "clashing with the target. Default 2.5 -- below any real contact distance "
        "(a C-C van der Waals contact is ~3.4 A, the shortest heavy-atom H-bond "
        "~2.6 A), so a count above zero means genuine interpenetration rather than "
        "a tight interface. Same default as ringfit's, deliberately.",
    )
    parser.add_argument(
        "--occlusion-cutoff",
        type=float,
        default=5.0,
        help="Heavy-atom distance (A) within which a TARGET epitope residue counts "
        "as occluded by protomer B. Default 5.0 -- the first-shell contact distance "
        "(ringfit uses the same for target contacts), so 'occluded' means 'protomer "
        "B is in contact with it', not merely 'nearby'.",
    )
    parser.add_argument(
        "--dimer-contact-cutoff",
        type=float,
        default=5.0,
        help="Heavy-atom distance (A) within which a protomer-A residue counts as "
        "part of the C2 interface. Default 5.0, matching --occlusion-cutoff so the "
        "two residue sets are comparable (frac_overlap compares them directly).",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=100,
        help="Designs measured per task (default 100). The work is well under a "
        "second per design, so this exists to amortise the container cold start, "
        "not the computation. Lower it for more parallelism.",
    )


def _split(value: str) -> list[str]:
    """Comma-joined list -> stripped, non-empty tokens."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


def _lineage_value(ctx: ManifestCtx[DimerfitArgs], name: str, column: str) -> Any:
    """One design's value for ``column``, walking up the lineage.

    Raises rather than defaulting: every column read here (the reference complex, the
    two epitope lists) is load-bearing, and a row that silently dropped out would
    look like a filter problem rather than like the missing column it is.
    """
    value = ctx.lookup(name, column)
    if value is None or (isinstance(value, float) and pd.isna(value)):
        raise ValueError(
            f"{name}: column {column!r} is not set on this row or any ancestor. "
            f"dimerfit needs the reference complex and both epitope lists; check the "
            f"--ref-column / --epitope-column / --target-epitope-column names "
            f"against the parent table (they normally live on the table ifacegeom "
            f"annotated)."
        )
    text = str(value).strip()
    if not text or text.lower() in ("nan", "na"):
        raise ValueError(f"{name}: column {column!r} is empty.")
    return text


def build_dimerfit_manifest(ctx: ManifestCtx[DimerfitArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool

    if ctx.args.designs_per_task < 1:
        raise ValueError("--designs-per-task must be >= 1.")
    for flag, value in (
        ("--clash-cutoff", ctx.args.clash_cutoff),
        ("--occlusion-cutoff", ctx.args.occlusion_cutoff),
        ("--dimer-contact-cutoff", ctx.args.dimer_contact_cutoff),
    ):
        if value <= 0:
            raise ValueError(f"{flag} must be > 0 (got {value}).")

    binder_chain = ctx.args.binder_chain_in_ref.strip()
    if not binder_chain:
        raise ValueError("--binder-chain-in-ref must name one chain.")
    target_spec = ctx.args.target_chains_in_ref.strip()
    if target_spec.lower() != "auto":
        targets = _split(target_spec)
        if not targets:
            raise ValueError(
                "--target-chains-in-ref must name at least one chain, or 'auto'."
            )
        if binder_chain in targets:
            raise ValueError(
                f"chain {binder_chain!r} is both --binder-chain-in-ref and a "
                f"--target-chains-in-ref chain; the two sides must be disjoint."
            )
        target_spec = ",".join(targets)

    dock_spec = ctx.args.dock_chains.strip()
    if dock_spec.lower() != "auto":
        dock = _split(dock_spec)
        if len(dock) != 2:
            raise ValueError(
                f"--dock-chains must name exactly 2 chains (protomer A then "
                f"protomer B) or 'auto', got {dock or 'none'}."
            )
        if dock[0] == dock[1]:
            raise ValueError(f"--dock-chains names the same chain twice: {dock}.")
        dock_spec = ",".join(dock)

    column = ctx.args.input_column
    if column == NO_DEFAULT_COLUMN or column not in ctx.df.columns:
        available = ", ".join(str(c) for c in ctx.df.columns if str(c).endswith("_path"))
        raise ValueError(
            f"dimerfit has no default input column: pass -i/--input-column with the "
            f"column holding the C2 DOCK (got {column!r}, which table "
            f"'{ctx.args.table}' does not have). Structure columns available: "
            f"{available or '(none)'}."
        )

    ready = ctx.ready
    members: list[tuple[str, str, str, str, str]] = []
    for name in ready.index:
        name = cast(str, name)
        dock_path = Path(str(ready.at[name, column]))
        if not dock_path.exists():
            print(f"{name}: MISSING {dock_path} (skipping)")
            continue
        # Reference and epitopes come from the parent complex table through lineage;
        # a missing one raises here, at submit time, which is the cheap place to fail.
        ref = Path(str(_lineage_value(ctx, name, ctx.args.ref_column)))
        if not ref.exists():
            raise FileNotFoundError(
                f"{name}: --ref-column {ctx.args.ref_column!r} points at {ref}, which "
                f"does not exist. Under --executor modal it must be a path on the "
                f"runs volume (the workstation's /runs)."
            )
        epitope = str(_lineage_value(ctx, name, ctx.args.epitope_column))
        target_epitope = str(
            _lineage_value(ctx, name, ctx.args.target_epitope_column)
        )
        members.append(
            (
                name,
                str(volume_path(dock_path)),
                str(volume_path(ref)),
                epitope,
                target_epitope,
            )
        )

    # One sub-manifest per task, and a top-level row per task pointing at it. The row
    # carries the task's design COUNT so the worker can refuse a half-visible task
    # file instead of quietly running short (see the guard in dimerfit.sh).
    tasks_dir = ctx.out_dir / "dimerfit_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    x = ctx.args.designs_per_task
    manifest_rows: list[tuple[str, ...]] = []
    for t, i in enumerate(range(0, len(members), x)):
        batch = members[i : i + x]
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            for row in batch:
                f.write("\t".join(row) + "\n")
        manifest_rows.append(
            (
                str(volume_path(task_file)),
                str(len(batch)),
                dock_spec,
                binder_chain,
                target_spec,
                ctx.args.resnum_match,
                str(ctx.args.clash_cutoff),
                str(ctx.args.occlusion_cutoff),
                # Kept last: always non-empty, so the manifest line never ends on an
                # empty field.
                str(ctx.args.dimer_contact_cutoff),
            )
        )

    return manifest_rows
