#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) asking, per design: DID THE PREDICTED BINDER KEEP THE
REFERENCE BINDING MODE, judged by the standard docking metric?

It wraps the PUBLISHED DockQ v2 package (Mirabello & Wallner, pip ``DockQ``, pinned to
2.1.3 in modal_image.py) and only selects, maps and reports its numbers. Premise, scope
limits, the chain-mapping and residue-matching choices, every flag and how to read the
columns: the tool's skill, .claude/skills/dockq/SKILL.md (and dockq_worker.py's
docstring for the mapping contract).

``action: update`` -- the binding-mode agreement is a property of a design that already
exists. CPU only: the builder forces ``gpus_per_task = 0``.

``default_input_column`` is the ``"not applicable"`` sentinel (as in hisbond,
epitope, cms, chainsel): there is no honest default MODEL column, so the builder raises
unless ``-i/--input-column`` is given, listing the table's ``_path`` columns.

BOTH STRUCTURE COLUMNS ARE RESOLVED UP THE LINEAGE, exactly as in hisbond: ``-i``
(the MODEL complex) and ``--ref-col`` (the REFERENCE, default ``structure_path``). A
column this table HAS is read from the row's own cell (a blank cell is an error row,
never an ancestor's file); a column this table LACKS is walked up parent_name /
parent_table. A design that resolves nothing is submitted anyway with an empty field
and collects as ``error: could not resolve ...``.
"""

import math
from argparse import ArgumentParser
from pathlib import Path
from typing import cast

import pandas as pd

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.core.naming import status_column
from prosapia.utils import expand_chain_spec

# The sentinel default_input_column (as in hisbond, epitope, cms, chainsel).
NO_DEFAULT_COLUMN = "not applicable"


class DockQArgs(CommonArgs):
    ref_col: str
    binder_chain: str
    target_chains: str
    ref_binder_chain: str
    ref_target_chains: str
    allowed_mismatches: int
    residue_match: str
    min_binder_identity: float
    min_native_contacts: int
    designs_per_task: int


def _clean(value: object) -> str | None:
    """A table cell as a non-empty string, or None for null/blank/absent."""
    if value is None:
        return None
    try:
        if pd.isna(value):  # type: ignore[arg-type]
            return None
    except (TypeError, ValueError):
        return None
    text = str(value).strip()
    return text or None


def _resolve(ctx: ManifestCtx[DockQArgs], name: str, column: str) -> str | None:
    """The row's own cell if this table HAS ``column`` (blank -> None, never an
    ancestor's value); otherwise the nearest ancestor's (the hisbond rule)."""
    if column in ctx.df.columns:
        if name not in ctx.df.index:
            return None
        return _clean(ctx.df.at[name, column])
    return _clean(ctx.lookup(name, column))


def _chains(flag: str, spec: str) -> list[str]:
    chains = expand_chain_spec(spec)
    if not chains:
        raise ValueError(f"{flag} must name at least one chain, got {spec!r}.")
    if len(set(chains)) != len(chains):
        raise ValueError(f"{flag} lists a chain twice: {chains}.")
    return chains


def add_run_dockq_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--ref-col",
        type=str,
        default="structure_path",
        help="Table column holding the REFERENCE complex, resolved up the lineage. "
        "Default 'structure_path' (the campaign's input complex on the root table). "
        "Must differ from -i. Recorded per row as dockq_ref_path.",
    )
    parser.add_argument(
        "--binder-chain",
        type=str,
        required=True,
        help="The binder chain in the MODEL, exactly one chain ID (auth_asym_id). NO "
        "DEFAULT: a defaulted letter would score the wrong molecule.",
    )
    parser.add_argument(
        "--target-chains",
        type=str,
        required=True,
        help="The target chains in the MODEL, chain mini-language ('A,B,C' or 'A:C'). "
        "NO DEFAULT. Must not contain --binder-chain.",
    )
    parser.add_argument(
        "--ref-binder-chain",
        type=str,
        required=True,
        help="The binder chain in the REFERENCE, exactly one chain ID. NO DEFAULT.",
    )
    parser.add_argument(
        "--ref-target-chains",
        type=str,
        required=True,
        help="The target chains in the REFERENCE, same count as --target-chains. NO "
        "DEFAULT. Which model copy maps onto which reference copy is SEARCHED, not "
        "taken from the order given here.",
    )
    parser.add_argument(
        "--allowed-mismatches",
        type=int,
        default=0,
        help="DockQ's allowed_mismatches, applied to the TARGET chains only (the binder "
        "is mapped explicitly). Default 0: the predicted target is the reference target "
        "sequence. A target pair differing at more positions is an error.",
    )
    parser.add_argument(
        "--residue-match",
        choices=("number", "sequence"),
        default="number",
        help="How residues correspond inside a mapped chain pair. 'number' (DEFAULT; "
        "DockQ's --no_align): by residue number, every position compared, numbering "
        "must be identical or the row errors. 'sequence' (DockQ's own default): by "
        "sequence alignment, and DockQ then DROPS every substituted position -- a "
        "design's mutated interface residues leave iRMSD/LRMSD and can never count "
        "as reproduced contacts, while still counting in Fnat's denominator.",
    )
    parser.add_argument(
        "--min-binder-identity",
        type=float,
        default=0.5,
        help="Minimum fraction of binder positions identical between model and "
        "reference (default 0.5). Below it the row errors: the model binder is not a "
        "variant of the reference binder (wrong reference or wrong chain). Designs "
        "with 30 %% substitutions sit at 0.7.",
    )
    parser.add_argument(
        "--min-native-contacts",
        type=int,
        default=1,
        help="Native binder interfaces with fewer native contacts (DockQ's nat_total: "
        "residue pairs within 5 A in the reference) are EXCLUDED from dockq_best, "
        "dockq_mean_natw and n_interfaces_acceptable, and kept in per_interface "
        "marked used=0. Default 1 = DockQ's own behaviour (every native interface "
        "counts). On a 1-2 contact interface Fnat can only be 0, 0.5 or 1, so a wrong "
        "pose can pass the 0.23 gate there; choose N from the per-interface nat "
        "counts of the parents. A row with no interface at or above N errors.",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=20,
        help="Designs per task (default 20; about 0.4 s each).",
    )


def build_dockq_manifest(ctx: ManifestCtx[DockQArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool: callers never need -g 0

    binder = _chains("--binder-chain", ctx.args.binder_chain)
    ref_binder = _chains("--ref-binder-chain", ctx.args.ref_binder_chain)
    for flag, got in (("--binder-chain", binder), ("--ref-binder-chain", ref_binder)):
        if len(got) != 1:
            raise ValueError(f"{flag} must name exactly one chain, got {got}.")
    targets = _chains("--target-chains", ctx.args.target_chains)
    ref_targets = _chains("--ref-target-chains", ctx.args.ref_target_chains)
    if binder[0] in targets:
        raise ValueError(
            f"--binder-chain {binder[0]} is also in --target-chains {targets}: the "
            f"binder cannot be its own target."
        )
    if ref_binder[0] in ref_targets:
        raise ValueError(
            f"--ref-binder-chain {ref_binder[0]} is also in --ref-target-chains "
            f"{ref_targets}: the binder cannot be its own target."
        )
    if len(targets) != len(ref_targets):
        raise ValueError(
            f"--target-chains has {len(targets)} chains and --ref-target-chains "
            f"{len(ref_targets)}: the model and the reference must hold the same target."
        )
    if ctx.args.allowed_mismatches < 0:
        raise ValueError("--allowed-mismatches must be >= 0.")
    identity = float(ctx.args.min_binder_identity)
    if not math.isfinite(identity) or not 0.0 <= identity <= 1.0:
        raise ValueError(f"--min-binder-identity must be in [0, 1], got {identity}.")
    if ctx.args.designs_per_task < 1:
        raise ValueError("--designs-per-task must be >= 1.")
    if ctx.args.min_native_contacts < 1:
        raise ValueError("--min-native-contacts must be >= 1 (1 = DockQ's own rule).")

    column = ctx.args.input_column
    if column == NO_DEFAULT_COLUMN:
        available = ", ".join(str(c) for c in ctx.df.columns if str(c).endswith("_path"))
        raise ValueError(
            f"dockq has no default input column: pass -i/--input-column with the MODEL "
            f"complex column (e.g. boltz_cof_path). Structure columns in table "
            f"'{ctx.args.table}': {available or '(none)'} -- a parent table's column is "
            f"also valid and is resolved up the lineage."
        )
    ref_col = ctx.args.ref_col.strip()
    if not ref_col:
        raise ValueError("--ref-col must name a table column.")
    if ref_col == column:
        raise ValueError(
            f"-i and --ref-col are both {column!r}: every design would be compared "
            f"with itself."
        )

    # NOT ctx.ready (it drops rows with a blank input cell; those must become error
    # rows). The resume skip is copied from ManifestCtx.ready so --force matches.
    ready = ctx.df
    status_col = status_column(ctx.out_dir.name)
    if not ctx.args.force and status_col in ready.columns:
        ready = ready[ready[status_col] != "OK"]

    members: list[tuple[str, str, str]] = []
    n_no_model = n_no_ref = 0
    refs: set[str] = set()
    for name in ready.index:
        name = cast(str, name)
        model = _resolve(ctx, name, column)
        ref = _resolve(ctx, name, ref_col)
        if model is None:
            n_no_model += 1
            print(f"{name}: could not resolve {column} -- submitting so the task records an error")
        if ref is None:
            n_no_ref += 1
            print(f"{name}: could not resolve {ref_col} -- submitting so the task records an error")
        else:
            refs.add(ref)
        # CIF is passed through as-is: DockQ's own mmCIF parser reads it.
        members.append(
            (
                name,
                "" if model is None else str(volume_path(Path(model))),
                "" if ref is None else str(volume_path(Path(ref))),
            )
        )

    if members and n_no_model == len(members):
        raise ValueError(
            f"dockq resolved no model structure for any of the {len(members)} ready "
            f"designs via {column!r}. Check -i/--input-column."
        )
    if members and n_no_ref == len(members):
        raise ValueError(
            f"dockq resolved no reference structure for any of the {len(members)} "
            f"ready designs via {ref_col!r}. Check --ref-col."
        )
    print(
        f"dockq: {len(members)} design(s) against {len(refs)} distinct reference "
        f"structure(s) from {ref_col!r} -- expect one per parent complex. "
        f"Unresolved: {n_no_model} model, {n_no_ref} reference. A row whose OWN "
        f"{column} cell is blank is skipped by collect and keeps an EMPTY status."
    )

    # Stale per-design results go only after every refusal above.
    for name, _m, _r in members:
        (ctx.out_dir / f"{name}.tsv").unlink(missing_ok=True)

    tasks_dir = ctx.out_dir / "dockq_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    per_task = ctx.args.designs_per_task
    rows: list[tuple[str, ...]] = []
    for t, i in enumerate(range(0, len(members), per_task)):
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            for name, model, ref in members[i : i + per_task]:
                f.write(f"{name}\t{model}\t{ref}\n")
        rows.append(
            (
                # Every field is non-empty, so none can vanish from `cut`.
                str(volume_path(task_file)),
                binder[0],
                ",".join(targets),
                ref_binder[0],
                ",".join(ref_targets),
                str(int(ctx.args.allowed_mismatches)),
                ctx.args.residue_match,
                str(identity),
                column,
                ref_col,
                # Appended as field 11; dockq.sh defaults a missing one to 1.
                str(int(ctx.args.min_native_contacts)),
            )
        )
    return rows
