#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) asking, per design: how confident is Boltz about the
BINDER and its interface, not the whole complex? Flags, columns and how to read them:
``.claude/skills/boltzconf/SKILL.md``.

``action: update``: confidence is a property of an existing design, so the table is
annotated in place; ``-t`` is required. CPU only (the builder forces
``gpus_per_task = 0``).

``default_input_column`` is the ``"not applicable"`` sentinel (as in epitope, cms,
chainsel): a run_dir holds several boltz leaves (``boltz_path``, ``boltz_full8b_path``,
...), none a defensible default, so the builder raises unless ``-i/--input-column``
names the column holding the predicted STRUCTURE; the sibling files are found from it.

The lineage walk is OFF by default, which is where boltzconf DIFFERS from epitope: a
Boltz PREDICTION belongs to one sequence, so a child table reading its parent's
``boltz_path`` would score the PARENT's prediction under the child's name. An absent
column raises unless ``--allow-lineage`` is given, and even then ``input_source`` and
``stem_matches_name`` record what happened. A blank cell in a table that HAS the column
never borrows an ancestor's file; ``filter_ready`` drops such a row at collect, so its
``boltzconf_status`` stays EMPTY rather than ``error:`` (count non-empty statuses
against the row count).

Before submitting, the builder deletes ``<out_dir>/<name>.tsv`` for every design, so a
killed task collects as ``missing``, never as the previous run's ``OK``.
"""

from argparse import ArgumentParser
from pathlib import Path
from typing import cast

import pandas as pd

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.core.naming import status_column
from prosapia.utils import expand_chain_spec

# The sentinel default_input_column (see the module docstring).
NO_DEFAULT_COLUMN = "not applicable"


class BoltzconfArgs(CommonArgs):
    binder_chains: str
    target_chains: str
    contact_cutoff: float
    ipsae_pae_cutoff: float
    designs_per_task: int
    allow_lineage: bool


def _clean(value: object) -> str | None:
    """A table cell as a non-empty string, or None for null/blank/absent."""
    if value is None:
        return None
    try:
        if pd.isna(value):  # type: ignore[arg-type]
            return None
    except (TypeError, ValueError):  # arrays and other non-scalars
        return None
    text = str(value).strip()
    return text or None


def _resolve(
    ctx: ManifestCtx[BoltzconfArgs], name: str, column: str, allow_lineage: bool
) -> tuple[str | None, str]:
    """(structure path or None, the source it came from).

    A table that HAS ``column`` answers from the row's own cell, and a blank never falls
    back to an ancestor. The lineage is walked only when the column is ABSENT and
    ``--allow-lineage`` was given (a deliberate difference from epitope).
    """
    if column in ctx.df.columns:
        if name not in ctx.df.index:
            return None, "own table"
        return _clean(ctx.df.at[name, column]), "own table"
    if not allow_lineage:
        return None, "lineage (refused)"
    return _clean(ctx.lookup(name, column)), "lineage"


def add_run_boltzconf_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--binder-chains",
        type=str,
        required=True,
        help="The binder chain, in the chain mini-language. EXACTLY ONE chain: "
        "chains_ptm and pair_chains_iptm are per chain, so a two-chain binder could "
        "only be reported as an average that reads like one chain's confidence. "
        "Named plural for consistency with epitope/cms; more than one is an error. "
        "The tool does not know which chain is the binder -- check the letter "
        "against the prediction's own chain order (after mkcomplex prepended the "
        "target, the binder is the LAST chain) and boltzconf_binder_n_res against "
        "the expected binder length.",
    )
    parser.add_argument(
        "--target-chains",
        type=str,
        required=True,
        help="The target chain(s), in the chain mini-language ('A,B,C', 'A:C'). "
        "Must be disjoint from --binder-chains. The ORDER is the order of the "
        "per-chain columns (iptm_a, iptm_b, ... / ipae_a_mean, ...), named after "
        "the chain letter in LOWER CASE.",
    )
    parser.add_argument(
        "--contact-cutoff",
        type=float,
        default=5.0,
        help="Heavy-atom distance in A defining an interface residue PAIR (default "
        "5.0). An ipae_mean is only comparable with another ipae_mean at the same "
        "cutoff, which is why the value is echoed into a column. It does NOT affect "
        "the ipsae_* columns, which use no distance at all.",
    )
    parser.add_argument(
        "--ipsae-pae-cutoff",
        type=float,
        default=15.0,
        help="PAE cutoff in A for the ipsae_* columns (default 15.0). The reference "
        "SCRIPT has no built-in default -- it takes the cutoff as a required "
        "positional. The ipSAE PAPER's only recommendation, closing its Figure 7 "
        "cutoff benchmark, is that 'cutoffs of 10 or 15 A may be most suitable'; it "
        "is not predictor-specific (the paper benchmarks no AF3 or boltz cutoff), "
        "and the paper's own worked examples and command lines use 15. This default "
        "is the upper end of that band. A residue pair counts towards ipSAE when "
        "pae[i][j] < this value, STRICTLY. Echoed into ipsae_pae_cutoff: two ipSAEs "
        "computed at different cutoffs are different quantities and must never be "
        "compared or pooled.",
    )
    parser.add_argument(
        "--allow-lineage",
        action="store_true",
        help="Allow the structure column to be read from an ANCESTOR table when "
        "this table does not have it. OFF by default, unlike epitope: a boltz "
        "prediction belongs to one sequence, so a child table would be scored on "
        "its parent's prediction. Turn it on only when the rows of -t are the same "
        "molecules as the rows that were predicted, and read the input_source and "
        "stem_matches_name columns afterwards.",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=20,
        help="Designs scored per task (default 20). The work is a small matrix per "
        "design, so one task per design would be nearly all container start-up.",
    )


def build_boltzconf_manifest(ctx: ManifestCtx[BoltzconfArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool: callers never need -g 0

    binder = expand_chain_spec(ctx.args.binder_chains)
    targets = expand_chain_spec(ctx.args.target_chains)
    if len(binder) != 1:
        raise ValueError(
            f"--binder-chains must name exactly one chain, got {len(binder)} "
            f"({ctx.args.binder_chains!r}): chains_ptm and pair_chains_iptm are per "
            f"chain, so two binder chains could only be averaged, and an average "
            f"that reads like one chain's confidence is exactly the silently-wrong "
            f"number this tool exists to avoid."
        )
    if not targets:
        raise ValueError("--target-chains must name at least one chain (e.g. 'A,B,C').")
    if len(set(targets)) != len(targets):
        raise ValueError(f"--target-chains lists a chain twice: {targets}.")
    if overlap := sorted(set(binder) & set(targets)):
        raise ValueError(
            f"chain(s) {','.join(overlap)} are named as both binder and target."
        )
    keys = [c.lower() for c in targets]
    if len(set(keys)) != len(keys):
        raise ValueError(
            f"--target-chains {targets} collide once lower-cased into column names "
            f"({keys}); rename the chains before scoring them."
        )
    # Column names need one-character chain keys, else a key could shadow a fixed name.
    if bad := sorted(k for k in keys if len(k) != 1 or not k.isalnum()):
        raise ValueError(
            f"--target-chains produced multi-character or non-alphanumeric column "
            f"key(s) {bad}: per-chain column names (ipsae_max_<x>, "
            f"ipsae_nres_binder_<x>, ...) assume a single-character chain key and "
            f"would otherwise be ambiguous."
        )
    if ctx.args.contact_cutoff <= 0:
        raise ValueError("--contact-cutoff must be > 0.")
    if ctx.args.ipsae_pae_cutoff <= 0:
        raise ValueError("--ipsae-pae-cutoff must be > 0.")
    if ctx.args.designs_per_task < 1:
        raise ValueError("--designs-per-task must be >= 1.")

    column = ctx.args.input_column
    if column == NO_DEFAULT_COLUMN:
        available = ", ".join(str(c) for c in ctx.df.columns if str(c).endswith("_path"))
        raise ValueError(
            f"boltzconf has no default input column: pass -i/--input-column with the "
            f"column holding the PREDICTED STRUCTURE (the confidence JSON and the "
            f"pae/plddt npz files are found beside it). Structure columns in table "
            f"'{ctx.args.table}': {available or '(none)'} -- a PARENT table's column "
            f"is also valid and is resolved up the lineage."
        )

    # Not ctx.ready: it drops blank-input rows, which must become error rows.
    ready = ctx.df
    status_col = status_column(ctx.out_dir.name)
    if not ctx.args.force and status_col in ready.columns:
        ready = ready[ready[status_col] != "OK"]

    if column not in ctx.df.columns and not ctx.args.allow_lineage:
        available = ", ".join(str(c) for c in ctx.df.columns if str(c).endswith("_path"))
        raise ValueError(
            f"table '{ctx.args.table}' has no column {column!r}. boltzconf does NOT "
            f"walk the lineage by default: a boltz prediction belongs to one "
            f"sequence, so reading an ancestor's prediction would score the PARENT's "
            f"structure under each child's name. Structure columns here: "
            f"{available or '(none)'}. Pass --allow-lineage if the rows of this "
            f"table really are the molecules that were predicted."
        )

    members: list[tuple[str, str, str]] = []
    n_unresolved = 0
    for name in ready.index:
        name = cast(str, name)
        # A stale result from an earlier run must not survive this submit.
        (ctx.out_dir / f"{name}.tsv").unlink(missing_ok=True)
        raw, source = _resolve(ctx, name, column, ctx.args.allow_lineage)
        if raw is None:
            n_unresolved += 1
            print(
                f"{name}: could not resolve {column} (blank own cell, or absent from "
                f"this table and every ancestor) -- submitting anyway so the task "
                f"records an error"
            )
            members.append((name, "", source))
            continue
        # As-is: the file must stay beside boltz's sibling JSON/npz.
        members.append((name, str(volume_path(Path(raw))), source))

    if n_unresolved:
        print(
            f"boltzconf: {n_unresolved}/{len(members)} design(s) could not resolve "
            f"{column}. Their tasks write 'error: could not resolve ...', but a row "
            f"whose OWN {column} cell is blank is skipped by collect (filter_ready) "
            f"and keeps an EMPTY boltzconf_status -- compare non-empty statuses with "
            f"the row count."
        )
    if members and n_unresolved == len(members):
        raise ValueError(
            f"boltzconf resolved no structure at all for any of the {len(members)} "
            f"ready designs via column {column!r} (searched table "
            f"'{ctx.args.table}' and its ancestors). Check -i/--input-column."
        )

    tasks_dir = ctx.out_dir / "boltzconf_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    per_task = ctx.args.designs_per_task
    manifest_rows: list[tuple[str, ...]] = []
    for t, i in enumerate(range(0, len(members), per_task)):
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            for name, src, source in members[i : i + per_task]:
                f.write(f"{name}\t{src}\t{source}\n")
        manifest_rows.append(
            (
                str(volume_path(task_file)),
                ",".join(binder),
                ",".join(targets),
                f"{ctx.args.contact_cutoff}",
                column,
                # Appended: the .sh cut -f1..5 indices are unchanged.
                f"{ctx.args.ipsae_pae_cutoff}",
            )
        )
    return manifest_rows
