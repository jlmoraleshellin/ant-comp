#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) asking, per design: how confident is AlphaFold3 about
the BINDER and its interface, not the whole complex? The AF3 counterpart of boltzconf.
Flags, columns and how to read them: ``.claude/skills/af3conf/SKILL.md``.

``action: update``: confidence is a property of an existing design, so the table is
annotated in place; ``-t`` is required. CPU only (the builder forces
``gpus_per_task = 0``), and it reads files an alphafold3 leaf already wrote, so there is
no re-prediction.

``default_input_column`` is the ``"not applicable"`` sentinel (as in boltzconf): a
run_dir can hold several alphafold3 leaves, so the builder raises unless
``-i/--input-column`` names the column holding AF3's predicted STRUCTURE
(``<prefix>_model.cif``); the sibling JSONs are found from that path.

The lineage walk is OFF by default: a prediction belongs to one sequence, so a child
table would score its PARENT's prediction under its own name. An absent column raises
unless ``--allow-lineage`` is given (``input_source`` and ``stem_matches_name`` then
record what happened). A blank cell in a table that HAS the column never borrows an
ancestor's file; since ``filter_ready`` drops such a row at collect, it KEEPS any
af3conf_* values it already had (a previous ``OK`` included), so gate filters on
``alphafold3_status == "OK"``.

``--sequence-column`` (optional, normally ``mkcomplex_sequence``) follows the same
own-table-first rule; the worker compares the CIF chain named by --binder-chains with
the SAME letter's segment and writes ``binder_seq_match``. It ties the file to the
sequence for that chain and cannot tell whether that chain is the binder.

Labels: no ``_``-separated segment of ``-l`` may equal or start with a per-chain column
stem (iptm, ipsae, ipae), because the collector resets per-chain columns by prefix and
would blank a sibling leaf's. Before submitting, the builder deletes
``<out_dir>/<name>.tsv`` for every design, so a killed task collects as ``missing``,
never as the previous run's ``OK``.
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

TOOL_NAME = "af3conf"
# No dir-label segment may start with these: the collector resets them by prefix.
RESERVED_LABEL_STEMS = ("iptm", "ipsae", "ipae")


class Af3confArgs(CommonArgs):
    binder_chains: str
    target_chains: str
    ipsae_pae_cutoff: float
    designs_per_task: int
    allow_lineage: bool
    sequence_column: str | None


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
    ctx: ManifestCtx[Af3confArgs], name: str, column: str, allow_lineage: bool
) -> tuple[str | None, str]:
    """(structure path or None, the source it came from).

    A table that HAS ``column`` answers from the row's own cell, and a blank never falls
    back to an ancestor. The lineage is walked only when the column is ABSENT and
    ``--allow-lineage`` was given.
    """
    if column in ctx.df.columns:
        if name not in ctx.df.index:
            return None, "own table"
        return _clean(ctx.df.at[name, column]), "own table"
    if not allow_lineage:
        return None, "lineage (refused)"
    return _clean(ctx.lookup(name, column)), "lineage"


def add_run_af3conf_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--binder-chains",
        type=str,
        required=True,
        help="The binder chain, in the chain mini-language. EXACTLY ONE chain: "
        "chain_ptm and chain_pair_iptm are per chain, so a two-chain binder could "
        "only be reported as an average that reads like one chain's confidence. "
        "Named plural for consistency with boltzconf/epitope/cms; more than one is "
        "an error. The tool does not know which chain is the binder -- AF3's chain "
        "letters are positional in the sequence it was given (after mkcomplex "
        "prepended the target, the binder is the LAST chain); check "
        "af3conf_binder_n_res against the expected binder length.",
    )
    parser.add_argument(
        "--target-chains",
        type=str,
        required=True,
        help="The target chain(s), in the chain mini-language ('A,B,C', 'A:C'). "
        "AF3 splits a multi-letter entity (ABC) into plain chains A, B, C in its "
        "output, so 'A:C' names them. Must be disjoint from --binder-chains. The "
        "ORDER is the order of the per-chain columns (iptm_a, iptm_b, ...), named "
        "after the chain letter in LOWER CASE.",
    )
    parser.add_argument(
        "--ipsae-pae-cutoff",
        type=float,
        default=15.0,
        help="PAE cutoff in A for the ipsae_* columns (default 15.0, as boltzconf). "
        "The reference script has no built-in default; the ipSAE paper's only "
        "recommendation is that 'cutoffs of 10 or 15 A may be most suitable', and "
        "its own examples use 15. A residue pair counts when pae[i][j] < this "
        "value, STRICTLY. Echoed into ipsae_pae_cutoff: two ipSAEs computed at "
        "different cutoffs are different quantities and must never be compared or "
        "pooled -- including an af3conf and a boltzconf column side by side.",
    )
    parser.add_argument(
        "--allow-lineage",
        action="store_true",
        help="Allow the structure column to be read from an ANCESTOR table when "
        "this table does not have it. OFF by default: an AF3 prediction belongs to "
        "one sequence, so a child table would be scored on its parent's "
        "prediction. Turn it on only when the rows of -t are the same molecules as "
        "the rows that were predicted, and read input_source and stem_matches_name "
        "afterwards.",
    )
    parser.add_argument(
        "--sequence-column",
        type=str,
        default=None,
        help="Optional '/'-chainbreak sequence column for the predicted complex "
        "(normally mkcomplex_sequence -- the column alphafold3 was run on). When "
        "given, the one-letter sequence of the CIF chain named by --binder-chains "
        "is compared with the SAME letter's segment of that column (segment i is "
        "chain chr(ord('A')+i), the rule alphafold3 letters its chains by) and "
        "binder_seq_match is written 1/0. It ties the file to the sequence; it does "
        "NOT check which chain is the binder (binder_n_res does). A "
        "mismatch is a trust value, not an error; an unresolvable sequence or a "
        "missing binder segment is an error for that design. Absent (the default): "
        "binder_seq_match is NA, which means UNCHECKED. Resolved like the structure "
        "column (own table; lineage only with --allow-lineage).",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=20,
        help="Designs scored per task (default 20). The work is a small matrix per "
        "design, so one task per design would be nearly all container start-up.",
    )


def build_af3conf_manifest(ctx: ManifestCtx[Af3confArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool: callers never need -g 0

    binder = expand_chain_spec(ctx.args.binder_chains)
    targets = expand_chain_spec(ctx.args.target_chains)
    if len(binder) != 1:
        raise ValueError(
            f"--binder-chains must name exactly one chain, got {len(binder)} "
            f"({ctx.args.binder_chains!r}): chain_ptm and chain_pair_iptm are per "
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
    if ctx.args.ipsae_pae_cutoff <= 0:
        raise ValueError("--ipsae-pae-cutoff must be > 0.")
    if ctx.args.designs_per_task < 1:
        raise ValueError("--designs-per-task must be >= 1.")
    leaf = ctx.out_dir.name
    label = leaf[len(TOOL_NAME) + 1 :] if leaf.startswith(f"{TOOL_NAME}_") else ""
    bad_segments = [
        s for s in label.split("_") if s.lower().startswith(RESERVED_LABEL_STEMS)
    ]
    if bad_segments:
        raise ValueError(
            f"-l/--dir-label {label!r} has the segment(s) {bad_segments} that equal or "
            f"start with a per-chain column stem ({'/'.join(RESERVED_LABEL_STEMS)}): "
            f"the collector resets per-chain columns by prefix, so collecting a "
            f"sibling af3conf leaf whose label is a prefix of this one (or the "
            f"unlabelled leaf) would read this leaf's {leaf}_* columns as its own "
            f"per-chain columns and blank them. Pick a label without those segments."
        )

    column = ctx.args.input_column
    if column == NO_DEFAULT_COLUMN:
        available = ", ".join(str(c) for c in ctx.df.columns if str(c).endswith("_path"))
        raise ValueError(
            f"af3conf has no default input column: pass -i/--input-column with the "
            f"column holding AF3's PREDICTED STRUCTURE (alphafold3_path or "
            f"alphafold3_<label>_path; the confidence JSONs are found beside it). "
            f"Structure columns in table '{ctx.args.table}': {available or '(none)'}."
        )

    # Not ctx.ready: it drops blank-input rows, which must become error rows.
    ready = ctx.df
    status_col = status_column(ctx.out_dir.name)
    if not ctx.args.force and status_col in ready.columns:
        ready = ready[ready[status_col] != "OK"]

    if column not in ctx.df.columns and not ctx.args.allow_lineage:
        available = ", ".join(str(c) for c in ctx.df.columns if str(c).endswith("_path"))
        raise ValueError(
            f"table '{ctx.args.table}' has no column {column!r}. af3conf does NOT "
            f"walk the lineage by default: an AF3 prediction belongs to one "
            f"sequence, so reading an ancestor's prediction would score the PARENT's "
            f"structure under each child's name. Structure columns here: "
            f"{available or '(none)'}. Pass --allow-lineage if the rows of this "
            f"table really are the molecules that were predicted."
        )

    seq_column = ctx.args.sequence_column
    if (
        seq_column is not None
        and seq_column not in ctx.df.columns
        and not ctx.args.allow_lineage
    ):
        raise ValueError(
            f"table '{ctx.args.table}' has no --sequence-column {seq_column!r}, and "
            f"af3conf does not walk the lineage by default: a parent's sequence is "
            f"not this design's. Pass --allow-lineage only if the rows of this table "
            f"really are the molecules that were predicted."
        )

    members: list[tuple[str, str, str, str]] = []
    n_unresolved = 0
    for name in ready.index:
        name = cast(str, name)
        # A stale result from an earlier run must not survive this submit.
        (ctx.out_dir / f"{name}.tsv").unlink(missing_ok=True)
        raw, source = _resolve(ctx, name, column, ctx.args.allow_lineage)
        sequence = ""
        if seq_column is not None:
            # epitope's normalisation; a blank stays blank (the worker errors on it).
            seq_raw, _ = _resolve(ctx, name, seq_column, ctx.args.allow_lineage)
            sequence = (seq_raw or "").replace(" ", "").upper()
        if raw is None:
            n_unresolved += 1
            print(
                f"{name}: could not resolve {column} (blank own cell, or absent from "
                f"this table and every ancestor) -- submitting anyway so the task "
                f"records an error"
            )
            members.append((name, "", source, sequence))
            continue
        # As-is: the file must stay beside AF3's sibling JSONs.
        members.append((name, str(volume_path(Path(raw))), source, sequence))

    if n_unresolved:
        print(
            f"af3conf: {n_unresolved}/{len(members)} design(s) could not resolve "
            f"{column}. Their tasks write 'error: could not resolve ...', but a row "
            f"whose OWN {column} cell is blank is skipped by collect (filter_ready) "
            f"and KEEPS whatever af3conf status it already had, including a stale "
            f"OK -- gate filters on alphafold3_status == 'OK'."
        )
    if members and n_unresolved == len(members):
        raise ValueError(
            f"af3conf resolved no structure at all for any of the {len(members)} "
            f"ready designs via column {column!r} (searched table "
            f"'{ctx.args.table}' and its ancestors). Check -i/--input-column."
        )

    tasks_dir = ctx.out_dir / "af3conf_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    per_task = ctx.args.designs_per_task
    manifest_rows: list[tuple[str, ...]] = []
    for t, i in enumerate(range(0, len(members), per_task)):
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            for name, src, source, sequence in members[i : i + per_task]:
                # The sequence rides LAST: it is the only field that may be empty.
                f.write(f"{name}\t{src}\t{source}\t{sequence}\n")
        manifest_rows.append(
            (
                str(volume_path(task_file)),
                ",".join(binder),
                ",".join(targets),
                f"{ctx.args.ipsae_pae_cutoff}",
                column,
                # Empty = no sequence check. Last, so cut -f6 reads '' when absent.
                seq_column or "",
            )
        )
    return manifest_rows
