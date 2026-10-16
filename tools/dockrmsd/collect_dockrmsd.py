#!/usr/bin/env python3
"""
Collect dockrmsd results into the table.

Reads ``<run_dir>/<table>/dockrmsd[_<label>]/<prefix>/<name>.tsv`` (written by the
task, or at submit time for a design that could not be submitted) and writes
``dockrmsd[_<label>]_<prefix>_<column>``. Like usalign and posrmsd the status is
PER-COMPARISON (``..._<prefix>_status``), carried in ``data`` with ``status=None``
so the framework's bare ``<leaf>_status`` is not stamped.

STATUS VALUES
    OK                 the target superposition converged within
                       --max-target-fit-rmsd and every number below is a docking
                       number
    warn: ...          numbers computed, but the target fit exceeded
                       --max-target-fit-rmsd, so the binder was moved into a frame
                       that is not the reference frame. Excluded by == "OK"
    error: ...         nothing computed (binder lengths unequal, a chain flag that
                       names no protein chain, no usable permutation, ...). Trust
                       columns already established are KEPT so the row shows why
    missing_input: ... not submitted: --col-b resolved to nothing on this row and
                       all its ancestors
    missing            no result file at all (not selected by --col-a, not yet
                       finished, or killed before it wrote)

TRUST COLUMNS -- read these before `rmsd`:
    binder_chain, ref_binder_chain, target_chains, ref_target_chains, fit_atoms
    target_match          'number' (author residue numbers) or 'sequence' (gapped
                          alignment) -- how the TARGET was paired. The BINDER is
                          always paired by file-order index, whatever this says
    pin_mapping           1 = no permutation search, the given chain order only
    model_chain_ids, ref_chain_ids
    binder_n_res, ref_binder_n_res    equal, or the row errors
    binder_seq_match_frac             NOT a gate: the two binders are the same
                                      backbone with different sequences
    target_fit_rmsd                   the winning permutation's target fit
    best_target_fit_rmsd              the best fit any permutation achieved
    n_target_res_used                 residue pairs the fit used
    n_target_res_used_per_chain       how they split per chain pair ('139/139/139')
                                      -- the floor is applied PER PAIR
    n_target_res_model, n_target_res_ref, target_res_used_frac,
    target_seq_match_frac
    target_resnum_range_model, target_resnum_range_ref
                          the author numbering each file spans ('1-149;1-149;
                          1-149' vs '1-139;1-139;1-139'): numbering provenance as
                          an observation, not an inference from what matched
    n_fit_atoms, n_binder_atoms
    n_permutations, n_candidates, permutation, margin, rmsd_all
SCIENCE:
    rmsd              the docking RMSD: minimum over candidate permutations
    max_dev           largest single binder-atom deviation under the winner
    complex_rmsd      EXACTLY the atom-count-weighted quadrature mean of
                      target_fit_rmsd and rmsd, so it carries no independent
                      information and is ALWAYS SMALLER than rmsd. Not "dominated
                      by the target" -- the terms are quadratic, so past ~2 A the
                      binder term leads. Never a substitute for rmsd
    sup_path          the model moved into the reference frame (mmCIF)

NA (empty), never 0, marks a value that was not measured.

Usage:
    sapia collect dockrmsd outputs/RUN -t table2 \\
        --col-a boltz_path --col-b rfdiffusion3_path
"""

import re
from argparse import ArgumentParser
from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectArgs, CollectCtx, CollectEach, Collected, DesignCtx

FLOAT_COLUMNS = [
    "binder_seq_match_frac",
    "target_fit_rmsd",
    "best_target_fit_rmsd",
    "target_res_used_frac",
    "target_seq_match_frac",
    "margin",
    "rmsd",
    "max_dev",
    "complex_rmsd",
]
INT_COLUMNS = [
    "pin_mapping",
    "binder_n_res",
    "ref_binder_n_res",
    "n_target_res_used",
    "n_target_res_model",
    "n_target_res_ref",
    "n_fit_atoms",
    "n_binder_atoms",
    "n_permutations",
    "n_candidates",
]
STR_COLUMNS = [
    "binder_chain",
    "ref_binder_chain",
    "target_chains",
    "ref_target_chains",
    "fit_atoms",
    "target_match",
    "model_chain_ids",
    "ref_chain_ids",
    "n_target_res_used_per_chain",
    "target_resnum_range_model",
    "target_resnum_range_ref",
    "permutation",
    "rmsd_all",
]
METRIC_COLUMNS = FLOAT_COLUMNS + INT_COLUMNS + STR_COLUMNS
# Never present on a non-OK/non-warn row: these are the numbers the measurement
# produces, and NA (not 0) is what "not measured" looks like.
SCIENCE_COLUMNS = ["rmsd", "max_dev", "complex_rmsd"]
PREFIX_RE = re.compile(r"^[A-Za-z0-9_]+$")


class DockRmsdCollectArgs(CollectArgs):
    col_a: str | None
    col_b: str | None
    output_prefix: str | None


def add_collect_dockrmsd_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--col-a",
        type=str,
        default=None,
        help="Same --col-a as the run (used to derive the prefix when "
        "--output-prefix is not set).",
    )
    parser.add_argument(
        "--col-b",
        type=str,
        default=None,
        help="Same --col-b as the run (used to derive the prefix).",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default=None,
        help="Same --output-prefix as the run. Defaults to "
        "'<col_a minus _path>_vs_<col_b minus _path>'.",
    )


def resolve_prefix(col_a: str | None, col_b: str | None, output_prefix: str | None) -> str:
    """Identical to run_dockrmsd.resolve_prefix (kept in step by a test)."""
    if output_prefix:
        prefix = output_prefix.strip()
    elif col_a and col_b:
        prefix = f"{col_a.strip().removesuffix('_path')}_vs_{col_b.strip().removesuffix('_path')}"
    else:
        raise ValueError("Provide --output-prefix, or both --col-a and --col-b.")
    if not PREFIX_RE.match(prefix):
        raise ValueError(
            f"comparison prefix {prefix!r} must be letters, digits and '_' only "
            f"(it names a directory and every column); pass --output-prefix"
        )
    return prefix


def _present(row: pd.Series, col: str) -> str | None:
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return str(value).strip()


def collect_dockrmsd(ctx: CollectCtx) -> CollectEach:
    """Per-design dockrmsd collector for one named comparison (prefix)."""
    prefix = resolve_prefix(
        getattr(ctx.args, "col_a", None),
        getattr(ctx.args, "col_b", None),
        getattr(ctx.args, "output_prefix", None),
    )
    results_dir = ctx.out_dir / prefix
    if not results_dir.is_dir():
        raise FileNotFoundError(
            f"Results dir not found: {results_dir}. collect needs the same "
            f"--col-a/--col-b or --output-prefix (and -l) as the run."
        )

    def key(col: str) -> str:
        return f"{prefix}_{col}"

    def empty() -> dict[str, Any]:
        data: dict[str, Any] = {key(c): None for c in METRIC_COLUMNS}
        data[key("sup_path")] = None
        return data

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = results_dir / f"{d.name}.tsv"
        data = empty()
        if not tsv_path.is_file():
            data[key("status")] = "missing"
            yield Collected(status=None, data=data)
            return

        # Everything as text: no NA-sniffing of a chain id that spells 'NA'.
        result_df = pd.read_csv(
            tsv_path, sep="\t", dtype=str, keep_default_na=False, na_values=[""]
        )
        if result_df.empty:
            data[key("status")] = "error: empty tsv"
            yield Collected(status=None, data=data)
            return

        row = result_df.iloc[0]
        status = str(row["status"])
        data[key("status")] = status
        measured = status == "OK" or status.startswith("warn:")

        for col in FLOAT_COLUMNS:
            if (v := _present(row, col)) is not None:
                data[key(col)] = float(v)
        for col in INT_COLUMNS:
            if (v := _present(row, col)) is not None:
                data[key(col)] = int(float(v))
        for col in STR_COLUMNS:
            if (v := _present(row, col)) is not None:
                data[key(col)] = v
        if measured:
            data[key("sup_path")] = _present(row, "sup_path")
        else:
            # An error row keeps the trust columns it managed to establish (chain
            # ids, binder lengths, the permutation audit) -- they are the
            # diagnosis -- but NEVER a science number.
            for col in SCIENCE_COLUMNS:
                data[key(col)] = None

        if measured:
            missing = [c for c in SCIENCE_COLUMNS if data[key(c)] is None]
            if missing:
                data = empty()
                data[key("status")] = (
                    f"error: status {status!r} but {', '.join(missing)} absent from "
                    f"{tsv_path.name}"
                )
        yield Collected(status=None, data=data)

    return one
