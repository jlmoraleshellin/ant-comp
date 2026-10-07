#!/usr/bin/env python3
"""
Collect polyv results into the table.

Reads the per-design <name>.tsv files that polyv_worker.py wrote under
<run_dir>/<table>/polyv[_<label>]/ and merges them back as <prefix>_status /
<prefix>_path (the packed poly-valine pose) plus one <prefix>_<field> column per
metric:

    MEASUREMENTS
    cms_target       ContactMolecularSurface of the TARGET surface weighted toward
                     the polyV binder, A^2
    cms_binder       the same with the sides swapped (binder surface), A^2
    dsasa            interface dSASA summed over BOTH partners, A^2
                     (InterfaceAnalyzerMover). NOTE this is about TWICE the
                     single-side buried surface area usually quoted in the binder
                     literature -- do not compare it with a one-sided BSA.
    dsasa_per_res    dsasa / nres_int; NA when nres_int is 0
    nres_int         interface residues, both sides
    nres_int_binder  interface residues on the binder
    nres_int_target  interface residues on the target
    sc               interface shape complementarity, 0-1; NA when Rosetta's SC
                     calculation failed (it can, on a tiny interface)

    TRUST METRICS -- read these before any of the above
    seed             the RNG seed used, derived from the design name alone. The
                     packer is stochastic and Rosetta's RNG state carries across the
                     designs of a task, so without a per-design seed a score would
                     depend on batch position. Two rows with the same seed and the
                     same input must agree to every decimal.
    n_mutated        binder residues that are VAL in the packed pose
    n_kept_gly_pro   binder residues left as GLY/PRO (--keep-gly-pro only; 0 otherwise)
    n_res_binder     binder residues SCORED (i.e. after --exclude-resnames)
    n_res_target     target residues SCORED (i.e. after --exclude-resnames)
    n_res_dropped_chain    residues deleted because their chain was in neither
                     selection
    n_res_dropped_resname  residues deleted because --exclude-resnames named them.
                     Kept separate from the chain count so that "a whole chain went"
                     and "the exclusion flag bit" are distinguishable
    chains_found     every chain ID present in the pose as loaded, comma-joined
    n_input_restypes distinct residue NAMES on the binder in the input. 1 == a uniform
                     backbone whatever letter the generator used; >5 == the pose
                     arrived carrying a real sequence, which polyv then erased
    seconds          wall time of this design

NA (an empty cell) rather than 0 marks a value that could not be computed. A dSASA,
CMS or nres_int of 0 is a real measurement -- the binder does not touch the target --
and is NOT a failure; the status stays OK.

Usage:
    sapia collect polyv outputs/20260928_100140_7ojg_binder --table table0 -l bb
"""

from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Bare column names; the driver leaf-prefixes them (polyv_<name>).
FLOAT_COLUMNS = [
    "cms_target",
    "cms_binder",
    "dsasa",
    "dsasa_per_res",
    "sc",
    "seconds",
]
INT_COLUMNS = [
    "nres_int",
    "nres_int_binder",
    "nres_int_target",
    "n_mutated",
    "n_kept_gly_pro",
    "n_res_binder",
    "n_res_target",
    "n_res_dropped_chain",
    "n_res_dropped_resname",
    "n_input_restypes",
    "seed",
]
STR_COLUMNS = ["chains_found"]
RESULT_COLUMNS = FLOAT_COLUMNS + INT_COLUMNS + STR_COLUMNS


def _empty() -> dict[str, Any]:
    return {c: None for c in RESULT_COLUMNS}


def _present(row: pd.Series, col: str) -> Any | None:
    """The cell's value, or None when the column is absent, NaN or blank."""
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def _fields(row: pd.Series) -> dict[str, Any]:
    data = _empty()
    for col in FLOAT_COLUMNS:
        if (value := _present(row, col)) is not None:
            data[col] = float(value)
    for col in INT_COLUMNS:
        if (value := _present(row, col)) is not None:
            data[col] = int(float(value))
    for col in STR_COLUMNS:
        if (value := _present(row, col)) is not None:
            data[col] = str(value)
    return data


def collect_polyv(ctx: CollectCtx) -> CollectEach:
    """Per-design polyv collector. Variants -- another structure column or another
    chain split off the same table -- are distinguished via --dir-label, matching the
    output dir."""

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"
        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty())
            return

        result_df = pd.read_csv(
            tsv_path, sep="\t", dtype={c: str for c in STR_COLUMNS}
        )
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty())
            return

        row = result_df.iloc[0]
        status = str(row["status"])
        if status != "OK":
            yield Collected(status=status, path="", data=_empty())
            return

        # Rebuild the path from out_dir rather than trusting the worker's absolute
        # one, so it is stored run-relative like every other tool's <leaf>_path.
        packed = ctx.out_dir / f"{d.name}.pdb"
        yield Collected(
            status=status,
            path=str(packed) if packed.is_file() else "",
            data=_fields(row),
        )

    return one
