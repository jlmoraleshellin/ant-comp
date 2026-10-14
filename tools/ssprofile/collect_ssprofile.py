#!/usr/bin/env python3
"""
Collect ssprofile results into the table.

Reads the per-design <name>.tsv files that ssprofile_worker.py wrote under
<run_dir>/<table>/ssprofile[_<label>]/ and merges them back as <prefix>_status /
<prefix>_path (the per-residue TSV: chain, resnum, resname, ss) plus one
<prefix>_<field> column per metric.

  TRUST COLUMNS -- read these first:
    chains_used         the chain(s) actually profiled, comma-joined
    sequence            one-letter sequence of the selection, chains joined by '/'
                        (modified residues lower-case, unknown 'X'). The wrong-chain
                        guard when binder and target have similar lengths.
    input_column        the table column the structure came from  } provenance, kept
    source_path         the file actually read (post CIF->PDB)   } on error rows too
    n_res               amino-acid residues with a CA in that selection: the CA
                        count behind Rg and the denominator of every fraction.
                        Compare it with the length you expect the chain to have.
    n_chain_breaks      C(i)-N(i+1) gaps > 2.5 A inside a selected chain. 0 on every
                        OK row; > 0 only on a 'warn:' row. NA when the backbone was
                        incomplete (the gap test needs N and C).

  SCIENCE COLUMNS:
    radius_of_gyration  A. CA-based, MASS-UNWEIGHTED: sqrt(mean(|r_i - centroid|^2))
                        over the selection's CA atoms, centroid = their plain mean.
    helix_fraction      DSSP H (incl. 3-10 G and pi I) / n_res
    beta_fraction       DSSP E (incl. isolated bridge B) / n_res
    loop_fraction       everything else / n_res
                        The three sum to 1.0 on an OK row and are NA together
                        otherwise.

Statuses:
    'OK'         every column measured.
    'warn: ...'  Rg / n_res / chains_used measured, the three fractions NA (missing
                 backbone atoms, a chain break, or a failed DSSP call).
    'error: ...' every measured column NA; input_column/source_path kept when the
                 worker wrote them (absent chain, no CA in a selected chain,
                 unreadable or unresolvable structure; a crashed worker writes none).
    'missing'    no TSV at all.

NA (an empty cell) is never 0: a 0.0 beta_fraction on an OK row is a measured
all-helix chain; an empty one was not measured.

Usage:
    sapia collect ssprofile outputs/20261009_binders --table table0 -l bb
"""

from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Bare column names; the driver leaf-prefixes them (ssprofile_<name>).
FLOAT_COLUMNS = [
    "radius_of_gyration",
    "helix_fraction",
    "beta_fraction",
    "loop_fraction",
]
INT_COLUMNS = ["n_res", "n_chain_breaks"]
STR_COLUMNS = ["chains_used", "sequence", "input_column", "source_path"]
RESULT_COLUMNS = FLOAT_COLUMNS + INT_COLUMNS + STR_COLUMNS
# Provenance: what the row was measured FROM. Kept on error rows too, because "which
# file failed" is exactly what an error row needs to say.
PROVENANCE = ["input_column", "source_path"]
FRACTIONS = ["helix_fraction", "beta_fraction", "loop_fraction"]


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


def _provenance(row: pd.Series) -> dict[str, Any]:
    data = _empty()
    for col in PROVENANCE:
        if (value := _present(row, col)) is not None:
            data[col] = str(value)
    return data


def _check(status: str, data: dict[str, Any]) -> str:
    """Re-assert the worker's contract at collect time, so a TSV that violates it
    (hand-edited, or written by a future buggy worker) lands as an error rather than
    as plausible numbers. Returns the status to record."""
    fracs = [data[c] for c in FRACTIONS]
    if status == "OK":
        if any(v is None for v in fracs) or data["radius_of_gyration"] is None:
            return "error: OK row with missing metrics in the result TSV"
        if data["n_res"] is None or data["n_res"] < 1:
            return "error: OK row with no residues in the result TSV"
        if abs(sum(fracs) - 1.0) > 1e-3:
            return f"error: SS fractions sum to {sum(fracs):.6f}, not 1.0"
    elif status.startswith("warn"):
        if any(v is not None for v in fracs):
            return "error: warn row carries SS fractions in the result TSV"
        if data["radius_of_gyration"] is None:
            return "error: warn row without radius_of_gyration in the result TSV"
    return status


def collect_ssprofile(ctx: CollectCtx) -> CollectEach:
    """Per-design ssprofile collector. Variants (the designed backbone and its
    prediction, or two different chain selections) are distinguished via
    --dir-label, matching the output dir."""

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"
        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty())
            return

        result_df = pd.read_csv(tsv_path, sep="\t", dtype={c: str for c in STR_COLUMNS})
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty())
            return

        row = result_df.iloc[0]
        status = str(row["status"])
        if status.startswith("error") or not (status == "OK" or status.startswith("warn")):
            # Absent chain, no CA, unreadable structure, crashed worker: no numbers,
            # but the provenance (when the worker got far enough to write it) stays.
            yield Collected(status=status, path="", data=_provenance(row))
            return

        data = _fields(row)
        checked = _check(status, data)
        if checked.startswith("error"):
            yield Collected(status=checked, path="", data=_provenance(row))
            return

        path = _present(row, "path")
        yield Collected(status=checked, path="" if path is None else str(path), data=data)

    return one
