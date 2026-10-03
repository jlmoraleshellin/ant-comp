#!/usr/bin/env python3
"""
Collect dimerfit results into the dock table.

Reads the per-design <name>.tsv files that dimerfit_worker.py wrote under
<run_dir>/<table>/dimerfit[_<label>]/ and merges them back as <prefix>_status /
<prefix>_path (the TRANSFORMED DIMER, no target -- the file the next step consumes)
plus one <prefix>_<field> column per metric:

    complex_path        the transformed dimer PLUS the reference target, for
                        inspection only

  mapping trust -- a wrong superposition yields plausible geometry, so these are
  the point of the tool, not an extra:
    align_rmsd          CA RMSD (A) of dock protomer A onto the reference binder
    n_align_atoms       how many CA pairs the fit used
    seq_match_frac      fraction of those pairs whose residue identity agrees
                        (expect 1.0: the dock IS the reference binder; below 0.9
                        the status becomes 'warn: ...' and the metrics still land)
    resnum_offset       modal ref_resnum - dock_resnum (expect 0: rpxdock renumbers
                        1..N and the reference binder is already numbered 1..N, so
                        anything else means a renumbering nobody accounted for ->
                        'warn: resnum_offset N')
    resnum_match        'ordinal' or 'resnum' -- how the pairs were formed
    dock_chains         what was actually used; 'auto' is resolved here. The first
    ref_binder_chain    dock chain is protomer A, the second protomer B
    ref_target_chains
    complex_target_chains  chain IDs the target got inside complex_path (renamed
                        only when they collided with a protomer's)
    n_protomers         2 -- any other count is an 'error:' row (the count is in
                        the error message, since every metric is then NA)
    n_res_a, n_res_b    amino acids per protomer; check them against the parent's
    ref_binder_len      ifacegeom_binder_len -- an invariant this tool does NOT
                        compute, so it is independent evidence
    n_atoms_b           protomer B heavy atoms = the denominator of clash_frac

  C2 interface geometry, relative to the epitope (we want overlap LOW):
    dimer_iface_res     the C2 interface residues on the protomer-A side, as
    n_dimer_iface_res   'A:12,A:15,...'
    dimer_iface_com     mass-weighted centres as 'x,y,z', in the REFERENCE frame
    epitope_com
    epitope_com_dist    |dimer_iface_com - epitope_com| (A)
    n_epitope_res       epitope residues mapped onto the dock = the denominator of
                        frac_overlap
    overlap_res         residues that are BOTH C2 interface and epitope, and the
    n_overlap_res       fraction of the epitope they are. LOW is what we want: the
    frac_overlap        dimer interface should sit ADJACENT to the epitope

  occlusion of the target binding site by protomer B (we want this HIGH, i.e.
  opposite-signed to frac_overlap -- do not collapse the two):
    n_clash             protomer B heavy atoms within --clash-cutoff of the target.
    clash_frac          A geometric COUNT, not an energy: not fa_rep
    min_dist_b_target   closest protomer-B-to-target heavy-atom distance (A)
    n_target_epitope_res  target binding-site residues read from
                        --target-epitope-column = the denominator below
    occluded_res        target epitope residues within --occlusion-cutoff of
    n_occluded_res      protomer B, and the fraction of them. HIGH = the partner
    occluded_frac       protomer blocks the target
    seconds             wall time of the design's measurement

NA (empty) rather than 0 marks a field that did not apply: the design failed, no TSV
was written, or the geometry was undefined. A blank must never be read as a passing
zero -- in particular a missing occluded_frac is not "nothing occluded".

A status of 'warn: ...' means the numbers were computed but something about the
mapping or the geometry is suspect (low seq_match_frac, a non-zero resnum_offset, or
no C2 interface at all). The metrics are kept so the suspicion can be judged from the
table, and the 'warn:' prefix keeps the row out of any `== "OK"` selection.

Usage:
    sapia collect dimerfit outputs/20261002_143419_dimer_phase2 --table table1
"""

from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Bare column names; the driver leaf-prefixes them (dimerfit_<name>).
FLOAT_COLUMNS = [
    "align_rmsd",
    "seq_match_frac",
    "epitope_com_dist",
    "frac_overlap",
    "clash_frac",
    "min_dist_b_target",
    "occluded_frac",
    "seconds",
]
INT_COLUMNS = [
    "n_align_atoms",
    "resnum_offset",
    "n_protomers",
    "n_res_a",
    "n_res_b",
    "ref_binder_len",
    "n_atoms_b",
    "n_dimer_iface_res",
    "n_epitope_res",
    "n_overlap_res",
    "n_clash",
    "n_target_epitope_res",
    "n_occluded_res",
]
STR_COLUMNS = [
    "complex_path",
    "resnum_match",
    "dock_chains",
    "ref_binder_chain",
    "ref_target_chains",
    "complex_target_chains",
    "dimer_iface_res",
    "dimer_iface_com",
    "epitope_com",
    "overlap_res",
    "occluded_res",
]
RESULT_COLUMNS = FLOAT_COLUMNS + INT_COLUMNS + STR_COLUMNS

# Every string column is read as a string: a residue list like 'A:12,A:15', a
# coordinate like '1.0,2.0,3.0' and a chain id like '1' must never be guessed at by
# pandas' type inference.
STR_DTYPES = {c: str for c in STR_COLUMNS}


def _empty() -> dict[str, Any]:
    return {c: None for c in RESULT_COLUMNS}


def _present(row: pd.Series, col: str) -> Any | None:
    """The cell's value, or None when the column is absent, NaN or blank."""
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def _fields(row: pd.Series) -> dict[str, Any]:
    """One result row -> collected columns, NA-safe (a missing or blank cell is a
    field that did not apply, not a zero)."""
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


def collect_dimerfit(ctx: CollectCtx) -> CollectEach:
    """Per-design dimerfit collector. The framework iterates ready designs and stamps
    status/path (keyed by the tool leaf); this reads one design's one-row <name>.tsv
    and adds the mapping-trust, geometry and occlusion columns. Variants -- another
    dock column, another reference, other cutoffs -- are distinguished via
    --dir-label, matching the output dir."""

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"
        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty())
            return

        result_df = pd.read_csv(tsv_path, sep="\t", dtype=STR_DTYPES)
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty())
            return

        row = result_df.iloc[0]
        status = str(row["status"])

        # An 'error:' row keeps every metric NA -- the message (which names the
        # reason, e.g. the protomer count) is the only thing worth carrying.
        if status != "OK" and not status.startswith("warn:"):
            yield Collected(status=status, path="", data=_empty())
            return

        path = _present(row, "path")
        yield Collected(
            status=status,
            path="" if path is None else str(path),
            data=_fields(row),
        )

    return one
