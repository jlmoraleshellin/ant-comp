#!/usr/bin/env python3
"""
Collect epitope results into the table.

Reads the per-design <name>.tsv files that epitope_worker.py wrote under
<run_dir>/<table>/epitope[_<label>]/ and merges them back as <prefix>_status /
<prefix>_path (the per-target-residue TSV: chain, resnum, resname, n_contacts,
min_dist, is_hotspot -- EVERY target residue, contacted or not) plus one
<prefix>_<field> column per metric:

  DID IT LAND ON THE EPITOPE I CHOSE:
    hotspot_recall        contacted hotspots / listed hotspots, 0.0-1.0. This is
                          the filter column. 0.0 is a REAL value (the binder bound
                          somewhere else); NA means nothing was measured.
    hotspot_hits          the hotspots contacted, comma-joined, or 'none'
    n_hotspots            how many were listed -- the denominator
    min_dist_hotspot      closest heavy-atom approach to ANY listed hotspot, A.
                          The graded signal when recall is 0: 6 A is a near miss,
                          40 A is the other side of the target.

  THE TRUST COLUMN -- read it before believing any of the above:
    hotspot_resnames      the residue NAME found at each listed hotspot position,
                          in the order listed: 'A96=LYS,A99=GLU,...'. If those are
                          not the residues you chose, the structure was renumbered
                          and every other column here answers a different question.
                          (A hotspot number that is ABSENT from the target is an
                          error status, not a recall of 0 -- this column is for the
                          nastier case where the number exists but means something
                          else.)

  WHICH EPITOPE IT ACTUALLY USED:
    n_iface_target_res    target residues with >= 1 contact
    iface_target_resnums  which ones, comma-joined ('A94,A96,A97,...'), or 'none'.
                          Compare two designs' lists to see whether a confident
                          binder found a different site.
    n_contacts            total heavy-atom contact pairs across the interface. A
                          crude size proxy only -- for interface AREA use cms
                          (sc, sc_area) or pyrosetta (if_dSASA).

  HOW CENTRAL THE LISTED HOTSPOTS ARE TO THAT INTERFACE -- the same matrix again,
  restricted to the listed residues (added 2026-10-02). hotspot_recall says
  'touched'; these say 'how much of the binding is there':
    hotspot_n_contacts    contact pairs on the LISTED hotspots only. Same cutoff,
                          same pair counting as n_contacts, fewer residues, so
                          hotspot_n_contacts <= n_contacts always.
    hotspot_contact_frac  hotspot_n_contacts / n_contacts, 0.0-1.0: the SHARE of
                          this design's interface carried by the listed residues.
                          THE FILTER COLUMN when the question is 'grazed it' vs
                          'built around it'. Read it as a share, not an area: an
                          interface of 300 pairs and one of 80 pairs can both sit
                          at 0.15. NA when n_contacts is 0 (the share is 0/0,
                          undefined) -- while hotspot_n_contacts is a real 0 there.
                          It counts ATOM pairs, so a big, deeply inserted side
                          chain scores higher than a small one at equal burial.
    hotspot_contacts      the per-hotspot pair counts in the order --hotspots was
                          given, 'A24=0,A36=3,A99=21'. Formatted like
                          hotspot_resnames and in the same order, so the two read
                          against each other. An explicit 0 is a measured zero.

  WHICH OF THE BINDER'S OWN RESIDUES DO THE BINDING -- the same matrix, read along
  its other axis (added 2026-09-30):
    n_iface_design_res    design residues with >= 1 contact
    iface_design_resnums  which ones, comma-joined ('B12,B15,B19,...') in ASCENDING
                          residue order, or 'none'
    iface_design_resnames THE BINDER-SIDE TRUST COLUMN: the same residues with the
                          identity at each one, 'B12=HIS,B15=TYR,B19=GLU,...',
                          formatted exactly like hotspot_resnames and in the same
                          order as iface_design_resnums. Eyeball it: a frame shift
                          between a designed sequence and its backbone is obvious
                          here and invisible everywhere else.

  WHERE THE DESIGN IDENTITIES CAME FROM -- read before believing the names above:
    seq_source            'structure' (identities from the coordinates -- the
                          default and the original behaviour) or the sequence
                          column they were taken from
    seq_len               length of the sequence actually used
    n_design_res_struct   amino-acid residues in the design chain(s) of the
                          STRUCTURE. A disagreement with seq_len is an 'error:',
                          never a truncated or padded mapping.

  STERIC OVERLAP (added 2026-10-08) -- the same heavy atoms, read for overlap:
    min_dist_any          smallest heavy-atom distance between ANY design atom and
                          ANY target atom, A (<= min_dist_hotspot always)
    n_clash               design/target heavy-atom PAIRS closer than clash_cutoff.
                          0 is a real zero on an OK row; NA on every error row.
                          Unit = atom PAIRS: never comparable with ringfit_n_clash,
                          which counts binder ATOMS.
    n_clash_design_res    distinct design residues in >= 1 clash pair (the same
                          quantity as ringfit_n_clash_res)
    clash_cutoff          the A cutoff used (--clash-cutoff, default 2.5)
    min_dist_bb           as min_dist_any with the DESIGN side restricted to its
                          backbone (N, CA, C, O, OXT); >= min_dist_any
    n_clash_bb            as n_clash over those pairs; <= n_clash. Both NA when the
                          design chains have no backbone-named atoms. The overlap
                          sequence design cannot remove. CB NOT included: the
                          FLOOR (lenient reading).
    min_dist_bbcb         as min_dist_bb with CB added (CB persists under every
    n_clash_bbcb          sequence except GLY): the honest persistent set and THE
                          GATE for a generator backbone. any <= bbcb <= bb for the
                          distance, bb <= bbcb <= any for the count.
    Geometry, not an energy; H-bond / salt-bridge distances (2.5-3.0 A) are not
    clashes. On a backbone-only binder there are no side chains to clash.

  THE ECHO OF WHAT WAS ASKED:
    design_chains         the binder chain(s) actually scored
    target_chains         the target chain(s) actually scored
    contact_cutoff        the A cutoff used -- a recall is only comparable with
                          another recall at the same cutoff

NA (empty) rather than 0 marks a metric that could not be measured: the design
failed, a chain was absent, a hotspot did not exist, or no TSV was written. A
hotspot_recall of 0.0 with status OK is the opposite -- a successful measurement of
a binder that ignored the epitope -- and the two must never be confused, which is
why nothing here defaults to zero.

Statuses: 'OK'; 'error: ...' for an absent design or target chain, a hotspot absent
from the target chain (the numbering contract), an unreadable structure, a structure
path that could not be resolved from the input column (on the row or up the
lineage), a --sequence-column whose length disagrees with the structure or which
holds a multi-chain '/' string, or a crashed worker; 'missing' when the design has
no TSV at all.

BACKWARD COMPATIBILITY. The columns added on 2026-09-30, 2026-10-02 and 2026-10-08
are appended, never substituted: every pre-existing column keeps its name, type,
formatting and meaning, and with --sequence-column omitted every pre-existing
number is bit-for-bit what it was (--clash-cutoff feeds only the 2026-10-08
columns). Rows collected before those dates simply have the newer columns empty --
a missing column in the per-design TSV reads back as NA, like any other unmeasured
metric. To backfill, re-run with --force under the same -l label AND collect with
--force (without it, already-OK rows are skipped and keep the old, empty cells).

Usage:
    sapia collect epitope outputs/20260930_egfr_binder --table table1 -l pred
"""

from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Bare column names; the driver leaf-prefixes them (epitope_<name>).
#
# APPEND-ONLY. Live tables carry epitope_bb_*, epitope_bbnterm_* and epitope_pred_*
# from earlier runs; no name here may be removed, renamed or change type.
FLOAT_COLUMNS = [
    "hotspot_recall",
    "min_dist_hotspot",
    "contact_cutoff",
    # --- added 2026-10-02 ---
    "hotspot_contact_frac",
]
INT_COLUMNS = [
    "n_hotspots",
    "n_iface_target_res",
    "n_contacts",
    # --- added 2026-09-30 ---
    "n_iface_design_res",
    "seq_len",
    "n_design_res_struct",
    # --- added 2026-10-02 ---
    "hotspot_n_contacts",
]
STR_COLUMNS = [
    "hotspot_hits",
    "hotspot_resnames",
    "iface_target_resnums",
    "design_chains",
    "target_chains",
    # --- added 2026-09-30 ---
    "iface_design_resnums",
    "iface_design_resnames",
    "seq_source",
    # --- added 2026-10-02 ---
    "hotspot_contacts",
]
# --- added 2026-10-08: steric overlap. Kept in their OWN typed list and appended
# at the END of RESULT_COLUMNS, so every pre-existing column keeps its position in a
# freshly collected table as well as its name and type.
CLASH_COLUMNS: dict[str, type] = {
    "min_dist_any": float,
    "n_clash": int,
    "n_clash_design_res": int,
    "clash_cutoff": float,
    # second 2026-10-08 pass: design BACKBONE vs every target heavy atom
    "min_dist_bb": float,
    "n_clash_bb": int,
    # third 2026-10-08 pass: design backbone PLUS CB vs every target heavy atom
    "min_dist_bbcb": float,
    "n_clash_bbcb": int,
}
RESULT_COLUMNS = FLOAT_COLUMNS + INT_COLUMNS + STR_COLUMNS + list(CLASH_COLUMNS)


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
    metric that could not be measured, not a zero)."""
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
    for col, kind in CLASH_COLUMNS.items():
        if (value := _present(row, col)) is not None:
            data[col] = int(float(value)) if kind is int else float(value)
    return data


def collect_epitope(ctx: CollectCtx) -> CollectEach:
    """Per-design epitope collector. Variants -- the designed backbone and its
    prediction, or two different epitopes scored on one table -- are distinguished
    via --dir-label, matching the output dir."""

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"
        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty())
            return

        result_df = pd.read_csv(
            tsv_path,
            sep="\t",
            dtype={c: str for c in STR_COLUMNS},
        )
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty())
            return

        row = result_df.iloc[0]
        status = str(row["status"])
        if status != "OK":
            # An absent chain, a hotspot the target does not have, a crashed
            # worker: no path, no numbers. A recall of 0.0 here would be read as
            # "bound elsewhere" -- the single most expensive confusion this tool
            # exists to prevent.
            yield Collected(status=status, path="", data=_empty())
            return

        path = _present(row, "path")
        yield Collected(
            status=status,
            path="" if path is None else str(path),
            data=_fields(row),
        )

    return one
