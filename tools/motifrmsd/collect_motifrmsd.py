#!/usr/bin/env python3
"""
Collect motifrmsd results into the table.

Scans ``<run_dir>/<table>/motifrmsd/`` for the per-design TSV files written by
motifrmsd.sh and merges them back as ``<prefix>_status`` / ``<prefix>_path`` (the
design superposed into the reference frame) plus one ``<prefix>_<metric>`` column
per metric.

TRUST COLUMNS -- read these before any RMSD:

    n_expected          residues listed in rfd3's diffused_index_map
    n_matched           pairs whose CA was found in BOTH structures. Equal to
                        n_expected on every OK row -- a shortfall is an error by
                        default, and a 'warn:' row only with --allow-partial
    seq_match_frac      fraction of matched pairs whose residue identity agrees
                        between reference and design. NOT gated on: rfd3's output
                        residue naming convention for scaffolded backbones is not
                        established in this repo, so a uniformly low value across
                        every design says something about rfd3, while a value that
                        varies per design says the mapping is wrong for those rows
    ref_source          'flag' (--ref-structure) or 'sidecar'
                        (specification.input)
    ref_matches_spec    1 = the reference measured against is byte-identical to
                        the one rfd3 recorded diffusing; 0 = a different file (the
                        status is 'warn:'); NA = the sidecar's input was not
                        readable, so the question is open -- never read NA as 1
    ref_used/spec_input the file actually read, and the one rfd3 recorded
    ref_chains          reference chains the map addressed, e.g. 'D'
    design_chains       design chains the map addressed, e.g. 'A'
    ref_motif_seq       one-letter motif sequence of the reference, segment order
    design_motif_seq    the same positions in the design -- the two together are
                        what makes a low seq_match_frac diagnosable
    sampled_contig      the resolved contig rfd3 actually built (per design)

SCIENCE:

    rmsd                CA RMSD (A) of the motif after a Kabsch fit ON THE MOTIF
    max_dev             largest single-CA deviation after that same fit -- an
                        RMSD hides one badly placed residue, this does not
    n_segments          contiguous runs of the motif in REFERENCE numbering
    worst_seg_rmsd      max over the per-segment RMSDs; the filter-safe scalar,
                        present on every OK row whatever the segment count
    seg_ranges          the segments, e.g. 'D24-37/D48-67', in reference order
    seg_rmsds           their RMSDs, '/'-joined, same order
    seg0_rmsd, seg1_rmsd, ...
                        one filterable column per segment, from the SAME single
                        global superposition -- never a per-segment refit. These
                        columns exist only for segment indices some design in the
                        table actually has.
    design_fixed_positions
                        THE SAME segments in the DESIGN's own output numbering,
                        formatted as a ProteinMPNN --fixed-positions value for a
                        single chain group: '28:35,50:59' -- inclusive start:end,
                        comma-joined, NO trailing slash (a trailing '/' would
                        declare a second, empty chain group). Fragment i is
                        segment i, same order as seg_ranges/seg_rmsds. Feed it as
                        --chains-to-design <design chain> --fixed-positions
                        '<this>' to redesign everything EXCEPT the graft. NA --
                        never a partial string -- when the motif spans more than
                        one design chain, or when that chain is not numbered
                        1..N in file order (ProteinMPNN positions are 1-indexed
                        CHAIN POSITIONS, not PDB residue numbers; they coincide
                        only for a gapless chain numbered from 1, which is what
                        rfd3 writes). Read fixed_positions_chain alongside it.

                        In a BINDER-SCAFFOLDING contig (a target held fixed plus
                        a motif grafted onto the binder chain) rfd3 maps the
                        whole target too, so the motif always spans several
                        design chains and this column is always NA. That is what
                        the run flag --motif-chain <reference chain letter> is
                        for: it restricts THIS COLUMN to the residues lifted from
                        one reference chain, leaving every other column measuring
                        the full mapping.
    fixed_positions_chain
                        the DESIGN-side chain letter design_fixed_positions
                        addresses -- the value to pass as --chains-to-design.
                        NA exactly when design_fixed_positions is NA. NOT
                        design_chains, which keeps listing every chain the RMSD
                        covers (all four of 'A,B,C,D' in a binder-scaffolding
                        run) and is therefore the wrong column to hand to
                        ProteinMPNN.

NA (empty) rather than 0 marks a metric that was not computed. Every error row
carries NA for every numeric column: a 0.0 RMSD would read as a perfect graft.

Usage:
    sapia collect motifrmsd outputs/20261008_085450_5M2J_vhh2_motif_graft \
        --table table0_seg2
"""

import re
from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Bare column names; the driver leaf-prefixes them (motifrmsd_<name>).
FLOAT_COLUMNS = [
    "seq_match_frac",
    "rmsd",
    "max_dev",
    "worst_seg_rmsd",
]
INT_COLUMNS = [
    "n_expected",
    "n_matched",
    "ref_matches_spec",
    "n_segments",
]
STR_COLUMNS = [
    "ref_source",
    "ref_used",
    "spec_input",
    "ref_chains",
    "design_chains",
    "ref_motif_seq",
    "design_motif_seq",
    "sampled_contig",
    "seg_ranges",
    "seg_rmsds",
    "design_fixed_positions",
    "fixed_positions_chain",
]
METRIC_COLUMNS = FLOAT_COLUMNS + INT_COLUMNS + STR_COLUMNS

# Per-segment columns are per design (a design with two motif segments writes
# seg0_rmsd and seg1_rmsd), so they are discovered rather than declared.
SEG_COLUMN = re.compile(r"^seg\d+_rmsd$")


def _empty() -> dict[str, Any]:
    """Every declared column as NA. Used for every non-OK row, so a failure is
    never indistinguishable from a measurement of zero."""
    return {c: None for c in METRIC_COLUMNS}


def _present(row: pd.Series, col: str) -> Any | None:
    """The cell's value, or None when the column is absent, NaN or blank."""
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def _metrics(row: pd.Series) -> dict[str, Any]:
    """One result row -> collected columns, NA-safe (a missing or blank cell is a
    metric that did not apply, not a zero)."""
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
    for col in row.index:
        if SEG_COLUMN.match(str(col)) and (value := _present(row, col)) is not None:
            data[str(col)] = float(value)
    return data


def collect_motifrmsd(ctx: CollectCtx) -> CollectEach:
    """Per-design motifrmsd collector. The framework iterates ready designs and
    stamps status/path (keyed by the tool leaf); this reads one design's one-row
    ``<name>.tsv`` and adds the motif-fidelity columns. Variants are distinguished
    via --dir-label, matching the output dir."""

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"

        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty())
            return

        result_df = pd.read_csv(tsv_path, sep="\t")
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty())
            return

        row_data = result_df.iloc[0]
        status = str(row_data["status"])

        # 'warn: ...' means the numbers were computed but something about the
        # correspondence or the reference is not what it should be (a partial
        # match under --allow-partial, or a reference that is not the file rfd3
        # diffused). Keep the numbers so the suspicion can be judged from the
        # table, while the status keeps the row out of any `== "OK"` selection.
        if status != "OK" and not status.startswith("warn:"):
            yield Collected(status=status, path="", data=_empty())
            return

        aligned_path = row_data.get("aligned_path")
        yield Collected(
            status=status,
            path=""
            if aligned_path is None or pd.isna(aligned_path)
            else str(aligned_path),
            data=_metrics(row_data),
        )

    return one
