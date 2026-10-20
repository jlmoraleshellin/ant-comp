#!/usr/bin/env python3
"""
Collect hbdesigner outputs into the child table it reserved.

Reads the one-row-per-rank TSV ``<run_dir>/<table>/hbdesigner/<design>.tsv`` written
by hbdesigner_worker.py and mints one child row per kept network, keyed
``<parent>_hb<rank>`` (rank 1 = HBDesigner's best) and carrying ``parent_name`` for
lineage. Nothing from the ancestors is copied -- resolve it with ``lookup`` when you
need it.

Columns (leaf-prefixed ``hbdesigner_``). The three TRUST columns come first:
``resnum_shift_max`` (is the mpnn mapping unshifted), ``pdb_parser`` (which reader
produced that mapping) and ``guide_atoms_removed`` (was the --guide_res marker
atom taken out of the file the row points at).

    resnum_shift_max      max |mpnn position - PDB resnum| over the network
                          residues. 0 = no renumbering happened, so `network`'s
                          resnums ARE the mpnn positions. Higher = more shift;
                          any non-zero value means read `resnum_offset` before
                          trusting a hand-off.
    pdb_parser            which reader produced the residue mapping: 'gemmi'
                          (Modal, and anywhere gemmi is importable) or
                          'biopython' (the fallback used in the vib conda env,
                          which has no gemmi). Not a quality score -- the two
                          agree on the standard 20 amino acids -- but they do NOT
                          agree on every residue name (39 names only gemmi calls
                          an amino acid, 942 only Biopython, and PCA is E to one
                          and Q to the other), so n_res_total can be reader-
                          dependent for a file containing one. Gate a
                          cross-backend comparison on this column.
    guide_atoms_removed   HETATM lines removed from this rank's PDB because they
                          were HBDesigner's --guide_res marker atom (resname ORI
                          AND atom name V1, element V). 0 = the file was checked
                          and had none (the normal case without --guide-res);
                          1 = the marker was taken out and the file as HBDesigner
                          wrote it kept under <design>/guide_raw/. Upstream adds
                          exactly one, so a value ABOVE 1 means the output does
                          not match the code this tool was built against -- it is
                          reported, not capped. NA = the strip pass did not run or
                          failed, so the file MAY still carry the marker -- NA is
                          not 0 here. Collect prints a count of the NA rows.

    fold_sequence        THE SEQUENCE A PREDICTOR MAY FOLD, and the reason this
                          tool is worth running on an already-designed binder. It
                          is the STAGED INPUT's binder-chain sequence with the
                          network's binder positions overwritten by the network
                          residues -- it keeps every residue the input already
                          had and changes only the network. `path` cannot be
                          folded (poly-glycine); this can. NA when the binder
                          chain cannot be identified or a check failed -- NEVER
                          an empty string, which would read as a binder with no
                          residues. Its quality is bounded by the staged input:
                          if that was poly-glycine, so is this.
    n_fold_changes       how many of the network's binder positions differ from
                          the staged input's residue there. 0 is legitimate (the
                          network re-chose what was already there); higher means
                          more of the binder moved. Not a quality score.
    fold_chain           the chain fold_sequence belongs to: the single network
                          chain NOT in --omit-chains. NA with the other two.
    centers               the designed HISTIDINES on the binder, in hisbond's
                          centres format: '<index0>:HIS' joined by ';', e.g.
                          '11:HIS;19:HIS'. index0 is 0-BASED INTO fold_sequence
                          (hence into the binder chain in file order), NOT a
                          residue number. The token is always a bare 'HIS':
                          hbdesigner designs the residue, not its protonation
                          state, so HIS-S / HIS-P would assert something it has
                          not measured. NA -- NEVER '' -- when the network put no
                          histidine on the binder; filter with .notna(), not
                          != ''. NO TOOL CONSUMES THIS CELL: hisbond's --centers
                          now takes 1-BASED positions as plain integers through
                          {expr} islands, and a packed string cannot ride through
                          one. Until hbdesigner publishes numeric position
                          columns, read the positions off here by hand and ADD 1.
    n_centers             how many such histidines. 0 IS A MEASUREMENT (this
                          network designed none) and pairs with centers NA. Both
                          are NA together only when fold_sequence is NA, i.e.
                          when there is no binder chain to index into at all.
    fold_seq_note        NA when the fold sequence was built; otherwise the
                          reason it was not. A row can be OK and still have no
                          fold sequence -- the commonest case is a network that
                          lands on two chains with neither omitted, so which one
                          is the binder is genuinely unknown. (A MONOMER run does
                          get a sequence: its one network chain is trivially the
                          only candidate.) This column says which case it was.

    path                  the designed structure. The GRAFTED file when chains were
                          grafted back, otherwise the raw rank PDB. Poly-glycine
                          everywhere except the network -- by design, not corruption.
                          DO NOT FOLD THIS; fold `fold_sequence`.
    rank                  1 = best by HBDesigner's own ordering (fewest buried
                          unsats, then highest saturation, then lowest HB_Score_full)
    hb_score_full         ref2015 energy of the network, per network residue,
                          relative to the same backbone as poly-glycine. MORE
                          NEGATIVE IS BETTER.
    hb_score_hb           the same difference under an hbond-only score function
                          (fa_rep 0.55 + hbond_sc + hbond_bb_sc). More negative is
                          better; this is the one that is about hydrogen bonds.
    avg_burial            mean sidechain-neighbour count over the network residues.
                          Higher = more buried. A buried network is the point.
    saturation            fraction of the network's polar groups that are satisfied
                          (range 0-2). Higher is better; --min-sat filtered on it.
    buried_heavy_unsats   buried unsatisfied heavy atoms left in the network
                          (--max-buns, default 0 = none tolerated). Lower is better.
    buried_unsat_hpol     buried unsatisfied polar hydrogens (--max-buphs, default
                          5). Lower is better.
    network               HBDesigner's own description of the network, as
                          ``<chain><resnum><aa>`` joined by ':', e.g. 'A12S:A16T'.
                          PDB NUMBERING. Includes --anchor-res residues.
    network_seq           the network's amino acids in fixed_positions order ('STN')
    n_network_res         how many residues the network has. Usually --n-res, but
                          anchors and symmetrization can change it -- read it, do
                          not assume it.
    network_chains        chains the network touches, comma-joined, in PDB order.
                          This is what --chains-to-design must list for
                          fixed_positions to line up.
    fixed_positions       the network in ProteinMPNN's --fixed-positions syntax
                          ('12,16/5': ',' between positions, '/' between chains, one
                          group per chain in network_chains order). 1-BASED WITHIN
                          EACH CHAIN -- mapped through the output PDB's residue
                          order, never copied from the resnums in `network`.
    fix1 .. fixK          the same positions one per column, so a per-row hand-off
                          can be written as --fixed-positions '{hbdesigner_fix1},
                          {hbdesigner_fix2}/' (prosapia's {expr} islands resolve
                          integers only, so a whole string column cannot be used).
    resnum_offset         per network chain, ``mpnn position - PDB resnum`` when it
                          is constant over that chain, else '<chain>:var'. All
                          zeroes means the PDB was numbered 1..L per chain and
                          `network`'s resnums ARE the mpnn positions.
    grafted             'yes' if the omitted chains were restored from the input
                          PDB (upstream graft_seq.py), 'no' otherwise.
    graft_identity        fraction of the grafted chains' residues that now match
                          the input. 1.0 = fully restored; NA when nothing was
                          grafted. Below 1.0 means the graft did not do its job --
                          treat the structure as wrong, not merely imperfect.
    n_res_total           residues in the output structure (checked equal to the
                          input's; a mismatch is an error status, not a column).

Two things to read carefully:

  * A ``<parent>_hb0`` row is a FAILURE RECORD, not a design: HBDesigner can find
    no network and still exit 0, writing neither a PDB nor a CSV. That row carries
    an ``error:`` status, an empty path and NA metrics, so the parent's failure is
    visible in the table instead of looking like a design that was never attempted.
    Never read it as a design; always gate on ``hbdesigner_status == 'OK'``.
  * The fan-out is NOT ``--top-k`` per parent. Upstream keeps
    ``min(top_k, n_surviving)``, so a parent may contribute 5 rows, 1 row, or only
    the ``_hb0`` failure record.

A design whose TSV is absent altogether (task still queued, killed, or never
submitted) is SKIPPED with a message -- no evidence is not the same as failure.

Usage:
    sapia collect hbdesigner outputs/<run> --table <the table the run reserved>
"""

import csv
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

FLOAT_COLUMNS = [
    "hb_score_full",
    "hb_score_hb",
    "avg_burial",
    "saturation",
    "graft_identity",
]
INT_COLUMNS = [
    "rank",
    "buried_heavy_unsats",
    "buried_unsat_hpol",
    "n_network_res",
    "resnum_shift_max",
    "n_res_total",
    "guide_atoms_removed",
    "n_fold_changes",
    "n_centers",
]
STR_COLUMNS = [
    "network",
    "network_seq",
    "network_chains",
    "fixed_positions",
    "resnum_offset",
    "grafted",
    "pdb_parser",
    "fold_sequence",
    "fold_chain",
    "fold_seq_note",
    "centers",
]
# `fixed_idx` is exploded into fix1..fixK instead of being collected as a column.
DATA_COLUMNS = FLOAT_COLUMNS + INT_COLUMNS + STR_COLUMNS


def _value(row: dict[str, str], key: str) -> str | None:
    """The cell, or None when absent or blank (NA, never a silent zero)."""
    value = row.get(key)
    if value is None or str(value).strip() == "":
        return None
    return str(value).strip()


def _fields(row: dict[str, str]) -> dict[str, Any]:
    """One TSV row -> collected columns. A blank stays NA rather than becoming 0."""
    data: dict[str, Any] = {c: None for c in DATA_COLUMNS}
    for col in FLOAT_COLUMNS:
        if (value := _value(row, col)) is not None:
            data[col] = float(value)
    for col in INT_COLUMNS:
        if (value := _value(row, col)) is not None:
            # Written through pandas upstream, so an int can arrive as '0.0'.
            data[col] = int(float(value))
    for col in STR_COLUMNS:
        if (value := _value(row, col)) is not None:
            data[col] = value
    if (idx := _value(row, "fixed_idx")) is not None:
        for i, position in enumerate(idx.split(","), start=1):
            data[f"fix{i}"] = int(position)
    return data


def _print_trust_counts(out_dir: Path) -> None:
    """One scan of every per-design TSV, two count lines about the trust columns.

    The counts exist so an NA is never read as a zero: ``guide_atoms_removed`` is
    NA only when the guide-atom strip pass did not run or failed, and a row read
    by Biopython was produced by a different parser from a row read by gemmi
    (Modal has gemmi, the vib conda env does not).
    """
    parsers: dict[str, int] = {}
    changes: dict[int, int] = {}
    centres: dict[int, int] = {}
    notes: dict[str, int] = {}
    ok_rows = no_strip = with_seq = with_centres = 0
    for tsv in sorted(out_dir.glob("*.tsv")):
        try:
            with open(tsv, newline="") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
        except OSError:
            continue
        for row in rows:
            if (_value(row, "status") or "") != "OK":
                continue
            ok_rows += 1
            parsers[_value(row, "pdb_parser") or "unrecorded"] = (
                parsers.get(_value(row, "pdb_parser") or "unrecorded", 0) + 1
            )
            if _value(row, "guide_atoms_removed") is None:
                no_strip += 1
            if _value(row, "fold_sequence") is not None:
                with_seq += 1
                n = _value(row, "n_fold_changes")
                if n is not None:
                    key = int(float(n))
                    changes[key] = changes.get(key, 0) + 1
                nc = _value(row, "n_centers")
                if nc is not None:
                    ckey = int(float(nc))
                    centres[ckey] = centres.get(ckey, 0) + 1
                if _value(row, "centers") is not None:
                    with_centres += 1
            else:
                note = _value(row, "fold_seq_note") or "no reason recorded"
                notes[note] = notes.get(note, 0) + 1
    summary = ", ".join(f"{n} {name}" for name, n in sorted(parsers.items())) or "none"
    print(f"hbdesigner: {ok_rows} OK row(s) on disk; PDB reader: {summary}")
    print(
        f"hbdesigner: {no_strip} of {ok_rows} OK row(s) carry NO guide-atom strip "
        f"record (guide_atoms_removed is NA, which is not 0)"
    )
    dist = ", ".join(f"{k}:{v}" for k, v in sorted(changes.items())) or "none"
    print(
        f"hbdesigner: {with_seq} of {ok_rows} OK row(s) have a fold_sequence "
        f"(the foldable one); n_fold_changes distribution: {dist}"
    )
    cdist = ", ".join(f"{k}:{v}" for k, v in sorted(centres.items())) or "none"
    print(
        f"hbdesigner: {with_centres} of {with_seq} row(s) with a fold_sequence "
        f"carry designed histidine centres (hisbond-ready); n_centers "
        f"distribution: {cdist}"
    )
    for note, n in sorted(notes.items(), key=lambda kv: -kv[1]):
        print(f"hbdesigner:   {n} OK row(s) have NO fold_sequence -- {note}")


def collect_hbdesigner(ctx: CollectCtx) -> CollectEach:
    """Per-parent hbdesigner collector. hbdesigner is a create tool: the framework
    iterates the ready parents and this mints one child row per kept network
    (``<parent>_hb<rank>``), carrying ``parent`` for lineage. The framework stamps
    status/path/parent_name from each Collected. Rows are rebuilt from disk on every
    collect, so re-collecting is idempotent."""

    _print_trust_counts(ctx.out_dir)

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"
        if not tsv_path.is_file():
            # No result file at all: the task did not run or did not finish. That is
            # not a design failure, so nothing is minted -- re-collect after the
            # task completes.
            print(f"{d.name}: no result TSV at {tsv_path.name} (skipping)")
            return

        with open(tsv_path, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))

        for row in rows:
            rank = int(float(_value(row, "rank") or 0))
            status = _value(row, "status") or "error: no status recorded"
            path = _value(row, "path") or ""
            yield Collected(
                # rank 0 never collides with a real design: upstream ranks from 1.
                name=f"{d.name}_hb{rank}",
                parent=d.name,
                path=path,
                status=status,
                # Whatever the worker managed to record is kept even on an error
                # row (blanks collect as NA); the status is the gate, not the data.
                data=_fields(row),
            )

    return one
