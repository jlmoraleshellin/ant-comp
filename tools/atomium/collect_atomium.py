"""Collect AtomiUM FASTA outputs into a table.

Scans the output directory for the batched subgroup dirs, parses the FASTAs and
writes one child row per designed sequence.

Lineage is derived from the directory structure:

    <output_dir>/grp_<digest>/seqs/<fasta_stem>.fa

  * grp_<digest> = a batched subgroup run by atomium.sh (several designs sharing
    params, run in one atomium.py call).
  * fasta_stem   = the staged input filename AtomiUM processed. The submitter
    symlinks each input as ``<design_name>.pdb``, so the stem IS the immediate
    parent-table row -- stamped as ``parent_name``.

WHAT THE FASTA CONTAINS -- MEASURED, 2026-10-14
-----------------------------------------------
AtomiUM's FASTA is NOT ProteinMPNN's, so this does not share proteinmpnn's parser
(verified against ``atomium.py``, which writes the file directly)::

    >_init_sequence, chains=D
    MKT...                                     <- the native input chain(s)
    >gen_sample_0, temperature=0.1, seq_rec=0.43
    MRE...

A 2-design probe on real output (run_dir
``20261009_130629_5M2J_vhh_redmotif_lowtemp_scaleup``, child table
``table1_atomium_probe``, ``--chains-to-design D`` with 17 pinned positions)
established, from the raw files:

  * the generated sequence is the **FULL designed chain, pinned residues
    included** -- chain-D length 144 and generated length 144 on one design, 120
    and 120 on the other; the 17 pinned characters read ``TSNYWYENTNLIKRSPF`` in
    the native chain and in all four generated sequences alike;
  * ``_init_sequence`` is **byte-identical to the native designed chain**, not to
    the native minus the pinned positions;
  * no ``/`` or ``:`` appeared in any entry (one designed chain).

So ``atomium_sequence`` IS the designed chain and nothing is reinserted here.
An earlier version of the skill claimed the opposite (fixed residues omitted, a
``sequence_full`` column carrying the whole chain); that claim was wrong and the
column never existed. **Do not add reinsertion** -- it would double-insert.

Like ProteinMPNN's, the FASTA covers the DESIGNED chains only; a chain that is
not designed never appears in it. **Above ONE designed chain the ``/`` join and the
order of chains within it are ASSUMED, not measured** -- the probe had a single
designed chain, and the convention is read off ProteinMPNN's printing code. So every
check here that could depend on chain order is written as a MULTISET comparison, and
``init_order_matches`` records the element-wise answer without ever failing a row.
The first real multi-chain run settles the question from the table.

TRUST COLUMNS -- WHY THEY EXIST
-------------------------------
The only reason the question above needed a GPU probe is that no column could
answer it. Every row now carries the evidence, taken from the task's own
intermediate files in the same ``grp_`` dir:

  * ``parsed_pdbs.jsonl``   -- one JSON line per design; ``seq_chain_<X>`` is that
    chain's sequence as ``parse_multiple_chains.py`` read it (confirmed
    byte-identical to the native chain D in the probe);
  * ``assigned_pdbs.jsonl`` -- one JSON object, ``{design: [[designed],[fixed]]}``;
    absent when no ``--chains-to-design`` was given (then every chain is designed);
  * ``fixed_pdbs.jsonl``    -- one JSON object, ``{design: {chain: [1-based pos]}}``;
    absent when nothing was pinned.

``init_matches_parsed`` and ``len_matches_parsed`` check the FRAME (is the row's
sequence the chain we think it is, at the length we think it has);
``fixed_held`` checks that every pinned position really carries the input
residue. ``input_pdb`` records WHICH structure the sequence was designed on, read
from the symlink the submitter staged -- note that the trust columns are computed
against the same submit's ``parsed_pdbs.jsonl``, so a FASTA left behind by an earlier
failed submit would be self-consistent and read as fully trusted; ``input_pdb`` does
not close that race, it makes it visible. A failure of any of them makes the row an **error row** rather than a
quietly wrong sequence: its status starts with ``error:`` and the trust column
reads ``0``, which therefore never appears on a row that is fine.

NO SCORE COLUMN
---------------
There is no ``score`` or ``global_score`` -- AtomiUM reports only sequence
recovery against the input, so a design cannot be ranked by model likelihood the
way an mpnn ``score`` column allows. Two consequences worth knowing:

  * ``seq_rec`` is similarity to the INPUT sequence, not a quality score.
  * ``sample`` is NOT unique per row. AtomiUM numbers samples
    ``n % num_seq_per_target`` and walks temperatures with
    ``n // num_seq_per_target``, so with several ``--sampling-temp`` values the
    same ``sample`` id recurs once per temperature. The row key therefore uses the
    running entry index (``<parent>_a<i>``), and (``sample``, ``temperature``)
    together identify a draw.

Each row carries only ``parent_name``; ancestor values are resolved on demand by
walking the lineage, so nothing is propagated here.

Usage:
    sapia collect atomium outputs/<run> --table <the table the run reserved>
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import pandas as pd

from prosapia.core import (
    Collected,
    CollectCtx,
    CollectEach,
    DesignCtx,
)

# The header AtomiUM gives the native sequence it echoes first in every FASTA.
# The REAL header carries trailing fields -- measured: ">_init_sequence, chains=D"
# -- so this is matched against the header's FIRST comma-separated field, never
# against the whole line. An exact-equality test on the whole header silently
# never fired, and the echo was collected as a design (`<parent>_a0`) at status
# OK carrying the native input sequence.
INIT_HEADER = "_init_sequence"

# Trailing ``key=value`` fields on a generated header.
HEADER_FIELDS: Dict[str, type] = {
    "temperature": float,
    "seq_rec": float,
}

SEQ_CHAIN_PREFIX = "seq_chain_"


def is_init_header(header: str) -> bool:
    """True for AtomiUM's echo of the native input sequence.

    Matched on the FIRST comma-separated field so any number of trailing
    ``key=value`` fields (``chains=D`` today, more tomorrow) is tolerated, while a
    generated header -- whose first field is always ``gen_sample_<k>`` -- can never
    match. Matching by header rather than by position is deliberate: a future extra
    preamble entry must not shift which row is dropped.
    """
    return header.split(",", 1)[0].strip() == INIT_HEADER


def parse_atomium_header(header: str) -> Dict[str, Any]:
    """Parse ``gen_sample_0, temperature=0.1, seq_rec=0.43`` into columns.

    The sample id is the ``gen_sample_<k>`` prefix rather than a ``key=value`` pair.
    Unparseable fields become NA rather than raising: a partially written header
    should still collect its sequence.
    """
    out: Dict[str, Any] = {"sample": pd.NA}
    out.update({k: pd.NA for k in HEADER_FIELDS})
    for part in header.split(","):
        part = part.strip()
        if part.startswith("gen_sample_"):
            try:
                out["sample"] = int(part[len("gen_sample_") :])
            except ValueError:
                pass
            continue
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        key, value = key.strip(), value.strip()
        if key in HEADER_FIELDS:
            try:
                out[key] = HEADER_FIELDS[key](value)
            except ValueError:
                out[key] = pd.NA
    return out


def parse_fasta(fasta_path: Path) -> List[Tuple[str, str]]:
    entries: List[Tuple[str, str]] = []
    header = ""
    seq_lines: List[str] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    entries.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header:
        entries.append((header, "".join(seq_lines)))
    return entries


# ---------------------------------------------------------------------------
# The task's own intermediate files, read back as evidence
# ---------------------------------------------------------------------------


@dataclass
class GroupFiles:
    """The three helper jsonls of one ``grp_`` dir, read once and cached."""

    parsed: Dict[str, dict] = field(default_factory=dict)
    assigned: Dict[str, list] = field(default_factory=dict)
    fixed: Dict[str, dict] = field(default_factory=dict)
    # Whether the optional files EXIST, which is not the same as whether they hold
    # this design. Absent = the flag was never passed (nothing pinned, every chain
    # designed); present but missing this design's key = an unexplained gap, which
    # is an error rather than a confident zero.
    has_assigned: bool = False
    has_fixed: bool = False
    error: str | None = None


def _read_json_lines(path: Path) -> List[Any]:
    objs: List[Any] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                objs.append(json.loads(line))
    return objs


def load_group(grp_dir: Path) -> GroupFiles:
    """Read ``parsed_``/``assigned_``/``fixed_pdbs.jsonl`` from one subgroup dir.

    ``parsed`` is required -- without it nothing about a row can be checked, so its
    absence or corruption is an error on every row of the group rather than a
    silently unverified sequence. The other two are optional by construction:
    ``assigned_pdbs.jsonl`` exists only with ``--chains-to-design``, and
    ``fixed_pdbs.jsonl`` only with ``--fixed-positions``.
    """
    group = GroupFiles()
    parsed_path = grp_dir / "parsed_pdbs.jsonl"
    if not parsed_path.is_file():
        group.error = f"no parsed_pdbs.jsonl in {grp_dir.name}"
        return group
    try:
        for entry in _read_json_lines(parsed_path):
            group.parsed[str(entry["name"])] = entry
    except (OSError, ValueError, KeyError) as exc:
        group.error = f"unreadable {grp_dir.name}/parsed_pdbs.jsonl ({exc})"
        return group

    for attr, fname in (("assigned", "assigned_pdbs.jsonl"), ("fixed", "fixed_pdbs.jsonl")):
        path = grp_dir / fname
        if not path.is_file():
            continue
        try:
            merged: Dict[str, Any] = {}
            for obj in _read_json_lines(path):
                merged.update(obj)
            setattr(group, attr, merged)
            setattr(group, f"has_{attr}", True)
        except (OSError, ValueError) as exc:
            group.error = f"unreadable {grp_dir.name}/{fname} ({exc})"
            return group
    return group


# Columns that describe the design's input, not one draw. NA until proven.
#
# NA, never "". An empty string does not survive the TSV round-trip -- pandas'
# read_csv returns NaN for it -- so a reader's `== ""` selects nothing. For the same
# family of reasons ``fixed_residues`` is written CHAIN-PREFIXED (``D:NA``): a
# two-position pin of Asn+Ala spells the bare string "NA", which is in read_csv's
# default na_values and would come back from the table as a null.
_NA_FACTS: Dict[str, Any] = {
    "input_pdb": pd.NA,
    "parsed_len": pd.NA,
    "n_chains_designed": pd.NA,
    "n_fixed": pd.NA,
    "init_matches_parsed": pd.NA,
    "init_order_matches": pd.NA,
    "fixed_positions": pd.NA,
    "fixed_residues": pd.NA,
}
_NA_DRAW: Dict[str, Any] = {
    "len_matches_parsed": pd.NA,
    "fixed_held": pd.NA,
}


def staged_input(grp_dir: Path, name: str) -> Any:
    """The structure this design was actually staged on, from the submitter's symlink.

    No column recorded this before, so a table could not answer "which backbone is
    this sequence for?" without re-deriving it from lineage. NA when the link is
    gone -- its absence is not made an error, because it is an audit column and the
    FASTA next to it is still what AtomiUM wrote.
    """
    link = grp_dir / "inputs" / f"{name}.pdb"
    try:
        if link.is_symlink():
            return str(link.readlink())
        if link.exists():
            return str(link)
    except OSError:
        pass
    return pd.NA


@dataclass
class Facts:
    """What the task's own files say about one design's input, plus a per-draw check."""

    columns: Dict[str, Any]
    error: str | None = None
    parsed_seqs: List[str] = field(default_factory=list)
    designed: List[str] = field(default_factory=list)
    fixed_pos: List[int] = field(default_factory=list)

    def check_draw(self, sequence: str) -> Tuple[Dict[str, Any], str | None]:
        """Per-draw trust columns: does this sequence fit the parsed input chain(s)?

        Length is compared as a MULTISET of per-chain lengths, so it does not depend
        on the order AtomiUM writes a multi-chain FASTA in (ProteinMPNN sorts the
        designed chains by their mask index, which is not necessarily the order
        ``assigned_pdbs.jsonl`` lists them in). ``fixed_held`` is only computed for a
        SINGLE designed chain: with several, the position-to-chain mapping of the
        ``/``-joined string is not established, and a guess there is exactly the
        class of error this column exists to catch.
        """
        cols = dict(_NA_DRAW)
        if self.error:
            return cols, None
        parts = sequence.split("/")
        want = sorted(len(s) for s in self.parsed_seqs)
        got = sorted(len(p) for p in parts)
        cols["len_matches_parsed"] = int(got == want)
        if got != want:
            return cols, (
                f"error: generated sequence has chain lengths {got} but the parsed "
                f"input chain(s) {','.join(self.designed)} have {want}"
            )
        if len(self.designed) != 1 or not self.fixed_pos:
            return cols, None
        seq, native = parts[0], self.parsed_seqs[0]
        wrong = [p for p in self.fixed_pos if seq[p - 1] != native[p - 1]]
        cols["fixed_held"] = int(not wrong)
        if wrong:
            shown = ", ".join(
                f"{p}:{native[p - 1]}->{seq[p - 1]}" for p in wrong[:6]
            )
            return cols, (
                f"error: {len(wrong)} of {len(self.fixed_pos)} fixed position(s) "
                f"were NOT held in chain {self.designed[0]}: {shown}"
            )
        return cols, None


def design_facts(
    group: GroupFiles, name: str, entries: List[Tuple[str, str]], grp_dir: Path
) -> Facts:
    """Everything checkable about one design's input, read from the group's jsonls."""
    cols: Dict[str, Any] = dict(_NA_FACTS)
    cols["input_pdb"] = staged_input(grp_dir, name)

    if group.error:
        return Facts(cols, f"error: {group.error}")

    entry = group.parsed.get(name)
    if entry is None:
        return Facts(
            cols,
            f"error: no entry for {name!r} in parsed_pdbs.jsonl "
            f"(the FASTA and the parsed inputs disagree about this design)",
        )

    chains_all = [k[len(SEQ_CHAIN_PREFIX) :] for k in entry if k.startswith(SEQ_CHAIN_PREFIX)]
    # No assigned_pdbs.jsonl at all means no --chains-to-design, i.e. every chain
    # designed. The FILE existing while this design is absent from it is a different
    # thing entirely, and is not survivable as a guess.
    assigned = group.assigned.get(name)
    if group.has_assigned and assigned is None:
        return Facts(
            cols,
            f"error: assigned_pdbs.jsonl exists but has no entry for {name!r}, so "
            f"which chains were designed is unknown (it is not safe to assume all)",
        )
    designed = [str(c) for c in assigned[0]] if assigned else list(chains_all)

    absent = [c for c in designed if f"{SEQ_CHAIN_PREFIX}{c}" not in entry]
    if absent:
        return Facts(
            cols,
            f"error: designed chain(s) {','.join(absent)} are not in the parsed "
            f"input (chains present: {','.join(chains_all) or '<none>'}) -- the "
            f"--chains-to-design letter does not exist in this structure",
        )

    parsed_seqs = [str(entry[f"{SEQ_CHAIN_PREFIX}{c}"]) for c in designed]
    cols["parsed_len"] = sum(len(s) for s in parsed_seqs)
    cols["n_chains_designed"] = len(designed)

    # A present fixed_pdbs.jsonl means positions WERE requested. A gap in it is
    # therefore an error, not a confident n_fixed of 0: zero is a measurement, and a
    # motif-graft row claiming nothing was pinned is exactly the silent wrong answer
    # these columns exist to prevent. (make_fixed_positions_dict.py writes every
    # chain of every parsed design, [] included.)
    if group.has_fixed and name not in group.fixed:
        return Facts(
            cols,
            f"error: fixed_pdbs.jsonl exists but has no entry for {name!r}, so the "
            f"pinned positions are unknown (n_fixed is NA, not 0)",
        )
    fixed_map = group.fixed.get(name, {})
    if group.has_fixed:
        unlisted = [c for c in designed if c not in fixed_map]
        if unlisted:
            return Facts(
                cols,
                f"error: fixed_pdbs.jsonl has no key for designed chain(s) "
                f"{','.join(unlisted)} of {name!r} (keys: "
                f"{','.join(sorted(map(str, fixed_map))) or '<none>'}) -- the pinned "
                f"positions of that chain are unknown",
            )
    fixed_by_chain = {c: [int(p) for p in fixed_map.get(c, [])] for c in designed}
    cols["n_fixed"] = sum(len(v) for v in fixed_by_chain.values())
    audit = ";".join(
        f"{c}:{','.join(str(p) for p in fixed_by_chain[c])}"
        for c in designed
        if fixed_by_chain[c]
    )
    cols["fixed_positions"] = audit or pd.NA

    # A pinned position outside the chain is an error, never a trimmed list: it
    # means the spec and the structure disagree about the numbering.
    for c, positions in fixed_by_chain.items():
        length = len(parsed_seqs[designed.index(c)])
        out_of_range = [p for p in positions if p < 1 or p > length]
        if out_of_range:
            return Facts(
                cols,
                f"error: fixed position(s) {out_of_range[:6]} are outside chain "
                f"{c} (1..{length})",
            )

    init = [seq for header, seq in entries if is_init_header(header)]
    if not init:
        return Facts(
            cols,
            "error: no _init_sequence entry in the FASTA -- the frame of the "
            "generated sequences cannot be checked against the parsed input",
        )
    init_parts = init[0].split("/")
    # Multiset: the pass/fail column must not depend on an unverified chain order.
    init_ok = sorted(init_parts) == sorted(parsed_seqs)
    cols["init_matches_parsed"] = int(init_ok)
    # Element-wise, in `designed` order. PURELY INFORMATIONAL -- it never changes a
    # status. 0 here with init_matches_parsed 1 means the FASTA's chain order is not
    # the order assigned_pdbs.jsonl lists, which is the open multi-chain question.
    cols["init_order_matches"] = int(init_parts == parsed_seqs)

    single = designed[0] if len(designed) == 1 else None
    if single is not None and fixed_by_chain[single]:
        cols["fixed_residues"] = f"{single}:" + "".join(
            parsed_seqs[0][p - 1] for p in fixed_by_chain[single]
        )

    if not init_ok:
        return Facts(
            cols,
            f"error: _init_sequence does not equal the parsed designed chain(s) "
            f"{','.join(designed)} (echo length {len(init[0])}, parsed "
            f"{cols['parsed_len']}) -- the FASTA is not the chain this row claims",
        )

    return Facts(
        cols,
        None,
        parsed_seqs=parsed_seqs,
        designed=designed,
        fixed_pos=fixed_by_chain[single] if single is not None else [],
    )


def collect_atomium(ctx: CollectCtx) -> CollectEach:
    """Per-parent AtomiUM collector. atomium is a create tool: the framework iterates
    the ready parents and this mints one child row (``<parent>_a<i>``) per sampled
    sequence, carrying ``parent`` for lineage. The framework stamps status/path/
    parent_name from each Collected; child rows are discovered on disk, so re-running
    rebuilds them (idempotent).

    ``i`` is the running index over ALL FASTA entries, the skipped ``_init_sequence``
    echo included -- so the first draw is ``_a1``, not ``_a0``. Row names are lineage
    keys and the documented naming of every existing atomium table, so they are left
    exactly as they were when the echo was (wrongly) being collected as ``_a0``.
    """
    # Each grp_<digest>/ output dir holds seqs/<design>.fa, where the staged input
    # was symlinked as <design>.pdb -- so the FASTA stem IS the parent-table row.
    fasta_by_parent: Dict[str, Path] = {}
    for subdir in sorted(p for p in ctx.out_dir.iterdir() if p.is_dir()):
        if subdir.name in ("atomium_logs", "atomium_tasks"):
            continue
        seqs_dir = subdir / "seqs"
        if not seqs_dir.is_dir():
            continue
        for fasta_path in seqs_dir.glob("*.fa"):
            fasta_by_parent[fasta_path.stem] = fasta_path

    group_cache: Dict[Path, GroupFiles] = {}

    def one(d: DesignCtx) -> Iterable[Collected]:
        fasta_path = fasta_by_parent.get(d.name)
        if fasta_path is None:
            print(f"{d.name}: no FASTA found (skipping)")
            return

        grp_dir = fasta_path.parent.parent
        if grp_dir not in group_cache:
            group_cache[grp_dir] = load_group(grp_dir)
        entries = parse_fasta(fasta_path)
        facts = design_facts(group_cache[grp_dir], d.name, entries, grp_dir)

        for i, (header, sequence) in enumerate(entries):
            # The first entry is AtomiUM's echo of the native input sequence.
            if is_init_header(header):
                continue

            # Bare column names; the driver leaf-prefixes them.
            data: Dict[str, Any] = {
                "iteration": i,
                "sequence": sequence,
                "seq_len": len(sequence.replace("/", "")),
            }
            data.update(parse_atomium_header(header))
            data.update(facts.columns)
            draw_cols, draw_error = facts.check_draw(sequence)
            data.update(draw_cols)
            yield Collected(
                name=f"{d.name}_a{i}",
                parent=d.name,
                path=fasta_path,
                status=facts.error or draw_error or "OK",
                data=data,
            )

    return one
