#!/usr/bin/env python3
"""
Build a forced-template mmCIF for Boltz -- and the ``templates:`` YAML that goes
with it -- from a subset of chains of a reference structure.

NOT a sapia tool (no spec.py, so `discover()` never picks it up). It is a
workstation-side utility, run once per template while staging inputs. In the
workstation it lives at ``/runs/utils/boltz_template.py`` (NOT under
``$PROSAPIA_TOOLS_DIR/_utils/``, which does not exist and exits 2). ``/runs`` is
the runs VOLUME, so the workstation runs whatever copy was last staged there:
an edit in the repo does not reach it until it is staged again (a campaign may
also stage its own copy inside its run_dir -- check the md5 of the copy you
call). So:

    sapia modal-shell --cmd 'python /runs/utils/boltz_template.py \\
        --fetch 7OJG --out /runs/inputs/7ojg_tmpl_kabc.cif \\
        --chains K,A,B,C --rename-to A,B,C,D --residues K:19-59+107-155 \\
        --yaml /runs/inputs/7ojg_tmpl_kabc.yaml'

Write outputs under ``/runs/`` -- every ``modal-shell`` call is a FRESH
container, so a file left in /tmp does not exist for the next one.

It exists because the template CIF is the one input in this workbench whose
failure modes are SILENT and EXPENSIVE, and because writing it ad hoc gets the
same four things wrong every time (see `.claude/skills/boltz/SKILL.md`):

1. **Boltz names template chains by ``label_asym_id``, not the auth chain ID.**
   A CIF written by gemmi's own ``setup_entities()`` labels subchains ``Axp``,
   ``Bxp``, ... and ``template_id: [A, ...]`` then dies with *"Template chain A
   is not one of the protein chains"*. Here every chain is given
   ``label_asym_id == auth_asym_id`` by construction.

2. **``_entity_poly_seq`` must equal the MODELLED residue count of every
   entity.** Boltz's ``parse_polymer`` uses it as the token list, and a declared
   residue absent from the model still consumes a token index -- so a
   too-long declaration shifts every residue index downstream and nothing
   errors. The classic way in: residues cloned from a source file keep their
   source ``subchain`` label, ``setup_entities()`` merges two chains into one
   entity, and ``full_sequence`` is taken from whichever came first (a real
   case: ``_entity_poly_seq`` = 49 for 41-residue chains). Here
   ``setup_entities()`` is never called -- one entity is built per output chain
   with ``full_sequence`` taken from the residues actually written, so the two
   cannot disagree.

3. **``_atom_site.label_seq_id`` must actually be written, 1..N per chain.**
   It is the index Boltz uses; a ``.`` leaves the residue with no token. gemmi's
   ``assign_label_seq_id()`` only numbers residues marked
   ``entity_type == Polymer``, and that flag is inherited from the SOURCE file --
   a plain PDB or a bare ``_atom_site`` mmCIF carries no entity annotation, so
   the clones arrive ``Unknown`` and every residue is skipped in silence.
   ``set_entities()`` therefore sets it, and raises if any residue still has no
   ``label_seq``.

4. **Boltz falls back to a sequence SEARCH unless ``chain_id`` and
   ``template_id`` are both given and equal length**, which with sequence-
   identical chains (any homo-oligomer) can assign the wrong template chain and
   destroy the geometry without a word. The emitted YAML always writes both.

5. Ligands, waters and a second model quietly become extra entities. Only the
   first model is used and only tabulated amino/nucleic acids are kept
   (``--keep-het`` to opt out).

A build writes ``<out>.tmp`` (and ``<yaml>.tmp``) first and renames them into
place only when ``verify()`` passes; on failure both temporaries are deleted.
So **an absent output file is the failure signal** -- which matters under
``sapia modal-shell --cmd``, which always exits 0.

Every write is followed by ``verify()``, which RE-READS the written file as raw
mmCIF categories -- not through gemmi's structure view, which would paper over
exactly the mismatch being looked for -- and counts ``_entity_poly_seq.mon_id``
rows per entity against the modelled residues of its subchains. A failed check
exits 1 and says which chain and which count. Run it on a CIF from anywhere with
``--check <file>``.

What it deliberately does NOT do: decide which chains to template. Do not
template the binder chain, and remember a forced template STEERS, it does not
constrain -- superpose the predicted target chains back onto this file
afterwards and drop the rows where the target did not land.

UNRESOLVED residues: ``--keep-gaps``
------------------------------------
By default the template declares exactly what is MODELLED, and the residues are
indexed 1..N in the order they appear. That is right when the source chain is
continuous, and WRONG when it is not: a crystal structure with an unbuilt loop
then has its declared sequence silently compacted, so position i of the template
is no longer position i of the query and every residue after the first gap is
off by the size of that gap. Boltz does not notice; the prediction just gets a
forced template that is subtly shifted.

``--keep-gaps`` (off by default, nothing changes without it) keeps the declared
sequence COMPLETE over the requested ``--residues <Src>:a-b`` range:

* ``_entity_poly_seq`` holds all ``b - a + 1`` positions in order, taking the
  residue name of an unresolved position from the SOURCE file's own
  ``_pdbx_poly_seq_scheme`` (which lists every declared residue, modelled or
  not). Positions are matched by AUTHOR number, never by position in the source
  entity -- the entity is usually longer than the range.
* ``_atom_site.label_seq_id`` of a modelled residue is its position in that
  range (``author number - a + 1``), so it has GAPS.
* an unresolved position gets NO atoms. Coordinates are never invented.

So the template is 1:1 with a query that carries the full range, which is the
only mapping Boltz can use. One range per chain, each chain still its own
entity. ``--check`` understands such a file (see ``--expect-declared``) and
prints how many positions are unresolved and where.

What Boltz actually does with it (read from jwohlwend/boltz
b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc), and why the checks are where they
are:

* ``parse_polymer`` aligns the MODELLED residue NAMES to the entity's
  ``full_sequence`` with ``gemmi.align_sequence_to_polymer`` -- it does NOT
  read ``label_seq_id``. A declared residue with no atoms simply becomes a
  token with ``is_present=False``, zero coordinates, masked off, and the force
  potential leaves it unrestrained. Gaps are therefore harmless; what is NOT
  harmless is a modelled residue whose name differs from the declared name at
  its place (the alignment treats it as absent without advancing the polymer
  index, desyncing everything after it) or a modelled residue not declared at
  all (an overrun, boltz issue 669). Both are refused at build time and by
  ``--check``.
* the query-to-template map, with ``chain_id`` and ``template_id`` both given,
  is a local BLOSUM62 alignment with gap penalty -1000 -- one gapless diagonal
  -- against the DECLARED sequence, then positional, with no identity
  threshold. A shifted or wrong-length declaration is used anyway, silently.
  ``--expect-sequence <query>`` compares the declared sequence with the query
  the prediction will carry and fails on a length difference or on more
  differing positions than ``--allow-mismatches`` (default 0) -- a shift by k
  differs almost everywhere, so it fails.
* ``verify()`` also runs the alignment Boltz runs
  (``gemmi.align_sequence_to_polymer`` of the declared sequence against the
  written polymer, default scoring) and fails unless every modelled residue
  lands on the declared position equal to its ``label_seq_id``.
* gaps are never filled. Invented coordinates would be restrained by the
  ``force``/``threshold`` potential like real ones.

Modes
-----
build   --out with --chains (and --in or --fetch)
check   --check <file.cif> [--expect-chains A,B,C] [--expect-declared N]
        -- verify only, writes nothing

Self-test (no network, no Modal): ``uv run python utils/test_boltz_template.py``.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import gemmi

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.cif.gz"


# --------------------------------------------------------------------------
# argument parsing
# --------------------------------------------------------------------------


def split_list(value: str) -> list[str]:
    """'A, B ,C' -> ['A', 'B', 'C']; '' -> []."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


def parse_residues(value: str) -> dict[str, list[tuple[int, int]]]:
    """'K:19-59+107-155,A:36-149' -> {'K': [(19,59), (107,155)], 'A': [(36,149)]}.

    Keyed by the SOURCE chain ID, applied before any rename. Bounds are
    inclusive and compared against the source residue numbers, so they read the
    same as the numbering in the reference file (and in the paper).
    """
    ranges: dict[str, list[tuple[int, int]]] = {}
    for spec in split_list(value):
        if spec.count(":") != 1:
            raise ValueError(
                f"--residues entry {spec!r} must be '<chain>:<first>-<last>[+<first>-<last>...]'"
            )
        chain, segments = (tok.strip() for tok in spec.split(":"))
        parsed: list[tuple[int, int]] = []
        for seg in segments.split("+"):
            first_str, _, last_str = seg.strip().partition("-")
            if not first_str or not last_str:
                raise ValueError(f"--residues segment {seg!r} must be '<first>-<last>'")
            first, last = int(first_str), int(last_str)
            if last < first:
                raise ValueError(f"--residues segment {seg!r} ends before it starts")
            parsed.append((first, last))
        ranges.setdefault(chain, []).extend(parsed)
    return ranges


def resolve_chains(chains: str, rename_to: str) -> list[tuple[str, str]]:
    """[(source_id, output_id), ...] in the order given; rename is positional."""
    kept = split_list(chains)
    renamed = split_list(rename_to)
    if not kept:
        raise ValueError("--chains is empty")
    if renamed and len(renamed) != len(kept):
        raise ValueError(
            f"--rename-to has {len(renamed)} IDs but --chains has {len(kept)}"
        )
    out_ids = renamed or kept
    if len(set(out_ids)) != len(out_ids):
        raise ValueError(f"output chain IDs are not unique: {','.join(out_ids)}")
    return list(zip(kept, out_ids))


# --------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------


def fetch_rcsb(pdb_id: str, dest: Path) -> Path:
    """Download an entry from RCSB to ``dest``, uncompressed.

    Done here rather than with ``modal volume put`` because this network
    intercepts TLS to Modal's blob storage and the upload fails; fetching from
    inside the workstation sidesteps it entirely.
    """
    url = RCSB_URL.format(pdb_id=pdb_id.upper())
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=120) as resp, gzip.GzipFile(
        fileobj=resp
    ) as gz, dest.open("wb") as out:
        shutil.copyfileobj(gz, out)
    print(f"fetched {url} -> {dest} ({dest.stat().st_size} bytes)")
    return dest


def load_structure(path: Path) -> gemmi.Structure:
    """Read a PDB or mmCIF file, single conformer, no hydrogens."""
    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    st = gemmi.read_structure(str(path))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    if len(st) == 0:
        raise ValueError(f"no model in {path}")
    return st


def is_polymer_residue(res: gemmi.Residue) -> bool:
    """True for a tabulated amino-acid or nucleic-acid residue (MSE included)."""
    info = gemmi.find_tabulated_residue(res.name)
    return bool(info and (info.is_amino_acid() or info.is_nucleic_acid()))


def one_letter(res: gemmi.Residue) -> str:
    info = gemmi.find_tabulated_residue(res.name)
    code = info.one_letter_code if info else ""
    return code.upper() if code else "X"


def one_letter_name(name: str) -> str:
    info = gemmi.find_tabulated_residue(name)
    code = info.one_letter_code if info else ""
    return code.upper() if code else "X"


def nonstandard(names: list[str]) -> list[tuple[int, str]]:
    """[(1-based position, name)] for every residue that is not one of the 20.

    gemmi's ``one_letter_code`` is LOWERCASE for a modified residue (MSE -> 'm')
    and this file upper-cases it, so an engineered residue silently reads as the
    parent amino acid in any sequence comparison. Boltz aligns the declared
    names, so it is not an error -- but it must be said out loud, never
    normalised in silence.
    """
    out: list[tuple[int, str]] = []
    for i, name in enumerate(names):
        info = gemmi.find_tabulated_residue(name)
        code = info.one_letter_code if info else ""
        if not code or not code.isupper():
            out.append((i + 1, name))
    return out


def condense(numbers: list[int]) -> str:
    """[1,2,3,7] -> '1-3,7'. For reporting unresolved stretches readably."""
    if not numbers:
        return "-"
    spans: list[tuple[int, int]] = []
    for n in sorted(numbers):
        if spans and n == spans[-1][1] + 1:
            spans[-1] = (spans[-1][0], n)
        else:
            spans.append((n, n))
    return ",".join(str(lo) if lo == hi else f"{lo}-{hi}" for lo, hi in spans)


def read_poly_seq_scheme(
    path: Path, chains: list[str] | None = None
) -> dict[str, list[tuple[int, str]]]:
    """``_pdbx_poly_seq_scheme`` as {auth chain: [(author number, mon_id), ...]}.

    This is the only place the DECLARED sequence of the source entity can be
    read: ``_atom_site`` has the modelled residues and nothing else, and
    gemmi's ``Structure`` view of a chain is the modelled residues too. The
    category lists every declared residue of every polymer chain, resolved or
    not, with its ``seq_id`` in the entity and its author number.

    Only the source ``chains`` asked for are read (all, if None), so an oddity
    in a chain nobody templates cannot stop the build.

    Mapping is by AUTHOR number (``pdb_seq_num``, falling back to
    ``auth_seq_num``) because ``--residues`` is in author numbering and the
    entity usually extends past the requested range. For an unresolved residue
    that has neither -- some depositions leave both ``?`` -- the number is
    derived by counting from the resolved neighbours in ``seq_id`` order and the
    derivation is printed. It is REFUSED when the neighbours on the two sides
    disagree about the offset (a numbering jump inside the unresolved stretch:
    which residue absorbs it cannot be known), and when the chain has no
    resolved neighbour at all.

    Two rows sharing one ``seq_id`` are MICROHETEROGENEITY (two residue types
    deposited at one position, ``hetero = y``). The position is kept with both
    names joined by ``/`` so that ``declared_range`` refuses it if, and only if,
    it falls inside a requested range: a template must declare ONE name there.
    """
    doc = gemmi.cif.read(str(path))
    block = doc.sole_block()
    table = block.find(
        "_pdbx_poly_seq_scheme.",
        [
            "pdb_strand_id", "seq_id", "mon_id", "pdb_seq_num", "auth_seq_num",
            "?hetero", "?pdb_ins_code",
        ],
    )
    if not len(table):
        raise ValueError(
            f"{path} has no _pdbx_poly_seq_scheme -- the declared sequence of the "
            "source chains is unknown, so --keep-gaps cannot name the unresolved "
            "residues. Use a source file from the PDB, or drop --keep-gaps"
        )
    wanted = set(chains) if chains is not None else None
    # chain -> seq_id -> [names], author number, hetero flags, insertion codes
    raw: dict[str, dict[int, dict]] = {}
    for row in table:
        chain = row[0]
        if wanted is not None and chain not in wanted:
            continue
        seq_id, mon_id = int(row[1]), row[2]
        num: int | None = None
        for candidate in (row[3], row[4]):
            text = candidate.strip()
            if text not in ("", ".", "?"):
                num = int(text)
                break
        entry = raw.setdefault(chain, {}).setdefault(
            seq_id, {"names": [], "num": None, "hetero": set(), "icode": set()}
        )
        entry["names"].append(mon_id)
        if entry["num"] is None:
            entry["num"] = num
        if row.has(5):
            entry["hetero"].add(row.str(5))
        if row.has(6) and row.str(6).strip() not in ("", ".", "?"):
            entry["icode"].add(row.str(6).strip())

    scheme: dict[str, list[tuple[int, str]]] = {}
    for chain, by_seq in raw.items():
        seqs = sorted(by_seq)
        anchors = [(seq, by_seq[seq]["num"]) for seq in seqs if by_seq[seq]["num"] is not None]
        if not anchors:
            raise ValueError(
                f"chain {chain}: _pdbx_poly_seq_scheme gives no author number for any "
                "residue, so the declared sequence cannot be placed in author "
                "numbering"
            )
        derived: list[int] = []
        filled: list[tuple[int, str]] = []
        for seq in seqs:
            entry = by_seq[seq]
            num = entry["num"]
            if num is None:
                before = [a for a in anchors if a[0] < seq]
                after = [a for a in anchors if a[0] > seq]
                offsets = set()
                if before:
                    offsets.add(before[-1][1] - before[-1][0])
                if after:
                    offsets.add(after[0][1] - after[0][0])
                if len(offsets) != 1:
                    raise ValueError(
                        f"chain {chain}: seq_id {seq} ({entry['names'][0]}) has no "
                        "author number, and its resolved neighbours disagree about "
                        f"the offset ({sorted(offsets)}) -- the numbering jumps "
                        "somewhere inside the unresolved stretch, so the author "
                        "number cannot be derived by counting"
                    )
                num = seq + offsets.pop()
                derived.append(num)
            names = sorted(set(entry["names"]))
            filled.append((num, "/".join(names)))
        numbers = [num for num, _ in filled]
        if len(set(numbers)) != len(numbers):
            dup = sorted({n for n in numbers if numbers.count(n) > 1})
            icodes = sorted({c for e in by_seq.values() for c in e["icode"]})
            raise ValueError(
                f"chain {chain}: two declared residues share author number(s) "
                f"{condense(dup)} in _pdbx_poly_seq_scheme"
                + (f" (insertion codes present: {','.join(icodes)})" if icodes else "")
                + ", so a residue range in author numbering is ambiguous"
            )
        if derived:
            print(
                f"note: chain {chain}: author number derived by counting from "
                f"resolved neighbours for {len(derived)} unresolved residue(s): "
                f"{condense(derived)}",
                file=sys.stderr,
            )
        scheme[chain] = filled
    return scheme


# --------------------------------------------------------------------------
# build
# --------------------------------------------------------------------------


@dataclass
class ChainReport:
    source: str
    out: str
    n_res: int
    first: str
    last: str
    sequence: str
    # --keep-gaps only; zero/empty in the default contiguous mode.
    n_declared: int = 0
    declared_sequence: str = ""
    unresolved: list[int] = field(default_factory=list)  # AUTHOR numbers
    offset: int = 0  # label_seq_id = author number - offset


@dataclass
class BuildResult:
    structure: gemmi.Structure
    chains: list[ChainReport] = field(default_factory=list)


def select(
    st: gemmi.Structure,
    groups: list[tuple[str, str]],
    ranges: dict[str, list[tuple[int, int]]],
    strip_het: bool,
    renumber: bool,
    model_index: int,
    keep_gaps: bool = False,
    scheme: dict[str, list[tuple[int, str]]] | None = None,
) -> BuildResult:
    """Copy the requested chains into a fresh structure, one entity per chain.

    A requested chain that is not in the file, or a residue range that selects
    nothing, is an ERROR -- never a shorter template. A silently partial
    template is well-formed, folds, and is wrong.

    With ``keep_gaps`` the entity declares the WHOLE ``--residues`` range of
    that chain (names for the unresolved positions come from ``scheme``, the
    source's ``_pdbx_poly_seq_scheme``) and each modelled residue keeps its
    position in that range as ``label_seq_id``. Nothing else changes: the same
    atoms are written, under the same chain IDs.
    """
    if model_index >= len(st):
        raise ValueError(f"--model {model_index} but the file has {len(st)} model(s)")
    model = st[model_index]
    present = {chain.name: chain for chain in model}
    unknown = [src for src, _ in groups if src not in present]
    if unknown:
        raise ValueError(
            f"chain(s) {','.join(unknown)} not in the structure "
            f"(have: {','.join(sorted(present))})"
        )

    out_st = gemmi.Structure()
    out_st.name = st.name or "template"
    out_st.spacegroup_hm = "P 1"
    out_st.cell = gemmi.UnitCell()
    out_model = gemmi.Model(1)
    result = BuildResult(structure=out_st)

    declared_per_chain: dict[str, list[str]] = {}

    for src, out_id in groups:
        segments = ranges.get(src)
        declared_names: list[str] = []
        first_pos = 0
        if keep_gaps:
            declared_names, first_pos = declared_range(src, segments, scheme)
            declared_per_chain[out_id] = declared_names
        out_chain = gemmi.Chain(out_id)
        for res in present[src]:
            if strip_het and not is_polymer_residue(res):
                continue
            num = res.seqid.num
            if num is None:
                raise ValueError(
                    f"chain {src}: residue {res.name} has no sequence number, so it "
                    "has no position in the declared range"
                )
            if segments and not any(lo <= num <= hi for lo, hi in segments):
                continue
            copy = res.clone()
            # The whole point: no source subchain label survives into the new
            # file, so nothing can merge two chains into one entity behind our
            # back. label_asym_id is set explicitly to the auth chain ID.
            copy.subchain = out_id
            copy.het_flag = "A"
            copy.label_seq = None
            if renumber:
                copy.seqid = gemmi.SeqId(len(out_chain) + 1, " ")
            if keep_gaps:
                copy.label_seq = position_in_range(
                    src, res, num, first_pos, declared_names
                )
            out_chain.add_residue(copy)
        if len(out_chain) == 0:
            raise ValueError(
                f"chain {src} selected no residues"
                + (f" in range {ranges_str(segments)}" if segments else "")
            )
        if keep_gaps:
            modelled_nums = {res.seqid.num for res in out_chain}
            unresolved = [
                first_pos + i
                for i in range(len(declared_names))
                if first_pos + i not in modelled_nums
            ]
            print(
                f"note: chain {src} -> {out_id}: declaring all {len(declared_names)} "
                f"residues of {first_pos}-{first_pos + len(declared_names) - 1}; "
                f"{len(unresolved)} unresolved and written with no atoms"
                + (f" (author numbers {condense(unresolved)})" if unresolved else ""),
                file=sys.stderr,
            )
            odd = nonstandard(declared_names)
            if odd:
                print(
                    f"note: chain {src} -> {out_id}: {len(odd)} declared residue(s) "
                    "are not one of the 20 standard amino acids and read as their "
                    "parent letter in any sequence comparison: "
                    + ", ".join(
                        f"position {p} (author {first_pos + p - 1}) {n}"
                        for p, n in odd[:8]
                    )
                    + (" ..." if len(odd) > 8 else ""),
                    file=sys.stderr,
                )
        if segments and not keep_gaps:
            want = sum(hi - lo + 1 for lo, hi in segments)
            if len(out_chain) != want:
                print(
                    f"note: chain {src} -> {out_id}: {len(out_chain)} residues for a "
                    f"range spanning {want} -- the source has gaps or unmodelled "
                    f"residues there. The template declares what is MODELLED, which "
                    f"is correct; check the gap is where you think it is.",
                    file=sys.stderr,
                )
        out_model.add_chain(out_chain)
        result.chains.append(
            ChainReport(
                source=src,
                out=out_id,
                n_res=len(out_chain),
                first=str(out_chain[0].seqid.num),
                last=str(out_chain[len(out_chain) - 1].seqid.num),
                sequence="".join(one_letter(r) for r in out_chain),
                n_declared=len(declared_names) if keep_gaps else 0,
                declared_sequence=(
                    "".join(one_letter_name(n) for n in declared_names)
                    if keep_gaps
                    else ""
                ),
                unresolved=unresolved if keep_gaps else [],
                offset=(first_pos - 1) if keep_gaps else 0,
            )
        )

    out_st.add_model(out_model)
    if keep_gaps:
        # label_seq_id was set per residue above and MUST NOT be reassigned:
        # assign_label_seq_id() would compact it back to 1..N and undo the gaps.
        set_entities(out_st, declared=declared_per_chain, assign=False)
    else:
        set_entities(out_st)
    return result


def declared_range(
    src: str,
    segments: list[tuple[int, int]] | None,
    scheme: dict[str, list[tuple[int, str]]] | None,
) -> tuple[list[str], int]:
    """The full declared residue list of ``src`` over its single --residues range.

    Returns (names in order, first author number). Raises rather than trimming:
    a range running past the declared sequence means the caller is asking for
    positions this structure cannot supply, and a shorter declaration would
    shift the query mapping -- the exact failure --keep-gaps exists to stop.
    """
    if scheme is None:
        raise ValueError("--keep-gaps needs the source's _pdbx_poly_seq_scheme")
    if not segments or len(segments) != 1:
        raise ValueError(
            f"--keep-gaps needs exactly one '{src}:<first>-<last>' range for chain "
            f"{src}; got {ranges_str(segments) if segments else 'none'}. The "
            "declared sequence must be one contiguous stretch of the query"
        )
    lo, hi = segments[0]
    chain_scheme = scheme.get(src)
    if chain_scheme is None:
        raise ValueError(
            f"chain {src} is not a polymer chain in _pdbx_poly_seq_scheme "
            f"(have: {','.join(sorted(scheme))})"
        )
    by_num = dict(chain_scheme)
    missing = [p for p in range(lo, hi + 1) if p not in by_num]
    if missing:
        declared_nums = [n for n, _ in chain_scheme]
        raise ValueError(
            f"chain {src}: --residues {lo}-{hi} is not inside the declared sequence "
            f"({min(declared_nums)}-{max(declared_nums)} in author numbering); "
            f"{len(missing)} position(s) are not declared at all: {condense(missing)}"
        )
    mixed = [(p, by_num[p]) for p in range(lo, hi + 1) if "/" in by_num[p]]
    if mixed:
        raise ValueError(
            f"chain {src}: microheterogeneity inside {lo}-{hi} -- "
            + ", ".join(f"author {p} is {names}" for p, names in mixed)
            + " (two _pdbx_poly_seq_scheme rows share one seq_id). A template "
            "must declare one residue per position; pick a source without it or "
            "a range that excludes it"
        )
    return [by_num[p] for p in range(lo, hi + 1)], lo


def position_in_range(
    src: str,
    res: gemmi.Residue,
    num: int,
    first_pos: int,
    declared_names: list[str],
) -> int:
    """label_seq_id of a modelled residue under --keep-gaps, fully checked.

    A mismatch here means the modelled residue is not the residue the declared
    sequence puts at that position -- i.e. the author numbering of the file and
    of its own scheme disagree. That is unrecoverable, so it raises.
    """
    if res.seqid.icode.strip():
        raise ValueError(
            f"chain {src}: residue {num}{res.seqid.icode}{res.name} has an insertion "
            "code, so it has no unique position in the declared range"
        )
    pos = num - first_pos + 1
    if not 1 <= pos <= len(declared_names):
        raise ValueError(
            f"chain {src}: modelled residue {num} lies outside the declared range "
            f"{first_pos}-{first_pos + len(declared_names) - 1}"
        )
    if declared_names[pos - 1] != res.name:
        raise ValueError(
            f"chain {src}: modelled residue {num} is {res.name} but the declared "
            f"sequence has {declared_names[pos - 1]} at position {pos} (author {num}) "
            "-- the coordinates and _pdbx_poly_seq_scheme disagree about the "
            "numbering"
        )
    return pos


def ranges_str(segments: list[tuple[int, int]]) -> str:
    return "+".join(f"{lo}-{hi}" for lo, hi in segments)


def set_entities(
    st: gemmi.Structure,
    declared: dict[str, list[str]] | None = None,
    assign: bool = True,
) -> None:
    """One polymer entity per chain, ``full_sequence`` = the modelled residues.

    Under ``--keep-gaps`` the caller passes ``declared`` (the full range, per
    output chain) and ``assign=False``, because ``label_seq_id`` has already
    been set to the position in that range and ``assign_label_seq_id()`` would
    compact it back to 1..N. Both arguments default to the contiguous
    behaviour, which is unchanged.

    ``setup_entities()`` is deliberately NOT used. It would relabel subchains
    (``Axp``) and deduplicate sequence-identical chains into a shared entity
    whose ``full_sequence`` comes from one of them -- the two ways this file
    goes wrong. Duplicating identical entities is redundant mmCIF, not invalid,
    and it makes the per-chain token count impossible to get wrong.

    Doing it by hand means the two things ``setup_entities()`` would have set on
    each RESIDUE must be set here, because ``assign_label_seq_id()`` resolves a
    residue's entity through them and silently skips any residue it cannot
    resolve, leaving ``label_seq_id`` as ``.`` in the written file:

    * ``res.subchain`` -- the key looked up in ``Entity.subchains``. Also set in
      ``select()``, so no source label survives.
    * ``res.entity_type`` -- gemmi only numbers residues marked ``Polymer``. It
      is inherited from the source file, and a plain PDB or a bare ``_atom_site``
      mmCIF (anything a predictor or a stripped-down target file emits) carries
      no entity annotation at all, so it arrives as ``Unknown`` and every
      residue is skipped. Under ``--keep-het`` this also marks kept
      ligands/waters as polymer residues, which is what ``full_sequence``
      already declares them to be.
    """
    st.entities.clear()
    for chain in st[0]:
        for res in chain:
            res.subchain = chain.name
            res.entity_type = gemmi.EntityType.Polymer
        ent = gemmi.Entity(chain.name)
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = gemmi.PolymerType.PeptideL
        ent.subchains = [chain.name]
        if declared is None:
            ent.full_sequence = [res.name for res in chain]
        else:
            if chain.name not in declared:
                raise ValueError(f"no declared sequence for output chain {chain.name}")
            ent.full_sequence = list(declared[chain.name])
        st.entities.append(ent)
    # label_seq_id then runs 1..N per chain, matching the token indices Boltz
    # builds from _entity_poly_seq.
    if assign:
        st.assign_label_seq_id(True)
    # Never write a `.` and hope: an unnumbered residue is the failure above.
    unnumbered = [
        f"{chain.name}:{res.seqid.num}{res.name}"
        for chain in st[0]
        for res in chain
        if res.label_seq is None
    ]
    if unnumbered:
        raise ValueError(
            f"{len(unnumbered)} residue(s) got no label_seq_id "
            f"(first: {', '.join(unnumbered[:4])}) -- they would be written as "
            "'.' and Boltz's token indices would not line up with "
            "_entity_poly_seq"
        )


def write_cif(st: gemmi.Structure, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    doc = st.make_mmcif_document()
    doc.write_file(str(out))


# --------------------------------------------------------------------------
# verify -- reads the written file back as raw mmCIF
# --------------------------------------------------------------------------


def gap_shaped(nums: list[str], modelled: list[str], mon_ids: list[str]) -> bool:
    """Does this chain look like a --keep-gaps chain, with no benefit of doubt?

    True only when the label_seq_ids are integers, strictly ascending, inside
    1..declared, the residue names match the declared sequence at exactly those
    positions, AND the numbering is NOT the plain contiguous 1..N. That last
    clause is what keeps the old failure intact: a template that declares more
    residues than it models while numbering them 1..M is the silent shift this
    file has always refused, and it must keep failing unless the caller states
    the declared length with --expect-declared.
    """
    if not nums or not all(n.isdigit() for n in nums):
        return False
    values = [int(n) for n in nums]
    if values == list(range(1, len(values) + 1)):
        return False
    if any(b <= a for a, b in zip(values, values[1:])):
        return False
    if values[0] < 1 or values[-1] > len(mon_ids):
        return False
    return all(mon_ids[v - 1] == comp for v, comp in zip(values, modelled))


def check_gap_chain(
    label_asym: str, entity: str, nums: list[str], modelled: list[str], mon_ids: list[str]
) -> tuple[list[str], list[str]]:
    """Validate one gapped chain. Returns (ok lines, failures)."""
    lines: list[str] = []
    failures: list[str] = []
    if not all(n.isdigit() for n in nums):
        failures.append(
            f"chain {label_asym} (entity {entity}): label_seq_id is not an integer "
            f"for every residue (starts {','.join(nums[:4])})"
        )
        return lines, failures
    values = [int(n) for n in nums]
    if any(b <= a for a, b in zip(values, values[1:])):
        failures.append(
            f"chain {label_asym} (entity {entity}): label_seq_id is not unique and "
            f"ascending (starts {','.join(nums[:6])})"
        )
    if values[0] < 1 or values[-1] > len(mon_ids):
        failures.append(
            f"chain {label_asym} (entity {entity}): label_seq_id runs "
            f"{values[0]}..{values[-1]}, outside 1..{len(mon_ids)} declared"
        )
    if len(modelled) > len(mon_ids):
        failures.append(
            f"chain {label_asym} (entity {entity}): {len(modelled)} residues are "
            f"modelled but only {len(mon_ids)} are declared"
        )
    mismatch = [
        (v, comp)
        for v, comp in zip(values, modelled)
        if 1 <= v <= len(mon_ids) and mon_ids[v - 1] != comp
    ]
    if mismatch:
        pos, comp = mismatch[0]
        failures.append(
            f"chain {label_asym} (entity {entity}): modelled residue at label_seq_id "
            f"{pos} is {comp} but _entity_poly_seq declares {mon_ids[pos - 1]} there "
            f"({len(mismatch)} such position(s)) -- the template would hand Boltz the "
            "wrong residue at that index"
        )
    if failures:
        return lines, failures
    unresolved = sorted(set(range(1, len(mon_ids) + 1)) - set(values))
    lines.append(
        f"OK  chain {label_asym}: _entity_poly_seq = {len(mon_ids)} residues, "
        f"{len(modelled)} modelled at label_seq_id {values[0]}..{values[-1]} "
        "(keep-gaps: declared complete, atoms only where resolved)"
    )
    lines.append(
        f"    chain {label_asym}: {len(unresolved)} unresolved position(s)"
        + (f": {condense(unresolved)}" if unresolved else "")
        + f"; {len(modelled)}/{len(mon_ids)} modelled "
        f"({100 * len(modelled) / len(mon_ids):.0f}%)"
    )
    return lines, failures


def check_boltz_alignment(
    st: gemmi.Structure, label_asym: str, mon_ids: list[str]
) -> tuple[list[str], list[str]]:
    """Place the modelled residues the way Boltz does, and compare with label_seq_id.

    Boltz's ``parse_polymer`` does NOT read ``label_seq_id``: it aligns the
    declared sequence to the modelled polymer with
    ``gemmi.align_sequence_to_polymer`` (default ``AlignmentScoring``) and walks
    ``match_string`` over the declared sequence, taking the next modelled
    residue only on an identity column ``|``. A mismatch column consumes a
    declared position WITHOUT consuming the residue, which desyncs everything
    after it. This repeats that walk and fails unless every modelled residue
    lands on the declared position equal to its own ``label_seq_id`` -- the
    property the template actually relies on, not a proxy for it.
    """
    polymer = st[0].get_subchain(label_asym)
    if len(polymer) == 0:
        return [], [f"chain {label_asym}: no polymer residues to align"]
    result = gemmi.align_sequence_to_polymer(
        mon_ids, polymer, gemmi.PolymerType.PeptideL, gemmi.AlignmentScoring()
    )
    match = result.match_string
    if len(match) != len(mon_ids):
        return [], [
            f"chain {label_asym}: the alignment Boltz runs has {len(match)} columns "
            f"for {len(mon_ids)} declared residues (cigar {result.cigar_str()}) -- "
            "some modelled residue has no declared position"
        ]
    placed: dict[int, int] = {}
    i = 0
    for j, column in enumerate(match):
        if i < len(polymer) and column == "|":
            placed[i] = j + 1
            i += 1
    wrong = [
        (k, polymer[k].label_seq, placed.get(k))
        for k in range(len(polymer))
        if placed.get(k) != polymer[k].label_seq
    ]
    if wrong:
        k, want, got = wrong[0]
        return [], [
            f"chain {label_asym}: Boltz's own alignment places {len(wrong)} of "
            f"{len(polymer)} modelled residue(s) away from their label_seq_id "
            f"(first: {polymer[k].name} label_seq_id {want} -> "
            f"{'no position' if got is None else f'position {got}'}; cigar "
            f"{result.cigar_str()}) -- the template would be read shifted"
        ]
    return [
        f"OK  chain {label_asym}: Boltz's alignment places all {len(polymer)} "
        "modelled residues at their label_seq_id"
    ], []


def check_declared_sequence(
    label_asym: str,
    mon_ids: list[str],
    query: str,
    allow_mismatches: int = 0,
    author_of=lambda pos: None,
) -> tuple[list[str], list[str]]:
    """Compare an entity's DECLARED sequence with the query, 1:1.

    Boltz maps query to template by a local BLOSUM62 alignment with a gap
    penalty of -1000 -- one gapless diagonal -- between the query chain and the
    DECLARED template sequence, with no identity threshold. So a declared
    sequence of a different length, or one shifted against the query, is mapped
    anyway and silently. A length difference FAILS, and so does any number of
    differing positions above ``allow_mismatches`` (default 0): a shift by k
    differs almost everywhere, while a genuine point difference (an engineered
    residue, a species variant) is a handful that the caller can allow
    explicitly after reading where they are.

    ``query`` must already be normalised (upper case, no whitespace).
    Positions are printed in BOTH numberings -- the 1-based position in the
    declared range and, where the file says, the author number -- so that a
    position is never read as a shift of the author number.
    """
    lines: list[str] = []
    failures: list[str] = []
    declared_seq = "".join(one_letter_name(m) for m in mon_ids)

    def where(pos: int) -> str:
        author = author_of(pos)
        return f"position {pos}" + (f" (author {author})" if author is not None else "")

    if len(declared_seq) != len(query):
        failures.append(
            f"chain {label_asym}: declares {len(declared_seq)} residues but the "
            f"query given with --expect-sequence has {len(query)} -- Boltz aligns "
            "these two gaplessly, so they cannot be made to correspond"
        )
        return lines, failures
    diffs = [i + 1 for i, (a, b) in enumerate(zip(declared_seq, query)) if a != b]
    detail = ", ".join(
        f"{where(p)}: {query[p - 1]}->{declared_seq[p - 1]}" for p in diffs[:8]
    ) + (" ..." if len(diffs) > 8 else "")
    if len(diffs) > allow_mismatches:
        failures.append(
            f"chain {label_asym}: declared sequence differs from the query at "
            f"{len(diffs)}/{len(query)} position(s), more than --allow-mismatches "
            f"{allow_mismatches}: {detail} -- a shifted or wrong template; if these "
            "are genuine point differences, allow exactly that many"
        )
    elif diffs:
        lines.append(
            f"NOTE chain {label_asym}: declared sequence differs from the query at "
            f"{len(diffs)}/{len(query)} position(s) (allowed {allow_mismatches}): "
            f"{detail}"
        )
    else:
        lines.append(
            f"OK  chain {label_asym}: declared sequence is the query, "
            f"all {len(query)} positions"
        )
    odd = nonstandard(mon_ids)
    if odd:
        lines.append(
            f"NOTE chain {label_asym}: {len(odd)} declared residue(s) are not "
            "standard amino acids and compare as their parent letter: "
            + ", ".join(f"{where(p)} {n}" for p, n in odd[:8])
            + (" ..." if len(odd) > 8 else "")
        )
    return lines, failures


def normalise_query(query: str) -> str:
    """Upper case, all whitespace removed -- a pasted FASTA line must not 'differ everywhere'."""
    return "".join(query.split()).upper()


def parse_expect_declared(value: str) -> int | dict[str, int] | None:
    """'152' -> 152 (every entity); 'A:152,B:140' -> per chain; '' -> None."""
    value = value.strip()
    if not value:
        return None
    if ":" not in value:
        return int(value)
    out: dict[str, int] = {}
    for item in split_list(value):
        chain, _, n = item.partition(":")
        out[chain.strip()] = int(n)
    return out


def verify(
    path: Path,
    expect_chains: list[str],
    expect_declared: int | dict[str, int] | None = None,
    expect_sequence: str = "",
    allow_mismatches: int = 0,
    expect_modelled_min: int | None = None,
) -> tuple[bool, list[str]]:
    """Check the invariants Boltz depends on. Returns (ok, report lines).

    Reads ``_atom_site``, ``_entity_poly_seq`` and ``_struct_asym`` directly:
    going through ``gemmi.read_structure`` would re-derive what is being
    checked and pass on a file Boltz rejects.

    ``expect_declared`` is the length of the range a ``--keep-gaps`` build
    declared. Given, every entity must declare exactly that many residues and
    the gapped checks are applied; not given, a chain is only treated as gapped
    when its own numbering proves it (see ``gap_shaped``), so a contiguous file
    that declares more than it models fails exactly as it always has.

    It may also be a {chain: length} dict -- what a build passes, so a chain
    whose range differs in length from its neighbours' is still checked against
    ITS OWN length (and a legitimate chain whose only unresolved residues are the
    C-terminal tail, numbered 1..M contiguously, is recognised). A chain absent
    from the dict falls back to detection.

    ``expect_sequence`` is the query sequence (one-letter, normalised) each
    entity is meant to template, compared position by position with the
    declared sequence; more than ``allow_mismatches`` differences fail.

    ``expect_modelled_min`` fails any chain modelling fewer residues than that
    (a 2-of-152 template is well-formed and nearly useless).

    Every chain is also put through the alignment Boltz itself runs
    (``check_boltz_alignment``). On a gapped chain it prints its own OK line; on
    a plain contiguous chain it is silent unless it fails, so the historic three
    OK lines of an ordinary template are unchanged.
    """
    doc = gemmi.cif.read(str(path))
    block = doc.sole_block()
    lines: list[str] = []
    failures: list[str] = []

    atoms = block.find(
        "_atom_site.",
        [
            "label_asym_id", "auth_asym_id", "label_entity_id", "label_seq_id",
            "label_comp_id", "?auth_seq_id",
        ],
    )
    if not len(atoms):
        return False, ["FAIL: no _atom_site rows"]
    # label_seq_id -> author number per chain, for printing BOTH numberings
    auth_by_label: dict[str, dict[int, int]] = {}
    for row in atoms:
        if row.has(5) and row[3].isdigit():
            text = row.str(5).strip()
            if text.lstrip("-").isdigit():
                auth_by_label.setdefault(row[0], {})[int(row[3])] = int(text)

    # modelled residues per label_asym_id, in file order
    residues: dict[str, list[tuple[str, str]]] = {}
    label_to_auth: dict[str, set[str]] = {}
    label_to_entity: dict[str, set[str]] = {}
    unnumbered: dict[str, int] = {}
    for row in atoms:
        label_asym, auth_asym, entity, label_seq, comp = (row[i] for i in range(5))
        label_to_auth.setdefault(label_asym, set()).add(auth_asym)
        label_to_entity.setdefault(label_asym, set()).add(entity)
        if label_seq in (".", "?"):
            unnumbered[label_asym] = unnumbered.get(label_asym, 0) + 1
        seen = residues.setdefault(label_asym, [])
        if not seen or seen[-1] != (label_seq, comp):
            seen.append((label_seq, comp))

    # 1. label_asym_id == auth_asym_id
    mismatched = {
        lab: sorted(auth) for lab, auth in label_to_auth.items() if auth != {lab}
    }
    if mismatched:
        failures.append(
            "label_asym_id != auth_asym_id for "
            + "; ".join(f"{lab} (auth {','.join(a)})" for lab, a in mismatched.items())
            + " -- Boltz addresses template chains by label_asym_id, so template_id "
            "must name these, or the run dies with 'not one of the protein chains'"
        )
    else:
        lines.append(f"OK  label_asym_id == auth_asym_id for all {len(residues)} chains")

    # 1b. every _atom_site row carries a label_seq_id. Reported before the count
    # check below because it is the CAUSE of a confusing symptom there: residues
    # are deduped on (label_seq_id, comp_id), so with label_seq_id == '.' two
    # ADJACENT residues of the same type collapse into one and the modelled
    # count comes out short by exactly the number of such pairs.
    if unnumbered:
        failures.append(
            "label_seq_id is '.' for "
            + "; ".join(
                f"{lab} ({n} atom rows)" for lab, n in sorted(unnumbered.items())
            )
            + " -- Boltz indexes template residues by label_seq_id, so an "
            "unnumbered residue has no token; the residue counts below are "
            "unreliable for these chains"
        )

    # 2. _entity_poly_seq rows per entity == modelled residues of its subchains
    poly = block.find("_entity_poly_seq.", ["entity_id", "num", "mon_id"])
    declared: dict[str, list[str]] = {}
    for row in poly:
        declared.setdefault(row[0], []).append(row[2])
    if not declared:
        failures.append(
            "no _entity_poly_seq rows at all -- Boltz's parse_polymer builds its "
            "token list from this category, so there is nothing to template"
        )
    try:
        aligned_st: gemmi.Structure | None = gemmi.make_structure_from_block(block)
    except Exception as exc:  # noqa: BLE001 -- reported, never swallowed
        aligned_st = None
        failures.append(f"gemmi could not build a structure from the file: {exc}")

    def author_lookup(label_asym: str):
        known = auth_by_label.get(label_asym, {})
        offsets = {a - l for l, a in known.items()}

        def author_of(pos: int) -> int | None:
            if pos in known:
                return known[pos]
            return pos + next(iter(offsets)) if len(offsets) == 1 else None

        return author_of

    gapped: set[str] = set()
    for label_asym, res_list in sorted(residues.items()):
        entities = label_to_entity[label_asym]
        if len(entities) != 1:
            failures.append(
                f"chain {label_asym} spans entities {','.join(sorted(entities))}"
            )
            continue
        entity = entities.pop()
        mon_ids = declared.get(entity, [])
        modelled = [comp for _, comp in res_list]
        nums = [n for n, _ in res_list]
        want = (
            expect_declared.get(label_asym)
            if isinstance(expect_declared, dict)
            else expect_declared
        )
        if expect_sequence:
            seq_lines, seq_failures = check_declared_sequence(
                label_asym,
                mon_ids,
                expect_sequence,
                allow_mismatches,
                author_lookup(label_asym),
            )
            lines.extend(seq_lines)
            failures.extend(seq_failures)
        if expect_modelled_min is not None and len(modelled) < expect_modelled_min:
            failures.append(
                f"chain {label_asym}: only {len(modelled)} residues are modelled, "
                f"fewer than --expect-modelled-min {expect_modelled_min}"
            )
        boltz_lines: list[str] = []
        if aligned_st is not None and label_asym not in unnumbered and mon_ids:
            boltz_lines, boltz_failures = check_boltz_alignment(
                aligned_st, label_asym, mon_ids
            )
            failures.extend(boltz_failures)
        if want is not None and len(mon_ids) != want:
            failures.append(
                f"chain {label_asym} (entity {entity}): _entity_poly_seq declares "
                f"{len(mon_ids)} residues but --expect-declared says {want} "
                "-- the declared sequence is not the range the query carries, so "
                "position i of the template is not position i of the query"
            )
            gapped.add(label_asym)
            continue
        if want is not None or gap_shaped(nums, modelled, mon_ids):
            gapped.add(label_asym)
            ok_lines, gap_failures = check_gap_chain(
                label_asym, entity, nums, modelled, mon_ids
            )
            lines.extend(ok_lines)
            if not gap_failures:
                lines.extend(boltz_lines)
            failures.extend(gap_failures)
            continue
        if len(mon_ids) != len(modelled):
            failures.append(
                f"chain {label_asym} (entity {entity}): _entity_poly_seq declares "
                f"{len(mon_ids)} residues but {len(modelled)} are modelled -- every "
                "declared-but-absent residue still consumes a token index, so all "
                "downstream residue indices shift silently. A shared entity between "
                "two sequence-identical chains is the usual cause"
            )
        elif mon_ids != modelled:
            first = next(
                i for i, (a, b) in enumerate(zip(mon_ids, modelled)) if a != b
            )
            failures.append(
                f"chain {label_asym} (entity {entity}): _entity_poly_seq and the "
                f"modelled residues differ at position {first + 1} "
                f"({mon_ids[first]} declared, {modelled[first]} modelled)"
            )
        else:
            lines.append(
                f"OK  chain {label_asym}: _entity_poly_seq = {len(mon_ids)} residues, "
                "identical to what is modelled"
            )

    # 3. label_seq_id contiguous from 1 (the token index Boltz builds).
    # Skipped for a keep-gaps chain, whose label_seq_id is checked against the
    # declared sequence in check_gap_chain() instead -- gaps are the point there.
    for label_asym, res_list in sorted(residues.items()):
        if label_asym in gapped:
            continue
        nums = [r[0] for r in res_list]
        if nums != [str(i + 1) for i in range(len(nums))]:
            failures.append(
                f"chain {label_asym}: label_seq_id is not 1..{len(nums)} "
                f"(starts {','.join(nums[:4])})"
            )

    # 4. the chains the YAML will name are actually there
    missing = [c for c in expect_chains if c not in residues]
    if missing:
        failures.append(
            f"template_id names {','.join(missing)} but the CIF has "
            f"{','.join(sorted(residues))}"
        )
    elif expect_chains:
        lines.append(f"OK  template chains present: {','.join(expect_chains)}")

    for f in failures:
        lines.append(f"FAIL {f}")
    return not failures, lines


# --------------------------------------------------------------------------
# YAML
# --------------------------------------------------------------------------


def template_yaml(
    cif: str, chain_ids: list[str], template_ids: list[str], force: bool, threshold: float
) -> str:
    """The ``templates:`` block, at column 0 -- boltz splices it verbatim.

    ``chain_id`` (prediction chains) and ``template_id`` (chains in this CIF)
    are ALWAYS both written and equal length: that is what makes Boltz zip them
    positionally instead of falling back to a sequence search.
    """
    if len(chain_ids) != len(template_ids):
        raise ValueError(
            f"--predict-chains has {len(chain_ids)} IDs but the template has "
            f"{len(template_ids)} chains -- they must zip 1:1"
        )
    body = [
        "templates:",
        f"  - cif: {cif}",
        f"    chain_id: [{', '.join(chain_ids)}]",
        f"    template_id: [{', '.join(template_ids)}]",
    ]
    if force:
        body.append("    force: true")
        body.append(f"    threshold: {threshold}")
    return "\n".join(body) + "\n"


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


class Args(argparse.Namespace):
    pass


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Build a Boltz forced-template mmCIF and its templates: YAML.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  # 8-chain template from two protomers of a fetched assembly\n"
            "  boltz_template.py --fetch 7OJG --out /runs/inputs/tmpl.cif \\\n"
            "      --chains K,A,B,C --residues K:19-59+107-155 \\\n"
            "      --yaml /runs/inputs/tmpl.yaml\n\n"
            "  # verify a CIF written by something else\n"
            "  boltz_template.py --check /runs/inputs/tmpl.cif --expect-chains A,B,C,D\n"
        ),
    )
    ap.add_argument("--in", dest="src", type=Path, help="Source PDB/mmCIF (or --fetch).")
    ap.add_argument("--fetch", type=str, default="", help="PDB ID to download from RCSB.")
    ap.add_argument(
        "--fetch-to",
        type=Path,
        help="Where to keep the fetched entry (default: <out>'s dir, <ID>.cif).",
    )
    ap.add_argument("--out", type=Path, help="Template mmCIF to write.")
    ap.add_argument("--chains", type=str, default="", help="Source chain IDs, in order.")
    ap.add_argument(
        "--rename-to", type=str, default="", help="Output chain IDs, positional."
    )
    ap.add_argument(
        "--residues",
        type=str,
        default="",
        help="Per-SOURCE-chain residue ranges, e.g. 'K:19-59+107-155,A:36-149'.",
    )
    ap.add_argument(
        "--predict-chains",
        type=str,
        default="",
        help="chain_id in the YAML (chains of the PREDICTION). Default: the "
        "template's own chain IDs.",
    )
    ap.add_argument("--yaml", dest="yaml_out", type=Path, help="Write the YAML block here.")
    ap.add_argument(
        "--cif-path",
        type=str,
        default="",
        help="Path to name in the YAML's cif: field (default: --out, which is "
        "already a workstation path).",
    )
    ap.add_argument(
        "--threshold", type=float, default=2.0, help="force threshold (default 2.0)."
    )
    ap.add_argument(
        "--no-force",
        action="store_true",
        help="Omit force/threshold (template as a hint, not forced).",
    )
    ap.add_argument("--renumber", action="store_true", help="Renumber each chain 1..N.")
    ap.add_argument(
        "--keep-gaps",
        action="store_true",
        help="Declare the WHOLE --residues range of each chain (names for "
        "unresolved positions from the source's _pdbx_poly_seq_scheme) and give "
        "each modelled residue its position in that range as label_seq_id, so "
        "the template stays 1:1 with a query that has no gaps. Needs exactly one "
        "range per chain. No coordinates are invented.",
    )
    ap.add_argument(
        "--keep-het",
        action="store_true",
        help="Keep ligands/ions/waters (default: polymer residues only).",
    )
    ap.add_argument("--model", type=int, default=0, help="Model index (default 0).")
    ap.add_argument("--check", type=Path, help="Verify an existing CIF and exit.")
    ap.add_argument(
        "--expect-chains", type=str, default="", help="With --check: chains that must exist."
    )
    ap.add_argument(
        "--expect-declared",
        type=str,
        default="",
        help="With --check: the number of residues every entity must DECLARE "
        "(the length of a --keep-gaps range), e.g. '152', or per chain "
        "'A:152,B:140'. Also tells the checker the file is a keep-gaps file, "
        "which is the only way to recognise one whose modelled residues happen "
        "to start at position 1 and run contiguously.",
    )
    ap.add_argument(
        "--expect-sequence",
        type=str,
        default="",
        help="The query sequence (one-letter) this template is for. Every entity's "
        "DECLARED sequence is compared with it position by position; a different "
        "length fails, and so do more differing positions than --allow-mismatches. "
        "Use it on a build or with --check. Boltz aligns query to declared "
        "template gaplessly, so this is the only check that a template is not "
        "silently shifted. Case and whitespace are normalised.",
    )
    ap.add_argument(
        "--allow-mismatches",
        type=int,
        default=0,
        help="With --expect-sequence: how many positions may legitimately differ "
        "from the query (default 0). Read the reported positions first -- e.g. "
        "8Z8M's TNF needs 1 for its real L143D.",
    )
    ap.add_argument(
        "--expect-modelled-min",
        type=int,
        default=None,
        help="With --check: fail any chain modelling fewer than this many residues "
        "(a template with almost no atoms is well-formed and useless).",
    )
    return ap


def main() -> None:
    args = build_parser().parse_args(namespace=Args())

    if args.check:
        if not args.check.exists():
            # Each modal-shell call is a FRESH container: anything written to
            # /tmp by an earlier call is gone. Templates belong on the Volume.
            raise SystemExit(
                f"no such file: {args.check} -- if a previous call wrote it to /tmp, "
                "note every `sapia modal-shell` call is a new container; write "
                "templates under /runs/"
            )
        ok, lines = verify(
            args.check,
            split_list(args.expect_chains),
            parse_expect_declared(args.expect_declared),
            normalise_query(args.expect_sequence),
            args.allow_mismatches,
            args.expect_modelled_min,
        )
        print(f"--- {args.check}")
        print("\n".join(lines))
        raise SystemExit(0 if ok else 1)

    if not args.out:
        raise SystemExit("--out is required (or use --check)")
    if not args.src and not args.fetch:
        raise SystemExit("one of --in / --fetch is required")

    src = args.src
    if args.fetch:
        dest = args.fetch_to or args.out.parent / f"{args.fetch.upper()}.cif"
        src = fetch_rcsb(args.fetch, dest)

    groups = resolve_chains(args.chains, args.rename_to)
    scheme = None
    if args.keep_gaps:
        if args.renumber:
            raise SystemExit(
                "--keep-gaps and --renumber contradict each other: --keep-gaps "
                "numbers each residue by its position in the declared range"
            )
        if args.keep_het:
            raise SystemExit(
                "--keep-gaps and --keep-het contradict each other: the declared "
                "range comes from _pdbx_poly_seq_scheme, which lists polymer "
                "residues only, so a kept ligand or water has no position in it"
            )
        if not args.residues:
            raise SystemExit(
                "--keep-gaps needs --residues '<chain>:<first>-<last>' per chain -- "
                "the declared range is what the query is mapped onto"
            )
        scheme = read_poly_seq_scheme(src, [source for source, _ in groups])
    built = select(
        load_structure(src),
        groups,
        parse_residues(args.residues),
        strip_het=not args.keep_het,
        renumber=args.renumber,
        model_index=args.model,
        keep_gaps=args.keep_gaps,
        scheme=scheme,
    )
    # Write to a temporary beside the target and rename only once verify()
    # passes: under `sapia modal-shell --cmd` the exit code is lost, so an
    # ABSENT output file has to be the failure signal.
    tmp_out = args.out.with_name(args.out.name + ".tmp")
    write_cif(built.structure, tmp_out)

    out_ids = [c.out for c in built.chains]
    predict = split_list(args.predict_chains) or out_ids
    block = template_yaml(
        cif=args.cif_path or str(args.out),
        chain_ids=predict,
        template_ids=out_ids,
        force=not args.no_force,
        threshold=args.threshold,
    )

    print(f"--- {args.out}")
    for c in built.chains:
        seq = c.sequence if len(c.sequence) <= 60 else c.sequence[:57] + "..."
        print(
            f"  {c.source} -> {c.out}: {c.n_res} res, {c.first}-{c.last}  {seq}"
        )
        if args.keep_gaps:
            print(
                f"      declared {c.n_declared} "
                f"({c.offset + 1}-{c.offset + c.n_declared}), "
                f"{len(c.unresolved)} unresolved"
                + (f" at author {condense(c.unresolved)}" if c.unresolved else "")
            )
    # A single --expect-declared only makes sense when every chain's range has
    # the same length; otherwise let verify() recognise the gaps per chain.
    # Per chain, so chains whose ranges differ in length are each checked
    # against their own declared length.
    expect_declared = (
        {c.out: c.n_declared for c in built.chains} if args.keep_gaps else None
    )
    ok, lines = verify(
        tmp_out,
        out_ids,
        expect_declared,
        normalise_query(args.expect_sequence),
        args.allow_mismatches,
        args.expect_modelled_min,
    )
    print("\n".join(f"  {line}" for line in lines))

    tmp_yaml = None
    if args.yaml_out:
        args.yaml_out.parent.mkdir(parents=True, exist_ok=True)
        tmp_yaml = args.yaml_out.with_name(args.yaml_out.name + ".tmp")
        tmp_yaml.write_text(block)
        print(f"--- {args.yaml_out}")
    else:
        print("--- templates: block (pass to boltz --template-yaml)")
    print(block, end="")

    if not ok:
        tmp_out.unlink(missing_ok=True)
        if tmp_yaml is not None:
            tmp_yaml.unlink(missing_ok=True)
        print(
            "\nNOT WRITTEN: the checks above failed, and neither "
            f"{args.out} nor the YAML was created. Fix before spending.",
            file=sys.stderr,
        )
        raise SystemExit(1)

    tmp_out.replace(args.out)
    if tmp_yaml is not None:
        tmp_yaml.replace(args.yaml_out)


if __name__ == "__main__":
    main()
