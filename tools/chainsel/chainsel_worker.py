#!/usr/bin/env python3
"""
Extract a named subset of chains from one structure into a new structure file,
optionally merging several source chains into one output chain.

This is the per-array-task step of the chainsel tool (and works standalone for a
single structure). It exists because every structure comparison in the workbench is
whole-file: ``usalign --mm 1`` on a binder complex is dominated by the large fixed
target and reads ~1.0 whatever the binder did, and ``--ter 2`` isolates only the
FIRST chain while a binder is the LAST. Writing the binder chain to its own file
makes binder-to-binder comparison (``usalign --mm 0``) possible.

The second job is the SPLIT PROTOMER: a target trimmed of a membrane belt is given
to a predictor as two entities per protomer, so the prediction holds one protomer as
two chains. ``--groups 'A+C:A,B+D:B'`` writes those back as one chain per protomer,
which is what ``ringfit --target-chains A,B`` needs before it can report
``bsa_t1``/``bsa_t2``/``bridge_ratio`` at all.

The one rule that matters: **a requested chain that is not in the structure is an
error, not a shorter file.** A silently partial extraction would produce a perfectly
well-formed PDB that every downstream number (TM-score, RMSD, BSA, energy) would then
be computed on -- wrong, plausible, and undetectable. So a missing chain is reported
as data (``status`` = ``error: chain E not in <file> (have: A,B,C)``) and no output
file is written. That holds for every chain named in a merge group too.

RESIDUE ORDER inside a merged chain is the order the group names its sources -- all
of A's residues, then all of C's -- and is never re-sorted by residue number: two
segments of a split protomer carry independent numbering from the predictor (each
1..N), so sorting would interleave them catastrophically.

NUMBERING is exactly one of three regimes:

* default: source residue numbers and insertion codes are preserved. If a MERGED
  chain then holds a duplicated (number, icode), that design is an ERROR -- a file
  whose residue numbers repeat inside one chain silently loses residues in any
  consumer that indexes by number (``ringfit --resnum-match resnum`` does);
* ``--renumber``: each OUTPUT chain is renumbered 1..N in written order, insertion
  codes dropped -- continuous across a merge seam;
* ``--renumber-from 'A:19-59,C:107-155'``: each SOURCE chain's segment is renumbered
  consecutively from its own first number, insertion codes dropped. When a ``-last``
  is given it is checked against the residues actually written, so a segment of the
  wrong length is an error rather than a quietly shifted numbering. It assumes each
  segment is internally gapless.

SPLIT MODE (``--split-groups``) is the inverse of a merge: it CUTS one source chain
into several output chains. RFdiffusion's symmetric mode writes one chain per subunit
that fuses the designed binder and the fixed target, with residue numbers that
REPEAT inside the chain (binder 1..60, then target 1..152). Downstream tools need the
binder and the target as separate chains. Segments are addressed by 1-BASED POSITION
in the chain's residue order (counted over the residues chainsel keeps -- polymer
residues only unless ``--keep-het``), never by residue number, because the numbers
repeat and cannot address anything. A position is a residue as delimited in the file
text (see below), so MICROHETEROGENEITY -- two alternative residues of different
names at one site, each with its own altloc -- counts as TWO positions, and both are
written (each keeps its first conformer). Hydrogens and alternative conformations are
dropped within each residue; no residue is ever removed by that, and one left with no
atoms (hydrogens only) is an error. The segments of each source chain must tile it
exactly (from position 1, no gap, no overlap, ending at the chain's real residue
count) -- a checksum: a mismatch is an error for that design, never a trimmed cut.
Output chains are written in ascending output-ID order, each renumbered 1..N
(insertion codes dropped). Split mode always renumbers. It also reports, per output
chain, the residue count and the mean B-factor (``chain_n_res``/``chain_mean_b``),
which is how a caller confirms which segment is the binder (RFdiffusion writes B 0 on
diffused residues and B 1 on fixed ones).

Split mode cannot use gemmi's ordinary reader as-is, and that is the whole hazard of
this mode. Measured with gemmi 0.7.5: the PDB and mmCIF readers merge a residue into
an EARLIER residue of the same chain with the same (number, insertion code, name) --
adjacent or not -- so two GLY 2 come back as ONE residue with 8 backbone atoms; and
``remove_alternative_conformations`` then DELETES every later residue whose number
repeats an earlier one (it treats them as microheterogeneity). On an RFdiffusion
symmetric chain that silently deletes target positions 61..120 and folds target GLYs
into binder GLYs. Split mode therefore first rewrites every residue's number in the
source TEXT to its position (a residue boundary is a change of chain / number /
insertion code / name from the previous atom record, OR an atom name+altloc already
seen in the current residue), so gemmi only ever sees unique numbers; and it checks
that gemmi returned exactly the residues it delimited, in order. The other modes
still use the ordinary reader, but now REFUSE a kept chain whose numbering repeats
(check_numbering_not_repeated: a residue holding the same atom twice, or a number
shared by residues one of which has an atom without altloc), with an error pointing
to --split-groups; true microheterogeneity is exempt. On every input the guard
accepts, their output is byte-for-byte what it was.

Reads PDB or mmCIF transparently (gemmi). Alternative conformations and hydrogens are
dropped on read, so the extraction is deterministic regardless of what the source
modelled. Only the FIRST model is used: a multi-model file (NMR ensemble, some
predictors' output) would otherwise write a chain several times over.

``--strip-het`` (the default; ``--keep-het`` disables it) keeps only residues gemmi
tabulates as amino acids or nucleic acids, which drops waters, ions and ligands while
keeping modified polymer residues such as MSE.

Writes a one-row TSV (name, status, path, n_chains, n_res, chains, n_atoms -- plus
chain_n_res, chain_mean_b in split mode ONLY, so the other modes' TSV is byte-for-byte
what it was before split mode existed) that collect_chainsel.py merges back into the
table. ``n_chains`` counts OUTPUT chains (so
a merge lowers it), while ``n_res``/``n_atoms`` count everything written (so a merge
does not change them). Errors are recorded as data (a status starting with 'error:')
rather than only crashing, so partial array runs still collect.

Normally invoked per array task by the chainsel tool -- run ``sapia run chainsel``
rather than calling this directly.

Usage (standalone, single structure):
    python tools/chainsel/chainsel_worker.py \\
        --name design_0 --src design_0.cif --chains E --rename-to A \\
        --out design_0_binder.pdb --result-tsv design_0.tsv

    python tools/chainsel/chainsel_worker.py \\
        --name design_0 --src design_0.cif \\
        --groups 'A+C:A,B+D:B,E:E' \\
        --renumber-from 'A:19-59,C:107-155,B:19-59,D:107-155' \\
        --out design_0_merged.pdb --result-tsv design_0.tsv

    python tools/chainsel/chainsel_worker.py \\
        --name design_0 --src design_0.pdb \\
        --split-groups 'A:1-60>D+61-212>A,B:1-60>E+61-212>B,C:1-60>F+61-212>C' \\
        --out design_0_split.pdb --result-tsv design_0.tsv
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

import gemmi

# Metric columns of the per-design TSV, in order (bare names: collect_chainsel.py
# hands them to the driver, which leaf-prefixes them to chainsel_<name>).
METRIC_COLUMNS = ["n_chains", "n_res", "chains", "n_atoms"]
RESULT_COLUMNS = ["name", "status", "path", *METRIC_COLUMNS]
# Written by split mode only (appended after METRIC_COLUMNS), so a non-split TSV is
# unchanged. Per output chain, in output-ID order: 'A=152,...,D=60' and
# 'A=1.00,...,D=0.00'.
SPLIT_METRIC_COLUMNS = ["chain_n_res", "chain_mean_b"]

# One output chain: source chain IDs in write order, and the ID written.
Group = tuple[list[str], str]
# One split source chain: its ID and its segments (first, last, out_id), 1-based
# inclusive POSITIONS in the chain's residue order, listed in order.
Segment = tuple[int, int, str]
Split = tuple[str, list[Segment]]


def split_list(value: str) -> list[str]:
    """Comma-joined list -> stripped, non-empty tokens."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


def parse_groups(value: str) -> list[Group]:
    """'A+C:A,B:B' -> [(['A', 'C'], 'A'), (['B'], 'B')]."""
    groups: list[Group] = []
    for spec in split_list(value):
        if spec.count(":") != 1:
            raise ValueError(
                f"--groups entry {spec!r} must be '<chain>[+<chain>...]:<out_id>'"
            )
        sources_str, out_id = (tok.strip() for tok in spec.split(":"))
        sources = [tok.strip() for tok in sources_str.split("+") if tok.strip()]
        if not sources or not out_id:
            raise ValueError(f"--groups entry {spec!r} is missing a chain or an id")
        groups.append((sources, out_id))
    if not groups:
        raise ValueError("--groups is empty")
    return groups


def resolve_groups(groups: str, chains: str, rename_to: str) -> list[Group]:
    """The selection as output groups.

    ``--groups`` is what chainsel.sh passes (the manifest carries the canonical
    form). ``--chains``/``--rename-to`` are kept for standalone use and become
    groups of one.
    """
    if groups.strip():
        return parse_groups(groups)
    kept = split_list(chains)
    renamed = split_list(rename_to)
    if not kept:
        raise ValueError("--chains is empty")
    if renamed and len(renamed) != len(kept):
        raise ValueError(
            f"--rename-to has {len(renamed)} IDs but --chains has {len(kept)}"
        )
    out_ids = renamed or kept
    return [([chain], out_id) for chain, out_id in zip(kept, out_ids)]


def parse_starts(value: str) -> dict[str, tuple[int, int]]:
    """'A:19-59,C:107' -> {'A': (19, 59), 'C': (107, 0)} (0 = no end given)."""
    starts: dict[str, tuple[int, int]] = {}
    for spec in split_list(value):
        if spec.count(":") != 1:
            raise ValueError(
                f"--renumber-from entry {spec!r} must be '<chain>:<first>[-<last>]'"
            )
        chain, rng = (tok.strip() for tok in spec.split(":"))
        first_str, _, last_str = rng.partition("-")
        first = int(first_str)
        last = int(last_str) if last_str else 0
        starts[chain] = (first, last)
    return starts


def parse_split_groups(value: str) -> list[Split]:
    """'A:1-60>D+61-212>A,B:1-60>E' -> [('A', [(1, 60, 'D'), (61, 212, 'A')]), ...].

    Comma separates SOURCE chains, '+' separates one chain's segments, and each
    segment is '<first>-<last>><out_id>' in 1-based POSITIONS (not residue numbers).
    Everything that can be checked without the structure is checked here: the
    grammar, a source chain named twice, a duplicated output ID, and the tiling --
    each chain's segments must start at position 1 and run on contiguously, in
    order, with no gap and no overlap. Whether the last segment ends exactly at the
    chain's real residue count is a per-structure check, made in extract_split().

    Kept in step with run_chainsel._parse_split_groups (the submit-side copy, which
    refuses the run before anything is queued); tests/test_chainsel.py checks that
    the two accept and reject the same specs.
    """
    splits: list[Split] = []
    seen_chains: set[str] = set()
    out_ids: list[str] = []
    entries = split_list(value)
    if not entries:
        raise ValueError("--split-groups is empty: give e.g. 'A:1-60>D+61-212>A'")
    for entry in entries:
        if entry.count(":") != 1:
            raise ValueError(
                f"--split-groups entry {entry!r} is malformed: write "
                f"'<chain>:<first>-<last>><out_id>[+<first>-<last>><out_id>...]'"
            )
        chain, segs_str = (tok.strip() for tok in entry.split(":"))
        if not chain:
            raise ValueError(f"--split-groups entry {entry!r} names no source chain")
        if chain in seen_chains:
            raise ValueError(f"--split-groups names source chain {chain} twice")
        seen_chains.add(chain)
        segments: list[Segment] = []
        for tok in segs_str.split("+"):
            tok = tok.strip()
            if tok.count(">") != 1:
                raise ValueError(
                    f"--split-groups segment {tok!r} of {entry!r} is malformed: "
                    f"write '<first>-<last>><out_id>', e.g. '1-60>D'"
                )
            rng, out_id = (part.strip() for part in tok.split(">"))
            if not out_id:
                raise ValueError(
                    f"--split-groups segment {tok!r} of {entry!r} names no output "
                    f"chain ID after '>'"
                )
            first_str, sep, last_str = rng.partition("-")
            if not sep:
                raise ValueError(
                    f"--split-groups segment {tok!r} of {entry!r} needs a position "
                    f"range '<first>-<last>' before '>'"
                )
            try:
                first, last = int(first_str), int(last_str)
            except ValueError:
                raise ValueError(
                    f"--split-groups segment {tok!r} of {entry!r} has a non-integer "
                    f"position"
                ) from None
            if first < 1:
                raise ValueError(
                    f"--split-groups segment {tok!r} of {entry!r}: positions are "
                    f"1-based"
                )
            if last < first:
                raise ValueError(
                    f"--split-groups segment {tok!r} of {entry!r} ends ({last}) "
                    f"before it starts ({first})"
                )
            segments.append((first, last, out_id))
            out_ids.append(out_id)
        expected = 1
        for first, last, out_id in segments:
            if first < expected:
                raise ValueError(
                    f"--split-groups chain {chain}: segment {first}-{last}>{out_id} "
                    f"overlaps the previous one (next free position is {expected})"
                )
            if first > expected:
                what = "the chain start" if expected == 1 else "the previous segment"
                raise ValueError(
                    f"--split-groups chain {chain}: gap -- segment "
                    f"{first}-{last}>{out_id} starts at {first} but {what} leaves "
                    f"position {expected} uncovered (segments must tile the chain "
                    f"from position 1, in order)"
                )
            expected = last + 1
        splits.append((chain, segments))
    dups = sorted({o for o in out_ids if out_ids.count(o) > 1})
    if dups:
        raise ValueError(
            f"--split-groups writes output chain ID(s) {dups} more than once; each "
            f"segment needs its own output chain ID"
        )
    return splits


def split_spec(splits: list[Split]) -> str:
    """Canonical text form of a parsed split spec (whitespace normalised)."""
    return ",".join(
        f"{chain}:" + "+".join(f"{a}-{b}>{o}" for a, b, o in segments)
        for chain, segments in splits
    )


def load_structure(
    path: Path, guard_chains: list[str] | None = None, strip_het: bool = True
) -> gemmi.Structure:
    """Read a PDB or mmCIF file into a single-conformer, hydrogen-free structure.

    With ``guard_chains`` (what the non-split modes pass: every source chain they
    keep), first refuse a chain whose numbering repeats -- see
    check_numbering_not_repeated(). The guard only ever raises; on any input it
    accepts, the structure returned is exactly what it was without it.
    """
    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    st = gemmi.read_structure(str(path))
    if guard_chains and len(st) > 0:
        check_numbering_not_repeated(st[0], guard_chains, strip_het, path)
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    if len(st) == 0:
        raise ValueError(f"no model in {path}")
    return st


def check_numbering_not_repeated(
    model: gemmi.Model, chain_ids: list[str], strip_het: bool, path: Path
) -> None:
    """Refuse (raise) a kept chain whose residue numbering repeats inside it.

    Measured with gemmi 0.7.5: the reader merges a residue into an earlier one of
    the same chain with the same (number, insertion code, name), and
    remove_alternative_conformations() then deletes every later residue whose
    number repeats an earlier one. The non-split modes would write such a chain
    short and say OK. Two signatures, checked on the residues that would be written
    (polymer only unless --keep-het), BEFORE alt-confs are dropped:

    (a) a POLYMER residue holding the same (heavy-atom name, altloc) twice -- two
        residues merged on read. Hydrogens and non-polymer residues are not looked
        at: a ligand whose carbons are all named 'C' (docking/OpenBabel output) or
        hydrogens with truncated, repeated names are legitimate and were handled
        before the guard existed. Two merged polymer residues always repeat their
        heavy backbone atoms, so this misses no fused-chain collision.
    (b) a (number, insertion code) shared by several residues where EITHER two or
        more of them have an atom with no altloc, OR two of them use the same altloc
        letter. Microheterogeneity -- alternative residues at one site, each with
        its own altloc letter(s), at most one of them carrying blank-altloc atoms
        (e.g. ARG 12 with a blank backbone and an altloc-A side chain, plus LYS 12
        with an altloc-B side chain) -- is exempt and keeps its old behaviour. A
        fused-chain collision always fires: binder and target residue both have
        blank-altloc atoms.

    Not caught: a STRIPPED residue (e.g. a water) that shares a number with a
    polymer residue and comes BEFORE it in the file -- the reader then drops the
    polymer residue. Polymer-first is the PDB/mmCIF convention.
    """
    for chain_id in chain_ids:
        by_seqid: dict[tuple[int, str], list[gemmi.Residue]] = {}
        for res in chain_residues(model, chain_id, strip_het):
            if is_polymer_residue(res):
                atoms = [
                    (atom.name, atom.altloc) for atom in res if not atom.is_hydrogen()
                ]
                if len(set(atoms)) != len(atoms):
                    raise ValueError(
                        f"numbering repeats inside chain {chain_id} of {path.name}: "
                        f"residue {res.name} "
                        f"{res.seqid.num}{res.seqid.icode.strip()} holds the same "
                        f"atom twice (two residues merged on read); this mode would "
                        f"write a wrong chain -- use --split-groups"
                    )
            key = (res.seqid.num, res.seqid.icode)
            by_seqid.setdefault(key, []).append(res)
        for (num, icode), group in by_seqid.items():
            if len(group) < 2:
                continue
            with_blank = sum(
                1 for res in group if any(not atom.has_altloc() for atom in res)
            )
            letter_sets = [
                {atom.altloc for atom in res if atom.has_altloc()} for res in group
            ]
            letters = [c for letter_set in letter_sets for c in letter_set]
            if with_blank >= 2 or len(set(letters)) != len(letters):
                names = "/".join(res.name for res in group)
                raise ValueError(
                    f"numbering repeats inside chain {chain_id} of {path.name}: "
                    f"{len(group)} residues ({names}) are numbered "
                    f"{num}{icode.strip()}; the reader would drop all but one -- "
                    f"use --split-groups"
                )


def _is_cif(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return name.endswith((".cif", ".mmcif"))


def _read_text(path: Path) -> str:
    if path.name.lower().endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return f.read()
    return path.read_text()


def _positional_pdb(text: str) -> tuple[str, dict[str, int]]:
    """Rewrite each residue's number (cols 23-27) to its 1-based position in its
    chain, in the FIRST model; return the text and the residue count per chain.

    A residue boundary is a change of (chain, name, number, insertion code) from the
    previous atom record, or an (atom name, altloc) already seen in the current
    residue -- so two adjacent 'GLY 2' records become two residues, not one with
    eight atoms. ANISOU records follow the atom before them. Records after the first
    ENDMDL are left alone (only the first model is used).
    """
    out: list[str] = []
    counts: dict[str, int] = {}
    prev_key: tuple[str, str, str] | None = None
    seen: set[tuple[str, str]] = set()
    current: dict[str, int] = {}
    done = False
    for line in text.splitlines(keepends=True):
        record = line[:6]
        if not done and record in ("ATOM  ", "HETATM", "ANISOU"):
            body = line.rstrip("\r\n")
            if len(body) < 27:
                raise ValueError(f"truncated {record.strip()} record: {body!r}")
            chain = line[20:22].strip()
            if record == "ANISOU":
                if chain not in current:
                    raise ValueError(f"ANISOU before any atom of chain {chain!r}")
                num = current[chain]
            else:
                key = (chain, line[17:21], line[22:27])
                atom = (line[12:16], line[16])
                if key != prev_key or atom in seen:
                    counts[chain] = counts.get(chain, 0) + 1
                    seen = set()
                prev_key = key
                seen.add(atom)
                num = current[chain] = counts[chain]
            if num > 9999:
                raise ValueError(
                    f"chain {chain} has more than 9999 residues; PDB columns cannot "
                    f"hold its positions -- use an mmCIF source"
                )
            line = f"{line[:22]}{num:4d} {line[27:]}"
        elif record == "ENDMDL" and counts:
            done = True
        out.append(line)
    return "".join(out), counts


def _positional_cif(doc: gemmi.cif.Document) -> dict[str, int]:
    """Same as _positional_pdb, in place on an mmCIF document's _atom_site loop:
    auth_seq_id (label_seq_id when there is no auth column) becomes the residue's
    position in its auth chain, insertion codes are cleared. Measured with gemmi
    0.7.5: residues are grouped on auth_seq_id, not label_seq_id."""
    block = doc.sole_block()

    def col(*tags: str):
        for tag in tags:
            values = block.find_values(f"_atom_site.{tag}")
            if len(values):
                return values
        return None

    seq = col("auth_seq_id", "label_seq_id")
    chain_col = col("auth_asym_id", "label_asym_id")
    comp = col("auth_comp_id", "label_comp_id")
    atom_col = col("auth_atom_id", "label_atom_id")
    if seq is None or chain_col is None or comp is None or atom_col is None:
        raise ValueError("mmCIF has no usable _atom_site loop")
    icode = col("pdbx_PDB_ins_code")
    alt = col("label_alt_id")
    model = col("pdbx_PDB_model_num")
    label_asym = col("label_asym_id")

    counts: dict[str, int] = {}
    prev_key: tuple[str, ...] | None = None
    seen: set[tuple[str, str]] = set()
    first_model = model[0] if model is not None else None
    for i in range(len(seq)):
        if model is not None and model[i] != first_model:
            continue
        chain = gemmi.cif.as_string(chain_col[i])
        key = (
            chain,
            comp[i],
            seq[i],
            icode[i] if icode is not None else "",
            label_asym[i] if label_asym is not None else "",
        )
        atom = (atom_col[i], alt[i] if alt is not None else "")
        if key != prev_key or atom in seen:
            counts[chain] = counts.get(chain, 0) + 1
            seen = set()
        prev_key = key
        seen.add(atom)
        seq[i] = str(counts[chain])
        if icode is not None:
            icode[i] = "?"
    return counts


def load_structure_positional(path: Path) -> gemmi.Structure:
    """load_structure() for split mode: every residue of the first model carries
    its 1-based POSITION in its chain as its number, so nothing downstream -- gemmi's
    reader, remove_alternative_conformations, a chain-part merge -- can merge,
    delete or reorder residues whose source numbers repeat.

    Then checks it: per chain, gemmi must hand back exactly positions 1..N in order.
    Any other answer means gemmi regrouped the residues, and the cut would be wrong.
    """
    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    if _is_cif(path):
        doc = gemmi.cif.read(str(path))
        counts = _positional_cif(doc)
        st = gemmi.make_structure_from_block(doc.sole_block())
    else:
        text, counts = _positional_pdb(_read_text(path))
        st = gemmi.read_pdb_string(text)
    st.merge_chain_parts()
    if len(st) == 0:
        raise ValueError(f"no model in {path}")

    def numbers(model: gemmi.Model) -> dict[str, list[int]]:
        by_chain: dict[str, list[int]] = {}
        for chain in model:
            by_chain.setdefault(chain.name, []).extend(r.seqid.num for r in chain)
        return by_chain

    for chain_id, nums in numbers(st[0]).items():
        if nums != list(range(1, counts.get(chain_id, 0) + 1)):
            raise ValueError(
                f"chain {chain_id} of {path.name}: the reader returned {len(nums)} "
                f"residues for {counts.get(chain_id, 0)} delimited, or out of order "
                f"-- refusing to cut by position"
            )

    st.remove_alternative_conformations()
    st.remove_hydrogens()
    # Neither call removes a residue (gemmi 0.7.5: a residue made only of
    # hydrogens stays in place with ZERO atoms; one made only of altloc-B atoms
    # keeps them) -- extract_split() refuses zero-atom residues. None may merge or
    # move either; check that here.
    for chain_id, nums in numbers(st[0]).items():
        if any(b <= a for a, b in zip(nums, nums[1:])):
            raise ValueError(
                f"chain {chain_id} of {path.name}: residue order changed while "
                f"dropping alternative conformations/hydrogens"
            )
    return st


def is_polymer_residue(res: gemmi.Residue) -> bool:
    """True for a tabulated amino-acid or nucleic-acid residue.

    Used as the --strip-het test instead of the HETATM flag, which would throw away
    modified polymer residues (MSE and friends are HETATM records but are part of
    the chain).
    """
    info = gemmi.find_tabulated_residue(res.name)
    return bool(info and (info.is_amino_acid() or info.is_nucleic_acid()))


def chain_residues(
    model: gemmi.Model, chain_id: str, strip_het: bool
) -> list[gemmi.Residue]:
    """Residues of every chain object named ``chain_id``, in file order.

    A PDB split by TER, or an mmCIF whose waters hang off the polymer's auth chain,
    can present one chain name as several chain objects; they are concatenated here
    rather than the first one silently winning.
    """
    residues: list[gemmi.Residue] = []
    for chain in model:
        if chain.name != chain_id:
            continue
        for res in chain:
            if strip_het and not is_polymer_residue(res):
                continue
            residues.append(res)
    return residues


def check_unique_numbering(chain: gemmi.Chain, sources: list[str], src: Path) -> None:
    """A merged chain must not repeat a (residue number, insertion code).

    Two segments of a split protomer each arrive numbered 1..N, so preserving the
    source numbering across a merge collides. The collision is invisible in the
    file itself but silently drops residues in any consumer that indexes residues
    by number -- ringfit's ``--resnum-match resnum`` builds exactly such an index.
    """
    keys = [f"{res.seqid.num}{str(res.seqid.icode).strip()}" for res in chain]
    seen: set[str] = set()
    dups: list[str] = []
    for key in keys:
        if key in seen and key not in dups:
            dups.append(key)
        seen.add(key)
    if dups:
        shown = ",".join(dups[:5]) + ("..." if len(dups) > 5 else "")
        raise ValueError(
            f"merged chain {chain.name} of {src.name} ({'+'.join(sources)}) has "
            f"duplicate residue numbers ({shown}): the segments carry independent "
            f"numbering. Pass --renumber (1..N across the merged chain) or "
            f"--renumber-from '<chain>:<first>[-<last>],...' (auth numbering per "
            f"segment)"
        )


def extract(
    src: Path,
    groups: list[Group],
    strip_het: bool,
    renumber: bool,
    starts: dict[str, tuple[int, int]],
) -> tuple[gemmi.Structure, list[str], int, int]:
    """One structure -> (extracted structure, written chain IDs, n_res, n_atoms).

    ``n_res``/``n_atoms`` count everything written, over all output chains; the
    number of output chains is ``len(groups)``, which a merge makes smaller than
    the number of source chains.

    Raises on anything unusable (the caller records the message as data): a chain
    that is absent, one that survives HET stripping with no residues left, a
    segment whose length contradicts its --renumber-from range, or a merged chain
    whose residue numbers collide.
    """
    st = load_structure(
        src,
        guard_chains=[chain for sources, _ in groups for chain in sources],
        strip_het=strip_het,
    )
    model = st[0]
    present = [chain.name for chain in model]

    out_model = gemmi.Model(1)
    out_ids: list[str] = []
    n_res = 0
    n_atoms = 0

    for sources, out_id in groups:
        out_chain = gemmi.Chain(out_id)
        for chain_id in sources:
            if chain_id not in present:
                raise ValueError(
                    f"chain {chain_id} not in {src.name} (have: {','.join(present)})"
                )
            residues = chain_residues(model, chain_id, strip_het)
            if not residues:
                hint = (
                    " after stripping HETATM/ligands/waters (--keep-het)"
                    if strip_het
                    else ""
                )
                raise ValueError(
                    f"chain {chain_id} of {src.name} has no residues{hint}"
                )

            first = None
            if starts:
                if chain_id not in starts:
                    raise ValueError(
                        f"--renumber-from gives no start for kept chain {chain_id}"
                    )
                first, last = starts[chain_id]
                if last and len(residues) != last - first + 1:
                    raise ValueError(
                        f"chain {chain_id} of {src.name} has {len(residues)} "
                        f"residues but --renumber-from says {first}-{last} "
                        f"({last - first + 1})"
                    )

            # Written in the order the group lists its sources, never re-sorted
            # by residue number.
            for i, res in enumerate(residues):
                copy = res.clone()
                if first is not None:
                    copy.seqid = gemmi.SeqId(first + i, " ")
                elif renumber:
                    # Continuous across the whole output chain, so a merged
                    # A+C comes out 1..90 rather than 1..41 then 1..49.
                    copy.seqid = gemmi.SeqId(len(out_chain) + 1, " ")
                out_chain.add_residue(copy)
                n_atoms += len(copy)
            n_res += len(residues)

        if len(sources) > 1:
            check_unique_numbering(out_chain, sources, src)
        out_model.add_chain(out_chain)
        out_ids.append(out_id)

    if n_res == 0:
        raise ValueError(f"nothing extracted from {src.name}")

    out_st = gemmi.Structure()
    out_st.name = src.stem
    out_st.spacegroup_hm = st.spacegroup_hm
    out_st.cell = st.cell
    out_st.add_model(out_model)
    out_st.setup_entities()
    return out_st, out_ids, n_res, n_atoms


def extract_split(
    src: Path, splits: list[Split], strip_het: bool
) -> tuple[gemmi.Structure, list[str], int, int, dict[str, int], dict[str, float]]:
    """Cut source chains into segments by POSITION; one output chain per segment.

    Returns (structure, written chain IDs in output order, n_res, n_atoms,
    residues per output chain, mean B-factor per output chain).

    Positions count the residues chainsel keeps (polymer only unless --keep-het),
    in file order, as delimited in the text (microheterogeneity = two positions). Raises (the caller records it as data, no file is written) on a
    source chain that is absent or empty, or whose segments do not end exactly at
    its residue count -- the checksum. Output chains are written in ascending
    output-ID order, each renumbered 1..N; atoms are cloned unchanged.
    """
    st = load_structure_positional(src)
    model = st[0]
    present = [chain.name for chain in model]

    cut: dict[str, list[gemmi.Residue]] = {}
    for chain_id, segments in splits:
        if chain_id not in present:
            raise ValueError(
                f"chain {chain_id} not in {src.name} (have: {','.join(present)})"
            )
        # Iterated by position -- never looked up by residue number.
        residues = chain_residues(model, chain_id, strip_het)
        if not residues:
            hint = (
                " after stripping HETATM/ligands/waters (--keep-het)"
                if strip_het
                else ""
            )
            raise ValueError(f"chain {chain_id} of {src.name} has no residues{hint}")
        # remove_hydrogens() leaves a hydrogen-only residue in place with no
        # atoms: it would count as a position and in chain_n_res but write nothing.
        hollow = [i + 1 for i, res in enumerate(residues) if len(res) == 0]
        if hollow:
            raise ValueError(
                f"chain {chain_id} of {src.name} has residue(s) with no atoms left "
                f"after dropping hydrogens/alternative conformations, at "
                f"position(s) {','.join(map(str, hollow[:5]))}"
                f"{'...' if len(hollow) > 5 else ''}"
            )
        covered = segments[-1][1]
        if covered != len(residues):
            shown = "+".join(f"{a}-{b}>{o}" for a, b, o in segments)
            raise ValueError(
                f"chain {chain_id} of {src.name} has {len(residues)} residues but "
                f"--split-groups '{chain_id}:{shown}' covers {covered}"
            )
        for first, last, out_id in segments:
            cut[out_id] = residues[first - 1 : last]

    out_model = gemmi.Model(1)
    out_ids = sorted(cut)
    n_res = 0
    n_atoms = 0
    chain_n_res: dict[str, int] = {}
    chain_mean_b: dict[str, float] = {}
    for out_id in out_ids:
        out_chain = gemmi.Chain(out_id)
        b_sum = 0.0
        b_n = 0
        for i, res in enumerate(cut[out_id]):
            copy = res.clone()
            copy.seqid = gemmi.SeqId(i + 1, " ")
            if copy.label_seq is not None:
                copy.label_seq = i + 1
            # A clone keeps the SOURCE chain's subchain (label_asym_id) and
            # entity_id (verified, gemmi 0.7.5), so binder D and target A cut from
            # one source chain would share label_asym_id 'A' and collide on
            # (label_asym_id, label_seq_id) in mmCIF output. Clear them so
            # setup_entities() assigns fresh ones per output chain.
            copy.subchain = ""
            copy.entity_id = ""
            out_chain.add_residue(copy)
            for atom in copy:
                b_sum += atom.b_iso
                b_n += 1
        if b_n == 0:
            raise ValueError(f"output chain {out_id} of {src.name} has no atoms")
        out_model.add_chain(out_chain)
        chain_n_res[out_id] = len(cut[out_id])
        chain_mean_b[out_id] = b_sum / b_n
        n_res += len(cut[out_id])
        n_atoms += b_n

    out_st = gemmi.Structure()
    out_st.name = src.stem
    out_st.spacegroup_hm = st.spacegroup_hm
    out_st.cell = st.cell
    out_st.add_model(out_model)
    out_st.setup_entities()
    return out_st, out_ids, n_res, n_atoms, chain_n_res, chain_mean_b


def per_chain(values: dict, fmt: str = "{}") -> str:
    """{'A': 152, 'D': 60} -> 'A=152,D=60' (dict order = output order)."""
    return ",".join(f"{k}={fmt.format(v)}" for k, v in values.items())


def write_structure(st: gemmi.Structure, out: Path, out_format: str) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out_format == "cif":
        st.make_mmcif_document().write_file(str(out))
    else:
        st.write_pdb(str(out))


def write_result(
    result_tsv: Path,
    name: str,
    status: str,
    path: str,
    metrics: dict,
    extra_columns: list[str] | None = None,
) -> None:
    """One-row TSV. A metric that did not apply is written as NA (an empty cell),
    never as 0 -- 0 chains and 0 residues are real, alarming values.

    ``extra_columns`` (split mode's SPLIT_METRIC_COLUMNS) are appended after the
    standard ones; without them the file is exactly what it always was."""
    metric_columns = METRIC_COLUMNS + (extra_columns or [])
    result_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(result_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(RESULT_COLUMNS + (extra_columns or []))
        writer.writerow(
            [name, status, path]
            + [
                "" if metrics.get(c) is None else str(metrics[c])
                for c in metric_columns
            ]
        )


class Args(argparse.Namespace):
    name: str
    src: Path
    chains: str
    rename_to: str
    groups: str
    out_format: str
    keep_het: bool
    renumber: bool
    renumber_from: str
    split_groups: str
    out: Path
    result_tsv: Path


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Extract (and optionally merge) a subset of chains."
    )
    ap.add_argument("--name", required=True)
    ap.add_argument(
        "--src", type=Path, required=True, help="Source structure (PDB or mmCIF)."
    )
    ap.add_argument(
        "--groups",
        default="",
        help="Canonical selection '<chain>[+<chain>...]:<out_id>', comma-joined, "
        "e.g. 'A+C:A,B+D:B'. What chainsel.sh passes; supersedes "
        "--chains/--rename-to.",
    )
    ap.add_argument(
        "--chains",
        default="",
        help="Chain IDs to keep, comma-joined and ordered (standalone use).",
    )
    ap.add_argument(
        "--rename-to",
        default="",
        help="New chain IDs, same order and number as --chains. Empty keeps them.",
    )
    ap.add_argument("--out-format", choices=("pdb", "cif"), default="pdb")
    ap.add_argument(
        "--keep-het",
        action="store_true",
        help="Keep HETATM/ligands/waters. Default strips them.",
    )
    ap.add_argument(
        "--renumber",
        action="store_true",
        help="Renumber each OUTPUT chain from 1, continuous across a merge seam. "
        "Default preserves the numbering.",
    )
    ap.add_argument(
        "--renumber-from",
        default="",
        help="Per-source-chain auth numbering '<chain>:<first>[-<last>]', "
        "comma-joined, e.g. 'A:19-59,C:107-155'. Excludes --renumber.",
    )
    ap.add_argument(
        "--split-groups",
        default="",
        help="Cut source chains by 1-based POSITION: "
        "'<chain>:<first>-<last>><out_id>[+...]', comma-joined per source chain, "
        "e.g. 'A:1-60>D+61-212>A'. Segments must tile each chain exactly. Output "
        "chains in ascending ID order, each renumbered 1..N. Excludes --groups, "
        "--chains, --rename-to and --renumber-from; --renumber is a no-op.",
    )
    ap.add_argument(
        "--out", type=Path, required=True, help="Extracted structure to write."
    )
    ap.add_argument(
        "--result-tsv", type=Path, required=True, help="Per-design result TSV to write."
    )
    args = ap.parse_args(namespace=Args())

    if args.split_groups.strip():
        run_split(args)
        return

    try:
        if args.renumber and args.renumber_from.strip():
            raise ValueError("--renumber and --renumber-from are mutually exclusive")
        groups = resolve_groups(args.groups, args.chains, args.rename_to)
        starts = parse_starts(args.renumber_from)
        st, out_ids, n_res, n_atoms = extract(
            args.src, groups, not args.keep_het, args.renumber, starts
        )
        write_structure(st, args.out, args.out_format)
        metrics = {
            "n_chains": len(out_ids),
            "n_res": n_res,
            "chains": ",".join(out_ids),
            "n_atoms": n_atoms,
        }
        print(
            f"{args.name}: wrote {args.out.name} -- {len(out_ids)} chain(s) "
            f"{','.join(out_ids)}, {n_res} residues, {n_atoms} atoms"
        )
        write_result(args.result_tsv, args.name, "OK", str(args.out), metrics)
    except Exception as e:  # noqa: BLE001 - errors are recorded as data
        print(f"{args.name}: ERROR {e}", file=sys.stderr)
        # No partial file is left behind: a half-written extraction would be read
        # by the next tool as if it were complete.
        try:
            args.out.unlink(missing_ok=True)
        except OSError:
            pass
        write_result(args.result_tsv, args.name, f"error: {e}", "", {})


def run_split(args: Args) -> None:
    """Split mode: its own path so the other modes' code is untouched."""
    try:
        clashing = [
            flag
            for flag, value in (
                ("--groups", args.groups),
                ("--chains", args.chains),
                ("--rename-to", args.rename_to),
                ("--renumber-from", args.renumber_from),
            )
            if value.strip()
        ]
        if clashing:
            raise ValueError(
                f"--split-groups cannot be combined with {', '.join(clashing)}"
            )
        # --renumber is accepted and ignored: split mode always renumbers 1..N.
        splits = parse_split_groups(args.split_groups)
        st, out_ids, n_res, n_atoms, chain_n_res, chain_mean_b = extract_split(
            args.src, splits, not args.keep_het
        )
        write_structure(st, args.out, args.out_format)
        metrics = {
            "n_chains": len(out_ids),
            "n_res": n_res,
            "chains": ",".join(out_ids),
            "n_atoms": n_atoms,
            "chain_n_res": per_chain(chain_n_res),
            "chain_mean_b": per_chain(chain_mean_b, "{:.2f}"),
        }
        print(
            f"{args.name}: wrote {args.out.name} -- split into {len(out_ids)} "
            f"chain(s) {metrics['chain_n_res']}, mean B {metrics['chain_mean_b']}, "
            f"{n_atoms} atoms"
        )
        write_result(
            args.result_tsv, args.name, "OK", str(args.out), metrics,
            SPLIT_METRIC_COLUMNS,
        )
    except Exception as e:  # noqa: BLE001 - errors are recorded as data
        print(f"{args.name}: ERROR {e}", file=sys.stderr)
        try:
            args.out.unlink(missing_ok=True)
        except OSError:
            pass
        write_result(
            args.result_tsv, args.name, f"error: {e}", "", {}, SPLIT_METRIC_COLUMNS
        )


if __name__ == "__main__":
    main()
