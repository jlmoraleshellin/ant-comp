#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) that extracts a named subset of chains from each
design's structure into a new structure file, optionally MERGING several source
chains into one output chain.

Every structure comparison we have is whole-file. ``usalign --mm 1`` on a binder
complex is dominated by the large, fixed target and reads ~1.0 whatever the binder
did; ``--ter 2`` isolates only the FIRST chain, and the binder is the LAST. So a
binder-vs-binder comparison is simply not expressible today. ``chainsel`` makes it
expressible: it writes ``<name>.pdb`` holding only the chains you name, and collects
its path as ``chainsel_path``, which ``usalign --col-a/--col-b`` (with ``--mm 0``)
then consumes like any other structure column.

That splits the two questions a binder campaign must keep apart:

* *did it fold?* -- binder chain of the prediction vs. binder chain of the design,
  each extracted alone, compared with ``--mm 0``;
* *did it stay?* -- the whole complex, compared in the target's frame (that is
  ``usalign --mm 1`` or ``ringfit``, not this tool).

The second reason the tool exists is the SPLIT PROTOMER. A target trimmed of a
membrane belt (SlyB residues 19-59 + 107-155, with a real 13.4 A gap between the
two segments) is handed to a structure predictor as two entities, so the predictor
returns one protomer as TWO chains. ``ringfit`` takes exactly two
``--target-chains`` -- one per protomer -- and cannot be told "t1 is the union of
A and C", so ``bsa_t1``/``bsa_t2``/``bridge_ratio`` are not computable at all on
such a file. ``--merge-groups 'A+C:A,B+D:B'`` writes the same coordinates back as
two chains, and ringfit works again.

The third reason is the inverse: the FUSED subunit. RFdiffusion's symmetric mode
writes one chain per subunit that holds the diffused binder AND the fixed target,
with residue numbers that repeat inside the chain (binder 1..60, then target
1..152). ``epitope``, ``atomium`` and ``usalign`` need binder and target as separate
chains with unique numbering. ``--split-groups 'A:1-60>D+61-212>A,...'`` cuts each
source chain by 1-based POSITION (numbers cannot address a residue when they
repeat), checks that the segments tile the chain exactly (a checksum: any mismatch
is an error for that design, never a trimmed cut), writes the output chains in
ascending ID order and ALWAYS renumbers each 1..N. It also collects, per output
chain, the residue count and mean B-factor (``chain_n_res``/``chain_mean_b``) -- the
trust columns that say which segment is binder (RFdiffusion: B 0) and which target
(B 1). Those two columns are NA in every other mode.

Merging is ONLY ever done when asked for explicitly. ``--rename-to`` still REFUSES
to map two kept chains onto one ID, because there it would be an accident; a merge
group says so in the flag itself.

``action: update`` -- an extracted sub-structure is a *property* of a design that
already exists (like a relaxed structure or a metric), not a new entity, so it
annotates the same table in place and mints no child. Run it twice with different
``-l/--dir-label`` to extract two different chain sets onto one table.

``default_input_column`` is the literal string ``"not applicable"`` -- the same
sentinel ``usalign`` uses. There is no honest default here: the tool is equally at
home on ``rfdiffusion3_path``, ``boltz_path`` or ``alphafold3_path``, and picking one
would make the other two fail *silently* (a column that is absent from the frame is a
no-op for ``filter_ready``, so every row would look ready and then read a missing
cell). The sentinel matches no column, and the builder raises a named error at submit
time instead -- so **-i/--input-column is effectively required**.

The chain lettering of a designed complex and of its prediction need not agree;
``--rename-to`` normalises the extracted chain to a fixed letter so the two files are
directly comparable.

NUMBERING (the whole reason a merge can go wrong quietly) -- exactly one of three
regimes, and the output always carries one of them:

* neither flag (default): every residue keeps the source residue number and
  insertion code. Two merged segments each numbered 1..N therefore COLLIDE, and the
  worker refuses that design with an error rather than writing a file whose residue
  numbers repeat (a ``{resnum: residue}`` index downstream would silently drop half
  of it);
* ``--renumber``: each OUTPUT chain is renumbered 1..N in written order, insertion
  codes dropped. A merged ``A+C:A`` comes out 1..90 continuous ACROSS the seam;
* ``--renumber-from 'A:19-59,C:107-155'``: each SOURCE chain's segment is renumbered
  consecutively from its own first number, insertion codes dropped -- so the output
  carries the target's real auth numbering and ``ringfit --resnum-match resnum`` can
  pair it against the reference assembly. The optional ``-<last>`` is checked against
  the residues actually written, so a segment of the wrong length is an error.

Usage:
    sapia run chainsel outputs/20260928_100140_7ojg_binder \
        --table table1 \
        --input-column boltz_path \
        --chains E --rename-to A \
        --dir-label binder_pred

    sapia run chainsel outputs/20260928_100140_7ojg_binder \
        --table table1 \
        --input-column boltz_path \
        --merge-groups 'A+C:A,B+D:B,E:E' \
        --renumber-from 'A:19-59,C:107-155,B:19-59,D:107-155' \
        --dir-label protomers

    sapia run chainsel outputs/<run_dir> \
        --table table0 \
        --input-column rfdiffusion_path \
        --split-groups 'A:1-60>D+61-212>A,B:1-60>E+61-212>B,C:1-60>F+61-212>C' \
        --dir-label split
"""

from argparse import ArgumentParser
from pathlib import Path
from typing import cast

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.utils import resolve_template

OUT_FORMATS = ("pdb", "cif")
# The sentinel default_input_column (see the module docstring): matches no column,
# so the builder can refuse the run with a real message instead of silently
# submitting rows whose input cell does not exist.
NO_DEFAULT_COLUMN = "not applicable"

# One output chain: the source chain IDs that feed it, in write order, and the ID
# it is written under. A plain --chains selection is groups of one.
Group = tuple[list[str], str]
# One split source chain: its ID and its segments (first, last, out_id) in 1-based
# inclusive POSITIONS of the chain's residue order.
Segment = tuple[int, int, str]
Split = tuple[str, list[Segment]]


class ChainselArgs(CommonArgs):
    chains: str
    rename_to: str
    merge_groups: str
    out_format: str
    keep_het: bool
    renumber: bool
    renumber_from: str
    split_groups: str


def add_run_chainsel_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--chains",
        type=str,
        default="",
        help="Chain IDs to KEEP, comma-joined and ordered, e.g. 'E' or 'A,B'. The "
        "output file writes them in the order given, one output chain each. A "
        "chain named here but ABSENT from a design is an error for that design "
        "(status 'error: chain X not in ...'), never a silently shorter file. "
        "Required unless --merge-groups is given; the two are alternatives.",
    )
    parser.add_argument(
        "--rename-to",
        type=str,
        default="",
        help="New chain IDs for the kept chains, comma-joined and in the same order "
        "and number as --chains. Empty (default) keeps the original IDs. Use it to "
        "normalise a binder chain to one fixed letter across tables whose chain "
        "lettering differs (rfdiffusion3's binder may be 'B' while the predictor "
        "calls it 'E'), so the two extracted files are directly comparable. Two "
        "kept chains may NOT be renamed onto the same ID -- that would merge them "
        "by accident; say so explicitly with --merge-groups instead. Cannot be "
        "combined with --merge-groups (a group states its output ID after ':').",
    )
    parser.add_argument(
        "--merge-groups",
        type=str,
        default="",
        help="Groups of source chains to MERGE into one output chain: "
        "'<chain>[+<chain>...]:<out_id>', groups comma-joined, e.g. "
        "'A+C:A,B+D:B'. Replaces --chains/--rename-to (it says both which chains "
        "to keep and what to call them), and is the ONLY way to get two source "
        "chains into one output chain. Residues are written in the order the "
        "group lists them -- A's residues then C's -- and are NEVER re-sorted by "
        "residue number, because two segments of a split protomer carry "
        "independent numbering from the predictor and sorting would interleave "
        "them. Use it when a target protomer was predicted as two chains (a "
        "membrane belt trimmed out of the middle) and a downstream tool such as "
        "ringfit needs exactly one chain per protomer. Merged numbering is "
        "governed by --renumber / --renumber-from; colliding numbers in a merged "
        "chain are an error for that design.",
    )
    parser.add_argument(
        "--out-format",
        type=str,
        choices=OUT_FORMATS,
        default="pdb",
        help="Format of the extracted structure: 'pdb' (default, what USalign and "
        "most Rosetta tooling want) or 'cif'. PDB format has one-character chain "
        "IDs, so a longer (or renamed-to-longer) ID is refused at submit time.",
    )
    parser.add_argument(
        "--keep-het",
        action="store_true",
        help="Keep HETATM records -- ligands, ions, waters, anything not a "
        "tabulated amino-acid or nucleic-acid residue. The default STRIPS them, "
        "because a waters-and-ions tail inflates the atom counts and makes "
        "structure comparison depend on what the predictor happened to model.",
    )
    parser.add_argument(
        "--renumber",
        action="store_true",
        help="Renumber each OUTPUT chain from 1 (and drop insertion codes). With "
        "--merge-groups that means the merged chain is numbered continuously "
        "ACROSS the seam: 'A+C:A' with a 41- and a 49-residue segment comes out "
        "1..90, not 1..41 then 1..49. Default off, which preserves the source "
        "numbering. Turn it on when a downstream step pairs residues BY NUMBER "
        "across two extracted files whose numbering differs; leave it off "
        "whenever the original numbering is what you will look up (hotspots, "
        "mutations). Mutually exclusive with --renumber-from.",
    )
    parser.add_argument(
        "--renumber-from",
        type=str,
        default="",
        help="Restore auth numbering per SOURCE chain: "
        "'<chain>:<first>[-<last>]', comma-joined, e.g. "
        "'A:19-59,C:107-155'. Each source chain's residues are renumbered "
        "consecutively from <first> in written order (insertion codes dropped), "
        "so a merged protomer carries the target's real residue numbers and "
        "'ringfit --resnum-match resnum' can pair it against the reference "
        "assembly. The optional '-<last>' is CHECKED against the residues "
        "actually written: a segment of a different length is an error for that "
        "design, not a shifted numbering. Every kept source chain must be listed "
        "exactly once (a half-renumbered output is refused at submit time). "
        "Mutually exclusive with --renumber. Note this assumes each segment is "
        "internally gapless -- it writes <first>, <first>+1, ... and does not "
        "reproduce numbering gaps inside a segment.",
    )
    parser.add_argument(
        "--split-groups",
        type=str,
        default="",
        help="CUT source chains into several output chains, the inverse of "
        "--merge-groups: '<chain>:<first>-<last>><out_id>[+<first>-<last>><out_id>"
        "...]', source chains comma-joined, e.g. "
        "'A:1-60>D+61-212>A,B:1-60>E+61-212>B,C:1-60>F+61-212>C'. <first>/<last> "
        "are 1-based POSITIONS in the chain's residue order in the file (counted "
        "over the residues chainsel keeps: polymer only unless --keep-het; "
        "alternative residues at one site count as two), NOT "
        "residue numbers -- RFdiffusion's symmetric mode numbers binder and target "
        "from 1 inside one chain, so numbers repeat and cannot address anything. "
        "CHECKSUM: each source chain's segments must tile it exactly -- from "
        "position 1, in order, no gap, no overlap (refused at submit time) and "
        "ending at the chain's real residue count (else that design is an error, "
        "no file written). A source chain not named is DROPPED; a named one that is "
        "absent is an error. Output chain IDs must be unique. Output chains are "
        "written in ascending ID order (A,B,C,D,E,F) whatever the listing order. "
        "NUMBERING: split mode ALWAYS renumbers each output chain 1..N (insertion "
        "codes dropped), because the source numbering repeats; --renumber is "
        "accepted as a redundant no-op, --renumber-from is an error. Replaces "
        "--chains/--rename-to/--merge-groups (any combination is an error). Accepts "
        "{expr} placeholders resolved per design, e.g. 'A:1-{binder_len}>D+...'. "
        "Collects chain_n_res and chain_mean_b per output chain -- read them to "
        "confirm which segment is binder and which target.",
    )


def _split(value: str) -> list[str]:
    """Comma-joined list -> stripped, non-empty tokens."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


def _parse_merge_groups(value: str) -> list[Group]:
    """'A+C:A,B+D:B' -> [(['A', 'C'], 'A'), (['B', 'D'], 'B')]."""
    groups: list[Group] = []
    for spec in _split(value):
        if spec.count(":") != 1:
            raise ValueError(
                f"--merge-groups entry {spec!r} is malformed: write "
                f"'<chain>[+<chain>...]:<out_id>', e.g. 'A+C:A', with the groups "
                f"comma-joined."
            )
        sources_str, out_id = spec.split(":")
        sources = [tok.strip() for tok in sources_str.split("+") if tok.strip()]
        out_id = out_id.strip()
        if not sources:
            raise ValueError(
                f"--merge-groups entry {spec!r} names no source chain before ':'."
            )
        if not out_id:
            raise ValueError(
                f"--merge-groups entry {spec!r} names no output chain ID after ':'."
            )
        if len(set(sources)) != len(sources):
            raise ValueError(
                f"--merge-groups entry {spec!r} names a source chain twice: {sources}."
            )
        groups.append((sources, out_id))
    if not groups:
        raise ValueError(
            "--merge-groups is empty: give at least one group, e.g. 'A+C:A'."
        )
    return groups


def resolve_groups(args: ChainselArgs) -> list[Group]:
    """The chain selection as output groups, whichever flag spelled it.

    Guards the --chains / --rename-to / --merge-groups interaction: --chains and
    --merge-groups are alternatives (exactly one), and --rename-to belongs to
    --chains only. Merging happens only through --merge-groups -- the --rename-to
    collision check below stays a hard error on purpose.
    """
    merge = args.merge_groups.strip()
    chains = _split(args.chains)
    rename_to = _split(args.rename_to)

    if merge and chains:
        raise ValueError(
            "--chains and --merge-groups are alternatives: --merge-groups already "
            "says which chains to keep (before ':') and what to call them (after "
            "':'). Drop --chains, or drop --merge-groups."
        )
    if merge and rename_to:
        raise ValueError(
            "--rename-to cannot be combined with --merge-groups: a group's output "
            "chain ID is the part after ':' (e.g. 'A+C:A'). Drop --rename-to."
        )
    if not merge and not chains:
        raise ValueError(
            "chainsel needs a chain selection: pass --chains (one output chain per "
            "ID, e.g. --chains E) or --merge-groups (several source chains joined "
            "into one output chain, e.g. --merge-groups 'A+C:A,B+D:B') or "
            "--split-groups (one source chain cut by position into several, e.g. "
            "--split-groups 'A:1-60>D+61-212>A')."
        )

    if merge:
        groups = _parse_merge_groups(merge)
    else:
        if len(set(chains)) != len(chains):
            raise ValueError(f"--chains lists a chain twice: {chains}.")
        if rename_to:
            if len(rename_to) != len(chains):
                raise ValueError(
                    f"--rename-to must have the same number of IDs as --chains "
                    f"({len(chains)}), got {len(rename_to)}: {rename_to}."
                )
            if len(set(rename_to)) != len(rename_to):
                raise ValueError(
                    f"--rename-to lists an ID twice: {rename_to}. Two kept chains "
                    f"cannot share an ID -- they would merge into one. If merging "
                    f"is what you want, say so: --merge-groups 'A+C:A'."
                )
        out_ids = rename_to or chains
        groups = [([chain], out_id) for chain, out_id in zip(chains, out_ids)]

    sources = [chain for group_sources, _ in groups for chain in group_sources]
    if len(set(sources)) != len(sources):
        raise ValueError(
            f"--merge-groups writes a source chain more than once: {sources}. "
            f"Each source chain belongs to exactly one group."
        )
    out_ids = [out_id for _, out_id in groups]
    if len(set(out_ids)) != len(out_ids):
        raise ValueError(
            f"--merge-groups gives two groups the same output chain ID: {out_ids}. "
            f"That would merge them implicitly -- write them as one group instead "
            f"(e.g. 'A+C:A')."
        )
    return groups


def _parse_renumber_from(value: str, kept: list[str]) -> dict[str, tuple[int, int]]:
    """'A:19-59,C:107' -> {'A': (19, 59), 'C': (107, 0)} (0 = no end given).

    Every kept source chain must appear exactly once: a partially renumbered
    output would mix two numbering regimes in one file with nothing to mark it.
    """
    starts: dict[str, tuple[int, int]] = {}
    for spec in _split(value):
        if spec.count(":") != 1:
            raise ValueError(
                f"--renumber-from entry {spec!r} is malformed: write "
                f"'<chain>:<first>' or '<chain>:<first>-<last>', e.g. 'A:19-59'."
            )
        chain, rng = (tok.strip() for tok in spec.split(":"))
        if not chain:
            raise ValueError(f"--renumber-from entry {spec!r} names no chain.")
        if chain in starts:
            raise ValueError(f"--renumber-from names chain {chain} twice.")
        first_str, _, last_str = rng.partition("-")
        try:
            first = int(first_str)
            last = int(last_str) if last_str else 0
        except ValueError:
            raise ValueError(
                f"--renumber-from entry {spec!r} has a non-integer residue number."
            ) from None
        if last_str and last < first:
            raise ValueError(
                f"--renumber-from entry {spec!r} ends ({last}) before it starts "
                f"({first})."
            )
        starts[chain] = (first, last)

    unknown = sorted(set(starts) - set(kept))
    if unknown:
        raise ValueError(
            f"--renumber-from names chain(s) {unknown} that this run does not keep. "
            f"Kept source chains: {kept}."
        )
    unlisted = [chain for chain in kept if chain not in starts]
    if unlisted:
        raise ValueError(
            f"--renumber-from must give a start for EVERY kept source chain; "
            f"{unlisted} are missing. A half-renumbered file mixes two numbering "
            f"regimes with nothing in the output to mark which is which."
        )
    return starts


def _starts_spec(starts: dict[str, tuple[int, int]]) -> str:
    """Normalised '<chain>:<first>[-<last>]' spec for the manifest/worker."""
    return ",".join(
        f"{chain}:{first}" if not last else f"{chain}:{first}-{last}"
        for chain, (first, last) in starts.items()
    )


def _groups_spec(groups: list[Group]) -> str:
    """Canonical 'A+C:A,B:B' spec: what the worker is actually handed, so a plain
    --chains selection and a merge take the same code path downstream."""
    return ",".join(f"{'+'.join(sources)}:{out_id}" for sources, out_id in groups)


def _parse_split_groups(value: str) -> list[Split]:
    """'A:1-60>D+61-212>A,B:1-60>E' -> [('A', [(1, 60, 'D'), (61, 212, 'A')]), ...].

    The submit-side copy of chainsel_worker.parse_split_groups (the worker must run
    without prosapia, so the two cannot share a module); tests/test_chainsel.py
    checks that they accept and reject the same specs. Checks everything that does
    not need the structure: the grammar, a source chain named twice, a duplicated
    output ID, and the tiling (from position 1, in order, no gap, no overlap).
    Whether a chain really has <last> residues is checked per design by the worker.
    """
    splits: list[Split] = []
    seen_chains: set[str] = set()
    out_ids: list[str] = []
    entries = _split(value)
    if not entries:
        raise ValueError(
            "--split-groups is empty: give at least one source chain, e.g. "
            "'A:1-60>D+61-212>A'."
        )
    for entry in entries:
        if entry.count(":") != 1:
            raise ValueError(
                f"--split-groups entry {entry!r} is malformed: write "
                f"'<chain>:<first>-<last>><out_id>[+<first>-<last>><out_id>...]', "
                f"e.g. 'A:1-60>D+61-212>A', with source chains comma-joined."
            )
        chain, segs_str = (tok.strip() for tok in entry.split(":"))
        if not chain:
            raise ValueError(f"--split-groups entry {entry!r} names no source chain.")
        if chain in seen_chains:
            raise ValueError(
                f"--split-groups names source chain {chain} twice: give all of its "
                f"segments in one entry, joined with '+'."
            )
        seen_chains.add(chain)
        segments: list[Segment] = []
        for tok in segs_str.split("+"):
            tok = tok.strip()
            if tok.count(">") != 1:
                raise ValueError(
                    f"--split-groups segment {tok!r} of {entry!r} is malformed: "
                    f"write '<first>-<last>><out_id>', e.g. '1-60>D'."
                )
            rng, out_id = (part.strip() for part in tok.split(">"))
            if not out_id:
                raise ValueError(
                    f"--split-groups segment {tok!r} of {entry!r} names no output "
                    f"chain ID after '>'."
                )
            first_str, sep, last_str = rng.partition("-")
            if not sep:
                raise ValueError(
                    f"--split-groups segment {tok!r} of {entry!r} needs a position "
                    f"range '<first>-<last>' before '>'."
                )
            try:
                first, last = int(first_str), int(last_str)
            except ValueError:
                raise ValueError(
                    f"--split-groups segment {tok!r} of {entry!r} has a non-integer "
                    f"position (positions are 1-based integers; wrap a table column "
                    f"in braces, e.g. '1-{{binder_len}}>D')."
                ) from None
            if first < 1:
                raise ValueError(
                    f"--split-groups segment {tok!r} of {entry!r}: positions are "
                    f"1-based."
                )
            if last < first:
                raise ValueError(
                    f"--split-groups segment {tok!r} of {entry!r} ends ({last}) "
                    f"before it starts ({first})."
                )
            segments.append((first, last, out_id))
            out_ids.append(out_id)
        expected = 1
        for first, last, out_id in segments:
            if first < expected:
                raise ValueError(
                    f"--split-groups chain {chain}: segment {first}-{last}>{out_id} "
                    f"overlaps the previous one (next free position is {expected}). "
                    f"Segments must tile the chain in order, each residue once."
                )
            if first > expected:
                what = "the chain start" if expected == 1 else "the previous segment"
                raise ValueError(
                    f"--split-groups chain {chain}: gap -- segment "
                    f"{first}-{last}>{out_id} starts at {first} but {what} leaves "
                    f"position {expected} uncovered. Segments must tile the chain "
                    f"from position 1, in order; to drop residues, there is no "
                    f"split syntax -- that is deliberate."
                )
            expected = last + 1
        splits.append((chain, segments))
    dups = sorted({o for o in out_ids if out_ids.count(o) > 1})
    if dups:
        raise ValueError(
            f"--split-groups writes output chain ID(s) {dups} more than once. Each "
            f"segment is its own output chain and needs its own ID (split mode never "
            f"merges; use --merge-groups for that)."
        )
    return splits


def _split_guards(args: ChainselArgs) -> None:
    """--split-groups replaces every other selection/numbering flag except
    --renumber, which it makes redundant (split mode always renumbers)."""
    clashing = [
        flag
        for flag, value in (
            ("--chains", args.chains),
            ("--rename-to", args.rename_to),
            ("--merge-groups", args.merge_groups),
            ("--renumber-from", args.renumber_from),
        )
        if value.strip()
    ]
    if clashing:
        raise ValueError(
            f"--split-groups cannot be combined with {', '.join(clashing)}: it says "
            f"which chains to keep (before ':'), how to cut them and what to call "
            f"each piece (after '>'), and it always renumbers each output chain "
            f"1..N because the source numbering repeats. Drop "
            f"{', '.join(clashing)}."
        )


def _canonical_split(value: str, out_format: str) -> str:
    """Parse and check one (already {expr}-resolved) split spec; return the
    canonical text the worker is handed."""
    splits = _parse_split_groups(value)
    if out_format == "pdb":
        too_long = [o for _, segs in splits for _, _, o in segs if len(o) != 1]
        if too_long:
            raise ValueError(
                f"--out-format pdb needs one-character chain IDs, but the "
                f"--split-groups output IDs {too_long} are longer. Shorten them "
                f"(the part after '>'), or use --out-format cif."
            )
    return ",".join(
        f"{chain}:" + "+".join(f"{a}-{b}>{o}" for a, b, o in segments)
        for chain, segments in splits
    )


def _build_split_manifest(ctx: ManifestCtx[ChainselArgs]) -> list[tuple[str, ...]]:
    """Split-mode rows: the canonical split spec in field 5, the group and
    --renumber-from fields empty (the worker renumbers 1..N in split mode whatever
    the numbering field says)."""
    _split_guards(ctx.args)
    raw = ctx.args.split_groups.strip()
    # Without placeholders the spec is the same for every design: check it once,
    # here, before the input column, so a bad spec refuses the run even on an
    # empty ready set.
    static = None if "{" in raw else _canonical_split(raw, ctx.args.out_format)
    column = _check_input_column(ctx)

    manifest_rows: list[tuple[str, ...]] = []
    for name in ctx.ready.index:
        name = cast(str, name)
        src = Path(str(ctx.ready.at[name, column]))
        if not src.exists():
            print(f"{name}: MISSING {src} (skipping)")
            continue
        if static is not None:
            spec = static
        else:
            try:
                spec = _canonical_split(
                    resolve_template(raw, ctx.lookup, name), ctx.args.out_format
                )
            except ValueError as e:
                raise ValueError(f"design {name!r}: {e}") from None
        manifest_rows.append(
            (
                name,
                str(volume_path(src)),
                "",  # no group spec in split mode
                "",  # no --renumber-from in split mode
                spec,
                ctx.args.out_format,
                "keep" if ctx.args.keep_het else "strip",
                "renumber" if ctx.args.renumber else "preserve",
            )
        )
    return manifest_rows


def _check_input_column(ctx: ManifestCtx[ChainselArgs]) -> str:
    # No honest default input column (see module docstring): refuse loudly rather
    # than submit rows whose input cell does not exist.
    column = ctx.args.input_column
    if column == NO_DEFAULT_COLUMN or column not in ctx.df.columns:
        available = ", ".join(
            str(c) for c in ctx.df.columns if str(c).endswith("_path")
        )
        raise ValueError(
            f"chainsel has no default input column: pass -i/--input-column with the "
            f"structure column to extract from (got {column!r}, which table "
            f"'{ctx.args.table}' does not have). Structure columns available: "
            f"{available or '(none)'}."
        )
    return column


def build_chainsel_manifest(ctx: ManifestCtx[ChainselArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool

    if ctx.args.split_groups.strip():
        return _build_split_manifest(ctx)

    groups = resolve_groups(ctx.args)
    kept = [chain for sources, _ in groups for chain in sources]
    out_ids = [out_id for _, out_id in groups]

    if ctx.args.out_format == "pdb":
        # PDB format has a single chain-ID column; a longer id would be truncated
        # (silently merging chains) by the writer.
        too_long = [c for c in out_ids if len(c) != 1]
        if too_long:
            raise ValueError(
                f"--out-format pdb needs one-character chain IDs, but the output "
                f"IDs {too_long} are longer. Rename them with --rename-to (or the "
                f"':<out_id>' of a merge group), or use --out-format cif."
            )

    renumber_from = ctx.args.renumber_from.strip()
    if renumber_from and ctx.args.renumber:
        raise ValueError(
            "--renumber and --renumber-from both set the output numbering and "
            "cannot be combined: --renumber numbers each OUTPUT chain 1..N "
            "(continuous across a merge seam), --renumber-from numbers each SOURCE "
            "chain from the number you give it (the target's auth numbering). Pick "
            "one."
        )
    starts_spec = (
        _starts_spec(_parse_renumber_from(renumber_from, kept)) if renumber_from else ""
    )

    column = _check_input_column(ctx)

    ready = ctx.ready

    groups_spec = _groups_spec(groups)

    manifest_rows: list[tuple[str, ...]] = []
    for name in ready.index:
        name = cast(str, name)
        # No CIF->PDB staging: the worker reads either format with gemmi.
        src = Path(str(ready.at[name, column]))
        if not src.exists():
            print(f"{name}: MISSING {src} (skipping)")
            continue
        manifest_rows.append(
            (
                name,
                str(volume_path(src)),
                # Canonical selection: '<src>[+<src>...]:<out_id>' per output
                # chain, so --chains and --merge-groups are one thing downstream.
                groups_spec,
                # May be empty, so never last.
                starts_spec,
                # The --split-groups spec: always empty here (split mode builds
                # its rows in _build_split_manifest).
                "",
                # Kept last: always non-empty, so the manifest line never ends on
                # an empty field.
                ctx.args.out_format,
                "keep" if ctx.args.keep_het else "strip",
                "renumber" if ctx.args.renumber else "preserve",
            )
        )

    return manifest_rows
