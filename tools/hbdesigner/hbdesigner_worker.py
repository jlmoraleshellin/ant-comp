#!/usr/bin/env python3
"""
Per-design post-step for hbdesigner: turn one design dir into the one-row-per-rank
TSV that collect_hbdesigner.py folds into the table.

Three jobs, in order of how much they matter:

1. RECORD A RUN THAT PRODUCED NOTHING. ``rank_and_save`` returns early -- writing
   no PDB and no CSV -- when no network survives scoring or symmetrization, and
   the process still exits 0. Left alone that reads downstream as a design with no
   children, i.e. as nothing at all. This writes a single rank-0 row with an
   ``error:`` status instead, so the table says so.

2. TRANSLATE THE NETWORK INTO PROTEINMPNN'S COORDINATES. HBDesigner reports its
   network as PDB chain+resnum+aa (``A12S:A16T``). ProteinMPNN's
   ``--fixed-positions`` counts 1..L WITHIN each parsed chain. Those agree only
   when the PDB happens to be numbered 1..L per chain. The mapping is therefore
   computed from the output PDB's own residue order, and reported with it:

     * ``resnum_offset``   per network chain, ``index - resnum`` when constant
                           over that chain, else ``<chain>:var``
     * ``resnum_shift_max`` max |index - resnum| over the network residues; 0 means
                           the resnums in ``network`` ARE the mpnn positions
     * a hard check that the residue sitting at each mapped position really is the
       amino acid HBDesigner said it was -- a mismatch is an ``error:`` status, not
       a plausible-looking number.

3. VERIFY THE GRAFT. When chains were omitted from design they come back
   poly-glycine and the task script runs upstream's graft_seq.py. This reports
   whether that ran (``grafted``) and how much of the grafted chains now matches
   the input (``graft_identity``; 1.0 = fully restored).

Residues counted are ATOM-record amino acids, in file order, merged across repeated
chain blocks -- the same set ProteinMPNN's parser indexes. HETATM (ligands, waters,
and the virtual guide atom HBDesigner adds with --guide_res) is ignored throughout.

HOW THE PDB IS READ, AND WHY THERE ARE TWO WAYS.
Every file this worker opens is a PDB: the staged input (``ensure_pdb`` converted
any CIF at SUBMIT time, in the prosapia venv) and HBDesigner's own rank outputs.
gemmi is the preferred reader, but it is NOT present in every environment this
worker runs in -- the Modal image has it, the vib conda env ``hbdesigner_gpu``
does not -- so gemmi is imported LAZILY and Biopython (which that env does have,
and which HBDesigner itself declares as a dependency) takes over when the import
fails. Which reader ran is reported per row as ``pdb_parser`` (``gemmi`` or
``biopython``) -- read it before comparing rows produced on different backends.
If NEITHER is importable the design gets an ``error:`` row naming both failures;
the task is never crashed.

The two readers agree on the standard 20 amino acids, which is the whole of this
tool's domain, but NOT on every residue name -- see ``_residues_biopython`` for
the exact list of the 39 names only gemmi accepts, the 942 only Biopython
accepts, and the one (PCA) they read as different letters.

Usage (called by hbdesigner.sh; standalone for debugging):
    python hbdesigner_worker.py --name design_0 --design-dir out/design_0 \\
        --ref-pdb out/design_0/design_0.pdb --graft-chains B \\
        --graft-dir out/design_0/grafted --run-rc 0 --result-tsv out/design_0.tsv
"""

import argparse
import csv
import re
import sys
from pathlib import Path

# Columns of the per-design TSV. `name` and `rank` key the row; everything after
# `path` is collected as a leaf-prefixed column (collect_hbdesigner.py).
RESULT_COLUMNS = [
    "name",
    "rank",
    "status",
    "path",
    "hb_score_full",
    "hb_score_hb",
    "avg_burial",
    "saturation",
    "buried_heavy_unsats",
    "buried_unsat_hpol",
    "network",
    "network_seq",
    "n_network_res",
    "network_chains",
    "fixed_positions",
    "fixed_idx",
    "resnum_offset",
    "resnum_shift_max",
    "grafted",
    "graft_identity",
    "n_res_total",
    "pdb_parser",
    "guide_atoms_removed",
    "fold_sequence",
    "n_fold_changes",
    "fold_chain",
    "fold_seq_note",
    "centers",
    "n_centers",
]

# Upstream's stats CSV header -> our bare column name. Verified against
# `rank_and_save`, which writes exactly these columns (plus a pandas index).
CSV_FIELDS = {
    "HB_Score_full": "hb_score_full",
    "HB_Score_hb": "hb_score_hb",
    "Avg_Burial": "avg_burial",
    "saturation": "saturation",
    "buried_heavy_unsats": "buried_heavy_unsats",
    "buried_unsat_Hpol": "buried_unsat_hpol",
    "network": "network",
}

# One entry of the colon-separated `network` string: chain, residue number, and
# the one-letter code of the residue placed there (get_network_res).
NETWORK_RE = re.compile(r"^([A-Za-z])(-?\d+)([A-Z])$")

# Which reader produced a ResidueMap. Reported as the `pdb_parser` column.
PARSER_GEMMI = "gemmi"
PARSER_BIOPYTHON = "biopython"

_GEMMI_MODULE: object | None = None
_GEMMI_IMPORT_ERROR: str = ""


def load_gemmi():
    """Import gemmi on first use; return None (and remember why) if it is absent.

    Lazy on purpose: this worker runs under whatever interpreter the tool's
    activation script selects, and that environment is not guaranteed to have
    gemmi (the vib conda env ``hbdesigner_gpu`` does not). A module-level import
    would kill the whole task before a single design could be recorded.
    """
    global _GEMMI_MODULE, _GEMMI_IMPORT_ERROR
    if _GEMMI_MODULE is None and not _GEMMI_IMPORT_ERROR:
        try:
            import gemmi as _gemmi
        except Exception as e:  # noqa: BLE001 - any import failure means "no gemmi"
            _GEMMI_IMPORT_ERROR = f"{type(e).__name__}: {e}"
        else:
            _GEMMI_MODULE = _gemmi
    return _GEMMI_MODULE


def _residues_gemmi(path: Path) -> list[tuple[str, int, str]]:
    """(chain, resnum, one-letter) for every ATOM-record amino acid, via gemmi."""
    gemmi = load_gemmi()
    if gemmi is None:  # pragma: no cover - guarded by the caller
        raise ValueError(f"gemmi is not importable ({_GEMMI_IMPORT_ERROR})")
    st = gemmi.read_structure(str(path))
    if len(st) == 0:
        raise ValueError(f"{path.name}: no model")
    out: list[tuple[str, int, str]] = []
    for chain in st[0]:
        for res in chain:
            if res.het_flag != "A":
                continue
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid():
                continue
            out.append((chain.name, res.seqid.num, (info.one_letter_code or "x").upper()))
    return out


def _residues_biopython(path: Path) -> list[tuple[str, int, str]]:
    """The same list via Biopython, for environments without gemmi.

    ``PDBParser(QUIET=True)``, first model, chains and residues in file order,
    residues whose hetero flag is ``' '`` (Biopython's marker for an ATOM record;
    HETATM becomes ``H_<resname>`` and water ``W``) and which ``is_aa`` accepts.
    The one-letter code comes from the same table ``is_aa`` tests against,
    ``protein_letters_3to1_extended``.

    THE TWO READERS DO NOT AGREE ON EVERY RESIDUE NAME, and this is not papered
    over. gemmi's amino-acid table has 129 entries; Biopython's has 1032. On the
    standard 20 they are identical (checked entry by entry), which is the whole
    of this tool's domain -- rfd3 backbones and HBDesigner output, which is
    written with ``unk_to_gly=True``. Outside it:

      * 39 names gemmi calls amino acids and Biopython does not, so gemmi counts
        them and Biopython skips them. Notably UNK, SEC, PYL, the ambiguity codes
        ASX/GLX, every D-amino acid (DAL, DAR, DAS, DCY, DGL, DGN, DHI, DIL, DLE,
        DLY, DPN, DPR, DSG, DSN, DTH, DTR, DTY, DVA) and the chromophores
        (CRO, CRQ, CR2, CR8, CRF, GYS, NRQ, MDO), plus 3FG, CSH, FGA, GHP, IAS,
        MED, MLU, OMZ.
      * 942 names Biopython calls amino acids and gemmi does not, all of them
        three characters or fewer and so all able to appear in a PDB.
      * exactly one shared name with a different letter: PCA, which gemmi reads
        as E and Biopython as Q.

    Checked, not assumed. On the real crystal-derived model and six well-formed
    synthetic layouts (two chains with a TER, a repeated chain block, HETATM
    water + MSE + the guide atom, non-colliding insertion codes, two MODELs,
    SEP/PTR/MLY/unknown-name residues) the two readers build an IDENTICAL
    ResidueMap. The observed differences, all on files outside that set:

      * UNK -> gemmi counts it as 'X', Biopython skips it (20 vs 18 residues on
        the probe). Same for every D-amino acid, e.g. DAL.
      * PCA -> both count it, with DIFFERENT letters (gemmi E, Biopython Q). The
        residue count matches, so only the identity check downstream catches it.
      * one residue's atoms split across a TER (not something a real writer
        produces) -> gemmi makes two residues and ResidueMap rejects the
        ambiguous key; Biopython merges them and succeeds.

    What that costs, concretely: the count does NOT simply cancel between the
    staged input and the output, because HBDesigner REWRITES the sequence. An
    input containing UNK (or a D-amino acid) comes back from HBDesigner as
    glycine, which Biopython does count -- so on vib the output has MORE
    residues than the input and the row becomes an ``error: output has N
    residues, input has N-k`` row. That is the real consequence on vib: an error
    row, loud, not a quietly reader-dependent ``n_res_total``. Under gemmi the
    same file passes, because gemmi counts UNK on both sides. PCA is the one
    that changes a letter without changing a count, and it surfaces through the
    residue-identity check in ``network_columns``. ``pdb_parser`` is in the table
    so a cross-backend comparison can be gated on it.
    """
    try:
        from Bio.Data.PDBData import protein_letters_3to1_extended
        from Bio.PDB import PDBParser
        from Bio.PDB.Polypeptide import is_aa
    except Exception as e:  # noqa: BLE001 - recorded as data by the caller
        raise ValueError(
            f"neither gemmi nor Biopython is importable, so no PDB can be read "
            f"here: gemmi gave {_GEMMI_IMPORT_ERROR or 'no error'}, Biopython "
            f"gave {type(e).__name__}: {e}"
        ) from None

    # PERMISSIVE=False, not for strictness in general but for ONE behaviour: in
    # its default permissive mode Biopython MERGES two residues that share
    # (hetero flag, resseq, icode) into one and carries on, so a duplicated
    # residue number silently becomes a shorter chain and every position after it
    # shifts. gemmi keeps both, which ResidueMap then rejects as an ambiguous
    # key. Strict mode makes Biopython raise there too, so the duplicate is an
    # error row under BOTH readers rather than a quiet renumbering under one.
    # Verified not to reject any well-formed file in the check set.
    parser = PDBParser(QUIET=True, PERMISSIVE=False)
    try:
        structure = parser.get_structure(path.stem, str(path))
    except Exception as e:  # noqa: BLE001 - recorded as data by the caller
        raise ValueError(
            f"{path.name}: Biopython refused to parse this PDB in strict mode "
            f"({type(e).__name__}: {e}). Strict mode rejects more than just "
            f"duplicates -- malformed or truncated records and bad numbering too "
            f"-- so read the message. A 'duplicate residue' specifically means "
            f"the same chain+resnum appears twice, which makes the ProteinMPNN "
            f"position mapping ambiguous."
        ) from None
    models = list(structure)
    if not models:
        raise ValueError(f"{path.name}: no model")
    out: list[tuple[str, int, str]] = []
    for chain in models[0]:
        for res in chain:
            # ' ' is Biopython's hetero flag for an ATOM record; HETATM (ligands,
            # waters, and the HBDesigner --guide_res marker atom) is skipped.
            if res.id[0] != " ":
                continue
            resname = res.get_resname().strip().upper()
            if not is_aa(resname):
                continue
            one = protein_letters_3to1_extended.get(f"{resname:<3s}")
            if one is None:
                one = protein_letters_3to1_extended.get(resname, "X")
            out.append((chain.id, int(res.id[1]), one.upper()))
    return out


def read_residues(path: Path) -> tuple[list[tuple[str, int, str]], str]:
    """Read *path* with gemmi if it is importable, else with Biopython."""
    if load_gemmi() is not None:
        return _residues_gemmi(path), PARSER_GEMMI
    return _residues_biopython(path), PARSER_BIOPYTHON


class ResidueMap:
    """A structure's amino-acid residues, per chain, in file order.

    ``index[(chain, resnum)]`` is the 1-based position of that residue within its
    chain -- what ProteinMPNN calls position N -- and ``seq[chain]`` its one-letter
    sequence. Repeated chain blocks are merged under one chain id, HETATM is
    skipped, and insertion codes are NOT part of the key (HBDesigner folds them
    into the residue number when it writes a PDB, so they cannot be matched back).

    A ``(chain, resnum)`` SEEN TWICE IS AN ERROR, not a last-one-wins overwrite.
    It happens when insertion codes collapse onto one number, when an altloc is a
    point mutation parsed as two residues, or when a file is simply numbered
    twice. The index is the whole basis of ``fixed_positions``, and when a key is
    ambiguous it cannot say which residue the network meant -- so the design
    becomes an error row rather than a plausible position. Raised by BOTH readers
    alike, since both feed this constructor.

    ``parser`` records which reader produced it (see ``read_residues``).
    """

    def __init__(self, path: Path):
        residues, self.parser = read_residues(path)
        self.chain_order: list[str] = []
        self.seq: dict[str, str] = {}
        self.resnums: dict[str, list[int]] = {}
        self.index: dict[tuple[str, int], int] = {}
        for name, resnum, one in residues:
            if name not in self.seq:
                self.chain_order.append(name)
                self.seq[name] = ""
                self.resnums[name] = []
            if (name, resnum) in self.index:
                raise ValueError(
                    f"{path.name}: residue {name}{resnum} appears more than once "
                    f"(first as position {self.index[(name, resnum)]} of chain "
                    f"{name}, again as position {len(self.resnums[name]) + 1}). "
                    f"Chain+resnum is how HBDesigner names a network residue, so "
                    f"an ambiguous one cannot be mapped to a ProteinMPNN position. "
                    f"Usual causes: insertion codes folded onto one number, an "
                    f"altloc point mutation, or duplicated numbering."
                )
            self.seq[name] += one
            self.resnums[name].append(resnum)
            self.index[(name, resnum)] = len(self.resnums[name])

    @property
    def n_res(self) -> int:
        return sum(len(v) for v in self.resnums.values())

    def offset(self, chain: str) -> str:
        """``index - resnum`` for this chain when it is constant, else 'var'."""
        offsets = {i + 1 - num for i, num in enumerate(self.resnums[chain])}
        return str(offsets.pop()) if len(offsets) == 1 else "var"


def parse_network(network: str) -> list[tuple[str, int, str]]:
    """'A16T:A12S' -> [('A', 16, 'T'), ('A', 12, 'S')] (order as written)."""
    entries = []
    for token in network.split(":"):
        token = token.strip()
        if not token:
            continue
        m = NETWORK_RE.match(token)
        if m is None:
            raise ValueError(
                f"unparseable network entry {token!r} (expected '<chain><resnum><aa>')"
            )
        entries.append((m.group(1), int(m.group(2)), m.group(3)))
    if not entries:
        raise ValueError("empty network string")
    return entries


def network_columns(
    network: str, rmap: ResidueMap, per_chain_out: dict | None = None
) -> dict[str, object]:
    """Map the network onto ProteinMPNN positions and report the mapping's trust.

    Sorted by chain (in PDB order) then by position, so ``fixed_positions`` reads
    the way a human writes one. Sorting is safe here: ProteinMPNN preserves the
    order it is given but only ``--tied-positions`` is order-sensitive.

    ``per_chain_out``, when given, is filled with ``{chain: [(index, aa), ...]}``
    for the graft-sequence step; it is working data, not a column.
    """
    entries = parse_network(network)
    per_chain: dict[str, list[tuple[int, str, int]]] = {}
    shifts: list[int] = []
    for chain, resnum, aa in entries:
        idx = rmap.index.get((chain, resnum))
        if idx is None:
            raise ValueError(
                f"network residue {chain}{resnum} is not in the output PDB "
                f"(chains present: {','.join(rmap.chain_order) or 'none'})"
            )
        actual = rmap.seq[chain][idx - 1]
        if actual != aa:
            raise ValueError(
                f"network says {chain}{resnum} is {aa} but the output PDB has "
                f"{actual} at that position -- the residue mapping is wrong, so "
                f"fixed_positions would fix the wrong residues"
            )
        per_chain.setdefault(chain, []).append((idx, aa, resnum))
        shifts.append(abs(idx - resnum))

    chains = [c for c in rmap.chain_order if c in per_chain]
    for chain in chains:
        per_chain[chain].sort()

    if per_chain_out is not None:
        per_chain_out.clear()
        per_chain_out.update(
            {c: [(idx, aa) for idx, aa, _ in per_chain[c]] for c in chains}
        )

    flat = [(chain, idx, aa) for chain in chains for idx, aa, _ in per_chain[chain]]
    return {
        "network_seq": "".join(aa for _, _, aa in flat),
        "n_network_res": len(flat),
        "network_chains": ",".join(chains),
        # ProteinMPNN --fixed-positions syntax: '/' between chains (one group per
        # chain in --chains-to-design order), ',' between positions.
        "fixed_positions": "/".join(
            ",".join(str(idx) for idx, _, _ in per_chain[chain]) for chain in chains
        ),
        # The same positions flat, for the per-row {expr} hand-off (fix1..fixK).
        "fixed_idx": ",".join(str(idx) for _, idx, _ in flat),
        "resnum_offset": ",".join(f"{c}:{rmap.offset(c)}" for c in chains),
        "resnum_shift_max": max(shifts),
    }


# Columns written when the fold sequence cannot be built. NEVER an empty string:
# an empty sequence would read downstream as "a binder with no residues".
FOLD_SEQ_NA: dict[str, object] = {
    "fold_sequence": None,
    "n_fold_changes": None,
    "fold_chain": None,
    # Both NA together: with no fold sequence there is no binder chain to index
    # into, so "no centres" would be a claim this row cannot make.
    "centers": None,
    "n_centers": None,
}

# The centre token. Only 'HIS' is ever emitted: hbdesigner designs the RESIDUE,
# not its protonation state, so 'HIS-S' / 'HIS-P' (which assert a state) would be
# a claim this tool has not made. hisbond accepts a bare 'HIS'.
CENTER_TOKEN = "HIS"


def fold_sequence_columns(
    per_chain: dict[str, list[tuple[int, str]]],
    out_map: ResidueMap,
    ref_map: ResidueMap,
    omit_chains: list[str],
) -> tuple[dict[str, object], str]:
    """The designed binder's real sequence with the network residues written in.

    HBDesigner hands back poly-glycine outside the network, which no predictor
    can fold. The staged input -- the same file this worker already reads as the
    residue-count reference -- still carries the binder's real sequence, and the
    network is a handful of substitutions on it. So the sequence a predictor
    should actually fold is: the STAGED INPUT's binder chain, with the network's
    binder positions overwritten by the network residues. That is this column.

    Which chain is the binder: THE SINGLE NETWORK CHAIN THAT SURVIVES THE OMIT
    FILTER. With ``--omit-chains X`` and an anchor on X the network chains are
    ``X,D``, X is omitted, D survives, and D is the binder. A MONOMER RUN DOES
    GET A SEQUENCE: its one network chain is trivially the only survivor, and
    that chain is the designed molecule. The columns are NA when the filter
    leaves zero or two or more chains -- a network on two chains with neither
    omitted is the common case, and there which one is the binder is genuinely
    unknown. A guess would be a sequence for the wrong molecule.

    HOW IT IS CHECKED. Three of the checks below (length, network letters
    present, nothing else changed) only re-read the loop that built the string,
    so on their own they prove nothing. The load-bearing ones are independent of
    it:

      * the output and the staged input must have the SAME ORDERED LIST OF
        RESIDUE NUMBERS for the binder chain. The network positions are indices
        into the OUTPUT; they only address the INPUT's sequence if the two agree
        residue for residue. This is the guard that catches a renumbered or
        re-length-ed chain.
      * every emitted CENTRE index must carry an 'H' in the fold sequence. The
        centres are derived from the same ``wanted`` map the sequence was built
        from, so this is a guard against an indexing slip (0-based vs 1-based)
        rather than an independent fact -- but an off-by-one here would point
        hisbond at the wrong residue, which is exactly the class of error that
        produces a plausible number.
      * every binder position that is NOT a network position must be GLYCINE in
        the rank PDB while the fold sequence carries the staged input's residue
        there. That is a fact about a different file (HBDesigner's output, which
        is poly-glycine by construction) than the one the sequence was built
        from, so it fails if the two files are not the pair they are assumed to
        be -- a stale rank PDB, a mismatched staged input, or a binder chain
        that was never actually redesigned.

    Returns ``(columns, note)``; ``note`` is '' when the sequence was built.
    """
    candidates = [c for c in per_chain if c not in omit_chains]
    if len(candidates) != 1:
        if not omit_chains:
            why = (
                "no --omit-chains was given, so no chain is marked as the target "
                "and the binder cannot be identified"
            )
        else:
            why = (
                f"the network touches {len(candidates)} non-omitted chain(s) "
                f"{candidates or '[]'} (omitted: {','.join(omit_chains)}); "
                f"exactly one is needed"
            )
        return dict(FOLD_SEQ_NA), f"no fold sequence: {why}"
    binder = candidates[0]

    if binder not in ref_map.seq:
        return dict(FOLD_SEQ_NA), (
            f"no fold sequence: binder chain {binder} is not in the staged input"
        )
    # The network positions are indices into the OUTPUT chain; they only address
    # the INPUT's sequence if the two chains are the same residues in the same
    # order. Checked rather than assumed.
    if out_map.resnums[binder] != ref_map.resnums[binder]:
        return dict(FOLD_SEQ_NA), (
            f"no fold sequence: chain {binder} has {len(out_map.resnums[binder])} "
            f"residues numbered {_span(out_map.resnums[binder])} in the output but "
            f"{len(ref_map.resnums[binder])} numbered {_span(ref_map.resnums[binder])} "
            f"in the staged input, so the network positions do not address the "
            f"input's sequence"
        )

    ref_seq = ref_map.seq[binder]
    grafted = list(ref_seq)
    for idx, aa in per_chain[binder]:
        if not 1 <= idx <= len(grafted):
            return dict(FOLD_SEQ_NA), (
                f"no fold sequence: network position {idx} is outside chain "
                f"{binder} (length {len(grafted)})"
            )
        grafted[idx - 1] = aa
    grafted_seq = "".join(grafted)

    # Verify before writing. The first three re-read the loop above and are kept
    # only as cheap guards against a coding slip; the GLYCINE check that follows
    # them is the independent one.
    wanted = {idx: aa for idx, aa in per_chain[binder]}
    if len(grafted_seq) != len(ref_seq):
        return dict(FOLD_SEQ_NA), (
            f"no fold sequence: built {len(grafted_seq)} residues for chain "
            f"{binder} but the staged input has {len(ref_seq)}"
        )
    bad = [i for i, aa in wanted.items() if grafted_seq[i - 1] != aa]
    if bad:
        return dict(FOLD_SEQ_NA), (
            f"no fold sequence: position(s) {bad} of chain {binder} do not carry "
            f"the network residue after grafting"
        )
    stray = [
        i + 1
        for i in range(len(ref_seq))
        if (i + 1) not in wanted and grafted_seq[i] != ref_seq[i]
    ]
    if stray:
        return dict(FOLD_SEQ_NA), (
            f"no fold sequence: position(s) {stray} of chain {binder} changed "
            f"although they are not network positions"
        )

    # THE INDEPENDENT CHECK. HBDesigner writes everything it did not design as
    # glycine (`to_pdb(unk_to_gly=True)`), so in the RANK PDB every non-network
    # position of the binder must be G -- a fact about the output file, not about
    # the string built above. If some position is not glycine there, this pair of
    # files is not what the column assumes (a stale rank PDB left from an earlier
    # run, a staged input that is not the one HBDesigner read, or a chain that was
    # never cleared for design), and the sequence would be a plausible fake.
    # The network positions themselves are excluded: those are exactly the ones
    # HBDesigner did design, and they were already matched residue by residue.
    out_seq = out_map.seq[binder]
    not_gly = [
        i + 1
        for i in range(len(out_seq))
        if (i + 1) not in wanted and out_seq[i] != "G"
    ]
    if not_gly:
        shown = ", ".join(
            f"{binder}{out_map.resnums[binder][i - 1]}={out_seq[i - 1]}"
            for i in not_gly[:6]
        )
        return dict(FOLD_SEQ_NA), (
            f"no fold sequence: {len(not_gly)} non-network position(s) of chain "
            f"{binder} are not glycine in the rank PDB ({shown}"
            f"{', ...' if len(not_gly) > 6 else ''}). HBDesigner writes every "
            f"position it did not design as glycine, so this output and the "
            f"staged input are not the pair this column assumes -- suspect a "
            f"stale rank PDB from an earlier run, or an input whose chain "
            f"{binder} was not designable. TREAT THE WHOLE ROW AS SUSPECT, not "
            f"just the fold columns: if the rank PDB is from an earlier run then "
            f"network, the scores, fixed_positions and path describe that run, "
            f"not this one. Re-run the design with --force under a new "
            f"--table-label rather than reading these numbers"
        )

    # The designed histidines on the binder, in hisbond's centres format. The
    # index is 0-BASED INTO THE BINDER SEQUENCE (hisbond's contract), i.e. the
    # 1-based position within the chain minus 1 -- the same indexing as
    # fold_sequence, which is what makes the check below meaningful.
    centers = [(idx - 1, aa) for idx, aa in sorted(wanted.items()) if aa == "H"]
    for index0, _ in centers:
        if not 0 <= index0 < len(grafted_seq) or grafted_seq[index0] != "H":
            found = (
                grafted_seq[index0]
                if 0 <= index0 < len(grafted_seq)
                else "out of range"
            )
            return dict(FOLD_SEQ_NA), (
                f"no fold sequence: centre index {index0} does not carry a "
                f"histidine in the fold sequence of chain {binder} (found "
                f"{found!r}), so the centres column would point hisbond at the "
                f"wrong residue"
            )

    return (
        {
            "fold_sequence": grafted_seq,
            "n_fold_changes": sum(
                1 for i, aa in wanted.items() if ref_seq[i - 1] != aa
            ),
            "fold_chain": binder,
            # NA, never '': hisbond raises on an empty centres cell, and an
            # empty string here would be a row that looks ready for it and is
            # not. n_centers stays 0 -- that IS the measurement for a network
            # with no designed histidine.
            "centers": (
                ";".join(f"{i}:{CENTER_TOKEN}" for i, _ in centers) or None
            ),
            "n_centers": len(centers),
        },
        "",
    )


def _span(resnums: list[int]) -> str:
    return f"{resnums[0]}-{resnums[-1]}" if resnums else "(none)"


def graft_columns(
    out_map: ResidueMap, ref_map: ResidueMap, graft_chains: list[str]
) -> dict[str, object]:
    """Did the graft put the omitted chains back, and how completely?"""
    if not graft_chains:
        return {"grafted": "no", "graft_identity": None}
    matched = total = 0
    for chain in graft_chains:
        if chain not in out_map.seq:
            raise ValueError(f"grafted chain {chain} is missing from the output PDB")
        if chain not in ref_map.seq:
            raise ValueError(f"grafted chain {chain} is missing from the input PDB")
        out_seq, ref_seq = out_map.seq[chain], ref_map.seq[chain]
        if len(out_seq) != len(ref_seq):
            raise ValueError(
                f"grafted chain {chain} has {len(out_seq)} residues but the input "
                f"has {len(ref_seq)}"
            )
        matched += sum(a == b for a, b in zip(out_seq, ref_seq))
        total += len(ref_seq)
    return {
        "grafted": "yes",
        "graft_identity": round(matched / total, 4) if total else None,
    }


# --- the --guide_res guide atom -------------------------------------------
#
# With --guide_res the installed HBDesigner appends a VIRTUAL GUIDE ATOM to every
# output PDB (`add_guide_atom`: atom name 'V1', residue name 'ORI', element 'V',
# `chain_index = chain_index.max() + 1`, written because rank_and_save calls
# `to_pdb(unk_to_gly=True, no_hetatm=False)`). It is a marker of where the network
# was steered, not part of the structure, and it is in the output of EVERY rank
# while the staged input has none. Two things downstream count residues and would
# see the difference: upstream's own graft_seq.py asserts
# `p_target.n_res == p_ref.n_res`, and anything that reads the file after us gets
# a one-atom chain of element V in a poly-glycine backbone. So it is removed from
# the rank PDBs BEFORE grafting, the untouched original is kept under
# `guide_raw/`, and the number of atoms removed is reported per rank as
# `guide_atoms_removed` (0 = the file was checked and had none; NA = the strip
# step did not run or failed, which is NOT the same as zero).
#
# A line must match BOTH the residue name and the atom name, so an unrelated ORI
# ligand is not silently deleted. `add_guide_atom` writes exactly ONE such atom,
# so the expected counts are 0 and 1 -- but the count is reported as written
# rather than capped, because a value above 1 means the output does not match
# the upstream code this tool was built against and is worth seeing, not hiding.
GUIDE_ATOM_RESNAME = "ORI"
GUIDE_ATOM_NAME = "V1"
GUIDE_STRIP_TSV = "guide_strip.tsv"
GUIDE_RAW_DIR = "guide_raw"


def _is_guide_atom_line(line: str) -> bool:
    return (
        line.startswith("HETATM")
        and line[17:20].strip().upper() == GUIDE_ATOM_RESNAME
        and line[12:16].strip().upper() == GUIDE_ATOM_NAME
    )


def strip_guide_atoms(design_dir: Path, name: str) -> list[tuple[str, int]]:
    """Remove the guide-atom HETATM from each rank PDB; keep the raw file.

    Returns (pdb file name, atoms removed) per rank PDB, and records the same in
    ``<design_dir>/guide_strip.tsv`` so the summarise pass can report it. A rank
    PDB with no guide atom is left untouched and reported as 0.

    ALWAYS READS THE RANK PDB ITSELF and refreshes ``guide_raw/`` from it
    whenever something is removed, so the raw copy always belongs to THIS
    invocation. An earlier version preferred an existing ``guide_raw/`` copy to
    make the count idempotent; on a ``--force`` rerun in the same leaf that
    copied the PREVIOUS run's coordinates over this run's output. Re-running the
    strip alone on an already-stripped file now reports 0, which is the honest
    answer for that file as it stands.

    Note on bytes: the file is round-tripped through ``read_text``/``write_text``,
    which normalise line endings to ``\\n``. For a CRLF input the kept raw copy is
    therefore content-identical but not byte-identical, and so is the rewritten
    rank PDB. A file with nothing to remove is not rewritten at all.
    """
    raw_dir = design_dir / GUIDE_RAW_DIR
    tsv = design_dir / GUIDE_STRIP_TSV
    # Remove any record from an earlier invocation FIRST: if this pass dies part
    # way, the summarise step must see no record (NA) rather than the previous
    # run's counts attached to this run's files.
    tsv.unlink(missing_ok=True)
    report: list[tuple[str, int]] = []
    try:
        for pdb in sorted(design_dir.glob(f"{name}_HBDes_rank_*.pdb")):
            lines = pdb.read_text().splitlines(keepends=True)
            kept = [line for line in lines if not _is_guide_atom_line(line)]
            removed = len(lines) - len(kept)
            if removed:
                raw_dir.mkdir(parents=True, exist_ok=True)
                (raw_dir / pdb.name).write_text("".join(lines))
                pdb.write_text("".join(kept))
            report.append((pdb.name, removed))
    finally:
        # Written even on failure, so the files already handled are recorded and
        # the rest simply have no row (NA).
        with open(tsv, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["pdb", "guide_atoms_removed"])
            writer.writerows(report)
    return report


def read_guide_strip(design_dir: Path) -> dict[str, int]:
    """``{rank pdb name: atoms removed}``; empty when the strip pass did not run."""
    tsv = design_dir / GUIDE_STRIP_TSV
    if not tsv.is_file():
        return {}
    with open(tsv, newline="") as f:
        return {
            row["pdb"]: int(row["guide_atoms_removed"])
            for row in csv.DictReader(f, delimiter="\t")
            if row.get("pdb")
        }


def write_result(result_tsv: Path, rows: list[dict]) -> None:
    # DictWriter(extrasaction="ignore") would silently drop a key that is not in
    # RESULT_COLUMNS, which is exactly how a newly added measurement goes missing
    # from the table. Catch it here instead.
    for row in rows:
        unknown = set(row) - set(RESULT_COLUMNS)
        assert not unknown, (
            f"hbdesigner_worker: {sorted(unknown)} would be dropped by the TSV "
            f"writer -- add them to RESULT_COLUMNS, and to the typed column lists "
            f"in collect_hbdesigner.py"
        )
    result_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(result_tsv, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=RESULT_COLUMNS, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in RESULT_COLUMNS})


def read_stats(csv_path: Path) -> list[dict[str, str]]:
    with open(csv_path, newline="") as f:
        return list(csv.DictReader(f))


def collect_design(args: argparse.Namespace) -> list[dict]:
    design_dir = Path(args.design_dir)
    stats_csv = design_dir / f"{args.name}_HBDes_stats.csv"
    graft_chains = [c.strip() for c in args.graft_chains.split(",") if c.strip()]
    omit_chains = [c.strip() for c in args.omit_chains.split(",") if c.strip()]

    if args.run_rc != 0:
        return [
            {
                "name": args.name,
                "rank": 0,
                "status": f"error: run_hbdesigner exited {args.run_rc}",
            }
        ]

    # The graft is not optional once it was asked for: without it the omitted
    # chains are poly-glycine and every downstream interface number is computed
    # against a chain that is not there. A failed graft is therefore the design's
    # failure, recorded here, rather than a dead task with no row at all.
    if args.graft_rc != 0:
        return [
            {
                "name": args.name,
                "rank": 0,
                "status": (
                    f"error: graft failed (hbdesigner.scripts.graft_seq exited "
                    f"{args.graft_rc} for chain(s) {args.graft_chains or '?'}; "
                    f"without the graft those chains are poly-glycine, so the "
                    f"structure is wrong rather than imperfect)"
                ),
            }
        ]

    # The zero-output success: upstream returned early without writing anything.
    if not stats_csv.is_file():
        return [
            {
                "name": args.name,
                "rank": 0,
                "status": (
                    "error: no networks (hbdesigner exited 0 without writing "
                    "a stats CSV -- nothing passed scoring or symmetrization)"
                ),
            }
        ]

    # The reference parse is what every row is checked against, so a reader that
    # cannot be loaded at all (no gemmi, no Biopython) or a PDB that neither can
    # read has to be a recorded error, not a crashed task.
    has_gemmi = load_gemmi() is not None
    try:
        ref_map = ResidueMap(Path(args.ref_pdb))
    except Exception as e:  # noqa: BLE001 - errors are recorded as data
        # Do not repeat the gemmi diagnosis when the error already carries it.
        hint = "" if has_gemmi or "gemmi" in str(e) else (
            f" (gemmi is not importable here: {_GEMMI_IMPORT_ERROR}; Biopython "
            f"was the reader)"
        )
        print(f"{args.name}: ERROR reading {args.ref_pdb}: {e}", file=sys.stderr)
        return [
            {
                "name": args.name,
                "rank": 0,
                "status": f"error: cannot read the staged input PDB: {e}{hint}",
                "pdb_parser": PARSER_GEMMI if has_gemmi else PARSER_BIOPYTHON,
            }
        ]

    guide_strip = read_guide_strip(design_dir)
    rows: list[dict] = []
    for raw in read_stats(stats_csv):
        rank = int(raw["Rank"])
        row: dict = {
            "name": args.name,
            "rank": rank,
            "status": "OK",
            "pdb_parser": ref_map.parser,
        }
        pdb_name = raw.get("Output_PDB") or f"{args.name}_HBDes_rank_{rank}.pdb"
        # When grafting ran, the grafted file is the one the next tool must use.
        pdb_path = (Path(args.graft_dir) if graft_chains else design_dir) / pdb_name
        # Guide-atom removal happens on the RANK pdb, before any graft, so it is
        # keyed by the rank file's name whichever directory the row points at.
        if pdb_name in guide_strip:
            row["guide_atoms_removed"] = guide_strip[pdb_name]
        try:
            for src, dst in CSV_FIELDS.items():
                if src not in raw:
                    raise ValueError(f"stats CSV has no column {src!r}")
                row[dst] = raw[src]
            if not pdb_path.is_file():
                raise ValueError(f"missing output structure {pdb_path.name}")
            out_map = ResidueMap(pdb_path)
            if out_map.n_res != ref_map.n_res:
                raise ValueError(
                    f"output has {out_map.n_res} residues, input has "
                    f"{ref_map.n_res} -- the scaffold was not preserved "
                    f"(both read with {out_map.parser}). If the counts differ by "
                    f"a few and the parser is {PARSER_BIOPYTHON}, suspect residue "
                    f"names Biopython does not call amino acids but HBDesigner "
                    f"rewrote: 'grep -c UNK' the staged input, and check for "
                    f"D-amino acids. gemmi counts those, Biopython does not, so "
                    f"the input comes up short against an output HBDesigner "
                    f"wrote as glycine"
                )
            row["n_res_total"] = out_map.n_res
            net_per_chain: dict[str, list[tuple[int, str]]] = {}
            row.update(network_columns(raw["network"], out_map, net_per_chain))
            row.update(graft_columns(out_map, ref_map, graft_chains))
            # Non-fatal by construction: it returns NA + a note, never raises,
            # so a row keeps its science columns when the binder cannot be named.
            graft_seq, note = fold_sequence_columns(
                net_per_chain, out_map, ref_map, omit_chains
            )
            row.update(graft_seq)
            row["fold_seq_note"] = note or None
            row["path"] = str(pdb_path)
        except Exception as e:  # noqa: BLE001 - errors are recorded as data
            print(f"{args.name} rank {rank}: ERROR {e}", file=sys.stderr)
            row["status"] = f"error: {e}"
            row["path"] = ""
        rows.append(row)

    if not rows:
        return [
            {
                "name": args.name,
                "rank": 0,
                "status": "error: no networks (stats CSV is empty)",
            }
        ]
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description="Summarise one hbdesigner design dir.")
    ap.add_argument(
        "--mode",
        choices=("summarise", "strip-guide"),
        default="summarise",
        help="'summarise' (default) writes the per-design TSV. 'strip-guide' "
        "removes the --guide_res guide-atom HETATM from the rank PDBs, keeping "
        "the originals under guide_raw/; it must run BEFORE graft_seq.py, whose "
        "n_res assertion the extra entry can trip.",
    )
    ap.add_argument("--name", required=True, help="Design (parent row) name.")
    ap.add_argument("--design-dir", required=True, help="This design's output dir.")
    ap.add_argument(
        "--ref-pdb",
        help="The PDB handed to HBDesigner: the graft reference and the residue-"
        "count reference. Required for --mode summarise.",
    )
    ap.add_argument(
        "--graft-chains",
        default="",
        help="Comma-joined chains that were grafted back ('' = none grafted).",
    )
    ap.add_argument("--graft-dir", default="", help="Dir holding the grafted PDBs.")
    ap.add_argument(
        "--omit-chains",
        default="",
        help="Comma-joined chains that were omitted from design ('' = none). Used "
        "to identify the designed binder chain for the fold sequence: it is the "
        "single network chain that is NOT omitted.",
    )
    ap.add_argument(
        "--run-rc", type=int, default=0, help="Exit code of run_hbdesigner."
    )
    ap.add_argument(
        "--graft-rc",
        type=int,
        default=0,
        help="Exit code of the graft_seq.py loop (0 when no graft was asked for).",
    )
    ap.add_argument("--result-tsv", help="Per-design TSV to write (--mode summarise).")
    args = ap.parse_args()

    if args.mode == "strip-guide":
        report = strip_guide_atoms(Path(args.design_dir), args.name)
        total = sum(n for _, n in report)
        print(
            f"{args.name}: guide-atom strip checked {len(report)} rank PDB(s), "
            f"removed {total} HETATM {GUIDE_ATOM_RESNAME} line(s)"
        )
        return

    for required in ("ref_pdb", "result_tsv"):
        if not getattr(args, required):
            ap.error(f"--{required.replace('_', '-')} is required for --mode summarise")

    rows = collect_design(args)
    write_result(Path(args.result_tsv), rows)
    ok = sum(1 for r in rows if r["status"] == "OK")
    print(f"{args.name}: wrote {len(rows)} row(s), {ok} OK -> {args.result_tsv}")


if __name__ == "__main__":
    main()
