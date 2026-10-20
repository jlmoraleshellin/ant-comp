#!/usr/bin/env python3
"""
Per-array-task step of the af3conf tool: the binder-specific confidence of ONE
AlphaFold3 prediction, computed from the files AF3 wrote beside the predicted structure
(nothing is re-predicted). The AF3 counterpart of boltzconf; how to run it and how to
read its columns: ``.claude/skills/af3conf/SKILL.md``.

A design is an ERROR ROW, never a default or a trimmed fit, unless:

  * there is exactly ONE binder chain and at least one disjoint target chain
    (``chain_ptm`` and ``chain_pair_iptm`` are per chain, so two binder chains could
    only be an average that reads like one chain's value);
  * the structure is AF3's own ``<prefix>_model.cif``, with both sibling JSONs;
  * there is one token per residue: a multi-atom ligand or modified residue is tokenised
    per atom and breaks the count, and a one-atom ligand in a NAMED chain has no
    label_seq_id (``check_tokens``);
  * the named chains are protein (``reject_nucleic_acid``, by component name): the ipSAE
    d0 floor here is the reference's PROTEIN 1.0 (2.0 for nucleic acids), so scoring one
    would inflate every ipsae_* column while the trust columns stayed clean;
  * every number is finite.

SIBLINGS. The input column names ``<dir>/<prefix>_model.cif``; AF3 writes beside it
``<prefix>_confidences.json`` (pae, atom_plddts, atom_chain_ids, token_chain_ids,
token_res_ids) and ``<prefix>_summary_confidences.json`` (chain_ptm, chain_pair_iptm and
chain_pair_pae_min as POSITIONAL chain arrays, ranking_score, ...). ``<prefix>`` is
taken from the FILE, never the design name: they differ when AF3 sanitised the name, and
a sibling picked by design name would belong to another design.

THE MAPPING, AND WHAT CHECKS IT. The confidences JSON names the chain of every token and
atom; the summary's per-chain arrays do not (index k is the k-th chain by first
appearance in token_chain_ids). Three bridges, each checked:

  token -> residue       token_chain_ids / token_res_ids must equal the CIF's residues
                         (chain, label_seq_id) element by element.
  atom -> residue        atom_chain_ids must equal the CIF's atom chains and
                         plddt_bfac_maxdiff = max |atom_plddts - CIF B-factor| must be
                         <= PLDDT_BFAC_TOL (AF3 writes the same pLDDT into both).
  chain index -> letter  chain_pae_min_maxdiff = max |summary chain_pair_pae_min - the
                         same minimum recomputed from pae over the token_chain_ids
                         blocks| must be <= CHAIN_PAE_MIN_TOL. The ONLY thing tying
                         the summary to a letter, it cannot see a permutation the matrix
                         is invariant to (identical copies, equal minima), and that
                         chain_ptm / chain_pair_iptm share the order is UNCHECKED.
                         chain_pae_min_margin is the smallest maxdiff over every swap of
                         two chains: at or below the tolerance the check could NOT have
                         caught a swap on that row (not an error, or symmetric rows go).

``stem_matches_name`` is 1 when <prefix> is the design name (a 1 does not prove
identity). ``binder_seq_match`` (with --sequence-column) is 1 when the CIF chain named
by --binder-chains has the sequence of the SAME letter's segment of that '/'-chainbreak
column (segment i is chain chr(ord('A') + i)), else 0; NA means UNCHECKED. It ties the
FILE to the column, not the named chain to the binder (``binder_n_res`` does), and is
the one numeric value kept on an error row. ``ranking_score`` is echoed so a filter can
require it to equal alphafold3_ranking_score (weak evidence: AF3 writes it to 0.01).

MEASURED on real AF3 3.0.2 output (one 11-chain homo-oligomer, 2026-10-11; not a binder
complex, and not the 3.0.4 this workspace runs): atom_plddts == CIF B-factors exactly;
chain_pair_pae_min recomputed from pae agrees to <= 0.05 A (rounding) and is 6.3-12.6 A
off with the chain order rotated, on a prediction with no confident interface, where
chain blocks differ most; chain_pair_iptm is exactly symmetric, so the symmetrisation is
a no-op kept for boltzconf's definition.

ipSAE is PORTED, NOT REIMPLEMENTED, and copied rather than imported from boltzconf (a
worker importing a sibling tool dir is not verified on Modal): ``ipsae_d0``,
``ipsae_asym``, ``reject_nucleic_acid`` and their constants are verbatim copies of
tools/boltzconf/boltzconf_worker.py, and test_af3conf.py asserts they are
source-identical. Reference:
https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py (md5, fetch date and
commit in boltzconf_worker.py's docstring).

Normally run per array task by af3conf.sh; run ``sapia run af3conf`` rather than this.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import gemmi
import numpy as np

# Max atom_plddts vs B-factor gap (measured identical; boltzconf's tolerance).
PLDDT_BFAC_TOL = 0.05

# Rounding allowance (A) for chain_pair_pae_min vs its recomputation from pae.
CHAIN_PAE_MIN_TOL = 0.1

# Bare metric names in write order; names shared with boltzconf are identical.
TRUST_COLUMNS = [
    "n_tokens",
    "n_res_struct",
    "chains",
    "binder_n_res",
    "n_target_res",
    "plddt_bfac_maxdiff",
    "chain_pae_min_maxdiff",
    "chain_pae_min_margin",
    "stem_matches_name",
    "binder_seq_match",
    "ranking_score",
    "input_source",
]

# One-letter codes; anything else reads 'X', so it can only make the comparison fail.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
SCIENCE_COLUMNS = [
    "binder_plddt",
    "binder_plddt_min",
    "chains_ptm_binder",
    "iptm_binder_mean",
    "iptm_binder_max",
    "iptm_binder_top2",
    # ipSAE, binder level: over binder-vs-target pairs, from the per-pair max direction.
    "ipsae_binder_mean",
    "ipsae_binder_max",
    "ipsae_binder_top2",
    "ipsae_pae_cutoff",
]

# Per target chain, the ipSAE block -- boltzconf's templates, unchanged.
IPSAE_PER_CHAIN_TEMPLATES = [
    "ipsae_binder_to_{k}",       # frame on the binder, chain k scored
    "ipsae_{k}_to_binder",       # frame on chain k, binder scored
    "ipsae_max_{k}",             # the reference's "max" row for this pair
    "ipsae_n0res_binder_to_{k}",  # support of ipsae_binder_to_{k}
    "ipsae_n0res_{k}_to_binder",  # support of ipsae_{k}_to_binder
    "ipsae_d0res_max_{k}",       # d0 of whichever direction won ipsae_max_{k}
    "ipsae_nres_binder_{k}",     # binder residues with any sub-cutoff PAE to k
    "ipsae_nres_{k}",            # chain-k residues reached from the binder
]


class Args(argparse.Namespace):
    task_file: Path
    binder_chains: str
    target_chains: str
    ipsae_pae_cutoff: float
    input_column: str
    sequence_column: str
    out_dir: Path


class StructureError(ValueError):
    """The file was read but holds nothing usable, or could not be parsed."""


def split_list(value: str) -> list[str]:
    return [tok.strip() for tok in value.split(",") if tok.strip()]


def chain_key(chain: str) -> str:
    """Chain letter -> the suffix used in a per-chain column name ('A' -> 'a')."""
    return "".join(c if c.isalnum() else "_" for c in chain).lower()


def per_chain_columns(target_chains: list[str]) -> list[str]:
    """Every per-target-chain column, grouped by metric, in the flag's order."""
    keys = [chain_key(c) for c in target_chains]
    columns = [f"iptm_{k}" for k in keys]
    for template in IPSAE_PER_CHAIN_TEMPLATES:
        columns += [template.format(k=k) for k in keys]
    return columns


def metric_columns(target_chains: list[str]) -> list[str]:
    return TRUST_COLUMNS + SCIENCE_COLUMNS + per_chain_columns(target_chains)


# --- reading AF3's mmCIF: residues in file order, plus every atom ---


class Residues:
    """Residues of the first model in file order, plus per-atom chain and B-factor.

    File order is the token order AF3 used and the order of atom_plddts. A residue is a
    run of atom records sharing (chain, seq, icode, comp) that also ends when an atom
    name repeats (boltzconf's rule), so repeated numbering never merges residues. """

    def __init__(self) -> None:
        self.chain: list[str] = []            # per residue
        self.comp: list[str] = []             # per residue (ALA, DA, ...)
        self.label_seq: list[int | None] = []  # per residue, None when '.'
        self.ca_atom: list[int | None] = []   # per residue, index into the atoms
        self.chain_order: list[str] = []
        self.atom_chain: list[str] = []       # per atom
        self.atom_bfac: list[float | None] = []  # per atom
        self._key: tuple | None = None
        self._seen: set[tuple[str, str]] = set()

    def add_atom(
        self,
        chain: str,
        comp: str,
        seq: str,
        label_seq: str,
        icode: str,
        atom: str,
        bfac: float | None,
        altloc: str = "",
    ) -> None:
        key = (chain, seq, icode, comp)
        if key != self._key or (atom, altloc) in self._seen:
            self._key = key
            self._seen = set()
            self.chain.append(chain)
            self.comp.append(comp)
            try:
                self.label_seq.append(int(label_seq) if label_seq else None)
            except ValueError:
                self.label_seq.append(None)
            self.ca_atom.append(None)
            if chain not in self.chain_order:
                self.chain_order.append(chain)
            elif self.chain_order[-1] != chain:
                raise StructureError(
                    f"chain {chain} appears in more than one block of the file; "
                    "token order cannot be taken from file order"
                )
        self._seen.add((atom, altloc))
        if atom == "CA" and self.ca_atom[-1] is None:
            self.ca_atom[-1] = len(self.atom_chain)
        self.atom_chain.append(chain)
        self.atom_bfac.append(bfac)

    def __len__(self) -> int:
        return len(self.chain)


def read_residues(path: Path) -> Residues:
    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    if not path.name.lower().endswith(".cif"):
        raise StructureError(
            f"{path.name} is not an mmCIF -- af3conf reads AF3's own "
            f"<prefix>_model.cif, beside which the confidence JSONs sit"
        )
    try:
        block = gemmi.cif.read(str(path)).sole_block()
    except Exception as e:  # noqa: BLE001 - any parse failure is "unreadable"
        raise StructureError(f"unreadable structure {path.name}: {e}") from e

    def col(*tags: str):
        for tag in tags:
            values = block.find_values(f"_atom_site.{tag}")
            if len(values):
                return [gemmi.cif.as_string(v) for v in values]
        return None

    chain_col = col("auth_asym_id", "label_asym_id")
    atom_col = col("label_atom_id", "auth_atom_id")
    comp_col = col("label_comp_id", "auth_comp_id")
    if None in (chain_col, atom_col, comp_col):
        raise StructureError(f"no usable _atom_site loop in {path.name}")
    seq_col = col("auth_seq_id", "label_seq_id")
    label_seq_col = col("label_seq_id")
    icode_col = col("pdbx_PDB_ins_code")
    alt_col = col("label_alt_id")
    bfac_col = col("B_iso_or_equiv")
    model_col = col("pdbx_PDB_model_num")

    def blank(v: str) -> str:
        return "" if v in ("", ".", "?") else v

    out = Residues()
    first_model = model_col[0] if model_col else None
    for i in range(len(chain_col)):
        if model_col is not None and model_col[i] != first_model:
            continue
        bfac = None
        if bfac_col is not None and blank(bfac_col[i]):
            try:
                bfac = float(bfac_col[i])
            except ValueError:
                bfac = None
        out.add_atom(
            chain=chain_col[i],
            comp=comp_col[i],
            seq=blank(seq_col[i]) if seq_col else "",
            label_seq=blank(label_seq_col[i]) if label_seq_col else "",
            icode=blank(icode_col[i]) if icode_col else "",
            atom=atom_col[i],
            bfac=bfac,
            altloc=blank(alt_col[i]) if alt_col else "",
        )
    if not len(out):
        raise StructureError(f"no atoms in the _atom_site loop of {path.name}")
    return out


# --- the sibling files AF3 wrote ---


def prediction_prefix(path: Path) -> str:
    """'x/d1_model.cif' -> 'd1'; 'x/d1_seed-42_sample-0_model.cif' -> 'd1_seed-42_sample-0'."""
    stem = Path(path.name).stem
    if not stem.endswith("_model") or stem == "_model":
        raise ValueError(
            f"{path.name} is not named <prefix>_model.cif, so AF3's sibling "
            f"<prefix>_confidences.json / <prefix>_summary_confidences.json cannot "
            f"be located -- the input column must hold AF3's own predicted structure"
        )
    return stem[: -len("_model")]


def sibling(path: Path, suffix: str) -> Path:
    return path.parent / f"{prediction_prefix(path)}{suffix}"


def load_json(path: Path, what: str) -> dict:
    if not path.exists():
        raise FileNotFoundError(
            f"{path.name} is missing beside the structure -- {what}; nothing is "
            f"guessed from another model or design"
        )
    with open(path) as f:
        payload = json.load(f)
    if not isinstance(payload, dict):
        raise ValueError(f"{path.name} does not hold a JSON object")
    return payload


def load_confidences(path: Path) -> tuple[dict, np.ndarray]:
    """(the confidences JSON without its matrices, the pae as a float array).

    json.load holds pae and contact_probs as Python lists (~32 B per entry), the
    dominant memory cost. contact_probs is dropped at once and pae converted before
    anything else is read; the peak is unchanged (see the skill for the task's memory
    limit). """
    conf = load_json(path, "the PAE matrix and per-atom pLDDT of this model live only there")
    conf.pop("contact_probs", None)
    pae = finite_array(require_key(conf, "pae", path.name), f"{path.name}:pae")
    del conf["pae"]
    return conf, pae


def binder_segment(sequence: str, binder_chain: str, column: str) -> str:
    """The binder's segment of a '/'-chainbreak sequence (segment i is chain 'A' + i).

    prosapia.utils.select_chains' rule, copied not imported (the task image has no
    prosapia); test_af3conf.py asserts the two agree.
    """
    segments = [seg.strip() for seg in sequence.split("/")]
    index = ord(binder_chain.upper()) - ord("A") if len(binder_chain) == 1 else -1
    if not 0 <= index < len(segments):
        raise ValueError(
            f"--sequence-column {column!r} has {len(segments)} '/'-segment(s), so it "
            f"has no segment for binder chain {binder_chain} (segment i is chain "
            f"chr(ord('A') + i), as alphafold3 letters its chains)"
        )
    if not segments[index]:
        raise ValueError(
            f"--sequence-column {column!r} has an EMPTY segment for binder chain "
            f"{binder_chain}"
        )
    return segments[index]


def chain_sequence(residues: "Residues", chain: str) -> str:
    """One-letter sequence of a chain in file order; non-standard residues read 'X'."""
    return "".join(
        THREE_TO_ONE.get(comp.strip().upper(), "X")
        for c, comp in zip(residues.chain, residues.comp)
        if c == chain
    )


def require_key(payload: dict, key: str, source: str):
    if key not in payload:
        raise ValueError(f"{source} has no '{key}'")
    return payload[key]


def finite_array(value, what: str) -> np.ndarray:
    """A numeric array with every entry finite.

    A JSON null becomes NaN and is refused: it would reach every mean and silently
    defeat the mapping checks, since ``nan > tol`` is False. """
    try:
        arr = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as e:
        raise ValueError(f"{what} is not a numeric array: {e}") from e
    if not np.isfinite(arr).all():
        n_bad = int((~np.isfinite(arr)).sum())
        raise ValueError(
            f"{what} holds {n_bad} non-finite or null value(s); every mean over "
            f"them would be NaN and the mapping checks would pass without checking "
            f"anything"
        )
    return arr


def chain_vector(summary: dict, key: str, n_chains: int) -> np.ndarray:
    """A positional per-chain array from the summary, checked against the chain count."""
    arr = finite_array(require_key(summary, key, "the summary JSON"), f"'{key}'")
    if arr.shape != (n_chains,):
        raise ValueError(
            f"'{key}' has shape {arr.shape} but the prediction has {n_chains} chain(s) "
            f"in token_chain_ids -- the summary and the confidences are not the same "
            f"complex, so its positional chain indices cannot be trusted"
        )
    return arr


def chain_matrix(summary: dict, key: str, n_chains: int) -> np.ndarray:
    """A positional chain x chain array from the summary, checked against the count."""
    raw = require_key(summary, key, "the summary JSON")
    if not isinstance(raw, list) or any(not isinstance(r, list) for r in raw):
        raise ValueError(f"'{key}' in the summary JSON is not a list of lists")
    if len(raw) != n_chains or any(len(r) != n_chains for r in raw):
        shape = f"{len(raw)} x {','.join(sorted({str(len(r)) for r in raw})) or '0'}"
        raise ValueError(
            f"'{key}' is {shape} but the prediction has {n_chains} chain(s) in "
            f"token_chain_ids -- the summary and the confidences are not the same "
            f"complex, so its positional chain indices cannot be trusted"
        )
    return finite_array(raw, f"'{key}'")


def pair_value(pairs: np.ndarray, a: int, b: int) -> float:
    """The SYMMETRISED chain_pair_iptm value (pair[a][b] + pair[b][a]) / 2.

    Exactly symmetric on real AF3 3.0.2 output (a no-op there); kept so iptm_<x> has
    boltzconf's definition whatever a later AF3 writes.
    """
    return float((pairs[a, b] + pairs[b, a]) / 2.0)


# --- ipSAE: VERBATIM copies of boltzconf_worker.py (test_af3conf.py asserts parity) ---

# calc_d0_array: L = max(26, L), floored at 1.0 (2.0 only for nucleic acids, refused).
IPSAE_L_MIN = 26.0
IPSAE_D0_MIN = 1.0


# The reference's nucleic-acid residue set; here it drives an outright refusal.
NUCLEIC_COMPS = frozenset({"DA", "DC", "DT", "DG", "A", "C", "U", "G"})


def reject_nucleic_acid(residues: "Residues", chains: list[str], path: Path) -> None:
    """Refuse a named chain that is a nucleic acid, so IPSAE_D0_MIN stays honest.

    IPSAE_D0_MIN is hard-coded to the reference's PROTEIN floor of 1.0, and this is
    the check that makes that an asserted invariant rather than an assumption. The
    reference uses 2.0 for any pair involving a nucleic acid, so a protein/DNA
    complex scored with the 1.0 floor would come back with a systematically
    INFLATED ipSAE -- and nothing else here would notice.

    The token-count check does NOT cover this: Boltz tokenises a standard
    nucleotide as ONE token, exactly like an amino acid, so n_tokens == n_res holds
    for a protein/nucleic-acid complex. (That check fires on ligands and modified
    residues, which are tokenised per atom.) Only the component names distinguish
    them, so only this function can.

    Checked for the NAMED chains only (binder + targets): an unnamed chain enters
    no ipSAE pair and no column here.
    """
    bad: dict[str, str] = {}
    for chain, comp in zip(residues.chain, residues.comp):
        if chain in chains and chain not in bad and comp.strip().upper() in NUCLEIC_COMPS:
            bad[chain] = comp.strip().upper()
    if bad:
        detail = ", ".join(f"{c} (has {bad[c]})" for c in sorted(bad))
        raise ValueError(
            f"chain(s) {detail} in {path.name} look like a nucleic acid. This tool "
            f"computes ipSAE with the reference's PROTEIN d0 floor of "
            f"{IPSAE_D0_MIN}; the reference uses 2.0 for any pair involving a "
            f"nucleic acid, so every ipsae_* column would be inflated and every "
            f"trust column would still read clean. Scoring a protein/nucleic-acid "
            f"interface needs the pair-type branch ported first"
        )


def ipsae_d0(n0: np.ndarray) -> np.ndarray:
    """d0 from a residue count -- a port of calc_d0_array() in ipsae.py.

    At L = 26 the Yang & Skolnick formula gives 0.9577, so the 1.0 floor bites for
    every n0res <= 26: a 2-residue and a 20-residue interface share d0 = 1.0 and
    differ only in how many pairs the mean runs over. That is the reference's
    behaviour, and it is why n0res is a column.
    """
    L = np.maximum(IPSAE_L_MIN, np.asarray(n0, dtype=float))
    return np.maximum(IPSAE_D0_MIN, 1.24 * (L - 15.0) ** (1.0 / 3.0) - 1.8)


def ipsae_asym(
    pae: np.ndarray,
    align_tokens: list[int],
    scored_tokens: list[int],
    pae_cutoff: float,
) -> dict:
    """ipSAE for ONE ORDERED chain pair: frame on `align_tokens`, score the others.

    The port of the reference's ``ipsae_d0res_asym`` branch:

        valid[i][j] = pae[i][j] < cutoff          (RAW pae, strict '<')
        n0res(i)    = number of valid j for that i
        d0(i)       = calc_d0_array(n0res(i))
        byres(i)    = mean over valid j of 1 / (1 + (pae[i][j]/d0(i))**2)
        ipSAE       = max over i of byres(i)      (a MAX OF MEANS)

    `pae` must be the raw asymmetric matrix. Feeding the symmetrised
    ``(pae + pae.T) / 2`` would produce a different, undefined quantity that still
    looks like an ipSAE -- the single most likely way to get this silently wrong.

    Returns ``ipsae``/``n0res``/``d0res`` as None when NO pair is under the cutoff
    (the reference reports 0.0 there; see the module docstring for why this does
    not). ``nres_align``/``nres_scored`` are always real counts, 0 included.
    """
    sub = pae[np.ix_(align_tokens, scored_tokens)]
    valid = sub < pae_cutoff
    n0 = valid.sum(axis=1)
    d0 = ipsae_d0(n0)
    ptm = 1.0 / (1.0 + (sub / d0[:, None]) ** 2.0)
    totals = np.where(valid, ptm, 0.0).sum(axis=1)
    byres = np.divide(
        totals, n0, out=np.zeros(len(align_tokens), dtype=float), where=n0 > 0
    )
    nres_align = int((n0 > 0).sum())
    # Unique residue counts by TOKEN INDEX, not residue number (module docstring).
    nres_scored = int(valid.any(axis=0).sum())
    if nres_align == 0:
        return {
            "ipsae": None,
            "n0res": None,
            "d0res": None,
            "nres_align": 0,
            "nres_scored": 0,
        }
    # np.argmax takes the FIRST maximum, as the reference does. The reference scans
    # the whole-complex array whose non-aligned rows are 0, so restricting to the
    # aligned chain picks the same residue whenever any value is positive -- and
    # when none is, this function has already returned None above.
    best = int(np.argmax(byres))
    return {
        "ipsae": float(byres[best]),
        "n0res": int(n0[best]),
        "d0res": float(d0[best]),
        "nres_align": nres_align,
        "nres_scored": nres_scored,
    }


# --- the measurement ---


def check_tokens(
    residues: Residues,
    token_chains: list,
    token_res: list,
    path: Path,
    named: list[str],
) -> None:
    """The token -> residue map: token i must BE residue i of the CIF.

    Element by element on the chain and the residue number. Every residue of a NAMED
    chain needs an integer label_seq_id (which also refuses a one-atom ligand named as
    binder or target); only an UNNAMED chain's residue without one (a lone ion) is
    compared on chain alone. """
    n_res = len(residues)
    n_tokens = len(token_chains)
    if n_tokens != n_res:
        raise ValueError(
            f"{n_tokens} tokens but {n_res} residues in {path.name} -- one token per "
            f"residue is this tool's premise (every chain a protein; a multi-atom "
            f"ligand or a modified residue is tokenised per atom), so the token "
            f"indices cannot be residue indices here"
        )
    for i, (t_chain, r_chain) in enumerate(zip(token_chains, residues.chain)):
        if str(t_chain) != r_chain:
            raise ValueError(
                f"token {i} is chain {t_chain} in token_chain_ids but residue {i} of "
                f"{path.name} is chain {r_chain} -- the confidences and the structure "
                f"do not list the same residues in the same order"
            )
    named_set = set(named)
    for i, (t_res, r_res) in enumerate(zip(token_res, residues.label_seq)):
        if r_res is None:
            if residues.chain[i] in named_set:
                raise ValueError(
                    f"residue {i} ({residues.comp[i]}) of NAMED chain "
                    f"{residues.chain[i]} in {path.name} has no integer label_seq_id "
                    f"-- it is not a polymer residue (a one-atom ligand such as an "
                    f"ion is one token, so only this catches it), and its token "
                    f"cannot be checked against token_res_ids"
                )
            continue
        if not isinstance(t_res, int) or t_res != r_res:
            raise ValueError(
                f"token {i} (chain {residues.chain[i]}) is residue {t_res!r} in "
                f"token_res_ids but {r_res} (label_seq_id) in {path.name} -- the "
                f"confidences and the structure are numbered differently"
            )


def block_minima(pae: np.ndarray, tokens: dict[str, list[int]], order: list[str]) -> np.ndarray:
    """R[a][b] = min of the pae block (chain order[a] aligned, order[b] scored)."""
    return np.array(
        [[float(pae[np.ix_(tokens[a], tokens[b])].min()) for b in order] for a in order]
    )


def chain_pae_min_maxdiff(recomputed: np.ndarray, reported: np.ndarray) -> float:
    """max |summary chain_pair_pae_min - recomputed block minimum| (identity order).

    The summary's [a][b] is taken to be chains order[a] x order[b]; recomputing it from
    the letter-keyed token blocks is what ties the positional index to a letter.
    """
    return float(np.abs(recomputed - reported).max())


def chain_pae_min_margin(recomputed: np.ndarray, reported: np.ndarray) -> float:
    """The smallest maxdiff over every transposition of two chains.

    How far the nearest single swap sits from the identity order: at or below
    CHAIN_PAE_MIN_TOL a summary with those chains swapped would have passed the order
    check. n(n-1)/2 swaps, each a permutation of the precomputed n x n matrix.
    """
    n = recomputed.shape[0]
    best = float("inf")
    for i in range(n):
        for j in range(i + 1, n):
            p = list(range(n))
            p[i], p[j] = j, i
            best = min(best, float(np.abs(recomputed[np.ix_(p, p)] - reported).max()))
    return best


def measure(
    path: Path,
    binder_chain: str,
    target_chains: list[str],
    ipsae_pae_cutoff: float,
    name: str = "",
    source: str = "",
    sequence: str | None = None,
    sequence_column: str = "",
    partial: dict | None = None,
) -> dict:
    """One prediction -> metrics; raises on anything that makes the numbers meaningless.

    main() records the exception as the status. ``sequence`` is the --sequence-column
    value (None when the flag is absent); ``partial`` receives binder_seq_match as soon
    as it is measured, so a later error can still report it. """
    prefix = prediction_prefix(path)  # fail on a wrong file before parsing it
    residues = read_residues(path)
    present = residues.chain_order
    wanted = [binder_chain, *target_chains]
    absent = [c for c in wanted if c not in present]
    if absent:
        raise ValueError(
            f"chain(s) {','.join(absent)} absent from {path.name} "
            f"(have: {','.join(present) or 'none'})"
        )

    # --- structure <-> sequence: measured first so it survives a later error ---
    seq_match: int | None = None
    if sequence is not None:
        if not sequence.strip():
            raise ValueError(
                f"--sequence-column {sequence_column!r} could not be resolved for this "
                f"design (blank own cell, or absent here and lineage not allowed), so "
                f"the binder_seq_match check that was asked for cannot be made"
            )
        expected = binder_segment(sequence.strip().upper(), binder_chain, sequence_column)
        seq_match = int(chain_sequence(residues, binder_chain) == expected)
        if partial is not None:
            partial["binder_seq_match"] = seq_match

    reject_nucleic_acid(residues, wanted, path)

    conf_name = f"{prefix}_confidences.json"
    conf, pae = load_confidences(sibling(path, "_confidences.json"))
    summary = load_json(
        sibling(path, "_summary_confidences.json"),
        "the per-chain confidence of this model lives only there",
    )
    ranking_score = finite_array(
        require_key(summary, "ranking_score", "the summary JSON"), "'ranking_score'"
    )
    if ranking_score.shape != ():
        raise ValueError("'ranking_score' in the summary JSON is not one number")

    # --- token -> residue ---------------------------------------------------
    token_chains = require_key(conf, "token_chain_ids", conf_name)
    token_res = require_key(conf, "token_res_ids", conf_name)
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise ValueError(f"the pae array is {pae.shape}, not a square matrix")
    n_tokens = int(pae.shape[0])
    if not isinstance(token_chains, list) or not isinstance(token_res, list):
        raise ValueError(f"token_chain_ids / token_res_ids in {conf_name} are not lists")
    if len(token_chains) != n_tokens or len(token_res) != n_tokens:
        raise ValueError(
            f"pae is {n_tokens} x {n_tokens} but token_chain_ids has "
            f"{len(token_chains)} entries and token_res_ids {len(token_res)}"
        )
    check_tokens(residues, token_chains, token_res, path, wanted)
    order = list(dict.fromkeys(str(c) for c in token_chains))  # == present, now

    # --- atom -> residue (the pLDDT path) ----------------------------------
    atom_chains = require_key(conf, "atom_chain_ids", conf_name)
    atom_plddts = finite_array(
        require_key(conf, "atom_plddts", conf_name), f"{conf_name}:atom_plddts"
    ).reshape(-1)
    n_atoms = len(residues.atom_chain)
    if not isinstance(atom_chains, list) or len(atom_chains) != n_atoms or atom_plddts.size != n_atoms:
        raise ValueError(
            f"{path.name} has {n_atoms} atoms but {conf_name} has "
            f"{len(atom_chains) if isinstance(atom_chains, list) else '?'} "
            f"atom_chain_ids and {atom_plddts.size} atom_plddts"
        )
    if [str(c) for c in atom_chains] != residues.atom_chain:
        raise ValueError(
            f"atom_chain_ids in {conf_name} do not follow the atom order of "
            f"{path.name} -- the per-atom pLDDT would be read off the wrong atoms"
        )
    if atom_plddts.min() < 0.0 or atom_plddts.max() > 100.0 or atom_plddts.max() <= 1.0:
        raise ValueError(
            f"atom_plddts span {atom_plddts.min():.3f}-{atom_plddts.max():.3f}; AF3 "
            f"writes pLDDT on a 0-100 scale, so this is not the file assumed"
        )
    if any(b is None for b in residues.atom_bfac):
        raise ValueError(
            f"{path.name} has atoms without a B-factor; AF3 writes pLDDT there, and "
            f"it is the independent check of the atom order"
        )
    bfac = np.asarray(residues.atom_bfac, dtype=float)
    plddt_maxdiff = float(np.abs(atom_plddts - bfac).max())
    if plddt_maxdiff > PLDDT_BFAC_TOL:
        raise ValueError(
            f"atom_plddts and the structure's B-factors disagree by up to "
            f"{plddt_maxdiff:.2f} (tolerance {PLDDT_BFAC_TOL}) -- AF3 writes the same "
            f"per-atom pLDDT into both, so the atom orders are not the same here"
        )

    # --- chain index -> letter ---------------------------------------------
    n_chains = len(order)
    chain_ptm = chain_vector(summary, "chain_ptm", n_chains)
    pairs = chain_matrix(summary, "chain_pair_iptm", n_chains)
    pae_min = chain_matrix(summary, "chain_pair_pae_min", n_chains)

    tokens: dict[str, list[int]] = {chain: [] for chain in order}
    for token, chain in enumerate(token_chains):
        tokens[str(chain)].append(token)
    recomputed = block_minima(pae, tokens, order)
    pae_min_diff = chain_pae_min_maxdiff(recomputed, pae_min)
    if pae_min_diff > CHAIN_PAE_MIN_TOL:
        raise ValueError(
            f"chain_pair_pae_min in the summary disagrees by up to {pae_min_diff:.2f} A "
            f"with the minimum recomputed from the pae matrix over the "
            f"token_chain_ids blocks (tolerance {CHAIN_PAE_MIN_TOL}) -- the summary's "
            f"positional chain indices are not the chains in token order, so every "
            f"chain_ptm / chain_pair_iptm value would belong to another chain. A "
            f"systematic order difference errors like this on non-confident rows "
            f"but can pass confident ones: distrust this batch's OK rows too"
        )
    pae_min_margin = chain_pae_min_margin(recomputed, pae_min)

    index_of = {chain: i for i, chain in enumerate(order)}
    binder_idx = index_of[binder_chain]
    binder_tokens = tokens[binder_chain]
    target_tokens = {c: tokens[c] for c in target_chains}
    all_target_tokens = [t for c in target_chains for t in target_tokens[c]]

    # Per-residue pLDDT is the CA atom's, ipsae.py's convention (boltzconf: per token).
    no_ca = [i for i in binder_tokens if residues.ca_atom[i] is None]
    if no_ca:
        raise ValueError(
            f"{len(no_ca)} residue(s) of binder chain {binder_chain} have no CA atom "
            f"(first: {residues.comp[no_ca[0]]} {residues.label_seq[no_ca[0]]}) -- "
            f"their pLDDT is undefined"
        )
    binder_plddt = atom_plddts[[residues.ca_atom[i] for i in binder_tokens]]

    per_target_iptm = [pair_value(pairs, binder_idx, index_of[c]) for c in target_chains]

    metrics: dict[str, object] = {
        "n_tokens": n_tokens,
        "n_res_struct": len(residues),
        "chains": f"{binder_chain}|{','.join(target_chains)}",
        "binder_n_res": len(binder_tokens),
        "n_target_res": len(all_target_tokens),
        "plddt_bfac_maxdiff": plddt_maxdiff,
        "chain_pae_min_maxdiff": pae_min_diff,
        "chain_pae_min_margin": pae_min_margin,
        # Siblings match by PATH; a 1 does not prove identity.
        "stem_matches_name": int(bool(name) and prefix == name),
        "binder_seq_match": seq_match,
        "ranking_score": float(ranking_score),
        "input_source": source,
        "binder_plddt": float(binder_plddt.mean()),
        "binder_plddt_min": float(binder_plddt.min()),
        "chains_ptm_binder": float(chain_ptm[binder_idx]),
        # Three numbers: an untouched target protomer drags the mean down.
        "iptm_binder_mean": float(np.mean(per_target_iptm)),
        "iptm_binder_max": float(np.max(per_target_iptm)),
        # With one target chain there is no second value, so top2 IS the max.
        "iptm_binder_top2": float(np.mean(sorted(per_target_iptm, reverse=True)[:2])),
        "ipsae_pae_cutoff": float(ipsae_pae_cutoff),
    }
    # Per target chain, in flag order; ipSAE takes the RAW pae, never symmetrised.
    ipsae_pair_max: list[float] = []  # reference convention: 0.0 for "no pair"
    any_ipsae_pair = False
    for chain in target_chains:
        key = chain_key(chain)
        metrics[f"iptm_{key}"] = per_target_iptm[target_chains.index(chain)]

        forward = ipsae_asym(pae, binder_tokens, target_tokens[chain], ipsae_pae_cutoff)
        reverse = ipsae_asym(pae, target_tokens[chain], binder_tokens, ipsae_pae_cutoff)
        available = [v for v in (forward["ipsae"], reverse["ipsae"]) if v is not None]
        if available:
            any_ipsae_pair = True
            pair_max = max(available)
            # On an exact tie the binder->target direction wins, as in boltzconf.
            won_forward = forward["ipsae"] is not None and forward["ipsae"] >= pair_max
            d0res_max = forward["d0res"] if won_forward else reverse["d0res"]
        else:
            pair_max = None
            d0res_max = None
        ipsae_pair_max.append(0.0 if pair_max is None else pair_max)

        metrics[f"ipsae_binder_to_{key}"] = forward["ipsae"]
        metrics[f"ipsae_{key}_to_binder"] = reverse["ipsae"]
        metrics[f"ipsae_max_{key}"] = pair_max
        metrics[f"ipsae_n0res_binder_to_{key}"] = forward["n0res"]
        metrics[f"ipsae_n0res_{key}_to_binder"] = reverse["n0res"]
        metrics[f"ipsae_d0res_max_{key}"] = d0res_max
        # Real counts: 0 means no binder residue is confident about this chain.
        metrics[f"ipsae_nres_binder_{key}"] = forward["nres_align"]
        metrics[f"ipsae_nres_{key}"] = forward["nres_scored"]

    # No-pair chains count 0.0, as in the reference; NA if no chain has a pair.
    metrics["ipsae_binder_mean"] = float(np.mean(ipsae_pair_max)) if any_ipsae_pair else None
    metrics["ipsae_binder_max"] = float(max(ipsae_pair_max)) if any_ipsae_pair else None
    metrics["ipsae_binder_top2"] = (
        float(np.mean(sorted(ipsae_pair_max, reverse=True)[:2])) if any_ipsae_pair else None
    )
    return metrics


# --- output ---


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def write_result(
    out_dir: Path,
    name: str,
    status: str,
    path: str,
    metrics: dict,
    columns: list[str],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    # Write to a temporary and rename: a killed task must not leave a half-row.
    tmp = out_dir / f".{name}.tsv.tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["name", "status", "path", *columns])
        w.writerow([name, status, path, *(_fmt(metrics.get(c)) for c in columns)])
    tmp.replace(out_dir / f"{name}.tsv")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Binder-specific AlphaFold3 confidence and ipSAE."
    )
    ap.add_argument("--task-file", type=Path, required=True,
                    help="Tab-separated 'name<TAB>structure[<TAB>source table"
                         "[<TAB>sequence]]' rows.")
    ap.add_argument("--binder-chains", required=True,
                    help="The binder chain. Exactly one (see the module docstring).")
    ap.add_argument("--target-chains", required=True,
                    help="Comma-joined target chain IDs, in column order.")
    ap.add_argument("--ipsae-pae-cutoff", type=float, default=15.0,
                    help="PAE cutoff in A for the ipsae_* columns (default 15.0, as "
                         "boltzconf). Two ipSAEs at different cutoffs are not "
                         "comparable.")
    ap.add_argument("--input-column", default="the input structure column",
                    help="Only used to name the column in the 'could not resolve' error.")
    ap.add_argument("--sequence-column", default="",
                    help="The --sequence-column the builder resolved into the task "
                         "file's 4th field. Empty (the default) = no sequence check, "
                         "binder_seq_match NA.")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(namespace=Args())

    binder = split_list(args.binder_chains)
    targets = split_list(args.target_chains)
    if len(binder) != 1:
        raise ValueError(
            f"--binder-chains must name exactly one chain, got {len(binder)} "
            f"({args.binder_chains!r}): chain_ptm and chain_pair_iptm are per "
            f"chain, so two binder chains could only be averaged, and an average "
            f"reading like one chain's confidence is the failure this refuses"
        )
    if not targets:
        raise ValueError("--target-chains is empty")
    overlap = sorted(set(binder) & set(targets))
    if overlap:
        raise ValueError(f"chain(s) {','.join(overlap)} are both binder and target")
    if args.ipsae_pae_cutoff <= 0:
        raise ValueError("--ipsae-pae-cutoff must be > 0")
    columns = metric_columns(targets)

    with open(args.task_file) as f:
        members = [line.rstrip("\n").split("\t") for line in f if line.strip()]

    for member in members:
        if len(member) < 2:
            raise ValueError(
                f"malformed task row {member!r}: expected 'name<TAB>structure'"
            )
        name, src = member[0], member[1]
        source = member[2] if len(member) > 2 else ""
        sequence = (member[3] if len(member) > 3 else "") if args.sequence_column else None
        partial: dict = {}
        try:
            if not src.strip():
                raise ValueError(
                    f"could not resolve a structure path from column "
                    f"{args.input_column!r} for this design -- neither on its own row "
                    f"nor on any ancestor (parent_name/parent_table)"
                )
            metrics = measure(
                Path(src), binder[0], targets, args.ipsae_pae_cutoff,
                name=name, source=source, sequence=sequence,
                sequence_column=args.sequence_column, partial=partial,
            )
            write_result(args.out_dir, name, "OK", src, metrics, columns)
            print(
                f"{name}: binder_plddt={metrics['binder_plddt']:.1f} "
                f"iptm_binder_max={metrics['iptm_binder_max']:.3f} "
                f"chain_pae_min_maxdiff={metrics['chain_pae_min_maxdiff']:.3f} "
                "ipsae_binder_max="
                + ("NA" if metrics["ipsae_binder_max"] is None
                   else f"{metrics['ipsae_binder_max']:.3f}")
            )
        except Exception as e:  # noqa: BLE001 - errors are recorded as data
            msg = " ".join(str(e).split())  # one line: it lands in a TSV cell
            print(f"{name}: ERROR {msg}", file=sys.stderr)
            # The echo and binder_seq_match survive; every other numeric column is NA.
            write_result(
                args.out_dir, name, f"error: {msg}", "",
                {"chains": f"{binder[0]}|{','.join(targets)}",
                 "input_source": source, **partial}, columns,
            )


if __name__ == "__main__":
    main()
