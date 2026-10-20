#!/usr/bin/env python3
"""
Per-array-task step of the boltzconf tool: the binder-specific confidence of ONE boltz
prediction, computed from the files boltz wrote beside the predicted structure (nothing
is re-read from the table or re-predicted). How to run it and how to read its columns:
``.claude/skills/boltzconf/SKILL.md``.

A design is an ERROR ROW, never a default or a trimmed fit, unless:

  * there is exactly ONE binder chain and at least one disjoint target chain
    (``chains_ptm`` and ``pair_chains_iptm`` are per chain, so two binder chains could
    only be an average that reads like one chain's value);
  * the three sibling files below are all present;
  * there is one token per residue (a ligand or modified residue is tokenised per atom);
  * the named chains are protein (``reject_nucleic_acid``, by component name): boltz
    gives a standard nucleotide ONE token, so the token check cannot see one, and the d0
    floor here is the reference's PROTEIN 1.0 (2.0 for nucleic acids), which would
    inflate every ipsae_* column while the trust columns stayed clean.

SIBLINGS. The input column names ``<...>/predictions/<name>/<name>_model_0.cif``;
everything else is derived from that file's own stem (``<name>_model_0``) in the same
directory, never from the design name (they differ when boltz sanitised the name):

    confidence_<stem>.json   chains_ptm, pair_chains_iptm (N x N, NOT symmetric)
    pae_<stem>.npz           key 'pae', n_tokens x n_tokens, Angstrom
    plddt_<stem>.npz         key 'plddt', n_tokens, 0-1 (rescaled here to 0-100)

THE MAPPING. The JSON keys chains by POSITIONAL INDEX in the input YAML and the PAE and
pLDDT arrays are indexed by TOKEN; neither carries a chain letter, so the only bridge is
the chain order of the structure file, assumed to be the order boltz was given. If that
is wrong every number is computed on the wrong residues and still looks reasonable.
Hence:

    n_tokens == n_res_struct, and len(pair_chains_iptm) == the structure's chain count,
        or it is an ERROR (never a trimmed fit);
    plddt_bfac_maxdiff = max |npz pLDDT*100 - CIF B-factor| over all tokens (boltz
        writes the same per-residue pLDDT into both): ~1e-5 when the token -> residue
        map is right, above PLDDT_BFAC_TOL an ERROR. It proves the npz and CIF orders
        agree, NOT the JSON's chain-index -> letter map (only the chain COUNT is checked
        there);
    stem_matches_name is 1 when the file is <name>_model_0 (a 0 means matched by PATH
        alone; not an error); input_source is the table the path was read from.

Interface pairs are binder residue i and target residue j whose closest HEAVY atoms are
within ``--contact-cutoff``; hydrogens are ignored.

ipSAE -- PORTED, NOT REIMPLEMENTED. The ``ipsae_*`` columns transcribe Roland Dunbrack's
reference implementation. This is the exact file ported, so a later editor re-reads it
rather than their memory of it:

    https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py
    md5 a48df7adc64afa36b4b475429b0a8aff (45140 bytes, 1010 lines), fetched
    2026-10-10; header "version 4 / January 3, 2026: Fixed Boltz2 issues";
    repo commit for that path 6174cf9e71cb1bd660cc805856a18c4871a6dec3
    paper: https://www.biorxiv.org/content/10.1101/2025.02.10.637595v2

Nothing is reconstructed from memory: the d0 clamp, the strict '<' on the RAW asymmetric
PAE, the per-residue normalisation and the max-of-means reduction are read off ipsae.py,
and changing any silently yields a number that looks like an ipSAE and is not one.
``ipsae_d0``, ``ipsae_asym``, ``reject_nucleic_acid`` and their constants are copied
verbatim into tools/af3conf/af3conf_worker.py, whose test asserts they are
source-identical: edit both or neither.

Normally run per array task by boltzconf.sh; run ``sapia run boltzconf`` rather than
this. """

import argparse
import csv
import gzip
import json
import statistics
import sys
from pathlib import Path

import gemmi
import numpy as np

# Max pLDDT vs B-factor gap: tight on purpose, so a shift inside a flat helix shows.
PLDDT_BFAC_TOL = 0.05

# A pLDDT max at or below this is on the 0-1 scale (the B-factor check proves it).
PLDDT_UNIT_MAX = 1.5

_CIF_SUFFIXES = (".cif", ".mmcif", ".cif.gz", ".mmcif.gz")

# Bare metric names in write order; collect_boltzconf.py hands them to the driver.
TRUST_COLUMNS = [
    "n_tokens",
    "n_res_struct",
    "chains",
    "binder_n_res",
    "n_target_res",
    "plddt_scale",
    "plddt_bfac_maxdiff",
    "stem_matches_name",
    "input_source",
]
SCIENCE_COLUMNS = [
    "binder_plddt",
    "binder_plddt_min",
    "chains_ptm_binder",
    "iptm_binder_mean",
    "iptm_binder_max",
    "iptm_binder_top2",
    "ipae_mean",
    "ipae_median",
    "n_if_pairs",
    "n_if_res_binder",
    "pae_binder_self",
    "contact_cutoff",
    # ipSAE, binder level: over binder-vs-target pairs, from the per-pair max direction.
    "ipsae_binder_mean",
    "ipsae_binder_max",
    "ipsae_binder_top2",
    "ipsae_pae_cutoff",
]

# Per target chain, the ipSAE block; both directions are kept (it is asymmetric).
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
    contact_cutoff: float
    ipsae_pae_cutoff: float
    input_column: str
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
    columns = [f"iptm_{k}" for k in keys] + [f"ipae_{k}_mean" for k in keys]
    for template in IPSAE_PER_CHAIN_TEMPLATES:
        columns += [template.format(k=k) for k in keys]
    return columns


def metric_columns(target_chains: list[str]) -> list[str]:
    return TRUST_COLUMNS + SCIENCE_COLUMNS + per_chain_columns(target_chains)


# --- reading the structure: residues in FILE ORDER, heavy atoms only ---


class Residues:
    """Residues of the first model in file order -- the token order Boltz used.

    A residue is a run of consecutive atom records sharing (chain, seq, icode, comp),
    and a run also ends when an atom name repeats. Runs, because a generator that
    restarts numbering inside a chain (RFdiffusion's symmetric mode writes 1,2,3,1,2,3)
    would otherwise collapse residues; the repeated-name rule, because three GLY all
    numbered 1 are indistinguishable by key and the token count would come out short.
    """

    def __init__(self) -> None:
        self.chain: list[str] = []          # per residue
        self.comp: list[str] = []           # per residue, the component (ALA, DA, ...)
        self.bfactor: list[float | None] = []
        self.coords: list[list[tuple[float, float, float]]] = []  # heavy atoms
        self.chain_order: list[str] = []
        self._key: tuple | None = None
        self._seen: set[tuple[str, str]] = set()

    def add_atom(
        self,
        chain: str,
        comp: str,
        seq: str,
        icode: str,
        atom: str,
        element: str,
        bfac: float | None,
        xyz: tuple[float, float, float],
        altloc: str = "",
    ) -> None:
        key = (chain, seq, icode, comp)
        if key != self._key or (atom, altloc) in self._seen:
            self._key = key
            self._seen = set()
            self.chain.append(chain)
            self.comp.append(comp)
            self.bfactor.append(None)
            self.coords.append([])
            if chain not in self.chain_order:
                self.chain_order.append(chain)
            elif self.chain_order[-1] != chain:
                raise StructureError(
                    f"chain {chain} appears in more than one block of the file; "
                    "token order cannot be taken from file order"
                )
        self._seen.add((atom, altloc))
        is_hydrogen = (
            element.upper() in ("H", "D")
            if element
            else atom[:1] in ("H", "D") and not atom[:1].isdigit()
        )
        if not is_hydrogen:
            self.coords[-1].append(xyz)
        # The B-factor holds per-residue pLDDT; CA is preferred so one atom is compared.
        if bfac is not None and (self.bfactor[-1] is None or atom == "CA"):
            self.bfactor[-1] = bfac

    def __len__(self) -> int:
        return len(self.chain)


def _read_text(path: Path) -> str:
    if path.name.lower().endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return f.read()
    return path.read_text()


def _is_cif(path: Path) -> bool:
    return path.name.lower().endswith(_CIF_SUFFIXES)


def _read_pdb(path: Path) -> Residues:
    out = Residues()
    n_atoms = 0
    for line in _read_text(path).splitlines():
        record = line[:6]
        if record == "ENDMDL" and n_atoms:
            break  # first model only
        if record not in ("ATOM  ", "HETATM"):
            continue
        if len(line) < 54:
            raise StructureError(f"truncated {record.strip()} record: {line!r}")
        n_atoms += 1
        bfac = None
        if len(line) >= 66 and line[60:66].strip():
            try:
                bfac = float(line[60:66])
            except ValueError:
                bfac = None
        out.add_atom(
            chain=line[20:22].strip(),
            comp=line[17:21].strip(),
            seq=line[22:26].strip(),
            icode=line[26].strip(),
            atom=line[12:16].strip(),
            element=line[76:78].strip() if len(line) >= 78 else "",
            bfac=bfac,
            xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
            altloc=line[16].strip(),
        )
    if n_atoms == 0:
        raise StructureError(f"no ATOM/HETATM records in {path.name}")
    return out


def _read_cif(path: Path) -> Residues:
    block = gemmi.cif.read(str(path)).sole_block()

    def col(*tags: str):
        for tag in tags:
            values = block.find_values(f"_atom_site.{tag}")
            if len(values):
                return [gemmi.cif.as_string(v) for v in values]
        return None

    x, y, z = col("Cartn_x"), col("Cartn_y"), col("Cartn_z")
    chain_col = col("auth_asym_id", "label_asym_id")
    atom_col = col("auth_atom_id", "label_atom_id")
    comp_col = col("auth_comp_id", "label_comp_id")
    if None in (x, y, z, chain_col, atom_col, comp_col):
        raise StructureError(f"no usable _atom_site loop in {path.name}")
    seq_col = col("auth_seq_id", "label_seq_id")
    icode_col = col("pdbx_PDB_ins_code")
    elem_col = col("type_symbol")
    alt_col = col("label_alt_id")
    bfac_col = col("B_iso_or_equiv")
    model_col = col("pdbx_PDB_model_num")

    def blank(v: str) -> str:
        return "" if v in ("", ".", "?") else v

    out = Residues()
    first_model = model_col[0] if model_col else None
    for i in range(len(x)):
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
            icode=blank(icode_col[i]) if icode_col else "",
            atom=atom_col[i],
            element=blank(elem_col[i]) if elem_col else "",
            bfac=bfac,
            xyz=(float(x[i]), float(y[i]), float(z[i])),
            altloc=blank(alt_col[i]) if alt_col else "",
        )
    if not len(out):
        raise StructureError(f"no atoms in the _atom_site loop of {path.name}")
    return out


def read_residues(path: Path) -> Residues:
    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    try:
        return _read_cif(path) if _is_cif(path) else _read_pdb(path)
    except StructureError:
        raise
    except Exception as e:  # noqa: BLE001 - any parse failure is "unreadable"
        raise StructureError(f"unreadable structure {path.name}: {e}") from e


# --- the sibling files Boltz wrote ---


def structure_stem(path: Path) -> str:
    """'x/d1_model_0.cif' -> 'd1_model_0' (also '.cif.gz')."""
    name = path.name
    for suffix in (".gz",):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    return Path(name).stem


def sibling(path: Path, prefix: str, suffix: str) -> Path:
    return path.parent / f"{prefix}_{structure_stem(path)}{suffix}"


def require_finite(values: np.ndarray, what: str) -> np.ndarray:
    """Refuse NaN/inf.

    A NaN propagates into every mean and silently defeats the mapping check, because
    ``nan > tol`` is False.
    """
    if not np.isfinite(values).all():
        n_bad = int((~np.isfinite(values)).sum())
        raise ValueError(
            f"{what} holds {n_bad} non-finite value(s) (NaN/inf); every mean over "
            f"them would be NaN and the pLDDT mapping check would pass without "
            f"checking anything"
        )
    return values


def load_npz(path: Path, key: str) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(
            f"{path.name} is missing beside the structure -- boltz writes it next to "
            f"the predicted CIF; nothing is guessed from another model or design"
        )
    with np.load(str(path)) as data:
        if key not in data:
            raise ValueError(
                f"{path.name} has no '{key}' array (keys: {','.join(data.files) or 'none'})"
            )
        return require_finite(np.asarray(data[key], dtype=float), f"{path.name}:{key}")


def load_confidence(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(
            f"{path.name} is missing beside the structure -- the per-chain confidence "
            f"of this model lives only there"
        )
    with open(path) as f:
        payload = json.load(f)
    if not isinstance(payload, dict):
        raise ValueError(f"{path.name} does not hold a JSON object")
    return payload


def chain_index_table(payload: dict, key: str, n_chains: int) -> dict[int, object]:
    """A positional-chain-index map from the JSON, checked against the structure.

    The JSON keys chains "0", "1", ... in the order of boltz's input YAML and names
    none, so a count that disagrees with the structure means a different complex. """
    if key not in payload:
        raise ValueError(f"the confidence JSON has no '{key}'")
    table = payload[key]
    if not isinstance(table, dict):
        raise ValueError(f"'{key}' in the confidence JSON is not an object")
    if len(table) != n_chains:
        raise ValueError(
            f"'{key}' covers {len(table)} chain(s) but the structure has {n_chains} "
            f"-- the JSON and the CIF are not the same complex, so the positional "
            f"chain indices cannot be trusted"
        )
    out: dict[int, object] = {}
    for raw_key, value in table.items():
        try:
            out[int(raw_key)] = value
        except (TypeError, ValueError) as e:
            raise ValueError(f"'{key}' has a non-integer chain key {raw_key!r}") from e
        for leaf in (value.values() if isinstance(value, dict) else [value]):
            try:
                finite = np.isfinite(float(leaf))  # type: ignore[arg-type]
            except (TypeError, ValueError) as e:
                raise ValueError(
                    f"'{key}' chain {raw_key} holds a non-numeric value {leaf!r}"
                ) from e
            if not finite:
                raise ValueError(
                    f"'{key}' chain {raw_key} holds a non-finite value {leaf!r}"
                )
    missing = [i for i in range(n_chains) if i not in out]
    if missing:
        raise ValueError(
            f"'{key}' is missing chain index/indices {','.join(map(str, missing))}"
        )
    return out


def pair_value(pairs: dict[int, object], a: int, b: int) -> float:
    """The SYMMETRISED pair_chains_iptm value (pair[a][b] + pair[b][a]) / 2.

    Boltz's matrix is not symmetric (up to ~0.08 apart) and neither direction is more
    correct, so both are used.
    """
    def one(i: int, j: int) -> float:
        row = pairs[i]
        if not isinstance(row, dict):
            raise ValueError("'pair_chains_iptm' rows are not objects")
        for candidate in (j, str(j)):
            if candidate in row:
                return float(row[candidate])  # type: ignore[index]
        raise ValueError(f"'pair_chains_iptm' has no entry [{i}][{j}]")

    return (one(a, b) + one(b, a)) / 2.0


# --- ipSAE: a transcription of the reference (source and md5: module docstring) ---

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


def residue_min_distances(
    residues: Residues, rows: list[int], cols: list[int]
) -> np.ndarray:
    """Minimum heavy-atom distance for every (rows[i], cols[j]) residue pair."""
    row_atoms = [np.asarray(residues.coords[i], dtype=float) for i in rows]
    col_atoms = [np.asarray(residues.coords[j], dtype=float) for j in cols]
    empty_rows = [rows[i] for i, a in enumerate(row_atoms) if a.size == 0]
    empty_cols = [cols[j] for j, a in enumerate(col_atoms) if a.size == 0]
    if empty_rows or empty_cols:
        raise ValueError(
            f"{len(empty_rows) + len(empty_cols)} residue(s) have no heavy atoms "
            f"(token indices {','.join(map(str, (empty_rows + empty_cols)[:5]))}) -- "
            "a distance to them is undefined"
        )
    col_counts = [len(a) for a in col_atoms]
    b = np.concatenate(col_atoms, axis=0)
    col_starts = np.cumsum([0, *col_counts[:-1]])
    out = np.empty((len(rows), len(cols)), dtype=float)
    # Blocked over whole ROW RESIDUES (never atoms): no residue straddles two blocks.
    block_size = 64
    for start in range(0, len(rows), block_size):
        chunk = row_atoms[start : start + block_size]
        counts = [len(a) for a in chunk]
        a = np.concatenate(chunk, axis=0)
        d2 = ((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=-1)
        per_col = np.minimum.reduceat(d2, col_starts, axis=1)
        row_starts = np.cumsum([0, *counts[:-1]])
        out[start : start + len(chunk)] = np.minimum.reduceat(
            per_col, row_starts, axis=0
        )
    return np.sqrt(out)


def measure(
    path: Path,
    binder_chain: str,
    target_chains: list[str],
    cutoff: float,
    ipsae_pae_cutoff: float,
    name: str = "",
    source: str = "",
) -> dict:
    """One prediction -> metrics; raises on anything that makes the numbers meaningless.

    main() records the exception as the design's status.
    """
    residues = read_residues(path)
    present = residues.chain_order
    wanted = [binder_chain, *target_chains]
    absent = [c for c in wanted if c not in present]
    if absent:
        raise ValueError(
            f"chain(s) {','.join(absent)} absent from {path.name} "
            f"(have: {','.join(present) or 'none'})"
        )
    reject_nucleic_acid(residues, wanted, path)

    confidence = load_confidence(sibling(path, "confidence", ".json"))
    pae = load_npz(sibling(path, "pae", ".npz"), "pae")
    plddt_raw = load_npz(sibling(path, "plddt", ".npz"), "plddt")

    n_res = len(residues)
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise ValueError(f"the pae array is {pae.shape}, not a square matrix")
    n_tokens = int(pae.shape[0])
    if n_tokens != n_res:
        raise ValueError(
            f"{n_tokens} pae tokens but {n_res} residues in {path.name} -- one token "
            f"per residue is this tool's premise (every chain a protein), so the "
            f"token indices cannot be residue indices here"
        )
    plddt = plddt_raw.reshape(-1)
    if plddt.size != n_tokens:
        raise ValueError(
            f"the plddt array has {plddt.size} entries but there are {n_tokens} tokens"
        )
    scale = "0-1" if float(plddt.max(initial=0.0)) <= PLDDT_UNIT_MAX else "0-100"
    plddt = plddt * 100.0 if scale == "0-1" else plddt

    # The independent check of the token -> residue mapping.
    bfac = [b for b in residues.bfactor]
    maxdiff: float | None = None
    if all(b is not None for b in bfac) and any(abs(b) > 0 for b in bfac):  # type: ignore[arg-type]
        maxdiff = float(np.abs(plddt - np.asarray(bfac, dtype=float)).max())
        if maxdiff > PLDDT_BFAC_TOL:
            raise ValueError(
                f"the plddt npz and the structure's B-factors disagree by up to "
                f"{maxdiff:.2f} (tolerance {PLDDT_BFAC_TOL}) -- boltz writes the same "
                f"per-residue pLDDT into both, so the token order and the residue "
                f"order are not the same here and every number would be computed on "
                f"the wrong residues"
            )

    pairs = chain_index_table(confidence, "pair_chains_iptm", len(present))
    chains_ptm = chain_index_table(confidence, "chains_ptm", len(present))
    index_of = {chain: i for i, chain in enumerate(present)}
    binder_idx = index_of[binder_chain]

    tokens = {chain: [] for chain in present}
    for token, chain in enumerate(residues.chain):
        tokens[chain].append(token)
    binder_tokens = tokens[binder_chain]
    target_tokens = {c: tokens[c] for c in target_chains}
    all_target_tokens = [t for c in target_chains for t in target_tokens[c]]

    per_target_iptm = [
        pair_value(pairs, binder_idx, index_of[c]) for c in target_chains
    ]
    pae_sym = (pae + pae.T) / 2.0
    distances = residue_min_distances(residues, binder_tokens, all_target_tokens)
    contact = distances <= cutoff
    binder_rows, target_cols = np.nonzero(contact)
    pair_pae = [
        float(pae_sym[binder_tokens[i], all_target_tokens[j]])
        for i, j in zip(binder_rows.tolist(), target_cols.tolist())
    ]

    metrics: dict[str, object] = {
        "n_tokens": n_tokens,
        "n_res_struct": n_res,
        "chains": f"{binder_chain}|{','.join(target_chains)}",
        "binder_n_res": len(binder_tokens),
        "n_target_res": len(all_target_tokens),
        "plddt_scale": scale,
        "plddt_bfac_maxdiff": maxdiff,
        # Siblings match by PATH; boltz sanitises names, so a 0 is shown, not an error.
        "stem_matches_name": int(bool(name) and structure_stem(path) == f"{name}_model_0"),
        "input_source": source,
        "binder_plddt": float(plddt[binder_tokens].mean()),
        "binder_plddt_min": float(plddt[binder_tokens].min()),
        "chains_ptm_binder": float(chains_ptm[binder_idx]),  # type: ignore[arg-type]
        # Three numbers: an untouched target protomer drags the mean down.
        "iptm_binder_mean": float(np.mean(per_target_iptm)),
        "iptm_binder_max": float(np.max(per_target_iptm)),
        # One target chain: no second value, so top2 IS the max.
        "iptm_binder_top2": float(np.mean(sorted(per_target_iptm, reverse=True)[:2])),
        "ipae_mean": float(np.mean(pair_pae)) if pair_pae else None,
        "ipae_median": float(statistics.median(pair_pae)) if pair_pae else None,
        "n_if_pairs": len(pair_pae),
        "n_if_res_binder": int(len(set(binder_rows.tolist()))),
        "pae_binder_self": _self_pae(pae_sym, binder_tokens),
        "contact_cutoff": cutoff,
        "ipsae_pae_cutoff": float(ipsae_pae_cutoff),
    }
    # Per target chain, in flag order. ipSAE takes the RAW pae, ipae_* pae_sym.
    offset = 0
    ipsae_pair_max: list[float] = []  # reference convention: 0.0 for "no pair"
    any_ipsae_pair = False
    for chain in target_chains:
        key = chain_key(chain)
        n_here = len(target_tokens[chain])
        in_chain = (target_cols >= offset) & (target_cols < offset + n_here)
        here = [p for p, keep in zip(pair_pae, in_chain.tolist()) if keep]
        metrics[f"iptm_{key}"] = per_target_iptm[target_chains.index(chain)]
        metrics[f"ipae_{key}_mean"] = float(np.mean(here)) if here else None
        offset += n_here

        forward = ipsae_asym(pae, binder_tokens, target_tokens[chain], ipsae_pae_cutoff)
        reverse = ipsae_asym(pae, target_tokens[chain], binder_tokens, ipsae_pae_cutoff)
        available = [v for v in (forward["ipsae"], reverse["ipsae"]) if v is not None]
        if available:
            any_ipsae_pair = True
            pair_max = max(available)
            # The winning direction's n0res/d0res (binder->target wins a tie).
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
        # Real counts (reference nres1/nres2): 0 means no confident binder residue.
        metrics[f"ipsae_nres_binder_{key}"] = forward["nres_align"]
        metrics[f"ipsae_nres_{key}"] = forward["nres_scored"]

    # Binder-vs-target pairs only; no-pair chains count 0.0, NA if no chain has one.
    metrics["ipsae_binder_mean"] = (
        float(np.mean(ipsae_pair_max)) if any_ipsae_pair else None
    )
    metrics["ipsae_binder_max"] = (
        float(max(ipsae_pair_max)) if any_ipsae_pair else None
    )
    metrics["ipsae_binder_top2"] = (
        float(np.mean(sorted(ipsae_pair_max, reverse=True)[:2]))
        if any_ipsae_pair
        else None
    )
    return metrics


def _self_pae(pae_sym: np.ndarray, tokens: list[int]) -> float | None:
    """Mean PAE inside the binder, diagonal excluded. None for a 1-residue chain."""
    if len(tokens) < 2:
        return None
    block = pae_sym[np.ix_(tokens, tokens)]
    off = ~np.eye(len(tokens), dtype=bool)
    return float(block[off].mean())


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
        description="Binder-specific Boltz confidence and interface PAE."
    )
    ap.add_argument("--task-file", type=Path, required=True,
                    help="Tab-separated 'name<TAB>structure[<TAB>source table]' rows.")
    ap.add_argument("--binder-chains", required=True,
                    help="The binder chain. Exactly one (see the module docstring).")
    ap.add_argument("--target-chains", required=True,
                    help="Comma-joined target chain IDs, in column order.")
    ap.add_argument("--contact-cutoff", type=float, default=5.0,
                    help="Heavy-atom distance defining an interface residue pair "
                         "(the ipae_* columns only; ipSAE never uses a distance).")
    ap.add_argument("--ipsae-pae-cutoff", type=float, default=15.0,
                    help="PAE cutoff in A for the ipsae_* columns (default 15.0: "
                         "the upper end of the 10-15 A band the ipSAE paper "
                         "recommends, and the value the paper's own worked "
                         "examples and command lines use). Two ipSAEs at "
                         "different cutoffs are not comparable.")
    ap.add_argument("--input-column", default="the input structure column",
                    help="Only used to name the column in the 'could not resolve' error.")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(namespace=Args())

    binder = split_list(args.binder_chains)
    targets = split_list(args.target_chains)
    if len(binder) != 1:
        raise ValueError(
            f"--binder-chains must name exactly one chain, got {len(binder)} "
            f"({args.binder_chains!r}): chains_ptm and pair_chains_iptm are per "
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
        try:
            if not src.strip():
                raise ValueError(
                    f"could not resolve a structure path from column "
                    f"{args.input_column!r} for this design -- neither on its own row "
                    f"nor on any ancestor (parent_name/parent_table)"
                )
            metrics = measure(
                Path(src), binder[0], targets, args.contact_cutoff,
                args.ipsae_pae_cutoff, name=name, source=source,
            )
            write_result(args.out_dir, name, "OK", src, metrics, columns)
            print(
                f"{name}: binder_plddt={metrics['binder_plddt']:.1f} "
                f"iptm_binder_max={metrics['iptm_binder_max']:.3f} "
                f"n_if_pairs={metrics['n_if_pairs']} "
                f"ipae_mean="
                + ("NA" if metrics["ipae_mean"] is None
                   else f"{metrics['ipae_mean']:.2f}")
                + " ipsae_binder_max="
                + ("NA" if metrics["ipsae_binder_max"] is None
                   else f"{metrics['ipsae_binder_max']:.3f}")
            )
        except Exception as e:  # noqa: BLE001 - errors are recorded as data
            msg = " ".join(str(e).split())  # one line: it lands in a TSV cell
            print(f"{name}: ERROR {msg}", file=sys.stderr)
            # The chain echo survives an error row; numeric columns are NA, never 0.
            write_result(
                args.out_dir, name, f"error: {msg}", "",
                {"chains": f"{binder[0]}|{','.join(targets)}",
                 "input_source": source}, columns,
            )


if __name__ == "__main__":
    main()
