#!/usr/bin/env python3
"""
Place a C2 dock back into the binder/target frame, and measure whether the partner
protomer occludes the target binding site and whether the two protomers can be linked.

This is the per-array-task step of the dimerfit tool (and works standalone). The
premise, in one sentence: *a binder we already know binds the target was docked
against itself with C2 symmetry, and we now want to know whether the resulting dimer
would sterically block the target (the "off" state of a pH switch) and whether the
C-terminus of one protomer can reach the N-terminus of the other.*

What it does, per design:

1. Read the C2 dock -- **exactly two protomers of the same binder**.
2. Kabsch-superpose dock protomer **A** onto the reference complex's **binder** chain,
   on paired CA atoms.
3. Apply that **one** transform to the **whole dimer** (both protomers), so protomer A
   lands where the reference binder sat and protomer B lands wherever the C2 operator
   put it relative to the target.
4. Measure, in that frame:
   * **occlusion** -- how much of the target's binding site protomer B covers
     (``occluded_frac``, ``n_clash``, ``clash_frac``): we WANT this high.
   * **linkability** -- CA(C-term of A) to CA(N-term of B) (``link_dist``), against a
     ~20 aa linker budget. Recorded, never gated on.
   * **where the C2 interface sits** relative to the epitope (``dimer_iface_res``,
     ``epitope_com_dist``, ``n_overlap_res``, ``frac_overlap``): we want the dimer
     interface ADJACENT to the epitope, so overlap LOW.
5. Write the transformed dimer (no target) and the transformed dimer plus target.

Scope limits -- where these numbers mean nothing
------------------------------------------------
* **Exactly two protomers.** A monomer, or a C3+ assembly, is an error for that row:
  ``link_dist`` and "the partner" are undefined with any other count.
* **The dock must be the same protein as the reference binder**, numbered
  compatibly. The tool does not do sequence alignment: residues are paired ordinally
  (i-th CA to i-th CA) when the two chains hold the same number of CA atoms and by
  residue number otherwise. ``seq_match_frac`` and ``resnum_offset`` are reported so a
  wrong pairing is visible instead of producing plausible geometry. **Backbone-only
  docks are usable for the geometry but the clash/occlusion counts will be too low**
  (no side chains), so dump the dock with real coordinates.
* **The epitope lists are taken as given**, in the reference's numbering. This tool
  re-derives the dock<->reference residue mapping and applies it to the binder-side
  epitope; it does not re-detect any interface. If the epitope column was measured on
  a different structure, every epitope-derived number is wrong and nothing here can
  tell.
* **``n_clash``/``clash_frac`` are heavy-atom distance counts, not an energy.** They
  are not ``fa_rep`` and must not be read as one -- there is no Rosetta in this image.
* **Nothing is relaxed or repacked.** Rigid-body geometry of the inputs as given.

Writes a one-row TSV (name, status, path, complex_path, then the metrics) that
collect_dimerfit.py merges back into the table. Errors are recorded as data (a status
starting with ``error:``) rather than only crashing, so partial array runs still
collect and the batch still exits 0.

Normally invoked per array task by the dimerfit tool -- run ``sapia run dimerfit``
rather than calling this directly.

Usage (standalone, one task file of
``name<TAB>dock<TAB>ref<TAB>binder_epitope<TAB>target_epitope`` lines):
    python tools/dimerfit/dimerfit_worker.py \\
        --task-file task_0.tsv --out-dir out/ \\
        --dock-chains auto --binder-chain-in-ref B --target-chains-in-ref A
"""

import argparse
import csv
import re
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np

# A C2 dock has exactly two protomers; anything else makes "the partner",
# link_dist and the occlusion question undefined.
N_PROTOMERS = 2
# Below this fraction of identity-agreeing matched pairs the dock<->reference
# pairing is suspect (wrong chain, wrong design, a register shift): warn and
# downgrade the status, but still report every metric.
SEQ_MATCH_MIN = 0.9
# How dock residues are paired with reference-binder residues.
MATCH_MODES = ("auto", "ordinal", "resnum")
# Chain IDs preferred for the target in the inspection complex when its own ID
# collides with a dock protomer's (the dock keeps its IDs: that file is the one
# downstream tools consume).
PREFERRED_COMPLEX_IDS = "TUVWXYZ"

# Metric columns of the per-design TSV, in order (bare names: collect_dimerfit.py
# hands them to the driver, which leaf-prefixes them to dimerfit_<name>).
METRIC_COLUMNS = [
    # --- mapping trust: a wrong superposition gives plausible geometry ----------
    "align_rmsd",
    "n_align_atoms",
    "seq_match_frac",
    "resnum_offset",
    "resnum_match",
    "dock_chains",
    "ref_binder_chain",
    "ref_target_chains",
    "complex_target_chains",
    "n_protomers",
    "n_res_a",
    "n_res_b",
    "ref_binder_len",
    "n_atoms_b",
    # --- C2 interface geometry, relative to the epitope ------------------------
    "dimer_iface_res",
    "n_dimer_iface_res",
    "dimer_iface_com",
    "epitope_com",
    "epitope_com_dist",
    "n_epitope_res",
    "overlap_res",
    "n_overlap_res",
    "frac_overlap",
    "link_dist",
    "cterm_to_nterm_res",
    # --- occlusion of the target binding site by protomer B -------------------
    "n_clash",
    "clash_frac",
    "min_dist_b_target",
    "n_target_epitope_res",
    "occluded_res",
    "n_occluded_res",
    "occluded_frac",
    "seconds",
]
RESULT_COLUMNS = ["name", "status", "path", "complex_path", *METRIC_COLUMNS]


# ---------------------------------------------------------------------------
# Residue bookkeeping
# ---------------------------------------------------------------------------


@dataclass
class Res:
    """One amino acid, reduced to what the geometry needs."""

    label: str  # "B:12" (chain:resnum[icode]) -- the format every tool here parses
    chain: str
    seqid: int
    icode: str
    resname: str
    code: str  # one-letter, 'X' when unknown
    heavy: np.ndarray  # (M, 3) heavy-atom coordinates
    masses: np.ndarray  # (M,) atomic weights, for the centre of mass
    ca: np.ndarray | None

    @property
    def key(self) -> tuple[int, str]:
        """Residue identity within its chain: (resnum, insertion code)."""
        return (self.seqid, self.icode)


def _is_heavy(atom: gemmi.Atom) -> bool:
    """Non-hydrogen. Tested on the element name rather than ``Element.is_hydrogen``,
    which is a method in some gemmi versions and a property in others (0.7.x)."""
    return atom.element.name.upper() not in ("H", "D")


def _is_amino_acid(res: gemmi.Residue) -> bool:
    info = gemmi.find_tabulated_residue(res.name)
    return bool(info and info.is_amino_acid())


def one_letter(res: gemmi.Residue) -> str:
    """One-letter code of an amino acid, 'X' when gemmi does not know it."""
    info = gemmi.find_tabulated_residue(res.name)
    code = (info.one_letter_code if info else "").upper()
    return code if code.isalpha() else "X"


def _residue(chain_name: str, res: gemmi.Residue) -> Res | None:
    """One gemmi residue -> ``Res``, or None when it has no heavy atoms."""
    coords: list[list[float]] = []
    masses: list[float] = []
    ca: np.ndarray | None = None
    for atom in res:
        if not _is_heavy(atom):
            continue
        pos = [atom.pos.x, atom.pos.y, atom.pos.z]
        coords.append(pos)
        masses.append(float(atom.element.weight))
        if atom.name == "CA" and ca is None:
            ca = np.asarray(pos, dtype=float)
    if not coords:
        return None
    icode = res.seqid.icode.strip()
    return Res(
        label=f"{chain_name}:{res.seqid.num}{icode}",
        chain=chain_name,
        seqid=res.seqid.num,
        icode=icode,
        resname=res.name.upper(),
        code=one_letter(res),
        heavy=np.asarray(coords, dtype=float),
        masses=np.asarray(masses, dtype=float),
        ca=ca,
    )


def chain_residues(model: gemmi.Model, chain_names: list[str]) -> list[Res]:
    """Amino acids of the named chains, in the order the chains were named."""
    out: list[Res] = []
    for name in chain_names:
        for chain in model:
            if chain.name != name:
                continue
            for res in chain:
                if _is_amino_acid(res) and (r := _residue(chain.name, res)) is not None:
                    out.append(r)
    return out


def protein_chain_names(model: gemmi.Model) -> list[str]:
    """Names of every chain holding at least one amino acid, in FILE ORDER.

    File order is what ``--dock-chains auto`` means: first = protomer A, second =
    protomer B.
    """
    seen: list[str] = []
    for chain in model:
        if chain.name not in seen and any(_is_amino_acid(r) for r in chain):
            seen.append(chain.name)
    return seen


def stack(residues: list[Res]) -> np.ndarray:
    """All heavy-atom coordinates of a residue list as one (N, 3) array."""
    if not residues:
        return np.empty((0, 3))
    return np.concatenate([r.heavy for r in residues])


def centre_of_mass(residues: list[Res]) -> np.ndarray | None:
    """Mass-weighted centre of the residues' heavy atoms."""
    if not residues:
        return None
    coords = stack(residues)
    masses = np.concatenate([r.masses for r in residues])
    total = float(masses.sum())
    if total <= 0:
        return None
    return (coords * masses[:, None]).sum(axis=0) / total


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def kabsch(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Optimal rigid superposition of ``p`` onto ``q`` (both (N, 3), paired).

    Returns ``(rot, trans, rmsd)`` with ``p @ rot.T + trans`` the moved ``p``.
    Reflections are excluded the usual way (flip the last singular vector when the
    naive rotation has determinant -1).
    """
    if p.shape != q.shape or p.shape[0] < 3:
        raise ValueError(f"need >=3 paired points for a superposition, got {p.shape[0]}")
    p_mean, q_mean = p.mean(axis=0), q.mean(axis=0)
    pc, qc = p - p_mean, q - q_mean
    u, _, vt = np.linalg.svd(pc.T @ qc)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    trans = q_mean - rot @ p_mean
    rmsd = float(np.sqrt((((pc @ rot.T) - qc) ** 2).sum(axis=1).mean()))
    return rot, trans, rmsd


def pair_stats(
    a: np.ndarray, b: np.ndarray, cutoff: float, chunk: int = 256
) -> tuple[float | None, np.ndarray]:
    """``(min distance a<->b, mask of a-atoms within cutoff of any b-atom)``.

    The minimum is None when either set is empty. Chunked over ``a`` so the pairwise
    block stays small.
    """
    mask = np.zeros(a.shape[0], dtype=bool)
    if a.size == 0 or b.size == 0:
        return None, mask
    best = np.inf
    for start in range(0, a.shape[0], chunk):
        block = a[start : start + chunk]
        delta = block[:, None, :] - b[None, :, :]
        dist = np.sqrt(np.einsum("ijk,ijk->ij", delta, delta))
        nearest = dist.min(axis=1)
        best = min(best, float(nearest.min()))
        mask[start : start + chunk] = nearest < cutoff
    return best, mask


def min_dist_to(residue: Res, coords: np.ndarray) -> float | None:
    """Smallest heavy-atom distance from one residue to a point set."""
    if coords.size == 0 or residue.heavy.size == 0:
        return None
    delta = residue.heavy[:, None, :] - coords[None, :, :]
    return float(np.sqrt(np.einsum("ijk,ijk->ij", delta, delta)).min())


# ---------------------------------------------------------------------------
# Residue matching (the trust layer)
# ---------------------------------------------------------------------------


def match_residues(
    dock: list[Res], ref: list[Res], mode: str
) -> tuple[list[tuple[int, int]], str]:
    """Pair dock residues with reference residues -> ``(index pairs, mode used)``.

    ``ordinal`` pairs the i-th with the i-th (right for a generator that renumbered
    1..N, which is what rpxdock does to its dump); ``resnum`` pairs equal residue
    numbers (tolerant of gaps); ``auto`` takes ordinal when the two chains hold the
    same number of CA atoms and resnum otherwise. Either way the pairing is CHECKED
    by residue identity and reported as ``seq_match_frac``.
    """
    if mode not in MATCH_MODES:
        raise ValueError(f"--resnum-match must be one of {MATCH_MODES}, got {mode!r}")
    if mode == "auto":
        mode = "ordinal" if len(dock) == len(ref) and len(dock) else "resnum"

    if mode == "ordinal":
        if len(dock) != len(ref):
            raise ValueError(
                f"ordinal matching needs an equal number of CA atoms, got {len(dock)} "
                f"(dock protomer A) vs {len(ref)} (reference binder); use "
                f"--resnum-match auto or resnum"
            )
        return [(i, i) for i in range(len(dock))], mode

    d_index = {r.key: i for i, r in enumerate(dock)}
    r_index = {r.key: i for i, r in enumerate(ref)}
    shared = sorted(set(d_index) & set(r_index))
    return [(d_index[k], r_index[k]) for k in shared], mode


def parse_res_labels(value: str, what: str) -> list[tuple[str, int, str]]:
    """``'B:12,B:15A'`` -> ``[('B', 12, ''), ('B', 15, 'A')]``.

    The ``chain:resnum[icode]`` form is what ifacegeom writes and what every tool in
    this campaign parses. The bare ``B12`` form (ringfit --hotspots) is rejected
    loudly rather than guessed at: the two conventions are not interchangeable.
    """
    out: list[tuple[str, int, str]] = []
    for token in (t.strip() for t in value.split(",")):
        if not token:
            continue
        m = re.fullmatch(r"([A-Za-z0-9]{1,4}):(-?\d+)([A-Za-z]?)", token)
        if not m:
            raise ValueError(
                f"cannot parse {what} residue {token!r}: expected 'chain:resnum' "
                f"(e.g. 'B:12'), not 'B12'"
            )
        out.append((m.group(1), int(m.group(2)), m.group(3).strip()))
    if not out:
        raise ValueError(f"{what} is empty: nothing to measure against")
    return out


# ---------------------------------------------------------------------------
# Chain resolution
# ---------------------------------------------------------------------------


def split_list(value: str) -> list[str]:
    """Comma-joined CLI list -> items, dropping blanks."""
    return [v.strip() for v in value.split(",") if v.strip()]


def resolve_dock_chains(model: gemmi.Model, spec: str, src: str) -> list[str]:
    """``--dock-chains`` -> exactly two chain IDs, [protomer A, protomer B].

    The dock file must hold **exactly two** protein chains, whatever ``--dock-chains``
    says: with one protomer there is no partner to occlude anything and no N-terminus
    to link to, and with three the C2 premise is simply false. ``auto`` takes the
    file's own chain order; an explicit list only decides which of the two is
    protomer A.
    """
    present = protein_chain_names(model)
    if len(present) != N_PROTOMERS:
        raise ValueError(
            f"{src} holds {len(present)} protein chain(s) "
            f"({','.join(present) or 'none'}), expected exactly {N_PROTOMERS} "
            f"protomers of a C2 dock: n_protomers={len(present)}"
        )
    if spec.strip().lower() == "auto":
        return present
    chains = split_list(spec)
    if len(chains) != N_PROTOMERS:
        raise ValueError(
            f"--dock-chains must name exactly {N_PROTOMERS} chains "
            f"(protomer A then protomer B), got {chains or 'none'}"
        )
    if chains[0] == chains[1]:
        raise ValueError(f"--dock-chains names the same chain twice: {chains}")
    if missing := [c for c in chains if c not in present]:
        raise ValueError(
            f"--dock-chains {missing} absent from {src} "
            f"(protein chains present: {','.join(present) or 'none'})"
        )
    return chains


def resolve_ref_chains(
    model: gemmi.Model, binder_spec: str, target_spec: str, src: str
) -> tuple[str, list[str]]:
    """``--binder-chain-in-ref``/``--target-chains-in-ref`` -> (binder, targets).

    A chain named but absent is an error for this design, never a smaller selection:
    a silently narrower target would understate every occlusion number, i.e. fail in
    the direction that looks like a passing design.
    """
    present = protein_chain_names(model)
    binder = binder_spec.strip()
    if not binder:
        raise ValueError("--binder-chain-in-ref must name one chain")
    if binder not in present:
        raise ValueError(
            f"--binder-chain-in-ref {binder!r} absent from {src} "
            f"(protein chains present: {','.join(present) or 'none'}). Note the "
            f"inherited hEGFR complexes have the BINDER on chain B and the TARGET on "
            f"chain A -- the opposite of the usual convention."
        )
    if target_spec.strip().lower() == "auto":
        targets = [c for c in present if c != binder]
        if not targets:
            raise ValueError(
                f"--target-chains-in-ref auto found no protein chain outside the "
                f"binder chain {binder!r} in {src}; this file holds no target, so "
                f"there is no binding site to occlude"
            )
    else:
        targets = split_list(target_spec)
        if not targets:
            raise ValueError(
                "--target-chains-in-ref must name at least one chain, or 'auto'"
            )
        if missing := [c for c in targets if c not in present]:
            raise ValueError(
                f"--target-chains-in-ref {missing} absent from {src} "
                f"(protein chains present: {','.join(present) or 'none'})"
            )
    if binder in targets:
        raise ValueError(
            f"chain {binder!r} is both the reference binder and a reference target "
            f"chain; the two sides must be disjoint"
        )
    return binder, targets


def complex_chain_ids(dock_ids: list[str], target_ids: list[str]) -> list[str]:
    """Chain IDs the reference target chains get in the inspection complex.

    The dock keeps its own IDs -- that file is the one downstream tools consume, so
    its chain names must not move. A target chain whose ID collides with a protomer's
    is renamed (preferring T, U, V ...), and the result is reported as
    ``complex_target_chains`` so the inspection file stays readable by hand.
    """
    used = set(dock_ids)
    out: list[str] = []
    for tid in target_ids:
        if tid not in used:
            used.add(tid)
            out.append(tid)
            continue
        candidates = PREFERRED_COMPLEX_IDS + "ABCDEFGHIJKLMNOPQRS" + "0123456789"
        for cand in candidates:
            if cand not in used:
                used.add(cand)
                out.append(cand)
                break
        else:  # pragma: no cover - 36 free IDs is far beyond any real input
            raise ValueError("no free chain ID left for the inspection complex")
    return out


# ---------------------------------------------------------------------------
# Per-design computation
# ---------------------------------------------------------------------------


def load_structure(path: Path) -> gemmi.Structure:
    """Read a PDB or mmCIF file into a single-model, single-conformer structure."""
    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    if len(st) == 0:
        raise ValueError(f"no model in {path}")
    return st


def _fmt(value: object) -> str:
    """TSV cell: NA (None) becomes empty, floats get 3 decimals."""
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def _fmt_point(point: np.ndarray | None) -> str | None:
    """A coordinate as ``x,y,z`` with 3 decimals (NA when absent)."""
    if point is None:
        return None
    return ",".join(f"{v:.3f}" for v in point)


@dataclass
class Design:
    """One sub-manifest line."""

    name: str
    dock: Path
    ref: Path
    epitope: str  # binder-side, in the REFERENCE binder's numbering
    target_epitope: str  # target-side, in the reference target's numbering


def build_complex(
    dock_model: gemmi.Model,
    dock_chains: list[str],
    ref_model: gemmi.Model,
    target_chains: list[str],
    new_target_ids: list[str],
    name: str,
) -> gemmi.Structure:
    """The transformed dimer plus the reference target, for inspection only."""
    out = gemmi.Structure()
    out.name = name
    model = gemmi.Model("1")
    for cid in dock_chains:
        for chain in dock_model:
            if chain.name == cid:
                model.add_chain(chain)
                break
    for cid, new_id in zip(target_chains, new_target_ids):
        for chain in ref_model:
            if chain.name != cid:
                continue
            new = gemmi.Chain(new_id)
            # Everything on the target chain, ligands included -- this file exists to
            # be looked at, so nothing is quietly dropped.
            for res in chain:
                new.add_residue(res)
            model.add_chain(new)
            break
    out.add_model(model)
    out.setup_entities()
    return out


def compute(
    args: "Args", design: Design
) -> tuple[str, dict[str, object], gemmi.Structure, gemmi.Structure]:
    """Every dimerfit metric for one dock, plus the two structures to write.

    Raises on anything unusable (the caller records the message as data).
    """
    started = time.time()

    dock_st = load_structure(design.dock)
    ref_st = load_structure(design.ref)
    ref_st.remove_waters()
    dock_model, ref_model = dock_st[0], ref_st[0]

    dock_chains = resolve_dock_chains(dock_model, args.dock_chains, design.dock.name)
    ref_binder_chain, ref_target_chains = resolve_ref_chains(
        ref_model,
        args.binder_chain_in_ref,
        args.target_chains_in_ref,
        design.ref.name,
    )

    prot_a = chain_residues(dock_model, dock_chains[:1])
    prot_b = chain_residues(dock_model, dock_chains[1:2])
    ref_binder = chain_residues(ref_model, [ref_binder_chain])
    ref_targets = chain_residues(ref_model, ref_target_chains)
    if not prot_a or not prot_b:
        raise ValueError(
            f"a protomer has no amino acids (A {len(prot_a)}, B {len(prot_b)})"
        )
    if not ref_binder:
        raise ValueError(f"reference binder chain {ref_binder_chain} has no amino acids")
    if not ref_targets:
        raise ValueError(
            f"reference target chain(s) {','.join(ref_target_chains)} have no amino "
            f"acids; there is no binding site to occlude"
        )

    # --- 1. the mapping, and its evidence ----------------------------------
    a_ca = [r for r in prot_a if r.ca is not None]
    ref_ca = [r for r in ref_binder if r.ca is not None]
    pairs, used_mode = match_residues(a_ca, ref_ca, args.resnum_match)
    if not pairs:
        raise ValueError(
            f"no residue matched between dock protomer {dock_chains[0]} "
            f"({len(a_ca)} CA) and reference binder {ref_binder_chain} "
            f"({len(ref_ca)} CA) under {used_mode} matching"
        )
    n_same = sum(a_ca[i].code == ref_ca[j].code for i, j in pairs)
    seq_match_frac = n_same / len(pairs)
    offsets = [ref_ca[j].seqid - a_ca[i].seqid for i, j in pairs]
    resnum_offset = Counter(offsets).most_common(1)[0][0]
    # Reference residue -> the dock protomer-A residue it corresponds to. This, not
    # an assumption about numbering, is how the epitope reaches the dock.
    ref_to_dock: dict[tuple[int, str], Res] = {
        ref_ca[j].key: a_ca[i] for i, j in pairs
    }
    if seq_match_frac < SEQ_MATCH_MIN:
        print(
            f"WARNING: {design.name}: only {seq_match_frac:.2f} of the {len(pairs)} "
            f"matched residues agree in identity ({used_mode} matching) between dock "
            f"chain {dock_chains[0]} and reference binder chain {ref_binder_chain}. "
            f"The superposition -- and every number under it -- is probably wrong: "
            f"check --binder-chain-in-ref and --dock-chains.",
            file=sys.stderr,
        )

    # --- 2. one transform, applied to the WHOLE dimer -----------------------
    rot, trans, align_rmsd = kabsch(
        np.asarray([a_ca[i].ca for i, _ in pairs]),
        np.asarray([ref_ca[j].ca for _, j in pairs]),
    )
    transform = gemmi.Transform()
    transform.mat.fromlist(rot.tolist())
    transform.vec.fromlist(trans.tolist())
    # transform_pos_and_adp on the MODEL moves every chain, which is the whole point:
    # protomer A lands on the reference binder and protomer B goes wherever the C2
    # operator put it relative to the target.
    dock_model.transform_pos_and_adp(transform)

    # Re-read the moved coordinates (the Res objects above hold pre-transform copies).
    prot_a = chain_residues(dock_model, dock_chains[:1])
    prot_b = chain_residues(dock_model, dock_chains[1:2])
    a_ca = [r for r in prot_a if r.ca is not None]
    ref_to_dock = {ref_ca[j].key: a_ca[i] for i, j in pairs}

    a_coords, b_coords = stack(prot_a), stack(prot_b)
    target_coords = stack(ref_targets)

    # --- 3. the C2 interface, on the protomer-A side ------------------------
    dimer_iface = [
        r
        for r in prot_a
        if (d := min_dist_to(r, b_coords)) is not None
        and d <= args.dimer_contact_cutoff
    ]
    dimer_iface_keys = {r.key for r in dimer_iface}

    # --- 4. the epitope, mapped onto protomer A ----------------------------
    epitope_labels = parse_res_labels(design.epitope, "--epitope-column")
    if wrong := sorted({c for c, _, _ in epitope_labels if c != ref_binder_chain}):
        raise ValueError(
            f"--epitope-column lists chain(s) {wrong} but --binder-chain-in-ref is "
            f"{ref_binder_chain!r}: the epitope was measured on a different selection, "
            f"so mapping it onto the dock would be meaningless"
        )
    epitope: list[Res] = []
    unmapped: list[str] = []
    for chain, num, icode in epitope_labels:
        res = ref_to_dock.get((num, icode))
        if res is None:
            unmapped.append(f"{chain}:{num}{icode}")
        else:
            epitope.append(res)
    if unmapped:
        raise ValueError(
            f"epitope residue(s) {','.join(unmapped)} are not among the "
            f"{len(pairs)} matched dock/reference residues, so the epitope cannot be "
            f"placed on the dock (resnum_offset {resnum_offset}, "
            f"seq_match_frac {seq_match_frac:.2f})"
        )

    overlap = [r for r in epitope if r.key in dimer_iface_keys]
    dimer_iface_com = centre_of_mass(dimer_iface)
    epitope_com = centre_of_mass(epitope)
    epitope_com_dist = (
        float(np.linalg.norm(dimer_iface_com - epitope_com))
        if dimer_iface_com is not None and epitope_com is not None
        else None
    )

    # --- 5. linkability: C-term of A -> N-term of B ------------------------
    cterm_a, nterm_b = prot_a[-1], prot_b[0]
    link_dist = (
        float(np.linalg.norm(cterm_a.ca - nterm_b.ca))
        if cterm_a.ca is not None and nterm_b.ca is not None
        else None
    )

    # --- 6. occlusion of the target binding site by protomer B -------------
    min_dist_b_target, clash_mask = pair_stats(
        b_coords, target_coords, args.clash_cutoff
    )
    n_clash = int(clash_mask.sum())
    clash_frac = n_clash / len(b_coords) if len(b_coords) else None

    target_epitope_labels = parse_res_labels(
        design.target_epitope, "--target-epitope-column"
    )
    by_key = {(r.chain, r.seqid, r.icode): r for r in ref_targets}
    target_epitope: list[Res] = []
    missing_target: list[str] = []
    for chain, num, icode in target_epitope_labels:
        res = by_key.get((chain, num, icode))
        if res is None:
            missing_target.append(f"{chain}:{num}{icode}")
        else:
            target_epitope.append(res)
    if missing_target:
        raise ValueError(
            f"--target-epitope-column residue(s) {','.join(missing_target)} are not "
            f"in the reference target chain(s) {','.join(ref_target_chains)}: the "
            f"epitope column and this reference do not describe the same structure"
        )
    occluded = [
        r
        for r in target_epitope
        if (d := min_dist_to(r, b_coords)) is not None and d <= args.occlusion_cutoff
    ]
    occluded_frac = len(occluded) / len(target_epitope) if target_epitope else None

    metrics: dict[str, object] = {
        "align_rmsd": align_rmsd,
        "n_align_atoms": len(pairs),
        "seq_match_frac": seq_match_frac,
        "resnum_offset": resnum_offset,
        "resnum_match": used_mode,
        "dock_chains": ",".join(dock_chains),
        "ref_binder_chain": ref_binder_chain,
        "ref_target_chains": ",".join(ref_target_chains),
        # Counted, not assumed -- resolve_dock_chains already refused anything but 2,
        # so this column is the evidence for that check rather than a restatement.
        "n_protomers": len(protein_chain_names(dock_model)),
        "n_res_a": len(prot_a),
        "n_res_b": len(prot_b),
        "ref_binder_len": len(ref_binder),
        "n_atoms_b": len(b_coords),
        "dimer_iface_res": ",".join(r.label for r in dimer_iface),
        "n_dimer_iface_res": len(dimer_iface),
        "dimer_iface_com": _fmt_point(dimer_iface_com),
        "epitope_com": _fmt_point(epitope_com),
        "epitope_com_dist": epitope_com_dist,
        "n_epitope_res": len(epitope),
        "overlap_res": ",".join(r.label for r in overlap),
        "n_overlap_res": len(overlap),
        "frac_overlap": len(overlap) / len(epitope) if epitope else None,
        "link_dist": link_dist,
        "cterm_to_nterm_res": f"{cterm_a.label}->{nterm_b.label}",
        "n_clash": n_clash,
        "clash_frac": clash_frac,
        "min_dist_b_target": min_dist_b_target,
        "n_target_epitope_res": len(target_epitope),
        "occluded_res": ",".join(r.label for r in occluded),
        "n_occluded_res": len(occluded),
        "occluded_frac": occluded_frac,
        "seconds": time.time() - started,
    }

    new_target_ids = complex_chain_ids(dock_chains, ref_target_chains)
    metrics["complex_target_chains"] = ",".join(new_target_ids)
    complex_st = build_complex(
        dock_model,
        dock_chains,
        ref_model,
        ref_target_chains,
        new_target_ids,
        design.name,
    )

    # A suspect mapping, and a numbering that is not the identity map this campaign
    # expects, are reported rather than swallowed: every metric is still written, but
    # the status keeps the row out of any `== "OK"` selection.
    if seq_match_frac < SEQ_MATCH_MIN:
        status = f"warn: seq_match_frac {seq_match_frac:.2f}"
    elif resnum_offset != 0:
        status = f"warn: resnum_offset {resnum_offset}"
    elif not dimer_iface:
        status = "warn: no dimer interface residues"
    else:
        status = "OK"
    return status, metrics, dock_st, complex_st


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def write_result(
    result_tsv: Path,
    name: str,
    status: str,
    path: str,
    complex_path: str,
    metrics: dict,
) -> None:
    result_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(result_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(RESULT_COLUMNS)
        writer.writerow(
            [name, status, path, complex_path]
            + [_fmt(metrics.get(c)) for c in METRIC_COLUMNS]
        )


def read_task_file(path: Path, expect: int | None) -> list[Design]:
    """``name<TAB>dock<TAB>ref<TAB>epitope<TAB>target_epitope`` lines -> designs.

    ``expect`` is the design count the manifest said this task carries. A mismatch is
    FATAL: a half-visible sub-manifest on a Volume reads as a short task rather than
    as an error, and that is exactly how 18 of 37 rpxdock tasks once ran with empty
    inputs (see the guard in dimerfit.sh).
    """
    if not path.is_file():
        raise FileNotFoundError(f"task file missing: {path}")
    designs: list[Design] = []
    for lineno, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 5:
            raise ValueError(
                f"malformed task line {lineno} in {path}: expected 5 tab-separated "
                f"fields (name, dock, ref, epitope, target_epitope), got "
                f"{len(fields)}: {line!r}"
            )
        name, dock, ref, epitope, target_epitope = fields
        if not all((name, dock, ref, epitope, target_epitope)):
            raise ValueError(
                f"task line {lineno} in {path} has an empty field: {line!r}"
            )
        designs.append(Design(name, Path(dock), Path(ref), epitope, target_epitope))
    if expect is not None and len(designs) != expect:
        raise ValueError(
            f"{path} holds {len(designs)} design(s) but the manifest says this task "
            f"carries {expect}. Refusing to run on a partially visible task file."
        )
    return designs


class Args(argparse.Namespace):
    task_file: Path
    n_designs: int
    dock_chains: str
    binder_chain_in_ref: str
    target_chains_in_ref: str
    resnum_match: str
    clash_cutoff: float
    occlusion_cutoff: float
    dimer_contact_cutoff: float
    out_dir: Path


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Place C2 docks back into the binder/target frame and measure "
        "occlusion and linkability for a batch of designs."
    )
    ap.add_argument(
        "--task-file",
        type=Path,
        required=True,
        help="Sub-manifest: one 'name<TAB>dock<TAB>ref<TAB>epitope<TAB>"
        "target_epitope' line per design.",
    )
    ap.add_argument(
        "--n-designs",
        type=int,
        default=-1,
        help="Design count the manifest says this task carries; -1 skips the check. "
        "A mismatch is fatal (a half-visible task file would otherwise run short).",
    )
    ap.add_argument("--dock-chains", default="auto")
    ap.add_argument("--binder-chain-in-ref", default="B")
    ap.add_argument("--target-chains-in-ref", default="A")
    ap.add_argument("--resnum-match", choices=MATCH_MODES, default="auto")
    ap.add_argument("--clash-cutoff", type=float, default=2.5)
    ap.add_argument("--occlusion-cutoff", type=float, default=5.0)
    ap.add_argument("--dimer-contact-cutoff", type=float, default=5.0)
    ap.add_argument(
        "--out-dir", type=Path, required=True, help="Where the per-design files go."
    )
    args = ap.parse_args(namespace=Args())

    args.out_dir.mkdir(parents=True, exist_ok=True)
    designs = read_task_file(
        args.task_file, None if args.n_designs < 0 else args.n_designs
    )

    failures = 0
    for design in designs:
        result_tsv = args.out_dir / f"{design.name}.tsv"
        dimer_pdb = args.out_dir / f"{design.name}_dimer.pdb"
        complex_pdb = args.out_dir / f"{design.name}_complex.pdb"
        try:
            status, metrics, dimer_st, complex_st = compute(args, design)
            dimer_st.setup_entities()
            dimer_st.write_pdb(str(dimer_pdb))
            complex_st.write_pdb(str(complex_pdb))
            print(
                f"{design.name}: {status} -- align_rmsd "
                f"{_fmt(metrics['align_rmsd'])} A, seq_match_frac "
                f"{_fmt(metrics['seq_match_frac'])}, resnum_offset "
                f"{metrics['resnum_offset']}, link_dist {_fmt(metrics['link_dist'])} "
                f"A, overlap {metrics['n_overlap_res']}/{metrics['n_epitope_res']}, "
                f"occluded {metrics['n_occluded_res']}/"
                f"{metrics['n_target_epitope_res']}, n_clash {metrics['n_clash']}"
            )
            write_result(
                result_tsv,
                design.name,
                status,
                str(dimer_pdb),
                str(complex_pdb),
                metrics,
            )
        except Exception as e:  # noqa: BLE001 - errors are recorded as data
            failures += 1
            print(f"{design.name}: ERROR {e}", file=sys.stderr)
            # No partial files left behind: a half-written structure would be worse
            # than none, since the column would point at something unusable.
            for stale in (dimer_pdb, complex_pdb):
                stale.unlink(missing_ok=True)
            write_result(result_tsv, design.name, f"error: {e}", "", "", {})

    # Every design has a result either way; exiting non-zero would only hide that
    # behind a task-level failure, so the batch succeeds and the table carries the
    # errors.
    print(f"dimerfit: {failures} failed of {len(designs)}")


if __name__ == "__main__":
    main()
