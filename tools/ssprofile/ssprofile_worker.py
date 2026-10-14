#!/usr/bin/env python3
"""
The shape and secondary-structure make-up of one chain selection, per design:
CA radius of gyration plus 3-state DSSP helix / strand / loop fractions.

This is the per-array-task step of the ssprofile tool (and works standalone). One task
scores several designs. For each design it reads ONE structure file, keeps only the
chains named by ``--design-chains``, and describes that selection's OWN backbone:

    radius_of_gyration  CA-based and MASS-UNWEIGHTED: centre the selection's CA
                        coordinates on their plain mean, then
                        Rg = sqrt(mean(|r_i - centroid|^2)), in A. One CA per
                        amino-acid residue; every CA weighs the same. With several
                        chains selected it is the Rg of their UNION.
    helix_fraction      residues DSSP calls H / n_res  (H includes 3-10 'G' and pi 'I')
    beta_fraction       residues DSSP calls E / n_res  (E includes isolated bridges 'B')
    loop_fraction       everything else / n_res        (turn 'T', bend 'S', coil ' ')
    n_res               amino-acid residues with a CA in the selection -- the
                        denominator of the three fractions and the CA count behind Rg
    n_chain_breaks      C(i)-N(i+1) gaps > 2.5 A inside a selected chain
    chains_used         the selected chains, post-validation, comma-joined
    sequence            one-letter sequence of the selection, chains joined by '/'
                        (modified residues lower-case, unknown 'X') -- the wrong-chain
                        guard when binder and target have similar lengths
    input_column        the table column the structure path came from   } provenance,
    source_path         the file actually read (post CIF->PDB staging)  } on EVERY row

THE SS ENGINE, AND WHY IT IS THIS ONE. rpxdock's ``frac_helix``/``frac_sheet``/
``frac_loop`` come from willutil's ``wu.dssp``. Read at willutil HEAD (ded19f6,
``willutil/chem/dssp.py``), that function is NOT a DSSP implementation: it writes the
backbone N/CA/C/O to a temporary PDB as poly-GLY residues and calls
``mdtraj.compute_dssp`` (simplified, 3-state), then renames 'C' to 'L'. This worker
makes the SAME call on the SAME input -- an mdtraj topology of N/CA/C/O atoms, every
residue named GLY, simplified=True -- built in memory instead of through a temporary
file. It does not install willutil itself: willutil is git-only (not on PyPI) and pulls
numba, xarray, netcdf4, seaborn, Bio, cppimport, pytest and more, all for a 15-line
wrapper around mdtraj.

Two consequences of copying that convention exactly, both deliberate:

  * Every residue is GLY to the SS engine, so PROLINE is treated as an H-bond donor,
    which strict DSSP does not do. Measured on 1UBQ, 2GB1, 1CRN and 4HHB (319
    residues, 15 prolines): one residue differs in total (1CRN), against the same
    call with real residue names. Kept for parity with rpxdock's columns.
  * The selection is assigned IN ISOLATION. H-bonds to unselected chains do not
    exist for DSSP here, so a binder strand that pairs only with a target strand
    reads as LOOP, not strand.

FAILURE CONTRACT -- status first, then which columns are trustworthy:

    'OK'        all of Rg, the three fractions, n_res, n_chain_breaks and
                chains_used are measured. The fractions sum to 1.0.
    'warn: ...' Rg, n_res, chains_used and sequence are measured; the three
                fractions are NA; n_chain_breaks is measured on a break-warn and NA
                on a missing-atom warn (the gap test needs N and C). Raised when
                the CA trace is usable but the SS assignment is not: an amino-acid
                residue missing N, CA, C or O (never guessed -- willutil would
                guess a missing O, this tool refuses; a CA-less residue is also
                absent from n_res and Rg), a backbone break inside a selected
                chain (the premise is ONE continuous polypeptide per chain), or
                any failure of the DSSP call itself.
    'error: ...' every measured column NA (input_column and source_path are still
                written). A chain in --design-chains absent from the structure
                (never a fallback to some other chain), ANY selected chain with
                zero CA atoms (empty, HETATM-only -- checked per chain, so 'A,L'
                with a ligand-only L is an error, not a smaller selection), a
                structure that cannot be read or could not be resolved.

A fraction of 0.0 on an OK row is a measurement (this binder has no strand); an empty
cell is the absence of one. Nothing here defaults to zero.

Writes, per design, ``<name>.tsv`` (one row: name, status, path, metrics...) and, on
OK and warn, ``<name>_residues.tsv`` -- one row per residue of the selection (chain,
resnum, resname, ss) in the order Rg and DSSP saw them; ``ss`` is NA on a warn row.

Normally invoked per array task by the ssprofile tool -- run ``sapia run ssprofile``
rather than calling this directly.

Usage (standalone):
    printf 'design_0\\tcomplex.pdb\\n' > task.tsv
    python tools/ssprofile/ssprofile_worker.py --task-file task.tsv \\
        --design-chains B --out-dir out/
"""

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np

# Metric columns of the per-design TSV, in order (bare names: collect_ssprofile.py
# hands them to the driver, which leaf-prefixes them to ssprofile_<name>).
METRIC_COLUMNS = [
    "radius_of_gyration",
    "helix_fraction",
    "beta_fraction",
    "loop_fraction",
    "n_res",
    "n_chain_breaks",
    "chains_used",
    # --- added 2026-10-09 after tool-reviewer: provenance + wrong-chain guard ---
    "sequence",
    "input_column",
    "source_path",
]
RESULT_COLUMNS = ["name", "status", "path", *METRIC_COLUMNS]
RESIDUE_COLUMNS = ["chain", "resnum", "resname", "ss"]

BACKBONE = ("N", "CA", "C", "O")
# DSSP's own kMaxPeptideBondLength: a C(i)-N(i+1) distance above this is a break.
# mdtraj applies the same idea internally (measured: a gapped chain given as one chain
# and the same chain split in two give identical assignments).
MAX_PEPTIDE_BOND = 2.5
# Simplified (3-state) mdtraj codes -> the column each one counts toward.
SS_CODES = {"H": "H", "E": "E", "C": "L"}


def split_list(value: str) -> list[str]:
    """Comma-joined list -> stripped, non-empty tokens."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


@dataclass
class Residue:
    chain: str
    resnum: str  # seqid number + insertion code, e.g. "96" or "96A"
    resname: str
    ca: np.ndarray  # (3,)
    backbone: np.ndarray | None  # (4, 3) N/CA/C/O, or None if any is missing
    missing: tuple[str, ...]  # backbone atom names that are absent
    letter: str  # one-letter code: upper for standard, gemmi's lower for modified, X unknown


def _is_amino_acid(res: gemmi.Residue) -> bool:
    """An amino-acid residue: tabulated as one by gemmi (covers MSE and other
    modified residues written as HETATM), or written as an ATOM record (covers
    names gemmi does not tabulate). Waters, ions and ligands are neither."""
    info = gemmi.find_tabulated_residue(res.name)
    if info is not None and info.is_amino_acid():
        return True
    return res.het_flag == "A"


def one_letter(name: str) -> str:
    """One-letter code for the sequence column. Standard residues upper-case; a
    modified residue keeps gemmi's LOWER-case parent letter (MSE -> 'm'), so it is
    visible rather than silently normalised; anything gemmi has no letter for is 'X'."""
    info = gemmi.find_tabulated_residue(name)
    code = info.one_letter_code if info is not None else " "
    return code if code.strip() else "X"


def _ca_atom(res: gemmi.Residue) -> gemmi.Atom | None:
    """The residue's alpha carbon: an atom NAMED CA whose element is CARBON. The
    element test is what keeps a calcium ion (also named CA) out of the trace."""
    for atom in res:
        if atom.name == "CA" and atom.element.name.upper() == "C":
            return atom
    return None


def read_model(path: Path) -> gemmi.Model:
    """First model of a PDB or mmCIF file, one altloc per atom, no hydrogens."""
    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    st = gemmi.read_structure(str(path))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    if len(st) == 0:
        raise ValueError(f"no model in {path.name}")
    return st[0]


def select(
    model: gemmi.Model, chains: list[str], src_name: str
) -> tuple[dict[str, list[Residue]], list[str]]:
    """Amino-acid residues WITH a CA atom, per selected chain, in file order, plus the
    labels of amino-acid residues that have NO CA (so the caller can flag them --
    without this, a CA-less residue at a chain END would vanish from n_res with no
    status, since the break test only sees interior gaps).

    Raises if a named chain is absent -- never a fallback to another chain, and never
    a quietly smaller selection. A chain that is present but contributes no CA is
    left empty here; the caller decides (zero CA over the whole selection is an
    error, see score())."""
    present = list(dict.fromkeys(ch.name for ch in model))
    missing = [c for c in chains if c not in present]
    if missing:
        raise ValueError(
            f"chain(s) {','.join(missing)} absent from {src_name} "
            f"(have: {','.join(present) or 'none'})"
        )

    out: dict[str, list[Residue]] = {c: [] for c in chains}
    no_ca: list[str] = []
    for name in chains:
        # gemmi can hold two chain objects with one name (a polymer and its
        # ligands after a TER). Walk all of them, in file order.
        for chain in model:
            if chain.name != name:
                continue
            for res in chain:
                if not _is_amino_acid(res):
                    continue
                ca = _ca_atom(res)
                if ca is None:
                    no_ca.append(f"{chain.name}{res.seqid.num}{res.seqid.icode.strip()}")
                    continue
                atoms = {a.name: a for a in res}
                missing_bb = tuple(a for a in BACKBONE if a not in atoms)
                backbone = (
                    None
                    if missing_bb
                    else np.array(
                        [[atoms[a].pos.x, atoms[a].pos.y, atoms[a].pos.z] for a in BACKBONE],
                        dtype=float,
                    )
                )
                out[name].append(
                    Residue(
                        chain=chain.name,
                        resnum=f"{res.seqid.num}{res.seqid.icode.strip()}",
                        resname=res.name.strip(),
                        ca=np.array([ca.pos.x, ca.pos.y, ca.pos.z], dtype=float),
                        backbone=backbone,
                        missing=missing_bb,
                        letter=one_letter(res.name),
                    )
                )
    return out, no_ca


def radius_of_gyration(ca: np.ndarray) -> float:
    """Mass-unweighted Rg of a (N, 3) coordinate set:
    sqrt(mean(|r_i - mean(r)|^2))."""
    if ca.ndim != 2 or ca.shape[1] != 3 or len(ca) == 0:
        raise ValueError(f"radius_of_gyration needs a non-empty (N, 3) array, got {ca.shape}")
    centred = ca - ca.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum(centred * centred, axis=1))))


def chain_breaks(residues: list[Residue]) -> list[str]:
    """Labels of residues followed by a C(i)-N(i+1) gap > MAX_PEPTIDE_BOND, within
    one chain. Needs the full backbone; the caller checks that first."""
    breaks: list[str] = []
    for a, b in zip(residues, residues[1:]):
        assert a.backbone is not None and b.backbone is not None
        if float(np.linalg.norm(b.backbone[0] - a.backbone[2])) > MAX_PEPTIDE_BOND:
            breaks.append(f"{a.chain}{a.resnum}")
    return breaks


def assign_ss(per_chain: list[list[Residue]]) -> list[str]:
    """3-state SS ('H'/'E'/'L') per residue, in order, chains kept separate.

    The same call wu.dssp makes: mdtraj.compute_dssp(simplified=True) over N/CA/C/O
    atoms of residues all named GLY, 'C' renamed to 'L'. Built in memory (nm), not
    through a temporary PDB. Raises if mdtraj returns anything but H/E/C for a
    residue -- an 'NA' code would otherwise be counted as loop."""
    import mdtraj as md  # imported here so the module's pure helpers load without it

    top = md.Topology()
    xyz: list[np.ndarray] = []
    for residues in per_chain:
        if not residues:
            continue
        chain = top.add_chain()
        for i, res in enumerate(residues):
            assert res.backbone is not None
            r = top.add_residue("GLY", chain, resSeq=i)
            for atom_name, element in zip(BACKBONE, ("N", "C", "C", "O")):
                top.add_atom(atom_name, md.element.get_by_symbol(element), r)
            xyz.append(res.backbone)
    coords = np.concatenate(xyz, axis=0)[None, :, :] / 10.0  # A -> nm
    codes = md.compute_dssp(md.Trajectory(coords, top), simplified=True)[0]

    n_expected = sum(len(r) for r in per_chain)
    if len(codes) != n_expected:
        raise ValueError(
            f"DSSP returned {len(codes)} assignments for {n_expected} residues"
        )
    bad = sorted({str(c) for c in codes if str(c) not in SS_CODES})
    if bad:
        raise ValueError(f"DSSP returned non-protein code(s) {bad}")
    return [SS_CODES[str(c)] for c in codes]


def score(src: Path, chains: list[str]) -> tuple[str, dict, list[list]]:
    """One structure -> (status, metrics, per-residue rows). Raises for an 'error:'
    (nothing measurable); returns a 'warn: ...' status when Rg is measured but the SS
    fractions are not."""
    if not chains:
        raise ValueError("no chains selected")
    model = read_model(src)
    per_chain_map, no_ca = select(model, chains, src.name)
    per_chain = [per_chain_map[c] for c in chains]
    residues = [r for rs in per_chain for r in rs]
    # Per chain, not just overall: with 'A,L' where L is ligand-only, an overall check
    # would pass and n_res / Rg would silently describe A alone while chains_used
    # claimed A,L -- a quietly smaller selection.
    empty = [c for c, rs in zip(chains, per_chain) if not rs]
    if empty:
        raise ValueError(
            f"chain(s) {','.join(empty)} hold no amino-acid residue with a CA atom "
            f"in {src.name}"
        )

    ca = np.stack([r.ca for r in residues])
    metrics: dict = {
        "radius_of_gyration": round(radius_of_gyration(ca), 4),
        "helix_fraction": None,
        "beta_fraction": None,
        "loop_fraction": None,
        "n_res": len(residues),
        "n_chain_breaks": None,
        "chains_used": ",".join(chains),
        # The wrong-chain guard n_res cannot be when binder and target have similar
        # lengths: the selection's own sequence, chains joined by '/'. On a backbone
        # generator's output this is whatever residue names it wrote (often poly-G).
        "sequence": "/".join("".join(r.letter for r in rs) for rs in per_chain),
    }
    rows = [[r.chain, r.resnum, r.resname, ""] for r in residues]

    incomplete = [f"{lab}(CA)" for lab in no_ca] + [
        f"{r.chain}{r.resnum}({'/'.join(r.missing)})" for r in residues if r.missing
    ]
    if incomplete:
        shown = ",".join(incomplete[:5]) + (",..." if len(incomplete) > 5 else "")
        return (
            f"warn: {len(incomplete)} residue(s) missing backbone atoms, SS not "
            f"assigned: {shown}",
            metrics,
            rows,
        )

    breaks = [b for rs in per_chain for b in chain_breaks(rs)]
    metrics["n_chain_breaks"] = len(breaks)
    if breaks:
        shown = ",".join(breaks[:5]) + (",..." if len(breaks) > 5 else "")
        return (
            f"warn: {len(breaks)} backbone break(s) (C-N > {MAX_PEPTIDE_BOND} A after "
            f"{shown}), SS not assigned",
            metrics,
            rows,
        )

    try:
        ss = assign_ss(per_chain)
    except Exception as e:  # noqa: BLE001 - a DSSP failure is a warn, Rg stands
        return f"warn: SS assignment failed: {e}", metrics, rows

    n = len(ss)
    metrics["helix_fraction"] = round(ss.count("H") / n, 6)
    metrics["beta_fraction"] = round(ss.count("E") / n, 6)
    metrics["loop_fraction"] = round(ss.count("L") / n, 6)
    for row, code in zip(rows, ss):
        row[3] = code
    return "OK", metrics, rows


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def write_result(out_dir: Path, name: str, status: str, path: str, metrics: dict) -> None:
    """One-row TSV. A metric that was not measured is an empty cell (NA), never 0."""
    write_tsv(
        out_dir / f"{name}.tsv",
        RESULT_COLUMNS,
        [
            [name, status, path]
            + ["" if metrics.get(c) is None else str(metrics[c]) for c in METRIC_COLUMNS]
        ],
    )


class Args(argparse.Namespace):
    task_file: Path
    design_chains: str
    input_column: str
    out_dir: Path


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Per-design CA radius of gyration and 3-state DSSP fractions of a chain selection."
    )
    ap.add_argument(
        "--task-file",
        type=Path,
        required=True,
        help="Tab-separated 'name<TAB>structure' rows, one per design. An empty "
        "structure field means the manifest builder could not resolve a path.",
    )
    ap.add_argument(
        "--design-chains",
        required=True,
        help="Comma-joined chain IDs, already expanded by the manifest builder.",
    )
    ap.add_argument(
        "--input-column",
        default="the input structure column",
        help="Name of the table column the paths came from; used only in the "
        "'could not resolve' error message.",
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(namespace=Args())

    chains = split_list(args.design_chains)
    if not chains:
        raise ValueError("--design-chains is empty")

    # Fail the whole task (not each row) if the SS engine is not importable: that is
    # an environment fault, and N rows of 'warn: SS assignment failed' would read
    # like a property of the designs. The .sh wrapper turns the crash into error rows.
    import mdtraj  # noqa: F401

    with open(args.task_file) as f:
        members = [line.rstrip("\n").split("\t") for line in f if line.strip()]

    for member in members:
        if len(member) < 2:
            raise ValueError(
                f"malformed task row {member!r}: expected 'name<TAB>structure'"
            )
        name, src = member[0], member[1]
        res_path = args.out_dir / f"{name}_residues.tsv"
        # Provenance, written on EVERY row including errors: which column and which
        # file this row was measured from, so two runs sharing a label (or a stale
        # file) are distinguishable from the table alone.
        provenance = {"input_column": args.input_column, "source_path": src or None}
        try:
            if not src.strip():
                raise ValueError(
                    f"could not resolve a structure path from column "
                    f"{args.input_column!r} for this design -- neither on its own "
                    f"row nor on any ancestor (parent_name/parent_table). Check "
                    f"-i/--input-column"
                )
            status, metrics, rows = score(Path(src), chains)
            metrics.update(provenance)
            write_tsv(res_path, RESIDUE_COLUMNS, rows)
            write_result(args.out_dir, name, status, str(res_path), metrics)
            print(
                f"{name}: rg={metrics['radius_of_gyration']} "
                f"H={metrics['helix_fraction']} E={metrics['beta_fraction']} "
                f"L={metrics['loop_fraction']} n_res={metrics['n_res']} "
                f"chains={metrics['chains_used']} {status}"
            )
        except Exception as e:  # noqa: BLE001 - errors are recorded as data
            print(f"{name}: ERROR {e}", file=sys.stderr)
            res_path.unlink(missing_ok=True)
            write_result(args.out_dir, name, f"error: {e}", "", provenance)


if __name__ == "__main__":
    main()
