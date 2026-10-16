#!/usr/bin/env python3
"""
Self-test for dockrmsd: synthetic C3 complexes with a hand-checkable answer, the
manifest builder against a real prosapia DataManager lineage walk, the collector,
and the task script executed through the REAL prosapia task prelude.

WHAT THIS PROVES, AND WHAT IT CANNOT.

It proves the arithmetic and the contracts: that a rigidly moved complex gives
docking RMSD 0; that the transform is computed on the TARGET ONLY (moving the
binder of the model does not change target_fit_rmsd); that the binder is moved
ONCE and not re-fitted (a translated binder keeps its translation); that the
permutation search finds the symmetry copy when the binder sits at an equivalent
site; that on a C3 trimer exactly 3 of the 6 permutations survive the fit
tolerance; that the target correspondence is by AUTHOR NUMBER, so extra modelled
residues outside the reference's numbering never enter the fit; that unequal
binder lengths are an error and never a trimmed comparison; that a binder LONGER
than a target protomer does not invert the roles; that a wrong chain flag errors
instead of returning a number; that an independent quaternion (Horn) fit agrees;
that a mirror image is not fitted; that every refusal path writes a non-OK status
with every SCIENCE column NA; and that the manifest fields reach the worker in the
order dockrmsd.sh cuts them.

It CANNOT prove the premise or the conventions of real files. The structures here
are written by this file and agree with it by construction. UNKNOWN until real
data:

  * the ten real docking RMSDs themselves. The pre-flight established the file
    shapes (model target chains 1-149, reference 1-139, binder 118 both sides,
    one model, no insertion codes, single-letter chain ids), and the hand-run
    script has since been read: it uses CA only, ONE Kabsch over all three
    protomers pooled, and the minimum over its TOP 3 permutations BY TARGET FIT
    -- the same basis as this tool, since on a true C3 "best three" and "within
    --fit-tolerance of the best" are the same set. So EXACT agreement (<0.05 A)
    is the acceptance bar. Nothing here has compared an output with any of the
    ten;
  * NOTE the one place the two rules could part company: a target where four or
    more mappings fit comparably. The script's [:3] is a hard-coded assumption;
    this tool's tolerance window is not. test_identity_and_warn covers the
    tolerance behaviour synthetically, but no real target has exercised it;
  * whether --min-target-identity 0.8 is the right floor UNDER SEQUENCE MODE,
    where it now measures "is this the same protein" rather than doubling as a
    numbering check. The pilot pair's 0.698 under `number` became the reason the
    mode exists; its value under `sequence` has not been measured on real files;
  * whether 1.0 A is the right --fit-tolerance on real oligomers: it was chosen
    from one measured C3 case (valid permutations ~1.05 A, impostors ~22.4 A), a
    gap of 20 A, and has not been exercised on a near-symmetric target where two
    inequivalent mappings fit almost as well (--pin-mapping is the escape hatch,
    and it too is untested on real data).

Run:  uv run python tools/dockrmsd/test_dockrmsd.py
"""

import importlib.util
import os
import subprocess
import sys
import tempfile
from argparse import Namespace
from functools import partial
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from prosapia.core import DataManager

HERE = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


worker = _load("dockrmsd_worker")
run_dockrmsd = _load("run_dockrmsd")
collect_dockrmsd = _load("collect_dockrmsd")

FAILURES: list[str] = []
CHECKS = 0


def check(label: str, got, want, tol: float | None = None) -> None:
    global CHECKS
    CHECKS += 1
    try:
        ok = abs(float(got) - float(want)) <= tol if tol is not None else got == want
    except (TypeError, ValueError):
        ok = False
    if not ok:
        FAILURES.append(f"{label}: got {got!r}, want {want!r}")


def check_true(label: str, value) -> None:
    global CHECKS
    CHECKS += 1
    if not value:
        FAILURES.append(f"{label}: expected truthy, got {value!r}")


def check_raises(label: str, fn, needle: str, exc=ValueError) -> None:
    global CHECKS
    CHECKS += 1
    try:
        fn()
    except exc as e:
        if needle not in str(e):
            FAILURES.append(f"{label}: message {str(e)!r} lacks {needle!r}")
        return
    except Exception as e:  # noqa: BLE001
        FAILURES.append(f"{label}: raised {type(e).__name__}: {e}")
        return
    FAILURES.append(f"{label}: did not raise")


# ---------------------------------------------------------------------------
# Synthetic structures
# ---------------------------------------------------------------------------

AA3 = {
    "A": "ALA", "G": "GLY", "S": "SER", "V": "VAL", "L": "LEU", "W": "TRP",
    "K": "LYS", "E": "GLU", "D": "ASP", "F": "PHE", "I": "ILE", "T": "THR",
    "N": "ASN", "Q": "GLN", "Y": "TYR", "R": "ARG", "H": "HIS", "P": "PRO",
    "M": "MET", "C": "CYS",
}


def helix(n: int, radius: float = 2.3, pitch: float = 1.5) -> np.ndarray:
    t = np.arange(n)
    return np.stack(
        [radius * np.cos(t * 1.745), radius * np.sin(t * 1.745), pitch * t], axis=1
    )


def rotation(axis, angle: float) -> np.ndarray:
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    x, y, z = axis
    c, s = np.cos(angle), np.sin(angle)
    C = 1 - c
    return np.array([
        [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


R120 = rotation([0, 0, 1], 2 * np.pi / 3)
R240 = rotation([0, 0, 1], 4 * np.pi / 3)
RIGID_R = rotation([1, 2, 3], 0.7)
RIGID_T = np.array([5.0, -3.0, 12.0])
# PDB files store coordinates with 3 decimals, so a write/read round trip of a
# rotated structure costs ~1e-3 A. Every "lands exactly" check uses this.
PDB_TOL = 3e-3

PROTOMER_SEQ = "AGSVLWKEAGSVLWKEAGSVLWKEAGSVLW"  # 30
BINDER_SEQ = "DFITNQYRHPMCDFITNQYRHPMCDF"  # 26
# The binder's "designed" sequence: the SAME backbone, a different sequence --
# the relationship dockrmsd must not require.
BINDER_SEQ_REDESIGNED = "KKITNQYRAAMCDFITAAYRHPMCEE"

# One protomer, pushed off the symmetry axis so the three copies are distinct.
PROTOMER = helix(len(PROTOMER_SEQ)) + np.array([14.0, 0.0, 0.0])
# A binder sitting against protomer A (NOT on the axis, so the three sites differ).
BINDER = helix(len(BINDER_SEQ), radius=1.8, pitch=1.9) + np.array([24.0, 6.0, 4.0])


OFFSETS = {
    "CA": np.zeros(3),
    "N": np.array([0.6, 0.2, -0.3]),
    "C": np.array([-0.5, 0.4, 0.3]),
    "O": np.array([-0.2, -0.7, 0.5]),
}


def residues(coords, seq, start: int = 1, atoms=("CA",), rot=None, trans=None):
    """[(resnum, one-letter, {atom: xyz})], one residue per coordinate.

    Backbone atoms are fixed offsets from the CA IN THE LOCAL FRAME, and the
    chain transform is applied to the finished residue -- so every atom is
    rigidly attached to its CA. Adding the offsets after the transform instead
    would make a "rigidly moved" model non-rigid, and the C3 copies
    inequivalent, which no backbone-mode test could then pass.
    """
    out = []
    for i, (xyz, code) in enumerate(zip(coords, seq)):
        pos = {a: move(np.atleast_2d(xyz + OFFSETS[a]), rot, trans)[0] for a in atoms}
        out.append((start + i, code, pos))
    return out


def write_structure(path: Path, chains: dict) -> None:
    """``chains`` = {chain: [(resnum, one-letter, {atom_name: xyz})]}."""
    st = gemmi.Structure()
    model = gemmi.Model("1")
    for chain_name, res_list in chains.items():
        chain = gemmi.Chain(chain_name)
        for resnum, code, atoms in res_list:
            res = gemmi.Residue()
            res.name = AA3[code]
            res.seqid = gemmi.SeqId(int(resnum), " ")
            res.het_flag = "A"
            for atom_name, xyz in atoms.items():
                atom = gemmi.Atom()
                atom.name = atom_name
                atom.element = gemmi.Element("N" if atom_name == "N" else
                                             "O" if atom_name == "O" else "C")
                atom.pos = gemmi.Position(*map(float, xyz))
                res.add_atom(atom)
            chain.add_residue(res)
        model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".cif":
        st.make_mmcif_document().write_file(str(path))
    else:
        st.write_pdb(str(path))


def move(coords: np.ndarray, rot=None, trans=None) -> np.ndarray:
    out = coords if rot is None else coords @ rot.T
    return out if trans is None else out + trans


def trimer(binder_coords, binder_seq=BINDER_SEQ, start=1, protomer=PROTOMER,
           rot=None, trans=None, target_extra=0, atoms=("CA",),
           chain_names=("A", "B", "C"), binder_chain="D"):
    """A C3 target (3 protomers 120 deg apart) plus one binder chain, all moved
    rigidly by (rot, trans). ``target_extra`` appends N residues NUMBERED BEYOND
    the reference's range, at absurd coordinates -- they must never enter a fit."""
    chains = {}
    for i, name in enumerate(chain_names):
        # One transform per chain: the C3 operation, then the rigid move.
        sym = np.linalg.matrix_power(R120, i)
        total = sym if rot is None else rot @ sym
        res = residues(protomer, PROTOMER_SEQ, start=start, atoms=atoms,
                       rot=total, trans=trans)
        if target_extra:
            junk = np.tile(np.array([900.0, 900.0, 900.0]), (target_extra, 1))
            res += residues(junk, "G" * target_extra,
                            start=start + len(PROTOMER_SEQ), atoms=atoms,
                            rot=total, trans=trans)
        chains[name] = res
    chains[binder_chain] = residues(
        binder_coords, binder_seq, start=1, atoms=atoms, rot=rot, trans=trans
    )
    return chains


# ---------------------------------------------------------------------------
# Worker helpers
# ---------------------------------------------------------------------------


def options(**over):
    fields = dict(
        binder_chain="D", target_chains=["A", "B", "C"],
        ref_binder_chain="D", ref_target_chains=["A", "B", "C"],
        fit_atoms="ca", target_match="number", min_target_res=10,
        min_target_res_frac=0.75, min_target_identity=0.8,
        fit_tolerance=1.0, max_target_fit_rmsd=5.0, pin_mapping=False,
    )
    fields.update(over)
    return worker.Options(**fields)


def measure(model: Path, ref: Path, **over):
    return worker.compute(model, ref, options(**over))


# ---------------------------------------------------------------------------
# 1. A known-good case: the model is the reference, rigidly moved
# ---------------------------------------------------------------------------


def test_known_good(tmp: Path) -> None:
    ref, model = tmp / "kg_ref.pdb", tmp / "kg_model.cif"
    write_structure(ref, trimer(BINDER))
    write_structure(model, trimer(BINDER, rot=RIGID_R, trans=RIGID_T))

    res = measure(model, ref)
    m = res.metrics
    check("known-good: status OK", res.status, "OK")
    check("known-good: docking rmsd ~0", m["rmsd"], 0.0, tol=PDB_TOL)
    check("known-good: target fit ~0", m["target_fit_rmsd"], 0.0, tol=PDB_TOL)
    check("known-good: max_dev ~0", m["max_dev"], 0.0, tol=PDB_TOL)
    check("known-good: complex rmsd ~0", m["complex_rmsd"], 0.0, tol=PDB_TOL)
    check("known-good: identity permutation wins", m["permutation"], "A>A,B>B,C>C")
    check("known-good: 6 permutations enumerated", m["n_permutations"], 6)
    check("known-good: C3 admits exactly 3 candidates", m["n_candidates"], 3)
    check("known-good: residues used", m["n_target_res_used"], 3 * len(PROTOMER_SEQ))
    check("known-good: binder lengths", m["binder_n_res"], len(BINDER_SEQ))
    check("known-good: target identity 1.0", m["target_seq_match_frac"], 1.0, tol=1e-9)
    check("known-good: used frac 1.0", m["target_res_used_frac"], 1.0, tol=1e-9)
    check("known-good: fit atoms = 1 per residue", m["n_fit_atoms"], 3 * len(PROTOMER_SEQ))
    check("known-good: chain ids echoed", m["model_chain_ids"], "A,B,C,D")
    check_true("known-good: margin is a real gap", float(m["margin"]) > 5.0)
    check_true("known-good: rmsd_all lists all 6",
               m["rmsd_all"].count(";") == 5 and "*A>A,B>B,C>C:fit=" in m["rmsd_all"])
    check("known-good: exactly the candidates are starred",
          m["rmsd_all"].count("*"), int(m["n_candidates"]))
    check("known-good: per-chain split is in the audit",
          m["rmsd_all"].count("per_chain=30/30/30"), 6)
    check("known-good: per-chain split is a column",
          m["n_target_res_used_per_chain"], "30/30/30")
    check("known-good: numbering provenance, model", m["target_resnum_range_model"],
          "1-30;1-30;1-30")
    check("known-good: numbering provenance, reference", m["target_resnum_range_ref"],
          "1-30;1-30;1-30")
    check("known-good: target match echoed", m["target_match"], "number")
    check("known-good: pin_mapping echoed", m["pin_mapping"], 0)
    # The three impostor permutations must be visible AND excluded, not silently
    # dropped: they are in rmsd_all with a fit of tens of angstroms.
    fits = [float(p.split(":fit=")[1].split(":")[0])
            for p in m["rmsd_all"].split(";") if ":fit=" in p]
    check("known-good: every permutation produced a fit", len(fits), 6)
    check("known-good: 3 good fits", sum(f < 1.0 for f in fits), 3)
    check_true("known-good: 3 impostor fits are far worse",
               sum(f > 5.0 for f in fits) == 3)

    # The superposed model really is in the reference's frame.
    sup = gemmi.read_structure(str(model))
    moved = res.moved
    check("known-good: moved structure keeps every chain", len(moved[0]), len(sup[0]))
    ref_st = gemmi.read_structure(str(ref))
    d = max(
        (moved[0][ci][ri][0].pos - ref_st[0][ci][ri][0].pos).length()
        for ci in range(4) for ri in range(len(moved[0][ci]))
    )
    check("known-good: sup structure lands on the reference", d, 0.0, tol=PDB_TOL)


# ---------------------------------------------------------------------------
# 2. The transform is the TARGET's, and the binder is moved ONCE
# ---------------------------------------------------------------------------


def test_binder_not_fitted(tmp: Path) -> None:
    ref, model = tmp / "nf_ref.pdb", tmp / "nf_model.pdb"
    shift = np.array([7.0, 0.0, 0.0])
    write_structure(ref, trimer(BINDER))
    # Same target, binder translated by |shift| -- a rigid-body displacement of
    # the binder alone, which is exactly what a docking error is.
    write_structure(model, trimer(BINDER + shift))

    m = measure(model, ref).metrics
    check("no-refit: target fit unaffected by the binder", m["target_fit_rmsd"], 0.0, tol=PDB_TOL)
    check("no-refit: the displacement is reported in full",
          m["rmsd"], float(np.linalg.norm(shift)), tol=PDB_TOL)
    check("no-refit: max_dev equals the displacement",
          m["max_dev"], float(np.linalg.norm(shift)), tol=PDB_TOL)
    # A tool that re-fitted the binder would report ~0 here.
    check_true("no-refit: not re-fitted (would be ~0)", float(m["rmsd"]) > 6.9)

    # A binder rotated about its own centroid: again no refit, so a real number.
    centroid = BINDER.mean(axis=0)
    spun = (BINDER - centroid) @ rotation([0, 1, 0], 0.9).T + centroid
    write_structure(model, trimer(spun))
    m2 = measure(model, ref).metrics
    check("no-refit: spun binder keeps the target fit", m2["target_fit_rmsd"], 0.0, tol=PDB_TOL)
    check_true("no-refit: spun binder has a large rmsd", float(m2["rmsd"]) > 3.0)


# ---------------------------------------------------------------------------
# 3. The symmetry search: an equivalent site must be found
# ---------------------------------------------------------------------------


def test_symmetry_search(tmp: Path) -> None:
    ref, model = tmp / "sym_ref.pdb", tmp / "sym_model.pdb"
    write_structure(ref, trimer(BINDER))
    # Model: identical target, binder moved to the NEXT equivalent site (R120).
    write_structure(model, trimer(BINDER @ R120.T))

    m = measure(model, ref).metrics
    check("symmetry: found an equivalent site (rmsd ~0)", m["rmsd"], 0.0, tol=PDB_TOL)
    # T must undo R120, i.e. model A lands on ref C.
    check("symmetry: the winning permutation (REFERENCE order)", m["permutation"],
          "B>A,C>B,A>C")
    check("symmetry: the winner's target fit is still ~0",
          m["target_fit_rmsd"], 0.0, tol=PDB_TOL)
    check("symmetry: 3 candidates", m["n_candidates"], 3)
    check_true("symmetry: the identity permutation was NOT the answer",
               "A>A,B>B,C>C:fit=" in m["rmsd_all"])
    identity_rmsd = float(
        [p for p in m["rmsd_all"].split(";") if p.lstrip("*").startswith("A>A,B>B,C>C")][0]
        .split(":binder=")[1].split(":")[0]
    )
    check_true("symmetry: identity permutation would have been badly wrong",
               identity_rmsd > 10.0)
    check_true("symmetry: margin is decisive", float(m["margin"]) > 5.0)

    # The same model scored WITHOUT the symmetry search (one target chain only)
    # must give the big number -- the search is what makes the difference.
    m1 = measure(model, ref, target_chains=["A"], ref_target_chains=["A"]).metrics
    check("symmetry: single-chain target has 1 permutation", m1["n_permutations"], 1)
    check("symmetry: a single candidate leaves margin NA (not 0)", m1["margin"], None)
    check("symmetry: without the search the answer is the wrong one",
          m1["rmsd"], identity_rmsd, tol=1e-2)


# ---------------------------------------------------------------------------
# 4. Equal-length binders are enforced
# ---------------------------------------------------------------------------


def test_binder_length_enforced(tmp: Path) -> None:
    ref, model = tmp / "bl_ref.pdb", tmp / "bl_model.pdb"
    write_structure(ref, trimer(BINDER))
    write_structure(model, trimer(BINDER[:-3], binder_seq=BINDER_SEQ[:-3]))

    try:
        measure(model, ref)
        FAILURES.append("binder length: did not raise")
        globals()["CHECKS"] = CHECKS + 1
    except worker.MeasureError as e:
        check_true("binder length: refused", "binder length mismatch" in str(e))
        check_true("binder length: names both counts",
                   "26" in str(e) and "23" in str(e))
        check("binder length: trust columns survive the error",
              e.metrics.get("binder_n_res"), 23)
        check("binder length: reference count too", e.metrics.get("ref_binder_n_res"), 26)
        check_true("binder length: no science on the error",
                   all(c not in e.metrics for c in ("rmsd", "max_dev", "complex_rmsd")))


# ---------------------------------------------------------------------------
# 5. The role-flip case: a binder longer than a target protomer
# ---------------------------------------------------------------------------


def test_role_flip(tmp: Path) -> None:
    """DockQ assigns receptor/ligand by LENGTH, so a binder longer than a target
    protomer inverts its measurement. dockrmsd takes the roles from the flags."""
    long_seq = (PROTOMER_SEQ + BINDER_SEQ)[: len(PROTOMER_SEQ) + 14]  # 44 > 30
    long_binder = helix(len(long_seq), radius=1.8, pitch=1.9) + np.array([24.0, 6.0, 4.0])
    ref, model = tmp / "rf_ref.pdb", tmp / "rf_model.pdb"
    write_structure(ref, trimer(long_binder, binder_seq=long_seq))
    shift = np.array([0.0, 4.0, 0.0])
    write_structure(
        model, trimer(long_binder + shift, binder_seq=long_seq, rot=RIGID_R, trans=RIGID_T)
    )

    m = measure(model, ref).metrics
    check_true("role flip: binder IS longer than a protomer",
               int(m["binder_n_res"]) > len(PROTOMER_SEQ))
    check("role flip: the binder's displacement is what is reported",
          m["rmsd"], float(np.linalg.norm(shift)), tol=PDB_TOL)
    check("role flip: the TARGET still set the frame", m["target_fit_rmsd"], 0.0, tol=PDB_TOL)
    check("role flip: the fit used target residues only",
          m["n_fit_atoms"], 3 * len(PROTOMER_SEQ))
    check("role flip: binder atoms are the binder's",
          m["n_binder_atoms"], len(long_seq))


# ---------------------------------------------------------------------------
# 6. Wrong chain arguments error instead of returning a number
# ---------------------------------------------------------------------------


def test_wrong_chains(tmp: Path) -> None:
    ref, model = tmp / "wc_ref.pdb", tmp / "wc_model.pdb"
    write_structure(ref, trimer(BINDER))
    write_structure(model, trimer(BINDER, rot=RIGID_R, trans=RIGID_T))

    cases = [
        (dict(binder_chain="Z"), "chain(s) Z (--binder-chain)"),
        (dict(ref_binder_chain="Z"), "chain(s) Z (--ref-binder-chain)"),
        (dict(target_chains=["A", "B", "Z"]), "chain(s) Z (--target-chains)"),
        (dict(ref_target_chains=["A", "B", "Z"]), "chain(s) Z (--ref-target-chains)"),
    ]
    for over, needle in cases:
        check_raises(f"wrong chain {over}: errors, never a number",
                     lambda o=over: measure(model, ref, **o), needle,
                     worker.MeasureError)

    # Naming the BINDER as one of the target chains: refused at submit (see the
    # refusals table in test_manifest), and if it ever got past, the identity
    # floor refuses it rather than fitting the binder onto a protomer.
    check_raises("wrong chain: binder as a target is caught by the identity floor",
                 lambda: measure(model, ref, target_chains=["A", "B", "D"],
                                 ref_target_chains=["A", "B", "C"]),
                 "no chain permutation gave a usable target superposition",
                 worker.MeasureError)


# ---------------------------------------------------------------------------
# 7. Target correspondence is by author number
# ---------------------------------------------------------------------------


def test_author_numbering(tmp: Path) -> None:
    ref, model = tmp / "an_ref.pdb", tmp / "an_model.pdb"
    write_structure(ref, trimer(BINDER))
    # The model protomers carry 8 extra residues, numbered past the reference's
    # range, at absurd coordinates. They are outside the intersection, so they
    # must not touch the fit.
    write_structure(
        model, trimer(BINDER, rot=RIGID_R, trans=RIGID_T, target_extra=8)
    )
    m = measure(model, ref).metrics
    check("numbering: extra modelled residues never enter the fit",
          m["target_fit_rmsd"], 0.0, tol=PDB_TOL)
    check("numbering: docking rmsd unaffected", m["rmsd"], 0.0, tol=PDB_TOL)
    check("numbering: used = the shared residues", m["n_target_res_used"], 90)
    check("numbering: model total counts the extras", m["n_target_res_model"], 114)
    check("numbering: reference total", m["n_target_res_ref"], 90)
    check("numbering: used frac vs the shorter side", m["target_res_used_frac"], 1.0, tol=1e-9)

    # Disjoint numbering (the rfd3 "renumbered from 1" trap) is an error naming it.
    write_structure(model, trimer(BINDER, start=501, rot=RIGID_R, trans=RIGID_T))
    check_raises("numbering: disjoint numbering errors, never a smaller fit",
                 lambda: measure(model, ref),
                 "no chain permutation gave a usable target superposition",
                 worker.MeasureError)
    try:
        measure(model, ref)
    except worker.MeasureError as e:
        check_true("numbering: the error row keeps the permutation audit",
                   "skipped_few_shared_residues" in str(e.metrics.get("rmsd_all", "")))


# ---------------------------------------------------------------------------
# 7b. The campaign's real geometry: an inserted loop, and the two match modes
# ---------------------------------------------------------------------------


def insert_loop(coords, seq, at: int, length: int):
    """A modelled loop the other file lacks: `length` extra residues spliced in
    at 0-based position `at`, renumbered 1..N+length as a predictor would."""
    bulge = coords[at] + np.stack([
        np.array([3.0, 0.0, 0.0]) * (i + 1) for i in range(length)
    ])
    return (
        np.vstack([coords[:at], bulge, coords[at:]]),
        seq[:at] + "RETPEGAEAK"[:length] + seq[at:],
    )


def test_inserted_loop(tmp: Path) -> None:
    """Measured on the real pilot pair: the Boltz model numbers its target chains
    1-149 and the rfd3 reference numbers its 1-139, BOTH FROM 1, and the model
    carries a 10-residue loop the reference does not. Author-number matching
    cannot work -- not because of a numbering convention but because the model
    models residues the reference lacks."""
    ref, model = tmp / "il_ref.pdb", tmp / "il_model.pdb"
    at, extra = 20, 10  # 20/30 is the pilot's 95/139 in miniature
    long_coords, long_seq = insert_loop(PROTOMER, PROTOMER_SEQ, at, extra)
    write_structure(ref, trimer(BINDER))
    write_structure(
        model, trimer(BINDER, protomer=long_coords, rot=RIGID_R, trans=RIGID_T)
    )
    # The inserted loop carries its own sequence, so rebuild the model with it.
    chains = trimer(BINDER, rot=RIGID_R, trans=RIGID_T)
    for i, name in enumerate(("A", "B", "C")):
        sym = np.linalg.matrix_power(R120, i)
        chains[name] = residues(long_coords, long_seq, start=1,
                                rot=RIGID_R @ sym, trans=RIGID_T)
    write_structure(model, chains)

    # --target-match number: the shared numbers 1..30 pair positions 1-20
    # correctly and 21-30 against the shifted tail, so identity collapses.
    try:
        measure(model, ref, target_match="number")
        FAILURES.append("inserted loop: author numbering did NOT error")
        globals()["CHECKS"] = CHECKS + 1
    except worker.MeasureError as e:
        check_true("inserted loop: author numbering errors loudly",
                   "no chain permutation gave a usable" in str(e))
        check_true("inserted loop: the audit names the identity floor",
                   "skipped_identity" in str(e.metrics.get("rmsd_all", "")))
        check_true("inserted loop: and the two numbering ranges are recorded",
                   e.metrics.get("target_resnum_range_model") == "1-40;1-40;1-40"
                   and e.metrics.get("target_resnum_range_ref") == "1-30;1-30;1-30")

    # --target-match sequence: the alignment places the loop in a gap and the
    # measurement works, which is the whole point of the mode.
    res = measure(model, ref, target_match="sequence")
    m = res.metrics
    check("inserted loop: sequence mode succeeds", res.status, "OK")
    check("inserted loop: target identity recovers", m["target_seq_match_frac"],
          1.0, tol=1e-9)
    check("inserted loop: every reference residue is paired",
          m["n_target_res_used"], 3 * len(PROTOMER_SEQ))
    check("inserted loop: the fit is clean", m["target_fit_rmsd"], 0.0, tol=PDB_TOL)
    check("inserted loop: and so is the docking rmsd", m["rmsd"], 0.0, tol=PDB_TOL)
    check("inserted loop: the mode is recorded", m["target_match"], "sequence")
    check("inserted loop: model total counts the loop", m["n_target_res_model"], 120)
    check("inserted loop: used frac measured against the shorter side",
          m["target_res_used_frac"], 1.0, tol=1e-9)
    check("inserted loop: the symmetry search still works", m["n_candidates"], 3)

    # The alignment itself, checked directly against the known splice. This is
    # also where the aligner is pinned to the REFERENCE IMPLEMENTATION's
    # correspondence: that script renumbered the prediction uniformly (+8) and
    # the design piecewise (+8 up to 94, +18 beyond), i.e. prefix 1:1 and suffix
    # shifted by exactly the loop length. If the aligner ever stopped producing
    # that shape, these three checks are where it would show.
    pairs = worker.align_sequences(long_seq, PROTOMER_SEQ)
    check("align: every reference residue is aligned", len(pairs), len(PROTOMER_SEQ))
    check("align: the prefix maps 1:1", [b for a, b in pairs[:at]], list(range(at)))
    check("align: the suffix is shifted by the loop length",
          [a - b for a, b in pairs[at:]], [extra] * (len(PROTOMER_SEQ) - at))
    check("align: identical sequences align 1:1",
          worker.align_sequences(PROTOMER_SEQ, PROTOMER_SEQ),
          [(i, i) for i in range(len(PROTOMER_SEQ))])

    # The BINDER is never sequence-aligned, whatever --target-match says: a
    # shorter binder is still a length error, not an alignment.
    short = dict(chains)
    short["D"] = residues(BINDER[:-3], BINDER_SEQ[:-3], rot=RIGID_R, trans=RIGID_T)
    write_structure(model, short)
    check_raises("inserted loop: sequence mode does NOT loosen the binder rule",
                 lambda: measure(model, ref, target_match="sequence"),
                 "binder length mismatch", worker.MeasureError)


def test_partial_misregistration(tmp: Path) -> None:
    """The dangerous case the identity floor only half-catches: the two files
    share author numbers across the WHOLE chain but are mis-registered from the
    midpoint. Caught here (0.67); a break nearer the C-terminus would score ~0.86
    and pass -- which is why the premise says to read the range columns."""
    ref, model = tmp / "mr_ref.pdb", tmp / "mr_model.pdb"
    write_structure(ref, trimer(BINDER))
    # Same 30 author numbers, but from position 15 the model's residues are the
    # reference's shifted by 10 -- a modelled loop that was renumbered away.
    shifted_seq = PROTOMER_SEQ[:15] + PROTOMER_SEQ[25:] + PROTOMER_SEQ[:10]
    shifted = np.vstack([PROTOMER[:15], PROTOMER[25:], PROTOMER[:10]])
    chains = {}
    for i, name in enumerate(("A", "B", "C")):
        sym = np.linalg.matrix_power(R120, i)
        chains[name] = residues(shifted, shifted_seq, start=1,
                                rot=RIGID_R @ sym, trans=RIGID_T)
    chains["D"] = residues(BINDER, BINDER_SEQ, rot=RIGID_R, trans=RIGID_T)
    write_structure(model, chains)

    try:
        measure(model, ref, target_match="number")
        FAILURES.append("mis-registration: did not error")
        globals()["CHECKS"] = CHECKS + 1
    except worker.MeasureError as e:
        check_true("mis-registration: errors rather than fitting wrong pairs",
                   "no chain permutation gave a usable" in str(e))
        check_true("mis-registration: the audit names the identity floor",
                   "skipped_identity" in str(e.metrics.get("rmsd_all", "")))
        check_true("mis-registration: the numbering ranges look IDENTICAL, which "
                   "is exactly why identity has to do the work here",
                   e.metrics.get("target_resnum_range_model")
                   == e.metrics.get("target_resnum_range_ref"))


def test_per_chain_pair_floor(tmp: Path) -> None:
    """F1: only ONE chain pair intersects, so a pooled floor would pass and the
    'whole-target' fit would silently become a single-protomer fit."""
    ref, model = tmp / "pc_ref.pdb", tmp / "pc_model.pdb"
    write_structure(ref, trimer(BINDER))
    # A predictor numbering the complex CONTINUOUSLY: A 1-30, B 31-60, C 61-90,
    # against a reference whose chains each start at 1. Only chain A overlaps.
    chains = {}
    for i, name in enumerate(("A", "B", "C")):
        sym = np.linalg.matrix_power(R120, i)
        chains[name] = residues(PROTOMER, PROTOMER_SEQ, start=1 + 30 * i,
                                rot=RIGID_R @ sym, trans=RIGID_T)
    chains["D"] = residues(BINDER, BINDER_SEQ, rot=RIGID_R, trans=RIGID_T)
    write_structure(model, chains)

    try:
        measure(model, ref, target_match="number")
        FAILURES.append("per-chain floor: a single-protomer fit was accepted")
        globals()["CHECKS"] = CHECKS + 1
    except worker.MeasureError as e:
        check_true("per-chain floor: refused, not a single-protomer fit",
                   "no chain permutation gave a usable" in str(e))
        check_true("per-chain floor: the message names the per-pair rule",
                   "PER CHAIN PAIR" in str(e))
        audit = str(e.metrics.get("rmsd_all", ""))
        check_true("per-chain floor: every permutation skipped on min_pair",
                   audit.count("skipped_few_shared_residues(min_pair=0") == 6)
        check_true("per-chain floor: the audit shows the 30/0/0 split",
                   "per_chain=30/0/0" in audit)
    # Pooled, this would have passed: 30 >= --min-target-res 10.
    check_true("per-chain floor: the pooled total would have cleared the floor",
               len(PROTOMER_SEQ) >= 10)


def test_pin_mapping(tmp: Path) -> None:
    """F3: on a pseudo-symmetric target the minimum reports the best of several
    INEQUIVALENT sites. --pin-mapping asks the honest question instead."""
    ref, model = tmp / "pm_ref.pdb", tmp / "pm_model.pdb"
    write_structure(ref, trimer(BINDER))
    write_structure(model, trimer(BINDER @ R120.T))  # binder at the next site

    free = measure(model, ref).metrics
    check("pin: the search finds the equivalent site", free["rmsd"], 0.0, tol=PDB_TOL)
    pinned = measure(model, ref, pin_mapping=True).metrics
    check("pin: only one mapping is evaluated", pinned["n_permutations"], 1)
    check("pin: it is the order as given", pinned["permutation"], "A>A,B>B,C>C")
    check("pin: margin NA with one candidate", pinned["margin"], None)
    check("pin: and the answer is the site-specific one", pinned["pin_mapping"], 1)
    check_true("pin: which is the LARGE number the search was hiding",
               float(pinned["rmsd"]) > 10.0)


# ---------------------------------------------------------------------------
# 8. The identity floor, the fit tolerance and the warn status
# ---------------------------------------------------------------------------


def test_identity_and_warn(tmp: Path) -> None:
    ref, model = tmp / "iw_ref.pdb", tmp / "iw_model.pdb"
    write_structure(ref, trimer(BINDER))
    write_structure(model, trimer(BINDER, rot=RIGID_R, trans=RIGID_T))

    # A HETERO target: chain C of the reference is given a different sequence, so
    # only the pairings that respect it clear --min-target-identity.
    hetero_ref = trimer(BINDER)
    hetero_ref["C"] = residues(
        PROTOMER, "P" * len(PROTOMER_SEQ), rot=np.linalg.matrix_power(R120, 2)
    )
    hetero_model = trimer(BINDER, rot=RIGID_R, trans=RIGID_T)
    hetero_model["C"] = residues(
        PROTOMER, "P" * len(PROTOMER_SEQ),
        rot=RIGID_R @ np.linalg.matrix_power(R120, 2), trans=RIGID_T,
    )
    write_structure(ref, hetero_ref)
    write_structure(model, hetero_model)
    m = measure(model, ref).metrics
    check("hetero: only identity-consistent pairings survive", m["n_candidates"], 1)
    check("hetero: the right one", m["permutation"], "A>A,B>B,C>C")
    check("hetero: single candidate -> margin NA", m["margin"], None)
    check_true("hetero: the rejected pairings are recorded, not hidden",
               m["rmsd_all"].count("skipped_identity") == 4)
    check("hetero: rmsd still ~0", m["rmsd"], 0.0, tol=PDB_TOL)

    # A model target that does not reproduce the reference: warn, numbers kept.
    write_structure(ref, trimer(BINDER))
    distorted = PROTOMER * np.array([1.25, 1.25, 1.35])  # commutes with R120
    write_structure(model, trimer(BINDER, protomer=distorted))
    res = measure(model, ref)
    check_true("warn: bad target fit downgrades the row",
               res.status.startswith("warn: target_fit_rmsd"))
    check_true("warn: the number is still there", res.metrics["rmsd"] is not None)
    check_true("warn: and the fit value is above the cutoff",
               float(res.metrics["target_fit_rmsd"]) > 5.0)
    res_ok = measure(model, ref, max_target_fit_rmsd=50.0)
    check("warn: the cutoff is a flag", res_ok.status, "OK")

    # --fit-tolerance is what keeps the minimum over genuine symmetry copies: a
    # huge tolerance lets the impostor permutations compete.
    write_structure(model, trimer(BINDER @ R120.T))
    write_structure(ref, trimer(BINDER))
    loose = measure(model, ref, fit_tolerance=1e6).metrics
    check("tolerance: a huge tolerance admits every permutation",
          loose["n_candidates"], 6)
    tight = measure(model, ref).metrics
    check("tolerance: the default admits only the symmetry copies",
          tight["n_candidates"], 3)
    check_true("tolerance: and only then is the winner a real fit",
               float(tight["target_fit_rmsd"]) < 1e-2)


# ---------------------------------------------------------------------------
# 9. Backbone atoms, sequence trust, geometry cross-checks
# ---------------------------------------------------------------------------


def horn_quaternion_fit(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Independent superposition via the quaternion (Horn 1987) method, used only
    to cross-check kabsch(). Returns the rotation."""
    pc, qc = p - p.mean(0), q - q.mean(0)
    h = pc.T @ qc
    trace = np.trace(h)
    k = np.array([
        [trace, h[1, 2] - h[2, 1], h[2, 0] - h[0, 2], h[0, 1] - h[1, 0]],
        [h[1, 2] - h[2, 1], h[0, 0] - h[1, 1] - h[2, 2], h[0, 1] + h[1, 0], h[0, 2] + h[2, 0]],
        [h[2, 0] - h[0, 2], h[0, 1] + h[1, 0], -h[0, 0] + h[1, 1] - h[2, 2], h[1, 2] + h[2, 1]],
        [h[0, 1] - h[1, 0], h[0, 2] + h[2, 0], h[1, 2] + h[2, 1], -h[0, 0] - h[1, 1] + h[2, 2]],
    ])
    w, v = np.linalg.eigh(k)
    q0, q1, q2, q3 = v[:, int(np.argmax(w))]
    return np.array([
        [q0*q0+q1*q1-q2*q2-q3*q3, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
        [2*(q1*q2+q0*q3), q0*q0-q1*q1+q2*q2-q3*q3, 2*(q2*q3-q0*q1)],
        [2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), q0*q0-q1*q1-q2*q2+q3*q3],
    ])


def test_geometry_and_atoms(tmp: Path) -> None:
    p = helix(40) + np.array([1.0, 2.0, 3.0])
    q = p @ RIGID_R.T + RIGID_T
    rot, trans = worker.kabsch(p, q)
    check("kabsch: recovers the rotation", float(np.abs(rot - RIGID_R).max()), 0.0, tol=1e-8)
    check("kabsch: moved points land", worker.rmsd(p @ rot.T + trans, q), 0.0, tol=1e-8)
    check("kabsch: agrees with an independent quaternion fit",
          float(np.abs(rot - horn_quaternion_fit(p, q)).max()), 0.0, tol=PDB_TOL)
    mirrored = p * np.array([1.0, 1.0, -1.0])
    rot_m, _ = worker.kabsch(mirrored, p)
    check("kabsch: no reflection", float(np.linalg.det(rot_m)), 1.0, tol=1e-9)
    check_raises("kabsch: refuses fewer than 3 points",
                 lambda: worker.kabsch(p[:2], q[:2]), "paired points",
                 worker.MeasureError)

    # --fit-atoms backbone uses 4 atoms per residue, on both sides of the measure.
    ref, model = tmp / "bb_ref.pdb", tmp / "bb_model.pdb"
    atoms = ("N", "CA", "C", "O")
    write_structure(ref, trimer(BINDER, atoms=atoms))
    write_structure(model, trimer(BINDER, rot=RIGID_R, trans=RIGID_T, atoms=atoms))
    m = measure(model, ref, fit_atoms="backbone").metrics
    check("backbone: 4 atoms per target residue", m["n_fit_atoms"], 4 * 90)
    check("backbone: 4 atoms per binder residue", m["n_binder_atoms"], 4 * len(BINDER_SEQ))
    check("backbone: rmsd still ~0", m["rmsd"], 0.0, tol=PDB_TOL)
    check("backbone: echoed", m["fit_atoms"], "backbone")
    # CA-only files cannot serve a backbone run: that is an error, not a silent
    # fallback to the atoms that happen to be there.
    write_structure(model, trimer(BINDER, rot=RIGID_R, trans=RIGID_T))
    check_raises("backbone: a CA-only model errors rather than falling back",
                 lambda: measure(model, ref, fit_atoms="backbone"),
                 "lack one of the N/CA/C/O atoms", worker.MeasureError)

    # Sequence identity of the BINDER is reported, never required: the two
    # binders are the same backbone with different sequences.
    write_structure(ref, trimer(BINDER))
    write_structure(model, trimer(BINDER, binder_seq=BINDER_SEQ_REDESIGNED,
                                  rot=RIGID_R, trans=RIGID_T))
    m = measure(model, ref).metrics
    expected = sum(a == b for a, b in zip(BINDER_SEQ, BINDER_SEQ_REDESIGNED)) / len(BINDER_SEQ)
    check("binder identity: reported", m["binder_seq_match_frac"], expected, tol=1e-9)
    check_true("binder identity: well below 1.0 and NOT a gate",
               float(m["binder_seq_match_frac"]) < 0.8)
    check("binder identity: the measurement still runs", m["rmsd"], 0.0, tol=PDB_TOL)


# ---------------------------------------------------------------------------
# 10. Manifest builder
# ---------------------------------------------------------------------------


def _run_dir(tmp: Path) -> tuple[Path, dict]:
    rd = tmp / "rundir"
    rd.mkdir()
    s = tmp / "structs"
    s.mkdir()
    paths = {}
    for k in ("p1", "p2", "p3", "c1", "c2", "c3", "c4", "own4"):
        paths[k] = s / f"{k}.pdb"
        paths[k].write_text("")
    t0 = pd.DataFrame(
        {"design_path": [str(paths["p1"]), str(paths["p2"]), str(paths["p3"])]},
        index=pd.Index(["p1", "p2", "p3"], name="name"),
    )
    t0.to_csv(rd / "table0.tsv", sep="\t")
    t1 = pd.DataFrame(
        {
            "boltz_path": [str(paths["c1"]), str(paths["c2"]), None, str(paths["c4"])],
            "design_path": [None, None, None, str(paths["own4"])],
            "parent_name": ["p1", "p2", "p1", "p3"],
            "parent_table": ["table0"] * 4,
            "gen": [1] * 4,
        },
        index=pd.Index(["c1", "c2", "c3", "c4"], name="name"),
    )
    # c2's parent p2 has a design_path; blank it to exercise missing_input.
    t1.to_csv(rd / "table1.tsv", sep="\t")
    return rd, paths


def _build(rd: Path, out_dir: Path, df=None, **overrides):
    dm = DataManager(rd)
    df = dm.read_frame("table1") if df is None else df
    fields = dict(
        gpus_per_task=1, col_a="boltz_path", col_b="design_path",
        binder_chain="D", target_chains="A:C",
        ref_binder_chain="D", ref_target_chains="A:C",
        fit_atoms="ca", target_match="number", min_target_res=10,
        min_target_res_frac=0.75, min_target_identity=0.8,
        fit_tolerance=1.0, max_target_fit_rmsd=5.0, pin_mapping=False,
        max_permutations=run_dockrmsd.DEFAULT_MAX_PERMUTATIONS,
        designs_per_task=25, output_prefix=None, force=False,
        input_column="not applicable", run_dir=rd, table="table1",
    )
    fields.update(overrides)
    args = Namespace(**fields)
    ctx = Namespace(df=df, args=args, out_dir=out_dir, lookup=partial(dm.lookup, df),
                    write_meta=lambda **k: None)
    return run_dockrmsd.build_dockrmsd_manifest(ctx), args  # type: ignore[arg-type]


def test_manifest(tmp: Path) -> None:
    rd, paths = _run_dir(tmp)
    out_dir = rd / "table1" / "dockrmsd"
    out_dir.mkdir(parents=True)
    prefix = "boltz_vs_design"
    res_dir = out_dir / prefix
    res_dir.mkdir()
    (res_dir / "c1.tsv").write_text("name\tstatus\tsup_path\nc1\tOK\t\n")  # stale

    rows, args = _build(rd, out_dir)
    check("manifest: gpus forced to 0", args.gpus_per_task, 0)
    check("manifest: one task for 3 designs", len(rows), 1)
    r = rows[0]
    check("manifest: fourteen fields", len(r), 14)
    check("manifest: chain fields", (r[1], r[2], r[3], r[4]), ("D", "A,B,C", "D", "A,B,C"))
    check("manifest: fit atoms", r[5], "ca")
    check("manifest: target match", r[6], "number")
    check("manifest: pin mapping off", r[12], "0")
    check("manifest: prefix last", r[13], prefix)
    rows_pin, _ = _build(rd, out_dir, pin_mapping=True, output_prefix="pinned")
    check("manifest: --pin-mapping reaches the task", rows_pin[0][12], "1")
    rows_seq, _ = _build(rd, out_dir, target_match="sequence", output_prefix="seqmode")
    check("manifest: --target-match reaches the task", rows_seq[0][6], "sequence")
    check_true("manifest: every field non-empty", all(f != "" for f in r))
    task_rows = [
        line.split("\t") for line in Path(r[0]).read_text().splitlines() if line
    ]
    by_name = {t[0]: t for t in task_rows}
    check("manifest: submitted designs", sorted(by_name), ["c1", "c2", "c4"])
    check("manifest: model from the row itself", by_name["c1"][1], str(paths["c1"].resolve()))
    check("manifest: reference resolved from the PARENT", by_name["c1"][2],
          str(paths["p1"].resolve()))
    check("manifest: a child's own col-b shadows the parent's",
          by_name["c4"][2], str(paths["own4"].resolve()))
    check_true("manifest: row without --col-a is not touched",
               not (res_dir / "c3.tsv").exists())
    check_true("manifest: stale result removed for a resubmitted design",
               not (res_dir / "c1.tsv").exists())

    # Batching.
    rows, _ = _build(rd, out_dir, designs_per_task=2)
    check("manifest: batching splits tasks", len(rows), 2)

    # Resume / --force.
    dm = DataManager(rd)
    df = dm.read_frame("table1")
    df[f"dockrmsd_{prefix}_status"] = ["OK", pd.NA, pd.NA, pd.NA]
    rows, _ = _build(rd, out_dir, df=df)
    names = [line.split("\t")[0] for line in Path(rows[0][0]).read_text().splitlines() if line]
    check("resume: OK row skipped", sorted(names), ["c2", "c4"])
    rows, _ = _build(rd, out_dir, df=df, force=True)
    names = [line.split("\t")[0] for line in Path(rows[0][0]).read_text().splitlines() if line]
    check("resume: --force resubmits", sorted(names), ["c1", "c2", "c4"])

    # A row whose --col-b resolves nowhere: missing_input, written here.
    t0 = dm.read_frame("table0")
    t0.loc["p2", "design_path"] = pd.NA
    t0.to_csv(rd / "table0.tsv", sep="\t")
    rows, _ = _build(rd, out_dir, output_prefix="nob")
    st = pd.read_csv(out_dir / "nob" / "c2.tsv", sep="\t").iloc[0]["status"]
    check_true("col-b: unresolved -> missing_input naming --col-b",
               st.startswith("missing_input:") and "--col-b 'design_path'" in st)
    names = [line.split("\t")[0] for line in Path(rows[0][0]).read_text().splitlines() if line]
    check("col-b: the rest are still submitted", sorted(names), ["c1", "c4"])

    refusals = [
        (dict(col_a="nope_path"), "not a column of the -t table"),
        (dict(col_b="nope_path"), "has no value on any"),
        (dict(col_a="boltz_path", col_b="boltz_path"), "compared with itself"),
        (dict(input_column="boltz_path"), "-i/--input-column"),
        (dict(binder_chain="D,E"), "exactly one chain"),
        (dict(ref_binder_chain=""), "at least one chain"),
        (dict(target_chains="A,B,D"), "is also in --target-chains"),
        (dict(ref_target_chains="A,B"), "same number of target chains"),
        (dict(target_chains="E:K", ref_target_chains="E:K"), "over --max-permutations"),
        (dict(min_target_res=2), "--min-target-res must be >= 3"),
        (dict(min_target_identity=1.5), "--min-target-identity must be in"),
        (dict(fit_tolerance=-1.0), "--fit-tolerance must be >= 0"),
        (dict(max_target_fit_rmsd=0.0), "--max-target-fit-rmsd must be > 0"),
        (dict(designs_per_task=0), "--designs-per-task must be >= 1"),
        (dict(min_target_res_frac=1.5), "--min-target-res-frac must be in"),
        (dict(output_prefix="bad prefix"), "letters, digits"),
    ]
    for over, needle in refusals:
        check_raises(f"refuse at submit: {over}", lambda o=over: _build(rd, out_dir, **o), needle)

    for a in (("x_path", "y_path", None), ("x", "y_path", "custom")):
        check(f"prefix: run == collect for {a}",
              run_dockrmsd.resolve_prefix(*a), collect_dockrmsd.resolve_prefix(*a))
    check("prefix: default", run_dockrmsd.resolve_prefix("boltz_path", "rfdiffusion3_path", None),
          "boltz_vs_rfdiffusion3")


# ---------------------------------------------------------------------------
# 11. Collector
# ---------------------------------------------------------------------------


def test_collector(tmp: Path) -> None:
    rd = tmp / "crun"
    out_dir = rd / "table1" / "dockrmsd"
    prefix = "boltz_vs_design"
    res_dir = out_dir / prefix
    res_dir.mkdir(parents=True)

    ref, model = tmp / "co_ref.pdb", tmp / "co_model.pdb"
    write_structure(ref, trimer(BINDER))
    write_structure(model, trimer(BINDER, rot=RIGID_R, trans=RIGID_T))
    worker.one_design("ok", str(model), str(ref), options(), res_dir)
    worker.one_design("warn", str(model), str(ref),
                      options(max_target_fit_rmsd=1e-12), res_dir)
    worker.one_design("err", str(model), str(ref), options(binder_chain="Z"), res_dir)
    # A binder-length error, to prove the trust columns survive it.
    bad = tmp / "co_short.pdb"
    write_structure(bad, trimer(BINDER[:-3], binder_seq=BINDER_SEQ[:-3]))
    worker.one_design("short", str(bad), str(ref), options(), res_dir)
    (res_dir / "mi.tsv").write_text(
        "name\tstatus\tsup_path\nmi\tmissing_input: --col-b 'design_path' has no value\t\n"
    )

    df = pd.DataFrame(index=pd.Index(["ok", "warn", "err", "short", "mi", "gone"], name="name"))
    ctx = Namespace(out_dir=out_dir, df=df,
                    args=Namespace(col_a="boltz_path", col_b="design_path",
                                   output_prefix=None))
    one = collect_dockrmsd.collect_dockrmsd(ctx)  # type: ignore[arg-type]

    def got(name):
        items = list(one(Namespace(name=name, out_dir=out_dir, lookup=None)))  # type: ignore[arg-type]
        check(f"collect {name}: exactly one Collected", len(items), 1)
        check(f"collect {name}: leaf status suppressed", items[0].status, None)
        return items[0].data

    k = lambda c: f"{prefix}_{c}"  # noqa: E731
    d = got("ok")
    check("collect ok: status", d[k("status")], "OK")
    check("collect ok: rmsd float", d[k("rmsd")], 0.0, tol=1e-3)
    check("collect ok: int column", d[k("n_permutations")], 6)
    check("collect ok: str column", d[k("permutation")], "A>A,B>B,C>C")
    check("collect ok: chain echo", d[k("binder_chain")], "D")
    check_true("collect ok: sup_path", bool(d[k("sup_path")]))
    check_true("collect ok: margin present", d[k("margin")] is not None)

    d = got("warn")
    check_true("collect warn: status kept", str(d[k("status")]).startswith("warn:"))
    check_true("collect warn: numbers kept", d[k("rmsd")] is not None)
    check_true("collect warn: sup_path kept", bool(d[k("sup_path")]))

    d = got("err")
    check_true("collect err: status", str(d[k("status")]).startswith("error:"))
    for c in ("rmsd", "max_dev", "complex_rmsd", "sup_path", "target_fit_rmsd"):
        check_true(f"collect err: {c} is NA, never 0", d[k(c)] is None)
    check("collect err: the chain echo is still there", d[k("binder_chain")], "Z")

    d = got("short")
    check_true("collect short: status", "binder length mismatch" in str(d[k("status")]))
    check("collect short: the two lengths are collected", d[k("binder_n_res")], 23)
    check("collect short: and the reference's", d[k("ref_binder_n_res")], 26)
    check_true("collect short: science NA", d[k("rmsd")] is None)

    d = got("mi")
    check_true("collect mi: missing_input kept",
               str(d[k("status")]).startswith("missing_input:"))
    nonnull = [c for c, v in d.items() if v is not None and c != k("status")]
    check("collect mi: every other column NA", nonnull, [])

    d = got("gone")
    check("collect gone: exactly 'missing'", d[k("status")], "missing")
    nonnull = [c for c, v in d.items() if v is not None and c != k("status")]
    check("collect gone: every other column NA", nonnull, [])

    # An OK row whose science is absent is an error, never a half-trusted OK.
    (res_dir / "liar.tsv").write_text("name\tstatus\tsup_path\nliar\tOK\t/x.cif\n")
    d = got("liar")
    check_true("collect: OK without science becomes an error",
               str(d[k("status")]).startswith("error: status 'OK'"))
    check_true("collect: and carries nothing else", d[k("rmsd")] is None)

    check("collector/worker column sets agree",
          set(collect_dockrmsd.METRIC_COLUMNS), set(worker.METRIC_COLUMNS))
    check("collector science set matches the worker's",
          collect_dockrmsd.SCIENCE_COLUMNS, worker.SCIENCE_COLUMNS)
    check_raises("collect: wrong prefix -> results dir not found",
                 lambda: collect_dockrmsd.collect_dockrmsd(Namespace(
                     out_dir=out_dir, df=df,
                     args=Namespace(col_a=None, col_b=None, output_prefix="other"))),  # type: ignore[arg-type]
                 "Results dir not found", FileNotFoundError)


# ---------------------------------------------------------------------------
# 12. The task script, through the real prosapia prelude
# ---------------------------------------------------------------------------


def test_task_script(tmp: Path) -> None:
    import prosapia.core as core

    prelude = Path(core.__file__).parent / "scripts" / "sapia_task_prelude.sh"
    work = tmp / "task"
    work.mkdir()
    ref, model = work / "t_ref.pdb", work / "t_model.pdb"
    write_structure(ref, trimer(BINDER))
    write_structure(model, trimer(BINDER @ R120.T, rot=RIGID_R, trans=RIGID_T))
    out_dir = work / "out"
    out_dir.mkdir()
    task_file = work / "task_0.tsv"
    task_file.write_text(f"t1\t{model}\t{ref}\nt2\t{model}\t{work / 'nope.pdb'}\n")
    manifest = work / "manifest.txt"
    manifest.write_text(
        "\t".join([str(task_file), "D", "A,B,C", "D", "A,B,C", "ca", "number",
                   "10", "0.75", "0.8", "1.0", "5.0", "0", "pfx"]) + "\n"
    )
    env = dict(os.environ, SAPIA_PRELUDE=str(prelude), SAPIA_SCHEDULER="modal",
               SAPIA_TASK_ID="1", SAPIA_TOOL_DIR=str(HERE), PIPELINE_PYTHON=sys.executable)
    proc = subprocess.run(["bash", str(HERE / "dockrmsd.sh"), str(manifest), str(out_dir)],
                          cwd=work, env=env, capture_output=True, text=True)
    check("task: exit 0", proc.returncode, 0)
    if proc.returncode != 0:
        FAILURES.append(f"task stderr: {proc.stderr[-600:]}")
        return
    read = lambda n: pd.read_csv(out_dir / "pfx" / f"{n}.tsv", sep="\t", dtype=str,  # noqa: E731
                                 keep_default_na=False, na_values=[""]).iloc[0]
    r1 = read("t1")
    check("task t1: status OK (every manifest field reached the worker)", r1["status"], "OK")
    check("task t1: the symmetry copy was found", float(r1["rmsd"]), 0.0, tol=1e-3)
    check("task t1: permutation", r1["permutation"], "B>A,C>B,A>C")
    check("task t1: fit atoms field", r1["fit_atoms"], "ca")
    check("task t1: target match field", r1["target_match"], "number")
    check("task t1: pin mapping field", int(r1["pin_mapping"]), 0)
    check_true("task t1: sup written", (out_dir / "pfx" / "t1_sup.cif").is_file())
    r2 = read("t2")
    check_true("task t2: a missing reference is an error row, not a crash",
               r2["status"].startswith("error:") and "not found" in r2["status"])
    check_true("task t2: science NA, never 0", pd.isna(r2["rmsd"]))


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        test_geometry_and_atoms(tmp)
        test_known_good(tmp)
        test_binder_not_fitted(tmp)
        test_symmetry_search(tmp)
        test_binder_length_enforced(tmp)
        test_role_flip(tmp)
        test_author_numbering(tmp)
        test_inserted_loop(tmp)
        test_partial_misregistration(tmp)
        test_per_chain_pair_floor(tmp)
        test_pin_mapping(tmp)
        test_identity_and_warn(tmp)
        test_manifest(tmp)
        test_wrong_chains(tmp)
        test_collector(tmp)
        test_task_script(tmp)

    if FAILURES:
        print(f"FAILED {len(FAILURES)} of {CHECKS} checks:")
        for f in FAILURES:
            print(f"  - {f}")
        return 1
    print(f"OK: {CHECKS} checks passed.")
    print(
        "NOTE: this proves the arithmetic and the contracts only. It CANNOT prove "
        "that --target-match sequence reproduces the ten hand-measured values on "
        "real files, nor that --fit-tolerance 1.0 separates symmetry copies on a "
        "real near-symmetric target -- see the module docstring."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
