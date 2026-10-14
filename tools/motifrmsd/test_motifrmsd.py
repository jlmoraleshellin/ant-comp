#!/usr/bin/env python3
"""
Self-test for motifrmsd: synthetic structures with a hand-computable answer.

WHAT THIS PROVES, AND WHAT IT CANNOT.

It proves the *geometry and the contracts*: that a rigid-body-moved design gives
RMSD 0, that a known per-residue displacement gives the RMSD you can compute by
hand, that a reflected ("mirror image") design is NOT fitted onto its reference,
that the segment grouping finds the runs, that per-segment RMSDs come from ONE
global fit rather than per-segment refits, and that every refusal path writes a
non-OK status with NA rather than a plausible number.

It CANNOT prove the premise or the mapping. The structures here are written by
this file, so they agree with this file by construction. In particular these
remain UNKNOWN until the tool is run on real rfd3 output:

  * whether real rfd3 motif-scaffolding output names the design chain the way
    diffused_index_map's values say it does;
  * whether rfd3 writes the motif's real residue identities into the output
    backbone or a poly-glycine, i.e. what seq_match_frac is *expected* to be;
  * whether every design in a real table resolves a sidecar through the
    .cif_to_pdb/*.src route (verified by directory listing for one run_dir, not
    by running the tool);
  * whether `specification.input` is readable from inside a task container, i.e.
    whether ref_matches_spec comes back 1 or NA in practice;
  * for --motif-chain: whether the letter a caller passes is really the one rfd3
    writes on the KEY side of diffused_index_map for the grafted chain in a given
    campaign's contig. The fixtures here spell the keys themselves, so they prove
    the filtering arithmetic and not the chain convention of real output.

Run:  uv run python tools/motifrmsd/test_motifrmsd.py
"""

import importlib.util
import json
import sys
import tempfile
from argparse import Namespace
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


worker = _load("motifrmsd_worker")
run_motifrmsd = _load("run_motifrmsd")
collect_motifrmsd = _load("collect_motifrmsd")

FAILURES: list[str] = []
CHECKS = 0


def check(label: str, got, want, tol: float | None = None) -> None:
    global CHECKS
    CHECKS += 1
    ok = (
        abs(float(got) - float(want)) <= tol
        if tol is not None
        else got == want
    )
    if not ok:
        FAILURES.append(f"{label}: got {got!r}, want {want!r}")


def check_true(label: str, value: bool) -> None:
    global CHECKS
    CHECKS += 1
    if not value:
        FAILURES.append(f"{label}: expected True")


# ---------------------------------------------------------------------------
# Synthetic structures
# ---------------------------------------------------------------------------

AA3 = {"A": "ALA", "G": "GLY", "S": "SER", "V": "VAL", "L": "LEU", "W": "TRP"}


def write_structure(
    path: Path, chains: dict[str, list[tuple[int, str, np.ndarray]]]
) -> None:
    """Write a PDB of CA-only residues: {chain: [(resnum, one_letter, xyz), ...]}."""
    st = gemmi.Structure()
    model = gemmi.Model("1")
    for chain_name, residues in chains.items():
        chain = gemmi.Chain(chain_name)
        for resnum, code, xyz in residues:
            res = gemmi.Residue()
            res.name = AA3[code]
            res.seqid = gemmi.SeqId(resnum, " ")
            atom = gemmi.Atom()
            atom.name = "CA"
            atom.element = gemmi.Element("C")
            atom.pos = gemmi.Position(*map(float, xyz))
            res.add_atom(atom)
            chain.add_residue(res)
        model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()
    path.parent.mkdir(parents=True, exist_ok=True)
    st.write_pdb(str(path))


def helix(n: int, start: int, offset: np.ndarray) -> list[np.ndarray]:
    """n points on a gentle helix -- non-collinear and non-degenerate, so Kabsch
    is well conditioned and a reflection is genuinely distinguishable."""
    out = []
    for i in range(n):
        t = (start + i) * 0.6
        out.append(offset + np.array([np.cos(t) * 5.0, np.sin(t) * 5.0, (start + i) * 1.5]))
    return out


def rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    k = np.array(
        [
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0],
        ]
    )
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * (k @ k)


# Reference: chain D, two motif segments 24-31 and 48-55, plus non-motif padding.
SEG_A = list(range(24, 32))
SEG_B = list(range(48, 56))
REF_CODES = {n: "AGSVLWAG"[i % 8] for i, n in enumerate(SEG_A + SEG_B)}

REF_COORDS: dict[int, np.ndarray] = {}
for i, n in enumerate(SEG_A):
    REF_COORDS[n] = helix(8, 0, np.zeros(3))[i]
for i, n in enumerate(SEG_B):
    REF_COORDS[n] = helix(8, 20, np.array([18.0, -7.0, 3.0]))[i]

# diffused_index_map: reference D24..D31 -> design A1..A8, D48..D55 -> A30..A37
# (rfd3 renumbers the output chain from 1, so the design numbering is unrelated
# to the reference's -- which is the whole reason the map has to be read.)
INDEX_MAP = {}
for i, n in enumerate(SEG_A):
    INDEX_MAP[f"D{n}"] = f"A{1 + i}"
for i, n in enumerate(SEG_B):
    INDEX_MAP[f"D{n}"] = f"A{30 + i}"

REF_ORDER = SEG_A + SEG_B  # reference (chain, resnum) sort order


def build_ref(path: Path) -> None:
    write_structure(
        path,
        {
            "D": [(n, REF_CODES[n], REF_COORDS[n]) for n in REF_ORDER],
            # A decoy chain, to prove the map's chain id is honoured.
            "E": [(n, "G", REF_COORDS[n] + np.array([100.0, 0, 0])) for n in REF_ORDER],
        },
    )


def build_design(
    path: Path,
    transform=None,
    perturb: dict[int, np.ndarray] | None = None,
    codes: dict[int, str] | None = None,
    drop: set[int] | None = None,
) -> None:
    """A design chain A whose motif residues are the reference's, optionally moved
    rigidly, perturbed per residue (keyed by REFERENCE resnum), relabelled, or
    with some residues removed. Plus scaffold residues the map never mentions."""
    drop = drop or set()
    perturb = perturb or {}
    codes = codes or {}
    residues: list[tuple[int, str, np.ndarray]] = []
    for ref_num in REF_ORDER:
        if ref_num in drop:
            continue
        design_num = int(INDEX_MAP[f"D{ref_num}"][1:])
        xyz = REF_COORDS[ref_num] + perturb.get(ref_num, np.zeros(3))
        if transform is not None:
            xyz = transform(xyz)
        residues.append((design_num, codes.get(ref_num, REF_CODES[ref_num]), xyz))
    # Scaffold: residues the motif map never names, far away. They must not enter
    # the fit at all.
    for k in range(10, 20):
        residues.append((k, "G", np.array([200.0 + k, 0.0, 0.0])))
    residues.sort()
    write_structure(path, {"A": residues})


def write_sidecar(
    path: Path,
    ref_path: Path,
    index_map: dict | None = None,
    contig: str = "10-25,D24-31,10-30,D48-55,10-30",
) -> None:
    data = {
        "metrics": {"max_ca_deviation": 99.0, "n_chainbreaks": 0},
        "specification": {
            "input": str(ref_path),
            "contig": contig,
            "extra": {"sampled_contig": "22P,D24,D31,28P", "example_id": "x"},
        },
        "ckpt_path": "/somewhere/rfd3_latest.ckpt",
        "seed": None,
    }
    if index_map is not None:
        data["diffused_index_map"] = index_map
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


def run_worker(tmp: Path, name: str, design: Path, sidecar: Path, ref: str = "",
               allow_partial: bool = False, motif_chain: str = ""):
    args = Namespace(
        name=name,
        design=design,
        sidecar=sidecar,
        ref=ref,
        allow_partial=allow_partial,
        motif_chain=motif_chain,
        out_pdb=tmp / f"{name}_onmotif.pdb",
        result_tsv=tmp / f"{name}.tsv",
    )
    try:
        result = worker.compute(args)
        assert result.structure is not None
        result.structure.setup_entities()
        result.structure.write_pdb(str(args.out_pdb))
        worker.write_result(
            args.result_tsv,
            name,
            result.status,
            str(args.out_pdb),
            result.metrics,
            worker.seg_columns(result.seg_rmsds),
        )
    except Exception as e:  # noqa: BLE001
        worker.write_result(args.result_tsv, name, f"error: {e}", "", {})
    return pd.read_csv(args.result_tsv, sep="\t").iloc[0]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_kabsch_identity_and_reflection() -> None:
    """A rigidly moved point set fits back exactly; a reflected one does not."""
    rng = np.random.default_rng(0)
    p = rng.normal(size=(12, 3)) * 8.0
    rot = rotation(np.array([0.3, -0.7, 0.5]), 1.1)
    trans = np.array([4.0, -9.0, 2.5])
    q = p @ rot.T + trans
    r, t = worker.kabsch(p, q)
    check("kabsch: rot recovered", float(np.abs(r - rot).max()), 0.0, tol=1e-9)
    check("kabsch: trans recovered", float(np.abs(t - trans).max()), 0.0, tol=1e-9)
    check("kabsch: rmsd of exact fit", worker.rmsd(p @ r.T + t, q), 0.0, tol=1e-9)

    # Determinant +1 always: a mirror image must not be flattered into a fit.
    mirrored = p * np.array([1.0, 1.0, -1.0])
    r2, t2 = worker.kabsch(mirrored, p)
    check("kabsch: no reflection", float(np.linalg.det(r2)), 1.0, tol=1e-9)
    check_true(
        "kabsch: mirrored set does not fit",
        worker.rmsd(mirrored @ r2.T + t2, p) > 1.0,
    )

    # Independent cross-check on a genuinely imperfect fit: the closed-form
    # residual from the singular values, which never calls worker.kabsch.
    noisy = q + rng.normal(size=q.shape) * 0.7
    rn, tn = worker.kabsch(p, noisy)
    pc, qc = p - p.mean(0), noisy - noisy.mean(0)
    sv = np.linalg.svd(pc.T @ qc, compute_uv=False)
    residual = (pc**2).sum() + (qc**2).sum() - 2 * sv.sum()
    check(
        "kabsch: rmsd matches singular-value residual",
        worker.rmsd(p @ rn.T + tn, noisy),
        float(np.sqrt(max(residual, 0.0) / len(p))),
        tol=1e-9,
    )
    check_true(
        "kabsch: noisy fit is a real residual",
        worker.rmsd(p @ rn.T + tn, noisy) > 0.3,
    )

    try:
        worker.kabsch(p[:2], q[:2])
        FAILURES.append("kabsch: 2 points should raise")
    except ValueError:
        pass
    global CHECKS
    CHECKS += 1


def test_parse_and_segments() -> None:
    check("parse D24", worker.parse_res_key("D24", "x"), ("D", 24))
    check("parse A1", worker.parse_res_key("A1", "x"), ("A", 1))
    for bad in ("D24A", "24", "", "1D", "D 24"):
        try:
            worker.parse_res_key(bad, "x")
            FAILURES.append(f"parse: {bad!r} should raise")
        except ValueError:
            pass
        global CHECKS
        CHECKS += 1

    keys = [("D", 25), ("D", 24), ("D", 48), ("D", 26), ("D", 49)]
    # Design side contiguous wherever the reference side is (D24-26 -> A1-3,
    # D48-49 -> A30-31), so only the reference gaps can split runs here.
    dkeys = [("A", 2), ("A", 1), ("A", 30), ("A", 3), ("A", 31)]
    runs = worker.segments(keys, dkeys)
    check("segments: count", len(runs), 2)
    check("segments: first run is 24,25,26", [keys[i][1] for i in runs[0]], [24, 25, 26])
    check("segments: second run is 48,49", [keys[i][1] for i in runs[1]], [48, 49])
    # A chain change breaks a run even when the numbers are consecutive.
    check(
        "segments: chain break",
        len(worker.segments([("D", 24), ("E", 25)], [("A", 1), ("A", 2)])),
        2,
    )

    # --- design-side contiguity (fix 2, 2026-10-08) ---
    # Reference-contiguous D26..D29, but the design has a gap after A2: two runs.
    ref4 = [("D", 26), ("D", 27), ("D", 28), ("D", 29)]
    check(
        "segments: design-side gap splits a reference-contiguous run",
        [[ref4[i][1] for i in r] for r in worker.segments(ref4, [("A", 1), ("A", 2), ("A", 9), ("A", 10)])],
        [[26, 27], [28, 29]],
    )
    # ... a design-side chain change does too ...
    check(
        "segments: design-side chain change splits",
        len(worker.segments(ref4, [("A", 1), ("A", 2), ("B", 3), ("B", 4)])),
        2,
    )
    # ... a design side running BACKWARDS is not contiguous ...
    check(
        "segments: design side reversed is not contiguous",
        len(worker.segments(ref4, [("A", 4), ("A", 3), ("A", 2), ("A", 1)])),
        4,
    )
    # ... and when both sides are contiguous nothing changes: one run.
    check(
        "segments: both sides contiguous stays one run",
        len(worker.segments(ref4, [("A", 5), ("A", 6), ("A", 7), ("A", 8)])),
        1,
    )
    try:
        worker.segments(ref4, [("A", 1)])
        FAILURES.append("segments: misaligned key lists should raise")
    except ValueError:
        pass
    CHECKS += 1


def test_perfect_graft(tmp: Path) -> None:
    """A design that is the reference motif rigidly moved: RMSD 0 everywhere."""
    ref = tmp / "ref.pdb"
    build_ref(ref)
    rot = rotation(np.array([1.0, 0.2, -0.3]), 0.9)
    trans = np.array([30.0, -12.0, 7.0])
    design = tmp / "perfect.pdb"
    build_design(design, transform=lambda v: rot @ v + trans)
    sidecar = tmp / "perfect.json"
    write_sidecar(sidecar, ref, INDEX_MAP)

    row = run_worker(tmp, "perfect", design, sidecar)
    check("perfect: status", str(row["status"]), "OK")
    check("perfect: n_expected", int(row["n_expected"]), 16)
    check("perfect: n_matched", int(row["n_matched"]), 16)
    check("perfect: rmsd", float(row["rmsd"]), 0.0, tol=1e-3)
    check("perfect: max_dev", float(row["max_dev"]), 0.0, tol=1e-3)
    check("perfect: n_segments", int(row["n_segments"]), 2)
    check("perfect: seg_ranges", str(row["seg_ranges"]), "D24-31/D48-55")
    check("perfect: seq_match_frac", float(row["seq_match_frac"]), 1.0, tol=1e-9)
    check("perfect: ref_chains", str(row["ref_chains"]), "D")
    check("perfect: design_chains", str(row["design_chains"]), "A")
    check("perfect: ref_source", str(row["ref_source"]), "sidecar")
    check("perfect: ref_matches_spec", int(row["ref_matches_spec"]), 1)
    check("perfect: seg0_rmsd", float(row["seg0_rmsd"]), 0.0, tol=1e-3)
    check("perfect: seg1_rmsd", float(row["seg1_rmsd"]), 0.0, tol=1e-3)
    check("perfect: sampled_contig", str(row["sampled_contig"]), "22P,D24,D31,28P")
    # The scaffold residues (never named by the map) must not have entered the fit.
    check("perfect: fit used only the motif", int(row["n_matched"]), 16)


def test_known_displacement(tmp: Path) -> None:
    """Every motif CA pushed 1 A along +z: the fit cannot absorb a shear, so the
    RMSD is a number written down before running anything.

    Residues 24-31 move +1 z, residues 48-55 move -1 z. The best rigid fit of
    that pattern is not the identity, so the RMSD is bounded above by 1.0 (the
    unfitted RMSD) and strictly positive. The hand-checkable fact is the simpler
    one below it: a UNIFORM shift is absorbed entirely.
    """
    ref = tmp / "ref.pdb"
    build_ref(ref)
    sidecar = tmp / "shift.json"
    write_sidecar(sidecar, ref, INDEX_MAP)

    # 1. A uniform shift of every motif atom is pure translation -> RMSD 0.
    design = tmp / "uniform.pdb"
    build_design(design, perturb={n: np.array([0.0, 0.0, 2.0]) for n in REF_ORDER})
    row = run_worker(tmp, "uniform", design, sidecar)
    check("uniform shift: status", str(row["status"]), "OK")
    check("uniform shift: rmsd absorbed", float(row["rmsd"]), 0.0, tol=1e-3)

    # 2. One residue displaced 3 A, the other 15 untouched. The fit redistributes
    #    a little, so RMSD <= 3/sqrt(16) = 0.75 and max_dev <= 3.0, both > 0.
    design = tmp / "one_off.pdb"
    build_design(design, perturb={SEG_A[3]: np.array([3.0, 0.0, 0.0])})
    row = run_worker(tmp, "one_off", design, sidecar)
    check("one_off: status", str(row["status"]), "OK")
    check_true("one_off: rmsd positive", float(row["rmsd"]) > 0.05)
    check_true("one_off: rmsd <= unfitted bound", float(row["rmsd"]) <= 0.75 + 1e-6)
    check_true("one_off: max_dev > rmsd", float(row["max_dev"]) > float(row["rmsd"]))
    check_true("one_off: max_dev <= 3", float(row["max_dev"]) <= 3.0 + 1e-6)


def test_global_fit_is_shared_across_segments(tmp: Path) -> None:
    """The headline contract: segment 1 is perfect in its own frame but rotated
    as a rigid body relative to segment 0. A per-segment refit would report ~0
    for BOTH segments and hide the dislocation; one global fit must not."""
    ref = tmp / "ref.pdb"
    build_ref(ref)
    sidecar = tmp / "dislocated.json"
    write_sidecar(sidecar, ref, INDEX_MAP)

    rot = rotation(np.array([0.0, 0.0, 1.0]), 0.35)
    centre = np.mean([REF_COORDS[n] for n in SEG_B], axis=0)
    perturb = {
        n: (rot @ (REF_COORDS[n] - centre) + centre) - REF_COORDS[n] for n in SEG_B
    }
    design = tmp / "dislocated.pdb"
    build_design(design, perturb=perturb)
    row = run_worker(tmp, "dislocated", design, sidecar)

    check("dislocated: status", str(row["status"]), "OK")
    seg0, seg1 = float(row["seg0_rmsd"]), float(row["seg1_rmsd"])
    check_true("dislocated: segment 0 is disturbed by the global fit", seg0 > 0.1)
    check_true("dislocated: segment 1 is disturbed by the global fit", seg1 > 0.1)
    check(
        "dislocated: worst_seg_rmsd is the max",
        float(row["worst_seg_rmsd"]),
        max(seg0, seg1),
        tol=1e-3,
    )
    check_true("dislocated: global rmsd positive", float(row["rmsd"]) > 0.1)

    # Prove the counterfactual explicitly: refitting segment 1 alone gives ~0,
    # which is the number a per-segment refit would have published.
    seg_b_ref = np.array([REF_COORDS[n] for n in SEG_B])
    seg_b_design = np.array([REF_COORDS[n] + perturb[n] for n in SEG_B])
    r, t = worker.kabsch(seg_b_design, seg_b_ref)
    check(
        "dislocated: a per-segment refit would have hidden it",
        worker.rmsd(seg_b_design @ r.T + t, seg_b_ref),
        0.0,
        tol=1e-6,
    )
    check_true(
        "dislocated: the shared fit does not hide it", seg1 > 0.5
    )


def test_error_contracts(tmp: Path) -> None:
    """Every refusal: non-OK status, and the collector turns it into NA."""
    ref = tmp / "ref.pdb"
    build_ref(ref)
    design = tmp / "ok.pdb"
    build_design(design)

    # 1. Sidecar missing entirely.
    row = run_worker(tmp, "nosidecar", design, tmp / "absent.json")
    check_true("no sidecar: error", str(row["status"]).startswith("error:"))

    # 2. Sidecar present but no diffused_index_map (the rc-foundry shape that
    #    this campaign actually found for ca_rmsd_to_input).
    sc = tmp / "nomap.json"
    write_sidecar(sc, ref, index_map=None)
    row = run_worker(tmp, "nomap", design, sc)
    check_true("no index map: error", str(row["status"]).startswith("error:"))
    check_true("no index map: names the key", "diffused_index_map" in str(row["status"]))

    # 3. A mapped residue absent from the DESIGN -> error by default, NA numbers.
    sc = tmp / "full.json"
    write_sidecar(sc, ref, INDEX_MAP)
    partial = tmp / "partial.pdb"
    build_design(partial, drop={SEG_A[2]})
    row = run_worker(tmp, "partial", partial, sc)
    check_true("partial: error by default", str(row["status"]).startswith("error:"))
    check_true("partial: counts are in the message", "15 of 16" in str(row["status"]))
    check_true("partial: rmsd is NA not 0", pd.isna(row.get("rmsd")))

    # 3b. ... and --allow-partial downgrades it to warn WITH numbers.
    row = run_worker(tmp, "partial_ok", partial, sc, allow_partial=True)
    check_true("partial+flag: warn", str(row["status"]).startswith("warn:"))
    check("partial+flag: n_expected", int(row["n_expected"]), 16)
    check("partial+flag: n_matched", int(row["n_matched"]), 15)
    check("partial+flag: rmsd", float(row["rmsd"]), 0.0, tol=1e-3)

    # 4. A reference that does not carry the mapped chain at all.
    other = tmp / "otherchain.pdb"
    write_structure(
        other, {"Z": [(n, REF_CODES[n], REF_COORDS[n]) for n in REF_ORDER]}
    )
    sc = tmp / "otherref.json"
    write_sidecar(sc, other, INDEX_MAP)
    row = run_worker(tmp, "wrongchain", design, sc)
    check_true("wrong ref chain: error", str(row["status"]).startswith("error:"))
    check_true("wrong ref chain: rmsd NA", pd.isna(row.get("rmsd")))

    # 5. Too few pairs for a fit.
    tiny_map = {f"D{SEG_A[0]}": INDEX_MAP[f"D{SEG_A[0]}"]}
    sc = tmp / "tiny.json"
    write_sidecar(sc, ref, tiny_map)
    row = run_worker(tmp, "tiny", design, sc)
    check_true("one pair: error", str(row["status"]).startswith("error:"))

    # 6. An unparseable residue key.
    sc = tmp / "badkey.json"
    write_sidecar(sc, ref, {**INDEX_MAP, "D24A": "A1"})
    row = run_worker(tmp, "badkey", design, sc)
    check_true("bad key: error", str(row["status"]).startswith("error:"))

    # 7. A reference file that is not there.
    sc = tmp / "noref.json"
    write_sidecar(sc, tmp / "does_not_exist.pdb", INDEX_MAP)
    row = run_worker(tmp, "noref", design, sc)
    check_true("missing ref: error", str(row["status"]).startswith("error:"))


def test_map_not_one_to_one(tmp: Path) -> None:
    """Fix 1 (2026-10-08): two reference residues mapped to the SAME design
    residue. Every count check passes (n_matched == n_expected), so without the
    one-to-one guard one design CA would silently enter the fit twice."""
    ref = tmp / "ref.pdb"
    build_ref(ref)
    design = tmp / "ok.pdb"
    build_design(design)

    dup_map = dict(INDEX_MAP)
    dup_map[f"D{SEG_A[1]}"] = INDEX_MAP[f"D{SEG_A[0]}"]  # D25 -> A1, same as D24
    sc = tmp / "dupvalue.json"
    write_sidecar(sc, ref, dup_map)

    row = run_worker(tmp, "dupvalue", design, sc)
    status = str(row["status"])
    check_true("dup value: error, not a number", status.startswith("error:"))
    check_true("dup value: says not one-to-one", "not one-to-one" in status)
    check_true("dup value: names the duplicated design key", "A1 <- D24, D25" in status)
    check_true("dup value: rmsd NA not 0", pd.isna(row.get("rmsd")))
    check_true("dup value: n_matched NA", pd.isna(row.get("n_matched")))
    check_true("dup value: worst_seg_rmsd NA", pd.isna(row.get("worst_seg_rmsd")))

    # --allow-partial is about MISSING residues; it must not launder a broken map.
    row = run_worker(tmp, "dupvalue_partial", design, sc, allow_partial=True)
    check_true(
        "dup value + --allow-partial: still error",
        str(row["status"]).startswith("error:"),
    )

    # Two key spellings that parse to the same reference residue ("D24"/"D024"),
    # each sent to a different design residue: the same CA-used-twice failure on
    # the reference side.
    dup_key_map = dict(INDEX_MAP)
    dup_key_map["D024"] = "A15"  # A15: present in the design (scaffold), unmapped
    sc = tmp / "dupkey.json"
    write_sidecar(sc, ref, dup_key_map)
    row = run_worker(tmp, "dupkey", design, sc)
    status = str(row["status"])
    check_true("dup parsed ref key: error", status.startswith("error:"))
    check_true("dup parsed ref key: names it", "D24 <- D24, D024" in status)
    check_true("dup parsed ref key: rmsd NA", pd.isna(row.get("rmsd")))


def test_reference_adjacent_segments_split_by_design(tmp: Path) -> None:
    """Fix 2 (2026-10-08): two motif pieces numbered consecutively in the
    reference (D26-32, D33-40) but placed apart in the design (A1-7, A20-27,
    scaffold between). They are two loops: n_segments 2, and worst_seg_rmsd
    isolates the bad one instead of blending it with the good one."""
    seg_p = list(range(26, 33))  # D26..D32, 7 residues
    seg_q = list(range(33, 41))  # D33..D40, 8 residues
    ref_order = seg_p + seg_q
    coords = dict(zip(ref_order, helix(len(ref_order), 0, np.zeros(3))))
    codes = {n: "AGSVLW"[i % 6] for i, n in enumerate(ref_order)}
    idx_map = {f"D{n}": f"A{1 + i}" for i, n in enumerate(seg_p)}
    idx_map.update({f"D{n}": f"A{20 + i}" for i, n in enumerate(seg_q)})

    ref = tmp / "adj_ref.pdb"
    write_structure(ref, {"D": [(n, codes[n], coords[n]) for n in ref_order]})

    # Segment P exact; segment Q distorted non-rigidly (alternating +/-1.5 A
    # across x), so no rigid fit can absorb it.
    perturb = {
        n: np.array([1.5 if i % 2 == 0 else -1.5, 0.0, 0.0])
        for i, n in enumerate(seg_q)
    }
    residues = [
        (int(idx_map[f"D{n}"][1:]), codes[n], coords[n] + perturb.get(n, np.zeros(3)))
        for n in ref_order
    ]
    residues += [(k, "G", np.array([200.0 + k, 0.0, 0.0])) for k in range(8, 20)]
    residues.sort(key=lambda r: r[0])
    design = tmp / "adj_design.pdb"
    write_structure(design, {"A": residues})

    sc = tmp / "adj.json"
    write_sidecar(sc, ref, idx_map, contig="5-10,D26-32,12,D33-40,5-10")
    row = run_worker(tmp, "adjacent", design, sc)

    check("adjacent: status", str(row["status"]), "OK")
    check("adjacent: n_matched", int(row["n_matched"]), 15)
    check("adjacent: two segments, not one", int(row["n_segments"]), 2)
    check("adjacent: seg_ranges", str(row["seg_ranges"]), "D26-32/D33-40")
    seg0, seg1 = float(row["seg0_rmsd"]), float(row["seg1_rmsd"])
    check_true("adjacent: the distorted segment is the worse one", seg1 > seg0)
    check(
        "adjacent: worst_seg_rmsd is the worse segment",
        float(row["worst_seg_rmsd"]),
        seg1,
        tol=1e-3,
    )
    # The blended answer the old (reference-only) grouping would have given: one
    # segment, whose RMSD IS the global RMSD. worst_seg_rmsd must exceed it.
    check_true(
        "adjacent: worst_seg_rmsd is not the blended global rmsd",
        float(row["worst_seg_rmsd"]) > float(row["rmsd"]) + 0.05,
    )

    # Control: the same reference ranges placed CONTIGUOUSLY in the design
    # (D33 -> A8) really are one loop, and stay one segment.
    joined_map = {f"D{n}": f"A{1 + i}" for i, n in enumerate(ref_order)}
    joined = tmp / "adj_joined.pdb"
    write_structure(
        joined,
        {"A": [(1 + i, codes[n], coords[n]) for i, n in enumerate(ref_order)]},
    )
    sc2 = tmp / "adj_joined.json"
    write_sidecar(sc2, ref, joined_map, contig="5-10,D26-40,5-10")
    row = run_worker(tmp, "adjacent_joined", joined, sc2)
    check("adjacent joined: status", str(row["status"]), "OK")
    check("adjacent joined: one segment", int(row["n_segments"]), 1)
    check("adjacent joined: seg_ranges", str(row["seg_ranges"]), "D26-40")


def _contiguous_design(
    path: Path, index_map: dict[str, str], chain_lengths: dict[str, int]
) -> None:
    """A design whose chains are numbered 1..N with NO gaps -- i.e. what rfd3 is
    understood to write, and the only shape in which residue numbers ARE
    ProteinMPNN's 1-indexed chain positions.

    Motif residues take their reference coordinates (so the graft is perfect and
    the RMSD question is out of the way); every other position is filler placed
    far off, where it cannot enter a fit it was never mapped into.
    """
    placed = {
        (v[0], int(v[1:])): REF_COORDS[int(k[1:])] for k, v in index_map.items()
    }
    chains: dict[str, list[tuple[int, str, np.ndarray]]] = {}
    for chain_name, length in chain_lengths.items():
        residues = []
        for num in range(1, length + 1):
            xyz = placed.get(
                (chain_name, num),
                np.array([300.0 + num, 50.0 * (ord(chain_name) - 65), 0.0]),
            )
            residues.append((num, "G", xyz))
        chains[chain_name] = residues
    write_structure(path, chains)


def test_design_fixed_positions(tmp: Path) -> None:
    """design_fixed_positions: the SAME segments as seg_ranges, printed from the
    design side, in ProteinMPNN's own single-chain-group syntax -- and NA, never
    a partial string, wherever it cannot be made honestly."""
    ref = tmp / "fp_ref.pdb"
    build_ref(ref)

    # D24-31 -> A5-12 and D48-55 -> A30-37. Deliberately NOT the reference
    # numbers, and deliberately not starting at 1, so a value that accidentally
    # printed the reference side (or an index) cannot pass.
    idx_map = {f"D{n}": f"A{5 + i}" for i, n in enumerate(SEG_A)}
    idx_map.update({f"D{n}": f"A{30 + i}" for i, n in enumerate(SEG_B)})

    design = tmp / "fp_design.pdb"
    _contiguous_design(design, idx_map, {"A": 40})
    sc = tmp / "fp.json"
    write_sidecar(sc, ref, idx_map)
    row = run_worker(tmp, "fp_ok", design, sc)

    check("fixed positions: status", str(row["status"]), "OK")
    check("fixed positions: two segments", int(row["n_segments"]), 2)
    check("fixed positions: seg_ranges unchanged (reference side)",
          str(row["seg_ranges"]), "D24-31/D48-55")
    value = str(row["design_fixed_positions"])
    check("fixed positions: design-side ranges", value, "5:12,30:37")

    # The three ways this column could be quietly wrong.
    check_true("fixed positions: no trailing slash", not value.endswith("/"))
    check_true("fixed positions: no chain separator at all", "/" not in value)
    check_true(
        "fixed positions: not the reference numbering",
        "24:31" not in value and "48:55" not in value,
    )
    check(
        "fixed positions: fragment count == segment count",
        len(value.split(",")),
        int(row["n_segments"]),
    )

    # Round-trip through prosapia's OWN position parser -- the code ProteinMPNN's
    # --fixed-positions actually goes through -- rather than through this file's
    # idea of the syntax. One chain group, holding exactly the mapped design
    # residues and nothing else.
    from prosapia.utils import parse_positions

    groups = parse_positions(value, lambda *a, **k: None, "fp_ok")
    check("fixed positions: parses to ONE chain group", len(groups), 1)
    expected = sorted(int(v[1:]) for v in idx_map.values())
    check("fixed positions: covers exactly the mapped design residues",
          sorted(groups[0]), expected)

    # ... and a trailing slash really would declare a second, empty group, which
    # is why this column must not carry one for a single --chains-to-design chain.
    check("fixed positions: a trailing slash would add a group",
          len(parse_positions(value + "/", lambda *a, **k: None, "fp_ok")), 2)

    # Reference-adjacent but design-split: one reference run, two design
    # fragments. Same contiguity rule as seg_ranges, applied to the other side.
    # D24-31 and D48-55 are NOT adjacent in the reference, so placing them
    # contiguously in the design must still give two fragments -- the reference
    # side of the same contiguity rule.
    joined_map = {f"D{n}": f"A{5 + i}" for i, n in enumerate(REF_ORDER)}
    joined = tmp / "fp_joined.pdb"
    _contiguous_design(joined, joined_map, {"A": 40})
    sc_j = tmp / "fp_joined.json"
    write_sidecar(sc_j, ref, joined_map)
    row = run_worker(tmp, "fp_joined", joined, sc_j)
    check("fixed positions: design-contiguous but reference-split stays split",
          int(row["n_segments"]), 2)
    check("fixed positions: two fragments for two segments",
          str(row["design_fixed_positions"]), "5:12,13:20")

    # --- NA contracts: a string that cannot be honest is absent, not partial ---

    # (1) a gapped design chain. ProteinMPNN positions are 1-indexed CHAIN
    # POSITIONS, so residue numbers from a chain that is not 1..N would fix the
    # wrong residues at exit 0.
    gapped = tmp / "fp_gapped.pdb"
    build_design(gapped)  # chain A: 1-8, 10-19, 30-37
    sc_g = tmp / "fp_gapped.json"
    write_sidecar(sc_g, ref, INDEX_MAP)
    row = run_worker(tmp, "fp_gapped", gapped, sc_g)
    check("gapped chain: row still OK", str(row["status"]), "OK")
    check_true("gapped chain: rmsd still reported", not pd.isna(row["rmsd"]))
    check_true("gapped chain: fixed positions NA",
               pd.isna(row.get("design_fixed_positions")))

    # (2) a chain numbered from 2: contiguous, but offset -- also not positions.
    offset_map = {k: f"A{int(v[1:]) + 1}" for k, v in idx_map.items()}
    offset = tmp / "fp_offset.pdb"
    _contiguous_design(offset, offset_map, {"A": 41})
    # Renumber chain A from 2 by dropping residue 1.
    st = gemmi.read_structure(str(offset))
    del st[0]["A"][0]
    st.setup_entities()
    st.write_pdb(str(offset))
    sc_o = tmp / "fp_offset.json"
    write_sidecar(sc_o, ref, offset_map)
    row = run_worker(tmp, "fp_offset", offset, sc_o)
    check("offset chain: row still OK", str(row["status"]), "OK")
    check_true("offset chain: fixed positions NA",
               pd.isna(row.get("design_fixed_positions")))

    # (3) a motif spanning TWO design chains: no single chain group can say it.
    two_map = {f"D{n}": f"A{1 + i}" for i, n in enumerate(SEG_A)}
    two_map.update({f"D{n}": f"B{1 + i}" for i, n in enumerate(SEG_B)})
    two = tmp / "fp_two_chains.pdb"
    _contiguous_design(two, two_map, {"A": 12, "B": 12})
    sc_t = tmp / "fp_two.json"
    write_sidecar(sc_t, ref, two_map)
    row = run_worker(tmp, "fp_two", two, sc_t)
    check("two design chains: row still OK", str(row["status"]), "OK")
    check("two design chains: design_chains records both",
          str(row["design_chains"]), "A,B")
    check_true("two design chains: fixed positions NA",
               pd.isna(row.get("design_fixed_positions")))

    # (4) an error row NAs it like every other column.
    row = run_worker(tmp, "fp_err", design, tmp / "nope.json")
    check_true("error row: status is error", str(row["status"]).startswith("error:"))
    check_true("error row: fixed positions NA",
               pd.isna(row.get("design_fixed_positions")))

    # (5) an --allow-partial row with n_matched < n_expected. The other columns
    # are measured on the matched subset exactly as before, but a fixed-positions
    # spec over that subset would leave the unmatched motif residue free to be
    # redesigned, so this column alone is NA. The residue is dropped from the
    # REFERENCE, not the design, so the design chain stays 1..N and gapless:
    # without the complete-match rule this row WOULD get a non-NA string
    # ("5:6,8:12,30:37"), which is what makes this check able to fail.
    dropped = SEG_A[2]  # D26 -> A7
    holed_ref = tmp / "fp_holed_ref.pdb"
    write_structure(
        holed_ref,
        {"D": [(n, REF_CODES[n], REF_COORDS[n]) for n in REF_ORDER if n != dropped]},
    )
    sc_p = tmp / "fp_partial.json"
    write_sidecar(sc_p, holed_ref, idx_map)

    # Without the flag the count mismatch is still an error, as before.
    row = run_worker(tmp, "fp_partial_noflag", design, sc_p)
    check_true("partial, no flag: still error",
               str(row["status"]).startswith("error:"))
    check_true("partial, no flag: fixed positions NA",
               pd.isna(row.get("design_fixed_positions")))

    import contextlib
    import io

    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        row = run_worker(tmp, "fp_partial", design, sc_p, allow_partial=True)
    check_true("partial+flag: status is warn", str(row["status"]).startswith("warn:"))
    check("partial+flag: n_expected", int(row["n_expected"]), 16)
    check("partial+flag: n_matched", int(row["n_matched"]), 15)
    check_true("partial+flag: rmsd still populated", not pd.isna(row["rmsd"]))
    check("partial+flag: rmsd on the matched subset", float(row["rmsd"]), 0.0,
          tol=1e-3)
    check("partial+flag: seg_ranges still populated (matched subset)",
          str(row["seg_ranges"]), "D24-25/D27-31/D48-55")
    check_true("partial+flag: seg_rmsds still populated",
               not pd.isna(row["seg_rmsds"]))
    check_true("partial+flag: fixed positions NA",
               pd.isna(row.get("design_fixed_positions")))
    check_true("partial+flag: NA reason written to stderr",
               "design_fixed_positions not emitted (NA)" in err.getvalue()
               and "15 of 16" in err.getvalue())


# ---------------------------------------------------------------------------
# --motif-chain (2026-10-09)
# ---------------------------------------------------------------------------

# Result TSV rows captured from the worker as it stood BEFORE --motif-chain was
# added, one per fixture shape below, with the temp dir written as '{}'. The
# regression test re-runs each fixture with the flag UNSET and compares the whole
# row -- every pre-change column, in order, including the status string -- against
# these bytes. They were produced by running both workers side by side, so this is
# a recorded diff of the old behaviour, not a restatement of the new one.
PRE_CHANGE_ROWS = {
    "gapped_chain": "name\tstatus\taligned_path\tn_expected\tn_matched\tseq_match_frac\tref_source\tref_matches_spec\tref_used\tspec_input\tref_chains\tdesign_chains\tref_motif_seq\tdesign_motif_seq\tsampled_contig\trmsd\tmax_dev\tn_segments\tworst_seg_rmsd\tseg_ranges\tseg_rmsds\tdesign_fixed_positions\tseg0_rmsd\tseg1_rmsd\r\nNAME\tOK\tALIGNED\t16\t16\t1.000\tsidecar\t1\t{}/ref.pdb\t{}/ref.pdb\tD\tA\tAGSVLWAGAGSVLWAG\tAGSVLWAGAGSVLWAG\t22P,D24,D31,28P\t0.000\t0.000\t2\t0.000\tD24-31/D48-55\t0.000/0.000\t\t0.000\t0.000\r\n",
    "missing_sidecar": "name\tstatus\taligned_path\tn_expected\tn_matched\tseq_match_frac\tref_source\tref_matches_spec\tref_used\tspec_input\tref_chains\tdesign_chains\tref_motif_seq\tdesign_motif_seq\tsampled_contig\trmsd\tmax_dev\tn_segments\tworst_seg_rmsd\tseg_ranges\tseg_rmsds\tdesign_fixed_positions\r\nNAME\terror: rfd3 sidecar JSON not found: {}/nope.json\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\r\n",
    "no_map": "name\tstatus\taligned_path\tn_expected\tn_matched\tseq_match_frac\tref_source\tref_matches_spec\tref_used\tspec_input\tref_chains\tdesign_chains\tref_motif_seq\tdesign_motif_seq\tsampled_contig\trmsd\tmax_dev\tn_segments\tworst_seg_rmsd\tseg_ranges\tseg_rmsds\tdesign_fixed_positions\r\nNAME\terror: rfd3 sidecar {}/sc_nomap.json has no 'diffused_index_map' (top-level keys: ['ckpt_path', 'metrics', 'seed', 'specification']). This is not motif-scaffolding output.\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\r\n",
    "partial_flag": "name\tstatus\taligned_path\tn_expected\tn_matched\tseq_match_frac\tref_source\tref_matches_spec\tref_used\tspec_input\tref_chains\tdesign_chains\tref_motif_seq\tdesign_motif_seq\tsampled_contig\trmsd\tmax_dev\tn_segments\tworst_seg_rmsd\tseg_ranges\tseg_rmsds\tdesign_fixed_positions\tseg0_rmsd\tseg1_rmsd\tseg2_rmsd\r\nNAME\twarn: matched 15 of 16 motif residues: 1 absent from reference holed_ref.pdb (D26)\tALIGNED\t16\t15\t0.267\tsidecar\t1\t{}/holed_ref.pdb\t{}/holed_ref.pdb\tD\tA\tAGVLWAGAGSVLWAG\tGGGGGGGGGGGGGGG\t22P,D24,D31,28P\t0.000\t0.000\t3\t0.000\tD24-25/D27-31/D48-55\t0.000/0.000/0.000\t\t0.000\t0.000\t0.000\r\n",
    "partial_noflag": "name\tstatus\taligned_path\tn_expected\tn_matched\tseq_match_frac\tref_source\tref_matches_spec\tref_used\tspec_input\tref_chains\tdesign_chains\tref_motif_seq\tdesign_motif_seq\tsampled_contig\trmsd\tmax_dev\tn_segments\tworst_seg_rmsd\tseg_ranges\tseg_rmsds\tdesign_fixed_positions\r\nNAME\terror: matched 15 of 16 motif residues: 1 absent from reference holed_ref.pdb (D26) (pass --allow-partial to measure the subset)\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\r\n",
    "perfect": "name\tstatus\taligned_path\tn_expected\tn_matched\tseq_match_frac\tref_source\tref_matches_spec\tref_used\tspec_input\tref_chains\tdesign_chains\tref_motif_seq\tdesign_motif_seq\tsampled_contig\trmsd\tmax_dev\tn_segments\tworst_seg_rmsd\tseg_ranges\tseg_rmsds\tdesign_fixed_positions\tseg0_rmsd\tseg1_rmsd\r\nNAME\tOK\tALIGNED\t16\t16\t0.250\tsidecar\t1\t{}/ref.pdb\t{}/ref.pdb\tD\tA\tAGSVLWAGAGSVLWAG\tGGGGGGGGGGGGGGGG\t22P,D24,D31,28P\t0.000\t0.000\t2\t0.000\tD24-31/D48-55\t0.000/0.000\t5:12,30:37\t0.000\t0.000\r\n",
    "perturbed": "name\tstatus\taligned_path\tn_expected\tn_matched\tseq_match_frac\tref_source\tref_matches_spec\tref_used\tspec_input\tref_chains\tdesign_chains\tref_motif_seq\tdesign_motif_seq\tsampled_contig\trmsd\tmax_dev\tn_segments\tworst_seg_rmsd\tseg_ranges\tseg_rmsds\tdesign_fixed_positions\tseg0_rmsd\tseg1_rmsd\r\nNAME\tOK\tALIGNED\t16\t16\t1.000\tsidecar\t1\t{}/ref.pdb\t{}/ref.pdb\tD\tA\tAGSVLWAGAGSVLWAG\tAGSVLWAGAGSVLWAG\t22P,D24,D31,28P\t0.519\t1.737\t2\t0.659\tD24-31/D48-55\t0.324/0.659\t\t0.324\t0.659\r\n",
    "ref_flag": "name\tstatus\taligned_path\tn_expected\tn_matched\tseq_match_frac\tref_source\tref_matches_spec\tref_used\tspec_input\tref_chains\tdesign_chains\tref_motif_seq\tdesign_motif_seq\tsampled_contig\trmsd\tmax_dev\tn_segments\tworst_seg_rmsd\tseg_ranges\tseg_rmsds\tdesign_fixed_positions\tseg0_rmsd\tseg1_rmsd\r\nNAME\tOK\tALIGNED\t16\t16\t0.250\tflag\t1\t{}/ref.pdb\t{}/ref.pdb\tD\tA\tAGSVLWAGAGSVLWAG\tGGGGGGGGGGGGGGGG\t22P,D24,D31,28P\t0.000\t0.000\t2\t0.000\tD24-31/D48-55\t0.000/0.000\t5:12,30:37\t0.000\t0.000\r\n",
    "rigid_moved": "name\tstatus\taligned_path\tn_expected\tn_matched\tseq_match_frac\tref_source\tref_matches_spec\tref_used\tspec_input\tref_chains\tdesign_chains\tref_motif_seq\tdesign_motif_seq\tsampled_contig\trmsd\tmax_dev\tn_segments\tworst_seg_rmsd\tseg_ranges\tseg_rmsds\tdesign_fixed_positions\tseg0_rmsd\tseg1_rmsd\r\nNAME\tOK\tALIGNED\t16\t16\t1.000\tsidecar\t1\t{}/ref.pdb\t{}/ref.pdb\tD\tA\tAGSVLWAGAGSVLWAG\tAGSVLWAGAGSVLWAG\t22P,D24,D31,28P\t0.001\t0.001\t2\t0.001\tD24-31/D48-55\t0.001/0.001\t\t0.001\t0.001\r\n",
    "two_chains": "name\tstatus\taligned_path\tn_expected\tn_matched\tseq_match_frac\tref_source\tref_matches_spec\tref_used\tspec_input\tref_chains\tdesign_chains\tref_motif_seq\tdesign_motif_seq\tsampled_contig\trmsd\tmax_dev\tn_segments\tworst_seg_rmsd\tseg_ranges\tseg_rmsds\tdesign_fixed_positions\tseg0_rmsd\tseg1_rmsd\r\nNAME\tOK\tALIGNED\t16\t16\t0.250\tsidecar\t1\t{}/ref.pdb\t{}/ref.pdb\tD\tA,B\tAGSVLWAGAGSVLWAG\tGGGGGGGGGGGGGGGG\t22P,D24,D31,28P\t0.000\t0.000\t2\t0.000\tD24-31/D48-55\t0.000/0.000\t\t0.000\t0.000\r\n"
}


def _regression_fixtures(tmp: Path) -> list[tuple[str, dict]]:
    """One fixture per shape the pre-change worker could land in: a clean graft, a
    gapped design chain, a two-design-chain motif, a rigid move, a per-residue
    perturbation, --allow-partial with and without the flag, a sidecar with no
    map, a missing sidecar, and an explicit --ref-structure."""
    ref = tmp / "ref.pdb"
    build_ref(ref)
    sc = tmp / "sc.json"
    write_sidecar(sc, ref, INDEX_MAP)

    idx = {f"D{n}": f"A{5 + i}" for i, n in enumerate(SEG_A)}
    idx.update({f"D{n}": f"A{30 + i}" for i, n in enumerate(SEG_B)})
    contig = tmp / "contig.pdb"
    _contiguous_design(contig, idx, {"A": 40})
    sc_c = tmp / "sc_c.json"
    write_sidecar(sc_c, ref, idx)

    two = {f"D{n}": f"A{1 + i}" for i, n in enumerate(SEG_A)}
    two.update({f"D{n}": f"B{1 + i}" for i, n in enumerate(SEG_B)})
    twop = tmp / "two.pdb"
    _contiguous_design(twop, two, {"A": 12, "B": 12})
    sc_t = tmp / "sc_t.json"
    write_sidecar(sc_t, ref, two)

    gapped = tmp / "gapped.pdb"
    build_design(gapped)

    moved = tmp / "moved.pdb"
    rot = rotation(np.array([0.2, 0.5, -0.8]), 0.9)
    build_design(moved, transform=lambda v: rot @ v + np.array([3.0, -4.0, 5.0]))
    pert = tmp / "pert.pdb"
    build_design(pert, perturb={SEG_A[0]: np.array([1.0, 0, 0]),
                                SEG_B[3]: np.array([0, 2.0, 0])})

    holed = tmp / "holed_ref.pdb"
    write_structure(holed, {"D": [(n, REF_CODES[n], REF_COORDS[n])
                                  for n in REF_ORDER if n != SEG_A[2]]})
    sc_p = tmp / "sc_p.json"
    write_sidecar(sc_p, holed, idx)
    sc_nomap = tmp / "sc_nomap.json"
    write_sidecar(sc_nomap, ref, None)

    return [
        ("perfect", dict(design=contig, sidecar=sc_c)),
        ("gapped_chain", dict(design=gapped, sidecar=sc)),
        ("two_chains", dict(design=twop, sidecar=sc_t)),
        ("rigid_moved", dict(design=moved, sidecar=sc)),
        ("perturbed", dict(design=pert, sidecar=sc)),
        ("partial_flag", dict(design=contig, sidecar=sc_p, allow_partial=True)),
        ("partial_noflag", dict(design=contig, sidecar=sc_p)),
        ("no_map", dict(design=contig, sidecar=sc_nomap)),
        ("missing_sidecar", dict(design=contig, sidecar=tmp / "nope.json")),
        ("ref_flag", dict(design=contig, sidecar=sc_c, ref=str(ref))),
    ]


def test_no_flag_regression(tmp: Path) -> None:
    """--motif-chain UNSET reproduces the pre-change worker byte for byte.

    The flag is additive by construction, but "by construction" is what every
    regression was before it happened. This compares the WHOLE result row -- all
    22 pre-change columns plus the per-segment ones, in order, with the status
    string and the float formatting -- against rows recorded from the worker as it
    stood before the change. The only permitted difference is the new
    fixed_positions_chain column, which is removed before the comparison (and
    checked separately below).
    """
    import contextlib
    import io

    work = tmp / "regress"
    work.mkdir()
    for label, kw in _regression_fixtures(work):
        name = f"rg_{label}"
        with contextlib.redirect_stderr(io.StringIO()):
            run_worker(work, name, **kw)  # flag unset
        # read_bytes, not read_text: text mode would translate the CRLF that
        # csv.writer emits into LF and the comparison would be against a
        # different file than the one on disk.
        raw = (work / f"{name}.tsv").read_bytes().decode()
        lines = raw.rstrip("\r\n").split("\r\n")
        header = lines[0].split("\t")
        rows = [line.split("\t") for line in lines[1:]]

        check_true(f"regression {label}: new column present",
                   "fixed_positions_chain" in header)
        drop = header.index("fixed_positions_chain")
        check(f"regression {label}: new column is NA with the flag unset "
              f"unless a spec was emitted",
              rows[0][drop] == "",
              rows[0][header.index("design_fixed_positions")] == "")

        kept = "\r\n".join(
            "\t".join(c for i, c in enumerate(cells) if i != drop)
            for cells in [header, *rows]
        ) + "\r\n"
        want = PRE_CHANGE_ROWS[label].replace("{}", str(work)).replace(
            "NAME", name).replace("ALIGNED", "")
        # aligned_path differs only because this harness does not write the PDB;
        # compare it as the pre-change capture did, with that one cell blanked.
        got = kept
        want_cells = [r.split("\t") for r in want.rstrip("\r\n").split("\r\n")]
        got_cells = [r.split("\t") for r in got.rstrip("\r\n").split("\r\n")]
        ap = want_cells[0].index("aligned_path")
        for row in (*want_cells[1:], *got_cells[1:]):
            row[ap] = ""
        check(f"regression {label}: row byte-identical to pre-change worker",
              got_cells, want_cells)


def _two_ref_chain_fixture(
    tmp: Path, stem: str, design_b_gap: bool = False
) -> tuple[Path, Path, dict[str, str]]:
    """The real binder-scaffolding shape, in miniature.

    Reference chain A (100-109) is a TARGET held fixed by the contig; reference
    chain B (30-33, 50-52, 95-99) carries the grafted CDR-like motif. rfd3 records
    BOTH in one diffused_index_map -- the whole point of the flag -- so the motif
    lands on two design chains: A1-10 (the target copy) and B11-14 / B45-47 /
    B87-91 (the graft).

    With *design_b_gap*, design chain B is missing residue 5, so it is no longer
    numbered 1..N in file order: restricting to B then STILL cannot produce a
    spec, which is the "restriction is not a free pass" case.
    """
    ref_a = list(range(100, 110))
    ref_b = list(range(30, 34)) + list(range(50, 53)) + list(range(95, 100))
    design_b = [11, 12, 13, 14, 45, 46, 47, 87, 88, 89, 90, 91]

    coords: dict[tuple[str, int], np.ndarray] = {}
    for i, n in enumerate(ref_a):
        coords[("A", n)] = helix(10, 0, np.array([0.0, 0.0, 0.0]))[i]
    for i, n in enumerate(ref_b):
        coords[("B", n)] = helix(12, 30, np.array([25.0, -12.0, 4.0]))[i]

    ref_path = tmp / f"{stem}_ref.pdb"
    write_structure(
        ref_path,
        {
            "A": [(n, "A", coords[("A", n)]) for n in ref_a],
            "B": [(n, "L", coords[("B", n)]) for n in ref_b],
        },
    )

    index_map = {f"A{n}": f"A{i + 1}" for i, n in enumerate(ref_a)}
    index_map.update({f"B{n}": f"B{design_b[i]}" for i, n in enumerate(ref_b)})

    placed = {
        (v[0], int(v[1:])): coords[(k[0], int(k[1:]))] for k, v in index_map.items()
    }
    chains: dict[str, list[tuple[int, str, np.ndarray]]] = {}
    for chain_name, length in (("A", 10), ("B", 100)):
        residues = []
        for num in range(1, length + 1):
            if design_b_gap and chain_name == "B" and num == 5:
                continue  # the hole that breaks 1..N
            residues.append(
                (
                    num,
                    "G",
                    placed.get(
                        (chain_name, num),
                        np.array([400.0 + num, 80.0 * (ord(chain_name) - 65), 0.0]),
                    ),
                )
            )
        chains[chain_name] = residues
    design_path = tmp / f"{stem}_design.pdb"
    write_structure(design_path, chains)

    sidecar = tmp / f"{stem}.json"
    write_sidecar(sidecar, ref_path, index_map)
    return design_path, sidecar, index_map


def test_motif_chain(tmp: Path) -> None:
    """--motif-chain restricts design_fixed_positions to ONE reference chain's
    residues -- and changes nothing else about the row.

    Models the campaign finding of 2026-10-09: in a binder-scaffolding contig that
    holds a target fixed, diffused_index_map names the whole target as well as the
    graft, so design_chains is multi-letter on every design and the unrestricted
    column is always NA.
    """
    work = tmp / "mc"
    work.mkdir()
    design, sidecar, index_map = _two_ref_chain_fixture(work, "mc")

    import contextlib
    import io

    def run(name: str, chain: str = "", **kw):
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            row = run_worker(work, name, design, sidecar, motif_chain=chain, **kw)
        return row, err.getvalue()

    # --- the situation the flag exists for: no flag, multi-chain, NA ---
    none_row, none_err = run("mc_noflag")
    check("no flag: status OK", str(none_row["status"]), "OK")
    check("no flag: design_chains lists both", str(none_row["design_chains"]), "A,B")
    check("no flag: ref_chains lists both", str(none_row["ref_chains"]), "A,B")
    check_true("no flag: fixed positions NA",
               pd.isna(none_row.get("design_fixed_positions")))
    check_true("no flag: NA reason is the multi-chain one",
               "maps into 2 design chains" in none_err)
    check_true("no flag: the reason points at the flag",
               "--motif-chain" in none_err)
    check_true("no flag: fixed_positions_chain NA too",
               pd.isna(none_row.get("fixed_positions_chain")))
    check("no flag: four segments over the whole map",
          int(none_row["n_segments"]), 4)
    check("no flag: seg_ranges over the whole map", str(none_row["seg_ranges"]),
          "A100-109/B30-33/B50-52/B95-99")

    # --- the fix: restrict to reference chain B, the grafted one ---
    row_b, err_b = run("mc_B", "B")
    check("motif-chain B: status still OK", str(row_b["status"]), "OK")
    check("motif-chain B: design-side spec for B only",
          str(row_b["design_fixed_positions"]), "11:14,45:47,87:91")
    check("motif-chain B: names the design chain",
          str(row_b["fixed_positions_chain"]), "B")
    check_true("motif-chain B: no NA reason on stderr",
               "design_fixed_positions not emitted" not in err_b)
    check_true("motif-chain B: spec holds no target residues",
               "1:10" not in str(row_b["design_fixed_positions"]))

    # The fragment count is the RESTRICTED segment count, which is NOT n_segments
    # (that still counts the whole map). Reading one for the other would be wrong.
    check("motif-chain B: three fragments, four segments",
          (len(str(row_b["design_fixed_positions"]).split(",")),
           int(row_b["n_segments"])), (3, 4))

    # Round-trip through prosapia's own parser: one chain group, exactly B's
    # mapped design residues.
    from prosapia.utils import parse_positions

    groups = parse_positions(
        str(row_b["design_fixed_positions"]), lambda *a, **k: None, "mc_B"
    )
    check("motif-chain B: parses to ONE chain group", len(groups), 1)
    check(
        "motif-chain B: covers exactly B's mapped design residues",
        sorted(groups[0]),
        sorted(int(v[1:]) for k, v in index_map.items() if k[0] == "B"),
    )

    # --- EVERY other column is untouched by the flag ---
    untouched = [c for c in none_row.index
                 if c not in ("name", "design_fixed_positions",
                              "fixed_positions_chain")]
    differing = [c for c in untouched
                 if str(none_row[c]).replace("mc_noflag", "") !=
                 str(row_b[c]).replace("mc_B", "")]
    check("motif-chain B: no other column changes", differing, [])
    check_true("motif-chain B: design_chains STILL lists both (honesty about "
               "what the RMSD covers)", str(row_b["design_chains"]) == "A,B")

    # --- the other reference chain works the same way ---
    row_a, _ = run("mc_A", "A")
    check("motif-chain A: the target chain's own spec",
          str(row_a["design_fixed_positions"]), "1:10")
    check("motif-chain A: names design chain A",
          str(row_a["fixed_positions_chain"]), "A")

    # --- a letter that is in no reference key: NA with the flag's own reason ---
    row_c, err_c = run("mc_C", "C")
    check("absent chain: status still OK", str(row_c["status"]), "OK")
    check_true("absent chain: fixed positions NA",
               pd.isna(row_c.get("design_fixed_positions")))
    check_true("absent chain: fixed_positions_chain NA",
               pd.isna(row_c.get("fixed_positions_chain")))
    check_true("absent chain: reason names the flag and the chain",
               "--motif-chain C given but no reference residue in chain C" in err_c)
    check_true("absent chain: reason lists the chains that ARE in the map",
               "reference chain(s) A, B" in err_c)
    check_true("absent chain: rmsd still measured",
               not pd.isna(row_c["rmsd"]))
    check("absent chain: rmsd unchanged by the bad flag",
          float(row_c["rmsd"]), float(none_row["rmsd"]), tol=1e-9)

    # --- restriction is not a free pass: the 1..N check still applies, to
    #     whichever chain the restricted residues landed on ---
    gap_design, gap_sidecar, _ = _two_ref_chain_fixture(work, "mcgap",
                                                        design_b_gap=True)
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        row_gap = run_worker(work, "mc_gap_B", gap_design, gap_sidecar,
                             motif_chain="B")
    check("gapped B: status still OK", str(row_gap["status"]), "OK")
    check_true("gapped B: rmsd still measured", not pd.isna(row_gap["rmsd"]))
    check_true("gapped B: fixed positions NA",
               pd.isna(row_gap.get("design_fixed_positions")))
    check_true("gapped B: fixed_positions_chain NA",
               pd.isna(row_gap.get("fixed_positions_chain")))
    check_true("gapped B: reason is the 1..N one, for chain B",
               "design chain 'B'" in err.getvalue()
               and "1-indexed chain positions" in err.getvalue())

    # ... while chain A of the SAME design is still fine, because the check is
    # about the chain the residues landed on, not about the whole file.
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        row_gap_a = run_worker(work, "mc_gap_A", gap_design, gap_sidecar,
                               motif_chain="A")
    check("gapped B, chain A: spec still emitted",
          str(row_gap_a["design_fixed_positions"]), "1:10")
    check("gapped B, chain A: names design chain A",
          str(row_gap_a["fixed_positions_chain"]), "A")

    # --- an error row NAs the new column like every other ---
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        row_err = run_worker(work, "mc_err", design, work / "nope.json",
                             motif_chain="B")
    check_true("error row: status is error", str(row_err["status"]).startswith("error:"))
    check_true("error row: fixed_positions_chain NA",
               pd.isna(row_err.get("fixed_positions_chain")))


def test_motif_chain_manifest(tmp: Path) -> None:
    """The flag reaches the task script in the field motifrmsd.sh cuts, and a
    malformed letter dies at submit rather than per design."""
    results = tmp / "mfrun" / "table0" / "rfd3" / "results_shard_0"
    results.mkdir(parents=True)
    design = results / "d_model_0.cif.gz"
    design.write_bytes(b"")
    (results / "d_model_0.json").write_text("{}")
    ready = pd.DataFrame({"rfdiffusion3_path": [str(design)]}, index=["d0"])

    def build(motif_chain: str):
        args = Namespace(
            gpus_per_task=1, ref_structure="", allow_partial=False,
            motif_chain=motif_chain, input_column="rfdiffusion3_path",
        )
        ctx = Namespace(args=args, ready=ready)
        return run_motifrmsd.build_motifrmsd_manifest(ctx)  # type: ignore[arg-type]

    # Field order is the contract with motifrmsd.sh: f5 is the chain, f6 the
    # allow-partial flag. An empty f5 must NOT shift f6.
    row = build("")[0]
    check("manifest: six fields", len(row), 6)
    check("manifest: empty motif chain in f5", row[4], "")
    check("manifest: allow_partial stays in f6", row[5], "no")
    row = build("D")[0]
    check("manifest: motif chain in f5", row[4], "D")
    check("manifest: allow_partial still f6", row[5], "no")
    check("manifest: whitespace stripped", build("  D  ")[0][4], "D")

    for bad in ("D1", "A,B", "3", "D-E", "D24-31"):
        try:
            build(bad)
        except ValueError as e:
            check_true(f"manifest: {bad!r} refused at submit",
                       "--motif-chain" in str(e))
        else:
            FAILURES.append(f"manifest: {bad!r} should not have built a manifest")
            global CHECKS
            CHECKS += 1


def test_reference_policy(tmp: Path) -> None:
    """--ref-structure wins; agreement is judged by content, not by path string."""
    ref = tmp / "ref.pdb"
    build_ref(ref)
    design = tmp / "ok.pdb"
    build_design(design)
    sc = tmp / "policy.json"
    write_sidecar(sc, ref, INDEX_MAP)

    # Same bytes under a different path -> matches_spec 1, status OK.
    copy = tmp / "copy_of_ref.pdb"
    copy.write_bytes(ref.read_bytes())
    row = run_worker(tmp, "samebytes", design, sc, ref=str(copy))
    check("same bytes: status", str(row["status"]), "OK")
    check("same bytes: ref_source", str(row["ref_source"]), "flag")
    check("same bytes: ref_matches_spec", int(row["ref_matches_spec"]), 1)
    check("same bytes: ref_used", str(row["ref_used"]), str(copy))

    # Genuinely different file -> matches_spec 0 and a warn status, numbers kept.
    moved = tmp / "moved_ref.pdb"
    write_structure(
        moved,
        {
            "D": [
                (n, REF_CODES[n], REF_COORDS[n] + np.array([0.0, 0.0, 0.4]))
                for n in REF_ORDER
            ]
        },
    )
    row = run_worker(tmp, "diffbytes", design, sc, ref=str(moved))
    check_true("different ref: warn", str(row["status"]).startswith("warn:"))
    check("different ref: ref_matches_spec", int(row["ref_matches_spec"]), 0)
    check_true("different ref: rmsd still written", not pd.isna(row.get("rmsd")))

    # Spec path unreadable -> matches_spec NA, never an implied 1.
    sc2 = tmp / "gone.json"
    write_sidecar(sc2, tmp / "vanished.pdb", INDEX_MAP)
    row = run_worker(tmp, "specgone", design, sc2, ref=str(ref))
    check("spec gone: status", str(row["status"]), "OK")
    check_true("spec gone: ref_matches_spec NA", pd.isna(row.get("ref_matches_spec")))
    check("spec gone: spec_input recorded", str(row["spec_input"]),
          str(tmp / "vanished.pdb"))


def test_sequence_trust(tmp: Path) -> None:
    """seq_match_frac reports, and never gates: a poly-glycine design is OK."""
    ref = tmp / "ref.pdb"
    build_ref(ref)
    sc = tmp / "seq.json"
    write_sidecar(sc, ref, INDEX_MAP)
    polyg = tmp / "polyg.pdb"
    build_design(polyg, codes={n: "G" for n in REF_ORDER})
    row = run_worker(tmp, "polyg", polyg, sc)
    check("polyG: status still OK", str(row["status"]), "OK")
    check("polyG: design_motif_seq", str(row["design_motif_seq"]), "G" * 16)
    expected = sum(REF_CODES[n] == "G" for n in REF_ORDER) / 16
    check("polyG: seq_match_frac", float(row["seq_match_frac"]), expected, tol=1e-9)
    check("polyG: rmsd unaffected", float(row["rmsd"]), 0.0, tol=1e-3)
    check(
        "polyG: ref_motif_seq in segment order",
        str(row["ref_motif_seq"]),
        "".join(REF_CODES[n] for n in REF_ORDER),
    )


def test_sidecar_resolution(tmp: Path) -> None:
    """The .src route beats a stem glob, which is why it exists."""
    cache = tmp / "run" / ".cif_to_pdb"
    cache.mkdir(parents=True)
    # Two runs, identical file stems, different output dirs -- the real campaign's
    # situation (table0_seg2 and table0_seg3 both hold ..._model_0.cif.gz).
    for label, digest in (("seg2", "8ecf9c29"), ("seg3", "ac344f9c")):
        results = tmp / "run" / f"table0_{label}" / f"rfd3_{label}" / "results_shard_0"
        results.mkdir(parents=True)
        source = results / "shard_0_x_diff_0_model_0.cif.gz"
        source.write_bytes(b"")
        (results / "shard_0_x_diff_0_model_0.json").write_text("{}")
        cached = cache / f"shard_0_x_diff_0_model_0.cif.{digest}.pdb"
        cached.write_text("")
        cached.with_suffix(".src").write_text(str(source))

        got, reason = run_motifrmsd.resolve_sidecar(cached)
        check(f"resolve {label}: reason", reason, "")
        check(
            f"resolve {label}: follows .src into the right table",
            str(got),
            str(results / "shard_0_x_diff_0_model_0.json"),
        )

    # A design path straight at an rfd3 .cif.gz resolves by sibling.
    direct = tmp / "run" / "table0_seg2" / "rfd3_seg2" / "results_shard_0"
    got, reason = run_motifrmsd.resolve_sidecar(direct / "shard_0_x_diff_0_model_0.cif.gz")
    check("resolve direct .cif.gz", str(got), str(direct / "shard_0_x_diff_0_model_0.json"))

    # Nothing to find -> (None, reason), never a guess.
    orphan = cache / "orphan.cif.deadbeef.pdb"
    orphan.write_text("")
    got, reason = run_motifrmsd.resolve_sidecar(orphan)
    check("resolve orphan: None", got, None)
    check_true("resolve orphan: reason given", bool(reason))


def test_bad_src_record_costs_one_row(tmp: Path) -> None:
    """Item 2 (2026-10-08): an empty / undecodable / path-less .src cache record
    is reported per design, never raised -- before this, an empty .src made
    _json_beside(Path('.')) raise ValueError and abort the WHOLE submit."""
    cache = tmp / "srcrun" / ".cif_to_pdb"
    cache.mkdir(parents=True)
    results = tmp / "srcrun" / "table0" / "rfd3" / "results_shard_0"
    results.mkdir(parents=True)
    good_src = results / "good_model_0.cif.gz"
    good_src.write_bytes(b"")
    (results / "good_model_0.json").write_text("{}")

    def cached(stem: str, record: bytes) -> Path:
        pdb = cache / f"{stem}.cif.0000.pdb"
        pdb.write_text("")
        pdb.with_suffix(".src").write_bytes(record)
        return pdb

    good = cached("good", str(good_src).encode())
    bad = {
        "empty": cached("empty", b""),
        "whitespace": cached("whitespace", b"  \n"),
        "not utf-8": cached("notutf8", b"\xff\xfe\x00garbage"),
        "root dir": cached("rootdir", b"/"),
    }
    for label, pdb in bad.items():
        try:
            got, reason = run_motifrmsd.resolve_sidecar(pdb)
        except Exception as e:  # noqa: BLE001
            FAILURES.append(f"bad .src ({label}): raised {type(e).__name__}: {e}")
            continue
        check(f"bad .src ({label}): no sidecar", got, None)
        check_true(f"bad .src ({label}): reason names the record", ".src" in reason)

    # The builder: one good design and every bad one -> ALL submitted, the bad
    # ones with an empty sidecar field (which the worker turns into error:).
    names = ["d_good", *[f"d_{k.replace(' ', '_')}" for k in bad]]
    paths = [str(good), *[str(v) for v in bad.values()]]
    ready = pd.DataFrame({"rfdiffusion3_path": paths}, index=names)
    args = Namespace(
        gpus_per_task=1, ref_structure="", allow_partial=False, motif_chain="",
        input_column="rfdiffusion3_path",
    )
    ctx = Namespace(args=args, ready=ready)
    try:
        rows = run_motifrmsd.build_motifrmsd_manifest(ctx)  # type: ignore[arg-type]
    except Exception as e:  # noqa: BLE001
        FAILURES.append(f"builder with bad .src: raised {type(e).__name__}: {e}")
        return
    check("builder: every design submitted", len(rows), len(names))
    check("builder: CPU-only forced", args.gpus_per_task, 0)
    by_name = {r[0]: r for r in rows}
    check(
        "builder: good design keeps its sidecar",
        by_name["d_good"][2],
        str(results / "good_model_0.json"),
    )
    for n in names[1:]:
        check(f"builder: {n} has empty sidecar field", by_name[n][2], "")

    # ... and that empty field becomes a per-design error row, NA numbers.
    ref = tmp / "ref.pdb"
    build_ref(ref)
    design = tmp / "ok.pdb"
    build_design(design)
    row = run_worker(tmp, "emptysidecar", design, Path(""))
    check_true("empty sidecar field: error", str(row["status"]).startswith("error:"))
    check_true(
        "empty sidecar field: says it was not located",
        "could not be located" in str(row["status"]),
    )
    check_true("empty sidecar field: rmsd NA", pd.isna(row.get("rmsd")))


def test_collector(tmp: Path) -> None:
    """The collector keeps OK/warn numbers and NAs every error column."""
    out_dir = tmp / "collect_out"
    out_dir.mkdir()
    ref = tmp / "ref.pdb"
    build_ref(ref)
    sc = tmp / "coll.json"
    write_sidecar(sc, ref, INDEX_MAP)
    design = tmp / "ok.pdb"
    build_design(design)

    ok_row = run_worker(tmp, "d_ok", design, sc)
    (out_dir / "d_ok.tsv").write_text((tmp / "d_ok.tsv").read_text())
    err_row = run_worker(tmp, "d_err", design, tmp / "absent.json")
    (out_dir / "d_err.tsv").write_text((tmp / "d_err.tsv").read_text())
    assert str(ok_row["status"]) == "OK" and str(err_row["status"]).startswith("error:")

    ctx = Namespace(out_dir=out_dir)
    one = collect_motifrmsd.collect_motifrmsd(ctx)  # type: ignore[arg-type]

    got = list(one(Namespace(name="d_ok", out_dir=out_dir, lookup=None)))  # type: ignore[arg-type]
    check("collect ok: one row", len(got), 1)
    check("collect ok: status", got[0].status, "OK")
    check("collect ok: rmsd", float(got[0].data["rmsd"]), 0.0, tol=1e-3)
    check("collect ok: seg columns present", "seg1_rmsd" in got[0].data, True)
    check_true("collect ok: path set", bool(got[0].path))

    got = list(one(Namespace(name="d_err", out_dir=out_dir, lookup=None)))  # type: ignore[arg-type]
    check("collect err: one row", len(got), 1)
    check_true("collect err: status", str(got[0].status).startswith("error:"))
    check_true("collect err: rmsd NA", got[0].data["rmsd"] is None)
    check_true("collect err: n_matched NA", got[0].data["n_matched"] is None)
    check_true(
        "collect err: no seg columns",
        not any(k.startswith("seg") and k.endswith("_rmsd") for k in got[0].data),
    )
    check("collect err: path empty", got[0].path, "")

    got = list(one(Namespace(name="absent_design", out_dir=out_dir, lookup=None)))  # type: ignore[arg-type]
    check("collect missing: status", got[0].status, "missing")
    check_true("collect missing: rmsd NA", got[0].data["rmsd"] is None)

    # Every column the collector declares must be one the worker can write.
    worker_columns = set(worker.METRIC_COLUMNS)
    declared = set(collect_motifrmsd.METRIC_COLUMNS)
    check("collector/worker column sets agree", declared, worker_columns)


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        test_kabsch_identity_and_reflection()
        test_parse_and_segments()
        test_perfect_graft(tmp)
        test_known_displacement(tmp)
        test_global_fit_is_shared_across_segments(tmp)
        test_error_contracts(tmp)
        test_map_not_one_to_one(tmp)
        test_reference_adjacent_segments_split_by_design(tmp)
        test_design_fixed_positions(tmp)
        test_no_flag_regression(tmp)
        test_motif_chain(tmp)
        test_motif_chain_manifest(tmp)
        test_reference_policy(tmp)
        test_sequence_trust(tmp)
        test_sidecar_resolution(tmp)
        test_bad_src_record_costs_one_row(tmp)
        test_collector(tmp)

    if FAILURES:
        print(f"FAILED {len(FAILURES)} of {CHECKS} checks:")
        for f in FAILURES:
            print(f"  - {f}")
        return 1
    print(f"OK: {CHECKS} checks passed.")
    print(
        "NOTE: this proves the geometry and the error contracts only. It cannot "
        "prove the chain conventions, the residue-naming convention or the "
        "sidecar layout of REAL rfd3 output -- see the module docstring."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
