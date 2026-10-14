#!/usr/bin/env python3
"""
Self-test for posrmsd: synthetic structures with a hand-checkable answer, the
manifest builder against a real prosapia DataManager lineage walk, the collector,
and the task script executed through the REAL prosapia task prelude.

WHAT THIS PROVES, AND WHAT IT CANNOT.

It proves the arithmetic and the contracts: that a rigidly moved A gives RMSD 0;
that the fit agrees with an independent quaternion (Horn) superposition; that a
mirror image is not fitted; that per-segment RMSDs come from ONE global fit; that
positions index the chain's polymer in FILE ORDER (not by residue number, not
counting waters/ligands); that the segment split is identical to prosapia's own
parse_positions; that every refusal path writes a non-OK status with NA; that the
lineage-resolved --col-b / --positions-column behave like usalign's --col-b; and
that the manifest fields reach the worker in the order posrmsd.sh cuts them.

It CANNOT prove the premise or the conventions of real files. The structures here
are written by this file and agree with it by construction. UNKNOWN until real
data:

  * whether a real Boltz prediction and its rfd3/chainsel parent really have the
    same polymer length and ordering on the chosen chains (a_len == b_len,
    seq_match_frac ~1.0 is the evidence);
  * which chain letter the binder carries in each real file (gemmi's chain name =
    auth_asym_id for mmCIF), i.e. what --chain-a/--chain-b must be;
  * whether motifrmsd_design_fixed_positions, as actually collected, is in the
    file-order numbering of the structure posrmsd reads as B (motifrmsd only emits
    it for a chain numbered 1..N, which makes the two coincide -- on the rfd3
    output; a chainsel-extracted chain with --renumber-from would break that);
  * whether rfd3 writes real residue identities at motif positions (if it wrote
    UNK there, seq_match_frac would be 0 by the X rule, not 1).

Run:  uv run python tools/posrmsd/test_posrmsd.py
"""

import importlib.util
import json
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
from prosapia.utils import parse_positions

HERE = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


worker = _load("posrmsd_worker")
run_posrmsd = _load("run_posrmsd")
collect_posrmsd = _load("collect_posrmsd")

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


def check_true(label: str, value: bool) -> None:
    global CHECKS
    CHECKS += 1
    if not value:
        FAILURES.append(f"{label}: expected True")


def check_raises(label: str, fn, needle: str, exc=ValueError) -> None:
    global CHECKS
    CHECKS += 1
    try:
        fn()
    except exc as e:
        if needle not in str(e):
            FAILURES.append(f"{label}: raised but message lacks {needle!r}: {e}")
        return
    except Exception as e:  # noqa: BLE001
        FAILURES.append(f"{label}: raised {type(e).__name__} not {exc.__name__}: {e}")
        return
    FAILURES.append(f"{label}: did not raise")


# ---------------------------------------------------------------------------
# Synthetic structures
# ---------------------------------------------------------------------------

AA3 = {"A": "ALA", "G": "GLY", "S": "SER", "V": "VAL", "L": "LEU", "W": "TRP",
       "K": "LYS", "E": "GLU", "X": "UNK"}


def write_structure(path: Path, chains: dict, het: dict | None = None) -> None:
    """CA-only structure. ``chains`` = {chain: [(resnum, code, xyz[, has_ca])]};
    ``het`` = {chain: [(resnum, resname, xyz)]} appended AFTER the polymer (waters,
    ligands). Written as PDB or mmCIF by suffix."""
    st = gemmi.Structure()
    model = gemmi.Model("1")
    het = het or {}
    for chain_name, residues in chains.items():
        chain = gemmi.Chain(chain_name)
        for item in residues:
            resnum, code, xyz = item[:3]
            has_ca = item[3] if len(item) > 3 else True
            res = gemmi.Residue()
            res.name = AA3[code]
            res.seqid = gemmi.SeqId(resnum, " ")
            res.het_flag = "A"
            atom = gemmi.Atom()
            atom.name = "CA" if has_ca else "CB"
            atom.element = gemmi.Element("C")
            atom.pos = gemmi.Position(*map(float, xyz))
            res.add_atom(atom)
            chain.add_residue(res)
        for resnum, resname, xyz in het.get(chain_name, []):
            res = gemmi.Residue()
            res.name = resname
            res.seqid = gemmi.SeqId(resnum, " ")
            res.het_flag = "H"
            atom = gemmi.Atom()
            atom.name = "O" if resname == "HOH" else "C1"
            atom.element = gemmi.Element("O" if resname == "HOH" else "C")
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


def helix(n: int) -> np.ndarray:
    """n points on an alpha-helix-like curve: non-collinear, well conditioned."""
    t = np.arange(n)
    return np.stack([2.3 * np.cos(t * 1.745), 2.3 * np.sin(t * 1.745), 1.5 * t], axis=1)


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


SEQ = "AGSVLWKEAGSVLWKEAGSVLWKEAGSVLW"  # 30 residues
N = len(SEQ)
COORDS = helix(N)
R0 = rotation([1, 2, 3], 0.7)
T0 = np.array([5.0, -3.0, 12.0])


def residues(coords=COORDS, seq=SEQ, start=1, gaps=()):
    """[(resnum, code, xyz)] numbered from ``start`` in file order, skipping any
    residue NUMBER in ``gaps`` (the file still lists every residue)."""
    out, num = [], start
    for code, xyz in zip(seq, coords):
        while num in gaps:
            num += 1
        out.append((num, code, xyz))
        num += 1
    return out


def segments_json(labels, segs) -> str:
    return json.dumps({"labels": labels, "segments": segs})


SEGS = [[3, 4, 5, 6], [12, 13, 14, 15, 16], [25, 26, 27]]
LABELS = ["3:6", "12:16", "25:27"]


def run_worker(tmp: Path, name: str, a: Path, b: Path, chain_a="", chain_b="",
               labels=LABELS, segs=SEGS, payload: str | None = None) -> pd.Series:
    args = Namespace(
        name=name, a=a, b=b, chain_a=chain_a, chain_b=chain_b,
        segments=payload if payload is not None else segments_json(labels, segs),
        out_struct=tmp / f"{name}_sup.cif", result_tsv=tmp / f"{name}.tsv",
    )
    worker.run(args)
    return pd.read_csv(args.result_tsv, sep="\t", dtype=str,
                       keep_default_na=False, na_values=[""]).iloc[0]


def horn_rmsd(p: np.ndarray, q: np.ndarray) -> float:
    """Independent superposition RMSD via Horn's quaternion method."""
    pc, qc = p - p.mean(0), q - q.mean(0)
    S = pc.T @ qc
    Sxx, Sxy, Sxz = S[0]
    Syx, Syy, Syz = S[1]
    Szx, Szy, Szz = S[2]
    K = np.array([
        [Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx],
        [Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz],
        [Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy],
        [Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz],
    ])
    lam = np.linalg.eigvalsh(K)[-1]
    e = (pc ** 2).sum() + (qc ** 2).sum() - 2 * lam
    return float(np.sqrt(max(e, 0.0) / len(p)))


# ---------------------------------------------------------------------------
# Position parsing: must be prosapia's, byte for byte
# ---------------------------------------------------------------------------


def test_parse_segments() -> None:
    lookup_vals = {"s": 11, "e": 18}

    def lookup(name, col):
        return lookup_vals.get(col, pd.NA)

    cases = {
        "11:18,45:54,87:92": ["11:18", "45:54", "87:92"],
        "5,7,9": ["5", "7", "9"],
        "[1:3,10:12]": ["1:3", "10:12"],
        " 11 : 18 , 20:22 ": ["11:18", "20:22"],
        "1:3,": ["1:3"],
        "{s}:{e},30:32": ["11:18", "30:32"],
        "{s}:{e+2}": ["11:20"],
    }
    for spec, want_labels in cases.items():
        labels, segs = run_posrmsd.parse_segments(spec, lookup, "d")
        check(f"parse {spec!r}: labels", labels, want_labels)
        flat = [p for s in segs for p in s]
        check(f"parse {spec!r}: == prosapia.parse_positions",
              [flat], parse_positions(spec, lookup, "d"))
        check(f"parse {spec!r}: one segment per fragment", len(segs), len(want_labels))
    _, segs = run_posrmsd.parse_segments("11:18,45:54", lookup, "d")
    check("parse: inclusive ends", (segs[0][0], segs[0][-1], len(segs[0])), (11, 18, 8))

    refusals = {
        "1:3/5:7": "'/' chain break",
        "10:": "open-ended",
        ":10": "open-ended",
        "18:11,1:5": "expands to no positions",
        "0:3": "< 1",
        "1:5,4:8": "more than once",
        "1:2": "needs at least 3",
        "a:b": "non-integer",
        "   ": "empty position spec",
        "1:2:3": "malformed",
        "{missing}:20": "not found",
    }
    for spec, needle in refusals.items():
        check_raises(f"parse refuses {spec!r}",
                     lambda spec=spec: run_posrmsd.parse_segments(spec, lookup, "d"),
                     needle)


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def test_kabsch() -> None:
    p = COORDS[:10]
    rot, trans = worker.kabsch(p @ R0.T + T0, p)
    moved = (p @ R0.T + T0) @ rot.T + trans
    check("kabsch: rigid motion recovered", worker.rmsd(moved, p), 0.0, tol=1e-9)
    check("kabsch: det +1", float(np.linalg.det(rot)), 1.0, tol=1e-9)
    mirror = p * np.array([1, 1, -1])
    rot, trans = worker.kabsch(mirror, p)
    check("kabsch: reflection never used", float(np.linalg.det(rot)), 1.0, tol=1e-9)
    check_true("kabsch: mirror not fitted to 0",
               worker.rmsd(mirror @ rot.T + trans, p) > 0.5)
    rng = np.random.default_rng(0)
    q = p + rng.normal(scale=0.8, size=p.shape)
    rot, trans = worker.kabsch(q, p)
    check("kabsch: agrees with independent Horn quaternion fit",
          worker.rmsd(q @ rot.T + trans, p), horn_rmsd(q, p), tol=1e-6)
    check_raises("kabsch: <3 points", lambda: worker.kabsch(p[:2], p[:2]), ">=3")


# ---------------------------------------------------------------------------
# Worker end to end
# ---------------------------------------------------------------------------


def test_perfect_and_frame(tmp: Path) -> None:
    b = tmp / "b.pdb"
    a = tmp / "a.cif"
    write_structure(b, {"A": residues()})
    write_structure(a, {"A": residues(COORDS @ R0.T + T0)})
    row = run_worker(tmp, "perfect", a, b)
    check("perfect: status", row["status"], "OK")
    check("perfect: rmsd", float(row["rmsd"]), 0.0, tol=1e-3)
    check("perfect: max_dev", float(row["max_dev"]), 0.0, tol=1e-3)
    check("perfect: worst_seg", float(row["worst_seg_rmsd"]), 0.0, tol=1e-3)
    check("perfect: n_positions", int(row["n_positions"]), 12)
    check("perfect: seq_match_frac", float(row["seq_match_frac"]), 1.0, tol=1e-9)
    check("perfect: a_motif_seq comma-joined segments", row["a_motif_seq"],
          ",".join("".join(SEQ[p - 1] for p in s) for s in SEGS))
    check("perfect: seqs equal", row["a_motif_seq"], row["b_motif_seq"])
    check("perfect: seg_ranges", row["seg_ranges"], "3:6/12:16/25:27")
    check("perfect: n_segments", int(row["n_segments"]), 3)
    check("perfect: seg2 present", float(row["seg2_rmsd"]), 0.0, tol=1e-3)
    check("perfect: offsets 0", (int(row["a_resnum_offset"]), int(row["b_resnum_offset"])), (0, 0))
    check("perfect: gapless", (row["a_gapless"], row["b_gapless"]), ("1", "1"))
    check("perfect: chains", (row["a_chain"], row["b_chain"]), ("A", "A"))
    check("perfect: lens", (int(row["a_len"]), int(row["b_len"])), (N, N))
    # A was MOVED into B's frame: the written superposition sits on B.
    sup = gemmi.read_structure(str(row["sup_path"]))
    xyz = np.array([[r[0].pos.x, r[0].pos.y, r[0].pos.z] for r in sup[0]["A"]])
    check("perfect: sup file is A moved onto B (whole chain)",
          worker.rmsd(xyz, COORDS), 0.0, tol=1e-2)


def test_global_fit_and_independent_check(tmp: Path) -> None:
    """Segment 1 of A is shifted 3 A: a per-segment refit would give 0 for it;
    the shared global fit must not. RMSD must equal Horn's, and the per-segment
    RMSDs must recombine into the global one."""
    shifted = COORDS.copy()
    idx = [p - 1 for p in SEGS[1]]
    shifted[idx] += np.array([3.0, 0.0, 0.0])
    b = tmp / "gb.pdb"
    a = tmp / "ga.pdb"
    write_structure(b, {"A": residues()})
    write_structure(a, {"A": residues(shifted @ R0.T + T0)})
    row = run_worker(tmp, "global", a, b)
    sel = [p - 1 for s in SEGS for p in s]
    want = horn_rmsd((shifted @ R0.T + T0)[sel], COORDS[sel])
    check("global: rmsd == Horn", float(row["rmsd"]), want, tol=2e-3)
    segs = [float(row[f"seg{i}_rmsd"]) for i in range(3)]
    check_true("global: shifted segment not refit to 0", segs[1] > 1.0)
    check_true("global: other segments pay for it too", segs[0] > 0.05 and segs[2] > 0.05)
    recombined = np.sqrt(sum(r * r * len(s) for r, s in zip(segs, SEGS)) / 12)
    check("global: segment RMSDs recombine to global", recombined, float(row["rmsd"]), tol=5e-3)
    check("global: worst_seg is the max", float(row["worst_seg_rmsd"]), max(segs), tol=1e-3)
    check_true("global: max_dev >= rmsd", float(row["max_dev"]) >= float(row["rmsd"]))
    check("global: seg_rmsds string", row["seg_rmsds"], "/".join(f"{v:.3f}" for v in segs))

    mirror = tmp / "mirror.pdb"
    write_structure(mirror, {"A": residues(COORDS * np.array([1, 1, -1]))})
    row = run_worker(tmp, "mirror", mirror, b)
    check("mirror: status OK", row["status"], "OK")
    check_true("mirror: not fitted to 0", float(row["rmsd"]) > 0.5)


def test_file_order_not_resnum(tmp: Path) -> None:
    """B numbered 101.. and A numbered 1..: position 3 is the 3rd residue of each."""
    b = tmp / "fb.pdb"
    a = tmp / "fa.pdb"
    write_structure(b, {"A": residues(start=101)})
    write_structure(a, {"A": residues(COORDS @ R0.T)})
    row = run_worker(tmp, "fileorder", a, b)
    check("file order: status", row["status"], "OK")
    check("file order: rmsd 0", float(row["rmsd"]), 0.0, tol=1e-3)
    check("file order: seq match 1", float(row["seq_match_frac"]), 1.0, tol=1e-9)
    check("file order: b offset 100", int(row["b_resnum_offset"]), 100)
    check("file order: a offset 0", int(row["a_resnum_offset"]), 0)
    check("file order: b still gapless", row["b_gapless"], "1")

    # A numbering gap at residue NUMBERS 10-12 (the file still lists 30 residues).
    g = tmp / "gap.pdb"
    write_structure(g, {"A": residues(gaps=(10, 11, 12))})
    row = run_worker(tmp, "gap", a, g)
    check("gap: status OK (file order is still defined)", row["status"], "OK")
    check("gap: rmsd 0 (same residues by file order)", float(row["rmsd"]), 0.0, tol=1e-3)
    check("gap: gapless 0", row["b_gapless"], "0")
    check_true("gap: offset NA when not constant over the picks",
               pd.isna(row["b_resnum_offset"]))
    row = run_worker(tmp, "gap_after", a, g, labels=["15:20"], segs=[list(range(15, 21))])
    check("gap: constant offset after the gap", int(row["b_resnum_offset"]), 3)


def test_hetero_residues_do_not_count(tmp: Path) -> None:
    b = tmp / "hb.pdb"
    a = tmp / "ha.cif"
    het = {"A": [(201, "HOH", (40, 40, 40)), (202, "HOH", (41, 40, 40)),
                 (203, "LIG", (42, 40, 40))]}
    write_structure(b, {"A": residues()}, het=het)
    write_structure(a, {"A": residues(COORDS @ R0.T)}, het=het)
    row = run_worker(tmp, "het", a, b)
    check("het: status", row["status"], "OK")
    check("het: a_len counts polymer only", int(row["a_len"]), N)
    check("het: b_len counts polymer only", int(row["b_len"]), N)
    row = run_worker(tmp, "het_beyond", a, b, labels=["29:32"], segs=[[29, 30, 31, 32]])
    check_true("het: a position past the polymer is an error even with waters",
               row["status"].startswith("error:") and "has 30 polymer residues" in row["status"])


def test_chain_selection(tmp: Path) -> None:
    decoy = residues(COORDS[::-1] + 20.0)  # a different shape, listed FIRST
    b2 = tmp / "two_b.pdb"
    write_structure(b2, {"A": decoy, "D": residues()})
    a1 = tmp / "one_a.pdb"
    write_structure(a1, {"B": residues(COORDS @ R0.T)})

    row = run_worker(tmp, "ambig", a1, b2)
    check_true("chain: >1 polymer chain without flag is an error",
               row["status"].startswith("error:"))
    check_true("chain: error names the chains and the flag",
               "A, D" in row["status"] and "--chain-b" in row["status"]
               and "structure B" in row["status"])
    check_true("chain: error row has NA rmsd", pd.isna(row["rmsd"]))

    row = run_worker(tmp, "pick", a1, b2, chain_b="D")
    check("chain: explicit --chain-b D", row["status"], "OK")
    check("chain: picked D, not the first chain", float(row["rmsd"]), 0.0, tol=1e-3)
    check("chain: reports chains used", (row["a_chain"], row["b_chain"]), ("B", "D"))

    row = run_worker(tmp, "absent", a1, b2, chain_b="E")
    check_true("chain: absent chain is an error naming the chains present",
               row["status"].startswith("error:") and "chain E" in row["status"]
               and "A, D" in row["status"])

    wateronly = tmp / "wo.pdb"
    write_structure(wateronly, {"A": residues(), "W": []},
                    het={"W": [(1, "HOH", (50, 50, 50))]})
    row = run_worker(tmp, "water_chain", wateronly, b2, chain_b="D")
    check("chain: a water-only chain does not make A ambiguous", row["status"], "OK")


def test_refusals(tmp: Path) -> None:
    b = tmp / "rb.pdb"
    write_structure(b, {"A": residues()})
    short = tmp / "short.pdb"
    write_structure(short, {"A": residues(COORDS[:20], SEQ[:20])})
    row = run_worker(tmp, "short", b, short)
    check_true("refuse: B shorter than max index is an error",
               row["status"].startswith("error:") and "structure B" in row["status"]
               and "25" in row["status"] and "chain A" in row["status"])
    for col in ("rmsd", "max_dev", "seq_match_frac", "n_positions", "worst_seg_rmsd", "a_len"):
        check_true(f"refuse: {col} NA on error row", pd.isna(row[col]))

    noca = residues()
    noca[13] = (*noca[13], False)  # position 14 has no CA
    nca = tmp / "noca.pdb"
    write_structure(nca, {"A": noca})
    row = run_worker(tmp, "noca", nca, b)
    check_true("refuse: missing CA is an error naming position and structure A",
               row["status"].startswith("error:") and "14" in row["status"]
               and "structure A" in row["status"] and "no CA" in row["status"])

    row = run_worker(tmp, "nofile", tmp / "absent.pdb", b)
    check_true("refuse: absent file", row["status"].startswith("error:")
               and "not found" in row["status"])
    row = run_worker(tmp, "badjson", b, b, payload="{not json")
    check_true("refuse: malformed payload", row["status"].startswith("error:"))
    row = run_worker(tmp, "dupjson", b, b, labels=["1:3", "2:4"], segs=[[1, 2, 3], [2, 3, 4]])
    check_true("refuse: duplicated positions in payload",
               row["status"].startswith("error:") and "more than once" in row["status"])


def test_sequence_trust(tmp: Path) -> None:
    b = tmp / "sb.pdb"
    write_structure(b, {"A": residues()})
    mutated = SEQ[:12] + ("K" if SEQ[12] != "K" else "E") + SEQ[13:]  # position 13
    a = tmp / "sa.pdb"
    write_structure(a, {"A": residues(COORDS, mutated)})
    row = run_worker(tmp, "mut", a, b)
    check("seq: status still OK (identity is reported, not gated)", row["status"], "OK")
    check("seq: one mismatch of 12", float(row["seq_match_frac"]), 11 / 12, tol=1e-3)
    a_segs, b_segs = row["a_motif_seq"].split(","), row["b_motif_seq"].split(",")
    check("seq: mismatch visible in segment 1 only",
          [x == y for x, y in zip(a_segs, b_segs)], [True, False, True])

    unk = tmp / "unk.pdb"
    write_structure(unk, {"A": residues(COORDS, "X" * N)})
    row = run_worker(tmp, "unk", unk, unk)
    check("seq: X never counts as a match", float(row["seq_match_frac"]), 0.0, tol=1e-9)
    check("seq: X shown", row["a_motif_seq"].replace(",", ""), "X" * 12)


def test_length_mismatch_warns(tmp: Path) -> None:
    b = tmp / "lb.pdb"
    write_structure(b, {"A": residues()})
    longer = tmp / "longer.pdb"
    extra = np.vstack([COORDS, COORDS[-1:] + [0, 0, 3.8]])
    write_structure(longer, {"A": residues(extra, SEQ + "A")})
    row = run_worker(tmp, "lenwarn", longer, b)
    check_true("length: warn status", row["status"].startswith("warn:")
               and "lengths differ" in row["status"])
    check("length: numbers still computed", float(row["rmsd"]), 0.0, tol=1e-3)


# ---------------------------------------------------------------------------
# Manifest builder, against prosapia's real lineage lookup
# ---------------------------------------------------------------------------


def _run_dir(tmp: Path) -> tuple[Path, dict]:
    rd = tmp / "rundir"
    rd.mkdir()
    s = tmp / "structs"
    s.mkdir()
    paths = {}
    for k in ("p1", "p2", "p4", "p5", "c1", "c2", "c4", "c5", "own5"):
        paths[k] = s / f"{k}.pdb"
        paths[k].write_text("")
    t0 = pd.DataFrame(
        {
            "backbone_path": [str(paths["p1"]), str(paths["p2"]), str(paths["p4"]), str(paths["p5"])],
            "motif_pos": ["3:6,12:16,25:27", None, "18:11,1:5", "1:3,5:7"],
            "mstart": [3, 3, 3, 3],
        },
        index=pd.Index(["p1", "p2", "p4", "p5"], name="name"),
    )
    t0.to_csv(rd / "table0.tsv", sep="\t")
    t1 = pd.DataFrame(
        {
            "boltz_path": [str(paths["c1"]), str(paths["c2"]), None, str(paths["c4"]), str(paths["c5"])],
            "backbone_path": [None, None, None, None, str(paths["own5"])],
            "parent_name": ["p1", "p2", "p1", "p4", "p5"],
            "parent_table": ["table0"] * 5,
            "gen": [1] * 5,
        },
        index=pd.Index(["c1", "c2", "c3", "c4", "c5"], name="name"),
    )
    t1.to_csv(rd / "table1.tsv", sep="\t")
    return rd, paths


def _build(rd: Path, out_dir: Path, df=None, **overrides):
    dm = DataManager(rd)
    df = dm.read_frame("table1") if df is None else df
    fields = dict(
        gpus_per_task=1, col_a="boltz_path", col_b="backbone_path",
        positions="", positions_column="motif_pos", chain_a="", chain_b="",
        output_prefix=None, force=False, input_column="not applicable",
        run_dir=rd,
    )
    fields.update(overrides)
    args = Namespace(**fields)
    ctx = Namespace(df=df, args=args, out_dir=out_dir, lookup=partial(dm.lookup, df),
                    write_meta=lambda **k: None)
    return run_posrmsd.build_posrmsd_manifest(ctx), args  # type: ignore[arg-type]


def test_manifest(tmp: Path) -> None:
    rd, paths = _run_dir(tmp)
    out_dir = rd / "table1" / "posrmsd"
    out_dir.mkdir(parents=True)
    prefix = "boltz_vs_backbone"
    res_dir = out_dir / prefix
    res_dir.mkdir()
    (res_dir / "c1.tsv").write_text("name\tstatus\tsup_path\nc1\tOK\t\n")  # stale

    rows, args = _build(rd, out_dir)
    check("manifest: gpus forced to 0", args.gpus_per_task, 0)
    by_name = {r[0]: r for r in rows}
    check("manifest: submitted rows", sorted(by_name), ["c1", "c5"])
    r = by_name["c1"]
    check("manifest: seven fields", len(r), 7)
    check("manifest: f2 is A (own row)", r[1], str(paths["c1"].resolve()))
    check("manifest: f3 is B resolved from the PARENT", r[2], str(paths["p1"].resolve()))
    check("manifest: empty chains keep their fields", (r[3], r[4]), ("", ""))
    payload = json.loads(r[5])
    check("manifest: segments from the parent's spec", payload["segments"], SEGS)
    check("manifest: labels", payload["labels"], LABELS)
    check("manifest: prefix last", r[6], prefix)
    check("manifest: child's own col-b shadows the parent's",
          by_name["c5"][2], str(paths["own5"].resolve()))
    check_true("manifest: stale result removed for a resubmitted design",
               not (res_dir / "c1.tsv").exists())

    st = pd.read_csv(res_dir / "c2.tsv", sep="\t").iloc[0]["status"]
    check_true("manifest: unresolved positions column -> missing_input row",
               st.startswith("missing_input:") and "--positions-column 'motif_pos'" in st)
    st = pd.read_csv(res_dir / "c4.tsv", sep="\t").iloc[0]["status"]
    check_true("manifest: bad per-design spec -> error row, not an aborted submit",
               st.startswith("error: position spec") and "expands to no positions" in st)
    check_true("manifest: row without --col-a is not touched",
               not (res_dir / "c3.tsv").exists())

    # Resume: an OK row is skipped unless --force.
    dm = DataManager(rd)
    df = dm.read_frame("table1")
    df[f"posrmsd_{prefix}_status"] = ["OK", pd.NA, pd.NA, pd.NA, pd.NA]
    rows, _ = _build(rd, out_dir, df=df)
    check("resume: OK row skipped", sorted(r[0] for r in rows), ["c5"])
    rows, _ = _build(rd, out_dir, df=df, force=True)
    check("resume: --force resubmits", sorted(r[0] for r in rows), ["c1", "c5"])

    # Literal positions, with an {expr} resolved up the lineage per row.
    rows, _ = _build(rd, out_dir, positions_column="", positions="{mstart}:6,10:12",
                     chain_a="B", chain_b="D", output_prefix="lit")
    r = {x[0]: x for x in rows}["c1"]
    check("literal: {expr} resolved from the parent",
          json.loads(r[5])["segments"], [[3, 4, 5, 6], [10, 11, 12]])
    check("literal: chains carried", (r[3], r[4]), ("B", "D"))
    check("literal: c2 now submitted (positions are literal)",
          sorted(x[0] for x in rows), ["c1", "c2", "c4", "c5"])
    check("literal: custom prefix", r[6], "lit")

    # Missing --col-b on one row.
    dm = DataManager(rd)
    df = dm.read_frame("table1")
    t0 = dm.read_frame("table0")
    t0.loc["p1", "backbone_path"] = pd.NA
    t0.to_csv(rd / "table0.tsv", sep="\t")
    rows, _ = _build(rd, out_dir, output_prefix="nob")
    st = pd.read_csv(out_dir / "nob" / "c1.tsv", sep="\t").iloc[0]["status"]
    check_true("col-b: unresolved -> missing_input naming --col-b",
               st.startswith("missing_input:") and "--col-b 'backbone_path'" in st)

    # Submit-time refusals.
    refusals = [
        (dict(positions="1:3"), "mutually exclusive"),
        (dict(positions_column=""), "exactly one of --positions"),
        (dict(input_column="boltz_path"), "-i/--input-column"),
        (dict(chain_a="A,B"), "--chain-a"),
        (dict(chain_b="D1-5"), "--chain-b"),
        (dict(positions_column="", positions="1:3/4:6"), "chain break"),
        (dict(positions_column="", positions="18:11,1:5"), "expands to no positions"),
        (dict(col_a="nope_path"), "not a column of the -t table"),
        (dict(col_b="nope_path"), "has no value on any"),
        (dict(positions_column="nope"), "has no value on any"),
        (dict(output_prefix="bad prefix"), "letters, digits"),
    ]
    for over, needle in refusals:
        check_raises(f"refuse at submit: {over}", lambda o=over: _build(rd, out_dir, **o), needle)


# ---------------------------------------------------------------------------
# Collector
# ---------------------------------------------------------------------------


def test_collector(tmp: Path) -> None:
    out_dir = tmp / "coll" / "posrmsd"
    prefix = "a_vs_b"
    res = out_dir / prefix
    res.mkdir(parents=True)
    b = tmp / "cb.pdb"
    a = tmp / "ca.pdb"
    write_structure(b, {"A": residues()})
    write_structure(a, {"A": residues(COORDS @ R0.T)})
    run_worker(res, "ok", a, b)
    longer = tmp / "clong.pdb"
    write_structure(longer, {"A": residues(np.vstack([COORDS, [[0, 0, 99]]]), SEQ + "A")})
    run_worker(res, "warn", longer, b)
    run_worker(res, "err", tmp / "nope.pdb", b)
    (res / "mi.tsv").write_text("name\tstatus\tsup_path\nmi\tmissing_input: --col-b 'x' has no value\t\n")

    df = pd.DataFrame(index=pd.Index(["ok", "warn", "err", "mi", "gone"], name="name"))
    df[f"posrmsd_{prefix}_seg5_rmsd"] = 1.0  # stale from an earlier collect
    ctx = Namespace(out_dir=out_dir, df=df,
                    args=Namespace(col_a="a_path", col_b="b_path", output_prefix=None))
    one = collect_posrmsd.collect_posrmsd(ctx)  # type: ignore[arg-type]

    def got(name):
        items = list(one(Namespace(name=name, out_dir=out_dir, lookup=None)))  # type: ignore[arg-type]
        check(f"collect {name}: exactly one Collected", len(items), 1)
        check(f"collect {name}: leaf status suppressed", items[0].status, None)
        return items[0].data

    k = lambda c: f"{prefix}_{c}"  # noqa: E731
    d = got("ok")
    check("collect ok: status", d[k("status")], "OK")
    check("collect ok: rmsd float", d[k("rmsd")], 0.0, tol=1e-3)
    check("collect ok: int column", d[k("n_positions")], 12)
    check("collect ok: seg column", d[k("seg2_rmsd")], 0.0, tol=1e-3)
    check_true("collect ok: stale seg5 NA'd", d[k("seg5_rmsd")] is None)
    check_true("collect ok: sup_path", bool(d[k("sup_path")]))
    check("collect ok: chain string", d[k("a_chain")], "A")

    d = got("warn")
    check_true("collect warn: status kept", d[k("status")].startswith("warn:"))
    check_true("collect warn: numbers kept", d[k("rmsd")] is not None)

    for name, start in (("err", "error:"), ("mi", "missing_input:"), ("gone", "missing")):
        d = got(name)
        check_true(f"collect {name}: status {start}", str(d[k("status")]).startswith(start))
        nonnull = [c for c, v in d.items() if v is not None and c != k("status")]
        check(f"collect {name}: every other column NA", nonnull, [])
    check("collect gone: exactly 'missing'", got("gone")[k("status")], "missing")

    check("collector/worker column sets agree",
          set(collect_posrmsd.METRIC_COLUMNS), set(worker.METRIC_COLUMNS))
    for args in (("x_path", "y_path", None), ("x", "y_path", "custom"), ("boltz_path", "rfd3_seg_path", None)):
        check(f"prefix: run == collect for {args}",
              run_posrmsd.resolve_prefix(*args), collect_posrmsd.resolve_prefix(*args))
    check("prefix: default", run_posrmsd.resolve_prefix("boltz_path", "chainsel_binder_path", None),
          "boltz_vs_chainsel_binder")
    check_raises("collect: wrong prefix -> results dir not found",
                 lambda: collect_posrmsd.collect_posrmsd(Namespace(
                     out_dir=out_dir, df=df,
                     args=Namespace(col_a=None, col_b=None, output_prefix="other"))),  # type: ignore[arg-type]
                 "Results dir not found", FileNotFoundError)


# ---------------------------------------------------------------------------
# The task script, through the real prosapia prelude
# ---------------------------------------------------------------------------


def test_task_script(tmp: Path) -> None:
    import prosapia.core as core

    prelude = Path(core.__file__).parent / "scripts" / "sapia_task_prelude.sh"
    work = tmp / "task"
    work.mkdir()
    b = work / "tb.pdb"
    a = work / "ta.pdb"
    write_structure(b, {"A": residues(), "D": residues(COORDS + 50)})
    write_structure(a, {"A": residues(COORDS @ R0.T)})
    out_dir = work / "out"
    out_dir.mkdir()
    manifest = work / "manifest.txt"
    rows = [
        ("t1", str(a), str(b), "", "A", segments_json(LABELS, SEGS), "pfx"),
        ("t2", str(a), str(b), "", "", segments_json(LABELS, SEGS), "pfx"),
    ]
    manifest.write_text("".join("\t".join(r) + "\n" for r in rows))
    for task_id in ("1", "2"):
        env = dict(os.environ, SAPIA_PRELUDE=str(prelude), SAPIA_SCHEDULER="modal",
                   SAPIA_TASK_ID=task_id, SAPIA_TOOL_DIR=str(HERE),
                   PIPELINE_PYTHON=sys.executable)
        proc = subprocess.run(["bash", str(HERE / "posrmsd.sh"), str(manifest), str(out_dir)],
                              cwd=work, env=env, capture_output=True, text=True)
        check(f"task {task_id}: exit 0", proc.returncode, 0)
        if proc.returncode != 0:
            FAILURES.append(f"task {task_id} stderr: {proc.stderr[-400:]}")
    read = lambda n: pd.read_csv(out_dir / "pfx" / f"{n}.tsv", sep="\t", dtype=str,  # noqa: E731
                                 keep_default_na=False, na_values=[""]).iloc[0]
    r1 = read("t1")
    check("task t1: status (chain-b field reached the worker)", r1["status"], "OK")
    check("task t1: rmsd", float(r1["rmsd"]), 0.0, tol=1e-3)
    check("task t1: b_chain", r1["b_chain"], "A")
    check_true("task t1: sup written", (out_dir / "pfx" / "t1_sup.cif").is_file())
    r2 = read("t2")
    check_true("task t2: empty chain-b field -> ambiguity error, not a guess",
               r2["status"].startswith("error:") and "--chain-b" in r2["status"])


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        test_parse_segments()
        test_kabsch()
        test_perfect_and_frame(tmp)
        test_global_fit_and_independent_check(tmp)
        test_file_order_not_resnum(tmp)
        test_hetero_residues_do_not_count(tmp)
        test_chain_selection(tmp)
        test_refusals(tmp)
        test_sequence_trust(tmp)
        test_length_mismatch_warns(tmp)
        test_manifest(tmp)
        test_collector(tmp)
        test_task_script(tmp)

    if FAILURES:
        print(f"FAILED {len(FAILURES)} of {CHECKS} checks:")
        for f in FAILURES:
            print(f"  - {f}")
        return 1
    print(f"OK: {CHECKS} checks passed.")
    print(
        "NOTE: this proves the arithmetic and the contracts only. It cannot prove "
        "that real A/B pairs share residue ordering, which chain letters real files "
        "carry, or the numbering of a real positions column -- see the docstring."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
