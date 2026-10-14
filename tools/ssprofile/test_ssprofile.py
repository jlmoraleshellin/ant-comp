#!/usr/bin/env python3
"""
Self-test for ssprofile: synthetic backbones with a hand-computable or independently
computed answer.

WHAT THIS PROVES, AND WHAT IT CANNOT.

It proves the arithmetic and the contracts: Rg against hand-computed values (and its
invariance to rigid motion and to residue identity, i.e. that it is unweighted), that
the SS fractions are the counts of the assignment divided by n_res, that the SS call
reproduces willutil's ``wu.dssp`` file round-trip character for character on the
fixtures, that an isolated strand reads as loop (the documented in-isolation
blindness), that every refusal path writes a non-OK status with NA rather than a zero,
that the builder refuses a missing input column / bad chain spec and submits
unresolvable rows as errors, that the collector re-checks the worker's contract, and
that the .sh task script cuts its manifest fields correctly and writes fallback rows
when the worker crashes.

It CANNOT prove the premise or the mapping. Every structure here is written by this
file, so it agrees with this file by construction. UNKNOWN until real data:

  * whether real rfd3 / boltz / bindcraft2 / boltzgen outputs put the binder on the
    chain letter callers will pass to --design-chains (n_res is the cross-check);
  * whether those outputs carry an O atom on every residue (a missing O, including on
    the C-terminal residue, is a 'warn:' here, never guessed);
  * whether those outputs have in-chain C-N gaps > 2.5 A (a 'warn:' here);
  * whether the mdtraj in rpxdock's image (unpinned there) is the same version as the
    1.11.0 pinned here -- parity with rpxdock's frac_* is "same code path", tested here
    against a reimplementation of willutil's wrapper, not against rpxdock itself;
  * whether DSSP (any implementation) is right on generated backbones -- nobody has
    cross-checked it against an independent SS assignment on this repo's data.

Run:  uv run --with mdtraj==1.11.0 python tools/ssprofile/test_ssprofile.py
"""

import importlib.util
import os
import subprocess
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


worker = _load("ssprofile_worker")
run_ssprofile = _load("run_ssprofile")
collect_ssprofile = _load("collect_ssprofile")

from prosapia.core import CollectCtx, DesignCtx, ManifestCtx  # noqa: E402

FAILURES: list[str] = []
CHECKS = 0


def check(cond: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not cond:
        FAILURES.append(label)


def close(a, b, tol=1e-6) -> bool:
    return a is not None and b is not None and abs(float(a) - float(b)) <= tol


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


def _place(a, b, c, bond, angle, torsion):
    """NeRF: the atom bonded to c at `bond`, angle b-c-d, torsion a-b-c-d (deg)."""
    angle, torsion = np.radians(angle), np.radians(torsion)
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d = np.array(
        [-bond * np.cos(angle), bond * np.sin(angle) * np.cos(torsion),
         bond * np.sin(angle) * np.sin(torsion)]
    )
    return c + d[0] * bc + d[1] * m + d[2] * n


def backbone(n: int, phi: float, psi: float) -> np.ndarray:
    """(n, 4, 3) ideal N/CA/C/O backbone with constant phi/psi, trans omega."""
    N = np.array([0.0, 0.0, 0.0])
    CA = np.array([1.458, 0.0, 0.0])
    t = np.radians(180 - 111.2)
    C = CA + 1.525 * np.array([np.cos(t), np.sin(t), 0.0])
    out = []
    for i in range(n):
        N_next = _place(N, CA, C, 1.329, 116.2, psi)
        O = _place(N, CA, C, 1.231, 120.5, psi + 180)
        out.append(np.array([N, CA, C, O]))
        if i == n - 1:
            break
        CA_next = _place(CA, C, N_next, 1.458, 121.7, 180.0)
        C_next = _place(C, N_next, CA_next, 1.525, 111.2, phi)
        N, CA, C = N_next, CA_next, C_next
    return np.array(out)


def _axis(v: np.ndarray) -> np.ndarray:
    """Unit vector with a deterministic sign (largest component positive), so the
    sheet fixture does not depend on SVD sign conventions."""
    v = v / np.linalg.norm(v)
    return v if v[np.argmax(np.abs(v))] > 0 else -v


def _rot(axis: np.ndarray, deg: float) -> np.ndarray:
    a = np.radians(deg)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(a) * k + (1 - np.cos(a)) * k @ k


def sheet(n: int = 12) -> tuple[np.ndarray, np.ndarray]:
    """Two antiparallel extended strands that mdtraj's DSSP pairs into a sheet
    (geometry found by a grid search during development; asserted below)."""
    s1 = backbone(n, -139.0, 135.0)
    ca = s1[:, 1]
    cen = ca.mean(axis=0)
    _, _, vt = np.linalg.svd(ca - cen)
    ax, p1, p2 = _axis(vt[0]), _axis(vt[1]), _axis(vt[2])
    flip = _rot(p1, 180.0)
    turn = _rot(ax, 150.0)
    s2 = ((s1 - cen) @ flip.T) @ turn.T + cen + 4.2 * p2 - 1.0 * ax
    return s1, s2


def write_pdb(path: Path, chains: dict[str, np.ndarray], resnames=None,
              drop: set[tuple[str, int, str]] = frozenset(), hetero: dict | None = None,
              start: int = 1) -> Path:
    """chains: letter -> (n, 4, 3). drop: {(chain, index, atom)} atoms to omit.
    hetero: letter -> list of (resname, atom name, element, xyz) HETATM residues."""
    st = gemmi.Structure()
    model = gemmi.Model("1")
    for letter, bb in chains.items():
        ch = gemmi.Chain(letter)
        for i, res_xyz in enumerate(bb):
            res = gemmi.Residue()
            res.name = (resnames or {}).get(letter, ["GLY"] * len(bb))[i]
            res.seqid = gemmi.SeqId(start + i, " ")
            res.entity_type = gemmi.EntityType.Polymer
            res.het_flag = "A"
            for name, el, xyz in zip(("N", "CA", "C", "O"), "NCCO", res_xyz):
                if (letter, i, name) in drop:
                    continue
                atom = gemmi.Atom()
                atom.name = name
                atom.element = gemmi.Element(el)
                atom.pos = gemmi.Position(*xyz)
                atom.occ = 1.0
                res.add_atom(atom)
            ch.add_residue(res)
        model.add_chain(ch)
    for letter, items in (hetero or {}).items():
        ch = gemmi.Chain(letter)
        for k, (rname, aname, el, xyz) in enumerate(items):
            res = gemmi.Residue()
            res.name = rname
            res.seqid = gemmi.SeqId(900 + k, " ")
            res.het_flag = "H"
            atom = gemmi.Atom()
            atom.name = aname
            atom.element = gemmi.Element(el)
            atom.pos = gemmi.Position(*xyz)
            atom.occ = 1.0
            res.add_atom(atom)
            ch.add_residue(res)
        model.add_chain(ch)
    st.add_model(model)
    st.setup_entities()
    path.parent.mkdir(parents=True, exist_ok=True)
    st.write_pdb(str(path))
    return path


def wu_dssp_roundtrip(ncaco: np.ndarray) -> str:
    """willutil's wu.dssp, re-implemented from its source (ded19f6): dump N/CA/C/O as
    GLY into ONE chain with willutil's own record format, mdtraj.load, simplified
    DSSP, 'C' -> 'L'. The independent path the worker's in-memory call must match."""
    import mdtraj as md

    fmt = (
        "ATOM  {ia:5d} {an:^4}{idx:^1}{rn:3s} {c:1}{ir:4d}{insert:1s}   "
        "{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}           {elem:1s}\n"
    )
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "tmp.pdb"
        ia = 0
        with open(path, "w") as f:
            for i, p in enumerate(ncaco):
                for j, an in enumerate(("N", "CA", "C", "O")):
                    f.write(fmt.format(ia=ia + j, an=an, idx=" ", rn="GLY", c="A", ir=i,
                                       insert=" ", x=p[j, 0], y=p[j, 1], z=p[j, 2],
                                       occ=1, b=1, elem=" "))
                ia += 4
        t = md.load(str(path))
        return "".join(md.compute_dssp(t)[0]).replace("C", "L")


def run_worker(tmp: Path, name: str, src: str, chains: str) -> pd.Series:
    out_dir = tmp / f"out_{name}"
    task = tmp / f"task_{name}.tsv"
    task.write_text(f"{name}\t{src}\n")
    sys.argv = ["ssprofile_worker.py", "--task-file", str(task), "--design-chains",
                chains, "--input-column", "rfdiffusion3_path", "--out-dir", str(out_dir)]
    worker.main()
    return pd.read_csv(out_dir / f"{name}.tsv", sep="\t", dtype=str).iloc[0]


def num(row: pd.Series, col: str):
    v = row.get(col)
    if v is None or pd.isna(v) or str(v).strip() == "":
        return None
    return float(v)


def all_na(row: pd.Series, cols) -> bool:
    return all(num(row, c) is None for c in cols)


NUMERIC = ["radius_of_gyration", "helix_fraction", "beta_fraction", "loop_fraction",
           "n_res", "n_chain_breaks"]
FRACS = ["helix_fraction", "beta_fraction", "loop_fraction"]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_rg_arithmetic() -> None:
    # Six points at distance 3 on the axes around (10, -4, 7): Rg = 3 exactly.
    pts = np.array([[3, 0, 0], [-3, 0, 0], [0, 3, 0], [0, -3, 0], [0, 0, 3], [0, 0, -3]],
                   dtype=float) + np.array([10.0, -4.0, 7.0])
    check(close(worker.radius_of_gyration(pts), 3.0, 1e-12), "Rg of an octahedron of radius 3 is 3")
    # Two points 10 apart: Rg = 5.
    check(close(worker.radius_of_gyration(np.array([[0, 0, 0], [10.0, 0, 0]])), 5.0, 1e-12),
          "Rg of two points 10 A apart is 5")
    # One point: 0 (a measurement, n_res says why it is meaningless).
    check(close(worker.radius_of_gyration(np.array([[1.0, 2, 3]])), 0.0, 1e-12), "Rg of one point is 0")
    # Rigid motion invariance.
    rng = np.random.default_rng(0)
    x = rng.normal(size=(50, 3)) * 7
    r = _rot(_axis(np.array([1.0, 2.0, 3.0])), 73.0)
    check(close(worker.radius_of_gyration(x), worker.radius_of_gyration(x @ r.T + 40.0), 1e-9),
          "Rg invariant under rotation + translation")
    # Matches the explicit formula.
    want = np.sqrt(((x - x.mean(0)) ** 2).sum(1).mean())
    check(close(worker.radius_of_gyration(x), want, 1e-12), "Rg equals sqrt(mean |r - centroid|^2)")
    try:
        worker.radius_of_gyration(np.empty((0, 3)))
        check(False, "Rg of an empty set must raise")
    except ValueError:
        check(True, "")


def test_helix(tmp: Path) -> None:
    bb = backbone(30, -57.0, -47.0)
    pdb = write_pdb(tmp / "helix.pdb", {"A": bb})
    row = run_worker(tmp, "helix", str(pdb), "A")
    check(row["status"] == "OK", f"helix OK ({row['status']})")
    ss = worker.assign_ss([worker.select(worker.read_model(pdb), ["A"], "helix.pdb")[0]["A"]])
    check(ss == list(wu_dssp_roundtrip(bb)), "helix: in-memory SS equals willutil's file round-trip")
    check(close(num(row, "helix_fraction"), round(ss.count("H") / 30, 6)), "helix_fraction = H/n_res")
    check(num(row, "helix_fraction") > 0.8, f"ideal helix is mostly H ({num(row, 'helix_fraction')})")
    check(close(sum(num(row, c) for c in FRACS), 1.0, 1e-5), "fractions sum to 1")
    check(num(row, "n_res") == 30, "n_res = 30")
    check(num(row, "n_chain_breaks") == 0, "n_chain_breaks = 0 on a continuous chain")
    check(row["chains_used"] == "A", "chains_used echoes A")
    check(row["sequence"] == "G" * 30, f"sequence is the selection's one-letter sequence ({row['sequence']})")
    check(row["input_column"] == "rfdiffusion3_path" and row["source_path"] == str(pdb),
          "provenance: input_column and source_path recorded on an OK row")
    ca = bb[:, 1]
    check(close(num(row, "radius_of_gyration"),
                np.sqrt(((ca - ca.mean(0)) ** 2).sum(1).mean()), 1e-3),
          "Rg column = CA Rg of the fixture")
    per_res = pd.read_csv(row["path"], sep="\t", dtype=str)
    check(len(per_res) == 30 and "".join(per_res["ss"]) == "".join(ss),
          "per-residue file carries the same assignment")


def test_residue_identity_does_not_matter(tmp: Path) -> None:
    """Unweighted Rg and GLY-only SS: renaming every residue changes nothing."""
    bb = backbone(20, -57.0, -47.0)
    names = (["TRP", "PRO", "ALA", "LYS", "GLU"] * 4)[:20]
    a = run_worker(tmp, "gly", str(write_pdb(tmp / "gly.pdb", {"A": bb})), "A")
    b = run_worker(tmp, "named", str(write_pdb(tmp / "named.pdb", {"A": bb},
                                               resnames={"A": names})), "A")
    for c in ["radius_of_gyration", *FRACS]:
        check(a[c] == b[c], f"{c} identical for GLY and named residues ({a[c]} vs {b[c]})")
    check(b["sequence"] == "WPAKE" * 4, f"sequence follows residue names ({b['sequence']})")
    # Modified and unknown residues stay visible in the sequence.
    c = run_worker(tmp, "odd", str(write_pdb(tmp / "odd.pdb", {"A": bb},
                   resnames={"A": ["MSE", "UNK"] + ["ALA"] * 18})), "A")
    check(c["sequence"] == "mX" + "A" * 18, f"MSE -> 'm', UNK -> 'X' ({c['sequence']})")


def test_sheet_and_isolation(tmp: Path) -> None:
    s1, s2 = sheet()
    pdb = write_pdb(tmp / "sheet.pdb", {"A": s1, "B": s2})
    # The worker takes the builder's EXPANDED list; an unexpanded range reaching it
    # is an absent chain, loudly.
    raw = run_worker(tmp, "sheet_raw", str(pdb), "A:B")
    check(raw["status"].startswith("error:"), "worker given an unexpanded 'A:B' errors, never guesses")
    both = run_worker(tmp, "sheet_ab", str(pdb), ",".join(run_ssprofile.parse_design_chains("A:B")))
    check(both["status"] == "OK", f"sheet OK ({both['status']})")
    check(both["chains_used"] == "A,B", f"A:B expands to A,B ({both['chains_used']})")
    check(num(both, "n_res") == 24, "n_res counts both chains")
    check(both["sequence"] == "G" * 12 + "/" + "G" * 12, f"multi-chain sequence joined by '/' ({both['sequence']})")
    model = worker.read_model(pdb)
    sel, _ = worker.select(model, ["A", "B"], "sheet.pdb")
    ss = worker.assign_ss([sel["A"], sel["B"]])
    check(ss.count("E") >= 10, f"fixture really is a sheet ({''.join(ss)})")
    check(close(num(both, "beta_fraction"), round(ss.count("E") / 24, 6)), "beta_fraction = E/n_res")
    allca = np.concatenate([s1[:, 1], s2[:, 1]])
    check(close(num(both, "radius_of_gyration"),
                np.sqrt(((allca - allca.mean(0)) ** 2).sum(1).mean()), 1e-3),
          "multi-chain Rg is the Rg of the union")
    # The documented blindness: strand A alone has no partner -> loop.
    alone = run_worker(tmp, "sheet_a", str(pdb), "A")
    check(alone["status"] == "OK" and num(alone, "beta_fraction") == 0.0,
          f"a strand whose partner is not selected reads as loop "
          f"(beta={alone['beta_fraction']})")
    check(num(alone, "loop_fraction") == 1.0, "...and all of it is loop")
    # Single-chain parity with willutil on a non-helical chain too.
    check(worker.assign_ss([sel["A"]]) == list(wu_dssp_roundtrip(s1)),
          "strand: in-memory SS equals willutil's file round-trip")


def test_error_contracts(tmp: Path) -> None:
    bb = backbone(12, -57.0, -47.0)
    pdb = write_pdb(tmp / "err.pdb", {"A": bb},
                    hetero={"L": [("HEM", "FE", "Fe", (30.0, 0.0, 0.0))]})

    row = run_worker(tmp, "absent", str(pdb), "D")
    check(row["status"].startswith("error:") and "absent" in row["status"],
          f"absent chain -> error ({row['status']})")
    check(all_na(row, NUMERIC) and pd.isna(row["chains_used"]) and pd.isna(row["sequence"]),
          "absent chain -> every measured column NA")
    check(row["input_column"] == "rfdiffusion3_path" and row["source_path"] == str(pdb),
          "absent chain -> provenance still recorded")

    row = run_worker(tmp, "partial_absent", str(pdb), "A,D")
    check(row["status"].startswith("error:"), "one absent chain of two -> error, not a smaller selection")
    check(all_na(row, NUMERIC), "...with every column NA")

    row = run_worker(tmp, "hetonly", str(pdb), "L")
    check(row["status"].startswith("error:") and "no amino-acid" in row["status"],
          f"ligand-only chain -> error ({row['status']})")
    check(all_na(row, NUMERIC), "ligand-only chain -> NA")

    row = run_worker(tmp, "mixed_empty", str(pdb), "A,L")
    check(row["status"].startswith("error:") and "L hold no" in row["status"],
          f"one selected chain with no CA -> error, not a smaller selection ({row['status']})")
    check(all_na(row, NUMERIC), "...with every column NA")

    row = run_worker(tmp, "nofile", str(tmp / "does_not_exist.pdb"), "A")
    check(row["status"].startswith("error:") and "missing" in row["status"], "missing file -> error")

    garbage = tmp / "garbage.pdb"
    garbage.write_text("this is not a structure\n")
    row = run_worker(tmp, "garbage", str(garbage), "A")
    check(row["status"].startswith("error:"), f"unparseable file -> error ({row['status']})")
    check(all_na(row, NUMERIC), "unparseable file -> NA")

    row = run_worker(tmp, "unresolved", "", "A")
    check(row["status"].startswith("error:") and "could not resolve" in row["status"]
          and "rfdiffusion3_path" in row["status"], "empty path -> error naming the column")


def test_calcium_is_not_a_residue(tmp: Path) -> None:
    """A Ca2+ ion is an atom NAMED CA; it must not enter the CA trace."""
    bb = backbone(12, -57.0, -47.0)
    pdb = write_pdb(tmp / "ca_ion.pdb", {"A": bb},
                    hetero={"A": [("CA", "CA", "Ca", (60.0, 60.0, 60.0))]})
    row = run_worker(tmp, "ca_ion", str(pdb), "A")
    check(row["status"] == "OK" and num(row, "n_res") == 12,
          f"calcium excluded: n_res 12 ({row['n_res']}, {row['status']})")
    ca = bb[:, 1]
    check(close(num(row, "radius_of_gyration"),
                np.sqrt(((ca - ca.mean(0)) ** 2).sum(1).mean()), 1e-3),
          "Rg not dragged by a far-away calcium")


def test_warn_contracts(tmp: Path) -> None:
    bb = backbone(16, -57.0, -47.0)
    pdb = write_pdb(tmp / "noO.pdb", {"A": bb}, drop={("A", 5, "O")})
    row = run_worker(tmp, "noO", str(pdb), "A")
    check(row["status"].startswith("warn:") and "missing backbone" in row["status"],
          f"missing O -> warn ({row['status']})")
    check(all_na(row, FRACS), "missing O -> fractions NA (never guessed, never 0)")
    check(num(row, "radius_of_gyration") is not None and num(row, "n_res") == 16,
          "missing O -> Rg and n_res still measured")
    check(num(row, "n_chain_breaks") is None, "missing O -> n_chain_breaks NA")

    # The C-terminal O too: strict.
    pdb = write_pdb(tmp / "noOterm.pdb", {"A": bb}, drop={("A", 15, "O")})
    row = run_worker(tmp, "noOterm", str(pdb), "A")
    check(row["status"].startswith("warn:"), "missing terminal O -> warn")

    # A residue with no CA at all is not in the trace; its neighbours still are.
    pdb = write_pdb(tmp / "noCA.pdb", {"A": bb}, drop={("A", 3, "CA")})
    row = run_worker(tmp, "noCA", str(pdb), "A")
    check(num(row, "n_res") == 15, f"a residue without CA is not counted ({row['n_res']})")
    check(row["status"].startswith("warn:") and "A4(CA)" in row["status"],
          f"an interior CA-less residue -> missing-backbone warn naming it ({row['status']})")
    pdb = write_pdb(tmp / "noCAterm.pdb", {"A": bb}, drop={("A", 15, "CA")})
    row = run_worker(tmp, "noCAterm", str(pdb), "A")
    check(row["status"].startswith("warn:") and "A16(CA)" in row["status"] and num(row, "n_res") == 15,
          f"a CA-less residue at the chain END is flagged, not silently dropped ({row['status']})")
    check(all_na(row, FRACS), "...with fractions NA")

    # A break: second half moved 10 A away.
    broken = bb.copy()
    broken[8:] += np.array([10.0, 0.0, 0.0])
    pdb = write_pdb(tmp / "broken.pdb", {"A": broken})
    row = run_worker(tmp, "broken", str(pdb), "A")
    check(row["status"].startswith("warn:") and "break" in row["status"],
          f"chain break -> warn ({row['status']})")
    check(num(row, "n_chain_breaks") == 1, f"n_chain_breaks = 1 ({row['n_chain_breaks']})")
    check(all_na(row, FRACS), "chain break -> fractions NA")
    check(num(row, "radius_of_gyration") is not None, "chain break -> Rg measured")

    # Two chains whose ends are far apart are NOT a break.
    pdb = write_pdb(tmp / "twochains.pdb", {"A": bb[:8], "B": bb[8:] + 30.0})
    row = run_worker(tmp, "twochains", str(pdb), "A,B")
    check(row["status"] == "OK" and num(row, "n_chain_breaks") == 0,
          f"a gap BETWEEN selected chains is not a break ({row['status']})")


def build_manifest(tmp: Path, df: pd.DataFrame, **kw):
    out_dir = tmp / "manifest_out" / kw.pop("leaf", "ssprofile")
    out_dir.mkdir(parents=True, exist_ok=True)
    args = Namespace(
        run_dir=tmp, table="table1", input_column=kw.pop("input_column", "rfdiffusion3_path"),
        force=kw.pop("force", False), gpus_per_task=1,
        design_chains=kw.pop("design_chains", "B"),
        designs_per_task=kw.pop("designs_per_task", 20),
    )
    lineage = kw.pop("lineage", {})
    ctx = ManifestCtx(df=df, args=args, out_dir=out_dir,
                      lookup=lambda name, column: lineage.get((name, column)),
                      write_meta=lambda **_: None)
    rows = run_ssprofile.build_ssprofile_manifest(ctx)
    sub = []
    for row in rows:
        sub.extend(l.split("\t") for l in Path(row[0]).read_text().splitlines() if l.strip())
    return rows, sub, args


def raises(fn, needle: str) -> bool:
    try:
        fn()
    except ValueError as e:
        return needle in str(e)
    return False


def test_builder(tmp: Path) -> None:
    pdb = write_pdb(tmp / "b.pdb", {"A": backbone(5, -57, -47), "B": backbone(5, -57, -47)})
    names = [f"d{i}" for i in range(5)]
    df = pd.DataFrame({"rfdiffusion3_path": [str(pdb)] * 5}, index=pd.Index(names, name="name"))

    rows, sub, args = build_manifest(tmp, df, design_chains="A:B", designs_per_task=2)
    check(args.gpus_per_task == 0, "builder forces gpus_per_task = 0")
    check(len(rows) == 3 and len(sub) == 5, f"5 designs in 3 tasks of <=2 ({len(rows)}, {len(sub)})")
    check(all(r[1] == "A,B" and r[2] == "rfdiffusion3_path" for r in rows),
          "manifest carries expanded chains and the input column")

    check(raises(lambda: build_manifest(tmp, df, input_column="not applicable"),
                 "no default input column"), "sentinel input column raises")
    check(raises(lambda: build_manifest(tmp, df, design_chains=""), "at least one chain"),
          "empty --design-chains raises")
    check(raises(lambda: build_manifest(tmp, df, design_chains="A,A"), "twice"),
          "duplicated chain raises")
    check(raises(lambda: build_manifest(tmp, df, design_chains="A:"), ""),
          "malformed chain range raises")

    # Lineage: column absent from this table, present on the parent.
    child = pd.DataFrame({"atomium_sequence": ["AAAA", "CCCC"]}, index=pd.Index(["c0", "c1"], name="name"))
    rows, sub, _ = build_manifest(tmp, child, lineage={("c0", "rfdiffusion3_path"): str(pdb)})
    got = {s[0]: s[1] for s in sub}
    check(got.get("c0", "").endswith("b.pdb"), "lineage-resolved path is used")
    check(got.get("c1") == "", "unresolvable row is submitted with an EMPTY path, not dropped")
    check(raises(lambda: build_manifest(tmp, child), "resolved no structure"),
          "all rows unresolvable -> the submit raises")

    # Resume: an OK row is skipped unless --force.
    done = df.copy()
    done["ssprofile_bb_status"] = ["OK", "error: x", "OK", "warn: y", None]
    _, sub, _ = build_manifest(tmp, done, leaf="ssprofile_bb")
    check(sorted(s[0] for s in sub) == ["d1", "d3", "d4"], f"OK rows skipped on rerun ({sub})")
    _, sub, _ = build_manifest(tmp, done, leaf="ssprofile_bb", force=True)
    check(len(sub) == 5, "--force resubmits OK rows")

    # mmCIF input is staged to PDB.
    st = gemmi.read_structure(str(pdb))
    cif = tmp / "b_cif.cif"
    st.make_mmcif_document().write_file(str(cif))
    cdf = pd.DataFrame({"rfdiffusion3_path": [str(cif)]}, index=pd.Index(["cx"], name="name"))
    _, sub, _ = build_manifest(tmp, cdf)
    check(sub[0][1].endswith(".pdb") and Path(sub[0][1]).exists(), f"CIF staged to PDB ({sub[0][1]})")


def collect_one(tmp: Path, name: str, tsv_text: str | None) -> dict:
    out_dir = tmp / "collect_out"
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv = out_dir / f"{name}.tsv"
    tsv.unlink(missing_ok=True)
    if tsv_text is not None:
        tsv.write_text(tsv_text)
    ctx = CollectCtx(
        df=pd.DataFrame(index=pd.Index([name])),
        args=Namespace(run_dir=tmp, table="table1", dir_label="", force=False),
        table_name="table1", out_dir=out_dir, status_col="ssprofile_status",
        path_col="ssprofile_path", parent_table=None, parent_df=pd.DataFrame(),
        lookup=lambda n, c: None, creates_table=False, default_input_column="not applicable",
    )
    emitted = list(collect_ssprofile.collect_ssprofile(ctx)(DesignCtx(name, out_dir, ctx.lookup)))
    assert len(emitted) == 1
    c = emitted[0]
    return {"status": c.status, "path": c.path, "data": dict(c.data)}


HEADER = "name\tstatus\tpath\t" + "\t".join(worker.METRIC_COLUMNS) + "\n"


def test_collector(tmp: Path) -> None:
    got = collect_one(tmp, "ok", HEADER + "ok\tOK\t/x.tsv\t12.5\t0.5\t0.25\t0.25\t40\t0\tB\n")
    check(got["status"] == "OK" and got["data"]["radius_of_gyration"] == 12.5
          and got["data"]["n_res"] == 40 and got["data"]["chains_used"] == "B"
          and got["data"]["beta_fraction"] == 0.25, f"OK row collects ({got})")
    check(set(got["data"]) == {"radius_of_gyration", "helix_fraction", "beta_fraction",
                               "loop_fraction", "n_res", "n_chain_breaks", "chains_used",
                               "sequence", "input_column", "source_path"},
          "collected column set")
    got = collect_one(tmp, "prov", HEADER + "prov\tOK\t/x.tsv\t12.5\t0.5\t0.25\t0.25\t4\t0\tB\tMKLV\tboltz_path\t/runs/a.pdb\n")
    check(got["data"]["sequence"] == "MKLV" and got["data"]["input_column"] == "boltz_path"
          and got["data"]["source_path"] == "/runs/a.pdb", f"new string columns collect ({got['data']})")
    got = collect_one(tmp, "errprov", HEADER + "errprov\terror: chain(s) D absent\t\t\t\t\t\t\t\t\t\tboltz_path\t/runs/a.pdb\n")
    check(got["status"].startswith("error:") and got["data"]["input_column"] == "boltz_path"
          and got["data"]["source_path"] == "/runs/a.pdb"
          and all(got["data"][c] is None for c in got["data"] if c not in ("input_column", "source_path")),
          f"error row keeps provenance, every measured column NA ({got['data']})")

    got = collect_one(tmp, "zero", HEADER + "zero\tOK\t/x.tsv\t9\t1.0\t0.0\t0.0\t30\t0\tA\n")
    check(got["data"]["beta_fraction"] == 0.0, "a measured 0.0 stays 0.0")

    got = collect_one(tmp, "warn", HEADER + "warn\twarn: 1 backbone break(s)\t/x.tsv\t9\t\t\t\t30\t1\tA\n")
    check(got["status"].startswith("warn:") and got["data"]["radius_of_gyration"] == 9.0
          and all(got["data"][c] is None for c in FRACS), f"warn row: Rg kept, fractions NA ({got})")

    got = collect_one(tmp, "err", HEADER + "err\terror: chain(s) D absent\t\t\t\t\t\t\t\t\n")
    check(got["status"].startswith("error:") and all(v is None for v in got["data"].values()),
          "error row: every column NA")

    got = collect_one(tmp, "crash", "name\tstatus\tpath\ncrash\terror: worker crashed\t\n")
    check(got["status"] == "error: worker crashed" and all(v is None for v in got["data"].values()),
          "crash fallback TSV (3 columns) collects as an error with NA")

    got = collect_one(tmp, "absent", None)
    check(got["status"] == "missing", "no TSV -> missing")

    got = collect_one(tmp, "badsum", HEADER + "badsum\tOK\t/x\t9\t0.5\t0.5\t0.5\t30\t0\tA\n")
    check(got["status"].startswith("error:") and got["data"]["radius_of_gyration"] is None,
          f"OK row whose fractions do not sum to 1 -> error, NA ({got['status']})")

    got = collect_one(tmp, "warnfrac", HEADER + "warnfrac\twarn: x\t/x\t9\t1\t0\t0\t30\t\tA\n")
    check(got["status"].startswith("error:"), "warn row carrying fractions -> error")

    got = collect_one(tmp, "okmissing", HEADER + "okmissing\tOK\t/x\t9\t\t\t\t30\t0\tA\n")
    check(got["status"].startswith("error:"), "OK row with missing fractions -> error")


def test_task_script(tmp: Path) -> None:
    """Run ssprofile.sh under a stub prelude: field cutting, and the crash fallback."""
    pdb = write_pdb(tmp / "sh.pdb", {"B": backbone(14, -57, -47)})
    out_dir = tmp / "sh_out"
    out_dir.mkdir()
    task = tmp / "sh_task.tsv"
    task.write_text(f"s0\t{pdb}\ns1\t\n")
    prelude = tmp / "prelude.sh"
    prelude.write_text(
        'OUT_DIR="$SH_OUT"; SAPIA_TASK_ID=1; SAPIA_LINE="$SH_LINE"\n'
        "sapia_activate() { :; }\n"
    )
    env = dict(os.environ, SAPIA_PRELUDE=str(prelude), SAPIA_TOOL_DIR=str(HERE),
               PIPELINE_PYTHON=sys.executable, SH_OUT=str(out_dir),
               SH_LINE=f"{task}\tB\trfdiffusion3_path")
    # A stale residues file from an "earlier run" for the unresolved design must not
    # survive this run.
    (out_dir / "s1_residues.tsv").write_text("chain\tresnum\tresname\tss\nB\t1\tGLY\tH\n")
    p = subprocess.run(["bash", str(HERE / "ssprofile.sh")], env=env, capture_output=True, text=True)
    check(not (out_dir / "s1_residues.tsv").exists(), "stale <name>_residues.tsv removed before the worker runs")
    check(p.returncode == 0, f"task script exits 0 ({p.returncode}: {p.stderr[-300:]})")
    s0 = pd.read_csv(out_dir / "s0.tsv", sep="\t", dtype=str).iloc[0]
    s1 = pd.read_csv(out_dir / "s1.tsv", sep="\t", dtype=str).iloc[0]
    check(s0["status"] == "OK" and s0["chains_used"] == "B", f"task script: s0 OK ({s0['status']})")
    check(s1["status"].startswith("error:") and "rfdiffusion3_path" in s1["status"],
          "task script: unresolved row -> error naming the column")

    # Crash: a worker that cannot import its engine -> every design gets an error row,
    # and the task exits non-zero.
    crash_out = tmp / "sh_crash"
    crash_out.mkdir()
    fake = tmp / "fakepy"
    fake.write_text("#!/bin/bash\nexit 3\n")
    fake.chmod(0o755)
    env.update(SH_OUT=str(crash_out), PIPELINE_PYTHON=str(fake))
    # RISK A: an OK result left by an EARLIER run must not survive a crash of this one.
    (crash_out / "s0.tsv").write_text((out_dir / "s0.tsv").read_text())
    (crash_out / "s0_residues.tsv").write_text("stale\n")
    p = subprocess.run(["bash", str(HERE / "ssprofile.sh")], env=env, capture_output=True, text=True)
    check(p.returncode == 3, f"crashed worker -> task exits with the worker's code ({p.returncode})")
    for n in ("s0", "s1"):
        r = pd.read_csv(crash_out / f"{n}.tsv", sep="\t", dtype=str).iloc[0]
        check(r["status"] == "error: worker crashed",
              f"crash fallback row for {n}, even over a stale OK file ({r['status']})")
    check(not (crash_out / "s0_residues.tsv").exists(), "stale residues file gone after a crash")


def main() -> None:
    tests = [test_rg_arithmetic]
    tmp_tests = [test_helix, test_residue_identity_does_not_matter, test_sheet_and_isolation,
                 test_error_contracts, test_calcium_is_not_a_residue, test_warn_contracts,
                 test_builder, test_collector, test_task_script]
    for t in tests:
        t()
    with tempfile.TemporaryDirectory() as d:
        for t in tmp_tests:
            sub = Path(d) / t.__name__
            sub.mkdir()
            try:
                t(sub)
            except Exception as e:  # noqa: BLE001
                FAILURES.append(f"{t.__name__} raised {type(e).__name__}: {e}")
    print("=" * 70)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S) of {CHECKS} checks:")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print(f"all {CHECKS} ssprofile checks passed")


if __name__ == "__main__":
    main()
