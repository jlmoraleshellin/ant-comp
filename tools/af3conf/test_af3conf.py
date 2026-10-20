#!/usr/bin/env python3
"""
Self-test for af3conf: synthetic AF3 predictions with a hand-countable answer.

Fixtures have the shape of real AlphaFold3 3.0.x output (``<name>_model.cif``, a
confidences JSON with pae / atom_plddts / atom_chain_ids / token_chain_ids /
token_res_ids, a summary JSON with POSITIONAL chain_ptm / chain_pair_iptm /
chain_pair_pae_min), read off one real AF3 3.0.2 prediction: the tests prove the
arithmetic and the error contracts on that shape, not that every AF3 output has it.

They pin: the ipSAE port is source-identical to boltzconf's and gives the same output;
a hand-computed ipSAE; the chain map comes from token_chain_ids (a summary in another
order is refused); a multi-letter entity scores as plain chains; and every error path
is an error row with NA, never 0.

Optional: test_matches_reference_on_real_af3 re-runs DunbrackLab's ipsae.py and af3conf
on the same REAL AF3 files when AF3CONF_REAL_CIF and IPSAE_REFERENCE are set.

Run:  uv run --no-sync --with pytest pytest tools/af3conf/test_af3conf.py
"""

import importlib
import importlib.util
import inspect
import json
import os
import re
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


worker = _load("af3conf_worker", HERE / "af3conf_worker.py")
run_af3conf = _load("run_af3conf", HERE / "run_af3conf.py")
collect_af3conf = _load("collect_af3conf", HERE / "collect_af3conf.py")
BOLTZCONF_WORKER = HERE.parent / "boltzconf" / "boltzconf_worker.py"
boltzconf_worker = _load("boltzconf_worker", BOLTZCONF_WORKER)

from prosapia.core import CollectCtx, DesignCtx, ManifestCtx  # noqa: E402

PROTEIN_ATOMS = ("N", "CA", "C", "O")
NUCLEOTIDE_ATOMS = ("P", "C1'", "N9")
LIGAND_ATOMS = ("C1", "C2", "C3", "O1", "N1")


# ---------------------------------------------------------------- fixtures


def make_prediction(
    root: Path,
    name: str,
    chains: list[tuple[str, list[str]]],
    pae: np.ndarray | None = None,
    chain_ptm: list[float] | None = None,
    pair_iptm: list[list[float]] | None = None,
    edit_conf=None,
    edit_summary=None,
    write_conf: bool = True,
    write_summary: bool = True,
    cif_name: str | None = None,
    blank_label_seq: str = "",
) -> Path:
    """An AF3-shaped prediction dir; returns the path of <name>_model.cif.

    ``chains`` is [(chain letter, [component per residue])] in FILE order. ALA & co. get
    4 atoms (one token each), DA 3 (one token), LIG a 5-atom HETATM ligand (one token
    PER ATOM), ZN a one-atom ion (one token, label_seq_id '.'). ``blank_label_seq``
    names a chain whose label_seq_id is '?'. chain_pair_pae_min is recomputed from
    ``pae`` as AF3's summary would hold it, unless edit_summary changes it. """
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    atom_lines: list[str] = []
    atom_chain_ids: list[str] = []
    atom_plddts: list[float] = []
    token_chain_ids: list[str] = []
    token_res_ids: list[int] = []
    serial = 0
    for entity, (chain, comps) in enumerate(chains, start=1):
        for resnum, comp in enumerate(comps, start=1):
            if comp == "LIG":
                atoms, group, label_seq = LIGAND_ATOMS, "HETATM", "."
            elif comp == "ZN":
                atoms, group, label_seq = ("ZN",), "HETATM", "."
            elif comp in ("DA", "DC", "DG", "DT"):
                atoms, group, label_seq = NUCLEOTIDE_ATOMS, "ATOM", str(resnum)
            else:
                atoms, group, label_seq = PROTEIN_ATOMS, "ATOM", str(resnum)
            if chain == blank_label_seq:
                label_seq = "?"
            for k, atom in enumerate(atoms):
                serial += 1
                # Distinct, 2-decimal pLDDTs: a one-atom shift always shows.
                plddt = round(40.0 + (serial * 7.31) % 55.0, 2)
                atom_tok = f'"{atom}"' if "'" in atom else atom  # as AF3 quotes C1'
                atom_lines.append(
                    f"{group} {serial} {atom[0]} {atom_tok} . {comp} "
                    f"{chain} {entity} {label_seq} ? {serial * 1.5:.3f} {k:.3f} 0.000 "
                    f"1.00 {plddt:.2f} {resnum} {chain} 1"
                )
                atom_chain_ids.append(chain)
                atom_plddts.append(plddt)
                if comp == "LIG":
                    token_chain_ids.append(chain)
                    token_res_ids.append(resnum)
            if comp != "LIG":
                token_chain_ids.append(chain)
                token_res_ids.append(resnum)
    header = "\n".join(
        [
            f"data_{name}",
            "#",
            "loop_",
            *(
                f"_atom_site.{c}"
                for c in (
                    "group_PDB id type_symbol label_atom_id label_alt_id "
                    "label_comp_id label_asym_id label_entity_id label_seq_id "
                    "pdbx_PDB_ins_code Cartn_x Cartn_y Cartn_z occupancy "
                    "B_iso_or_equiv auth_seq_id auth_asym_id pdbx_PDB_model_num"
                ).split()
            ),
        ]
    )
    cif = d / (cif_name or f"{name}_model.cif")
    cif.write_text(header + "\n" + "\n".join(atom_lines) + "\n#\n")

    n_tok = len(token_chain_ids)
    if pae is None:
        rng = np.random.default_rng(7)
        pae = np.round(rng.uniform(0.5, 30.0, size=(n_tok, n_tok)), 1)
    order = list(dict.fromkeys(token_chain_ids))
    n_ch = len(order)
    blocks = {c: [i for i, x in enumerate(token_chain_ids) if x == c] for c in order}
    pae_min = [
        [round(float(pae[np.ix_(blocks[a], blocks[b])].min()), 2) for b in order]
        for a in order
    ]
    if chain_ptm is None:
        chain_ptm = [round(0.1 + 0.07 * i, 2) for i in range(n_ch)]
    if pair_iptm is None:
        pair_iptm = [
            [chain_ptm[i] if i == j else round(0.2 + 0.05 * (i + j), 2) for j in range(n_ch)]
            for i in range(n_ch)
        ]
    conf = {
        "atom_chain_ids": atom_chain_ids,
        "atom_plddts": atom_plddts,
        "contact_probs": np.eye(n_tok).tolist(),
        "pae": np.asarray(pae).tolist(),
        "token_chain_ids": token_chain_ids,
        "token_res_ids": token_res_ids,
    }
    summary = {
        "chain_iptm": [0.5] * n_ch,
        "chain_pair_iptm": pair_iptm,
        "chain_pair_pae_min": pae_min,
        "chain_ptm": chain_ptm,
        "fraction_disordered": 0.0,
        "has_clash": 0.0,
        "iptm": 0.5,
        "ptm": 0.5,
        "ranking_score": 0.5,
    }
    if edit_conf:
        edit_conf(conf)
    if edit_summary:
        edit_summary(summary)
    prefix = name
    if write_conf:
        (d / f"{prefix}_confidences.json").write_text(json.dumps(conf))
    if write_summary:
        (d / f"{prefix}_summary_confidences.json").write_text(json.dumps(summary))
    return cif


def run_worker(tmp: Path, cif: Path, name: str, binder: str, targets: str,
               cutoff: float = 15.0, sequence: str = "",
               sequence_column: str = "") -> pd.Series:
    """One design through the worker's main(), as af3conf.sh would; its TSV row.
    ``sequence_column`` empty = the flag absent, as the builder writes it."""
    out_dir = tmp / f"out_{name}"
    task = tmp / f"task_{name}.tsv"
    task.write_text(f"{name}\t{cif}\town table\t{sequence}\n")
    sys.argv = [
        "af3conf_worker.py", "--task-file", str(task), "--binder-chains", binder,
        "--target-chains", targets, "--ipsae-pae-cutoff", str(cutoff),
        "--input-column", "alphafold3_path", "--sequence-column", sequence_column,
        "--out-dir", str(out_dir),
    ]
    worker.main()
    return pd.read_csv(out_dir / f"{name}.tsv", sep="\t", dtype=str,
                       keep_default_na=False).iloc[0]


def collect_one(tmp: Path, name: str, table_columns: dict | None = None) -> dict:
    """The same design through the collector, as the driver would call it."""
    out_dir = tmp / f"out_{name}"
    out_dir.mkdir(parents=True, exist_ok=True)
    # The collector derives the leaf from out_dir.name, so give it the real one.
    leaf_dir = tmp / "collect" / "af3conf"
    if leaf_dir.exists():
        for p in leaf_dir.iterdir():
            p.unlink()
    else:
        leaf_dir.mkdir(parents=True)
    for p in out_dir.iterdir():
        (leaf_dir / p.name).write_bytes(p.read_bytes())
    df = pd.DataFrame(index=pd.Index([name]))
    for col, value in (table_columns or {}).items():
        df[col] = value
    ctx = CollectCtx(
        df=df,
        args=Namespace(run_dir=tmp, table="table1", dir_label="", force=False),
        table_name="table1",
        out_dir=leaf_dir,
        status_col="af3conf_status",
        path_col="af3conf_path",
        parent_table=None,
        parent_df=pd.DataFrame(),
        lookup=lambda n, column: None,
        creates_table=False,
        default_input_column="alphafold3_path",
    )
    one = collect_af3conf.collect_af3conf(ctx)
    emitted = list(one(DesignCtx(name, leaf_dir, ctx.lookup)))
    assert len(emitted) == 1
    c = emitted[0]
    return {"data": dict(c.data), "path": c.path, "status": c.status}


def numeric_cells(row: pd.Series) -> dict:
    skip = {"name", "status", "path", "chains", "input_source"}
    return {k: v for k, v in row.items() if k not in skip}


# ---------------------------------------------------------------- the ipSAE port


def test_ipsae_source_identical_to_boltzconf():
    """The copy is the code. Any edit to one side without the other fails here."""
    for fn in ("ipsae_d0", "ipsae_asym", "reject_nucleic_acid"):
        assert inspect.getsource(getattr(worker, fn)) == inspect.getsource(
            getattr(boltzconf_worker, fn)
        ), f"{fn} drifted from boltzconf's copy"
    assert worker.IPSAE_L_MIN == boltzconf_worker.IPSAE_L_MIN
    assert worker.IPSAE_D0_MIN == boltzconf_worker.IPSAE_D0_MIN
    assert worker.NUCLEIC_COMPS == boltzconf_worker.NUCLEIC_COMPS
    assert worker.IPSAE_PER_CHAIN_TEMPLATES == boltzconf_worker.IPSAE_PER_CHAIN_TEMPLATES


def test_ipsae_parity_same_matrix_same_output():
    """The SAME raw PAE matrix through both copies gives identical dicts -- including
    values exactly AT the cutoff (strict '<'), a row with nothing under it, and a
    pair with nothing under it at all (None, not 0.0)."""
    rng = np.random.default_rng(2026)
    pae = np.round(rng.uniform(0.0, 31.0, size=(60, 60)), 1)
    pae[5, 40:50] = 15.0       # exactly at the default cutoff: must NOT count
    pae[50:60, 0:10] = 31.0    # chain pair with nothing under any cutoff below 31
    binder = list(range(0, 10))
    target_a = list(range(10, 40))
    target_b = list(range(40, 60))
    cases = [
        (binder, target_a), (target_a, binder), (binder, target_b),
        (target_b, binder), (list(range(50, 60)), binder),
    ]
    for cutoff in (5.0, 10.0, 15.0, 31.0):
        for align, scored in cases:
            ours = worker.ipsae_asym(pae, align, scored, cutoff)
            theirs = boltzconf_worker.ipsae_asym(pae, align, scored, cutoff)
            assert ours == theirs, (cutoff, align[0], scored[0])
    nothing = worker.ipsae_asym(pae, list(range(50, 60)), list(range(0, 10)), 15.0)
    assert nothing["ipsae"] is None and nothing["n0res"] is None
    assert nothing["nres_align"] == 0 and nothing["nres_scored"] == 0
    n0 = np.arange(0, 400)
    np.testing.assert_array_equal(worker.ipsae_d0(n0), boltzconf_worker.ipsae_d0(n0))


def test_ipsae_d0_clamp():
    # L raised to 26 first: 1.24*(26-15)^(1/3) - 1.8 = 0.9577 -> floored at 1.0.
    assert worker.ipsae_d0(np.array([0, 10, 26]))[0] == 1.0
    assert np.all(worker.ipsae_d0(np.array([0, 10, 26])) == 1.0)
    expected = 1.24 * (100 - 15) ** (1 / 3) - 1.8  # 3.652...
    assert worker.ipsae_d0(np.array([100]))[0] == pytest.approx(expected)
    assert expected == pytest.approx(3.6520, abs=1e-3)


def test_ipsae_hand_computed(tmp_path):
    """Target A (3 residues, tokens 0-2), binder B (2 residues, tokens 3-4).

    binder->A, cutoff 15, d0 = 1.0 (n0 <= 26):
      token 3: pae 1, 3, 20  -> 2 valid -> (1/(1+1) + 1/(1+9)) / 2 = 0.30
      token 4: pae 20,20,20  -> 0 valid
      ipsae_binder_to_a = 0.30, n0res 2, nres_binder_a 1, nres_a 2
    A->binder:
      token 0: pae 0, 0      -> 1.0 (n0 2)
      token 1: pae 2, 30     -> 1/(1+4) = 0.2 (n0 1)
      token 2: pae 30, 30    -> nothing
      ipsae_a_to_binder = 1.0, n0res 2
    max = 1.0 (reverse won) -> d0res_max_a = 1.0; one target chain, so the three
    binder-level aggregates are all 1.0.
    """
    pae = np.full((5, 5), 30.0)
    np.fill_diagonal(pae, 0.8)
    pae[3, 0:3] = [1.0, 3.0, 20.0]
    pae[4, 0:3] = [20.0, 20.0, 20.0]
    pae[0, 3:5] = [0.0, 0.0]
    pae[1, 3:5] = [2.0, 30.0]
    pae[2, 3:5] = [30.0, 30.0]
    cif = make_prediction(tmp_path, "hand", [("A", ["ALA"] * 3), ("B", ["GLY"] * 2)], pae=pae)
    m = worker.measure(cif, "B", ["A"], 15.0, name="hand", source="own table")
    assert m["ipsae_binder_to_a"] == pytest.approx(0.30)
    assert m["ipsae_n0res_binder_to_a"] == 2
    assert m["ipsae_nres_binder_a"] == 1
    assert m["ipsae_nres_a"] == 2
    assert m["ipsae_a_to_binder"] == pytest.approx(1.0)
    assert m["ipsae_n0res_a_to_binder"] == 2
    assert m["ipsae_max_a"] == pytest.approx(1.0)
    assert m["ipsae_d0res_max_a"] == 1.0
    for agg in ("ipsae_binder_mean", "ipsae_binder_max", "ipsae_binder_top2"):
        assert m[agg] == pytest.approx(1.0)
    # The RAW matrix, never the symmetrised one: the two directions differ.
    assert m["ipsae_binder_to_a"] != m["ipsae_a_to_binder"]
    # Trust columns on a clean fixture.
    assert m["n_tokens"] == m["n_res_struct"] == 5
    assert m["binder_n_res"] == 2 and m["n_target_res"] == 3
    assert m["plddt_bfac_maxdiff"] == 0.0
    assert m["chain_pae_min_maxdiff"] <= 0.005
    assert m["stem_matches_name"] == 1


def test_na_is_never_zero_and_aggregates_keep_reference_zero(tmp_path):
    """Binder C reaches A under the cutoff, never B. ipsae_*_b is NA with
    ipsae_nres_binder_b a measured 0; _top2 carries the reference's 0.0 for B, so
    it is HALF _max -- the one-protomer signature. With nothing reached anywhere,
    the aggregates themselves go NA."""
    pae = np.full((6, 6), 25.0)
    pae[4, 0] = 0.0   # C token 4 -> A token 0
    cif = make_prediction(
        tmp_path, "half", [("A", ["ALA"] * 2), ("B", ["ALA"] * 2), ("C", ["GLY"] * 2)],
        pae=pae,
    )
    m = worker.measure(cif, "C", ["A", "B"], 15.0, name="half")
    assert m["ipsae_binder_to_a"] == pytest.approx(1.0)
    assert m["ipsae_binder_to_b"] is None and m["ipsae_max_b"] is None
    assert m["ipsae_n0res_binder_to_b"] is None and m["ipsae_d0res_max_b"] is None
    assert m["ipsae_nres_binder_b"] == 0 and m["ipsae_nres_b"] == 0
    assert m["ipsae_binder_max"] == pytest.approx(1.0)
    assert m["ipsae_binder_top2"] == pytest.approx(0.5)
    assert m["ipsae_binder_mean"] == pytest.approx(0.5)

    row = run_worker(tmp_path, cif, "half", "C", "A,B")
    assert row["ipsae_binder_to_b"] == ""          # NA in the TSV, not 0
    assert row["ipsae_nres_binder_b"] == "0"       # a measured zero
    got = collect_one(tmp_path, "half")
    assert got["status"] == "OK"
    assert got["data"]["ipsae_binder_to_b"] is None
    assert got["data"]["ipsae_nres_binder_b"] == 0.0

    cif2 = make_prediction(
        tmp_path, "none", [("A", ["ALA"] * 2), ("B", ["GLY"] * 2)],
        pae=np.full((4, 4), 25.0),
    )
    m2 = worker.measure(cif2, "B", ["A"], 15.0, name="none")
    for agg in ("ipsae_binder_mean", "ipsae_binder_max", "ipsae_binder_top2"):
        assert m2[agg] is None
    row2 = run_worker(tmp_path, cif2, "none", "B", "A")
    assert row2["status"] == "OK"
    got2 = collect_one(tmp_path, "none")
    assert got2["status"] == "OK" and got2["data"]["ipsae_binder_max"] is None
    # iptm_a is never optional: it comes from the summary, not from the PAE cutoff.
    assert got2["data"]["iptm_a"] is not None


# ---------------------------------------------------------------- the chain map


def test_chain_map_from_token_chain_ids_out_of_order(tmp_path):
    """File order D, B, A -- binder FIRST and the targets reversed. Summary index 0
    is D, 1 is B, 2 is A, read from token_chain_ids; never from the letters'
    alphabetical order or from --target-chains' order. chain_pair_iptm is given
    ASYMMETRIC for D/A to pin the symmetrisation."""
    chain_ptm = [0.11, 0.22, 0.33]            # D, B, A
    pair = [
        [0.11, 0.45, 0.60],                   # D->D, D->B, D->A
        [0.45, 0.22, 0.90],
        [0.70, 0.90, 0.33],                   # A->D = 0.70: (0.60 + 0.70) / 2 = 0.65
    ]
    cif = make_prediction(
        tmp_path, "ooo",
        [("D", ["GLY"] * 3), ("B", ["ALA"] * 4), ("A", ["SER"] * 5)],
        chain_ptm=chain_ptm, pair_iptm=pair,
    )
    m = worker.measure(cif, "D", ["A", "B"], 15.0, name="ooo")
    assert m["chains_ptm_binder"] == pytest.approx(0.11)
    assert m["iptm_a"] == pytest.approx(0.65)
    assert m["iptm_b"] == pytest.approx(0.45)
    assert m["iptm_binder_max"] == pytest.approx(0.65)
    assert m["iptm_binder_mean"] == pytest.approx(0.55)
    assert m["iptm_binder_top2"] == pytest.approx(0.55)
    assert m["binder_n_res"] == 3 and m["n_target_res"] == 9
    assert m["chains"] == "D|A,B"
    # binder pLDDT: the CA atom (index 1 of N, CA, C, O) of each binder residue.
    conf = json.loads((cif.parent / "ooo_confidences.json").read_text())
    ca = [conf["atom_plddts"][4 * r + 1] for r in range(3)]
    assert m["binder_plddt"] == pytest.approx(np.mean(ca))
    assert m["binder_plddt_min"] == pytest.approx(min(ca))


def test_summary_in_a_different_chain_order_is_refused(tmp_path):
    """The summary's positional arrays describe chains in another order than the
    tokens: chain_pair_pae_min no longer matches the pae blocks, and that is an
    ERROR -- never per-chain numbers about the wrong chain."""
    n = 9
    pae = np.full((n, n), 28.0)
    pae[0:3, 3:6] = 4.0      # A->B min 4
    pae[3:6, 0:3] = 9.0      # B->A min 9
    pae[0:3, 6:9] = 17.0     # A->C min 17
    np.fill_diagonal(pae, 0.8)

    def permute(summary):
        p = [1, 0, 2]   # swap the first two chains
        for k in ("chain_ptm", "chain_iptm"):
            summary[k] = [summary[k][i] for i in p]
        for k in ("chain_pair_iptm", "chain_pair_pae_min"):
            summary[k] = [[summary[k][i][j] for j in p] for i in p]

    cif = make_prediction(
        tmp_path, "perm", [("A", ["ALA"] * 3), ("B", ["ALA"] * 3), ("C", ["GLY"] * 3)],
        pae=pae, edit_summary=permute,
    )
    with pytest.raises(ValueError, match="chain_pair_pae_min"):
        worker.measure(cif, "C", ["A", "B"], 15.0)
    row = run_worker(tmp_path, cif, "perm", "C", "A,B")
    assert row["status"].startswith("error:") and "positional chain" in row["status"]
    assert all(v == "" for v in numeric_cells(row).values())


def test_multi_letter_entity_is_plain_chains(tmp_path):
    """AF3 splits an ABC entity into chains A, B, C in token_chain_ids and in the
    CIF. 'A:C' on the submit side must expand to exactly those letters, and the
    per-chain columns come out a, b, c."""
    cif = make_prediction(
        tmp_path, "trimer",
        [("A", ["ALA", "LYS"]), ("B", ["ALA", "LYS"]), ("C", ["ALA", "LYS"]),
         ("D", ["GLY"] * 3)],
    )
    rows, sub, args = build_manifest(
        tmp_path, {"trimer": cif}, binder_chains="D", target_chains="A:C"
    )
    assert rows[0][2] == "A,B,C"
    row = run_worker(tmp_path, cif, "trimer", "D", rows[0][2])
    assert row["status"] == "OK"
    for k in ("a", "b", "c"):
        assert row[f"iptm_{k}"] != "" and f"ipsae_max_{k}" in row.index
    assert row["n_target_res"] == "6" and row["binder_n_res"] == "3"


# ---------------------------------------------------------------- error paths


def _error(tmp_path, cif, name, binder="B", targets="A") -> str:
    row = run_worker(tmp_path, cif, name, binder, targets)
    assert row["status"].startswith("error:"), row["status"]
    assert all(v == "" for v in numeric_cells(row).values()), "an error row is NA, never a number"
    assert row["chains"] == f"{binder}|{targets}"
    got = collect_one(tmp_path, name)
    assert got["status"] == row["status"]
    assert all(
        v is None for k, v in got["data"].items() if k not in ("chains", "input_source")
    )
    return row["status"]


TWO_CHAINS = [("A", ["ALA"] * 3), ("B", ["GLY"] * 2)]


def test_missing_confidences_sibling(tmp_path):
    cif = make_prediction(tmp_path, "noconf", TWO_CHAINS, write_conf=False)
    assert "noconf_confidences.json is missing" in _error(tmp_path, cif, "noconf")


def test_missing_summary_sibling(tmp_path):
    cif = make_prediction(tmp_path, "nosum", TWO_CHAINS, write_summary=False)
    assert "nosum_summary_confidences.json is missing" in _error(tmp_path, cif, "nosum")


def test_sibling_from_file_prefix_not_design_name(tmp_path):
    """The table row is called 'renamed'; the files are 'af3name_*'. The siblings are
    found from the file, and stem_matches_name says the name did not match."""
    cif = make_prediction(tmp_path, "af3name", TWO_CHAINS)
    row = run_worker(tmp_path, cif, "renamed", "B", "A")
    assert row["status"] == "OK" and row["stem_matches_name"] == "0"


def test_not_a_model_cif(tmp_path):
    cif = make_prediction(tmp_path, "odd", TWO_CHAINS, cif_name="odd.cif")
    assert "not named <prefix>_model.cif" in _error(tmp_path, cif, "odd")


def test_ligand_tokens_do_not_match_residues(tmp_path):
    """A ligand in an UNNAMED chain is still five tokens for one residue: the token
    indices stop being residue indices, so it is an error, not a skipped chain."""
    cif = make_prediction(tmp_path, "lig", [*TWO_CHAINS, ("C", ["LIG"])])
    assert "tokens but" in _error(tmp_path, cif, "lig")


def test_renumbered_tokens_are_refused(tmp_path):
    cif = make_prediction(
        tmp_path, "renum", TWO_CHAINS,
        edit_conf=lambda c: c.update(token_res_ids=[r + 1 for r in c["token_res_ids"]]),
    )
    assert "numbered differently" in _error(tmp_path, cif, "renum")


def test_token_chain_order_must_match_structure(tmp_path):
    def swap(c):
        c["token_chain_ids"] = ["B", "B", "A", "A", "A"]
    cif = make_prediction(tmp_path, "tswap", TWO_CHAINS, edit_conf=swap)
    assert "do not list the same residues" in _error(tmp_path, cif, "tswap")


def test_nucleic_acid_among_named_chains(tmp_path):
    """A standard nucleotide is ONE token, so every token check passes -- only the
    component name catches it."""
    cif = make_prediction(tmp_path, "dna", [("A", ["DA", "DC", "DG"]), ("B", ["GLY"] * 2)])
    assert "look like a nucleic acid" in _error(tmp_path, cif, "dna")


def test_nucleic_acid_in_unnamed_chain_is_fine(tmp_path):
    cif = make_prediction(tmp_path, "dna_off", [*TWO_CHAINS, ("C", ["DA", "DT"])])
    row = run_worker(tmp_path, cif, "dna_off", "B", "A")
    assert row["status"] == "OK"


def test_wrong_chain_count_in_chain_pair_iptm(tmp_path):
    def drop(s):
        s["chain_pair_iptm"] = [r[:1] for r in s["chain_pair_iptm"][:1]]
    cif = make_prediction(tmp_path, "cnt", TWO_CHAINS, edit_summary=drop)
    assert "'chain_pair_iptm' is 1 x 1" in _error(tmp_path, cif, "cnt")


def test_wrong_chain_count_in_chain_ptm(tmp_path):
    cif = make_prediction(
        tmp_path, "cnt2", TWO_CHAINS,
        edit_summary=lambda s: s.update(chain_ptm=s["chain_ptm"] + [0.9]),
    )
    assert "'chain_ptm' has shape (3,)" in _error(tmp_path, cif, "cnt2")


def test_nan_in_pae(tmp_path):
    # json.dumps writes NaN as the bare token NaN, which json.load reads back.
    cif = make_prediction(
        tmp_path, "nan", TWO_CHAINS,
        edit_conf=lambda c: c["pae"][3].__setitem__(0, float("nan")),
    )
    assert "non-finite" in _error(tmp_path, cif, "nan")


def test_null_in_chain_ptm(tmp_path):
    cif = make_prediction(
        tmp_path, "null", TWO_CHAINS,
        edit_summary=lambda s: s["chain_ptm"].__setitem__(1, None),
    )
    assert "non-finite or null" in _error(tmp_path, cif, "null")


def test_plddt_disagrees_with_bfactor(tmp_path):
    def shift(c):
        c["atom_plddts"] = c["atom_plddts"][1:] + c["atom_plddts"][:1]
    cif = make_prediction(tmp_path, "shift", TWO_CHAINS, edit_conf=shift)
    assert "B-factors disagree" in _error(tmp_path, cif, "shift")


def test_absent_chain(tmp_path):
    cif = make_prediction(tmp_path, "absent", TWO_CHAINS)
    assert "absent from" in _error(tmp_path, cif, "absent", binder="Z", targets="A")


def test_unresolvable_path_is_error_row(tmp_path):
    out_dir = tmp_path / "out_blank"
    task = tmp_path / "task_blank.tsv"
    task.write_text("blank\t\town table\n")
    sys.argv = ["w", "--task-file", str(task), "--binder-chains", "B",
                "--target-chains", "A", "--out-dir", str(out_dir)]
    worker.main()
    row = pd.read_csv(out_dir / "blank.tsv", sep="\t", dtype=str, keep_default_na=False).iloc[0]
    assert row["status"].startswith("error: could not resolve")


def test_named_one_atom_ion_is_refused(tmp_path):
    """One atom, one token, one residue: every count check passes. Only the
    named-chain label_seq_id rule catches it."""
    cif = make_prediction(tmp_path, "ion", [("A", ["ALA"] * 3), ("B", ["GLY"] * 2), ("C", ["ZN"])])
    status = _error(tmp_path, cif, "ion", binder="B", targets="A,C")
    assert "NAMED chain C" in status and "no integer label_seq_id" in status


def test_unnamed_one_atom_ion_is_fine(tmp_path):
    cif = make_prediction(tmp_path, "ion_off", [*TWO_CHAINS, ("C", ["ZN"])])
    row = run_worker(tmp_path, cif, "ion_off", "B", "A")
    assert row["status"] == "OK" and row["n_res_struct"] == "6"


def test_named_chain_with_unusable_label_seq_is_refused(tmp_path):
    cif = make_prediction(tmp_path, "noseqid", TWO_CHAINS, blank_label_seq="A")
    assert "no integer label_seq_id" in _error(tmp_path, cif, "noseqid")


# ---------------------------------------------------------------- new trust columns


def test_binder_segment_matches_prosapia_select_chains():
    """The worker's copy of the letter rule agrees with prosapia's own."""
    from prosapia.utils import select_chains
    for seq in ("AAA/BBB/CCC/DDDD", "MKV", "AA/G", "X/Y/Z/W/V"):
        for letter, segment in select_chains(seq, None).items():
            assert worker.binder_segment(seq, letter, "col") == segment
    with pytest.raises(ValueError, match="no segment for binder chain D"):
        worker.binder_segment("AAA/BBB", "D", "col")
    with pytest.raises(ValueError, match="EMPTY segment"):
        worker.binder_segment("AAA//CCC", "B", "col")


COMPLEX = [("A", ["ALA", "LYS", "GLU"]), ("B", ["GLY", "TRP"])]


def test_binder_seq_match_one_zero_and_na(tmp_path):
    cif = make_prediction(tmp_path, "seq", COMPLEX)
    ok = run_worker(tmp_path, cif, "seq", "B", "A", sequence="AKE/GW",
                    sequence_column="mkcomplex_sequence")
    assert ok["status"] == "OK" and ok["binder_seq_match"] == "1"
    # A mismatch is a TRUST VALUE, not an error: the row is still scored.
    bad = run_worker(tmp_path, cif, "seq", "B", "A", sequence="AKE/GY",
                     sequence_column="mkcomplex_sequence")
    assert bad["status"] == "OK" and bad["binder_seq_match"] == "0"
    got = collect_one(tmp_path, "seq")
    assert got["status"] == "OK" and got["data"]["binder_seq_match"] == 0
    # Flag absent: NA (unchecked), and an OK row without it is still OK.
    off = run_worker(tmp_path, cif, "seq", "B", "A")
    assert off["status"] == "OK" and off["binder_seq_match"] == ""
    got = collect_one(tmp_path, "seq")
    assert got["status"] == "OK" and got["data"]["binder_seq_match"] is None


def test_binder_seq_match_cannot_tell_which_chain_is_the_binder(tmp_path):
    """THE LIMIT. Both sides are keyed by the --binder-chains letter, so a TARGET
    letter passed as the binder still matches its own segment: binder_seq_match = 1
    on a row whose binder_* / iptm_* describe a target protomer. Only binder_n_res
    against the DESIGNED binder length (2 here) exposes it."""
    target = ["ALA", "LYS", "GLU", "SER"]
    cif = make_prediction(
        tmp_path, "wrongbinder",
        [("A", target), ("B", target), ("C", target), ("D", ["GLY", "TRP"])],
    )
    row = run_worker(tmp_path, cif, "wrongbinder", "A", "B,C,D", sequence="AKES/AKES/AKES/GW",
                     sequence_column="mkcomplex_sequence")
    assert row["status"] == "OK" and row["binder_seq_match"] == "1"
    assert row["binder_n_res"] == "4"   # a target protomer's length, not the binder's 2


def test_binder_seq_unresolvable_or_missing_segment_is_error(tmp_path):
    cif = make_prediction(tmp_path, "seqerr", COMPLEX)
    row = run_worker(tmp_path, cif, "seqerr", "B", "A", sequence="",
                     sequence_column="mkcomplex_sequence")
    assert "could not be resolved" in row["status"]
    row = run_worker(tmp_path, cif, "seqerr", "B", "A", sequence="AKE",
                     sequence_column="mkcomplex_sequence")
    assert "no segment for binder chain B" in row["status"]
    assert row["binder_seq_match"] == ""


def test_binder_seq_match_survives_a_later_error(tmp_path):
    """Measured before the error, so it stays on the error row -- '0 and the
    summary is missing' says the file is not this design's."""
    cif = make_prediction(tmp_path, "seqlate", COMPLEX, write_summary=False)
    row = run_worker(tmp_path, cif, "seqlate", "B", "A", sequence="AKE/GY",
                     sequence_column="mkcomplex_sequence")
    assert row["status"].startswith("error:") and "summary_confidences" in row["status"]
    assert row["binder_seq_match"] == "0"
    assert row["binder_plddt"] == "" and row["ranking_score"] == ""
    got = collect_one(tmp_path, "seqlate")
    assert got["data"]["binder_seq_match"] == 0
    assert all(v is None for k, v in got["data"].items()
               if k not in ("chains", "input_source", "binder_seq_match"))


def test_ranking_score_echoed_and_required(tmp_path):
    cif = make_prediction(tmp_path, "rank", TWO_CHAINS,
                          edit_summary=lambda s: s.update(ranking_score=0.83))
    row = run_worker(tmp_path, cif, "rank", "B", "A")
    assert float(row["ranking_score"]) == 0.83
    cif = make_prediction(tmp_path, "norank", TWO_CHAINS,
                          edit_summary=lambda s: s.pop("ranking_score"))
    assert "has no 'ranking_score'" in _error(tmp_path, cif, "norank")
    cif = make_prediction(tmp_path, "nanrank", TWO_CHAINS,
                          edit_summary=lambda s: s.update(ranking_score=None))
    assert "non-finite or null" in _error(tmp_path, cif, "nanrank")


def _block_pae(n_per: int, minima: np.ndarray) -> np.ndarray:
    """A pae whose chain-block minima are exactly ``minima`` (n_chains x n_chains)."""
    k = minima.shape[0]
    pae = np.full((k * n_per, k * n_per), 29.0)
    for a in range(k):
        for b in range(k):
            pae[a * n_per, b * n_per] = minima[a, b]
    return pae


def test_chain_pae_min_margin_separated_row(tmp_path):
    """Distinct block minima: the nearest swap is far from the identity order."""
    minima = np.array([[0.8, 4.0, 17.0], [9.0, 0.8, 12.0], [21.0, 6.0, 0.8]])
    cif = make_prediction(
        tmp_path, "sep", [("A", ["ALA"] * 2), ("B", ["ALA"] * 2), ("C", ["GLY"] * 2)],
        pae=_block_pae(2, minima),
    )
    m = worker.measure(cif, "C", ["A", "B"], 15.0)
    # The three swaps, permuting the KNOWN minima directly (independent of the tool).
    best = min(
        np.abs(minima[np.ix_(p, p)] - minima).max()
        for p in ([1, 0, 2], [2, 1, 0], [0, 2, 1])
    )
    assert m["chain_pae_min_maxdiff"] == pytest.approx(0.0)
    assert m["chain_pae_min_margin"] == pytest.approx(best)
    assert m["chain_pae_min_margin"] > worker.CHAIN_PAE_MIN_TOL


def test_chain_pae_min_margin_blind_row_is_ok_not_error(tmp_path):
    """Two target copies with identical block minima to everything: swapping them
    changes nothing, so margin is 0 -- the order check could not have caught that
    swap. The row is still OK: such rows may be right."""
    minima = np.array([[0.8, 5.0, 3.0], [5.0, 0.8, 3.0], [3.0, 3.0, 0.8]])
    cif = make_prediction(
        tmp_path, "blind", [("A", ["ALA"] * 2), ("B", ["ALA"] * 2), ("C", ["GLY"] * 2)],
        pae=_block_pae(2, minima),
    )
    row = run_worker(tmp_path, cif, "blind", "C", "A,B")
    assert row["status"] == "OK"
    assert float(row["chain_pae_min_margin"]) <= worker.CHAIN_PAE_MIN_TOL
    got = collect_one(tmp_path, "blind")
    assert got["data"]["chain_pae_min_margin"] == pytest.approx(0.0)


def test_confidences_matrices_dropped_after_load(tmp_path):
    cif = make_prediction(tmp_path, "mem", TWO_CHAINS)
    conf, pae = worker.load_confidences(cif.parent / "mem_confidences.json")
    assert "contact_probs" not in conf and "pae" not in conf
    assert pae.shape == (5, 5)


# ---------------------------------------------------------------- collect


def test_collect_missing_blanks_stale_per_chain_columns(tmp_path):
    got = collect_one(tmp_path, "gone", table_columns={
        "af3conf_iptm_z": 0.9, "af3conf_ipsae_max_z": 0.8, "af3conf_iptm_binder_max": 0.7,
    })
    assert got["status"] == "missing"
    assert got["data"]["iptm_z"] is None and got["data"]["ipsae_max_z"] is None
    assert all(v is None for v in got["data"].values())


def test_collect_ok_row_with_a_hole_is_an_error(tmp_path):
    cif = make_prediction(tmp_path, "hole", TWO_CHAINS)
    run_worker(tmp_path, cif, "hole", "B", "A")
    tsv = tmp_path / "out_hole" / "hole.tsv"
    df = pd.read_csv(tsv, sep="\t", dtype=str, keep_default_na=False)
    df.loc[0, "iptm_a"] = ""
    df.to_csv(tsv, sep="\t", index=False)
    got = collect_one(tmp_path, "hole")
    assert got["status"] == "error: OK row missing iptm_a"


def test_collect_types_on_ok_row(tmp_path):
    cif = make_prediction(tmp_path, "types", TWO_CHAINS)
    run_worker(tmp_path, cif, "types", "B", "A")
    got = collect_one(tmp_path, "types")
    d = got["data"]
    assert got["status"] == "OK" and got["path"] == str(cif)
    assert isinstance(d["n_tokens"], int) and isinstance(d["binder_plddt"], float)
    assert d["chains"] == "B|A" and d["input_source"] == "own table"


# ---------------------------------------------------------------- column names


def test_shared_column_names_identical_to_boltzconf():
    """Every af3conf column except its own new trust column exists in boltzconf
    under the same bare name, so the two leaves sit side by side."""
    ours = worker.metric_columns(["A", "B"])
    theirs = set(boltzconf_worker.metric_columns(["A", "B"]))
    assert [c for c in ours if c not in theirs] == [
        "chain_pae_min_maxdiff", "chain_pae_min_margin", "binder_seq_match",
        "ranking_score",
    ]
    assert set(collect_af3conf.RESULT_COLUMNS) == set(worker.TRUST_COLUMNS + worker.SCIENCE_COLUMNS)


# ---------------------------------------------------------------- manifest


def build_manifest(tmp: Path, structures: dict[str, Path], **overrides):
    out_dir = tmp / overrides.pop("leaf", "af3conf")
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        {"alphafold3_path": [str(p) for p in structures.values()]},
        index=pd.Index(list(structures), name="name"),
    )
    for key, value in overrides.pop("columns", {}).items():
        df[key] = value
    if overrides.pop("drop_column", False):
        df = df.drop(columns=["alphafold3_path"])
    args = Namespace(
        run_dir=tmp,
        table="table1",
        input_column=overrides.pop("input_column", "alphafold3_path"),
        force=False,
        gpus_per_task=1,
        binder_chains=overrides.pop("binder_chains", "B"),
        target_chains=overrides.pop("target_chains", "A"),
        ipsae_pae_cutoff=overrides.pop("ipsae_pae_cutoff", 15.0),
        designs_per_task=overrides.pop("designs_per_task", 20),
        allow_lineage=overrides.pop("allow_lineage", False),
        sequence_column=overrides.pop("sequence_column", None),
    )
    assert not overrides, overrides
    ctx = ManifestCtx(
        df=df, args=args, out_dir=out_dir,
        lookup=lambda name, column: None, write_meta=lambda **_: None,
    )
    rows = run_af3conf.build_af3conf_manifest(ctx)
    sub = []
    for row in rows:
        sub.extend(line.split("\t") for line in Path(row[0]).read_text().splitlines())
    return rows, sub, args


def test_manifest_shape_and_cpu_only(tmp_path):
    rows, sub, args = build_manifest(
        tmp_path, {"d1": tmp_path / "d1_model.cif"}, target_chains="A:C", binder_chains="D",
        ipsae_pae_cutoff=10.0,
    )
    assert args.gpus_per_task == 0
    assert rows[0][1:] == ("D", "A,B,C", "10.0", "alphafold3_path", "")
    assert sub == [["d1", os.path.abspath(tmp_path / "d1_model.cif"), "own table", ""]]


def test_manifest_sequence_column(tmp_path):
    """The sequence rides as the 4th sub-manifest field, normalised as epitope does
    (spaces dropped, upper-cased); the column name is manifest field 6. A blank cell
    stays blank (the worker errors on it); an absent column refuses at submit."""
    s = {"d1": tmp_path / "d1_model.cif", "d2": tmp_path / "d2_model.cif"}
    rows, sub, _ = build_manifest(
        tmp_path, s, binder_chains="D", target_chains="A:C",
        sequence_column="mkcomplex_sequence",
        columns={"mkcomplex_sequence": ["aaa/AAA/AAA/ g g g", ""]},
    )
    assert rows[0][5] == "mkcomplex_sequence"
    assert sub[0][3] == "AAA/AAA/AAA/GGG" and sub[1][3] == ""
    with pytest.raises(ValueError, match="no --sequence-column"):
        build_manifest(tmp_path, s, sequence_column="mkcomplex_sequence")


def test_manifest_refuses_labels_that_collide_with_per_chain_columns(tmp_path):
    """Any '_'-separated label segment equal to or starting with iptm/ipsae/ipae is
    refused -- not only a leading one: '-l x_iptm' collides with '-l x'."""
    s = {"d1": tmp_path / "d1_model.cif"}
    for bad in ("af3conf_iptm", "af3conf_ipsae2", "af3conf_IPAE_x", "af3conf_x_iptm",
                "af3conf_run1_IPSAE_b"):
        with pytest.raises(ValueError, match="per-chain column stem"):
            build_manifest(tmp_path, s, leaf=bad)
    for good in ("af3conf_tmpl", "af3conf_x_tip", "af3conf_noiptm"):
        rows, _, _ = build_manifest(tmp_path, s, leaf=good)
        assert rows


def _stale_for(leaf: str, table_columns: list[str]) -> list[str]:
    ctx = Namespace(out_dir=Path(leaf), df=pd.DataFrame(columns=table_columns))
    return collect_af3conf._stale_per_chain(ctx)


def test_prefix_collision_that_the_label_refusal_prevents():
    """WHY the refusal exists. Collecting leaf af3conf_x reads a sibling
    '-l x_iptm' leaf's columns as its OWN per-chain columns -- status, path and
    metrics included -- and would blank them. The unlabelled leaf does the same to
    '-l iptm'. A label without such a segment is never touched."""
    sibling = ["af3conf_x_iptm_status", "af3conf_x_iptm_path", "af3conf_x_iptm_ipsae_max_a"]
    assert _stale_for("af3conf_x", sibling) == ["iptm_status", "iptm_path", "iptm_ipsae_max_a"]
    assert _stale_for("af3conf", ["af3conf_iptm_status"]) == ["iptm_status"]
    assert _stale_for("af3conf_x", ["af3conf_x_tip_status", "af3conf_x_tip_iptm_a"]) == []


def test_manifest_refusals(tmp_path):
    s = {"d1": tmp_path / "d1_model.cif"}
    with pytest.raises(ValueError, match="no default input column"):
        build_manifest(tmp_path, s, input_column=run_af3conf.NO_DEFAULT_COLUMN)
    with pytest.raises(ValueError, match="exactly one chain"):
        build_manifest(tmp_path, s, binder_chains="B,C")
    with pytest.raises(ValueError, match="both binder and target"):
        build_manifest(tmp_path, s, binder_chains="A", target_chains="A,B")
    with pytest.raises(ValueError, match="does NOT walk the lineage"):
        build_manifest(tmp_path, s, drop_column=True)
    with pytest.raises(ValueError, match="must be > 0"):
        build_manifest(tmp_path, s, ipsae_pae_cutoff=0.0)


def test_manifest_blank_cell_never_borrows(tmp_path):
    rows, sub, _ = build_manifest(
        tmp_path, {"d1": tmp_path / "d1_model.cif", "d2": tmp_path / "d2_model.cif"},
        columns={"alphafold3_path": [str(tmp_path / "d1_model.cif"), ""]},
    )
    assert sub[1] == ["d2", "", "own table", ""]


# ---------------------------------------------------------------- wiring


def test_spec_wiring():
    sys.path.insert(0, str(HERE.parent))
    try:
        spec = importlib.import_module(f"{HERE.name}.spec")
    finally:
        sys.path.pop(0)
    tool = spec.TOOL
    assert tool.name == "af3conf" and tool.action == "update"
    assert tool.default_input_column == run_af3conf.NO_DEFAULT_COLUMN
    assert Path(tool.default_script).name == "af3conf.sh"


def test_task_script_parses_and_avoids_bash_specials():
    sh = HERE / "af3conf.sh"
    subprocess.run(["bash", "-n", str(sh)], check=True)
    assigned = set(re.findall(r"^([A-Z_]+)=", sh.read_text(), flags=re.M))
    specials = {"GROUPS", "UID", "EUID", "PPID", "PIPESTATUS", "SECONDS", "RANDOM",
                "LINENO", "HOSTNAME", "IFS", "REPLY", "PWD", "PATH", "HOME"}
    assert not assigned & specials


# ---------------------------------------------------------------- real data (opt-in)


@pytest.mark.skipif(
    not (os.environ.get("AF3CONF_REAL_CIF") and os.environ.get("IPSAE_REFERENCE")),
    reason="set AF3CONF_REAL_CIF=<AF3 <name>_model.cif> and IPSAE_REFERENCE=<ipsae.py>",
)
def test_matches_reference_on_real_af3(tmp_path):
    """DunbrackLab's ipsae.py and af3conf on the SAME real AF3 files must agree on
    every ipSAE number, the ipTM it reads from chain_pair_iptm, and the CA pLDDT.
    Optional env: AF3CONF_BINDER (default: last chain), AF3CONF_TARGETS (default:
    every other chain), AF3CONF_CUTOFF (default 15)."""
    import shutil

    src = Path(os.environ["AF3CONF_REAL_CIF"])
    prefix = worker.prediction_prefix(src)
    for suffix in ("_model.cif", "_confidences.json", "_summary_confidences.json"):
        shutil.copy(src.parent / f"{prefix}{suffix}", tmp_path / f"{prefix}{suffix}")
    cif = tmp_path / src.name
    cutoff = float(os.environ.get("AF3CONF_CUTOFF", "15"))
    order = worker.read_residues(cif).chain_order
    binder = os.environ.get("AF3CONF_BINDER", order[-1])
    targets = os.environ.get("AF3CONF_TARGETS", ",".join(c for c in order if c != binder))
    targets = [c for c in targets.split(",") if c]
    m = worker.measure(cif, binder, targets, cutoff, name=prefix)
    subprocess.run(
        [sys.executable, "-I", os.environ["IPSAE_REFERENCE"],
         str(tmp_path / f"{prefix}_confidences.json"), str(cif), f"{cutoff:g}", "10"],
        cwd=tmp_path, check=True, capture_output=True,
    )
    out = next(tmp_path.glob(f"*_{cutoff:g}_10.txt"))
    ref = {}
    for line in out.read_text().splitlines():
        f = line.split()
        if len(f) > 20 and f[4] in ("asym", "max"):
            ref[(f[0], f[1], f[4])] = f
    for t in targets:
        k = worker.chain_key(t)
        for (c1, c2), col in (((binder, t), f"ipsae_binder_to_{k}"),
                              ((t, binder), f"ipsae_{k}_to_binder")):
            # Reference columns: 5 ipSAE, 8 ipTM_af, 13 n0res, 19 nres1, 20 nres2.
            r = ref[(c1, c2, "asym")]
            if int(r[13]) == 0:
                assert m[col] is None, col
            else:
                assert m[col] == pytest.approx(float(r[5]), abs=1e-6), col
        fwd, rev = ref[(binder, t, "asym")], ref[(t, binder, "asym")]
        # The reference reads chain_pair_iptm per direction; af3conf symmetrises it.
        assert m[f"iptm_{k}"] == pytest.approx(
            (float(fwd[8]) + float(rev[8])) / 2, abs=1e-6
        )
        assert m[f"ipsae_nres_binder_{k}"] == int(fwd[19])
        assert m[f"ipsae_nres_{k}"] == int(fwd[20])
    assert m["plddt_bfac_maxdiff"] <= worker.PLDDT_BFAC_TOL
    assert m["chain_pae_min_maxdiff"] <= worker.CHAIN_PAE_MIN_TOL
    summary = json.loads((tmp_path / f"{prefix}_summary_confidences.json").read_text())
    assert m["ranking_score"] == summary["ranking_score"]
    # Reported, not asserted: on a confident complex it may legitimately be small.
    print(f"chain_pae_min_margin={m['chain_pae_min_margin']:.3f} "
          f"maxdiff={m['chain_pae_min_maxdiff']:.3f}")
