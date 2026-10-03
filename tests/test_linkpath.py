"""Tests for the linkpath tool.

The important one is `test_worker_and_collector_columns_agree`: the worker writes a
TSV whose header it owns, and the collector reads that TSV by column name. The two
lists live in different files, so they can drift silently — a column added to the
worker and not to the collector simply never reaches the table, with no error
anywhere. This is the mechanical check for that, and it also asserts the *order*
rule the campaign cares about (each rise constant beside the count it produced).

    .venv/bin/pytest tests/test_linkpath.py -v
"""

import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
TOOL = REPO / "tools" / "linkpath"


def _load(name: str):
    """Import a tools/linkpath module standalone (the package has no __init__)."""
    spec = importlib.util.spec_from_file_location(name, TOOL / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


worker = _load("linkpath_worker")
collector = _load("collect_linkpath")


def test_worker_and_collector_columns_agree():
    """Every metric the worker writes is read back, and nothing else is expected."""
    written = set(worker.METRIC_COLUMNS)
    read = set(collector.RESULT_COLUMNS)
    assert written == read, (
        f"worker-only (written, never collected): {sorted(written - read)}; "
        f"collector-only (expected, never written): {sorted(read - written)}"
    )


def test_collector_type_lists_partition_the_columns():
    """Each collected column is typed exactly once."""
    buckets = (
        collector.FLOAT_COLUMNS,
        collector.INT_COLUMNS,
        collector.BOOL_COLUMNS,
        collector.STR_COLUMNS,
    )
    flat = [c for b in buckets for c in b]
    assert len(flat) == len(set(flat)), "a column is in two type buckets"
    assert set(flat) == set(collector.RESULT_COLUMNS)


def test_worker_header_matches_metric_columns():
    """The TSV the worker actually writes has the header the collector assumes."""
    assert worker.RESULT_COLUMNS == ["name", "status", "path", *worker.METRIC_COLUMNS]


def test_rise_constants_sit_beside_the_counts_they_produced():
    """A reader must see the convention next to the number it generated."""
    order = worker.METRIC_COLUMNS
    assert order[order.index("n_res_min") + 1] == "taut_rise"
    assert order[order.index("n_res_relaxed") + 1] == "relaxed_rise"


@pytest.fixture(scope="module")
def structure(tmp_path_factory):
    """A tiny two-chain structure: two single-residue chains 30 A apart."""
    gemmi = pytest.importorskip("gemmi")
    path = tmp_path_factory.mktemp("linkpath") / "pair.pdb"
    st = gemmi.Structure()
    st.add_model(gemmi.Model("1"))
    for name, x in (("A", 0.137), ("B", 30.137)):
        chain = gemmi.Chain(name)
        res = gemmi.Residue()
        res.name = "GLY"
        res.seqid = gemmi.SeqId(1, " ")
        for aname, el, dx in (("N", "N", 0.0), ("CA", "C", 1.45), ("C", "C", 2.9)):
            atom = gemmi.Atom()
            atom.name = aname
            atom.element = gemmi.Element(el)
            atom.pos = gemmi.Position(x + dx, 0.211, 0.073)
            atom.occ = 1.0
            res.add_atom(atom)
        chain.add_residue(res)
        st[0].add_chain(chain)
    st.setup_entities()
    st.write_pdb(str(path))
    return path


def _run(tmp_path, structure, *extra) -> dict[str, str]:
    task = tmp_path / "task.tsv"
    task.write_text(f"d1\t{structure}\t-\t-\n")
    out = tmp_path / "out"
    subprocess.run(
        [sys.executable, str(TOOL / "linkpath_worker.py"),
         "--task-file", str(task), "--out-dir", str(out), *extra],
        check=True, capture_output=True, text=True,
    )
    rows = list(csv.reader((out / "d1.tsv").open(), delimiter="\t"))
    assert rows[0] == worker.RESULT_COLUMNS
    return dict(zip(rows[0], rows[1]))


def test_rises_are_recorded_with_the_counts(tmp_path, structure):
    row = _run(tmp_path, structure)
    assert row["status"] == "OK"
    assert row["n_res_min"] and row["n_res_relaxed"]
    assert float(row["taut_rise"]) == 3.5
    assert float(row["relaxed_rise"]) == 2.1
    # The counts are re-derivable from the table alone, which is the point.
    import math
    assert int(row["n_res_min"]) == math.ceil(
        float(row["path_dist"]) / float(row["taut_rise"])
    )
    assert int(row["n_res_relaxed"]) == math.ceil(
        float(row["path_dist"]) / float(row["relaxed_rise"])
    )


def test_non_default_rises_are_echoed(tmp_path, structure):
    row = _run(tmp_path, structure, "--taut-rise", "3.8", "--relaxed-rise", "2.4")
    assert float(row["taut_rise"]) == 3.8
    assert float(row["relaxed_rise"]) == 2.4


def test_no_residue_estimate_blanks_the_rises_too(tmp_path, structure):
    """Recording a rise that was not applied to anything is worse than nothing."""
    row = _run(tmp_path, structure, "--no-residue-estimate")
    assert row["status"] == "OK"
    assert row["path_dist"], "distances are still recorded"
    for col in ("n_res_min", "n_res_relaxed", "taut_rise", "relaxed_rise"):
        assert row[col] == "", f"{col} should be NA under --no-residue-estimate"


def test_error_rows_blank_every_metric(tmp_path, structure):
    task = tmp_path / "bad.tsv"
    task.write_text(f"d1\t{tmp_path / 'nope.pdb'}\t-\t-\n")
    out = tmp_path / "badout"
    p = subprocess.run(
        [sys.executable, str(TOOL / "linkpath_worker.py"),
         "--task-file", str(task), "--out-dir", str(out)],
        capture_output=True, text=True,
    )
    assert p.returncode == 0, "a per-design failure must not fail the batch"
    rows = list(csv.reader((out / "d1.tsv").open(), delimiter="\t"))
    row = dict(zip(rows[0], rows[1]))
    assert row["status"].startswith("error:")
    assert all(row[c] == "" for c in worker.METRIC_COLUMNS)
