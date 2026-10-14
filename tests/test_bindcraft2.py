"""Tests for the bindcraft2 tool.

Three layers, deliberately:

1. Pure unit tests of the manifest builder and the collector.
2. A cross-check of every settings file this tool generates against **BindCraft2's
   own validator** (``reject_unrecognized_settings`` / ``load_settings`` /
   ``resolve_prepared_states``), run in the upstream checkout's own venv. Asserting
   that our JSON "looks right" proves nothing; asserting that upstream accepts it,
   and resolves the targets we meant, is the check that would actually have failed.
3. An end-to-end run of ``bindcraft2.sh`` against a fake ``bindcraft`` binary that
   dumps its argv, which is the only way to prove the quoting fix really survives
   the manifest -> args file -> bash -> argv path.

Layer 2 and 3 skip themselves when their prerequisites are absent, so the suite
still runs on a machine without the upstream checkout.

    .venv/bin/pytest tests/test_bindcraft2.py -v
"""

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
TOOL_DIR = REPO / "tools" / "bindcraft2"

# The upstream checkout and the CPU venv built from it. Layer-2 tests need both.
UPSTREAM = Path(
    os.environ.get("BINDCRAFT2_SOURCE", REPO.parent / "BindCraft2")
).resolve()
UPSTREAM_PYTHON = UPSTREAM / ".venv-cpu" / "bin" / "python"
needs_upstream = pytest.mark.skipif(
    not UPSTREAM_PYTHON.is_file(),
    reason=(
        f"no BindCraft2 venv at {UPSTREAM_PYTHON}; clone PacesaLab/BindCraft2 at "
        f"v1.0.3 and `uv venv --python 3.12 .venv-cpu && uv pip install -e .`"
    ),
)


def _load(module_name: str):
    """Load one of the tool's modules straight from its file.

    The tool folder is not a package on disk (prosapia loads it as a synthetic one),
    so the tests import by path rather than by name.
    """
    spec = importlib.util.spec_from_file_location(
        f"_bc2_test_{module_name}", TOOL_DIR / f"{module_name}.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


run_bc2 = _load("run_bindcraft2")
collect_bc2 = _load("collect_bindcraft2")


class Args:
    """A stand-in for the parsed argparse namespace the builder reads."""

    def __init__(self, **overrides):
        defaults = dict(
            fetch_weights_only=False,
            trajectory_only=False,
            reuse_campaigns=None,
            target_pdb=None,
            shipped_target=[],
            target=[],
            extra_target=[],
            targets_file=None,
            chains=None,
            hotspots=None,
            coldspots=None,
            target_weight=None,
            target_objective=None,
            binder_lengths=None,
            num_designs=None,
            max_trajectories=None,
            modality=None,
            property=[],
            core=None,
            campaign_seed=None,
            metadata=None,
            design_workers=None,
            workers_per_gpu=None,
            save_monomers=False,
            no_resume=False,
            extra_settings=None,
            set=[],
            table=None,
        )
        defaults.update(overrides)
        for key, value in defaults.items():
            setattr(self, key, value)


def no_lookup(name, column):  # pragma: no cover - only hit by a stray {expr}
    raise AssertionError(f"unexpected lineage lookup for {column!r}")


def build(args, name="camp", target_path=None, extra=None, campaign_dir=None):
    return run_bc2._build_settings(
        name,
        target_path,
        args,
        no_lookup,
        campaign_dir or Path("/runs/out/campaigns/camp"),
        extra or {},
    )


# --------------------------------------------------------------------------
# 1. target specs
# --------------------------------------------------------------------------


def test_target_spec_parses_fields_and_aliases(tmp_path):
    entry = run_bc2.parse_target_spec(
        "name=hPDL1;path=t.pdb;chains=A,B;hotspots=A54,A56,A66-70;weight=-0.5",
        "--target #1",
    )
    assert entry == {
        "name": "hPDL1",
        "target_path": "t.pdb",
        "chains": "A,B",
        # the comma-separated hotspot list must survive: ';' is the field separator
        # precisely so that it can.
        "hotspots": "A54,A56,A66-70",
        "weight": "-0.5",
    }


def test_target_spec_defaults_name_to_file_stem():
    entry = run_bc2.parse_target_spec("path=targets/hPD1.pdb", "--target #1")
    assert entry["name"] == "hPD1"


@pytest.mark.parametrize(
    "spec, message",
    [
        ("name=X;wieght=1", "unknown target field"),
        ("name=X;name=Y", "given twice"),
        ("chains=A", "needs a `name`"),
        ("name=X;justawordhere", "not `key=value`"),
        ("", "empty target spec"),
    ],
)
def test_target_spec_refuses_bad_input(spec, message):
    with pytest.raises(run_bc2.SettingsConfigError, match=message):
        run_bc2.parse_target_spec(spec, "--target #1")


def test_reserved_target_names_are_refused():
    """A target called `mean` would land in the same column as the `__mean` summary."""
    for name in ("mean", "worst", "best", "spread", "selectivity"):
        with pytest.raises(run_bc2.SettingsConfigError, match="reserved"):
            run_bc2.parse_target_spec(f"name={name};path=x.pdb", "--target #1")


def test_shipped_and_described_targets_are_refused():
    """Upstream accepts this combination and silently drops the shipped targets.

    Verified against v1.0.3: campaign_over_presets({'target': ['hPDL1','mPDL1'],
    'targets': [...]}) returns ONLY the explicit list. A campaign that quietly
    designs against one target when it was told three is worse than a refused submit.
    """
    args = Args(shipped_target=["hPDL1"], target=["name=mine;path=x.pdb"])
    with pytest.raises(run_bc2.SettingsConfigError, match="REPLACE the shipped"):
        run_bc2.campaign_targets(args, "camp", None)


def test_duplicate_target_names_are_refused():
    args = Args(target=["name=A;path=x.pdb", "name=A;path=y.pdb"])
    with pytest.raises(run_bc2.SettingsConfigError, match="used more than once"):
        run_bc2.campaign_targets(args, "camp", None)


def test_shipped_targets_accumulate():
    args = Args(shipped_target=["hPDL1", "mPDL1"])
    targets, shipped = run_bc2.campaign_targets(args, "camp", None)
    assert targets is None
    assert shipped == ["hPDL1", "mPDL1"]


# --------------------------------------------------------------------------
# 2. settings assembly
# --------------------------------------------------------------------------


def test_single_target_child_run_shape(tmp_path):
    """The single-target path is the one every existing campaign used: keep it."""
    target = tmp_path / "t.pdb"
    target.write_text("ATOM\n")
    settings = build(Args(chains="A", hotspots="A54,A56"), target_path=target)
    assert settings["targets"] == [
        {
            "name": "camp",
            "target_path": str(target),
            "chains": "A",
            "hotspots": "A54,A56",
        }
    ]
    assert settings["campaign_name"] == "camp"
    assert settings["resume"] is True
    assert settings["save_design_trajectory"] is True


def test_multi_target_with_detarget(tmp_path):
    for name in ("hPDL1", "mPDL1", "hPD1"):
        (tmp_path / f"{name}.pdb").write_text("ATOM\n")
    args = Args(
        target=[
            f"name=hPDL1;path={tmp_path}/hPDL1.pdb;chains=A;hotspots=A54,A56",
            f"name=mPDL1;path={tmp_path}/mPDL1.pdb;chains=A",
            f"name=hPD1;path={tmp_path}/hPD1.pdb;objective=detarget;weight=-0.5",
        ]
    )
    settings = build(args)
    assert [t["name"] for t in settings["targets"]] == ["hPDL1", "mPDL1", "hPD1"]
    # weight has to reach the file as a NUMBER: upstream does float(...) on it and a
    # quoted "-0.5" would be the first thing to notice.
    assert settings["targets"][2]["weight"] == -0.5
    assert isinstance(settings["targets"][2]["weight"], float)
    assert settings["targets"][2]["objective"] == "detarget"


def test_extra_target_is_appended_after_the_table_target(tmp_path):
    primary = tmp_path / "primary.pdb"
    off = tmp_path / "off.pdb"
    primary.write_text("ATOM\n")
    off.write_text("ATOM\n")
    args = Args(extra_target=[f"name=off;path={off};weight=-1"], hotspots="A10")
    settings = build(args, name="row1", target_path=primary)
    assert [t["name"] for t in settings["targets"]] == ["row1", "off"]
    assert settings["targets"][0]["hotspots"] == "A10"
    assert settings["targets"][1]["weight"] == -1.0


def test_targets_file_is_accepted(tmp_path):
    target = tmp_path / "t.pdb"
    target.write_text("ATOM\n")
    spec = tmp_path / "targets.yaml"
    spec.write_text(
        json.dumps({"targets": [{"name": "T", "target_path": str(target)}]})
    )
    settings = build(Args(targets_file=str(spec)))
    assert settings["targets"][0]["name"] == "T"


def test_extra_settings_may_own_targets_when_no_target_flag_is_used(tmp_path):
    """The escape hatch the old tool closed: --extra-settings supplying `targets`.

    With no target flag in play the tool does not claim `targets`, so there is no
    collision and the settings file's own list flows straight through.
    """
    settings = build(
        Args(), extra={"targets": [{"name": "T", "target_path": "/x/t.pdb"}]}
    )
    assert settings["targets"] == [{"name": "T", "target_path": "/x/t.pdb"}]


def test_extra_settings_targets_collides_with_a_table_target(tmp_path):
    target = tmp_path / "t.pdb"
    target.write_text("ATOM\n")
    with pytest.raises(run_bc2.SettingsConfigError, match="targets"):
        build(
            Args(),
            target_path=target,
            extra={"targets": [{"name": "T", "target_path": "/x/t.pdb"}]},
        )


def test_a_campaign_with_no_target_is_refused():
    with pytest.raises(run_bc2.SettingsConfigError, match="declares no target"):
        build(Args())


def test_execution_and_monomer_flags_reach_the_settings(tmp_path):
    target = tmp_path / "t.pdb"
    target.write_text("ATOM\n")
    settings = build(
        Args(design_workers=4, workers_per_gpu="2", save_monomers=True),
        target_path=target,
    )
    assert settings["design_workers"] == 4
    assert settings["workers_per_gpu"] == "2"
    assert settings["save_binder_monomers"] is True


@pytest.mark.parametrize(
    "spec, expected",
    [("80", [80]), ("60-100", [60, 100]), ("60,80,100", [60, 80, 100])],
)
def test_binder_lengths(spec, expected):
    assert run_bc2._coerce_binder_lengths(spec) == expected


# --------------------------------------------------------------------------
# 3. CLI tokens
# --------------------------------------------------------------------------


def test_cli_tokens_keep_values_whole(tmp_path):
    meta = tmp_path / "meta.json"
    meta.write_text("{}")
    tokens = run_bc2._cli_tokens(
        Args(
            core="benchmark",
            modality="VHH,induced_fit",
            property=["humanize", "termini_accessible"],
            metadata=str(meta),
            # The case the old tool mangled: a value with spaces and braces.
            set=['binder_lengths=[70, 90]', 'aa_bias={"C": 0}'],
        )
    )
    assert tokens == [
        "--core",
        "benchmark",
        "--modality",
        "VHH,induced_fit",
        "--humanize",
        "--termini-accessible",
        "--metadata",
        str(meta),
        "--set",
        "binder_lengths=[70, 90]",
        "--set",
        'aa_bias={"C": 0}',
    ]


def test_cli_tokens_refuse_a_newline():
    with pytest.raises(run_bc2.SettingsConfigError, match="newline"):
        run_bc2._cli_tokens(Args(set=["a=b\nc=d"]))


def test_missing_metadata_file_is_refused(tmp_path):
    with pytest.raises(run_bc2.SettingsConfigError, match="--metadata file not found"):
        run_bc2._cli_tokens(Args(metadata=str(tmp_path / "nope.json")))


# --------------------------------------------------------------------------
# 4. the settings files upstream actually has to accept
# --------------------------------------------------------------------------

_UPSTREAM_CHECK = r"""
import json, sys
from bindcraft.settings import (
    load_settings, reject_unrecognized_settings, campaign_over_presets,
    resolve_prepared_states,
)
settings = json.load(open(sys.argv[1]))
merged = campaign_over_presets(settings)
reject_unrecognized_settings(merged)
load_settings(settings)
states = resolve_prepared_states(merged)
print(json.dumps([
    {"name": s.name, "weight": s.weight, "objective": s.objective,
     "hotspots": s.hotspots, "chains": s.chains}
    for s in states
]))
"""


def upstream_resolve(settings: dict, tmp_path: Path) -> list[dict]:
    """Validate a settings mapping with BindCraft2 itself, and return its targets.

    Raises with upstream's own message if it refuses the file -- which is the point:
    this is the check that catches a key we spelled wrong or a type we got wrong,
    rather than our own idea of what BindCraft2 accepts.
    """
    settings_path = tmp_path / "settings.json"
    settings_path.write_text(json.dumps(settings))
    check = tmp_path / "check.py"
    check.write_text(_UPSTREAM_CHECK)
    result = subprocess.run(
        [str(UPSTREAM_PYTHON), str(check), str(settings_path)],
        capture_output=True,
        text=True,
        cwd=UPSTREAM,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"BindCraft2 refused the settings this tool generated:\n{result.stderr}"
        )
    return json.loads(result.stdout.strip().splitlines()[-1])


@needs_upstream
def test_generated_single_target_settings_validate(tmp_path):
    target = tmp_path / "t.pdb"
    target.write_text("ATOM\n")
    settings = build(
        Args(chains="A", hotspots="A54,A56", binder_lengths="60-100", num_designs=5),
        target_path=target,
    )
    states = upstream_resolve(settings, tmp_path)
    assert [s["name"] for s in states] == ["camp"]
    assert states[0]["hotspots"] == "A54,A56"
    assert states[0]["objective"] == "target"


@needs_upstream
def test_generated_multi_target_settings_validate_and_detarget(tmp_path):
    """The campaign the old tool could not express at all, end to end through
    upstream's own resolver: two orthologs bound, one off-target avoided."""
    for name in ("hPDL1", "mPDL1", "hPD1"):
        (tmp_path / f"{name}.pdb").write_text("ATOM\n")
    args = Args(
        target=[
            f"name=hPDL1;path={tmp_path}/hPDL1.pdb;chains=A;hotspots=A54,A56",
            f"name=mPDL1;path={tmp_path}/mPDL1.pdb;chains=A;hotspots=A36,A38",
            f"name=hPD1;path={tmp_path}/hPD1.pdb;objective=detarget;weight=-0.5",
        ],
        num_designs=10,
    )
    states = upstream_resolve(build(args), tmp_path)
    by_name = {s["name"]: s for s in states}
    assert set(by_name) == {"hPDL1", "mPDL1", "hPD1"}
    assert by_name["hPDL1"]["objective"] == "target"
    assert by_name["mPDL1"]["hotspots"] == "A36,A38"
    # a negative weight selects detargeting, and upstream normalises the sign
    assert by_name["hPD1"]["objective"] == "detarget"
    assert by_name["hPD1"]["weight"] == -0.5


@needs_upstream
def test_generated_shipped_target_settings_validate(tmp_path):
    settings = build(Args(shipped_target=["hPDL1", "mPDL1"], num_designs=2))
    states = upstream_resolve(settings, tmp_path)
    assert [s["name"] for s in states] == ["hPDL1", "mPDL1"]


@needs_upstream
def test_execution_settings_are_names_upstream_knows(tmp_path):
    """design_workers / workers_per_gpu / save_binder_monomers are real setting
    names, not ones this tool invented. reject_unrecognized_settings is the judge."""
    target = tmp_path / "t.pdb"
    target.write_text("ATOM\n")
    settings = build(
        Args(design_workers=2, workers_per_gpu="auto", save_monomers=True),
        target_path=target,
    )
    upstream_resolve(settings, tmp_path)


# --------------------------------------------------------------------------
# 5. collector: multi-target
# --------------------------------------------------------------------------

import pandas as pd  # noqa: E402  (after the module loader above)


def row(**fields) -> pd.Series:
    return pd.Series(fields)


def test_row_targets_parses_names_and_weights():
    r = row(targets="hPDL1;mPDL1;hPD1", target_weights="1;1;-0.5")
    assert collect_bc2.row_targets(r) == [("hPDL1", 1.0), ("mPDL1", 1.0), ("hPD1", -0.5)]


def test_row_targets_is_empty_for_a_single_target_campaign():
    """BindCraft2 writes neither column below two targets, so there is nothing to
    split and the single-target table keeps exactly its old shape."""
    assert collect_bc2.row_targets(row(i_pTM=0.8)) == []


def test_row_targets_survives_a_weight_mismatch():
    r = row(targets="A;B", target_weights="1")
    assert collect_bc2.row_targets(r) == [("A", 1.0), ("B", 1.0)]


def test_split_target_columns():
    targets = [("hPDL1", 1.0), ("mPDL1", 1.0), ("hPD1", -0.5)]
    data = {"i_pTM": "0.82;0.79;0.21", "outcome": "passed"}
    derived = collect_bc2.split_target_columns(data, targets)
    assert derived["i_pTM__hPDL1"] == 0.82
    assert derived["i_pTM__mPDL1"] == 0.79
    assert derived["i_pTM__hPD1"] == 0.21
    # a cell with the wrong number of positions is not mapped onto the targets
    assert not any(k.startswith("outcome__") for k in derived)


def test_split_skips_a_cell_whose_positions_do_not_match():
    """The guard that stops a semicolon inside a residue list being read as a
    per-target packing. Same rule as campaign_output.per_target_readings."""
    targets = [("A", 1.0), ("B", 1.0), ("C", 1.0)]
    derived = collect_bc2.split_target_columns({"Interface_Residues": "7;9"}, targets)
    assert derived == {}


def test_summaries_use_binding_targets_only():
    targets = [("A", 1.0), ("B", 1.0), ("off", -1.0)]
    data = {"i_pTM": "0.90;0.70;0.20"}
    derived = collect_bc2.split_target_columns(data, targets)
    # mean/worst/best ignore the off-target: averaging it in would reward binding
    # the thing the campaign is avoiding.
    assert derived["i_pTM__mean"] == pytest.approx(0.80)
    assert derived["i_pTM__worst"] == pytest.approx(0.70)
    assert derived["i_pTM__best"] == pytest.approx(0.90)
    assert derived["i_pTM__spread"] == pytest.approx(0.20)
    # weakest binder minus strongest off-target, positive = the intended distinction
    assert derived["i_pTM__selectivity"] == pytest.approx(0.50)


def test_selectivity_flips_for_a_lower_is_better_metric():
    """i_pAE is lower-is-better, so the off-target should read HIGHER than the
    binding targets and the margin is computed the other way round."""
    targets = [("A", 1.0), ("off", -1.0)]
    derived = collect_bc2.split_target_columns({"i_pAE": "0.20;0.60"}, targets)
    assert derived["i_pAE__worst"] == pytest.approx(0.20)  # only one binder
    assert derived["i_pAE__selectivity"] == pytest.approx(0.40)


def test_no_summary_for_a_metric_with_no_known_direction():
    targets = [("A", 1.0), ("B", 1.0)]
    derived = collect_bc2.split_target_columns({"Made_Up_Metric": "1;2"}, targets)
    assert derived == {"Made_Up_Metric__A": 1.0, "Made_Up_Metric__B": 2.0}


# --------------------------------------------------------------------------
# 6. collector: structures
# --------------------------------------------------------------------------


def test_find_structures_picks_the_highest_weight_binding_target(tmp_path):
    for suffix in ("hPDL1", "mPDL1", "hPD1"):
        (tmp_path / f"d1_{suffix}.cif").write_text("x")
    (tmp_path / "d1_monomer.cif").write_text("x")
    # ordered as BindCraft2 orders them: by (-weight, name)
    targets = [("hPDL1", 2.0), ("mPDL1", 1.0), ("hPD1", -0.5)]
    found = collect_bc2.find_structures(tmp_path, "d1", targets)
    assert found.primary == tmp_path / "d1_hPDL1.cif"
    assert set(found.per_target) == {"hPDL1", "mPDL1", "hPD1"}
    assert found.monomer == tmp_path / "d1_monomer.cif"
    assert found.fell_back is False


def test_find_structures_single_target(tmp_path):
    (tmp_path / "d1.cif").write_text("x")
    found = collect_bc2.find_structures(tmp_path, "d1", [])
    assert found.primary == tmp_path / "d1.cif"
    assert found.per_target == {}
    assert found.fell_back is False


def test_find_structures_never_returns_a_monomer_as_the_complex(tmp_path):
    """The free binder is a different structure; handing it back as `<leaf>_path`
    would silently turn an interface metric into a monomer metric."""
    (tmp_path / "d1_monomer.cif").write_text("x")
    found = collect_bc2.find_structures(tmp_path, "d1", [])
    assert found.primary is None
    assert found.monomer == tmp_path / "d1_monomer.cif"


def test_find_structures_flags_a_glob_fallback(tmp_path):
    (tmp_path / "d1_unexpectedname.cif").write_text("x")
    found = collect_bc2.find_structures(tmp_path, "d1", [])
    assert found.primary == tmp_path / "d1_unexpectedname.cif"
    assert found.fell_back is True


# --------------------------------------------------------------------------
# 7. collector: archived trajectories
# --------------------------------------------------------------------------


def _make_archive(stage_dir: Path, design: str) -> Path:
    """A trajectory archive exactly as archive_trajectory_folder writes one:
    members relative to the trajectories folder, i.e. '<design>/<file>'."""
    stage_dir.mkdir(parents=True, exist_ok=True)
    archive_path = stage_dir / f"{design}.zip"
    with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(f"{design}/{design}_trajectory.cif", "data_x\n")
        archive.writestr(f"{design}/{design}_losses.csv", "phase,round\n")
    return archive_path


def test_archived_trajectory_is_read_out_of_the_zip(tmp_path):
    stage_dir = tmp_path / "1_Trajectories"
    _make_archive(stage_dir, "d1")
    dest = tmp_path / "unarchived"
    restored = collect_bc2.unarchived_dir(stage_dir, "d1", dest)
    assert restored == dest / "d1"
    assert (restored / "d1_trajectory.cif").is_file()
    # the campaign folder itself is untouched: the archive is still there
    assert (stage_dir / "d1.zip").is_file()
    found = collect_bc2.find_structures(restored, "d1_trajectory", [])
    assert found.primary == restored / "d1_trajectory.cif"


def test_unarchive_is_idempotent(tmp_path):
    stage_dir = tmp_path / "1_Trajectories"
    _make_archive(stage_dir, "d1")
    dest = tmp_path / "unarchived"
    first = collect_bc2.unarchived_dir(stage_dir, "d1", dest)
    second = collect_bc2.unarchived_dir(stage_dir, "d1", dest)
    assert first == second


def test_unarchive_refuses_a_traversing_member(tmp_path):
    stage_dir = tmp_path / "1_Trajectories"
    stage_dir.mkdir(parents=True)
    archive_path = stage_dir / "d1.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("d1/../../escaped.cif", "x")
        archive.writestr("d1/d1_trajectory.cif", "data_x\n")
    dest = tmp_path / "unarchived"
    collect_bc2.unarchived_dir(stage_dir, "d1", dest)
    assert not (tmp_path.parent / "escaped.cif").exists()
    assert not (tmp_path / "escaped.cif").exists()
    assert (dest / "d1" / "d1_trajectory.cif").is_file()


def test_no_archive_returns_none(tmp_path):
    stage_dir = tmp_path / "1_Trajectories"
    stage_dir.mkdir(parents=True)
    assert collect_bc2.unarchived_dir(stage_dir, "d1", tmp_path / "u") is None


# --------------------------------------------------------------------------
# 8. the task script, end to end
# --------------------------------------------------------------------------

PRELUDE = (
    REPO
    / ".venv/lib/python3.13/site-packages/prosapia/core/scripts/sapia_task_prelude.sh"
)
needs_prelude = pytest.mark.skipif(
    not PRELUDE.is_file() or shutil.which("bash") is None,
    reason="needs bash and the installed prosapia task prelude",
)


def _run_task(tmp_path: Path, line: str, args_tokens: list[str] | None = None):
    """Run bindcraft2.sh for one manifest line against a fake `bindcraft` that
    dumps its argv, one entry per line."""
    fake = tmp_path / "fake-bindcraft"
    fake.write_text('#!/bin/bash\nfor a in "$@"; do echo "ARG:$a"; done\n')
    fake.chmod(0o755)

    manifest = tmp_path / "manifest.txt"
    manifest.write_text(line + "\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir(exist_ok=True)

    if args_tokens is not None:
        (tmp_path / "design.args").write_text(
            "".join(f"{t}\n" for t in args_tokens)
        )

    return subprocess.run(
        ["bash", str(TOOL_DIR / "bindcraft2.sh"), str(manifest), str(out_dir)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        env={
            **os.environ,
            "SAPIA_PRELUDE": str(PRELUDE),
            "SAPIA_SCHEDULER": "modal",  # makes sapia_activate a no-op
            "SAPIA_TASK_ID": "1",
            "BINDCRAFT_BIN": str(fake),
        },
    )


@needs_prelude
def test_task_script_keeps_a_value_with_spaces_as_one_argument(tmp_path):
    """The regression this rewrite exists to prevent.

    The old script word-split an unquoted flags string, so `--set
    'binder_lengths=[70, 90]'` reached bindcraft as two arguments and the campaign
    either died or, worse, read a truncated value.
    """
    tokens = ["--set", "binder_lengths=[70, 90]", "--set", 'aa_bias={"C": 0}']
    settings = tmp_path / "camp.json"
    settings.write_text("{}")
    line = f"camp\tdesign\t{settings}\t{tmp_path / 'design.args'}"
    result = _run_task(tmp_path, line, tokens)
    assert result.returncode == 0, result.stderr
    passed = [
        l[len("ARG:") :] for l in result.stdout.splitlines() if l.startswith("ARG:")
    ]
    assert passed == ["design", str(settings), *tokens]


@needs_prelude
def test_task_script_handles_no_flags(tmp_path):
    settings = tmp_path / "camp.json"
    settings.write_text("{}")
    line = f"camp\tdesign\t{settings}\t{tmp_path / 'design.args'}"
    result = _run_task(tmp_path, line, [])
    assert result.returncode == 0, result.stderr
    passed = [
        l[len("ARG:") :] for l in result.stdout.splitlines() if l.startswith("ARG:")
    ]
    assert passed == ["design", str(settings)]


@needs_prelude
def test_task_script_fetch_weights_mode(tmp_path):
    result = _run_task(tmp_path, "bindcraft2_weights\tfetch-weights\t-\t-")
    assert result.returncode == 0, result.stderr
    passed = [
        l[len("ARG:") :] for l in result.stdout.splitlines() if l.startswith("ARG:")
    ]
    assert passed == ["fetch-weights"]


# --------------------------------------------------------------------------
# 9. the manifest builder itself, and its contract with the task script
# --------------------------------------------------------------------------


class Ctx:
    """A stand-in for prosapia's ManifestCtx, enough for a root run."""

    def __init__(self, args, out_dir: Path):
        self.args = args
        self.out_dir = out_dir
        self.df = pd.DataFrame()
        self.lookup = no_lookup
        self.meta: dict = {}

    def write_meta(self, **fields):
        self.meta.update(fields)

    @property
    def ready(self):
        return self.df


def test_manifest_builder_writes_settings_and_args(tmp_path):
    for name in ("hPDL1", "hPD1"):
        (tmp_path / f"{name}.pdb").write_text("ATOM\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    args = Args(
        run_dir=tmp_path,
        target=[
            f"name=hPDL1;path={tmp_path}/hPDL1.pdb;chains=A;hotspots=A54,A56",
            f"name=hPD1;path={tmp_path}/hPD1.pdb;objective=detarget;weight=-0.5",
        ],
        num_designs=3,
        set=["binder_lengths=[70, 90]"],
    )
    ctx = Ctx(args, out_dir)
    rows = run_bc2.build_bindcraft2_manifest(ctx)

    assert len(rows) == 1
    name, mode, settings_path, args_path = rows[0]
    assert name == "hPDL1_bc2"  # named after the BINDING target, not the off-target
    assert mode == run_bc2.MODE_DESIGN

    settings = json.loads(Path(settings_path).read_text())
    assert [t["name"] for t in settings["targets"]] == ["hPDL1", "hPD1"]
    assert settings["number_of_final_designs"] == 3

    # the args file is one argv token per line, so the spaced value stays whole
    assert Path(args_path).read_text() == "--set\nbinder_lengths=[70, 90]\n"

    # the sidecar is what collect reads to find the campaigns and pick the stage
    assert ctx.meta["trajectory_only"] is False
    assert ctx.meta["root_designs"] == ["hPDL1_bc2"]
    assert ctx.meta["campaigns_root"].endswith("campaigns")


def test_fetch_weights_only_builds_one_cpu_task(tmp_path):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    args = Args(run_dir=tmp_path, fetch_weights_only=True, gpus_per_task=1)
    ctx = Ctx(args, out_dir)
    rows = run_bc2.build_bindcraft2_manifest(ctx)
    assert rows == [("bindcraft2_weights", run_bc2.MODE_FETCH_WEIGHTS, "-", "-")]
    assert args.gpus_per_task == 0
    assert ctx.meta["weights_only"] is True


def test_fetch_weights_only_refuses_campaign_flags(tmp_path):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    args = Args(
        run_dir=tmp_path, fetch_weights_only=True, gpus_per_task=1, num_designs=5
    )
    with pytest.raises(ValueError, match="runs no campaign"):
        run_bc2.build_bindcraft2_manifest(Ctx(args, out_dir))


def test_manifest_arity_matches_the_task_script():
    """The manifest and bindcraft2.sh agree on how many fields a line has.

    A mismatch here is the worst failure signature there is: `cut -f4` on a 3-field
    line yields an empty string, the script runs anyway, and the campaign dies (or
    silently drops its flags) with nothing in the logs explaining why.
    """
    script = (TOOL_DIR / "bindcraft2.sh").read_text()
    cut_fields = {
        int(part)
        for part in __import__("re").findall(r"cut -f(\d+)", script)
    }
    assert cut_fields == {1, 2, 3, 4}


def test_collector_reports_nothing_for_a_weights_run(tmp_path, capsys):
    """A --fetch-weights-only run reserves a table it must not try to fill.

    The sidecar filename comes from prosapia, not from a guess here: an earlier
    version of this test invented one, so the branch never ran and the test passed
    while proving nothing.
    """
    from prosapia.core import RUN_META_FILENAME

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / RUN_META_FILENAME).write_text(json.dumps({"weights_only": True}))

    class Stub:
        class args:
            run_dir = tmp_path
            stage = collect_bc2.STAGE_AUTO
            metrics = ""
            all_metrics = False
            no_split_targets = False

    ctx = Stub()
    ctx.out_dir = out_dir

    class Design:
        name = "anything"

    one = collect_bc2.collect_bindcraft2(ctx)
    assert list(one(Design())) == []
    # both branches yield no rows, so the row count proves nothing. What
    # distinguishes them is that the weights run never looks for a stage at all.
    said = capsys.readouterr().out
    assert "nothing to collect" in said
    assert "Collecting" not in said


def test_collector_does_not_early_out_for_an_ordinary_run(tmp_path, capsys):
    """The other half of the branch: without weights_only it collects normally."""
    from prosapia.core import RUN_META_FILENAME

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / RUN_META_FILENAME).write_text(json.dumps({"trajectory_only": False}))

    class Stub:
        class args:
            run_dir = tmp_path
            stage = collect_bc2.STAGE_AUTO
            metrics = ""
            all_metrics = False
            no_split_targets = False

    ctx = Stub()
    ctx.out_dir = out_dir

    class Design:
        name = "missing_campaign"

    one = collect_bc2.collect_bindcraft2(ctx)
    # no campaign dir on disk -> no rows, but it got as far as looking
    assert list(one(Design())) == []
    said = capsys.readouterr().out
    assert "Collecting 'ranked' designs" in said
    assert "nothing to collect" not in said
