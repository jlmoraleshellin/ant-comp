#!/usr/bin/env python3
"""
Self-test for bindcraft2's --targets flag and its multi-target collector.

BindCraft2 itself is not installed where `sapia run` executes, so the two things this
test covers are exactly the two that can only be got wrong here:

    * the SUBMIT side -- that --targets validates against BindCraft2's seven per-target
      keys, refuses the flags that would give the epitope two sources, resolves {expr}
      and target_path per campaign, and leaves the single-target settings file
      BYTE-IDENTICAL to what this tool has always written;
    * the COLLECT side -- that a `;` cell is split by NAME (never by position or by
      file order), that a detarget's columns are kept apart from a binder's, that an
      empty position is NA and not 0, that a count mismatch is an ERROR status rather
      than a best-effort parse, that a missing per-target complex is NA and never a
      substituted sibling, and that a single-target campaign still collects exactly as
      before.

Everything is built from fabricated campaign folders and CSVs written in upstream's own
layout (campaign_output.py at v1.0.3) -- no GPU, no bindcraft, no network.

Run:  uv run python tools/bindcraft2/test_bindcraft2.py
"""

import importlib.util
import json
import sys
import tempfile
from argparse import Namespace
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


run_bc2 = _load("run_bindcraft2")
collect_bc2 = _load("collect_bindcraft2")

from prosapia.core import CollectCtx, DesignCtx  # noqa: E402

FAILURES: list[str] = []


def check(condition: bool, message: str) -> None:
    print(("  ok   " if condition else "  FAIL ") + message)
    if not condition:
        FAILURES.append(message)


def raises(fn, fragment: str, message: str) -> None:
    try:
        fn()
    except Exception as e:  # noqa: BLE001 - the point is what was raised
        check(fragment in str(e), f"{message} (said: {str(e)[:140]})")
        return
    check(False, f"{message} -- nothing raised")


# ---------------------------------------------------------------- fixtures


def write_targets(tmp: Path, entries, name: str = "egfr_pair.yaml") -> Path:
    path = tmp / name
    path.write_text(json.dumps(entries, indent=2))
    return path


def touch_pdb(tmp: Path, stem: str) -> Path:
    path = tmp / f"{stem}.pdb"
    path.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000\nEND\n")
    return path


MINIMAL_CIF = """data_x
_bindcraft.bindcraft_version 1.0.3
_bindcraft.redesigned_residues {spans}
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 CA ALA A 1 0.000 0.000 0.000
ATOM 2 CA ALA B 1 3.800 0.000 0.000
"""


def write_cif(path: Path, spans: str = "B1-20") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(MINIMAL_CIF.format(spans=spans))
    return path


def make_campaign(
    root: Path, group: str, rows: list[dict], complexes: list[str], spans: str = "B1-20"
) -> None:
    """A fabricated 3_Ranked stage: its CSV and the complexes it claims to have."""
    ranked = root / group / "3_Ranked"
    ranked.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(ranked / "!_Ranked.csv", index=False)
    for stem in complexes:
        write_cif(ranked / f"{stem}.cif", spans)


def collect_rows(out_dir: Path, groups: list[str], **args) -> dict[str, dict]:
    """Run the collector over ``groups`` and return the emitted rows by name."""
    ctx = CollectCtx(
        df=pd.DataFrame(),
        args=Namespace(
            run_dir=out_dir.parent,
            table="table1",
            dir_label="",
            force=False,
            stage=args.pop("stage", "ranked"),
            metrics=args.pop("metrics", ""),
            all_metrics=args.pop("all_metrics", False),
        ),
        table_name="table1",
        out_dir=out_dir,
        status_col="bindcraft2_status",
        path_col="bindcraft2_path",
        parent_table="table0",
        parent_df=pd.DataFrame(index=pd.Index(groups)),
        lookup=lambda name, column: None,
        creates_table=True,
        default_input_column="pdb_path",
    )
    one = collect_bc2.collect_bindcraft2(ctx)
    emitted: dict[str, dict] = {}
    for group in groups:
        for c in one(DesignCtx(group, out_dir, ctx.lookup)):
            emitted[c.name or group] = {
                "data": dict(c.data),
                "path": c.path,
                "status": c.status,
                "parent": c.parent,
            }
    return emitted


# ---------------------------------------------------------------- submit side


def t_schema_validation(tmp: Path) -> None:
    good = [
        {"name": "hEGFR", "target_path": "h.pdb", "chains": "A", "hotspots": "A355"},
        {"name": "mEGFR", "target_path": "m.pdb", "chains": "A", "weight": 1.0},
        {"name": "hERBB2", "target_path": "e.pdb", "weight": -0.5},
    ]
    for stem in ("h", "m", "e"):
        touch_pdb(tmp, stem)
    spec = run_bc2.load_targets(str(write_targets(tmp, good)))
    check(len(spec.entries) == 3, "a valid 3-target file parses")
    check(spec.stem == "egfr_pair", "the design group derives from the file stem")

    wrapped = write_targets(tmp, {"targets": good}, "wrapped.json")
    check(
        len(run_bc2.load_targets(str(wrapped)).entries) == 3,
        "the {targets: [...]} shape from upstream's own examples is accepted too",
    )

    raises(
        lambda: run_bc2.validate_target_entries(
            [{"name": "a", "target_path": "x.pdb", "hotspot": "A1"}]
        ),
        "did you mean 'hotspots'",
        "an unknown per-target key is refused at submit time, with a suggestion",
    )
    raises(
        lambda: run_bc2.validate_target_entries([{"target_path": "x.pdb"}]),
        "name is required",
        "a target with no name is refused (upstream would KeyError per container)",
    )
    raises(
        lambda: run_bc2.validate_target_entries([{"name": "a"}]),
        "target_path is required",
        "a target with no target_path is refused",
    )
    raises(
        lambda: run_bc2.validate_target_entries(
            [{"name": "h.EGFR", "target_path": "x.pdb"}]
        ),
        "metric suffix",
        "a '.' in a target name is refused (it is upstream's state separator)",
    )
    raises(
        lambda: run_bc2.validate_target_entries(
            [{"name": "h;EGFR", "target_path": "x.pdb"}]
        ),
        "must be letters/digits",
        "a ';' in a target name is refused (it is the value separator)",
    )
    raises(
        lambda: run_bc2.validate_target_entries(
            [{"name": "off_x", "target_path": "x.pdb"}]
        ),
        "off_",
        "a target named off_* is refused (it would read as a detarget column)",
    )
    raises(
        lambda: run_bc2.validate_target_entries(
            [
                {"name": "a", "target_path": "x.pdb"},
                {"name": "a", "target_path": "y.pdb"},
            ]
        ),
        "repeats",
        "duplicate target names are refused",
    )
    raises(
        lambda: run_bc2.validate_target_entries(
            [{"name": "a", "target_path": "x.pdb", "objective": "avoid"}]
        ),
        "objective",
        "an objective other than target/detarget is refused",
    )
    raises(
        lambda: run_bc2.validate_target_entries(
            [{"name": "a", "target_path": "x.pdb", "weight": -1}]
        ),
        "nothing to bind",
        "a list of nothing but detargets is refused",
    )
    raises(
        lambda: run_bc2.load_targets(str(tmp / "nope.yaml")),
        "not found",
        "a missing --targets file is refused",
    )


def t_exclusivity() -> None:
    def args(**kw):
        base = dict(
            targets="t.yaml",
            target_pdb=None,
            shipped_target=None,
            chains=None,
            hotspots=None,
            coldspots=None,
        )
        base.update(kw)
        return Namespace(**base)

    raises(
        lambda: run_bc2.check_targets_exclusivity(args(target_pdb=Path("x.pdb"))),
        "--target-pdb",
        "--targets and --target-pdb are mutually exclusive",
    )
    raises(
        lambda: run_bc2.check_targets_exclusivity(args(shipped_target="hPDL1")),
        "--shipped-target",
        "--targets and --shipped-target are mutually exclusive",
    )
    for flag, kw in (
        ("--hotspots", {"hotspots": "A54"}),
        ("--chains", {"chains": "A"}),
        ("--coldspots", {"coldspots": "A10"}),
    ):
        raises(
            lambda kw=kw: run_bc2.check_targets_exclusivity(args(**kw)),
            "PER-TARGET",
            f"{flag} is refused with --targets, naming it a per-target sub-key",
        )
    run_bc2.check_targets_exclusivity(args(targets=None, hotspots="A54"))
    check(True, "without --targets the top-level --hotspots is still fine")


def t_target_path_and_expr(tmp: Path) -> None:
    structures = tmp / "structures"
    structures.mkdir(exist_ok=True)
    touch_pdb(structures, "hEGFR")
    spec_dir = tmp / "spec"
    spec_dir.mkdir(exist_ok=True)
    touch_pdb(spec_dir, "beside")

    resolved = run_bc2.resolve_target_path(
        str(structures / "hEGFR.pdb"), spec_dir, "targets[0]"
    )
    check(resolved.is_absolute() and resolved.exists(), "an absolute target_path passes")
    check(
        run_bc2.resolve_target_path("beside.pdb", spec_dir, "targets[0]").exists(),
        "a relative target_path falls back to the --targets file's own directory",
    )
    raises(
        lambda: run_bc2.resolve_target_path("ghost.pdb", spec_dir, "targets[0]"),
        "does not exist",
        "a target_path that is nowhere fails at SUBMIT time, naming both attempts",
    )

    spec = run_bc2.TargetsSpec(
        [
            {
                "name": "hEGFR",
                "target_path": "beside.pdb",
                "hotspots": "A{epitope_start}-{epitope_end}",
            }
        ],
        spec_dir,
        "pair",
    )

    def lookup(name, column):
        return {"epitope_start": 355, "epitope_end": 360}[column]

    block = run_bc2.build_targets_block(spec, lookup, "row0")
    check(
        block[0]["hotspots"] == "A355-360",
        "{expr} in a --targets value resolves up the lineage, per campaign",
    )
    check(
        Path(block[0]["target_path"]).is_absolute(),
        "every target_path in the written settings is absolute",
    )


def t_settings_regression(tmp: Path) -> None:
    """The single-target settings file must not move a byte."""
    target = touch_pdb(tmp, "PDL1")
    args = Namespace(
        targets=None,
        target_pdb=target,
        shipped_target=None,
        chains="A",
        hotspots="A54,A56",
        coldspots=None,
        binder_lengths="60-100",
        num_designs=10,
        max_trajectories=200,
        campaign_seed=None,
        trajectory_only=False,
        no_resume=False,
    )
    built = run_bc2._build_settings(
        "PDL1_bc2", target, args, lambda n, c: None, tmp / "camp", {}
    )
    expected = {
        "save_design_trajectory": True,
        "campaign_name": "PDL1_bc2",
        "project_folder": str(tmp / "camp"),
        "resume": True,
        "targets": [
            {"name": "PDL1_bc2", "target_path": str(target), "chains": "A",
             "hotspots": "A54,A56"}
        ],
        "binder_lengths": [60, 100],
        "number_of_final_designs": 10,
        "max_trajectories": 200,
    }
    check(built == expected, "single-target settings are unchanged by this edit")

    multi = run_bc2._build_settings(
        "pair_bc2",
        None,
        Namespace(**{**vars(args), "target_pdb": None, "chains": None, "hotspots": None}),
        lambda n, c: None,
        tmp / "camp",
        {},
        run_bc2.TargetsSpec(
            [
                {"name": "hEGFR", "target_path": str(target), "weight": 1.0},
                {"name": "hERBB2", "target_path": str(target), "weight": -0.5},
            ],
            tmp,
            "pair",
        ),
    )
    check(
        [t["name"] for t in multi["targets"]] == ["hEGFR", "hERBB2"]
        and multi["campaign_name"] == "pair_bc2",
        "a --targets run writes BindCraft2's own targets list, in file order",
    )
    raises(
        lambda: run_bc2._build_settings(
            "pair_bc2",
            None,
            args,
            lambda n, c: None,
            tmp / "camp",
            {"targets": [{"name": "x", "target_path": "y"}]},
            run_bc2.TargetsSpec(
                [{"name": "hEGFR", "target_path": str(target)}], tmp, "pair"
            ),
        ),
        "--extra-settings",
        "a `targets:` in --extra-settings still collides with --targets",
    )


def _manifest_args(tmp: Path, **kw) -> "Namespace":
    base = dict(
        run_dir=tmp,
        table=None,
        input_column="pdb_path",
        force=False,
        targets=None,
        target_pdb=None,
        shipped_target=None,
        chains=None,
        hotspots=None,
        coldspots=None,
        binder_lengths=None,
        num_designs=None,
        max_trajectories=None,
        campaign_seed=None,
        modality=None,
        core=None,
        property=[],
        set=[],
        trajectory_only=False,
        no_resume=False,
        reuse_campaigns=None,
        extra_settings=None,
    )
    base.update(kw)
    return Namespace(**base)


def t_manifest_end_to_end(tmp: Path) -> None:
    """The whole submit path: groups, settings files, manifest rows, sidecar."""
    from prosapia.core import ManifestCtx

    work = tmp / "manifest"
    (work / "out").mkdir(parents=True, exist_ok=True)
    touch_pdb(work, "hEGFR")
    touch_pdb(work, "mEGFR")
    targets_file = write_targets(
        work,
        [
            {"name": "hEGFR", "target_path": str(work / "hEGFR.pdb"), "weight": 1.0},
            {"name": "mEGFR", "target_path": str(work / "mEGFR.pdb"), "weight": 1.0},
        ],
        "egfr_pair.yaml",
    )
    meta: dict = {}
    out_dir = work / "out"

    # root run: one campaign, group named for the file stem
    rows = run_bc2.build_bindcraft2_manifest(
        ManifestCtx(
            df=pd.DataFrame(),
            args=_manifest_args(work, targets=str(targets_file)),
            out_dir=out_dir,
            lookup=lambda n, c: None,
            write_meta=lambda **kw: meta.update(kw),
        )
    )
    check(len(rows) == 1 and rows[0][0] == "egfr_pair_bc2",
          "a root --targets run is ONE campaign, named `<file stem>_bc2`")
    settings = json.loads((out_dir / "settings" / "egfr_pair_bc2.json").read_text())
    check([t["name"] for t in settings["targets"]] == ["hEGFR", "mEGFR"],
          "its settings file carries both targets in one campaign")
    check(meta.get("submitted_targets") == ["hEGFR", "mEGFR"],
          "the run records the submitted target names for collect to cross-check")

    # child run: one campaign per ready row, same target list
    parent = pd.DataFrame(
        {"pdb_path": [str(work / "hEGFR.pdb"), str(work / "mEGFR.pdb")]},
        index=pd.Index(["row0", "row1"], name="name"),
    )
    rows = run_bc2.build_bindcraft2_manifest(
        ManifestCtx(
            df=parent,
            args=_manifest_args(work, table="table0", targets=str(targets_file)),
            out_dir=out_dir,
            lookup=lambda n, c: None,
            write_meta=lambda **kw: meta.update(kw),
        )
    )
    check([r[0] for r in rows] == ["row0", "row1"],
          "a child --targets run is one campaign per ready ROW, all on the same targets")
    check(
        [t["name"] for t in json.loads(
            (out_dir / "settings" / "row0.json").read_text())["targets"]]
        == ["hEGFR", "mEGFR"],
        "and the row's own --input-column is not used as a target",
    )

    # a placeholder in a root --targets run has nothing to resolve against
    placeholder = write_targets(
        work,
        [{"name": "hEGFR", "target_path": str(work / "hEGFR.pdb"),
          "hotspots": "A{epitope_start}"}],
        "ph.yaml",
    )
    raises(
        lambda: run_bc2.build_bindcraft2_manifest(
            ManifestCtx(
                df=pd.DataFrame(),
                args=_manifest_args(work, targets=str(placeholder)),
                out_dir=out_dir,
                lookup=lambda n, c: None,
                write_meta=lambda **kw: None,
            )
        ),
        "root run",
        "a {expr} in --targets on a ROOT run is refused, not resolved to nonsense",
    )


# ---------------------------------------------------------------- collect side


def t_single_target_collect(tmp: Path) -> None:
    out_dir = tmp / "table1" / "bindcraft2"
    root = out_dir / "campaigns"
    make_campaign(
        root,
        "PDL1_bc2",
        [
            {
                "design": "PDL1_bc2_l80_ab12_seq0",
                "Binder_Sequence": "AAAA",
                "i_pTM": 0.82,
                "i_pAE": 8.1,
                "Interface_Residues": 14,
                "rank": 1,
                "hash": "ab12",
            }
        ],
        ["PDL1_bc2_l80_ab12_seq0", "PDL1_bc2_l80_ab12_seq0_monomer"],
    )
    rows = collect_rows(out_dir, ["PDL1_bc2"])
    check(len(rows) == 1, "a single-target campaign yields one row per design")
    row = next(iter(rows.values()))
    data = row["data"]
    check(row["status"] == "OK", "its status is OK")
    check(data["i_pTM"] == 0.82 and data["Interface_Residues"] == 14,
          "metrics land in bare columns, unsplit and unrenamed")
    check("targets" not in data and data["n_targets"] == 1,
          "no `targets` column is invented for a single-target campaign")
    check(str(row["path"]).endswith(".pdb"), "the complex is converted to PDB")
    check(data["binder_chain"] == "B" and data["binder_chain_src"] == "stamp",
          "the binder chain letter is read off upstream's own span stamp")
    check("_monomer" not in str(data["cif_path"]),
          "the free-binder _monomer prediction is never the design's structure")


def t_multi_target_collect(tmp: Path) -> None:
    out_dir = tmp / "table1_multi" / "bindcraft2"
    root = out_dir / "campaigns"
    design = "pair_bc2_l80_cd34_seq0"
    # Upstream order: descending weight, then name. hEGFR 1.0, mEGFR 0.8, hERBB2 -0.5.
    make_campaign(
        root,
        "pair_bc2",
        [
            {
                "design": design,
                "Binder_Sequence": "AAAA",
                "targets": "hEGFR;mEGFR;hERBB2",
                "target_weights": "1;0.8;-0.5",
                "i_pTM": "0.82;0.79;0.21",
                "i_pAE": "8.1;;13.4",  # mEGFR not measured
                "Interface_Residues": "14;12;3",
                "Binder_Length": 80,  # shared: no `;`
                "failed_filters": "i_pTM.mEGFR,pLDDT.hEGFR",  # comma-joined, not split
                "Timing": "worker=0;start=1;reprediction=2",  # ';' but NOT per-target
                "hash": "cd34",
            }
        ],
        [f"{design}_hEGFR", f"{design}_hERBB2", f"{design}_monomer"],
        spans="C1-20",
    )
    rows = collect_rows(out_dir, ["pair_bc2"], metrics="Timing")
    check(len(rows) == 1, "a multi-target campaign is still ONE row per design")
    row = next(iter(rows.values()))
    data = row["data"]

    check(
        data.get("i_pTM_hEGFR") == 0.82 and data.get("i_pTM_mEGFR") == 0.79,
        "each ';' value becomes a column named for ITS target",
    )
    check(
        data.get("i_pTM_off_hERBB2") == 0.21 and "i_pTM_hERBB2" not in data,
        "a detarget's reading is kept apart under an off_ infix, never merged",
    )
    check("i_pTM" not in data, "the joined ';' cell itself is not left in the table")
    check(
        pd.isna(data.get("i_pAE_mEGFR")) and data.get("i_pAE_hEGFR") == 8.1,
        "an EMPTY position is NA, not 0",
    )
    check(data.get("Binder_Length") == 80,
          "a shared (non per-target) metric stays a single column")
    check(
        data.get("failed_filters") == "i_pTM.mEGFR,pLDDT.hEGFR",
        "failed_filters is comma-joined and is never split by target",
    )
    check(
        data.get("Timing") == "worker=0;start=1;reprediction=2",
        "Timing uses ';' for its own reasons and is never split by target",
    )
    check(
        data.get("targets") == "hEGFR;mEGFR;hERBB2"
        and data.get("target_weights") == "1;0.8;-0.5"
        and data.get("n_targets") == 3,
        "the split key is recorded in the table",
    )
    check(
        data.get("on_targets") == "hEGFR,mEGFR"
        and data.get("off_targets") == "hERBB2",
        "on- and off-targets are named explicitly",
    )
    check(
        str(data.get("path_hEGFR", "")).endswith(".pdb")
        and str(data.get("path_off_hERBB2", "")).endswith(".pdb"),
        "every target with a complex gets its own path column",
    )
    check(
        pd.isna(data.get("path_mEGFR")) and data.get("n_complexes") == 2,
        "a MISSING per-target complex is NA -- never a substituted sibling file",
    )
    check(
        Path(str(row["path"])).stem.startswith(f"{design}_hEGFR"),
        "`path` is the HIGHEST-WEIGHT on-target complex, not the alphabetically first",
    )
    check(
        data.get("binder_chain") == "C",
        "the binder chain letter follows the target's chain count (C behind 2 chains)",
    )
    check(row["status"] == "OK", "the row collects OK")


def t_count_mismatch_is_an_error(tmp: Path) -> None:
    out_dir = tmp / "table1_bad" / "bindcraft2"
    root = out_dir / "campaigns"
    design = "pair_bc2_l80_ef56_seq0"
    make_campaign(
        root,
        "pair_bc2",
        [
            {
                "design": design,
                "targets": "hEGFR;mEGFR;hERBB2",
                "target_weights": "1;0.8;-0.5",
                "i_pTM": "0.82;0.79",  # two values, three targets
            }
        ],
        [f"{design}_hEGFR"],
    )
    row = next(iter(collect_rows(out_dir, ["pair_bc2"]).values()))
    check(
        str(row["status"]).startswith("error:") and "i_pTM" in str(row["status"]),
        "a ';' count that disagrees with `targets` is an ERROR status, not a parse",
    )
    check(
        "i_pTM_hEGFR" not in row["data"] and row["path"] == "",
        "nothing of that design's per-target split is written, and no path",
    )

    out_dir2 = tmp / "table1_bad2" / "bindcraft2"
    make_campaign(
        out_dir2 / "campaigns",
        "pair_bc2",
        [{"design": design, "targets": "hEGFR;mEGFR", "target_weights": "1"}],
        [f"{design}_hEGFR"],
    )
    row = next(iter(collect_rows(out_dir2, ["pair_bc2"]).values()))
    check(
        "target_weights" in str(row["status"]) and str(row["status"]).startswith("error"),
        "targets/target_weights of different lengths is an error, not a guess",
    )


def t_missing_primary_complex(tmp: Path) -> None:
    out_dir = tmp / "table1_nop" / "bindcraft2"
    design = "pair_bc2_l80_gh78_seq0"
    make_campaign(
        out_dir / "campaigns",
        "pair_bc2",
        [
            {
                "design": design,
                "targets": "hEGFR;hERBB2",
                "target_weights": "1;-0.5",
                "i_pTM": "0.82;0.21",
            }
        ],
        [f"{design}_hERBB2"],  # only the OFF-target complex exists
    )
    row = next(iter(collect_rows(out_dir, ["pair_bc2"]).values()))
    check(
        "on-target" in str(row["status"]) and str(row["status"]).startswith("error"),
        "a detarget complex is NEVER promoted to the design's `path`",
    )
    check(row["path"] == "", "and the row carries no structure path at all")


def t_no_prefix_collision(tmp: Path) -> None:
    """`<design>_seq1*` also matches `<design>_seq10`: never pick one of several."""
    out_dir = tmp / "table1_pfx" / "bindcraft2"
    ranked = out_dir / "campaigns" / "g" / "3_Ranked"
    ranked.mkdir(parents=True)
    pd.DataFrame(
        [{"design": "g_seq1", "i_pTM": 0.5}, {"design": "g_seq10", "i_pTM": 0.6}]
    ).to_csv(ranked / "!_Ranked.csv", index=False)
    write_cif(ranked / "g_seq10.cif")  # g_seq1's own file is absent
    rows = collect_rows(out_dir, ["g"])
    check(
        "g_seq1" not in rows and "g_seq10" in rows,
        "a design with no file of its own is skipped, not given a sibling's structure",
    )


def main() -> None:
    with tempfile.TemporaryDirectory() as raw:
        tmp = Path(raw)
        for title, fn in (
            ("--targets schema validation", lambda: t_schema_validation(tmp)),
            ("--targets flag exclusivity", t_exclusivity),
            ("target_path resolution and {expr}", lambda: t_target_path_and_expr(tmp)),
            ("settings file regression", lambda: t_settings_regression(tmp)),
            ("manifest end to end", lambda: t_manifest_end_to_end(tmp)),
            ("single-target collect (regression)", lambda: t_single_target_collect(tmp)),
            ("multi-target collect", lambda: t_multi_target_collect(tmp)),
            ("error contract: ';' count mismatch", lambda: t_count_mismatch_is_an_error(tmp)),
            ("error contract: no on-target complex", lambda: t_missing_primary_complex(tmp)),
            ("no design-name prefix collision", lambda: t_no_prefix_collision(tmp)),
        ):
            print(f"\n{title}")
            fn()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILED:")
        for message in FAILURES:
            print(f"  - {message}")
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
