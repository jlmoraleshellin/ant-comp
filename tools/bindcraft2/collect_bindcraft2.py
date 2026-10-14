#!/usr/bin/env python3
"""
Collect BindCraft2 campaign results into the (child) binder table.

``run_bindcraft2.py`` points each campaign's ``project_folder`` at
``<out_dir>/campaigns/<name>/``, where BindCraft2 writes three stage folders, each
with a table this collector can read:

    campaigns/<name>/1_Trajectories/!_Trajectories.csv        one row per attempt
    campaigns/<name>/1_Trajectories/<design>/<design>_trajectory[_<target>].cif
    campaigns/<name>/2_Refolded/!_Refolded.csv       every scored candidate + outcome
    campaigns/<name>/2_Refolded/Complexes/<design>[_<target>].cif
    campaigns/<name>/3_Ranked/!_Ranked.csv          the accepted designs, best-first
    campaigns/<name>/3_Ranked/<design>[_<target>].cif

Which one to read is **decided by the run, not guessed here**: a `--trajectory-only`
run fills `1_Trajectories` and never accepts anything, so `3_Ranked` stays empty. The
run records that in the out_dir sidecar and `--stage auto` (the default) follows it,
the same contract that keeps `input_column` consistent between the two phases. Pass
`--stage` explicitly only to override.

* ``trajectories`` -- the hallucinated **backbones**, before any sequence redesign.
  The rows to hand to `atomium` / `proteinmpnn`. Note their `sequence` is the one the
  gradient happened to land on, which is exactly what you are about to replace.
* ``refolded`` -- every candidate the campaign scored, rejects included, with
  ``outcome`` and ``failed_filters``. What to read when a campaign accepted nothing.
* ``ranked`` -- only what the campaign accepted.

A trajectory row's ``terminated`` column names the stage the attempt stopped at;
blank means it ran to completion. They are collected rather than dropped, so filter
on it (``-f``) instead of assuming every row is a usable backbone.

MULTI-TARGET
------------
A campaign with TWO OR MORE targets does not write one column per target. BindCraft2
packs each metric into a single cell, semicolon-separated, in the order given by the
row's own ``targets`` column -- which is sorted by weight, so off-targets (negative
weight) come last (``campaign_output.target_ordered_row``; a campaign with fewer than
two targets is left completely untouched, so single-target tables are unaffected).

A packed cell is not a number: ``-f`` cannot threshold it and a sort would order it
as text. So this collector keeps the packed cell verbatim **and** splits it:

    bindcraft2_targets          hPDL1;mPDL1;hPD1
    bindcraft2_target_weights   1;1;-0.5
    bindcraft2_i_pTM            0.82;0.79;0.21     <- verbatim, as BC2 wrote it
    bindcraft2_i_pTM__hPDL1     0.82               <- filterable
    bindcraft2_i_pTM__mPDL1     0.79
    bindcraft2_i_pTM__hPD1      0.21
    bindcraft2_i_pTM__worst     0.79               <- binding targets only
    bindcraft2_i_pTM__selectivity  0.58            <- weakest binder vs best off-target

The summary columns cover the **binding** targets only, mirroring
``campaign_output.on_target_mean``: averaging an off-target in would reward binding
the thing you are trying to avoid. They are emitted only for metrics whose direction
is known (METRIC_DIRECTION below); a metric outside that map still gets its
per-target split, just no summary.

Structures are per-target too (``<design>_<target>.cif``). ``<leaf>_path`` is the
highest-weight binding target's complex -- the one a downstream predictor should
look at -- and every other target's complex is kept beside it as
``<leaf>_path__<target>``. The earlier version of this collector globbed and took
the first match, which silently picked one arbitrary target and dropped the rest.

Rows are keyed by BindCraft2's own ``design`` identity rather than by position:
``!_Ranked.csv`` is re-sorted by ``i_pDAE`` after every acceptance, so an index into
it is not stable across a re-collect, while the design name (campaign, modality,
length and recipe hash) is. The submitter sets ``campaign_name`` to the parent row's
name, so a design name already starts with its parent; the prefix is only added here
when a --set override has changed that.

Each complex is converted to PDB (downstream tools consume PDB) and becomes the row's
``<leaf>_path``; the source mmCIF stays available as ``<leaf>_cif_path``. Note the
structure is the **complex**, binder + target, not the binder alone -- unless the
campaign ran with ``--save-monomers``, in which case the free binder is recorded
separately as ``<leaf>_monomer_path``.

If ``archive_trajectories`` zipped the per-design folders, the structures are read
straight out of the zip (stdlib ``zipfile``; the archive is a plain deflate zip whose
members are relative to ``1_Trajectories``) into ``<out_dir>/.unarchived/``. Nothing
in the campaign folder is modified.

Safe to re-run: rows are rebuilt from the CSV and the structures on disk.

Usage:
    sapia collect bindcraft2 outputs/RUN --table <the table the run reserved>
    sapia collect bindcraft2 outputs/RUN --table <table> --stage refolded --force
"""

import json
import zipfile
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, Callable, Iterable, NamedTuple

import pandas as pd

from prosapia.core import (
    RUN_META_FILENAME,
    CollectArgs,
    CollectCtx,
    CollectEach,
    Collected,
    DesignCtx,
)
from prosapia.utils import ensure_pdb

# Layout run_bindcraft2.py imposes on the out_dir -- keep the two in step.
CAMPAIGNS_DIRNAME = "campaigns"

# --stage default: read what the run recorded in the sidecar instead of guessing.
STAGE_AUTO = "auto"

# BindCraft2's multi-target packing (campaign_output.py). The separator, and the two
# columns that say what the packed positions mean.
TARGET_VALUE_SEPARATOR = ";"
TARGET_NAME_COLUMN = "targets"
TARGET_WEIGHT_COLUMN = "target_weights"

# Suffix joining a metric (or a path) to the target it was measured against. Two
# underscores so it cannot be confused with BindCraft2's own single-underscore names
# (`i_pTM_detarget` is one metric; `i_pTM__hPD1` is one metric on one target).
TARGET_COLUMN_SEP = "__"

# Where unarchived trajectory structures are staged, under the collect out_dir. The
# campaign folder itself is never written to.
UNARCHIVED_DIRNAME = ".unarchived"


def _flat_structure(search_dir: Path, design: str) -> Path | None:
    """``<search_dir>/<design>.cif``, the un-suffixed single-target form."""
    exact = search_dir / f"{design}.cif"
    return exact if exact.is_file() else None


def _target_structure(search_dir: Path, design: str, target: str) -> Path | None:
    """``<search_dir>/<design>_<target>.cif``, as a multi-target campaign writes."""
    path = search_dir / f"{design}_{target}.cif"
    return path if path.is_file() else None


def _monomer_structure(search_dir: Path, design: str) -> Path | None:
    """The free binder, re-predicted alone (``save_binder_monomers``)."""
    path = search_dir / f"{design}_monomer.cif"
    return path if path.is_file() else None


def _any_structure(search_dir: Path, design: str) -> Path | None:
    """Last resort: any complex for this design, monomers excluded.

    Only reached when neither the plain nor the per-target name matched, which means
    BindCraft2 named the file in a way this collector does not know about. Returning
    something beats dropping the design, but the caller says so in the row's status.
    """
    matches = sorted(
        p for p in search_dir.glob(f"{design}*.cif") if "_monomer" not in p.name
    )
    return matches[0] if matches else None


class Stage(NamedTuple):
    """Where a stage keeps its table and its structures.

    ``search`` is the folder holding a design's structures and ``base`` the filename
    stem upstream builds them from; a multi-target complex is exactly
    ``<search>/<base>_<target name>.cif`` (``accepted_state_suffixes``), which is why
    the two are kept apart rather than baked into one finder.
    """

    folder: str  # under project_folder
    table: str  # its CSV, inside that folder
    # (stage_dir, design) -> the dir holding that design's structures, or None
    search: Callable[[Path, str], Path | None]
    # the structure filename stem, which for a trajectory is not just the design
    stem: Callable[[str], str]


def _ranked_dir(stage_dir: Path, design: str) -> Path | None:
    return stage_dir if stage_dir.is_dir() else None


def _refolded_dir(stage_dir: Path, design: str) -> Path | None:
    complexes = stage_dir / "Complexes"
    return complexes if complexes.is_dir() else None


def _trajectory_dir(stage_dir: Path, design: str) -> Path | None:
    """``1_Trajectories/<design>/`` -- a per-design subfolder, unlike the two later
    stages. When ``archive_trajectories`` packed it into ``<design>.zip`` the folder
    is gone and the structures live inside the archive; see _unarchived_dir."""
    design_dir = stage_dir / design
    return design_dir if design_dir.is_dir() else None


STAGES: dict[str, Stage] = {
    "trajectories": Stage(
        "1_Trajectories",
        "!_Trajectories.csv",
        _trajectory_dir,
        lambda design: f"{design}_trajectory",
    ),
    "refolded": Stage(
        "2_Refolded", "!_Refolded.csv", _refolded_dir, lambda design: design
    ),
    "ranked": Stage("3_Ranked", "!_Ranked.csv", _ranked_dir, lambda design: design),
}

# The column carrying the design's identity, and the one carrying its sequence.
DESIGN_COL = "design"
SEQUENCE_COL = "Binder_Sequence"
# 1_Trajectories only: the stage the attempt stopped at, blank when it completed.
TERMINATED_COL = "terminated"

# Column prefixes always collected whatever --metrics says: user metadata supplied
# with `--metadata` (meta_*), and the record of which presets actually ran
# (settings_core / settings_modality / settings_property / settings_target /
# settings_overrides). Provenance is cheap and a campaign cannot be re-read without
# it.
ALWAYS_PREFIXES = ("meta_", "settings_")

# Collected by default: the metrics a binder campaign is actually read on. Everything
# else in the CSV is reachable with --metrics or --all-metrics rather than being
# poured into the table by default -- BindCraft2 writes ~60 columns per design.
CORE_METRICS = (
    # confidence
    "i_pDAE",
    "i_pTM",
    "i_pAE",
    "pLDDT",
    "pTM",
    "Unbound_Binder_pLDDT",
    "Target_pLDDT",
    # interface
    "Interface_Residues",
    "Interface_BuriedArea",
    "Hotspot_Contact_Fraction",
    "Coldspot_Contact_Fraction",
    "Interface_Binder_Residues",
    "Interface_Target_Residues",
    # detargeting: only present when the campaign carried an off-target, and the
    # whole point of one. Interface confidence alone does not decide avoidance --
    # a peptide can read 0.27 i_pTM with its whole face on the off-target -- so the
    # residue count is collected beside it.
    "i_pTM_detarget",
    "i_pAE_detarget",
    "Interface_Residues_detarget",
    # geometry sanity
    "Backbone_Clashes",
    "All_Atom_Clashes",
    "Binder_Chain_Breaks",
    # developability
    "Surface_Hydrophobicity",
    "Binder_Length",
    "Binder_Net_Charge",
    "Binder_Free_Cysteines",
    # multi-target bookkeeping: what the packed cells mean, and in what order.
    TARGET_NAME_COLUMN,
    TARGET_WEIGHT_COLUMN,
    # provenance / triage. Which of these exist depends on the stage: `outcome` and
    # `failed_filters` on refolded, `rank` on ranked, `terminated` and `autotuned`
    # on trajectories. Absent ones are simply not collected.
    "rank",
    "hash",
    "trajectory",
    "outcome",
    "failed_filters",
    "terminated",
    "autotuned",
    "length",
)

# Direction of the metrics this collector will summarise across targets: True when a
# larger reading is the better one. A metric outside this map still gets its
# per-target split; it just gets no __mean/__worst/__best/__selectivity, because
# those three words are meaningless without a direction. Read off
# docs/source/outputs.md (BindCraft2 v1.0.3).
METRIC_DIRECTION: dict[str, bool] = {
    # higher is better
    "i_pDAE": True,
    "i_pTM": True,
    "pTM": True,
    "pLDDT": True,
    "Unbound_Binder_pLDDT": True,
    "Target_pLDDT": True,
    "SS_pLDDT": True,
    "Interface_Residues": True,
    "Interface_BuriedArea": True,
    "Hotspot_Contact_Fraction": True,
    "Epitope_Residues_Contacted": True,
    "Receptor_Chains_Contacted": True,
    "Framework_Packing_Fraction": True,
    "Domain_Separation_Ratio": True,
    "Scaffold_Sequence_Retained_Fraction": True,
    # lower is better
    "i_pAE": False,
    "Coldspot_Contact_Fraction": False,
    "Off_Paratope_Contact_Fraction": False,
    "Off_Epitope_Contact_Fraction": False,
    "Surface_Hydrophobicity": False,
    "Backbone_Clashes": False,
    "All_Atom_Clashes": False,
    "Binder_Chain_Breaks": False,
    "Binder_Free_Cysteines": False,
    "MHC_Anchor_Score": False,
    "Protease_Site_Score": False,
    "Exposed_Loop_Fraction": False,
    "Terminus_Exposure": False,
    "Interdomain_Contact_Fraction": False,
    "Scaffold_Framework_RMSD": False,
    "Oligomer_Symmetry_RMSD": False,
}


class BindCraft2CollectArgs(CollectArgs):
    stage: str
    metrics: str
    all_metrics: bool
    no_split_targets: bool


def add_collect_bindcraft2_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--stage",
        choices=(STAGE_AUTO, *STAGES),
        default=STAGE_AUTO,
        help="Which campaign stage to collect. 'auto' (the default) reads what the "
        "run recorded: 'trajectories' after a --trajectory-only run, 'ranked' "
        "otherwise. 'trajectories' takes the hallucinated backbones before any "
        "sequence redesign; 'refolded' takes every scored candidate, rejects "
        "included, with its `outcome` and `failed_filters`; 'ranked' takes only what "
        "the campaign accepted. NOTE -l/--dir-label does NOT separate two stages into "
        "separate columns: it only selects a leaf dir some run already created, and "
        "both stages write the same <leaf>_* column names. Collecting a second stage "
        "into the SAME leaf ADDS ROWS (stage row keys differ), it does not add "
        "columns. For two stages side by side, give the second one its own leaf via a "
        "--reuse-campaigns run with a different -l.",
    )
    parser.add_argument(
        "--metrics",
        type=str,
        default="",
        help="Comma-separated extra CSV columns to collect on top of the core set "
        "(e.g. 'Binder_pI,Binder_Helix_Fraction,Interface_W_Count'). Names are "
        "BindCraft2's own, as spelled in !_Ranked.csv. On a multi-target campaign a "
        "per-target column still splits into one column per target.",
    )
    parser.add_argument(
        "--all-metrics",
        action="store_true",
        help="Collect every column in the stage's CSV instead of the core set. "
        "Wide: BindCraft2 writes ~60 measurements per design, and a multi-target "
        "campaign then splits each of them per target as well.",
    )
    parser.add_argument(
        "--no-split-targets",
        action="store_true",
        help="Keep BindCraft2's semicolon-packed multi-target cells as they are, "
        "without adding the per-target and summary columns. The packed cell is a "
        "string, so -f cannot threshold it -- only pass this if the extra columns "
        "are in the way. No effect on a single-target campaign, which BindCraft2 "
        "never packs.",
    )


def _run_meta(out_dir: Path) -> dict:
    """The run's sidecar, or ``{}`` when absent or unreadable."""
    path = out_dir / RUN_META_FILENAME
    try:
        return json.loads(path.read_text()) if path.is_file() else {}
    except (OSError, json.JSONDecodeError):
        return {}


def resolve_stage(ctx: CollectCtx[BindCraft2CollectArgs], meta: dict) -> str:
    """The stage to collect: the flag, or what the run recorded in the sidecar.

    A run older than the sidecar key (or one whose settings set ``trajectory_only``
    through the verbatim ``--set``) reads as a full campaign, which is the safe
    default -- it looks in 3_Ranked and reports nothing rather than inventing rows.
    """
    if ctx.args.stage != STAGE_AUTO:
        return ctx.args.stage
    return "trajectories" if meta.get("trajectory_only") else "ranked"


def resolve_campaigns_root(ctx: CollectCtx[BindCraft2CollectArgs], meta: dict) -> Path:
    """Where the campaigns actually live.

    Normally this out_dir, but a ``--reuse-campaigns`` run reserves a table over an
    EARLIER run's campaigns so one campaign can be collected into two tables (its
    backbones into one, its accepted designs into another). The run records the path;
    collect follows it rather than assuming the two coincide.
    """
    recorded = meta.get("campaigns_root")
    return Path(recorded) if recorded else ctx.out_dir / CAMPAIGNS_DIRNAME


def _read_stage_table(campaign_dir: Path, stage: str) -> tuple[pd.DataFrame, Path]:
    """The stage's CSV as a frame, plus the stage dir its structures hang off.

    An absent stage folder or CSV yields an empty frame: a campaign that accepted
    nothing (or was killed before its first acceptance) is a normal outcome, not an
    error -- `2_Refolded` then still says why.

    Left on pandas' own dtype inference deliberately, so a single-target campaign
    collects byte-for-byte as it always has. A per-target column is a ``;``-joined
    string in *every* row (``target_ordered_row`` joins over every target name), so it
    arrives as text without being asked to.
    """
    folder, csv_name = STAGES[stage].folder, STAGES[stage].table
    stage_dir = campaign_dir / folder
    csv_path = stage_dir / csv_name
    if not csv_path.is_file():
        return pd.DataFrame(), stage_dir
    try:
        return pd.read_csv(csv_path), stage_dir
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
        print(f"  unreadable {csv_path}: {e}")
        return pd.DataFrame(), stage_dir


def _packed(value: Any) -> list[str]:
    """Split one of BindCraft2's packed multi-target cells into its positions."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    return str(value).split(TARGET_VALUE_SEPARATOR)


def row_targets(row: pd.Series) -> list[tuple[str, float]]:
    """``[(target name, weight), …]`` for a row, in BindCraft2's own packed order.

    Empty for a single-target campaign: BindCraft2 writes neither column then
    (``target_ordered_row`` returns the row untouched below two targets), and there
    is nothing to split. A row whose names and weights disagree in length is treated
    as unusable rather than zipped into a wrong mapping.
    """
    if TARGET_NAME_COLUMN not in row.index:
        return []
    names = [name for name in _packed(row.get(TARGET_NAME_COLUMN)) if name]
    weights = _packed(row.get(TARGET_WEIGHT_COLUMN))
    if not names:
        return []
    if len(weights) != len(names):
        # Weights missing or mismatched: keep the order, assume every target binds.
        return [(name, 1.0) for name in names]
    paired = []
    for name, weight in zip(names, weights):
        try:
            paired.append((name, float(weight)))
        except ValueError:
            paired.append((name, 1.0))
    return paired


def _number(value: str) -> float | str:
    """A packed position as a number when it is one, else verbatim (it may be a
    residue list, an empty position, or a filter name)."""
    text = value.strip()
    if not text:
        return ""
    try:
        return float(text)
    except ValueError:
        return text


def _summaries(
    metric: str, readings: dict[str, float], weights: dict[str, float]
) -> dict[str, float]:
    """Across-target summary columns for one metric, binding targets only.

    Mirrors ``campaign_output.on_target_mean``: an off-target reading is gated at
    acceptance, and averaging it in would reward binding the thing the campaign is
    trying to avoid. ``__selectivity`` is the one column that deliberately compares
    the two groups -- the weakest binding target against the strongest off-target, in
    the metric's own direction, so a positive margin always favours the intended
    distinction.
    """
    higher = METRIC_DIRECTION.get(metric)
    if higher is None:
        return {}
    binding = [v for name, v in readings.items() if weights.get(name, 1.0) > 0]
    against = [v for name, v in readings.items() if weights.get(name, 1.0) < 0]
    if len(binding) < 2 and not against:
        return {}
    out: dict[str, float] = {}
    if binding:
        out[f"{metric}{TARGET_COLUMN_SEP}mean"] = sum(binding) / len(binding)
        out[f"{metric}{TARGET_COLUMN_SEP}worst"] = min(binding) if higher else max(binding)
        out[f"{metric}{TARGET_COLUMN_SEP}best"] = max(binding) if higher else min(binding)
        out[f"{metric}{TARGET_COLUMN_SEP}spread"] = max(binding) - min(binding)
    if binding and against:
        out[f"{metric}{TARGET_COLUMN_SEP}selectivity"] = (
            min(binding) - max(against) if higher else min(against) - max(binding)
        )
    return out


def split_target_columns(
    data: dict[str, Any], targets: list[tuple[str, float]]
) -> dict[str, Any]:
    """Per-target and summary columns derived from the packed cells in ``data``.

    A cell is only split when it holds exactly as many positions as the row has
    targets -- the same guard ``campaign_output.per_target_readings`` applies, and
    the thing that stops a residue list that happens to contain a semicolon from
    being mapped onto the wrong targets.
    """
    names = [name for name, _ in targets]
    weights = dict(targets)
    derived: dict[str, Any] = {}
    for column, value in data.items():
        if column in (TARGET_NAME_COLUMN, TARGET_WEIGHT_COLUMN):
            continue
        positions = _packed(value)
        if len(positions) != len(names) or len(names) < 2:
            continue
        readings: dict[str, float] = {}
        for name, position in zip(names, positions):
            parsed = _number(position)
            derived[f"{column}{TARGET_COLUMN_SEP}{name}"] = parsed
            if isinstance(parsed, float):
                readings[name] = parsed
        if len(readings) == len(names):
            derived.update(_summaries(column, readings, weights))
    return derived


def _safe_members(archive: zipfile.ZipFile, design: str) -> list[str]:
    """Members of a trajectory archive belonging to ``design``, path-traversal safe.

    ``archive_trajectory_folder`` writes members relative to ``1_Trajectories``, so
    they all start ``<design>/``. Anything absolute, escaping upward, or outside that
    prefix is not ours and is skipped.
    """
    prefix = f"{design}/"
    members = []
    for member in archive.namelist():
        if member.startswith(("/", "\\")) or ".." in Path(member).parts:
            continue
        if member.startswith(prefix):
            members.append(member)
    return members


def unarchived_dir(stage_dir: Path, design: str, dest_root: Path) -> Path | None:
    """Restore one archived trajectory folder into ``dest_root`` and return it.

    ``archive_trajectories`` zips each ``1_Trajectories/<design>/`` into
    ``<design>.zip`` and deletes the folder, which used to leave this collector
    reporting a campaign's worth of designs as having no structure. The archive is a
    plain deflate zip, so it is read with the stdlib here rather than shelling out to
    `bindcraft unarchive` -- which is not installed in the workstation, and would
    rewrite someone's campaign folder. Nothing in the campaign is modified.
    """
    archive_path = stage_dir / f"{design}.zip"
    if not archive_path.is_file():
        return None
    destination = dest_root / design
    if destination.is_dir():
        return destination
    try:
        with zipfile.ZipFile(archive_path) as archive:
            members = _safe_members(archive, design)
            if not members:
                return None
            dest_root.mkdir(parents=True, exist_ok=True)
            archive.extractall(dest_root, members=members)
    except (OSError, zipfile.BadZipFile) as e:
        print(f"  unreadable archive {archive_path}: {e}")
        return None
    return destination if destination.is_dir() else None


class Structures(NamedTuple):
    """The structure files one design produced at one stage."""

    primary: Path | None  # the complex `<leaf>_path` points at
    per_target: dict[str, Path]  # target name -> its own complex
    monomer: Path | None  # the free binder, when save_binder_monomers was on
    fell_back: bool  # primary found by glob rather than by name


def find_structures(
    search_dir: Path, stem: str, targets: list[tuple[str, float]]
) -> Structures:
    """Locate every structure this design wrote, by name rather than by glob order.

    Single target: ``<stem>.cif``. Several: one ``<stem>_<target>.cif`` per target,
    and the primary is the highest-weight BINDING target -- the complex a downstream
    predictor or interface tool should be pointed at. ``targets`` arrives already
    ordered by ``(-weight, name)``, so the first binding entry is that one.
    """
    per_target: dict[str, Path] = {}
    for name, _ in targets:
        found = _target_structure(search_dir, stem, name)
        if found is not None:
            per_target[name] = found

    primary: Path | None = None
    for name, weight in targets:
        if weight > 0 and name in per_target:
            primary = per_target[name]
            break
    if primary is None and per_target:
        primary = next(iter(per_target.values()))

    fell_back = False
    if primary is None:
        primary = _flat_structure(search_dir, stem)
    if primary is None:
        primary = _any_structure(search_dir, stem)
        fell_back = primary is not None

    return Structures(
        primary=primary,
        per_target=per_target,
        monomer=_monomer_structure(search_dir, stem),
        fell_back=fell_back,
    )


def _row_data(row: pd.Series, columns: list[str]) -> dict[str, Any]:
    """The tool-specific columns for one design, as bare names (the driver prefixes)."""
    data: dict[str, Any] = {}
    if SEQUENCE_COL in row.index:
        # Named `sequence` like every other sequence-producing tool here, so a
        # downstream predictor takes `-i bindcraft2_sequence` and nothing else. One
        # binder sequence per design, multi-target or not -- that is the whole point.
        data["sequence"] = row[SEQUENCE_COL]
    for col in columns:
        if col in row.index:
            data[col] = row[col]
    return data


def _selected_columns(df: pd.DataFrame, ctx: CollectCtx[BindCraft2CollectArgs]) -> list[str]:
    """Which CSV columns to carry into the table."""
    extra = [m.strip() for m in ctx.args.metrics.split(",") if m.strip()]
    if ctx.args.all_metrics:
        chosen = [c for c in df.columns if c not in (DESIGN_COL, SEQUENCE_COL)]
    else:
        chosen = [c for c in (*CORE_METRICS, *extra) if c in df.columns]
    always = [
        c
        for c in df.columns
        if c.startswith(ALWAYS_PREFIXES) and c not in chosen
    ]
    return chosen + always


def collect_bindcraft2(ctx: CollectCtx[BindCraft2CollectArgs]) -> CollectEach:
    """Per-campaign BindCraft2 collector. A create tool: the framework iterates the
    ready parents (table rows, or the root design groups the run recorded) and this
    mints one child row per design the campaign produced, carrying ``parent`` for
    lineage. A campaign with no results yields no rows. The framework stamps
    status/path/parent_name from each Collected."""
    meta = _run_meta(ctx.out_dir)
    if meta.get("weights_only"):
        print(
            "This run only fetched the AlphaFold parameters into the cache volume; "
            "it ran no campaign, so there is nothing to collect."
        )
        return lambda d: ()

    stage = resolve_stage(ctx, meta)
    stage_spec = STAGES[stage]
    campaigns_root = resolve_campaigns_root(ctx, meta)
    run_dir = ctx.args.run_dir
    unarchive_root = ctx.out_dir / UNARCHIVED_DIRNAME
    how = "from --stage" if ctx.args.stage != STAGE_AUTO else "per the run's sidecar"
    print(f"Collecting '{stage_name}' designs ({how}) from campaigns in {campaigns_root}")
    if submitted:
        print(f"  run submitted {len(submitted)} target(s): {', '.join(submitted)}")
    warned_mismatch: set[str] = set()

    def as_pdb(cif: Path) -> tuple[Path, str]:
        """(structure for the table, status). A design whose mmCIF won't parse still
        has its metrics, and one bad file among hundreds shouldn't abort the collect."""
        try:
            return ensure_pdb(cif, run_dir), "OK"
        except (ValueError, RuntimeError, OSError) as e:
            return cif, f"error: unreadable mmCIF ({type(e).__name__})"

    def one(d: DesignCtx) -> Iterable[Collected]:
        campaign_dir = campaigns_root / d.name
        if not campaign_dir.is_dir():
            print(f"{d.name}: no campaign dir, skipping")
            return

        df, stage_dir = _read_stage_table(campaign_dir, stage_name)
        if df.empty or DESIGN_COL not in df.columns:
            print(f"{d.name}: no '{stage_name}' designs")
            return

        columns = _selected_columns(df, ctx)

        # Sorted by identity, not by rank: !_Ranked.csv is re-sorted after every
        # acceptance, so collecting in file order would be collecting in an order
        # that changes under a resumed campaign.
        n = 0
        restored = 0
        ordered = df.sort_values(DESIGN_COL).set_index(DESIGN_COL)
        for design, row in ordered.iterrows():
            design = str(design)
            search_dir = stage_spec.search(stage_dir, design)
            if search_dir is None and stage == "trajectories":
                # archive_trajectories packed the folder away; read it out of the zip
                # instead of reporting the design as structureless.
                search_dir = unarchived_dir(
                    stage_dir, design, unarchive_root / d.name
                )
                restored += search_dir is not None
            if search_dir is None:
                print(f"{d.name}: {design} has no structure dir under {stage_dir}")
                continue

            targets = row_targets(row)
            stem = stage_spec.stem(design)
            found = find_structures(search_dir, stem, targets)
            if found.primary is None:
                print(f"{d.name}: {design} has no structure in {search_dir}, skipping")
                continue

            data = _row_data(row, columns)
            if targets and not ctx.args.no_split_targets:
                data.update(split_target_columns(data, targets))
            data["cif_path"] = str(found.primary)
            if found.monomer is not None:
                data["monomer_path"] = str(found.monomer)
            if stage == "trajectories":
                # `terminated` is blank when the attempt ran to completion, so its
                # NA means success -- the opposite of NA everywhere else in a
                # prosapia table. Carry the polarity explicitly so a --filter reads
                # `bindcraft2_completed == True` instead of testing for a blank.
                stopped_at = row.get(TERMINATED_COL)
                data["completed"] = bool(
                    pd.isna(stopped_at) or not str(stopped_at).strip()
                )

            # A design whose mmCIF won't parse still has its metrics, and one bad
            # file among hundreds shouldn't abort the collect: keep the row, point
            # it at the mmCIF, and say so in its status.
            try:
                path, status = ensure_pdb(found.primary, run_dir), "OK"
            except (ValueError, RuntimeError, OSError) as e:
                print(f"{d.name}: {design} mmCIF unreadable ({e})")
                path, status = found.primary, f"error: unreadable mmCIF ({type(e).__name__})"

            # Every other target's complex, beside the primary one. Converted too:
            # a multi-target campaign is read by comparing the poses against each
            # other, and a half-converted set makes that a two-step job.
            for target, cif in found.per_target.items():
                if cif == found.primary:
                    continue
                try:
                    data[f"path{TARGET_COLUMN_SEP}{target}"] = str(
                        ensure_pdb(cif, run_dir)
                    )
                except (ValueError, RuntimeError, OSError):
                    data[f"path{TARGET_COLUMN_SEP}{target}"] = str(cif)

            if found.fell_back and status == "OK":
                # The name did not match either convention, so which target this
                # structure belongs to is a guess. Say so in the row rather than
                # letting it pass as a clean result.
                status = "OK: structure matched by glob, not by target name"

            yield Collected(
                name=row_name, parent=d.name, path=path, status=status, data=data
            )
            n += 1
        note = f", {restored} unarchived" if restored else ""
        print(f"{d.name}: OK ({n} {stage} design(s){note})")

    return one


def _single_target_row(
    row: pd.Series,
    columns: list[str],
    search_dir: Path,
    base: str,
    as_pdb: Callable[[Path], tuple[Path, str]],
) -> tuple[dict[str, Any], Path, str] | None:
    """One design of a single-target campaign -- the shape this tool has always had.

    Returns ``None`` when there is no structure at all (the design is then skipped, as
    before).
    """
    cif = _unambiguous_structure(search_dir, base)
    if cif is None:
        return None
    data = _shared_row_data(row, columns)
    data["cif_path"] = str(cif)
    data["n_targets"] = 1
    path, status = as_pdb(cif)
    binder_chain, source = binder_chain_letters(cif, path if status == "OK" else None)
    data["binder_chain"] = binder_chain
    data["binder_chain_src"] = source
    return data, path, status


def _multi_target_row(
    row: pd.Series,
    columns: list[str],
    targets: list[TargetRef],
    search_dir: Path,
    base: str,
    as_pdb: Callable[[Path], tuple[Path, str]],
) -> tuple[dict[str, Any], Path | str, str] | None:
    """One design of a multi-target campaign: one row, one column set per target.

    The design stays a single row -- a binder sequence is one entity. Every ``;`` cell
    becomes ``<metric>_<target>`` (or ``<metric>_off_<target>`` for a detarget), every
    target gets its own ``path_<target>``, and ``path`` points at the highest-weight
    on-target complex so the existing downstream steps keep working unchanged.
    """
    shared = [c for c in columns if c not in (TARGETS_COL, TARGET_WEIGHTS_COL)]
    data = _shared_row_data(row, shared)
    data[TARGETS_COL] = TARGET_SEP.join(t.name for t in targets)
    data[TARGET_WEIGHTS_COL] = TARGET_SEP.join(f"{t.weight:g}" for t in targets)
    data["n_targets"] = len(targets)
    on = [t for t in targets if not t.off]
    data["on_targets"] = ",".join(t.name for t in on) or pd.NA
    data["off_targets"] = ",".join(t.name for t in targets if t.off) or pd.NA

    # Split into a side dict first: a mismatch anywhere makes the WHOLE design an
    # error, and a half-applied split would leave the row looking collected.
    per_target: dict[str, Any] = {}
    try:
        for col in shared:
            if col not in row.index:
                continue
            parts = split_target_cell(col, row[col], len(targets))
            if parts is None:
                continue
            data.pop(col, None)  # a per-target metric has no single shared value
            for target, value in zip(targets, parts):
                per_target[f"{col}_{target.label}"] = value
    except TargetParseError as e:
        return _shared_row_data(row, shared) | {
            TARGETS_COL: data[TARGETS_COL],
            TARGET_WEIGHTS_COL: data[TARGET_WEIGHTS_COL],
            "n_targets": len(targets),
        }, "", f"error: {e}"
    data |= per_target

    # One complex per target, matched by NAME. Never by file order (alphabetical) and
    # never by falling back to a sibling: an absent complex is NA.
    found = 0
    for target in targets:
        cif = search_dir / f"{base}_{target.name}.cif"
        if not cif.is_file():
            data[f"path_{target.label}"] = pd.NA
            continue
        found += 1
        per_target_path, per_target_status = as_pdb(cif)
        data[f"path_{target.label}"] = (
            str(per_target_path) if per_target_status == "OK" else pd.NA
        )
    data["n_complexes"] = found
    if not found:
        return None

    # `path` = the highest-weight ON-target complex. `targets` is already ordered by
    # descending weight then name (campaign_output.weighted_target_order), so the first
    # attracting entry is the one. A detarget complex is never the design's `path`.
    if not on:
        return data, "", "error: every target is a detarget, so there is no on-target complex"
    primary = on[0]
    primary_cif = search_dir / f"{base}_{primary.name}.cif"
    if not primary_cif.is_file():
        missing = (
            f"error: no complex for the highest-weight on-target "
            f"{primary.name!r} at {primary_cif.name}"
        )
        return data, "", missing
    data["cif_path"] = str(primary_cif)
    path, status = as_pdb(primary_cif)
    binder_chain, source = binder_chain_letters(
        primary_cif, path if status == "OK" else None
    )
    data["binder_chain"] = binder_chain
    data["binder_chain_src"] = source
    return data, path, status
