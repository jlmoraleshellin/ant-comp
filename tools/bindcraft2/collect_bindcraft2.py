#!/usr/bin/env python3
"""
Collect BindCraft2 campaign results into the (child) binder table.

``run_bindcraft2.py`` points each campaign's ``project_folder`` at
``<out_dir>/campaigns/<name>/``, where BindCraft2 writes three stage folders, each
with a table this collector can read:

    campaigns/<name>/1_Trajectories/!_Trajectories.csv        one row per attempt
    campaigns/<name>/1_Trajectories/<design>/<design>_trajectory.cif
    campaigns/<name>/2_Refolded/!_Refolded.csv       every scored candidate + outcome
    campaigns/<name>/2_Refolded/Complexes/<design>.cif
    campaigns/<name>/3_Ranked/!_Ranked.csv          the accepted designs, best-first
    campaigns/<name>/3_Ranked/<design>.cif

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

Rows are keyed by BindCraft2's own ``design`` identity rather than by position:
``!_Ranked.csv`` is re-sorted by ``i_pDAE`` after every acceptance, so an index into
it is not stable across a re-collect, while the design name (campaign, modality,
length and recipe hash) is. The submitter sets ``campaign_name`` to the parent row's
name, so a design name already starts with its parent; the prefix is only added here
when a --set override has changed that.

Each accepted complex is converted to PDB (downstream tools consume PDB) and becomes
the row's ``<leaf>_path``; the source mmCIF stays available as ``<leaf>_cif_path``.
Note the structure is the **complex**, binder + target, not the binder alone.

Multi-target campaigns (``sapia run bindcraft2 --targets ...``)
--------------------------------------------------------------
A design is still **one row** -- one binder sequence is one entity, however many
targets it was optimised against. What changes is that every reading becomes several.
Two upstream behaviours, both silently mis-read by a single-target collector:

* **File fan-out.** With two or more prepared states the written filename gains a
  ``_<target name>`` suffix (``campaign_output.accepted_state_suffixes``), and that
  suffix is *empty* for a single target -- so the filename shape changes the moment a
  second target appears. One design writes one complex per target. File order is
  alphabetical while the CSV's target order is by descending weight, so "the first
  file" and "the first value in a cell" are not the same target.
* **Semicolon cells.** ``campaign_output.target_ordered_row`` collapses every
  per-target reading into one ``;``-joined string (``0.82;0.79;0.31``), ordered by
  descending weight then name, with **empty positions kept** for missing readings, and
  adds the row's own key: ``targets`` and ``target_weights``.

So this collector splits each ``;`` cell into one column per target, keyed on the
row's own ``targets`` cell -- never on position, never on alphabetical order:

    bindcraft2_i_pTM_hEGFR        an ON-target reading
    bindcraft2_i_pTM_off_hERBB2   an OFF-target (detarget) reading
    bindcraft2_path_hEGFR         that target's complex
    bindcraft2_path               the HIGHEST-WEIGHT on-target complex

An off-target's columns carry an ``off_`` infix because upstream deliberately excludes
detargets from its own ranking (``on_target_mean``), and averaging or comparing a
detarget reading against a binding reading inverts its meaning. ``run_bindcraft2.py``
refuses a target actually named ``off_*`` so the infix stays unambiguous.

**Error contracts.** A cell whose ``;`` value count disagrees with the number of names
in ``targets`` makes that design an ``error`` status -- never a best-effort parse. A
missing per-target complex is ``NA``, never a substituted sibling file. An empty
position in a ``;`` cell is ``NA``, never ``0``.

Safe to re-run: rows are rebuilt from the CSV and the structures on disk.

Usage:
    sapia collect bindcraft2 outputs/RUN --table <the table the run reserved>
    sapia collect bindcraft2 outputs/RUN --table <table> --stage refolded --force
"""

import json
import re
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

# BindCraft2's multi-target row format (campaign_output.py at v1.0.3):
# TARGET_VALUE_SEPARATOR / TARGET_NAME_COLUMN / TARGET_WEIGHT_COLUMN. Both key columns
# are written ONLY when the campaign prepared two or more target states, so their
# presence is itself the single/multi discriminator.
TARGET_SEP = ";"
TARGETS_COL = "targets"
TARGET_WEIGHTS_COL = "target_weights"

# The infix an off-target's columns carry. Must match OFF_TARGET_PREFIX in
# run_bindcraft2.py, which refuses a target whose own name starts with it.
OFF_TARGET_PREFIX = "off_"

# Columns that hold a `;` for reasons of their own and are NOT per-target readings.
# `Timing` is the dangerous one: `timing_stamp` joins `worker=0;start=...;design=...`
# with the same separator, so on a three-target campaign it would split cleanly into
# three plausible-looking columns of nonsense. The rest are upstream's own
# TRAILING_METADATA_COLUMNS / SEQUENCE_COLUMNS / TEXT_COLUMNS, which `target_ordered_row`
# leaves in the shared part of the row.
NEVER_SPLIT = frozenset(
    {
        "Timing",
        "failed_filters",
        "hash",
        "phase",
        "round",
        "terminated",
        "autotuned",
        "rank",
        "design",
        "trajectory",
        "length",
        "outcome",
        "Binder_Sequence",
        "Interface_Binder_Residues",
        "Interface_Target_Residues",
        "settings_core",
        "settings_modality",
        "settings_property",
        "settings_target",
        "settings_overrides",
        TARGETS_COL,
        TARGET_WEIGHTS_COL,
    }
)

# BindCraft2 stamps accepted structures with a `_bindcraft` mmCIF category; the two
# span fields below are written against the letters of the OUTPUT file, so they name
# the binder chain as a reader will see it. See campaign_output.designed_span_stamp.
_STAMP_SPANS = re.compile(
    r"^_bindcraft\.(?:redesigned_residues|paratope_residues)\s+(\S+)\s*$", re.MULTILINE
)
_SPAN_CHAIN = re.compile(r"([A-Za-z])\d")


class TargetParseError(ValueError):
    """A multi-target row this collector refuses to guess at. Becomes an error status
    on that design's row rather than a best-effort parse."""


class TargetRef(NamedTuple):
    """One target of a multi-target campaign, in the CSV's own order.

    ``off`` follows upstream exactly: ``resolve_prepared_states`` forces a detarget's
    weight negative, so a negative weight in ``target_weights`` IS a detarget. A weight
    of exactly 0 is an on-target that contributes nothing to the loss.
    """

    name: str
    weight: float
    off: bool

    @property
    def label(self) -> str:
        """The column-name fragment: ``off_<name>`` for a detarget, ``<name>`` else."""
        return f"{OFF_TARGET_PREFIX}{self.name}" if self.off else self.name


def parse_target_order(row: pd.Series) -> list[TargetRef]:
    """The row's targets, in its own order, or ``[]`` for a single-target campaign.

    The row is the key, not the run's flags and not the files on disk -- upstream can
    legitimately prepare states the submit never named (a FASTA target cropped into
    ``<name>_epitope_<i>``), and the ``;`` cells are ordered to match this list.
    """
    if TARGETS_COL not in row.index:
        return []
    raw_names = _cell(row, TARGETS_COL)
    if not raw_names:
        return []
    names = [n for n in raw_names.split(TARGET_SEP) if n]
    if len(names) < 2:
        return []
    raw_weights = _cell(row, TARGET_WEIGHTS_COL)
    weights = [w for w in raw_weights.split(TARGET_SEP) if w != ""] if raw_weights else []
    if len(weights) != len(names):
        raise TargetParseError(
            f"'{TARGETS_COL}' names {len(names)} target(s) ({raw_names!r}) but "
            f"'{TARGET_WEIGHTS_COL}' holds {len(weights)} weight(s) ({raw_weights!r}); "
            f"without the weights an on-target cannot be told from a detarget"
        )
    try:
        parsed = [float(w) for w in weights]
    except ValueError as e:
        raise TargetParseError(
            f"'{TARGET_WEIGHTS_COL}' {raw_weights!r} is not a list of numbers ({e})"
        ) from e
    return [TargetRef(n, w, w < 0) for n, w in zip(names, parsed)]


def _cell(row: pd.Series, col: str) -> str:
    """A CSV cell as a string, with an absent column or a NaN/None read as ''."""
    if col not in row.index:
        return ""
    value = row[col]
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):  # array-like: not an NA scalar
        pass
    return str(value)


def split_target_cell(col: str, value: Any, n_targets: int) -> list[Any] | None:
    """One CSV cell split into its per-target readings, or ``None`` if it is shared.

    Returns a list of length ``n_targets`` whose empty positions are ``pd.NA`` -- an
    unmeasured reading, **not** zero (upstream's own docs say so, and a 0 i_pTM reads
    as "does not bind" rather than "was not measured").

    Raises TargetParseError when a cell that must be per-target has the wrong number of
    values: ``target_ordered_row`` always joins over *every* target name, so a mismatch
    means the row and its ``targets`` cell disagree, and there is no safe way to say
    which value belongs to which target.
    """
    if col in NEVER_SPLIT or not isinstance(value, str) or TARGET_SEP not in value:
        return None
    parts = value.split(TARGET_SEP)
    if len(parts) != n_targets:
        raise TargetParseError(
            f"column {col!r} holds {len(parts)} ';'-separated value(s) ({value!r}) for "
            f"{n_targets} target(s)"
        )
    return [_number_or_na(p) for p in parts]


def _number_or_na(part: str) -> Any:
    """A per-target reading: a float when it is one, the raw text when it is not, and
    ``pd.NA`` for an empty position (never 0.0)."""
    text = part.strip()
    if not text:
        return pd.NA
    try:
        return float(text)
    except ValueError:
        return text


def _unambiguous_structure(search_dir: Path, base: str) -> Path | None:
    """``<search_dir>/<base>.cif`` -- the exact name and nothing else.

    A single-target campaign writes no suffix at all (``accepted_state_suffixes``
    returns ``''`` below two prepared states), so the exact name is the only correct
    one. An earlier version of this collector fell back to ``glob(f'{base}*.cif')`` and
    took the first hit, which silently attached ``<design>_seq10``'s structure to
    ``<design>_seq1`` whenever the latter's own file was missing, and attached an
    arbitrary target's complex in a multi-target campaign. There is no safe glob here:
    a design with no file of its own has no structure.
    """
    exact = search_dir / f"{base}.cif"
    return exact if exact.is_file() else None


class Stage(NamedTuple):
    """Where a stage keeps its table and its structures.

    ``search`` is the folder holding a design's structures and ``base`` the filename
    stem upstream builds them from; a multi-target complex is exactly
    ``<search>/<base>_<target name>.cif`` (``accepted_state_suffixes``), which is why
    the two are kept apart rather than baked into one finder.
    """

    folder: str  # under project_folder
    table: str  # its CSV, inside that folder
    search: Callable[[Path, str], Path]  # (stage_dir, design) -> structure folder
    base: Callable[[str], str]  # design -> filename stem


STAGES: dict[str, Stage] = {
    # 1_Trajectories keeps a folder per design, and the filename carries the design
    # name twice (campaign_output.trajectory_output_path prefixes it).
    "trajectories": Stage(
        "1_Trajectories",
        "!_Trajectories.csv",
        lambda stage_dir, design: stage_dir / design,
        lambda design: f"{design}_trajectory",
    ),
    "refolded": Stage(
        "2_Refolded",
        "!_Refolded.csv",
        lambda stage_dir, _design: stage_dir / "Complexes",
        lambda design: design,
    ),
    "ranked": Stage(
        "3_Ranked",
        "!_Ranked.csv",
        lambda stage_dir, _design: stage_dir,
        lambda design: design,
    ),
}

# The column carrying the design's identity, and the one carrying its sequence.
DESIGN_COL = "design"
SEQUENCE_COL = "Binder_Sequence"
# 1_Trajectories only: the stage the attempt stopped at, blank when it completed.
TERMINATED_COL = "terminated"

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
    # interface
    "Interface_Residues",
    "Interface_BuriedArea",
    "Hotspot_Contact_Fraction",
    # developability
    "Surface_Hydrophobicity",
    "Binder_Length",
    "Binder_Net_Charge",
    "Binder_Free_Cysteines",
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


class BindCraft2CollectArgs(CollectArgs):
    stage: str
    metrics: str
    all_metrics: bool


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
        "the campaign accepted. Override with -l/--dir-label to keep two stages of "
        "one run in separate columns.",
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
        "campaign multiplies the per-target ones by its target count.",
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


def binder_chain_letters(cif: Path, pdb: Path | None) -> tuple[Any, str]:
    """The chain letter(s) the binder carries in the written complex, and where from.

    Upstream writes chains sorted on ``(is_binder_chain(name), name)``
    (``protein.written_chains``), so the **binder is last** -- ``B`` behind a
    single-chain target, ``C`` behind a two-chain one. The ``binder_chain`` setting
    (default ``null`` -> ``"binder"``) is the campaign's internal name and is NOT that
    letter; reading it instead is how a downstream `--chains-to-design` ends up
    redesigning the target.

    Two sources, best first:

    * ``stamp`` -- ``_bindcraft.redesigned_residues`` / ``paratope_residues`` in the
      accepted mmCIF, whose spans upstream writes against the output letters. Exact,
      and correct for an oligomeric binder (several letters).
    * ``last_chain`` -- the last chain in the written structure. Right for a single
      binder chain, and an UNDER-count when ``copies > 1``.

    Returns ``(pd.NA, "none")`` when neither works, never a guess.
    """
    try:
        text = cif.read_text(errors="replace")
    except OSError:
        text = ""
    letters = {
        match.group(1)
        for value in _STAMP_SPANS.findall(text)
        for span in value.strip("'\"").split(",")
        for match in [_SPAN_CHAIN.match(span)]
        if match
    }
    if letters:
        return "".join(sorted(letters)), "stamp"
    last = _last_pdb_chain(pdb) if pdb is not None else None
    return (last, "last_chain") if last else (pd.NA, "none")


def _last_pdb_chain(pdb: Path) -> str | None:
    """The chain id of the last ATOM/HETATM record of a PDB file."""
    try:
        chains = [
            line[21]
            for line in pdb.read_text(errors="replace").splitlines()
            if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21
        ]
    except OSError:
        return None
    return chains[-1].strip() or None if chains else None


def _shared_row_data(row: pd.Series, columns: list[str]) -> dict[str, Any]:
    """The columns that are one value per design however many targets there are.

    Values are passed through exactly as pandas read them -- this is the single-target
    collector's behaviour, unchanged.
    """
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


def collect_bindcraft2(ctx: CollectCtx[BindCraft2CollectArgs]) -> CollectEach:
    """Per-campaign BindCraft2 collector. A create tool: the framework iterates the
    ready parents (table rows, or the root design groups the run recorded) and this
    mints one child row per design the campaign produced, carrying ``parent`` for
    lineage. A campaign with no results yields no rows. The framework stamps
    status/path/parent_name from each Collected."""
    meta = _run_meta(ctx.out_dir)
    stage_name = resolve_stage(ctx, meta)
    stage = STAGES[stage_name]
    campaigns_root = resolve_campaigns_root(ctx, meta)
    run_dir = ctx.args.run_dir
    extra_metrics = [m.strip() for m in ctx.args.metrics.split(",") if m.strip()]
    submitted = list(meta.get("submitted_targets") or [])
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

        # archive_trajectories packs each design folder into <design>.zip, and the
        # structures then exist only inside it. Say that, rather than reporting a
        # campaign's worth of designs as having no structure.
        if stage_name == "trajectories" and any(stage_dir.glob("*.zip")):
            print(
                f"{d.name}: {stage_dir} holds archived trajectories — run "
                f"`bindcraft unarchive {campaign_dir}` before collecting"
            )

        columns = (
            [
                c
                for c in df.columns
                if c not in (DESIGN_COL, SEQUENCE_COL, TARGETS_COL, TARGET_WEIGHTS_COL)
            ]
            if ctx.args.all_metrics
            else [c for c in (*CORE_METRICS, *extra_metrics) if c in df.columns]
        )

        # Sorted by identity, not by rank: !_Ranked.csv is re-sorted after every
        # acceptance, so collecting in file order would be collecting in an order
        # that changes under a resumed campaign.
        n = 0
        ordered = df.sort_values(DESIGN_COL).set_index(DESIGN_COL)
        for design, row in ordered.iterrows():
            design = str(design)
            row_name = design if design.startswith(d.name) else f"{d.name}_{design}"
            search_dir, base = stage.search(stage_dir, design), stage.base(design)

            try:
                targets = parse_target_order(row)
            except TargetParseError as e:
                print(f"{d.name}: {design} {e}")
                yield Collected(
                    name=row_name,
                    parent=d.name,
                    path="",
                    status=f"error: {e}",
                    data=_shared_row_data(row, columns),
                )
                n += 1
                continue

            if targets:
                if submitted and set(submitted) - {t.name for t in targets}:
                    missing = sorted(set(submitted) - {t.name for t in targets})
                    if d.name not in warned_mismatch:
                        warned_mismatch.add(d.name)
                        print(
                            f"{d.name}: the CSV's targets do not include {missing} "
                            f"which the run submitted; splitting on the CSV's own "
                            f"`{TARGETS_COL}` cell, which is what the values are "
                            f"ordered by"
                        )
                emitted = _multi_target_row(
                    row, columns, targets, search_dir, base, as_pdb
                )
            else:
                emitted = _single_target_row(row, columns, search_dir, base, as_pdb)

            if emitted is None:
                print(
                    f"{d.name}: {design} has no structure under {search_dir}, skipping"
                )
                continue
            data, path, status = emitted
            if stage_name == "trajectories":
                # `terminated` is blank when the attempt ran to completion, so its
                # NA means success -- the opposite of NA everywhere else in a
                # prosapia table. Carry the polarity explicitly so a --filter reads
                # `bindcraft2_completed == True` instead of testing for a blank.
                data["completed"] = not _cell(row, TERMINATED_COL).strip()
            if status != "OK":
                print(f"{d.name}: {design} {status}")
            yield Collected(
                name=row_name, parent=d.name, path=path, status=status, data=data
            )
            n += 1
        print(f"{d.name}: OK ({n} {stage_name} design(s))")

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
