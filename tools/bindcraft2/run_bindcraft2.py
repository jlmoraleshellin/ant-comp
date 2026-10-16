#!/usr/bin/env python3
"""
Submit BindCraft2 binder-design campaigns, one array task per target.

BindCraft2 is campaign-driven, not design-driven: one ``bindcraft design
<settings.json>`` process takes a target, hallucinates binder backbones with AF2,
redesigns their sequences with ProteinMPNN, refolds and filters the candidates,
and keeps going until ``number_of_final_designs`` have been accepted or
``max_trajectories`` attempts are spent. So the unit of work here is **one campaign
per target**, not one task per design -- the designs only exist once the campaign
has run, which is why this is a ``create`` tool and why the row count is known only
at collect time.

This tool writes each campaign's settings JSON itself (from the flags below, merged
with ``--extra-settings``) and points ``project_folder`` at
``<out_dir>/campaigns/<name>/``. Everything BindCraft2 exposes that has no dedicated
flag is reachable through ``--extra-settings`` (a file, merged into every campaign)
or ``--set KEY=VALUE`` (verbatim, per run).

Hotspots, coldspots and binder lengths are authored in BindCraft2's own syntax, with
``{expr}`` placeholders resolved per-design against the table lineage (integers, bare
column names, and + - * // arithmetic; see resolve_expr):

    --hotspots 'A54,A56,A66-70'            # literal, BindCraft2 syntax untouched
    --hotspots 'A{epitope_start}-{epitope_end}'   # resolved per row

Usage:
    # child run: one campaign per target structure already in a table
    sapia run bindcraft2 outputs/RUN --table table0 -i pdb_path \\
        --hotspots 'A54,A56,A66,A115' --binder-lengths 60-100 --num-designs 10

    # root run: one campaign against a target not in any table yet
    sapia run bindcraft2 outputs/RUN \\
        --target-pdb targets/PDL1.pdb --chains A \\
        --hotspots 'A54,A56' --modality VHH --property humanize

    # root run against a target BindCraft2 ships
    sapia run bindcraft2 outputs/RUN --shipped-target hPDL1 --num-designs 10

    # backbones only: stop before BindCraft2's own ProteinMPNN, to redesign the
    # sequences with this workspace's tools instead
    sapia run bindcraft2 outputs/RUN --target-pdb targets/PDL1.pdb \\
        --trajectory-only --max-trajectories 40 --hotspots 'A54,A56'

    # ONE binder optimised jointly against SEVERAL targets (multi-specificity),
    # optionally with an off-target to counter-select against
    sapia run bindcraft2 outputs/RUN --targets targets/egfr_pair.yaml --num-designs 10

Multi-specificity (``--targets``)
---------------------------------
BindCraft2 v1.0.3 designs one binder against a *list* of targets, not a loop over
targets: ``loss.py`` puts one loss instance per (loss x target) into a single weighted
sum, ``multitarget_merged_gradients`` (default true) refuses a sequence update until
every target has contributed, and ``multitarget_tied_redesign`` (default true) ties the
ProteinMPNN redesign across them. That is a different experiment from N campaigns, and
``--targets`` is the only way to ask for it here.

``--targets`` takes a YAML/JSON file: a list of target mappings (or ``{targets: [...]}``),
each carrying only BindCraft2's own seven per-target keys -- ``name``, ``target_path``,
``chains``, ``hotspots``, ``coldspots``, ``weight``, ``objective``. A negative ``weight``
or ``objective: detarget`` makes a target a **counter-selection**: the campaign is pushed
away from binding it.

    # targets/egfr_pair.yaml
    - name: hEGFR
      target_path: targets/hEGFR_d3.pdb
      chains: A
      hotspots: A355,A356,A440,A441
      weight: 1.0
    - name: mEGFR
      target_path: targets/mEGFR_d3.pdb
      chains: A
      hotspots: A355,A356,A440,A441
      weight: 1.0
    - name: hERBB2            # off-target: do NOT bind this
      target_path: targets/hERBB2.pdb
      chains: A
      weight: -0.5

What it does NOT do: it does not loop over targets (that is a child run with one row per
target), it does not merge two campaigns after the fact, and it does not make the
per-target numbers comparable across campaigns.
"""

import difflib
import json
import re
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, NamedTuple, cast

import yaml

from prosapia.core import CommonArgs, ManifestCtx, build_tool_leaf
from prosapia.core.data_manager import LookupFn
from prosapia.core.executors import volume_path
from prosapia.utils import resolve_template

TOOL_NAME = "bindcraft2"

# Layout this tool imposes on its out_dir. collect_bindcraft2.py reads the same
# convention -- keep the two in step.
CAMPAIGNS_DIRNAME = "campaigns"
SETTINGS_DIRNAME = "settings"

# Design-property presets BindCraft2 ships under settings/property/, each turned on
# by its own ``--<name>`` flag on the bindcraft CLI. Kept as a list so --property
# validates against it instead of forwarding a typo that bindcraft would swallow as
# a stray path argument.
PROPERTY_PRESETS = (
    "bigbang",
    "disulfide_staple",
    "forced_targeting",
    "humanize",
    "initial_guess",
    "mixed_topology",
    "protease_stable",
    "termini_accessible",
    "termini_together",
)

# BindCraft2's WHOLE per-target schema (settings.py:38, TARGET_SETTING_NAMES at v1.0.3).
# Copied rather than imported: bindcraft is not installed where `sapia run` executes, and
# a typo has to fail at submit time rather than once per container. Bump with the pin in
# modal_image.py.
TARGET_SETTING_NAMES = (
    "name",
    "target_path",
    "chains",
    "hotspots",
    "coldspots",
    "weight",
    "objective",
)
TARGET_OBJECTIVES = ("target", "detarget")

# A target's name is not free text here: upstream uses it as a metric-state suffix
# (`i_pTM.<name>`), as a written-filename suffix (`<design>_<name>.cif`) and as an entry
# in the `;`-joined `targets` cell, and collect_bindcraft2.py turns it into a column-name
# fragment. So anything that would break one of those is refused up front.
_TARGET_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")

# collect_bindcraft2.py names an off-target's columns `<metric>_off_<name>`, so a target
# actually CALLED `off_x` would make `i_pTM_off_x` ambiguous. Refuse the collision here
# rather than let it read as an off-target downstream.
OFF_TARGET_PREFIX = "off_"

# A {expr} placeholder island, resolved per-design up the table lineage. Meaningless
# in a root run (no table), so we reject it there.
_HAS_PLACEHOLDER = re.compile(r"\{[^}]*\}")


class SettingsConfigError(ValueError):
    """A run-wide settings misconfiguration that applies to every campaign (e.g. an
    --extra-settings key colliding with a dedicated flag). Unlike a per-row error it
    is not swallowed by the warn-and-skip loop -- it fails the whole submit up front."""


class BindCraft2Args(CommonArgs):
    trajectory_only: bool
    reuse_campaigns: str | None
    targets: str | None
    target_pdb: Path | None
    shipped_target: str | None
    chains: str | None
    hotspots: str | None
    coldspots: str | None
    binder_lengths: str | None
    num_designs: int | None
    max_trajectories: int | None
    modality: str | None
    property: list[str]
    core: str | None
    campaign_seed: int | None
    no_resume: bool
    extra_settings: str | None
    set: list[str]


def add_run_bindcraft2_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--trajectory-only",
        action="store_true",
        help="Stop each campaign after the AF2 hallucination stage: produce "
        "BACKBONES and nothing else, with no ProteinMPNN redesign, no refold, no "
        "filtering and no accepted design. Use it to hand the backbones to this "
        "workspace's own sequence designers (`atomium`, `proteinmpnn`) instead of "
        "BindCraft2's. Implies `save_design_trajectory` so the structures are "
        "written, and makes --max-trajectories the only budget (BindCraft2 defaults "
        "it to 100). `sapia collect` picks the matching stage up automatically.",
    )
    parser.add_argument(
        "--reuse-campaigns",
        type=str,
        default=None,
        metavar="TABLE[:LABEL]",
        help="Submit nothing; reserve a table over the campaigns an EARLIER "
        "bindcraft2 run already wrote, named by the table it collected into (and its "
        "-l label, if any). Use it to collect a second stage of one campaign into a "
        "second table -- e.g. the backbones into one and BindCraft2's own accepted "
        "designs into another, from the same GPU hours. Pass the SAME -t the "
        "original run used, so the design groups still line up, and pick the stage "
        "with `sapia collect --stage`. Incompatible with every target and campaign "
        "flag, since no campaign is run.",
    )
    parser.add_argument(
        "--targets",
        type=str,
        default=None,
        metavar="FILE",
        help="Path to a YAML/JSON file listing SEVERAL targets ONE binder is optimised "
        "against jointly (BindCraft2's own `targets` list -- true multi-specificity, "
        "not a loop). Either a bare list of mappings or `{targets: [...]}`. Each entry "
        "takes only BindCraft2's seven per-target keys: " + ", ".join(TARGET_SETTING_NAMES)
        + ". A negative `weight` (or `objective: detarget`) counter-selects against an "
        "off-target. Validated at submit time: unknown keys, duplicate or unusable "
        "names, bad objectives and missing target_path files all fail here rather than "
        "once per container. String values may embed {expr} placeholders resolved "
        "per-design up the lineage. Mutually exclusive with --target-pdb / "
        "--shipped-target, and with the top-level --chains/--hotspots/--coldspots (in "
        "this mode those are per-target sub-keys). Without --table the design group is "
        "named `<file stem>_bc2`; with --table one campaign per ready row is run "
        "against the SAME target list, and the row's --input-column is NOT used as a "
        "target -- it only decides which rows are ready and supplies the lineage that "
        "{expr} resolves against.",
    )
    parser.add_argument(
        "--target-pdb",
        type=Path,
        default=None,
        help="Single target structure to design binders against in a ROOT run (no "
        "--table): a PDB/mmCIF/FASTA not in any table yet. Only valid without "
        "--table (with a table, targets come from --input-column). The design group "
        "is named `<stem>_bc2`.",
    )
    parser.add_argument(
        "--shipped-target",
        type=str,
        default=None,
        help="Name of a target BindCraft2 ships (hPDL1, hPD1, mPDL1, hIL2R, hIL7RA, "
        "dynorphin_a; `bindcraft design --list-targets` is authoritative). ROOT runs "
        "only, and mutually exclusive with --target-pdb. The design group is named "
        "`<name>_bc2`.",
    )
    parser.add_argument(
        "--chains",
        type=str,
        default=None,
        help="Target chains to design against (per-target `chains`, e.g. 'A' or "
        "'A,B'). Omitted by default (BindCraft2 uses every chain in the file).",
    )
    parser.add_argument(
        "--hotspots",
        type=str,
        default=None,
        help="Target residues the binder should contact (per-target `hotspots`), in "
        "BindCraft2's own syntax: comma-separated residues and ranges, chain-prefixed "
        "(e.g. 'A54,A56,A66-70'). May embed {expr} placeholders resolved per-design "
        "up the lineage. Omitted by default (BindCraft2 picks the epitope itself).",
    )
    parser.add_argument(
        "--coldspots",
        type=str,
        default=None,
        help="Target regions to avoid contacting (per-target `coldspots`), same "
        "syntax as --hotspots. Omitted by default.",
    )
    parser.add_argument(
        "--binder-lengths",
        type=str,
        default=None,
        help="Binder size (`binder_lengths`): 'N' for one length, 'min-max' for a "
        "range drawn from per trajectory (e.g. '80' or '60-100'). May embed {expr} "
        "placeholders. Omitted by default (BindCraft2's own default, or the "
        "modality preset's).",
    )
    parser.add_argument(
        "--num-designs",
        type=int,
        default=None,
        help="Accepted designs to stop the campaign at (`number_of_final_designs`). "
        "This is the number of CHILD ROWS a campaign aims to produce, not a batch "
        "size -- BindCraft2 keeps spending trajectories until it has them. Omitted "
        "by default.",
    )
    parser.add_argument(
        "--max-trajectories",
        type=int,
        default=None,
        help="Design attempts to spend before giving up (`max_trajectories`). The "
        "real cost knob: a campaign runs until --num-designs are accepted OR this "
        "many attempts are spent. Omitted by default.",
    )
    parser.add_argument(
        "--modality",
        type=str,
        default=None,
        help="Binder format (-> `bindcraft design --modality`): binder, VHH, "
        "peptide, cyclic_peptide, ARP, scFv, Fab, large_binder, homo_oligomer, "
        "multidomain, induced_fit, fold_switch. Comma-separated to combine. "
        "Defaults to BindCraft2's own default (binder).",
    )
    parser.add_argument(
        "--property",
        action="append",
        default=[],
        choices=PROPERTY_PRESETS,
        metavar="NAME",
        help="Design-property preset to switch on (-> `bindcraft design --<name>`). "
        "Repeatable. One of: " + ", ".join(PROPERTY_PRESETS) + ".",
    )
    parser.add_argument(
        "--core",
        type=str,
        default=None,
        help="Core profile applied under every preset (-> `bindcraft design --core`), "
        "e.g. 'benchmark' for a reproducible run. Omitted by default.",
    )
    parser.add_argument(
        "--campaign-seed",
        type=int,
        default=None,
        help="Seed every trajectory is drawn from (`campaign_seed`). Set it with "
        "--core benchmark for a reproducible campaign. Omitted by default.",
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Start each campaign from scratch instead of carrying on into a "
        "project_folder already written. By default this tool sets `resume: true`, "
        "so re-running a submit continues the campaigns it already started (the "
        "framework's own resume filter can't help here -- a create tool's status "
        "column lives in the child table, not the one it reads).",
    )
    parser.add_argument(
        "--extra-settings",
        type=str,
        default=None,
        help="Path to a YAML or JSON file: a mapping of (extra) BindCraft2 campaign "
        "settings merged into EVERY campaign's settings file (e.g. objective, "
        "aa_bias, min_iptm_final, save_design_trajectory). String values may embed "
        "{expr} placeholders resolved per-design up the lineage. "
        "`bindcraft design --list-settings` names every setting it accepts.",
    )
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Extra `bindcraft design --set` override, appended verbatim and applied "
        "over the generated settings file. Repeatable. Escape hatch for settings "
        "without a dedicated flag. Must contain no spaces (the task script "
        "word-splits these tokens); use --extra-settings for anything richer.",
    )


def load_extra_settings(extra_settings: str | None) -> dict[str, Any]:
    """Parse the --extra-settings YAML/JSON file into a mapping of campaign settings.

    Returns ``{}`` when unset or empty. YAML is a JSON superset, so ``yaml.safe_load``
    parses both. Raises FileNotFoundError for a missing path and TypeError if the
    top level isn't a mapping.
    """
    if not extra_settings:
        return {}
    path = Path(extra_settings)
    if not path.is_file():
        raise FileNotFoundError(f"--extra-settings file not found: {extra_settings}")
    data = yaml.safe_load(path.read_text())
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise TypeError(
            f"--extra-settings must be a mapping of campaign settings, got "
            f"{type(data).__name__}"
        )
    return data


class TargetsSpec(NamedTuple):
    """A validated ``--targets`` file: its entries, the directory it lives in (used as
    the fallback base for relative ``target_path``s) and its stem (the root design-group
    name). ``entries`` are still un-``{expr}``-resolved -- resolution is per campaign."""

    entries: list[dict[str, Any]]
    base_dir: Path
    stem: str


def _unknown_target_key(key: str) -> str:
    suggestion = difflib.get_close_matches(key, TARGET_SETTING_NAMES, n=1)
    hint = f"; did you mean {suggestion[0]!r}?" if suggestion else ""
    return f"targets[].{key!r} is not a BindCraft2 per-target setting{hint}"


def validate_target_entries(entries: list[Any]) -> list[dict[str, Any]]:
    """Check a ``--targets`` list against BindCraft2's per-target schema, and return it.

    Every check here exists because the alternative is a failure one container-hour
    later, or -- worse -- a number that reads fine and means something else:

    * unknown keys: upstream rejects them too (``reject_unrecognized_settings``), but
      only after the image is built and the campaign has started.
    * ``name`` charset and uniqueness: the name is a metric-state suffix, a filename
      suffix, an entry in the ``;``-joined ``targets`` cell and a column-name fragment.
      A ``.``, ``;``, ``/`` or space in it silently breaks one of those.
    * at least one attracting target: a list of nothing but detargets is a campaign with
      no binding objective at all -- it would run, and produce binders of nothing.

    Raises ValueError listing every problem at once.
    """
    problems: list[str] = []
    seen: dict[str, int] = {}
    attracting = 0
    for index, entry in enumerate(entries):
        where = f"targets[{index}]"
        if not isinstance(entry, dict):
            problems.append(f"{where} is a {type(entry).__name__}, not a mapping")
            continue
        for key in entry:
            if key not in TARGET_SETTING_NAMES:
                problems.append(f"{where}: {_unknown_target_key(str(key))}")

        name = entry.get("name")
        if not isinstance(name, str) or not name:
            problems.append(f"{where}.name is required and must be a non-empty string")
        elif not _TARGET_NAME_RE.match(name):
            problems.append(
                f"{where}.name {name!r} must be letters/digits/_/- starting with a "
                f"letter or digit. The name becomes a metric suffix (`i_pTM.{name}`), a "
                f"filename suffix and a table column name, so '.', ';', '/' and spaces "
                f"would all break something downstream."
            )
        elif name.startswith(OFF_TARGET_PREFIX):
            problems.append(
                f"{where}.name {name!r} starts with {OFF_TARGET_PREFIX!r}, which the "
                f"collector uses to mark an off-target's columns "
                f"(`bindcraft2_i_pTM_off_<name>`). Rename it."
            )
        elif name in seen:
            problems.append(f"{where}.name {name!r} repeats targets[{seen[name]}].name")
        else:
            seen[name] = index

        if not isinstance(entry.get("target_path"), str) or not entry["target_path"]:
            problems.append(
                f"{where}.target_path is required and must be a path to the target "
                f"structure (BindCraft2 reads `target['name']` and `target_path` for "
                f"every entry)."
            )

        for key in ("chains", "hotspots", "coldspots"):
            if key in entry and not isinstance(entry[key], str):
                problems.append(
                    f"{where}.{key} must be a string in BindCraft2's own syntax "
                    f"(e.g. 'A54,A56,A66-70'), got {type(entry[key]).__name__}"
                )

        objective = entry.get("objective", "target")
        if objective not in TARGET_OBJECTIVES:
            problems.append(
                f"{where}.objective {objective!r} must be one of "
                f"{', '.join(TARGET_OBJECTIVES)}"
            )

        weight = entry.get("weight", 1.0)
        if isinstance(weight, bool) or not isinstance(weight, (int, float)):
            problems.append(
                f"{where}.weight must be a number, got {type(weight).__name__}"
            )
        elif objective != "detarget" and float(weight) > 0:
            attracting += 1

    if entries and not problems and not attracting:
        problems.append(
            "every target is a detarget (objective: detarget, or weight <= 0), so the "
            "campaign has nothing to bind. Give at least one target a positive weight."
        )

    if problems:
        raise ValueError(
            "--targets file is not a valid BindCraft2 target list:\n  "
            + "\n  ".join(problems)
        )
    return cast(list[dict[str, Any]], entries)


def load_targets(targets: str) -> TargetsSpec:
    """Parse and validate the ``--targets`` YAML/JSON file.

    Accepts either a bare list of target mappings or the full ``{targets: [...]}`` shape
    a BindCraft2 settings file uses, so a user can lift the block straight out of one of
    upstream's examples (``examples/pdl1_crossreactive_detarget.json``).
    """
    path = Path(targets)
    if not path.is_file():
        raise FileNotFoundError(f"--targets file not found: {targets}")
    data = yaml.safe_load(path.read_text())
    if isinstance(data, dict):
        if "targets" not in data:
            raise ValueError(
                f"--targets {targets} is a mapping with no 'targets' key. Write either "
                f"a bare list of target mappings or {{targets: [...]}}."
            )
        extra = sorted(set(data) - {"targets"})
        if extra:
            raise ValueError(
                f"--targets {targets} carries {extra} alongside 'targets'. This file "
                f"describes targets only -- campaign-wide settings go in "
                f"--extra-settings."
            )
        data = data["targets"]
    if not isinstance(data, list) or not data:
        raise ValueError(
            f"--targets {targets} must hold a non-empty list of target mappings, got "
            f"{type(data).__name__}"
        )
    return TargetsSpec(validate_target_entries(data), path.parent, path.stem)


def resolve_target_path(raw: str, base_dir: Path, where: str) -> Path:
    """One target's ``target_path`` as an absolute path the task container can open.

    A relative path is read against the submit cwd first (the convention every other
    path flag here follows -- `sapia` runs from `/runs`), then against the --targets
    file's own directory, so a self-contained targets file next to its structures also
    works. Both attempts are reported when neither exists; the chosen one is printed, so
    which rule fired is never a guess.
    """
    candidate = Path(raw)
    if candidate.is_absolute():
        resolved = volume_path(candidate)
        if not resolved.exists():
            raise ValueError(f"{where}.target_path {resolved} does not exist")
        return resolved
    tried = [volume_path(candidate), volume_path(base_dir / candidate)]
    for resolved in tried:
        if resolved.exists():
            return resolved
    raise ValueError(
        f"{where}.target_path {raw!r} does not exist. Tried "
        + " and ".join(str(p) for p in tried)
    )


def _tree_has_placeholder(obj: Any) -> bool:
    """True if any string key/value anywhere in ``obj`` embeds a ``{expr}`` island.

    Only string leaves are inspected -- structural dict/list braces don't count.
    """
    if isinstance(obj, dict):
        return any(
            _tree_has_placeholder(k) or _tree_has_placeholder(v) for k, v in obj.items()
        )
    if isinstance(obj, list):
        return any(_tree_has_placeholder(v) for v in obj)
    if isinstance(obj, str):
        return bool(_HAS_PLACEHOLDER.search(obj))
    return False


def _resolve_tree(obj: Any, lookup: LookupFn, name: str) -> Any:
    """Recursively resolve ``{expr}`` placeholders in a parsed settings structure.

    Strings (and dict keys) pass through ``resolve_template``; dicts and lists are
    walked; other scalars (int/float/bool/None) are returned unchanged. So native
    YAML/JSON types survive except where a string embeds ``{expr}``.
    """
    if isinstance(obj, dict):
        return {
            resolve_template(str(k), lookup, name): _resolve_tree(v, lookup, name)
            for k, v in obj.items()
        }
    if isinstance(obj, list):
        return [_resolve_tree(v, lookup, name) for v in obj]
    if isinstance(obj, str):
        return resolve_template(obj, lookup, name)
    return obj


def _coerce_binder_lengths(spec: str) -> list[int]:
    """'80' -> [80]; '60-100' / '60,100' -> [60, 100] (BindCraft2's [min, max])."""
    parts = [p.strip() for p in re.split(r"[-,]", spec) if p.strip()]
    try:
        return [int(p) for p in parts]
    except ValueError as e:
        raise ValueError(
            f"--binder-lengths {spec!r} is not a length or a 'min-max' range"
        ) from e


def resolve_reuse_root(spec: str, run_dir: Path) -> Path:
    """``TABLE[:LABEL]`` -> the campaigns dir of an earlier bindcraft2 run.

    A run's outputs live at ``run_dir/<the table it reserved>/<leaf>/campaigns``, so
    naming that table (and the ``-l`` label it used, if any) is enough to find them.
    """
    table, _, label = spec.partition(":")
    if not table:
        raise ValueError(
            f"--reuse-campaigns {spec!r} names no table; expected TABLE[:LABEL], "
            f"e.g. 'table1' or 'table1:traj'."
        )
    root = run_dir / table / build_tool_leaf(TOOL_NAME, label) / CAMPAIGNS_DIRNAME
    if not root.is_dir():
        raise FileNotFoundError(
            f"--reuse-campaigns {spec!r} points at {root}, which does not exist. "
            f"Name the table the earlier bindcraft2 run collected into, plus its "
            f"-l label after a colon if it used one."
        )
    return root


def _reject_campaign_flags(args: BindCraft2Args) -> None:
    """A reuse run submits nothing, so anything describing a campaign is a mistake."""
    named = [
        flag
        for flag, value in (
            ("--targets", args.targets),
            ("--target-pdb", args.target_pdb),
            ("--shipped-target", args.shipped_target),
            ("--trajectory-only", args.trajectory_only),
            ("--hotspots", args.hotspots),
            ("--coldspots", args.coldspots),
            ("--chains", args.chains),
            ("--binder-lengths", args.binder_lengths),
            ("--num-designs", args.num_designs),
            ("--max-trajectories", args.max_trajectories),
            ("--modality", args.modality),
            ("--core", args.core),
            ("--campaign-seed", args.campaign_seed),
            ("--extra-settings", args.extra_settings),
            ("--property", args.property),
            ("--set", args.set),
        )
        if value
    ]
    if named:
        raise ValueError(
            f"--reuse-campaigns runs no campaign, so {', '.join(named)} would have "
            f"no effect. Drop them; the campaign was already run with its own "
            f"settings, and only `sapia collect --stage` still applies."
        )


def trajectory_only_run(args: BindCraft2Args, extra_fields: dict[str, Any]) -> bool:
    """Whether this run stops at backbones.

    Reads the flag *and* --extra-settings, because a user may set ``trajectory_only``
    in the settings file instead. The verbatim ``--set`` tokens are not inspected --
    setting it that way leaves the sidecar saying otherwise, and ``sapia collect``
    will default to the wrong stage.
    """
    return bool(args.trajectory_only or extra_fields.get("trajectory_only"))


def build_targets_block(
    spec: TargetsSpec, lookup: LookupFn, name: str
) -> list[dict[str, Any]]:
    """One campaign's ``targets`` list from a --targets file.

    ``{expr}`` is resolved per campaign (so a child run can carry per-row hotspots) and
    every ``target_path`` is made absolute and checked to exist *here*, after resolution
    -- a placeholder can change which file a row points at, so the existence check has
    to happen per campaign rather than once at parse time.
    """
    block: list[dict[str, Any]] = []
    for index, entry in enumerate(spec.entries):
        resolved = cast(dict[str, Any], _resolve_tree(entry, lookup, name))
        resolved["target_path"] = str(
            resolve_target_path(
                str(resolved["target_path"]), spec.base_dir, f"targets[{index}]"
            )
        )
        block.append(resolved)
    return block


def _build_settings(
    name: str,
    target_path: Path | None,
    args: BindCraft2Args,
    lookup: LookupFn,
    campaign_dir: Path,
    extra_fields: dict[str, Any],
    targets_spec: TargetsSpec | None = None,
) -> dict[str, Any]:
    """Build one campaign's BindCraft2 settings mapping.

    ``target_path`` is the (absolute) target structure, or ``None`` when the target
    comes from ``--targets``, ``--shipped-target`` or --extra-settings.
    ``targets_spec`` is the parsed --targets file, when one was given. Raises ValueError
    on a per-row problem (an unresolvable {expr}, a target file that isn't there) and
    the run-wide ``SettingsConfigError`` on an extra-settings collision.
    """
    # Settings this tool merely defaults, so --extra-settings can still turn them
    # off. Unlike tool_fields these are not collision-checked.
    defaults: dict[str, Any] = {}

    # Bookkeeping this tool owns outright: the campaign's identity and where it
    # writes. Not negotiable via --extra-settings, hence checked for collisions.
    tool_fields: dict[str, Any] = {
        "campaign_name": name,
        "project_folder": str(campaign_dir),
        "resume": not args.no_resume,
    }

    # Keep the fold each trajectory ended on. Under --trajectory-only it is the
    # entire output; under a full campaign it is what lets the same GPU hours be
    # collected a second time as backbones (see --reuse-campaigns). One CIF per
    # trajectory is cheap, so it is a default rather than a second flag --
    # --extra-settings can still turn it off.
    defaults["save_design_trajectory"] = True
    if args.trajectory_only:
        tool_fields["trajectory_only"] = True

    if targets_spec is not None:
        # Multi-specificity: the file IS the target list, verbatim in BindCraft2's own
        # schema. Nothing here invents a per-target field -- chains/hotspots/coldspots
        # are refused at the top level in this mode precisely so there is one source.
        tool_fields["targets"] = build_targets_block(targets_spec, lookup, name)
    elif target_path is not None:
        target: dict[str, Any] = {"name": name, "target_path": str(target_path)}
        if args.chains is not None:
            target["chains"] = resolve_template(args.chains, lookup, name)
        if args.hotspots is not None:
            target["hotspots"] = resolve_template(args.hotspots, lookup, name)
        if args.coldspots is not None:
            target["coldspots"] = resolve_template(args.coldspots, lookup, name)
        tool_fields["targets"] = [target]
    elif args.shipped_target is not None:
        # A shipped target is named rather than described: its own preset carries the
        # path and chains, so per-target overrides would have nowhere to land.
        tool_fields["target"] = args.shipped_target
        for flag, value in (
            ("--chains", args.chains),
            ("--hotspots", args.hotspots),
            ("--coldspots", args.coldspots),
        ):
            if value is not None:
                raise SettingsConfigError(
                    f"{flag} describes a target file and cannot be combined with "
                    f"--shipped-target (its preset already names the epitope). Pass "
                    f"the structure with --target-pdb, or override the preset with "
                    f"--extra-settings."
                )

    if args.binder_lengths is not None:
        tool_fields["binder_lengths"] = _coerce_binder_lengths(
            resolve_template(args.binder_lengths, lookup, name)
        )
    if args.num_designs is not None:
        tool_fields["number_of_final_designs"] = args.num_designs
    if args.max_trajectories is not None:
        tool_fields["max_trajectories"] = args.max_trajectories
    if args.campaign_seed is not None:
        tool_fields["campaign_seed"] = args.campaign_seed

    resolved_extra = _resolve_tree(extra_fields, lookup, name)
    collisions = sorted(set(tool_fields) & set(resolved_extra))
    if collisions:
        raise SettingsConfigError(
            f"setting(s) {collisions} set by both a dedicated flag and "
            "--extra-settings; remove them from one source"
        )

    return {**defaults, **resolved_extra, **tool_fields}


def _cli_flags(args: BindCraft2Args) -> str:
    """The run-wide ``bindcraft design`` flags (identical for every campaign).

    Presets and --set live on the command line rather than in the settings file
    because that is the interface BindCraft2 documents for them: --modality/--core
    name preset files to layer under the campaign, and --set is applied over it.
    """
    flags: list[str] = []
    if args.core is not None:
        flags += ["--core", args.core]
    if args.modality is not None:
        flags += ["--modality", args.modality]
    for prop in args.property:
        flags.append("--" + prop.replace("_", "-"))
    for assignment in args.set:
        flags += ["--set", assignment]
    return " ".join(flags)


def check_targets_exclusivity(args: BindCraft2Args) -> None:
    """``--targets`` owns the whole target description, so nothing may describe one too.

    Refused up front and by name rather than merged: ``--chains/--hotspots/--coldspots``
    are top-level *shorthands* for the single-target case, but in BindCraft2's schema
    they are per-target sub-keys. Accepting both would leave two sources for the
    epitope, and silently applying one of them to all three targets is exactly the kind
    of plausible-looking wrong answer this tool exists to avoid.
    """
    if args.targets is None:
        return
    named = [
        flag
        for flag, value in (
            ("--target-pdb", args.target_pdb),
            ("--shipped-target", args.shipped_target),
        )
        if value
    ]
    if named:
        raise SettingsConfigError(
            f"--targets already lists every target, so {', '.join(named)} would name "
            f"another one. Pass exactly one of --targets / --target-pdb / "
            f"--shipped-target."
        )
    per_target = [
        flag
        for flag, value in (
            ("--chains", args.chains),
            ("--hotspots", args.hotspots),
            ("--coldspots", args.coldspots),
        )
        if value is not None
    ]
    if per_target:
        raise SettingsConfigError(
            f"{', '.join(per_target)} cannot be combined with --targets: with several "
            f"targets these are PER-TARGET settings, and BindCraft2 reads them from "
            f"each entry of the targets list. Move them into the --targets file as the "
            f"`chains` / `hotspots` / `coldspots` key of the target they describe."
        )


def _target_groups(
    ctx: ManifestCtx[BindCraft2Args], targets_spec: TargetsSpec | None = None
) -> list[tuple[str, Path | None]]:
    """The (name, target_path) groups this run designs against: one per ready table
    row for a child run, a single group for a root run.

    With ``--targets`` the target_path is always ``None`` -- the targets come from the
    file, not from the table or from ``--target-pdb``.
    """
    args = ctx.args

    if args.table is not None:
        if args.target_pdb is not None or args.shipped_target is not None:
            raise ValueError(
                "--target-pdb/--shipped-target are only valid for a root run (no "
                "--table); with --table, targets come from the table's "
                "--input-column. Drop one of them."
            )
        groups: list[tuple[str, Path | None]] = []
        for name in ctx.ready.index:
            name = cast(str, name)
            if targets_spec is not None:
                # The row parameterises the campaign (it decides readiness and supplies
                # the lineage {expr} resolves against); the targets are the file's.
                groups.append((name, None))
                continue
            target_path = volume_path(str(ctx.ready.at[name, args.input_column]))
            if not target_path.exists():
                print(f"{name}: MISSING {target_path} (skipping)")
                continue
            groups.append((name, target_path))
        if targets_spec is not None and groups:
            print(
                f"--targets: {len(groups)} campaign(s), each designing ONE binder "
                f"against all {len(targets_spec.entries)} target(s) in "
                f"{targets_spec.stem}. The --input-column ({args.input_column}) is NOT "
                f"used as a target here -- it only selects the ready rows."
            )
        return groups

    # Root run: {expr} placeholders resolve up a table lineage this run doesn't have.
    if any(
        s and _HAS_PLACEHOLDER.search(s)
        for s in (args.hotspots, args.coldspots, args.chains, args.binder_lengths)
    ):
        raise ValueError(
            "--hotspots/--coldspots/--chains/--binder-lengths contain a {expr} "
            "placeholder, but this is a root run (no --table) with no table lineage "
            "to resolve it against. Use literal values, or run with --table."
        )
    if targets_spec is not None:
        if _tree_has_placeholder(targets_spec.entries):
            raise ValueError(
                "--targets contains a {expr} placeholder, but this is a root run (no "
                "--table) with no table lineage to resolve it against. Use literal "
                "values, or run with --table."
            )
        return [(f"{targets_spec.stem}_bc2", None)]

    if args.target_pdb is not None and args.shipped_target is not None:
        raise ValueError(
            "--target-pdb and --shipped-target both name a target; pass exactly one."
        )
    if args.target_pdb is not None:
        target_path = volume_path(args.target_pdb)
        if not target_path.exists():
            raise FileNotFoundError(f"--target-pdb {target_path} does not exist.")
        return [(f"{target_path.stem}_bc2", target_path)]
    if args.shipped_target is not None:
        return [(f"{args.shipped_target}_bc2", None)]

    raise ValueError(
        "A root run (no --table) needs a target: pass --target-pdb <file>, "
        "--targets <file> or --shipped-target <name>. With --table, targets come from "
        "--input-column."
    )


def _reuse_existing_campaigns(ctx: ManifestCtx[BindCraft2Args]) -> list[tuple[str, ...]]:
    """Reserve this run's table over an earlier run's campaigns, and submit nothing.

    The table and its out_dir are reserved by the driver before the manifest is
    built, so returning no rows still leaves a table for `sapia collect` to fill --
    it just fills it from the recorded campaigns rather than from this out_dir.
    """
    _reject_campaign_flags(ctx.args)
    root = resolve_reuse_root(ctx.args.reuse_campaigns or "", ctx.args.run_dir)
    ctx.write_meta(campaigns_root=str(volume_path(root)), trajectory_only=False)

    if ctx.args.table is None:
        # Root reuse: no parent table to iterate, so the campaign dirs on disk ARE
        # the design groups -- the same record a root campaign run would write.
        groups = sorted(p.name for p in root.iterdir() if p.is_dir())
        ctx.write_meta(root_designs=groups)
    else:
        groups = sorted(map(str, ctx.ready.index))

    print(f"Reusing {len(groups)} campaign(s) under {root}")
    print("Nothing to submit. Collect the stage you want, e.g.:")
    print("  sapia collect bindcraft2 <run_dir> -t <this table> --stage ranked")
    return []


def build_bindcraft2_manifest(
    ctx: ManifestCtx[BindCraft2Args],
) -> list[tuple[str, ...]]:
    if ctx.args.reuse_campaigns:
        return _reuse_existing_campaigns(ctx)

    check_targets_exclusivity(ctx.args)
    targets_spec = load_targets(ctx.args.targets) if ctx.args.targets else None
    extra_fields = load_extra_settings(ctx.args.extra_settings)

    if ctx.args.table is None and _tree_has_placeholder(extra_fields):
        raise ValueError(
            "--extra-settings contains a {expr} placeholder, but this is a root run "
            "(no --table) with no table lineage to resolve it against."
        )

    trajectory_only = trajectory_only_run(ctx.args, extra_fields)
    if trajectory_only and ctx.args.num_designs is not None:
        print(
            "--num-designs is meaningless with --trajectory-only: no design is ever "
            "accepted, so the budget is --max-trajectories alone (BindCraft2 "
            "defaults it to 100)."
        )
    # The stage `sapia collect` should read is decided here, not guessed there:
    # a trajectory-only run fills 1_Trajectories and leaves 3_Ranked empty, and a
    # collect that looked in the wrong place would report zero rows rather than a
    # mismatch. Same contract as input_column -- the run records, collect reads.
    ctx.write_meta(
        trajectory_only=trajectory_only,
        campaigns_root=str(volume_path(ctx.out_dir) / CAMPAIGNS_DIRNAME),
    )
    if targets_spec is not None:
        # Recorded so collect can say "the CSV names targets this run never submitted"
        # instead of trusting the CSV blindly. It is a cross-check, not the key: the
        # key is always the row's own `targets` cell, because upstream can legitimately
        # add states the submit never named (a FASTA target cropped into
        # `<name>_epitope_<i>`).
        ctx.write_meta(
            submitted_targets=[str(t["name"]) for t in targets_spec.entries],
            submitted_target_weights=[
                float(t.get("weight", 1.0)) for t in targets_spec.entries
            ],
        )

    groups = _target_groups(ctx, targets_spec)
    if not groups:
        return []

    settings_dir = ctx.out_dir / SETTINGS_DIRNAME
    settings_dir.mkdir(parents=True, exist_ok=True)
    campaigns_root = volume_path(ctx.out_dir) / CAMPAIGNS_DIRNAME

    flags = _cli_flags(ctx.args)
    manifest_rows: list[tuple[str, ...]] = []
    submitted: list[str] = []
    for name, target_path in groups:
        try:
            settings = _build_settings(
                name,
                target_path,
                ctx.args,
                ctx.lookup,
                campaigns_root / name,
                extra_fields,
                targets_spec,
            )
        except SettingsConfigError:
            # A run-wide misconfiguration hits every row identically: fail fast
            # instead of silently skipping the entire table.
            raise
        except ValueError as e:
            # One bad row (e.g. an unresolvable {expr}) shouldn't sink the whole
            # array: warn and skip it.
            print(f"{name}: {e} (skipping)")
            continue
        settings_json = settings_dir / f"{name}.json"
        settings_json.write_text(json.dumps(settings, indent=2))
        manifest_rows.append((name, str(volume_path(settings_json)), flags))
        submitted.append(name)

    if ctx.args.table is None:
        # A root run has no parent table for collect to iterate; record the group
        # names so `sapia collect` can find their campaigns and rebuild their rows.
        ctx.write_meta(root_designs=submitted)

    return manifest_rows
