#!/usr/bin/env python3
"""
Submit BindCraft2 binder-design campaigns, one array task per campaign.

BindCraft2 is campaign-driven, not design-driven: one ``bindcraft design
<settings.json>`` process takes a target (or SEVERAL targets), hallucinates binder
backbones with AF2, redesigns their sequences with ProteinMPNN, refolds and filters
the candidates, and keeps going until ``number_of_final_designs`` have been accepted
or ``max_trajectories`` attempts are spent. So the unit of work here is **one
campaign per task**, not one task per design -- the designs only exist once the
campaign has run, which is why this is a ``create`` tool and why the row count is
known only at collect time.

This tool writes each campaign's settings JSON itself (from the flags below, merged
with ``--extra-settings``) and points ``project_folder`` at
``<out_dir>/campaigns/<name>/``. Everything BindCraft2 exposes that has no dedicated
flag is reachable through ``--extra-settings`` (a file, merged into every campaign)
or ``--set KEY=VALUE`` (verbatim, per run). ``bindcraft design --list-settings``
names all of them.

Targets, hotspots, coldspots and binder lengths are authored in BindCraft2's own
syntax, with ``{expr}`` placeholders resolved per-design against the table lineage
(integers, bare column names, and + - * // arithmetic; see resolve_expr):

    --hotspots 'A54,A56,A66-70'                   # literal, BC2 syntax untouched
    --hotspots 'A{epitope_start}-{epitope_end}'   # resolved per row

MULTI-TARGET AND DETARGETING
----------------------------
A campaign may carry several targets: orthologs to be bound cross-reactively, and
off-targets to be avoided. Each is one ``--target`` (root) or ``--extra-target``
(added beside the table's own target), written as ``;``-separated ``key=value``
fields in BindCraft2's own ``targets[]`` vocabulary -- ``;`` rather than ``,``
because hotspot lists already use commas:

    --target 'name=hPDL1;path=t/hPDL1.pdb;chains=A;hotspots=A54,A56;weight=1'
    --target 'name=hPD1;path=t/hPD1.pdb;objective=detarget;weight=-0.5'

A negative ``weight`` also selects detargeting, exactly as upstream does. Shipped
targets are named instead of described, and ``--shipped-target`` is repeatable.

Usage:
    # child run: one campaign per target structure already in a table
    sapia run bindcraft2 outputs/RUN --table table0 -i pdb_path \\
        --hotspots 'A54,A56,A66,A115' --binder-lengths 60-100 --num-designs 10

    # root run: one campaign against a target not in any table yet
    sapia run bindcraft2 outputs/RUN \\
        --target-pdb targets/PDL1.pdb --chains A \\
        --hotspots 'A54,A56' --modality VHH --property humanize

    # cross-reactive against two orthologs, detargeting a third protein
    sapia run bindcraft2 outputs/RUN \\
        --target 'name=hPDL1;path=t/hPDL1.pdb;chains=A;hotspots=A54,A56' \\
        --target 'name=mPDL1;path=t/mPDL1.pdb;chains=A;hotspots=A36,A38' \\
        --target 'name=hPD1;path=t/hPD1.pdb;objective=detarget;weight=-0.5' \\
        --num-designs 10

    # the table's target, plus an off-target to avoid
    sapia run bindcraft2 outputs/RUN -t table0 -i pdb_path \\
        --extra-target 'name=hPD1;path=t/hPD1.pdb;weight=-0.5'

    # backbones only: stop before BindCraft2's own ProteinMPNN, to redesign the
    # sequences with this workspace's tools instead
    sapia run bindcraft2 outputs/RUN --target-pdb targets/PDL1.pdb \\
        --trajectory-only --max-trajectories 40 --hotspots 'A54,A56'

    # warm the 5.3 GB AlphaFold parameter Volume, running no campaign
    sapia run bindcraft2 outputs/RUN --fetch-weights-only --table-label weights
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

# Manifest modes, read by bindcraft2.sh as field 2.
MODE_DESIGN = "design"
MODE_FETCH_WEIGHTS = "fetch-weights"

# Design-property presets BindCraft2 ships under settings/property/, each turned on
# by its own ``--<name>`` flag on the bindcraft CLI. Kept as a list so --property
# validates against it instead of forwarding a typo that bindcraft would swallow as
# a stray path argument. A property this list does not know (a newer upstream) is
# still reachable as `--set <name>=true`.
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

# BindCraft2's own `targets[]` vocabulary (settings.py: TARGET_SETTING_NAMES). A key
# outside this set is refused here rather than by the campaign, so a typo costs a
# submit instead of a GPU hour.
TARGET_SETTING_NAMES = frozenset(
    {"name", "target_path", "chains", "hotspots", "coldspots", "weight", "objective"}
)
# Friendly spellings accepted on the command line for `target_path`.
TARGET_FIELD_ALIASES = {"path": "target_path", "pdb": "target_path"}
# Fields inside one --target spec are separated by ';' because hotspot and chain
# lists already use ',' ("hotspots=A54,A56").
TARGET_FIELD_SEP = ";"
# Fields whose value is a number rather than a string.
TARGET_NUMERIC_FIELDS = {"weight"}

# Target names the collector cannot use, because it derives `<metric>__<target>`
# columns alongside `<metric>__mean` and friends. A target called `mean` would land
# both in the same column and the summary would win -- one number quietly standing
# in for another. Keep in step with collect_bindcraft2._summaries.
RESERVED_TARGET_NAMES = frozenset({"mean", "worst", "best", "spread", "selectivity"})

# A {expr} placeholder island, resolved per-design up the table lineage. Meaningless
# in a root run (no table), so we reject it there.
_HAS_PLACEHOLDER = re.compile(r"\{[^}]*\}")


class SettingsConfigError(ValueError):
    """A run-wide settings misconfiguration that applies to every campaign (e.g. an
    --extra-settings key colliding with a dedicated flag). Unlike a per-row error it
    is not swallowed by the warn-and-skip loop -- it fails the whole submit up front."""


class BindCraft2Args(CommonArgs):
    fetch_weights_only: bool
    trajectory_only: bool
    reuse_campaigns: str | None
    targets: str | None
    target_pdb: Path | None
    shipped_target: list[str]
    target: list[str]
    extra_target: list[str]
    targets_file: str | None
    chains: str | None
    hotspots: str | None
    coldspots: str | None
    target_weight: float | None
    target_objective: str | None
    binder_lengths: str | None
    copies: int | None
    num_designs: int | None
    max_trajectories: int | None
    modality: str | None
    property: list[str]
    core: str | None
    campaign_seed: int | None
    metadata: str | None
    design_workers: int | None
    workers_per_gpu: str | None
    save_monomers: bool
    no_resume: bool
    extra_settings: str | None
    set: list[str]


def add_run_bindcraft2_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--fetch-weights-only",
        action="store_true",
        help="Run NO campaign: submit a single task that runs `bindcraft "
        "fetch-weights`, downloading the ~5.3 GB AlphaFold parameters into the cache "
        "Volume so later campaigns start warm. It needs no GPU and asks for none "
        "(gpus-per-task is forced to 0). Collect reports nothing, by design -- give "
        "the run its own --table-label so the empty table it reserves does not "
        "shadow a real one.",
    )
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

    targets = parser.add_argument_group(
        "targets",
        "What each campaign designs against. A campaign may carry several targets: "
        "orthologs to bind cross-reactively, and off-targets to avoid.",
    )
    targets.add_argument(
        "--target-pdb",
        type=Path,
        default=None,
        help="Single target structure to design binders against in a ROOT run (no "
        "--table): a PDB, mmCIF or FASTA (a FASTA target is treated as disordered "
        "and cropped; see `crop_fasta_sequence`) not in any table yet. Only valid "
        "without --table (with a table, targets come from --input-column). The "
        "design group is named `<stem>_bc2`. For several targets use --target.",
    )
    targets.add_argument(
        "--target",
        action="append",
        default=[],
        metavar="SPEC",
        help="One target, described in BindCraft2's own `targets[]` vocabulary as "
        "`;`-separated key=value fields (`;` and not `,`, because hotspot lists "
        "already use commas). Repeat for a multi-target campaign. Keys: name, path "
        "(= target_path), chains, hotspots, coldspots, weight, objective. `name` "
        "defaults to the file stem; a negative `weight` or `objective=detarget` "
        "makes it an off-target to avoid. Values may embed {expr}. Example: "
        "--target 'name=hPD1;path=t/hPD1.pdb;chains=A;objective=detarget;weight=-0.5'",
    )
    targets.add_argument(
        "--extra-target",
        action="append",
        default=[],
        metavar="SPEC",
        help="An ADDITIONAL target appended after the primary one (the table row's "
        "structure, or --target-pdb). Same syntax as --target. This is the flag for "
        "'design against the target in my table, while avoiding this off-target'.",
    )
    targets.add_argument(
        "--shipped-target",
        action="append",
        default=[],
        metavar="NAME",
        help="Name of a target BindCraft2 ships (hPDL1, hPD1, mPDL1, hIL2R, hIL7RA, "
        "dynorphin_a; `bindcraft design --list-targets` is authoritative). ROOT runs "
        "only. Repeatable: several accumulate into one cross-reactive campaign, as "
        "`\"target\": [...]` does upstream. The design group is named after them. "
        "Cannot be mixed with a described target -- upstream would silently drop the "
        "shipped ones (verified), so this tool refuses it instead.",
    )
    targets.add_argument(
        "--targets-file",
        type=str,
        default=None,
        metavar="FILE",
        help="YAML or JSON holding the whole `targets` list verbatim (either a bare "
        "list of target objects, or a mapping with a `targets:` key). The escape "
        "hatch when a campaign's targets are easier to keep in a file than on the "
        "command line. String values may embed {expr}.",
    )
    targets.add_argument(
        "--chains",
        type=str,
        default=None,
        help="Target chains to design against (per-target `chains`, e.g. 'A' or "
        "'A,B'), applied to the PRIMARY target (the table row, or --target-pdb). "
        "Omitted by default (BindCraft2 uses every chain in the file).",
    )
    targets.add_argument(
        "--hotspots",
        type=str,
        default=None,
        help="Target residues the binder should contact on the PRIMARY target "
        "(per-target `hotspots`), in BindCraft2's own syntax: comma-separated "
        "residues and ranges, chain-prefixed (e.g. 'A54,A56,A66-70'). May embed "
        "{expr} placeholders resolved per-design up the lineage. Omitted by default "
        "(BindCraft2 picks the epitope itself).",
    )
    targets.add_argument(
        "--coldspots",
        type=str,
        default=None,
        help="Regions of the PRIMARY target to avoid contacting (per-target "
        "`coldspots`), same syntax as --hotspots. Omitted by default.",
    )
    targets.add_argument(
        "--target-weight",
        type=float,
        default=None,
        help="Relative importance of the PRIMARY target (`targets[].weight`, default "
        "1). Only meaningful alongside --target/--extra-target; a negative value "
        "would make the table's own target an off-target, which is almost certainly "
        "a mistake.",
    )
    targets.add_argument(
        "--target-objective",
        type=str,
        default=None,
        choices=("target", "detarget"),
        help="Whether the PRIMARY target is to be bound or avoided "
        "(`targets[].objective`, default `target`).",
    )

    campaign = parser.add_argument_group("campaign", "Size, budget and format.")
    campaign.add_argument(
        "--binder-lengths",
        type=str,
        default=None,
        help="Binder size (`binder_lengths`): 'N' for one length, 'min-max' for a "
        "range drawn from per trajectory, or 'a,b,c' for a discrete choice (e.g. "
        "'80', '60-100', '60,80,100'). May embed {expr} placeholders. Omitted by "
        "default (BindCraft2's own default, or the modality preset's).",
    )
    campaign.add_argument(
        "--num-designs",
        type=int,
        default=None,
        help="Accepted designs to stop the campaign at (`number_of_final_designs`). "
        "This is the number of CHILD ROWS a campaign aims to produce, not a batch "
        "size -- BindCraft2 keeps spending trajectories until it has them. Omitted "
        "by default.",
    )
    campaign.add_argument(
        "--max-trajectories",
        type=int,
        default=None,
        help="Design attempts to spend before giving up (`max_trajectories`). The "
        "real cost knob: a campaign runs until --num-designs are accepted OR this "
        "many attempts are spent. Omitted by default.",
    )
    campaign.add_argument(
        "--modality",
        type=str,
        default=None,
        help="Binder format (-> `bindcraft design --modality`): binder, VHH, "
        "peptide, cyclic_peptide, ARP, scFv, Fab, large_binder, homo_oligomer, "
        "multidomain, induced_fit, fold_switch. Comma-separated to combine. "
        "Defaults to BindCraft2's own default (binder).",
    )
    campaign.add_argument(
        "--property",
        action="append",
        default=[],
        choices=PROPERTY_PRESETS,
        metavar="NAME",
        help="Design-property preset to switch on (-> `bindcraft design --<name>`). "
        "Repeatable. One of: " + ", ".join(PROPERTY_PRESETS) + ". A property this "
        "list does not know is still reachable as `--set <name>=true`.",
    )
    campaign.add_argument(
        "--core",
        type=str,
        default=None,
        help="Core profile applied under every preset (-> `bindcraft design --core`), "
        "e.g. 'benchmark' for a reproducible run. Comma-separated to combine. "
        "Omitted by default.",
    )
    campaign.add_argument(
        "--campaign-seed",
        type=int,
        default=None,
        help="Seed every trajectory is drawn from (`campaign_seed`). Set it with "
        "--core benchmark for a reproducible campaign. Omitted by default.",
    )
    campaign.add_argument(
        "--metadata",
        type=str,
        default=None,
        metavar="FILE",
        help="JSON object of descriptive fields (author, project, note) recorded "
        "with the campaign and written into its tables as `meta_<name>` columns "
        "(-> `bindcraft design --metadata`). Provenance, not settings: it changes "
        "nothing about the design.",
    )
    campaign.add_argument(
        "--save-monomers",
        action="store_true",
        help="Also keep the binder re-predicted ALONE beside each complex "
        "(`save_binder_monomers`), which the collector then records as "
        "`<leaf>_monomer_path`. The free binder is what a self-consistency check "
        "wants, so keeping it here saves a later `chainsel` step.",
    )
    campaign.add_argument(
        "--no-resume",
        action="store_true",
        help="Start each campaign from scratch instead of carrying on into a "
        "project_folder already written. By default this tool sets `resume: true`, "
        "so re-running a submit continues the campaigns it already started (the "
        "framework's own resume filter can't help here -- a create tool's status "
        "column lives in the child table, not the one it reads).",
    )

    execution = parser.add_argument_group(
        "execution", "How one campaign uses the hardware it was given."
    )
    execution.add_argument(
        "--design-workers",
        type=int,
        default=None,
        help="Parallel design workers for ONE campaign (`design_workers`). A "
        "campaign is a long serial loop by default; with several GPUs on the task "
        "(`-g N`) this fans its trajectories across them. Leave unset for one "
        "worker. Raising it without raising -g just contends for the same card.",
    )
    execution.add_argument(
        "--workers-per-gpu",
        type=str,
        default=None,
        metavar="N|auto",
        help="Workers packed onto each GPU (`workers_per_gpu`, default `auto`). "
        "More than one only helps when a single trajectory leaves the card idle.",
    )

    advanced = parser.add_argument_group(
        "advanced", "Everything BindCraft2 exposes that has no dedicated flag."
    )
    advanced.add_argument(
        "--extra-settings",
        type=str,
        default=None,
        help="Path to a YAML or JSON file: a mapping of (extra) BindCraft2 campaign "
        "settings merged into EVERY campaign's settings file (e.g. objective, "
        "aa_bias, min_iptm_final, max_detarget_iptm_final, multitarget_steps, "
        "save_design_trajectory). String values may embed {expr} placeholders "
        "resolved per-design up the lineage. "
        "`bindcraft design --list-settings` names every setting it accepts.",
    )
    advanced.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Extra `bindcraft design --set` override, passed through verbatim and "
        "applied over the generated settings file. Repeatable. Escape hatch for "
        "settings without a dedicated flag; spaces and JSON are fine (each token is "
        "carried as its own argv entry), e.g. --set 'binder_lengths=[70, 90]'.",
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


def load_targets_file(targets_file: str | None) -> list[dict[str, Any]]:
    """Parse --targets-file into a list of BindCraft2 target objects.

    Accepts either a bare list of target objects or a mapping carrying a ``targets``
    key, so the same file works as a --targets-file and as the targets block of an
    --extra-settings file.
    """
    if not targets_file:
        return []
    path = Path(targets_file)
    if not path.is_file():
        raise FileNotFoundError(f"--targets-file not found: {targets_file}")
    data = yaml.safe_load(path.read_text())
    if isinstance(data, dict):
        data = data.get("targets")
    if not isinstance(data, list) or not data:
        raise ValueError(
            f"--targets-file {targets_file} must hold a non-empty list of target "
            f"objects, or a mapping with a 'targets:' list"
        )
    entries = []
    for position, entry in enumerate(data, start=1):
        if not isinstance(entry, dict):
            raise ValueError(
                f"--targets-file {targets_file}: entry {position} is a "
                f"{type(entry).__name__}, expected a mapping of target fields"
            )
        entries.append(_validated_target(entry, f"--targets-file entry {position}"))
    return entries


def _validated_target(entry: dict[str, Any], source: str) -> dict[str, Any]:
    """Check one target object against BindCraft2's own `targets[]` vocabulary.

    Upstream rejects an unknown key too, but only once the campaign starts; catching
    it here costs a submit instead of a queued GPU task.
    """
    unknown = sorted(set(entry) - TARGET_SETTING_NAMES)
    if unknown:
        raise SettingsConfigError(
            f"{source}: unknown target field(s) {unknown}; BindCraft2 accepts "
            f"{sorted(TARGET_SETTING_NAMES)}"
        )
    if not entry.get("name"):
        raise SettingsConfigError(
            f"{source}: a target needs a `name` -- it identifies the target in the "
            f"output tables and in every per-target structure filename"
        )
    if str(entry["name"]) in RESERVED_TARGET_NAMES:
        raise SettingsConfigError(
            f"{source}: target name {entry['name']!r} is reserved. The collector "
            f"writes `<metric>__<target>` beside `<metric>__mean`/`__worst`/"
            f"`__best`/`__spread`/`__selectivity`, so this name would collide with a "
            f"summary column and one number would silently stand in for another. "
            f"Pick another name."
        )
    return entry


def parse_target_spec(spec: str, source: str) -> dict[str, Any]:
    """``'name=hPD1;path=t/hPD1.pdb;weight=-0.5'`` -> a BindCraft2 target object.

    Fields are ``;``-separated because hotspot and chain lists already use ``,``.
    ``path``/``pdb`` are accepted as friendlier spellings of ``target_path``, and
    ``name`` defaults to the structure's file stem. Placeholders are NOT resolved
    here -- that happens per design row, in _resolve_target.
    """
    entry: dict[str, Any] = {}
    for field in spec.split(TARGET_FIELD_SEP):
        field = field.strip()
        if not field:
            continue
        key, sep, value = field.partition("=")
        key = key.strip()
        if not sep:
            raise SettingsConfigError(
                f"{source}: field {field!r} is not `key=value`. A target reads "
                f"`name=X{TARGET_FIELD_SEP}path=Y{TARGET_FIELD_SEP}chains=A`, with "
                f"`{TARGET_FIELD_SEP}` between fields because hotspots already use "
                f"commas."
            )
        key = TARGET_FIELD_ALIASES.get(key, key)
        if key in entry:
            raise SettingsConfigError(f"{source}: field {key!r} given twice")
        entry[key] = value.strip()
    if not entry:
        raise SettingsConfigError(f"{source}: empty target spec")
    if "name" not in entry and entry.get("target_path"):
        entry["name"] = Path(str(entry["target_path"])).stem
    return _validated_target(entry, source)


def _coerce_target_numbers(entry: dict[str, Any], source: str) -> dict[str, Any]:
    """Turn the numeric target fields into numbers after placeholder resolution.

    Everything arrives from the command line as a string; `weight` has to reach the
    settings file as a number or BindCraft2's ``float(target.get('weight', 1.0))``
    is the first thing that sees the mistake.
    """
    coerced = dict(entry)
    for field in TARGET_NUMERIC_FIELDS & set(coerced):
        value = coerced[field]
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            continue
        try:
            coerced[field] = float(str(value))
        except ValueError as e:
            raise SettingsConfigError(
                f"{source}: {field}={value!r} is not a number"
            ) from e
    return coerced


def _resolve_target(
    entry: dict[str, Any], lookup: LookupFn, design: str, source: str
) -> dict[str, Any]:
    """Resolve {expr} placeholders in one target and make its path absolute."""
    resolved = cast(dict[str, Any], _resolve_tree(entry, lookup, design))
    if resolved.get("target_path"):
        path = volume_path(str(resolved["target_path"]))
        if not path.exists():
            raise SettingsConfigError(f"{source}: target_path {path} does not exist")
        resolved["target_path"] = str(path)
    return _coerce_target_numbers(resolved, source)


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
    """'80' -> [80]; '60-100' / '60,100' -> [60, 100]; '60,80,100' -> discrete set.

    BindCraft2 reads [min, max] for two entries and a discrete choice for three or
    more, so the separator is preserved as given rather than interpreted here.
    """
    parts = [p.strip() for p in re.split(r"[-,]", spec) if p.strip()]
    try:
        return [int(p) for p in parts]
    except ValueError as e:
        raise ValueError(
            f"--binder-lengths {spec!r} is not a length, a 'min-max' range or a "
            f"comma-separated set of lengths"
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


# Every flag that describes a campaign. A --reuse-campaigns run submits nothing, so
# any of them is a mistake; the list is here (rather than inline) so a new flag has
# one obvious place to be registered.
def _campaign_flag_values(args: BindCraft2Args) -> list[tuple[str, Any]]:
    return [
        ("--target-pdb", args.target_pdb),
        ("--target", args.target),
        ("--extra-target", args.extra_target),
        ("--shipped-target", args.shipped_target),
        ("--targets-file", args.targets_file),
        ("--trajectory-only", args.trajectory_only),
        ("--fetch-weights-only", args.fetch_weights_only),
        ("--hotspots", args.hotspots),
        ("--coldspots", args.coldspots),
        ("--chains", args.chains),
        ("--target-weight", args.target_weight),
        ("--target-objective", args.target_objective),
        ("--binder-lengths", args.binder_lengths),
        ("--num-designs", args.num_designs),
        ("--max-trajectories", args.max_trajectories),
        ("--modality", args.modality),
        ("--core", args.core),
        ("--campaign-seed", args.campaign_seed),
        ("--metadata", args.metadata),
        ("--design-workers", args.design_workers),
        ("--workers-per-gpu", args.workers_per_gpu),
        ("--save-monomers", args.save_monomers),
        ("--extra-settings", args.extra_settings),
        ("--property", args.property),
        ("--set", args.set),
    ]


def _reject_campaign_flags(args: BindCraft2Args) -> None:
    """A reuse run submits nothing, so anything describing a campaign is a mistake."""
    named = [flag for flag, value in _campaign_flag_values(args) if value]
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


def _primary_target(args: BindCraft2Args, name: str, target_path: Path) -> dict:
    """The target a child row (or --target-pdb) contributes, as a target object."""
    entry: dict[str, Any] = {"name": name, "target_path": str(target_path)}
    for field, value in (
        ("chains", args.chains),
        ("hotspots", args.hotspots),
        ("coldspots", args.coldspots),
        ("objective", args.target_objective),
    ):
        if value is not None:
            entry[field] = value
    if args.target_weight is not None:
        entry["weight"] = args.target_weight
    return entry


def campaign_targets(
    args: BindCraft2Args, name: str, target_path: Path | None
) -> tuple[list[dict] | None, list[str] | None]:
    """The target entries, and the shipped target names, this campaign declares.

    Returns ``(targets, shipped)``, either of which may be None when this tool is not
    the one deciding it -- that is what lets --extra-settings own `targets` outright
    without tripping the collision check.

    Refuses the one combination upstream accepts and then silently ruins: a shipped
    target beside an explicit `targets` list. ``campaign_over_presets`` layers the
    accumulated preset targets *under* the request, and a list replaces rather than
    merges, so the shipped entries vanish without a word (verified against v1.0.3).
    """
    described: list[dict] = []
    if target_path is not None:
        described.append(_primary_target(args, name, target_path))
    for position, spec in enumerate(args.target, start=1):
        described.append(parse_target_spec(spec, f"--target #{position}"))
    for position, spec in enumerate(args.extra_target, start=1):
        described.append(parse_target_spec(spec, f"--extra-target #{position}"))
    described += load_targets_file(args.targets_file)

    shipped = list(args.shipped_target)
    if shipped and described:
        raise SettingsConfigError(
            f"--shipped-target {shipped} cannot be combined with a described target "
            f"({', '.join(sorted({str(t.get('name')) for t in described}))}). "
            f"BindCraft2 lets an explicit `targets` list REPLACE the shipped presets "
            f"rather than extend them, so the shipped targets would be dropped "
            f"silently. Write the shipped target out as its own --target (its "
            f"structure and hotspots are in settings/target/<name>.json), or use "
            f"--shipped-target alone."
        )

    if described:
        names = [str(entry.get("name")) for entry in described]
        duplicated = sorted({n for n in names if names.count(n) > 1})
        if duplicated:
            raise SettingsConfigError(
                f"target name(s) {duplicated} used more than once. Names identify a "
                f"target in the output tables and in every per-target structure "
                f"filename, so they have to be distinct."
            )
        return described, None
    if shipped:
        return None, shipped
    # Neither: --extra-settings may still own `targets`, checked after the merge.
    return None, None


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

    ``target_path`` is the (absolute) primary target structure, or ``None`` when the
    targets come from --target/--shipped-target/--targets-file/--extra-settings.
    Raises ValueError on a per-row problem (an unresolvable {expr}) and the run-wide
    ``SettingsConfigError`` on a misconfiguration that would hit every row.
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

    targets, shipped = campaign_targets(args, name, target_path)
    if targets is not None:
        tool_fields["targets"] = [
            _resolve_target(entry, lookup, name, f"target {entry.get('name')!r}")
            for entry in targets
        ]
    if shipped is not None:
        tool_fields["target"] = shipped
        # A shipped target is named rather than described: its own preset carries the
        # path, chains and hotspots, so per-target overrides would have nowhere to
        # land.
        for flag, value in (
            ("--chains", args.chains),
            ("--hotspots", args.hotspots),
            ("--coldspots", args.coldspots),
            ("--target-weight", args.target_weight),
            ("--target-objective", args.target_objective),
        ):
            if value is not None:
                raise SettingsConfigError(
                    f"{flag} describes a target file and cannot be combined with "
                    f"--shipped-target (its preset already names the epitope). Pass "
                    f"the structure with --target-pdb or --target, or override the "
                    f"preset with --extra-settings."
                )

    if args.binder_lengths is not None:
        tool_fields["binder_lengths"] = _coerce_binder_lengths(
            resolve_template(args.binder_lengths, lookup, name)
        )
    if args.copies is not None:
        if args.copies < 1:
            raise SettingsConfigError(
                f"--copies {args.copies} is not a chain count; it is the number of "
                f"identical chains the binder is built from, so it must be 1 or more "
                f"(1 = an ordinary single-chain binder)."
            )
        tool_fields["copies"] = args.copies
    if args.num_designs is not None:
        tool_fields["number_of_final_designs"] = args.num_designs
    if args.max_trajectories is not None:
        tool_fields["max_trajectories"] = args.max_trajectories
    if args.campaign_seed is not None:
        tool_fields["campaign_seed"] = args.campaign_seed
    if args.design_workers is not None:
        tool_fields["design_workers"] = args.design_workers
    if args.workers_per_gpu is not None:
        tool_fields["workers_per_gpu"] = args.workers_per_gpu
    if args.save_monomers:
        tool_fields["save_binder_monomers"] = True

    resolved_extra = _resolve_tree(extra_fields, lookup, name)
    collisions = sorted(set(tool_fields) & set(resolved_extra))
    if collisions:
        raise SettingsConfigError(
            f"setting(s) {collisions} set by both a dedicated flag and "
            "--extra-settings; remove them from one source"
        )

    settings = {**defaults, **resolved_extra, **tool_fields}
    if not settings.get("targets") and not settings.get("target"):
        raise SettingsConfigError(
            "this campaign declares no target. Give it one with --target-pdb, "
            "--target, --shipped-target or --targets-file, run with --table so the "
            "table's --input-column supplies it, or put a `targets:` list in "
            "--extra-settings."
        )
    return settings


def _cli_tokens(args: BindCraft2Args) -> list[str]:
    """The run-wide ``bindcraft design`` argv tokens (identical for every campaign).

    Presets and --set live on the command line rather than in the settings file
    because that is the interface BindCraft2 documents for them: --modality/--core
    name preset files to layer under the campaign, and --set is applied over it.

    Returned as a LIST of argv entries, not a joined string: each is written to the
    campaign's args file on its own line and read back with `mapfile`, so a value
    carrying spaces or JSON survives intact.
    """
    tokens: list[str] = []
    if args.core is not None:
        tokens += ["--core", args.core]
    if args.modality is not None:
        tokens += ["--modality", args.modality]
    for prop in args.property:
        tokens.append("--" + prop.replace("_", "-"))
    if args.metadata is not None:
        metadata = volume_path(args.metadata)
        if not metadata.is_file():
            raise SettingsConfigError(f"--metadata file not found: {args.metadata}")
        tokens += ["--metadata", str(metadata)]
    for assignment in args.set:
        tokens += ["--set", assignment]

    # The args file is one token per line, so a token containing a newline would be
    # read back as two. Nothing BindCraft2 accepts needs one.
    for token in tokens:
        if "\n" in token:
            raise SettingsConfigError(
                f"argument {token!r} contains a newline, which the per-campaign args "
                f"file cannot carry. Put the value in --extra-settings instead."
            )
    return tokens


def _write_args_file(path: Path, tokens: list[str]) -> Path:
    """One argv token per line, for `mapfile -t` in the task script."""
    path.write_text("".join(f"{token}\n" for token in tokens))
    return path


def _target_groups(ctx: ManifestCtx[BindCraft2Args]) -> list[tuple[str, Path | None]]:
    """The (name, primary target path) groups this run designs against: one per ready
    table row for a child run, a single group for a root run."""
    args = ctx.args

    if args.table is not None:
        if args.target_pdb is not None or args.shipped_target:
            raise ValueError(
                "--target-pdb/--shipped-target are only valid for a root run (no "
                "--table); with --table, the primary target comes from the table's "
                "--input-column. Add further targets with --extra-target."
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
    placeheld = [
        flag
        for flag, value in (
            ("--hotspots", args.hotspots),
            ("--coldspots", args.coldspots),
            ("--chains", args.chains),
            ("--binder-lengths", args.binder_lengths),
            *((f"--target {s!r}", s) for s in args.target),
            *((f"--extra-target {s!r}", s) for s in args.extra_target),
        )
        if value and _HAS_PLACEHOLDER.search(str(value))
    ]
    if placeheld:
        raise ValueError(
            f"{', '.join(placeheld)} contain a {{expr}} placeholder, but this is a "
            f"root run (no --table) with no table lineage to resolve it against. Use "
            f"literal values, or run with --table."
        )
    if targets_spec is not None:
        if _tree_has_placeholder(targets_spec.entries):
            raise ValueError(
                "--targets contains a {expr} placeholder, but this is a root run (no "
                "--table) with no table lineage to resolve it against. Use literal "
                "values, or run with --table."
            )
        return [(f"{targets_spec.stem}_bc2", None)]

    if args.target_pdb is not None and args.shipped_target:
        raise ValueError(
            "--target-pdb and --shipped-target both name a target; pass exactly one. "
            "To combine several described targets, repeat --target instead."
        )
    if args.target_pdb is not None:
        target_path = volume_path(args.target_pdb)
        if not target_path.exists():
            raise FileNotFoundError(f"--target-pdb {target_path} does not exist.")
        return [(f"{target_path.stem}_bc2", target_path)]
    if args.shipped_target:
        return [("_".join(args.shipped_target) + "_bc2", None)]
    if args.target or args.targets_file:
        # A described multi-target root campaign: name the group after its binding
        # targets so the campaign folder says what it was designed against.
        described, _ = campaign_targets(args, "", None)
        binding = [
            str(entry["name"])
            for entry in (described or [])
            if float(entry.get("weight", 1.0)) > 0
            and entry.get("objective") != "detarget"
        ]
        return [("_".join(binding or ["targets"]) + "_bc2", None)]

    raise ValueError(
        "A root run (no --table) needs a target: pass --target-pdb <file>, --target "
        "<spec>, --shipped-target <name> or --targets-file <file>. With --table, the "
        "primary target comes from --input-column."
    )


def _fetch_weights_only(ctx: ManifestCtx[BindCraft2Args]) -> list[tuple[str, ...]]:
    """Submit one task that downloads the AlphaFold parameters, and no campaign.

    The parameters are ~5.3 GB and land in the cache Volume on first use. Warming
    them once beats N cold containers each pulling the same archive, and this is the
    supported way to do it -- `bindcraft fetch-weights` verifies the checkpoints and
    exits non-zero if they are missing or unfinished, which a dummy campaign does not.
    """
    named = [flag for flag, value in _campaign_flag_values(ctx.args) if value]
    named = [flag for flag in named if flag != "--fetch-weights-only"]
    if named:
        raise ValueError(
            f"--fetch-weights-only runs no campaign, so {', '.join(named)} would "
            f"have no effect. Submit the warm-up on its own, then the campaign."
        )
    ctx.args.gpus_per_task = 0
    ctx.write_meta(weights_only=True, root_designs=[])
    print("Fetching the AlphaFold parameters into the cache volume; no campaign.")
    print("Nothing will be collected from this run.")
    return [("bindcraft2_weights", MODE_FETCH_WEIGHTS, "-", "-")]


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
    if ctx.args.fetch_weights_only:
        return _fetch_weights_only(ctx)
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

    # Run-wide, so a bad --metadata path or an unquotable token fails once, here,
    # rather than once per campaign.
    args_file = _write_args_file(
        settings_dir / "design.args", _cli_tokens(ctx.args)
    )

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
        manifest_rows.append(
            (
                name,
                MODE_DESIGN,
                str(volume_path(settings_json)),
                str(volume_path(args_file)),
            )
        )
        submitted.append(name)

    if ctx.args.table is None:
        # A root run has no parent table for collect to iterate; record the group
        # names so `sapia collect` can find their campaigns and rebuild their rows.
        ctx.write_meta(root_designs=submitted)

    return manifest_rows
