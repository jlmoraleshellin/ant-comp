"""Submit AtomiUM's ``design_ph.py`` on binder COMPLEXES, one child row per sequence.

``design_ph.py`` designs pH-selective binders: it pins His/Asp/Glu centres into a
protonation state and redesigns around them so that state is favoured in the complex
more than in the free binder (``ddE = dE(complex) - dE(apo binder)``). The input must
be a complex, so ``--binder-chains`` and ``--target-chains`` are required and disjoint,
``--center-mode`` has no default, and ``--lam`` must contain 0.0 (the control cannot be
recovered without a GPU rerun). Premise, columns and traps:
``.claude/skills/atomiumph/SKILL.md``.

``action="create"``: rows are ``<parent>_p<i>``, one per (lambda, sample) draw, so a
structure yields ``len(--lam) x --num-seq`` rows.

``--binder-alone`` runs the binder with no target in the graph, to get the binder's own
``core`` (design_ph counts burial over every chain present). ddE is then undefined (NA).
Refused in this mode: a target chain, any ``--lam`` but 0.0, ``--placement-region
interface``, ``--set`` of ``--placement_metric``/``--coupling_metric`` ``ddE``, a staged
structure holding any chain but the binder, and flipping the flag in a leaf that holds
output. ``--apo-weight`` and ``--charge-prior`` are sent, and recorded, as 0.0.

Re-submits: ``_check_no_stale_outputs`` refuses a submit that would leave old output in
the leaf, so retry failed parents with the whole same filter under the same label. The
driver writes the sidecar before this builder runs, so even a refused submit moves its
mtime (``atomiumph_output_predates_meta``); ``ctx.write_meta`` is therefore the
builder's last statement.

Batching mirrors ``tools/atomium``: one ``design_ph.py --pdb_dir`` call per subgroup of
symlinks ``grp_<g>/inputs/<name>.pdb``, bin-packed up to ``--designs-per-task``. The
symlink stem is the parent row name, which is the lineage contract (see the collector).
"""

import json
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path
from typing import cast

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.utils import ensure_pdb, expand_chain_spec

# The four 30-class pH checkpoints design_ph.py accepts (atomium's are 21-class).
CHECKPOINT_NOISES = ("n03", "n05", "n03_v2", "n05_v2")

CENTER_MODES = ("his_p", "his_s", "acid_p", "all_p")
PLACEMENT_REGIONS = ("interface", "core", "surface", "any")
REDESIGN_SCOPES = ("neighbourhood", "binder", "coupling")
COUPLING_PROBES = ("ala", "mask")
SEED_MODES = ("native", "generate")
OUTPUT_SCOPES = ("binder", "all")

# Sidecar keys read back by the collector as run-wide trust columns.
EXPECTED_KEY = "atomiumph_expected_per_structure"
BINDER_CHAINS_KEY = "atomiumph_binder_chains"
TARGET_CHAINS_KEY = "atomiumph_target_chains"
# Run-wide, and decisive for reading rows: center_mode sets which pH binds.
CENTER_MODE_KEY = "atomiumph_center_mode"
OUTPUT_SCOPE_KEY = "atomiumph_output_scope"
# The noise level decides WHICH centres are picked, so levels are not comparable.
MODEL_NOISE_KEY = "atomiumph_model_noise"
LAM_SWEEP_KEY = "atomiumph_lam_sweep"
# Changes which positions may move, so runs differing in it are not comparable.
PROTECT_TITRATABLE_KEY = "atomiumph_protect_titratable"
# Modelled because the upstream client default (0.5) differs from design_ph's (0.0).
APO_WEIGHT_KEY = "atomiumph_apo_weight"
# What design_ph was GIVEN (--binder-alone forces 0.0), not what the caller typed.
CHARGE_PRIOR_KEY = "atomiumph_charge_prior"
# --binder-alone: no target in the graph; the collector NAs target-dependent columns.
BINDER_ALONE_KEY = "atomiumph_binder_alone"

# Flags emitted here, mapped to their typed flag; argparse takes the LAST occurrence.
MANAGED_FLAGS = {
    "--binder_chains": "--binder-chains",
    "--target_chains": "--target-chains",
    "--lam": "--lam",
    "--center_mode": "--center-mode",
    "--centers_per_type": "--centers-per-type",
    "--placement_region": "--placement-region",
    "--checkpoint_model_noise": "--model-noise",
    "--protect_titratable": "--protect-titratable",
    "--apo_weight": "--apo-weight",
    "--temperature": "--temperature",
    "--num_seq": "--num-seq",
    "--seed_mode": "--seed-mode",
    "--seed_temperature": "--seed-temperature",
    "--seed": "--seed",
    "--redesign_scope": "--redesign-scope",
    "--neighbour_k": "--neighbour-k",
    "--max_mutations": "--max-mutations",
    "--coupling_radius": "--coupling-radius",
    "--coupling_probe": "--coupling-probe",
    "--ban_cys": "--ban-cys",
    "--charge_prior": "--charge-prior",
    "--charge_cutoff": "--charge-cutoff",
    "--interface_cutoff": "--interface-cutoff",
    "--burial_cutoff": "--burial-cutoff",
    "--output_scope": "--output-scope",
    # Emitted by atomiumph.sh, not by the tail, but just as fatal to override.
    "--pdb_dir": "(the staged inputs dir)",
    "--pdb_path": "(the staged inputs dir)",
    "--pattern": "(the staged inputs dir)",
    "--out_folder": "(the per-subgroup output dir)",
}


class AtomiumphArgs(CommonArgs):
    binder_chains: str
    target_chains: str
    binder_alone: bool
    center_mode: str
    centers_per_type: int
    placement_region: str
    model_noise: str
    lam: list[float]
    no_control: bool
    temperature: float
    num_seq: int
    seed_mode: str
    seed_temperature: float
    seed: int
    redesign_scope: str
    neighbour_k: int
    max_mutations: int
    coupling_radius: float
    coupling_probe: str
    ban_cys: int
    protect_titratable: int
    apo_weight: float
    charge_prior: float
    charge_cutoff: float
    interface_cutoff: float
    burial_cutoff: int
    output_scope: str
    designs_per_task: int
    set: list[str]


def _chain_list(spec: str, flag: str) -> list[str]:
    """Expand the chain mini-language (``A:C,E`` -> ``[A, B, C, E]``) for one flag."""
    try:
        return expand_chain_spec(spec)
    except ValueError as exc:
        raise ValueError(f"{flag}: {exc}")


def _check_chains(args: AtomiumphArgs) -> tuple[str, str]:
    """Binder/target chain lists as space-separated strings, validated up front.

    Both are required (upstream's binder default ``A`` is wrong for rfd3 complexes,
    where the binder is usually B). ``--binder-alone`` requires an empty target list.
    """
    binder = _chain_list(args.binder_chains, "--binder-chains")
    target = _chain_list(args.target_chains, "--target-chains")
    if not binder:
        raise ValueError(
            "--binder-chains is required: name the chain(s) to redesign, e.g. "
            "'B' (rfd3 binder complexes in this workspace are usually target=A, "
            "binder=B -- check the structure, do not assume)"
        )
    if args.binder_alone:
        if target:
            raise ValueError(
                f"--binder-alone and --target-chains {args.target_chains!r} contradict "
                f"each other: the mode exists to keep the target OUT of the graph, so "
                f"that 'core' (a Cb-neighbour count over every chain present) is the "
                f"binder's own core. Drop --target-chains, or drop --binder-alone and "
                f"run the complex."
            )
        return " ".join(binder), ""
    if not target:
        raise ValueError(
            "--target-chains is required: design_ph scores the protonation centre in "
            "the COMPLEX against the APO binder, so without a target chain every ddE "
            "is structurally zero and the run measures nothing"
        )
    overlap = [c for c in binder if c in target]
    if overlap:
        raise ValueError(
            f"--binder-chains and --target-chains overlap on {','.join(overlap)}: the "
            f"same chain cannot be both the redesigned binder and the target that the "
            f"apo reference removes"
        )
    return " ".join(binder), " ".join(target)


def _check_lam(
    lams: list[float], no_control: bool, binder_alone: bool = False
) -> list[str]:
    """``--lam`` must be a non-empty sweep with the 0.0 control; returns argv strings.

    Checks run on the rendered strings (``_fmt_lam`` rounds to 6 decimals), which is
    what design_ph receives. ``--binder-alone`` accepts only ``0.0``: with no apo
    graph every lambda arm would draw from the same distribution.
    """
    if not lams:
        raise ValueError("--lam: at least one lambda value is required")
    rendered = [_fmt_lam(x) for x in lams]
    if binder_alone:
        if rendered != ["0.0"]:
            raise ValueError(
                f"--binder-alone accepts only --lam 0.0, got {' '.join(rendered)}. "
                f"lambda scales the apo-vs-complex selective term and there is no apo "
                f"contrast without a target chain, so every arm would draw from the "
                f"same distribution and the extra rows would differ only in their "
                f"lambda column."
            )
        return rendered
    duplicates = sorted({x for x in rendered if rendered.count(x) > 1})
    if duplicates:
        raise ValueError(
            f"--lam: value(s) {', '.join(duplicates)} appear more than once once "
            f"rendered for argv (values are rounded to 6 decimals), which would produce "
            f"indistinguishable duplicate rows; list each lambda once"
        )
    if "0.0" in rendered or no_control:
        return rendered
    raise ValueError(
        "--lam must include 0.0, the internal control. lambda=0.0 switches off the "
        "selective term, so those designs are what this backbone produces WITHOUT "
        "pH selection pressure. Without them a good ddE cannot be distinguished from "
        "'this backbone was simply good', and the control cannot be recovered "
        "afterwards without re-running the whole sweep on GPU. Add 0.0, or pass "
        "--no-control to say deliberately that you already have it."
    )


def _check_choices(args: AtomiumphArgs) -> None:
    """Re-validate enumerated flags: argparse ``choices`` only guards the CLI."""
    for flag, value, allowed in (
        ("--model-noise", args.model_noise, CHECKPOINT_NOISES),
        ("--placement-region", args.placement_region, PLACEMENT_REGIONS),
        ("--redesign-scope", args.redesign_scope, REDESIGN_SCOPES),
        ("--coupling-probe", args.coupling_probe, COUPLING_PROBES),
        ("--seed-mode", args.seed_mode, SEED_MODES),
        ("--output-scope", args.output_scope, OUTPUT_SCOPES),
    ):
        if value not in allowed:
            extra = (
                " (these are the four fine-tuned 30-class pH checkpoints under "
                "ft_weights_for_pH/; tools/atomium's n00..n07 are the 21-class BASE "
                "weights and do not apply here. The _v2 pair are different weights, "
                "not newer builds of the other two)"
                if flag == "--model-noise"
                else ""
            )
            raise ValueError(
                f"{flag}: {value!r} is not one of {', '.join(allowed)}{extra}"
            )
    # Required, never defaulted: each value gives a healthy-looking, different molecule.
    if not args.center_mode:
        raise ValueError(
            "--center-mode is required: it decides which pH behaviour you are "
            "designing, and every value produces a full table that looks fine. "
            "his_s = binds at NEUTRAL pH, releases in acid (recycling / endosome "
            "release; TNF-alpha). his_p = binds in ACID, releases at neutral (EGFR). "
            "acid_p / all_p pin Asp/Glu, which only switch below ~5.8 -- diagnostic "
            "only, not usable for either His challenge. Upstream's acid_p default is "
            "deliberately not inherited."
        )
    if args.center_mode not in CENTER_MODES:
        raise ValueError(
            f"--center-mode: {args.center_mode!r} is not one of "
            f"{', '.join(CENTER_MODES)}"
        )
    if args.ban_cys not in (0, 1):
        raise ValueError(f"--ban-cys: expected 0 or 1, got {args.ban_cys}")
    if args.protect_titratable not in (0, 1):
        raise ValueError(
            f"--protect-titratable: expected 0 or 1, got {args.protect_titratable}"
        )
    if args.num_seq < 1:
        raise ValueError(f"--num-seq: must be at least 1, got {args.num_seq}")
    if args.max_mutations < 1:
        raise ValueError(f"--max-mutations: must be at least 1, got {args.max_mutations}")
    if args.centers_per_type < 1:
        raise ValueError(
            f"--centers-per-type: must be at least 1, got {args.centers_per_type}"
        )
    if args.neighbour_k < 1:
        raise ValueError(f"--neighbour-k: must be at least 1, got {args.neighbour_k}")
    # --binder-alone x interface provably fails (empty mask, design_ph.py:571-579).
    if args.binder_alone and args.placement_region == "interface":
        raise ValueError(
            "--binder-alone with --placement-region interface: there is no target to "
            "be at the interface of. region_masks returns an all-False interface mask "
            "without target chains and design_ph then raises \"Region 'interface' is "
            "empty\" for every structure. Use --placement-region core (the reason this "
            "mode exists), or surface/any."
        )


def _check_set(tokens: list[str], binder_alone: bool = False) -> None:
    """``--set`` is for unmodelled flags only and must survive a TSV round trip.

    Refused: a token repeating a flag this builder emits (argparse takes the last), one
    with a tab/newline (the sub-manifest is TSV), and under ``--binder-alone`` a
    ``--placement_metric``/``--coupling_metric`` of ``ddE`` (no apo graph, so upstream
    silently falls back to dE_complex).
    """
    for token in tokens:
        for bad, what in (("\t", "a tab"), ("\n", "a newline"), ("\r", "a carriage return")):
            if bad in token:
                raise ValueError(
                    f"--set: {token!r} contains {what}. The sub-manifest is TSV and the "
                    f"task script cuts the argv tail out of it, so that character would "
                    f"silently truncate or forge a row."
                )
        # Every word: a --set string may carry a flag and its value.
        for word in token.split():
            flag = word.split("=", 1)[0]
            if flag in MANAGED_FLAGS:
                owner = MANAGED_FLAGS[flag]
                raise ValueError(
                    f"--set: {flag} is set by this tool; passing it again (in "
                    f"{token!r}) would override the validated value (design_ph.py "
                    f"takes the last occurrence) without appearing in any column. "
                    f"Use {owner} instead."
                )

    if not binder_alone:
        return
    # Value as well as flag: `--placement_metric dE_complex` is honest and allowed.
    words: list[str] = []
    for token in tokens:
        words.extend(token.split())
    for i, word in enumerate(words):
        flag, _, inline = word.partition("=")
        if flag not in ("--placement_metric", "--coupling_metric"):
            continue
        value = inline if inline else (words[i + 1] if i + 1 < len(words) else "")
        if value == "ddE":
            raise ValueError(
                f"--set: {flag} ddE is impossible under --binder-alone. ddE needs the "
                f"apo graph, which design_ph builds only when there are target chains, "
                f"so the run would silently use dE_complex instead"
                + (
                    " -- and --coupling_metric reaches no column, so nothing in the "
                    "table would record that it did."
                    if flag == "--coupling_metric"
                    else "."
                )
            )


def _pdb_chain_counts(path: Path) -> dict[str, int]:
    """``{chain_id: n_atoms}`` over the ATOM/HETATM records of a PDB file.

    ``ensure_pdb`` has already converted anything else, so PDB suffices. Raises on a
    missing or empty file: a guard that cannot read its input must not pass.
    """
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} does not exist, so --binder-alone cannot check that it holds the "
            f"binder chain only"
        )
    counts: dict[str, int] = {}
    with open(path, errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
                chain = line[21].strip() or "_"
                counts[chain] = counts.get(chain, 0) + 1
    if not counts:
        raise ValueError(
            f"{path} holds no ATOM/HETATM records; --binder-alone cannot verify that "
            f"it is the binder chain alone"
        )
    return counts


def _check_binder_alone_input(name: str, path: Path, binder_chains: str) -> None:
    """Under ``--binder-alone``, refuse a staged structure holding any other chain.

    ``core`` counts Cb neighbours over every position in the graph (design_ph.py:582),
    so an extra chain, or the padding it may leave, inflates burial invisibly.
    """
    wanted = set(binder_chains.split())
    counts = _pdb_chain_counts(path)
    extra = sorted(c for c in counts if c not in wanted)
    missing = sorted(c for c in wanted if c not in counts)
    if extra:
        detail = ", ".join(f"{c} ({counts[c]} atoms)" for c in extra)
        raise ValueError(
            f"--binder-alone: {name} ({path}) contains chain(s) {detail} that are not "
            f"in --binder-chains '{binder_chains}'. This mode needs a structure holding "
            f"the binder chain ALONE: design_ph's 'core' is a Cb-neighbour count over "
            f"every position in the graph (region_masks), so anything else in the file "
            f"inflates the burial it measures, invisibly. Cut the binder out first "
            f"(chainsel) and point -i at that column."
        )
    if missing:
        raise ValueError(
            f"--binder-alone: {name} ({path}) has no atoms in chain(s) "
            f"{', '.join(missing)}, which --binder-chains names; the file holds "
            f"{', '.join(sorted(counts))}. Check --binder-chains against the structure."
        )


def _apo_weight_value(args: AtomiumphArgs) -> str:
    """``--apo_weight`` as design_ph receives it: 0.0 under ``--binder-alone``.

    One function for the argv tail and the recorded column, so they cannot drift.
    """
    return "0.0" if args.binder_alone else str(args.apo_weight)


def _charge_prior_value(args: AtomiumphArgs) -> str:
    """``--charge_prior`` as design_ph receives it: 0.0 under ``--binder-alone``.

    The prior sums TARGET charges, so with no target it is sent as zero.
    """
    return "0.0" if args.binder_alone else str(args.charge_prior)


def _design_ph_extra(args: AtomiumphArgs) -> str:
    """The run-wide ``design_ph.py`` argv tail, whitespace-free (expanded unquoted).

    Chains and ``--lam`` travel as their own manifest fields. Under ``--binder-alone``
    both target-dependent weights are sent as 0.0 instead of trusting upstream's guards.
    """
    parts = [
        "--center_mode", args.center_mode,
        "--centers_per_type", str(args.centers_per_type),
        "--placement_region", args.placement_region,
        "--checkpoint_model_noise", args.model_noise,
        "--temperature", str(args.temperature),
        "--num_seq", str(args.num_seq),
        "--seed_mode", args.seed_mode,
        "--seed_temperature", str(args.seed_temperature),
        "--seed", str(args.seed),
        "--redesign_scope", args.redesign_scope,
        "--neighbour_k", str(args.neighbour_k),
        "--max_mutations", str(args.max_mutations),
        "--coupling_radius", str(args.coupling_radius),
        "--coupling_probe", args.coupling_probe,
        "--ban_cys", str(args.ban_cys),
        "--protect_titratable", str(args.protect_titratable),
        "--apo_weight", _apo_weight_value(args),
        "--charge_prior", _charge_prior_value(args),
        "--charge_cutoff", str(args.charge_cutoff),
        "--interface_cutoff", str(args.interface_cutoff),
        "--burial_cutoff", str(args.burial_cutoff),
        "--output_scope", args.output_scope,
    ]
    parts.extend(args.set)
    return " ".join(parts)


# The signature is constant per run; kept so a per-design flag cannot mis-batch.
Member = tuple[str, Path]  # (design_name, staged_pdb_source)
Signature = tuple[str, str]  # (binder_chains, target_chains)
Subgroup = tuple[Signature, list[Member]]


def _pack_subgroups(chunks: list[Subgroup], max_per_task: int) -> list[list[Subgroup]]:
    """First-fit-pack subgroup chunks (each <= max_per_task) into tasks."""
    tasks: list[list[Subgroup]] = []
    remaining: list[int] = []
    for chunk in chunks:
        size = len(chunk[1])
        for t, cap in enumerate(remaining):
            if size <= cap:
                tasks[t].append(chunk)
                remaining[t] -= size
                break
        else:
            tasks.append([chunk])
            remaining.append(max_per_task - size)
    return tasks


# Suffixes design_ph.py writes per stem; mirrors the collector's SUFFIXES.
_OUTPUT_SUFFIXES = (".metrics.csv", ".couplings.csv", ".centers.json", ".ext.fa", ".fa")
# Dirs the run owns (no design_ph output); mirrors the collector's SKIP_DIRS.
_SKIP_DIRS = ("atomiumph_logs", "atomiumph_tasks", "inputs")


def _output_stems(out_dir: Path) -> set[str]:
    """Every design_ph output stem already in the leaf (one recursive walk).

    Shared by both submit-time guards; the run's own dirs are skipped.
    """
    stems: set[str] = set()
    if not out_dir.is_dir():
        return stems
    for path in out_dir.rglob("*"):
        if not path.is_file():
            continue
        if any(part in _SKIP_DIRS for part in path.relative_to(out_dir).parts[:-1]):
            continue
        for suffix in _OUTPUT_SUFFIXES:
            if path.name.endswith(suffix) and len(path.name) > len(suffix):
                stems.add(path.name[: -len(suffix)])
                break
    return stems


def _read_meta(out_dir: Path) -> dict:
    """The leaf's sidecar as a dict, or ``{}`` if absent, unreadable or not an object.

    Read before this submit's ``ctx.write_meta``, so it holds the previous run's values.
    """
    path = out_dir / ".meta.json"
    if not path.is_file():
        return {}
    try:
        meta = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return meta if isinstance(meta, dict) else {}


def _check_binder_alone_unchanged(
    out_dir: Path, existing: set[str], binder_alone: bool
) -> None:
    """Refuse to flip ``--binder-alone`` in a leaf that already holds output.

    Collect stamps the whole leaf from one sidecar, so mixed modes would mislabel every
    row of the other. A missing key means False; an empty leaf has nothing to mislabel.
    """
    if not existing:
        return
    previous = bool(_read_meta(out_dir).get(BINDER_ALONE_KEY, False))
    if previous == bool(binder_alone):
        return
    raise ValueError(
        f"{out_dir} already holds design_ph output from a run with "
        f"atomiumph_binder_alone={previous}, and this submit has "
        f"{bool(binder_alone)}. The two are different questions -- one has the target "
        f"in the graph, the other does not -- and collect stamps the whole leaf from "
        f"ONE sidecar, so the rows of one mode would be published with the other's "
        f"flag and its NA pattern. --binder-alone always gets its OWN --table-label "
        f"(and its own leaf); the two modes must never share one."
    )


def _check_no_stale_outputs(out_dir: Path, existing: set[str], wanted: set[str]) -> None:
    """Refuse to submit into a leaf holding results this submit will not regenerate.

    ``grp_<g>`` numbering restarts at 0 and nothing prunes old output, so collect would
    stamp older results with this run's parameters. Scanned leaf-wide, since a smaller
    rerun never revisits ``grp_1/out``. Fail-closed and non-destructive: it names the
    stems and deletes nothing.
    """
    stale = existing - wanted
    if stale:
        shown = ", ".join(sorted(stale)[:10]) + (" …" if len(stale) > 10 else "")
        raise ValueError(
            f"{out_dir} already holds design_ph output for {len(stale)} stem(s) this "
            f"submit will NOT regenerate: {shown}. grp_<g> numbering restarts at 0 on "
            f"every submit and the run sidecar is rewritten, so collect would stamp "
            f"those older results with THIS run's parameters (--model-noise, "
            f"--center-mode, the chains, --apo-weight, --protect-titratable). "
            f"RETRYING A FAILED SUBSET NEEDS A NEW -l/--dir-label: the two leaves "
            f"collect into the SAME table, so nothing is lost by splitting them. "
            f"Alternatively re-submit the full set, or remove that directory yourself "
            f"if the old results are finished with. Nothing is deleted automatically."
        )


def _stage_inputs(inputs_dir: Path, members: list[Member]) -> None:
    """Make ``inputs_dir`` hold exactly this subgroup's designs, as symlinks.

    ``grp_<gid>`` restarts at 0 each submit and ``--pdb_dir`` globs the whole directory,
    so leftover links would re-run stale designs under the new parameters.
    """
    wanted = {name for name, _ in members}
    for stale in inputs_dir.glob("*.pdb"):
        if stale.stem not in wanted:
            stale.unlink()
    for name, pdb_src in members:
        _symlink_member(pdb_src, inputs_dir, name)


def _symlink_member(pdb_src: Path, inputs_dir: Path, name: str) -> None:
    """Symlink a staged design as ``<name>.pdb`` (design_ph names outputs by stem)."""
    link = inputs_dir / f"{name}.pdb"
    if link.is_symlink():
        # A stale link from a rerun with a changed -i would re-design the old file.
        if link.readlink() == pdb_src:
            return
        link.unlink()
    elif link.exists():
        raise FileExistsError(
            f"{link} exists and is not a symlink; refusing to overwrite it (delete the "
            f"leaf dir, or rerun with a new -l/--dir-label)"
        )
    link.symlink_to(pdb_src)


def build_atomiumph_manifest(ctx: ManifestCtx[AtomiumphArgs]) -> list[tuple[str, ...]]:
    # Guardrails first: all cheap, all before a single container starts.
    binder_chains, target_chains = _check_chains(ctx.args)
    lam_values = _check_lam(
        list(ctx.args.lam), ctx.args.no_control, ctx.args.binder_alone
    )
    _check_choices(ctx.args)
    _check_set(list(ctx.args.set), ctx.args.binder_alone)

    extra = _design_ph_extra(ctx.args)
    lam_field = " ".join(lam_values)

    ready = ctx.ready
    groups: dict[Signature, list[Member]] = defaultdict(list)
    for design_name in sorted(cast(str, n) for n in ready.index):
        input_path = Path(str(ready.at[design_name, ctx.args.input_column]))
        # Local path first: volume_path gives the container's /runs form.
        local_pdb = ensure_pdb(input_path, ctx.args.run_dir)
        if ctx.args.binder_alone:
            _check_binder_alone_input(design_name, local_pdb, binder_chains)
        pdb_src = volume_path(local_pdb)
        groups[(binder_chains, target_chains)].append((design_name, pdb_src))

    x = ctx.args.designs_per_task
    chunks: list[Subgroup] = []
    for sig, members in groups.items():
        for i in range(0, len(members), x):
            chunks.append((sig, members[i : i + x]))
    tasks = _pack_subgroups(chunks, x)

    # Once, leaf-wide, before any grp_<g> dir is created or symlink moved.
    existing_stems = _output_stems(ctx.out_dir)
    _check_binder_alone_unchanged(ctx.out_dir, existing_stems, ctx.args.binder_alone)
    _check_no_stale_outputs(
        ctx.out_dir,
        existing_stems,
        {name for members in groups.values() for name, _ in members},
    )

    tasks_dir = ctx.out_dir / "atomiumph_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[tuple[str, ...]] = []
    gid = 0
    for t, subgroups in enumerate(tasks):
        sub_rows: list[tuple[str, ...]] = []
        for (binder_sig, target_sig), members in subgroups:
            grp_dir = ctx.out_dir / f"grp_{gid}"
            inputs_dir = grp_dir / "inputs"
            inputs_dir.mkdir(parents=True, exist_ok=True)
            _stage_inputs(inputs_dir, members)
            sub_rows.append((str(grp_dir), binder_sig, target_sig, lam_field, extra))
            gid += 1
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            f.writelines("\t".join(row) + "\n" for row in sub_rows)
        manifest_rows.append((str(task_file),))

    # LAST: a refused or aborted submit must not rewrite the recorded run parameters.
    ctx.write_meta(
        **{
            EXPECTED_KEY: len(lam_values) * ctx.args.num_seq,
            BINDER_CHAINS_KEY: binder_chains,
            TARGET_CHAINS_KEY: target_chains,
            CENTER_MODE_KEY: ctx.args.center_mode,
            OUTPUT_SCOPE_KEY: ctx.args.output_scope,
            MODEL_NOISE_KEY: ctx.args.model_noise,
            LAM_SWEEP_KEY: lam_field,
            PROTECT_TITRATABLE_KEY: str(ctx.args.protect_titratable),
            APO_WEIGHT_KEY: _apo_weight_value(ctx.args),
            CHARGE_PRIOR_KEY: _charge_prior_value(ctx.args),
            BINDER_ALONE_KEY: bool(ctx.args.binder_alone),
        }
    )

    return manifest_rows


def _fmt_lam(value: float) -> str:
    """Render a lambda for argv without exponent notation or a stray trailing dot."""
    text = f"{value:.6f}".rstrip("0")
    return text + "0" if text.endswith(".") else text


def add_run_atomiumph_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--binder-chains",
        type=str,
        default="",
        help="REQUIRED. Chain(s) of the binder to redesign, chain mini-language "
        "(':' is an inclusive letter range, ',' separates: 'A:C,E'). Passed to "
        "design_ph.py as a space-separated --binder_chains. No default on purpose: "
        "upstream's 'A' is wrong for the rfd3 binder complexes in this workspace, "
        "which usually carry the target on A and the binder on B.",
    )
    parser.add_argument(
        "--target-chains",
        type=str,
        default="",
        help="REQUIRED unless --binder-alone. Chain(s) of the target, same "
        "mini-language. These are the chains DELETED to build the apo reference, so "
        "ddE is defined only when they are given and disjoint from --binder-chains.",
    )
    parser.add_argument(
        "--binder-alone",
        action="store_true",
        help="Run design_ph.py on the BINDER ONLY: no target chain in the graph. This "
        "is a DIFFERENT QUESTION, not a cheaper version of the normal run. It exists "
        "for --placement-region core: design_ph defines core as a Cb-neighbour count "
        "within 10 A over EVERY chain present (region_masks), so with a target in the "
        "file the centre can be buried under the TARGET and be fully exposed in the "
        "free binder. With the target absent, core is the binder's own core. "
        "CONSEQUENCES, all structural: ddE is undefined (no apo contrast -- design_ph "
        "writes nan and the collector reports NA, never 0); every target-dependent "
        "column is NA and each row carries atomiumph_binder_alone True; "
        "--apo-weight and --charge-prior are FORCED TO 0.0 and recorded as 0.0. "
        "REFUSED in this mode: a non-empty --target-chains, any --lam but 0.0, "
        "--placement-region interface (the mask is all False without a target), "
        "--set '--placement_metric ddE' / '--coupling_metric ddE', a staged structure "
        "holding any chain other than --binder-chains (checked per file at submit -- "
        "cut the binder out with chainsel and point -i at that column), and flipping "
        "this flag in a leaf that already holds output. GIVE IT ITS OWN "
        "--table-label: the two modes must never share a leaf.",
    )
    parser.add_argument(
        "--center-mode",
        type=str,
        default="",
        choices=list(CENTER_MODES),
        help="REQUIRED, no default. Which protonation centres to pin, i.e. which pH "
        "behaviour you are designing. his_s: neutral His -- binds at NEUTRAL pH, "
        "releases in acid (the recycling / endosome-release direction; TNF-alpha). "
        "his_p: protonated His -- binds in ACID, releases at neutral (tumour or "
        "endosomal capture; EGFR). acid_p / all_p: Asp/Glu, which sit near pKa 4 and "
        "only switch below ~5.8 -- DIAGNOSTIC ONLY, not usable for either His "
        "challenge (README v2 l.30-32). Upstream defaults to acid_p; that default is "
        "deliberately NOT inherited, because it is wrong for both real use cases and "
        "produces a healthy-looking table of the wrong molecule. NOTE his_s reverses "
        "the charge prior -- it seeks a POSITIVE (Arg/Lys) target environment.",
    )
    parser.add_argument(
        "--centers-per-type",
        type=int,
        default=2,
        help="How many centres to pin per TOKEN TYPE (default 2, README v2's "
        "TNF-alpha recipe: one histidine gives a transition too shallow to separate "
        "6.5 from 7.4 when the pH window is <=1.5 units, l.201-203. Upstream's own "
        "default is 1). "
        "his_p/his_s have one token type, acid_p two (Asp, Glu) and all_p three, so 2 "
        "here is 2, 4 or 6 centres respectively. Every extra centre scales the "
        "--redesign-scope coupling scan linearly and adds a value to the ;-joined "
        "per-centre metrics cells. Only the first 12 centres get an "
        "atomiumph_centre<k>_idx1 position column; past that the count column "
        "atomiumph_n_centres_idx1 falls below atomiumph_n_centers and "
        "atomiumph_centres_idx1_basis reads 'truncated', so nothing is lost silently.",
    )
    parser.add_argument(
        "--placement-region",
        type=str,
        default="interface",
        choices=list(PLACEMENT_REGIONS),
        help="Where centres may be placed (default interface). 'interface' uses "
        "--interface-cutoff, which defaults to a strict 8.0 A.",
    )
    parser.add_argument(
        "--model-noise",
        "--noise",
        dest="model_noise",
        type=str,
        default="n03_v2",
        choices=list(CHECKPOINT_NOISES),
        help="Noise level of the fine-tuned 30-class (protonation-aware) checkpoint: "
        "n03, n05, n03_v2 or n05_v2 -- exactly what design_ph.py accepts, resolving "
        "to ft_weights_for_pH/ft_model_epoch10_<noise>.pt. tools/atomium's n00..n07 "
        "are the 21-class BASE weights and are refused here. The _v2 pair are "
        "DIFFERENT weights, not newer builds of the other two, so never rank rows "
        "from two levels together. Default n03_v2 (a user decision of 2026-10-09; "
        "upstream's client default is n05) -- the ONE deliberate exception to this "
        "tool's defaults-are-upstream's rule.",
    )
    parser.add_argument(
        "--lam",
        type=float,
        nargs="+",
        default=[0.0, 0.25, 0.5, 0.75],
        metavar="L",
        help="Lambda sweep: the weight of the selective term "
        "log p(target context) - lambda * log p(comparison context). 0.0 is the "
        "internal control (no selection pressure) and is REQUIRED unless "
        "--no-control is passed. --num-seq sequences are drawn at EACH value, so the "
        "row count per structure is len(lam) x num_seq.",
    )
    parser.add_argument(
        "--no-control",
        action="store_true",
        help="Deliberately allow a --lam sweep without 0.0. Only use this when the "
        "0.0 control for these same backbones already exists in another table.",
    )
    parser.add_argument(
        "--num-seq",
        type=int,
        default=8,
        help="Sequences PER LAMBDA VALUE (default 8). Total rows per structure is "
        "len(--lam) x this.",
    )
    parser.add_argument("--temperature", type=float, default=0.1,
                        help="Sampling temperature for the redesign (default 0.1).")
    parser.add_argument(
        "--seed-mode",
        type=str,
        default="native",
        choices=list(SEED_MODES),
        help="'native' (default) redesigns from the sequence in the input PDB -- "
        "correct when the binder is already designed. 'generate' first produces a "
        "full binder sequence (plain mask-predict, centres unpinned) and redesigns "
        "from that -- needed for a raw poly-Gly backbone.",
    )
    parser.add_argument("--seed-temperature", type=float, default=0.1,
                        help="Temperature of the --seed-mode generate pass (default 0.1).")
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default 42). Every structure is seeded identically, so a "
        "structure gives the same designs alone or inside a batch.",
    )
    parser.add_argument(
        "--redesign-scope",
        type=str,
        default="coupling",
        choices=list(REDESIGN_SCOPES),
        help="Which positions may change. 'coupling' (default: README v2 calls it "
        "the 'default choice for production', l.211; upstream's own default is "
        "neighbourhood). 'neighbourhood': "
        "graph neighbours of the centres, capped by --neighbour-k and --max-mutations. "
        "'binder': every designable binder position; for a raw backbone or an "
        "untrusted sequence. 'coupling': ranked by MEASURED influence on the centres' "
        "dE -- best quality, writes a <stem>.couplings.csv alanine scan, and costs a "
        "scan of n_candidates x n_centers forward passes per structure.",
    )
    parser.add_argument("--neighbour-k", type=int, default=12,
                        help="Neighbours considered per centre (default 12).")
    parser.add_argument("--max-mutations", type=int, default=20,
                        help="Cap on mutations per design (default 20).")
    parser.add_argument("--coupling-radius", type=float, default=12.0,
                        help="Radius of the coupling candidate set, A (default 12.0: "
                        "README v2 calls 12 'the recommended setting' -- 80%% of the "
                        "measured influence at 3.4x less compute than 20 A, l.213. "
                        "Upstream's own default is 15.0).")
    parser.add_argument(
        "--coupling-probe",
        type=str,
        default="ala",
        choices=list(COUPLING_PROBES),
        help="Probe substitution used to measure coupling: 'ala' (default, Gly where "
        "already Ala) or 'mask'.",
    )
    parser.add_argument(
        "--ban-cys",
        type=int,
        default=1,
        choices=[0, 1],
        help="1 (default) bans cysteine: a disulfide is redox- not pH-switched, an "
        "interchain one makes the binder covalent, and unpaired Cys break expression.",
    )
    parser.add_argument(
        "--apo-weight",
        type=float,
        default=0.5,
        help="Weight of the apo term in the sampling objective, forwarded as "
        "design_ph.py's --apo_weight. Default 0.5 = THE UPSTREAM CLIENT'S value; "
        "design_ph.py's own default is 0.0, so this flag is modelled precisely "
        "because omitting it would NOT reproduce upstream. Above 0 the apo graph "
        "enters SAMPLING and not only scoring, which is what pushes a design toward "
        "selectivity rather than toward being merely good in the complex; it also "
        "roughly doubles the infill forward passes. Recorded as atomiumph_apo_weight. UNDER --binder-alone THIS FLAG IS FORCED TO 0.0 -- there is no apo graph, and whether every consumer of apo_weight guards on it was not read, so the value is zeroed rather than trusted to be ignored; the recorded column says 0.0.",
    )
    parser.add_argument(
        "--protect-titratable",
        type=int,
        default=1,
        choices=[0, 1],
        help="Forwarded verbatim as design_ph.py's --protect_titratable (default 1, a "
        "user decision of 2026-10-09). SEMANTICS UNVERIFIED: the likely meaning is "
        "'keep titratable residues (His/Asp/Glu) that are NOT pinned centres from "
        "being redesigned away', but nothing read here states it. The build dumps "
        "upstream's --help to /opt/atomium_ph/DESIGN_PH_HELP.txt and the task script "
        "echoes it into the .out -- READ IT on the first run and correct the skill.",
    )
    parser.add_argument(
        "--charge-prior",
        type=float,
        default=0.5,
        help="Weight of the charge-complementarity prior on centre placement "
        "(default 0.5; 0 disables). Its SIGN follows --center-mode: his_p/acid_p seek "
        "a negative target environment, his_s a positive one. UNDER --binder-alone IT "
        "IS FORCED TO 0.0: the field sums TARGET charges only, so there is nothing for "
        "it to weigh. Recorded as atomiumph_charge_prior.",
    )
    parser.add_argument("--charge-cutoff", type=float, default=8.0,
                        help="Radius for the charge field, A (default 8.0).")
    parser.add_argument(
        "--interface-cutoff",
        type=float,
        default=8.0,
        help="Interface radius for --placement-region interface, A (default 8.0: "
        "README v2 gives 6-8 as the STRICT FIRST SHELL, l.235-236, and the TNF-alpha "
        "worked example uses 8.0. Upstream's own default is 12.5, which is permissive "
        "and admits a second shell, so an 'interface' centre need not touch the "
        "target at all).",
    )
    parser.add_argument("--burial-cutoff", type=int, default=16,
                        help="Neighbour count above which a position counts as buried "
                        "(default 16).")
    parser.add_argument(
        "--output-scope",
        type=str,
        default="binder",
        choices=list(OUTPUT_SCOPES),
        help="'binder' (default) writes the binder sequence only; 'all' writes every "
        "chain joined by '/'. Leave it at 'binder' unless you know the downstream "
        "consumer splits on '/'. 'all' ALSO SUPPRESSES THE NUMERIC CENTRE POSITIONS: "
        "atomiumph_centre<k>_idx1 is an index into the designed sequence, which under "
        "'all' is the whole complex and not the binder chain, so the columns collect "
        "as NA with atomiumph_centres_idx1_basis = 'output_scope_all' rather than as "
        "a position on the wrong molecule. The packed atomiumph_centers audit cell is "
        "written either way.",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=5,
        help="Max structures per task (default 5, lower than atomium's 10 because a "
        "design_ph run is much more expensive). A task that hits the timeout loses the "
        "WHOLE chunk, including structures already finished -- outputs only come back "
        "when the call completes.",
    )
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="FLAG",
        help="Raw flag forwarded verbatim to design_ph.py (repeatable), for flags this "
        "tool does not model -- e.g. --set '--contact_cutoff 6.0' (design_ph.py has "
        "NO --checkpoint; that flag is the upstream CLIENT's only). Tokens "
        "are split on whitespace, so a value containing a space cannot be passed this "
        "way, and a tab or newline is refused (the sub-manifest is TSV). A token "
        "repeating a flag the tool already sets is REFUSED at submit time: argparse "
        "takes the last occurrence, so it would override a validated value invisibly.",
    )
