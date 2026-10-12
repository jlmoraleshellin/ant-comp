"""Submit tasks running AtomiUM's ``design_ph.py`` on binder COMPLEXES, spawning a
child table with one row per designed sequence.

PREMISE -- read this before using any column it produces
--------------------------------------------------------
``design_ph.py`` designs **pH-selective** binder sequences: it pins titratable binder
residues (His/Asp/Glu) into a chosen protonation state and redesigns their environment
so that state is favoured *in the complex more than in the free binder*. Its selection
metric is therefore a DIFFERENCE OF DIFFERENCES::

    dE  = log p(target_state) - log p(comparison_state)   (one center, one context)
    ddE = dE(complex) - dE(apo binder)

The input must be a **complex**: a binder chain set AND a target chain set, in the same
file. Run it on a lone binder and there is no apo/complex contrast at all -- every ddE
is structurally zero and the numbers mean nothing. There is no way for this tool to
detect that for you beyond the chain checks below, so ``--binder-chains`` and
``--target-chains`` are both REQUIRED and must be disjoint.

``action="create"``: one child row per **(lambda, sample) draw**, named
``<parent>_p<i>``. ``--num-seq`` is per lambda value, so the row count per input
structure is ``len(--lam) x --num-seq``.

REFUSED AT SUBMIT
-----------------
* ``--binder-chains`` / ``--target-chains`` stay REQUIRED. The upstream client defaults
  the binder to ``"A"`` and the target to ``""``; here the binder is the LAST chain, a
  wrong-but-present letter redesigns the target with ``status == "OK"``, and an empty
  target makes every ddE structurally zero.
* ``--lam`` must contain ``0.0``. The upstream client only WARNS; this refuses at
  submit, because the control cannot be recovered without re-running the sweep on GPU.
  The duplicate / rendered-lambda checks are this tool's own, as are
  ``--designs-per-task`` and the other run-shape flags.
* ``--lam >= 1.0`` is NOT refused here, though the upstream client refuses it (sampling
  becomes uniform); ``design_ph.py`` does not.

BATCHING, AND THE LINEAGE CONTRACT IT RESTS ON
----------------------------------------------
Mirrors ``tools/atomium``: designs sharing a parameter signature are staged as symlinks
into ``grp_<g>/inputs/<design_name>.pdb`` and run as ONE batched
``design_ph.py --pdb_dir ... --pattern '*.pdb'`` call (checkpoint loaded once), and
subgroups are bin-packed onto tasks up to ``--designs-per-task`` designs each. The
symlink STEM is the parent-table row name, and ``design_ph.py`` names every output file
after the input stem -- which is the whole of the lineage contract (see
``collect_atomiumph.py``).
"""

from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path
from typing import cast

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.utils import ensure_pdb, expand_chain_spec

# The fine-tuned 30-class (protonation-aware) checkpoints. DELIBERATELY NOT shared with
# tools/atomium's MODEL_NOISES, which lists nine levels: the pH head exists only at
# these two, and importing the other tuple would let n07 through to a GPU container that
# then dies on a missing file.
# Exactly what design_ph.py accepts at the pin (design_ph.py:2132-2134). All four are
# fine-tuned 30-class pH checkpoints under ft_weights_for_pH/; the _v2 pair are
# DIFFERENT weights, not a newer build of the same ones, so a table must never mix them
# (atomiumph_model_noise is what says which). tools/atomium's n00..n07 are the 21-class
# base weights and are still refused here.
CHECKPOINT_NOISES = ("n03", "n05", "n03_v2", "n05_v2")

CENTER_MODES = ("his_p", "his_s", "acid_p", "all_p")
PLACEMENT_REGIONS = ("interface", "core", "surface", "any")
REDESIGN_SCOPES = ("neighbourhood", "binder", "coupling")
COUPLING_PROBES = ("ala", "mask")
SEED_MODES = ("native", "generate")
OUTPUT_SCOPES = ("binder", "all")

# Sidecar keys the collector reads back as trust columns. The chains are run-wide and
# would otherwise die in the sub-manifest: without them the table cannot say WHICH
# molecule these energies describe.
EXPECTED_KEY = "atomiumph_expected_per_structure"
BINDER_CHAINS_KEY = "atomiumph_binder_chains"
TARGET_CHAINS_KEY = "atomiumph_target_chains"
# Run-wide and decisive for how every row must be READ: --center-mode decides what a
# high ddE_mean even means (which pH binds), and --output-scope all makes `sequence` a
# slash-joined multi-chain string, at which point seq_len is no longer a binder length.
CENTER_MODE_KEY = "atomiumph_center_mode"
OUTPUT_SCOPE_KEY = "atomiumph_output_scope"
# MEASURED 2026-10-09: the model-noise level decides WHICH CENTRES design_ph picks --
# two runs over the same structures chose different centres at n03 and n05 -- so two
# noise levels are not comparable and the table must say which one a row came from.
# Until now the only record was the dir label. The lambda sweep is recorded with it
# because it is what makes `n_expected` a count rather than a number, and it says
# whether the lambda=0.0 control arm was in this run at all.
MODEL_NOISE_KEY = "atomiumph_model_noise"
LAM_SWEEP_KEY = "atomiumph_lam_sweep"
# Recorded for the same reason as the noise level: it changes WHICH positions the model
# was allowed to touch, so two runs that differ only in it are not comparable, and the
# table is otherwise silent about it. SEMANTICS UNVERIFIED -- see add_run_atomiumph_args.
PROTECT_TITRATABLE_KEY = "atomiumph_protect_titratable"
# Modelled ONLY because the upstream client's default (0.5) differs from design_ph.py's
# own (0.0): omitting the flag would silently run at 0.0 and not match upstream. It
# changes what is sampled, so it is recorded like the noise level.
APO_WEIGHT_KEY = "atomiumph_apo_weight"

# Flags this builder emits itself, mapped to the typed flag that owns them. design_ph.py
# is argparse-based, so a later occurrence WINS: a --set token repeating one of these
# would silently override a value the guardrails just validated, appear in no column,
# and make n_expected a lie. Refused at submit time instead.
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

    Both are REQUIRED (upstream defaults ``--binder_chains A --target_chains ''``; we
    refuse those because the premise is a complex, and because the binder is NOT
    reliably chain A here -- rfd3 binder output in this workspace commonly carries the
    target on A and the binder on B, so inheriting upstream's default would quietly
    design the target and score the binder as context).
    """
    binder = _chain_list(args.binder_chains, "--binder-chains")
    target = _chain_list(args.target_chains, "--target-chains")
    if not binder:
        raise ValueError(
            "--binder-chains is required: name the chain(s) to redesign, e.g. "
            "'B' (rfd3 binder complexes in this workspace are usually target=A, "
            "binder=B -- check the structure, do not assume)"
        )
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


def _check_lam(lams: list[float], no_control: bool) -> list[str]:
    """``--lam`` must be a non-empty sweep containing the 0.0 control; returns argv text.

    Both checks run on the RENDERED strings, which is what design_ph.py actually
    receives: ``_fmt_lam`` rounds to 6 decimals, so ``--lam 0.0 1e-7`` is two identical
    sweep arms however different the floats look here.
    """
    if not lams:
        raise ValueError("--lam: at least one lambda value is required")
    rendered = [_fmt_lam(x) for x in lams]
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
    """Re-validate every enumerated flag in the builder, not only in argparse.

    argparse ``choices`` catches a CLI typo; this catches a programmatic caller and
    makes the failure a submit-time ValueError naming the flag.
    """
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
    # Same treatment as --binder-chains/--target-chains: a REQUIRED choice, refused
    # at submit rather than defaulted, because every value of it produces a complete
    # and healthy-looking table -- of a different molecule's pH behaviour.
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


def _check_set(tokens: list[str]) -> None:
    """``--set`` is for UNMODELLED flags only, and must survive a TSV round trip.

    Two refusals, both at submit time:

    * a token repeating a flag this builder emits. design_ph.py's argparse takes the
      LAST occurrence, so ``--set '--checkpoint_model_noise n07'`` would override a
      value the guardrails just validated, show up in no column, and make the
      ``n_expected`` trust column a lie.
    * a token carrying a tab or a newline. The sub-manifest is TSV and ``atomiumph.sh``
      pulls the argv tail with ``cut -f5``: a tab truncates the tail at that point and a
      newline forges an extra subgroup row.
    """
    for token in tokens:
        for bad, what in (("\t", "a tab"), ("\n", "a newline"), ("\r", "a carriage return")):
            if bad in token:
                raise ValueError(
                    f"--set: {token!r} contains {what}. The sub-manifest is TSV and the "
                    f"task script cuts the argv tail out of it, so that character would "
                    f"silently truncate or forge a row."
                )
        # EVERY word, not just the first: a --set string is explicitly allowed to carry
        # a flag and its value, so '--contact_cutoff 6.0 --binder_chains A' would
        # otherwise sail past a first-word-only check and reassign the binder.
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


def _design_ph_extra(args: AtomiumphArgs) -> str:
    """The run-wide ``design_ph.py`` argv tail, shared by every subgroup.

    WHITESPACE-FREE BY CONSTRUCTION: the task script expands this unquoted so each
    token becomes its own argv entry. ``--binder_chains``/``--target_chains`` (single
    space-separated strings) and ``--lam`` (several values) are NOT here -- each
    travels as its own manifest field, see ``atomiumph.sh``.
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
        "--apo_weight", str(args.apo_weight),
        "--charge_prior", str(args.charge_prior),
        "--charge_cutoff", str(args.charge_cutoff),
        "--interface_cutoff", str(args.interface_cutoff),
        "--burial_cutoff", str(args.burial_cutoff),
        "--output_scope", args.output_scope,
    ]
    parts.extend(args.set)
    return " ".join(parts)


# A design ready to run, and the signature deciding which designs share one batched
# design_ph.py call. Every design_ph parameter is currently run-wide, so the signature
# is constant across a run; it is kept so a future per-design flag (a residue mask, a
# per-row centre list) cannot silently batch incompatible designs together.
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


# Suffixes design_ph.py writes per stem. Mirrors the collector's SUFFIXES; kept here so
# the submit-time stale check can name exactly the stems that would be collected.
_OUTPUT_SUFFIXES = (".metrics.csv", ".couplings.csv", ".centers.json", ".ext.fa", ".fa")
# Dirs the RUN owns, which never hold design_ph output. Mirrors the collector's
# SKIP_DIRS so the two agree on what counts as a result.
_SKIP_DIRS = ("atomiumph_logs", "atomiumph_tasks", "inputs")


def _check_no_stale_outputs(out_dir: Path, wanted: set[str]) -> None:
    """Refuse to submit into a leaf holding results this submit will not regenerate.

    THE PROBLEM: ``grp_<gid>`` numbering restarts at 0 on every submit, nothing prunes
    ``grp_<g>/out``, and ``.meta.json`` is rewritten each submit. So re-running the same
    ``-l/--dir-label`` with a different ``--model-noise``, ``--center-mode``, chain
    assignment, ``--apo-weight`` or ``--protect-titratable`` leaves the PREVIOUS run's
    outputs in place, and collect would stamp them with the NEW run's parameters -- a
    whole parent's numbers attributed to settings that never produced them.

    SCANNED LEAF-WIDE AND ONCE, before any subgroup dir is touched. A per-subgroup check
    is not enough and was the original bug: a rerun with a SMALLER ready set -- the
    normal case after a partial failure or a ``-f`` filter -- recreates only ``grp_0``
    and never looks at ``grp_1/out``, whose stems the collector would then mint as OK
    rows under the new run's parameters. ``wanted`` is therefore every member of the
    WHOLE submit, and the scan is a recursive glob over the same suffixes the collector
    indexes.

    FAIL CLOSED AND NON-DESTRUCTIVE: refuses, names the stems, deletes nothing. A stem
    this submit WILL regenerate is fine (its files are overwritten); only a stem that
    would survive untouched is refused.
    """
    if not out_dir.is_dir():
        return
    stale: set[str] = set()
    for path in out_dir.rglob("*"):
        if not path.is_file():
            continue
        # Never look inside the dirs the RUN owns: the staged symlinks are named
        # <stem>.pdb (not an output suffix) but the sub-manifests and logs live here
        # too, and a future log name must not be mistaken for a design.
        if any(part in _SKIP_DIRS for part in path.relative_to(out_dir).parts[:-1]):
            continue
        for suffix in _OUTPUT_SUFFIXES:
            if path.name.endswith(suffix) and len(path.name) > len(suffix):
                stem = path.name[: -len(suffix)]
                if stem not in wanted:
                    stale.add(stem)
                break
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
    """Make ``inputs_dir`` hold EXACTLY this subgroup's designs, as symlinks.

    The pruning matters because ``grp_<gid>`` numbering restarts at 0 on every submit:
    re-running the same label with a different ready set would otherwise leave the
    previous occupants of ``grp_0/inputs`` in place, and ``--pdb_dir`` globs the whole
    directory -- so stale designs would be re-run under the NEW parameters and their
    outputs would overwrite the earlier ones under the same stem.
    """
    wanted = {name for name, _ in members}
    for stale in inputs_dir.glob("*.pdb"):
        if stale.stem not in wanted:
            stale.unlink()
    for name, pdb_src in members:
        _symlink_member(pdb_src, inputs_dir, name)


def _symlink_member(pdb_src: Path, inputs_dir: Path, name: str) -> None:
    """Symlink a staged design into ``inputs_dir`` as ``<name>.pdb`` (no copy).

    The stem is what collect keys lineage off: design_ph.py names every output after
    the input file stem, so ``<name>.pdb`` -> ``<name>.metrics.csv`` -> parent row.
    """
    link = inputs_dir / f"{name}.pdb"
    if link.is_symlink():
        # A --force rerun after changing -i leaves a link pointing at the PREVIOUS
        # structure. Reusing it would re-design the old file under the new column's
        # name, with nothing in the table to show it.
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
    lam_values = _check_lam(list(ctx.args.lam), ctx.args.no_control)
    _check_choices(ctx.args)
    _check_set(list(ctx.args.set))

    extra = _design_ph_extra(ctx.args)
    lam_field = " ".join(lam_values)

    # Recorded for the collector's trust columns. `n_expected` because --num-seq is PER
    # LAMBDA, so the row count per structure is a product; the chains because they are
    # run-wide and would otherwise exist nowhere a reader of the table can see.
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
            APO_WEIGHT_KEY: str(ctx.args.apo_weight),
        }
    )

    ready = ctx.ready
    groups: dict[Signature, list[Member]] = defaultdict(list)
    for design_name in sorted(cast(str, n) for n in ready.index):
        input_path = Path(str(ready.at[design_name, ctx.args.input_column]))
        pdb_src = volume_path(ensure_pdb(input_path, ctx.args.run_dir))
        groups[(binder_chains, target_chains)].append((design_name, pdb_src))

    x = ctx.args.designs_per_task
    chunks: list[Subgroup] = []
    for sig, members in groups.items():
        for i in range(0, len(members), x):
            chunks.append((sig, members[i : i + x]))
    tasks = _pack_subgroups(chunks, x)

    # ONCE, leaf-wide, before any grp_<g> dir is created or any symlink moved: the
    # stems of the WHOLE submit against every output already in the leaf.
    _check_no_stale_outputs(
        ctx.out_dir, {name for members in groups.values() for name, _ in members}
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
        help="REQUIRED. Chain(s) of the target, same mini-language. These are the "
        "chains DELETED to build the apo reference, so ddE is defined only when they "
        "are given and disjoint from --binder-chains.",
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
        "per-centre metrics cells.",
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
        "roughly doubles the infill forward passes. Recorded as atomiumph_apo_weight.",
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
        "a negative target environment, his_s a positive one.",
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
        "consumer splits on '/'.",
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
