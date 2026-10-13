"""Submit tasks running AtomiUM on designs, spawning a new child table.

AtomiUM is a ProteinMPNN-like sequence designer (a PyG re-implementation with
noise-conditioned weights), and it ships ProteinMPNN's own helper scripts under
``proteinMPNN_helper_scripts/`` (on Modal; off Modal the activation script supplies
the helper directory, see ``atomium.sh``). Its jsonl interface is the same one
(``--jsonl_path``, ``--chain_id_jsonl``, ``--fixed_positions_jsonl``,
``--tied_positions_jsonl``, ``--bias_AA_jsonl``), so this builder is deliberately
the proteinmpnn builder with the same two mini-languages and the same batching:
designs whose params are identical -- the signature ``(chains, fixed_positions,
tie_mode, tied_positions)`` -- are grouped into a single batched ``atomium.py``
call (model loaded once), and groups are bin-packed onto tasks up to
``--designs-per-task`` designs each.

The pipeline (``parse_multiple_chains`` -> optional helper steps -> ``atomium.py``)
runs inline in ``atomium.sh``, once per group.

Chain mini-language (``--chains-to-design``): the chains you design, in ProteinMPNN
order. ``:`` is an inclusive letter range, ``,`` separates:

    --chains-to-design A:C,E     # -> "A B C E"
    --chains-to-design A,C        # -> "A C" (non-contiguous)
    --chains-to-design ''         # design all chains (default)

Position mini-language (``--fixed-positions`` / ``--tied-positions``): ``/`` breaks
chains (its groups map one-to-one onto ``--chains-to-design``, in order), ``,``
separates fragments within a chain, ``start:end`` expands inclusively, a single
position passes through; table column expressions live in ``{...}`` islands (resolved
up the lineage with ``+ - * //`` arithmetic), everything else is a literal integer:

    --fixed-positions 9:23/10,11,18:20,22    # -> "9 10 ... 23, 10 11 18 19 20 22"
    --tied-positions  1:8/1:8                # -> "1 2 3 4 5 6 7 8, 1 2 3 4 5 6 7 8"

Positions are 1-indexed within each parsed chain (every chain is renumbered to
1..L), NOT original PDB numbering.

Symmetry (``--symmetry auto|N``): ties all chains via ``make_tied_positions_dict
--homooligomer 1``, and makes the position mini-language describe a single
asymmetric unit (a lone group is broadcast across the designed chains). ``auto``
takes the order from the input's polymer chain count. Mutually exclusive with
``--tied-positions``.

What differs from the proteinmpnn tool, all verified against the checkout:

  * ``--model-noise`` picks the weights (``model_noised_=_<noise>.pt``) and replaces
    ProteinMPNN's model-name/backbone-noise knobs. AtomiUM has no ``--ca_only``,
    ``--use_soluble_model``, ``--save_score`` or ``--save_probs``.
  * This tool exposes no ``--batch-size``. CORRECTED 2026-10-11: ``atomium.py``'s own
    ``--batch_size`` was DEAD at commit 292e29c, which is what this note used to say,
    but it is LIVE at the pinned commit dee5a4d67ae02293cee59f86c86085c8fa11879d --
    default 0, meaning chunks are created dynamically up to the maximum node budget.
    It affects CHUNKING ONLY. The sampling count is unchanged and is still
    ``num_seq_per_target * len(temperatures)``.
  * ``--sampling-temp`` may carry SEVERAL temperatures, and each one produces a full
    ``--num-seq-per-target`` set (so N temps => N x num_seq sequences per design).
    It therefore travels as its own manifest field and is quoted in the task script;
    folding it into the unquoted argv tail would split "0.1 0.2" into two arguments.
  * ``--seed 0`` means "pick a random seed" (AtomiUM treats 0 as falsy).

Usage:
    sapia run atomium outputs/<run> --table table0 --table-label atomium \\
        --model-noise n05 --num-seq-per-target 8

    # scan two noise levels as two runs into two child tables
    sapia run atomium outputs/<run> -t table0 --table-label n03 --model-noise n03
    sapia run atomium outputs/<run> -t table0 --table-label n05 --model-noise n05
"""

import json
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path
from typing import cast

from prosapia.core import (
    CommonArgs,
    LookupFn,
    ManifestCtx,
)
from prosapia.core.executors import volume_path
from prosapia.utils import (
    ensure_pdb,
    expand_chain_spec,
    parse_positions,
    polymer_chain_names,
)

# The weight files shipped in the checkout: model_weights/model_noised_=_<noise>.pt.
# Guarded up front because a typo would otherwise burn a GPU container per design
# before torch.load raises.
MODEL_NOISES = ("n00", "n01", "n02", "n03", "n04", "n05", "n05_v2", "n06", "n07")


class AtomiumArgs(CommonArgs):
    chains_to_design: str
    fixed_positions: str
    tied_positions: str
    symmetry: str | None
    bias_aa: str
    set: list[str]
    model_noise: str
    num_seq_per_target: int
    sampling_temp: str
    seed: int
    designs_per_task: int


def _parse_chains(spec: str) -> str:
    """Expand the chain mini-language into a space-separated ``--chain_list`` string."""
    return " ".join(expand_chain_spec(spec))


def _parse_positions(
    spec: str,
    lookup: LookupFn,
    name: str,
    *,
    broadcast_to: int | None = None,
) -> str:
    """Expand the position mini-language into a ``--position_list``.

    Per-chain groups joined by ``", "``, positions within a group by spaces. Order is
    preserved and NOT de-duplicated (tied positions are index-parallel). Empty spec
    -> "". ``broadcast_to`` (set from ``--symmetry``) replicates a spec describing a
    SINGLE asymmetric unit across that many chains; any other group count passes
    through untouched.
    """
    groups = parse_positions(spec, lookup, name)
    if broadcast_to and len(groups) == 1:
        groups = groups * broadcast_to
    return ", ".join(" ".join(str(p) for p in group) for group in groups)


def _check_symmetry_spec(spec: str | None) -> None:
    """Validate ``--symmetry`` up front, before anything is submitted."""
    if spec is None or spec == "auto":
        return
    try:
        order = int(spec)
    except ValueError:
        raise ValueError(
            f"--symmetry: expected 'auto' or an integer chain count, got {spec!r} "
            f"(write the number of chains in the oligomer, e.g. --symmetry 12)"
        )
    if order < 2:
        raise ValueError(f"--symmetry: order must be at least 2, got {order}")


def _check_model_noise(noise: str) -> None:
    """Reject a noise level with no weight file, before any container starts."""
    if noise not in MODEL_NOISES:
        raise ValueError(
            f"--model-noise: {noise!r} has no weights in the AtomiUM checkout; "
            f"choose one of {', '.join(MODEL_NOISES)}"
        )


def _check_sampling_temp(spec: str) -> None:
    """Every ``--sampling-temp`` token must be a float; AtomiUM splits on whitespace."""
    tokens = spec.split()
    if not tokens:
        raise ValueError("--sampling-temp: at least one temperature is required")
    for token in tokens:
        try:
            float(token)
        except ValueError:
            raise ValueError(
                f"--sampling-temp: {token!r} is not a number (expected one or more "
                f"space-separated temperatures, e.g. '0.1' or '0.1 0.2')"
            )


def _resolve_chains(args: AtomiumArgs, pdb_src: Path) -> str:
    """One design's designed-chain list, as a ``--chain_list`` string.

    ``--chains-to-design`` wins when given. Otherwise ``--symmetry`` derives the list
    from the structure itself -- the first ``order`` polymer chains, where ``auto``
    takes every one. With neither flag, "" leaves every chain designed.
    """
    if args.chains_to_design:
        return _parse_chains(args.chains_to_design)
    if args.symmetry is None:
        return ""
    names = polymer_chain_names(pdb_src)
    order = len(names) if args.symmetry == "auto" else int(args.symmetry)
    return " ".join(names[:order])


def _write_bias_jsonl(spec: str, out_dir: Path) -> str:
    """Write a global AA-bias jsonl from an ``AA:bias`` spec; return its path.

    ``make_bias_AA.py`` is purely generative (no structure input), so the dict is
    written directly and the AtomiUM env is not needed on the submit node.
    """
    spec = spec.strip()
    if not spec:
        return ""
    bias: dict[str, float] = {}
    for token in spec.split():
        if ":" not in token:
            raise ValueError(
                f"--bias-aa: malformed pair {token!r} (expected 'AA:bias', "
                f"e.g. 'D:1.39')"
            )
        aa, value = token.split(":", 1)
        try:
            bias[aa] = float(value)
        except ValueError:
            raise ValueError(f"--bias-aa: bias for {aa!r} is not a number: {value!r}")
    bias_path = out_dir / "bias_AA.jsonl"
    with open(bias_path, "w") as f:
        f.write(json.dumps(bias) + "\n")
    return str(bias_path)


def _atomium_extra(args: AtomiumArgs, bias_path: str) -> str:
    """The run-wide atomium.py argv tail (same for every design).

    Whitespace-free by construction: the task script expands this unquoted so each
    token becomes its own argv entry. ``--sampling_temp`` is NOT here -- it may hold
    several space-separated temperatures and travels as its own quoted field.
    """
    parts = [
        "--model_noise",
        args.model_noise,
        "--num_seq_per_target",
        str(args.num_seq_per_target),
        "--seed",
        str(args.seed),
    ]
    if bias_path:
        parts += ["--bias_AA_jsonl", bias_path]
    parts.extend(args.set)
    return " ".join(parts)


# A design ready to run, and the (chains, fixed_positions, tie_mode, tied_positions)
# signature that decides which designs can share one batched atomium.py call.
Member = tuple[str, Path]  # (design_name, staged_pdb_source)
Signature = tuple[str, str, str, str]
Subgroup = tuple[Signature, list[Member]]


def _pack_subgroups(chunks: list[Subgroup], max_per_task: int) -> list[list[Subgroup]]:
    """First-fit-pack subgroup chunks (each <= max_per_task) into tasks.

    Each task holds at most ``max_per_task`` designs total, possibly spanning several
    signatures; every subgroup stays intact so batching is preserved.
    """
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


def _symlink_member(pdb_src: Path, inputs_dir: Path, name: str) -> None:
    """Symlink a staged design into ``inputs_dir`` as ``<name>.pdb`` (no copy).

    The stem is what collect keys lineage off: AtomiUM names each output FASTA after
    the parsed input, so ``<name>.pdb`` -> ``seqs/<name>.fa`` -> parent row ``name``.
    """
    link = inputs_dir / f"{name}.pdb"
    if not link.exists() and not link.is_symlink():
        link.symlink_to(pdb_src)


def build_atomium_manifest(ctx: ManifestCtx[AtomiumArgs]) -> list[tuple[str, ...]]:
    # A present input --input-column (ctx.ready) already implies the upstream step
    # succeeded (collect writes a _path only on status OK), so no extra prefilter.
    ready = ctx.ready

    # Run-wide pieces (all structure-independent): argv tail and the bias jsonl.
    bias_path = _write_bias_jsonl(ctx.args.bias_aa, ctx.out_dir)
    atomium_extra = _atomium_extra(ctx.args, bias_path)

    # Guardrails, all cheap and all before submission.
    _check_model_noise(ctx.args.model_noise)
    _check_sampling_temp(ctx.args.sampling_temp)
    _check_symmetry_spec(ctx.args.symmetry)
    if ctx.args.symmetry is not None and ctx.args.tied_positions:
        raise ValueError(
            "--symmetry (homo-oligomer tie) and --tied-positions (explicit tie) are "
            "mutually exclusive"
        )
    if (
        (ctx.args.fixed_positions or ctx.args.tied_positions)
        and not ctx.args.chains_to_design
        and ctx.args.symmetry is None
    ):
        raise ValueError(
            "--fixed-positions/--tied-positions require --chains-to-design (their "
            "groups map one-to-one onto those chains), or --symmetry to derive the "
            "chains and broadcast one group across them"
        )

    # Stage each design (CIF->PDB via the shared cache; PDBs returned as-is) and
    # compute its signature. Positions may resolve per-design via lineage, so grouping
    # by signature lets same-param designs share one batched atomium.py call.
    groups: dict[Signature, list[Member]] = defaultdict(list)
    for design_name in sorted(cast(str, n) for n in ready.index):
        input_path = Path(str(ready.at[design_name, ctx.args.input_column]))
        pdb_src = volume_path(ensure_pdb(input_path, ctx.args.run_dir))
        chains = _resolve_chains(ctx.args, pdb_src)
        # Under --symmetry a lone position group describes one asymmetric unit and is
        # replicated across the designed chains (the helper wants one group per chain).
        broadcast = len(chains.split()) if ctx.args.symmetry is not None else None
        fixed_pl = _parse_positions(
            ctx.args.fixed_positions,
            ctx.lookup,
            design_name,
            broadcast_to=broadcast,
        )
        if ctx.args.symmetry is not None:
            tie_mode, tied_pl = "homo", ""
        elif ctx.args.tied_positions:
            tie_mode, tied_pl = (
                "explicit",
                _parse_positions(ctx.args.tied_positions, ctx.lookup, design_name),
            )
        else:
            tie_mode, tied_pl = "", ""
        sig: Signature = (chains, fixed_pl, tie_mode, tied_pl)
        groups[sig].append((design_name, pdb_src))

    # Split each signature group into subgroups of <= X, then bin-pack subgroups into
    # tasks with a total budget of X designs each.
    x = ctx.args.designs_per_task
    chunks: list[Subgroup] = []
    for sig, members in groups.items():
        for i in range(0, len(members), x):
            chunks.append((sig, members[i : i + x]))
    tasks = _pack_subgroups(chunks, x)

    # Materialize: one grp_<g>/ output dir per subgroup (inputs/ holds symlinks), one
    # sub-manifest per task (a row per subgroup), and a 1-field top-level row per task.
    tasks_dir = ctx.out_dir / "atomium_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[tuple[str, ...]] = []
    gid = 0
    for t, subgroups in enumerate(tasks):
        sub_rows: list[tuple[str, ...]] = []
        for (chains_sig, fixed_pl, tie_mode, tied_pl), members in subgroups:
            grp_dir = ctx.out_dir / f"grp_{gid}"
            inputs_dir = grp_dir / "inputs"
            inputs_dir.mkdir(parents=True, exist_ok=True)
            for name, pdb_src in members:
                _symlink_member(pdb_src, inputs_dir, name)
            sub_rows.append(
                (
                    str(grp_dir),
                    chains_sig,
                    fixed_pl,
                    tie_mode,
                    tied_pl,
                    ctx.args.sampling_temp,
                    atomium_extra,
                )
            )
            gid += 1
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            f.writelines("\t".join(row) + "\n" for row in sub_rows)
        manifest_rows.append((str(task_file),))

    return manifest_rows


def add_run_atomium_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--model-noise",
        type=str,
        default="n05",
        metavar="N00..N07",
        help="Noise level of the AtomiUM weights to design with, i.e. the sigma the "
        "model was trained to denoise: n00=0.0A, n01=0.1A ... n07=0.7A (plus the "
        "n05_v2 retrain). Higher tolerates rougher backbones; lower trusts the input "
        "geometry more. Default n05. Selects model_weights/model_noised_=<noise>.pt, "
        "so an unknown value is rejected at submit time.",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=10,
        help="Max designs per task (default 10). Designs are auto-grouped by identical "
        "params (chains, fixed/tied positions, symmetry) so a group shares one batched "
        "atomium.py call, and groups are bin-packed up to this budget. Lower it for "
        "more parallelism, raise it for fewer tasks.",
    )
    parser.add_argument(
        "--chains-to-design",
        type=str,
        default="",
        help="Chains to design, in ProteinMPNN order, as a chain mini-language: ':' is "
        "an inclusive letter range and ',' separates (e.g. 'A:C,E' -> 'A B C E'). Passed "
        "as --chain_list to assign_fixed_chains and, for the position flags, to "
        "make_{fixed,tied}_positions_dict. Empty (default): design all chains, or the "
        "structure's first N when --symmetry says so.",
    )
    parser.add_argument(
        "--fixed-positions",
        type=str,
        default="",
        help="Positions to keep FIXED (not redesigned), as a position mini-language: "
        "'/' breaks chains (groups map one-to-one onto --chains-to-design, in order), "
        "',' separates fragments within a chain, 'start:end' is an inclusive range; table "
        "column expressions go in {...} islands (resolved up the lineage with + - * // "
        "arithmetic), everything else is a literal integer. 1-indexed within each parsed "
        "chain. E.g. '9:23/10,11,18:20,22'. One group per designed chain is expected, "
        "unless --symmetry broadcasts a single group across them. Empty (default): "
        "redesign everything.",
    )
    parser.add_argument(
        "--tied-positions",
        type=str,
        default="",
        help="Positions to TIE across chains for symmetric design, same position "
        "mini-language as --fixed-positions (groups map to --chains-to-design; tied "
        "groups must be equal length and are tied index-parallel). E.g. '1:8/1:8'. "
        "Mutually exclusive with --symmetry. Empty (default): no explicit tie.",
    )
    parser.add_argument(
        "--symmetry",
        type=str,
        default=None,
        metavar="AUTO|N",
        help="Homo-oligomer convenience. Ties all chains via make_tied_positions_dict "
        "--homooligomer 1, AND lets --fixed-positions describe a single asymmetric unit: "
        "a lone position group is broadcast across the designed chains. 'auto' takes the "
        "order from the input structure's polymer chain count; otherwise pass that count "
        "as a plain integer. With no --chains-to-design, the designed chains are the "
        "structure's first N. Omitted by default. Mutually exclusive with "
        "--tied-positions.",
    )
    parser.add_argument(
        "--bias-aa",
        type=str,
        default="",
        help="Global AA composition bias as space-separated AA:bias pairs "
        "(e.g. 'D:1.39 E:1.39 K:1.39'), applied via --bias_AA_jsonl. "
        "Empty (default): no bias.",
    )
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="FLAG",
        help="Raw flag forwarded verbatim to atomium.py (repeatable); also accepts "
        "pre-made jsonl paths. E.g. --set '--omit_AAs AC', --set '--max_length 5000', "
        "--set '--pssm_jsonl /path/pssm.jsonl'. Tokens are split on whitespace, so a "
        "value containing a space cannot be passed this way.",
    )
    parser.add_argument("--num-seq-per-target", type=int, default=2)
    parser.add_argument(
        "--sampling-temp",
        type=str,
        default="0.1",
        help="One or more space-separated sampling temperatures (e.g. '0.1' or "
        "'0.1 0.2'). AtomiUM generates a full --num-seq-per-target set at EACH "
        "temperature, so N temperatures give N x num_seq sequences per design. "
        "Default '0.1'.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=37,
        help="Random seed (default 37). 0 means pick a random seed per task.",
    )
