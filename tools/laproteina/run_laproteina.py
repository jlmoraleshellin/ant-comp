#!/usr/bin/env python3
"""
Submit La-Proteina Complexa (LPC) binder-generation jobs, one array task per job.

LPC is a fully atomistic binder generator: a partially latent flow-matching model
that samples backbone geometry, side chains and sequence jointly, conditioned on a
target structure and a hotspot list. Its unit of work is a **generation job**
(``complexa generate <config> --job-id N``), which draws ``nsamples`` binders
against one registered target. So the unit of work here is **one job**, not one
design -- the designs only exist once the job has run, which is why this is a
``create`` tool and why the row count is known only at collect time.

The target is **not** a file path. LPC resolves targets by name out of its own
registry (``configs/targets/targets_dict.yaml``), which carries the structure, the
chain/residue ranges the model sees, the hotspots and the binder-length range. So
this tool takes ``--task-name``, and the target must already be registered on the
LPC side. ``sapia run laproteina --help`` cannot list them; ``complexa target list``
can.

Only the **generate** stage is run. LPC's own ``design`` pipeline would continue
into filtering, its own SolubleMPNN sequence design and AF2/RF3 refolding; this tool
deliberately stops at backbones so the sequences can be designed by this workspace's
own tools (``atomium``, ``proteinmpnn``) and validated independently with ``boltz``.
That is the same division of labour as ``bindcraft2 --trajectory-only``.

Fan-out is by job, and it is a **split, not a multiplier**: ``--nsamples`` is the
total for the run and ``--num-jobs N`` divides it across N array tasks. LPC does the
dividing itself in ``split_by_job``, which is also why ``gen_njobs`` must be set to
match the task count -- see ``build_laproteina_manifest``.

Usage:
    # root run: 20 backbones against a registered LPC target, over 4 tasks
    sapia run laproteina outputs/RUN \\
        --task-name 90_TNFa_HUMAN_xreact --num-jobs 4 --nsamples 20

    # smaller, for a first smoke test: 2 backbones, one task
    sapia run laproteina outputs/RUN \\
        --task-name 90_TNFa_HUMAN_xreact --num-jobs 1 --nsamples 2

    # pass anything else straight through to Hydra
    sapia run laproteina outputs/RUN --task-name 90_TNFa_HUMAN_xreact \\
        --set generation.dataloader.batch_size=8
"""

from argparse import ArgumentParser
from pathlib import Path

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path

TOOL_NAME = "laproteina"

# Where each task parks the directory LPC wrote, under the tool's out_dir. The
# collector scans this and nothing else.
SAMPLES_DIRNAME = "jobs"

DEFAULT_CONFIG = "configs/search_binder_local_pipeline.yaml"


class LaProteinaArgs(CommonArgs):
    task_name: str | None
    num_jobs: int
    nsamples: int | None
    config: str
    seed: int | None
    design_prefix: str | None
    search: str
    rewards: bool
    set: list[str]


def add_run_laproteina_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--task-name",
        type=str,
        default=None,
        help="The LPC target to design against, by its key in LPC's "
        "configs/targets/targets_dict.yaml (e.g. 90_TNFa_HUMAN_xreact). The entry "
        "carries the structure, the chain/residue ranges, the hotspots and the "
        "binder-length range -- none of those are passed from here. Run "
        "`complexa target list` to see what is registered. Required.",
    )
    parser.add_argument(
        "--num-jobs",
        type=int,
        default=1,
        help="How many array tasks to SPLIT the run across, one generation job "
        "each. This is parallelism, not extra work: --nsamples is the total and "
        "LPC divides it between the jobs. Sets LPC's own `gen_njobs` to match, "
        "which is mandatory -- see the note on split_by_job below. Default 1.",
    )
    parser.add_argument(
        "--nsamples",
        type=int,
        default=None,
        help="TOTAL binders drawn by the whole run, divided across --num-jobs. "
        "Overrides generation.dataloader.dataset.nres.nsamples (shipped default: "
        "4). Leave unset to take the config's value. NOTE this is a total, not a "
        "per-job count: --num-jobs 4 --nsamples 20 is 20 designs, 5 per job.",
    )
    parser.add_argument(
        "--config",
        type=str,
        default=DEFAULT_CONFIG,
        help=f"LPC pipeline config, relative to the LPC checkout. Default "
        f"{DEFAULT_CONFIG}. Only the generation.* section is used -- the evaluate "
        "and analyze sections are never reached, because this tool stops at "
        "backbones.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Base seed. Job j runs with seed+j, so jobs don't redraw the same "
        "binders. Leave unset to take the config's seed, in which case the jobs are "
        "still separated by --job-id but the seed itself is shared.",
    )
    parser.add_argument(
        "--design-prefix",
        type=str,
        default=None,
        help="Prefix for the design-group names, and so for every row this run "
        "mints. Defaults to --task-name, which is descriptive but long; pass "
        "something shorter (e.g. 'tnfh') to keep row names readable.",
    )
    parser.add_argument(
        "--search",
        type=str,
        default="single-pass",
        choices=["single-pass", "best-of-n", "beam-search", "fk-steering", "mcts"],
        help="Inference-time search. Defaults to single-pass (no search), because "
        "every other algorithm scores candidates with a reward model and the "
        "--minimal image has none installed. The shipped config's default is "
        "best-of-n; turning it back on needs --rewards AND an image built with "
        "JAX/ColabDesign.",
    )
    parser.add_argument(
        "--rewards",
        action="store_true",
        help="Keep the config's reward model (AF2) instead of setting it to null. "
        "Only meaningful on a full (non---minimal) image with AF2_DIR set; "
        "otherwise Hydra fails resolving ${oc.env:AF2_DIR} before anything runs.",
    )
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="A raw Hydra override, passed through verbatim with a '++' prefix "
        "(e.g. --set generation.dataloader.batch_size=8). Repeatable. This is the "
        "escape hatch for everything LPC exposes that has no flag here.",
    )


def build_laproteina_manifest(
    ctx: ManifestCtx[LaProteinaArgs],
) -> list[tuple[str, ...]]:
    args = ctx.args

    if args.table is not None:
        raise ValueError(
            "laproteina is a root-only create: LPC resolves its target by name from "
            "its own registry, not from a table column, so there is nothing for "
            "--table/-t to supply. Drop -t and pass --task-name."
        )
    if not args.task_name:
        raise ValueError(
            "--task-name is required: it names the target in LPC's "
            "configs/targets/targets_dict.yaml. Run `complexa target list` to see "
            "what is registered."
        )
    if args.num_jobs < 1:
        raise ValueError(f"--num-jobs must be >= 1, got {args.num_jobs}")
    if args.nsamples is not None and args.nsamples < 1:
        raise ValueError(f"--nsamples must be >= 1, got {args.nsamples}")

    prefix = args.design_prefix or args.task_name

    # Hydra overrides shared by every job. Per-job values (seed, run token) are
    # appended in the loop, because they differ per task.
    base_overrides: list[str] = [f"generation.task_name={args.task_name}"]

    # gen_njobs MUST match the number of tasks. generate.py reads it as
    # `njobs = cfg.get("gen_njobs", 1)` and hands it to split_by_job(cfg, job_id,
    # njobs), which slices the sample budget: job j takes
    # ceil(nsamples/njobs) samples starting at j*that. Leave gen_njobs at the
    # config's default of 1 while submitting job_id>0 and the slice check
    #     nsamples_per_split * job_id >= nsamples
    # is true immediately, so the job logs "Job id N get 0 samples" and calls
    # sys.exit(0) -- exit 0, no structures, nothing to collect.
    base_overrides.append(f"gen_njobs={args.num_jobs}")

    if args.nsamples is not None:
        # The count lives on the nres sampler, not on dataset itself. Overriding
        # `dataset.nsamples` with ++ would quietly ADD a dead key and leave the
        # real value at its default.
        base_overrides.append(
            f"generation.dataloader.dataset.nres.nsamples={args.nsamples}"
        )

    # The shipped binder config turns on inference-time search (`best-of-n`) and
    # an AF2 reward model (`af_params_dir: ${oc.env:AF2_DIR}`). Neither can work
    # in the --minimal image: ColabDesign and JAX are not installed, and AF2_DIR
    # is unset, so Hydra fails at config resolution before the model ever loads.
    # Turn both off by default. `initialize_reward_model` handles a null cleanly
    # (logs "No reward model configured" and returns None).
    #
    # These go BEFORE --set so a caller running a full image can switch them back
    # on; with Hydra's ++, the later override wins.
    base_overrides.append(f"generation.search.algorithm={args.search}")
    if not args.rewards:
        base_overrides.append("generation.reward_model=null")

    base_overrides.extend(args.set)

    samples_root = volume_path(ctx.out_dir) / SAMPLES_DIRNAME

    # The stage contract, recorded for collect the same way bindcraft2 records its
    # own: the run decides, collect reads. Without this the collector would have to
    # guess where the task parked LPC's output.
    ctx.write_meta(
        task_name=args.task_name,
        samples_root=str(samples_root),
        nsamples=args.nsamples,
        config=args.config,
    )

    manifest_rows: list[tuple[str, ...]] = []
    group_names: list[str] = []
    for job_id in range(args.num_jobs):
        name = f"{prefix}_j{job_id}"
        overrides = list(base_overrides)
        if args.seed is not None:
            overrides.append(f"seed={args.seed + job_id}")
        # run_name is what makes LPC's output directory unique per task, and it is
        # how the task script finds that directory again afterwards. It must not
        # collide between jobs of the same run.
        overrides.append(f"run_name={name}")

        manifest_rows.append(
            (
                name,
                args.config,
                str(job_id),
                " ".join(f"++{o}" for o in overrides),
                str(samples_root / name),
            )
        )
        group_names.append(name)

    # A root run has no parent table for collect to iterate; record the group names
    # so `sapia collect` can find their outputs and rebuild their rows.
    ctx.write_meta(root_designs=group_names)

    if args.nsamples is not None:
        per_job = -(-args.nsamples // args.num_jobs)  # ceil, as split_by_job does
        detail = f" -> {args.nsamples} designs total (~{per_job}/job)"
        if args.num_jobs > args.nsamples:
            print(
                f"WARNING: --num-jobs ({args.num_jobs}) exceeds --nsamples "
                f"({args.nsamples}). LPC gives each job ceil(nsamples/njobs) "
                f"samples and exits 0 on any job whose slice starts past the end, "
                f"so some tasks will produce nothing and be reported as failures."
            )
    else:
        detail = " (nsamples from the config)"
    print(
        f"Submitting {args.num_jobs} LPC generation job(s) against "
        f"'{args.task_name}'{detail}"
    )
    return manifest_rows
