#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) scoring each design's backbone as a POLY-VALINE
binder against a fixed target, with PyRosetta only.

PREMISE -- sequence-blind backbone triage. This tool answers one question:
*does this backbone present a large, well-fitted surface to a fixed target, before
any sequence exists?* It erases the binder's sequence to poly-valine, packs only the
binder's side chains against a completely frozen target, and measures interface
geometry (contact molecular surface, dSASA, interface residue counts, shape
complementarity). Nothing is minimised and nothing is relaxed, anywhere: speed is the
point, and a moved backbone would no longer be the backbone that was triaged.

Two properties of the per-design work that the flags here control, and that anyone
editing those flags needs to know:

* **The result is deterministic per design name.** Rosetta's packer is stochastic and
  its RNG state carries ACROSS the designs of a task, so the worker reseeds from
  ``sha256(<design name>)`` before each design. One input file therefore gives one
  answer regardless of its position in the batch, of ``--designs-per-task`` or of row
  order, and ``--designs-per-task`` is a scheduling knob only. The seed is collected as
  ``seed``. Two DIFFERENT names on the same file get different seeds and so may differ
  in the last decimals.
* **Residues can be deleted before scoring.** Chains named in neither
  ``--binder-chains`` nor ``--target-chains`` are removed (InterfaceAnalyzerMover takes
  exactly two partners), and so is any NON-POLYMER residue named in
  ``--exclude-resnames``. Both counts are collected (``n_res_dropped_chain`` /
  ``n_res_dropped_resname``), and ``n_res_binder`` / ``n_res_target`` are the counts
  that were SCORED, i.e. after those deletions. A polymer residue may not be excluded
  (it would leave a hole and change every remaining residue's surface), and a requested
  resname that matches nothing is an error.

IT IS NOT an interface-quality score for a designed sequence, and it does NOT replace
the ``cms`` tool. Its numbers are **not comparable** to ``cms`` (different
implementation -- Rosetta's ContactMolecularSurface filter vs cms-cuda -- and polyV vs
the as-built sequence) and **not comparable** to ``pyrosetta``'s ``if_dSASA``
(different pose, no relax, no separated repack). Rank polyv columns only against other
polyv columns from the same run.

Its numbers are meaningless when:
  * the input is not a binder/target complex in a common frame (one chain, or a
    binder that is not docked onto the target);
  * the chain letters passed do not mean what you think -- read ``polyv_chains_found``
    before anything else;
  * you want chemistry. There is no energy, no dG, no hydrogen bond, no electrostatics
    here, by design. Poly-valine has no polar side chains at all.

``action: update`` -- a geometric score is a property of a backbone that already
exists, so it annotates the table in place. ``-t`` is required. Run it on whichever
structure column holds the backbone (``-i``), and label each run with ``-l`` so the
columns do not collide.

``default_input_column`` is the ``"not applicable"`` sentinel (as in cms, usalign and
chainsel): there is no honest default structure column, and a wrong one produces a
plausible wrong number rather than a failure, so the builder refuses the run unless
``-i`` names a column the table has.

Several designs are packed into one task (``--designs-per-task``, default 100):
``pyrosetta.init()`` costs ~3-4 s per container, so a task per design would be mostly
startup. ``n_tasks`` is therefore the number of CHUNKS, not the number of designs.

Usage:
    sapia run polyv outputs/20260928_100140_7ojg_binder \
        --table table0 \
        --input-column rfdiffusion3_path \
        --binder-chains B --target-chains A,C \
        --dir-label bb
"""

from argparse import ArgumentParser
from pathlib import Path
from typing import cast

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path

# The sentinel default_input_column (see the module docstring).
NO_DEFAULT_COLUMN = "not applicable"
EXTRA_ROTAMERS = ("none", "ex1", "ex1ex2")
# polyv renders a BLANK chain ID as "_" in chains_found; it is not selectable (it is
# also the separator in InterfaceAnalyzerMover's partner string).
BLANK_CHAIN = "_"

# Residue names refused in --exclude-resnames: deleting a POLYMER residue leaves a hole
# in the chain, exposing backbone and neighbour atoms it previously occluded, so every
# REMAINING residue's SASA, CMS and SC change for reasons unrelated to the interface.
# This list is fast submit-time feedback only and cannot be exhaustive -- the worker's
# authoritative test is Rosetta's own ``residue.is_polymer()``.
POLYMER_RESNAMES = frozenset(
    # the 20 canonical amino acids
    """ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR
       VAL""".split()
    # protonation / tautomer / termini variants seen in prepared structures
    + "HID HIE HIP HSD HSE HSP CYX CYM ASH GLH LYN".split()
    # modified polymer residues, and the unknown-polymer placeholder
    + "MSE SEP TPO PTR UNK".split()
    # nucleotides
    + "A C G U T DA DC DG DT RA RC RG RU ADE CYT GUA THY URA".split()
)


class PolyvArgs(CommonArgs):
    binder_chains: str
    target_chains: str
    keep_gly_pro: bool
    extra_rotamers: str
    distance_weight: float
    exclude_resnames: str
    designs_per_task: int


def add_run_polyv_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--binder-chains",
        type=str,
        required=True,
        help="The binder's chain IDs, comma-joined (e.g. 'B'). These are the chains "
        "whose sequence is erased to poly-valine and whose side chains are packed. A "
        "chain named here but ABSENT from a design is an error for that design, never "
        "a smaller selection.",
    )
    parser.add_argument(
        "--target-chains",
        type=str,
        required=True,
        help="The target's chain IDs, comma-joined (e.g. 'A,C' for a multimeric "
        "target). The target keeps its sequence and is completely frozen: no "
        "mutation, no repacking, no minimisation. Must be disjoint from "
        "--binder-chains (checked at submit time).",
    )
    parser.add_argument(
        "--keep-gly-pro",
        action="store_true",
        help="Keep GLY and PRO on the binder instead of mutating them to VAL. Off by "
        "default: the default polyV pose is a uniform 'what could this backbone "
        "present' probe. This does NOT weaken the trust contract -- the check becomes "
        "n_mutated + n_kept_gly_pro == n_res_binder, recounted from the pose and "
        "exactly as strong, so a residue that failed to mutate still cannot hide in "
        "the gap.",
    )
    parser.add_argument(
        "--extra-rotamers",
        type=str,
        choices=EXTRA_ROTAMERS,
        default="ex1",
        help="Extra rotamer sampling for the binder repack (default 'ex1'). 'none' is "
        "fastest and coarsest; 'ex1ex2' is finer and several times slower. Valine has "
        "a single chi, so ex2 buys little.",
    )
    parser.add_argument(
        "--distance-weight",
        type=float,
        default=0.5,
        help="distance_weight of Rosetta's ContactMolecularSurface filter (default "
        "0.5, the Baker-lab convention). Changing it makes the cms_* columns "
        "non-comparable with every run that used another value.",
    )
    parser.add_argument(
        "--exclude-resnames",
        type=str,
        default="",
        help="NON-POLYMER residue names dropped before scoring, comma-joined (e.g. "
        "'HEM,ZN' -- a ligand, ion, lipid or cofactor you do not want counted). "
        "Default: none, and polyv then scores EVERYTHING Rosetta accepted on the named "
        "chains, heteroatoms included. A polymer residue (any amino acid or "
        "nucleotide) is REFUSED: deleting one leaves a hole that changes every "
        "remaining residue's surface. A name that matches no residue is an error, not "
        "a smaller selection. Removals are collected as n_res_dropped_resname.",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=100,
        help="Designs scored per task (default 100). pyrosetta.init() costs ~3-4 s per "
        "container, so one task per design would be mostly startup. Lower it for more "
        "parallelism on very large tables.",
    )


def _split(value: str) -> list[str]:
    """Comma-joined list -> stripped, non-empty tokens."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


def build_polyv_manifest(ctx: ManifestCtx[PolyvArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool: no caller needs to pass -g 0

    binder = _split(ctx.args.binder_chains)
    target = _split(ctx.args.target_chains)
    if not binder or not target:
        raise ValueError(
            "--binder-chains and --target-chains must both name at least one chain."
        )
    overlap = sorted(set(binder) & set(target))
    if overlap:
        raise ValueError(
            f"chains {overlap} are in both --binder-chains and --target-chains; an "
            f"interface needs two disjoint sides."
        )
    if ctx.args.designs_per_task < 1:
        raise ValueError("--designs-per-task must be >= 1.")
    if not ctx.args.distance_weight > 0:
        raise ValueError("--distance-weight must be > 0.")
    if ctx.args.extra_rotamers not in EXTRA_ROTAMERS:
        raise ValueError(f"--extra-rotamers must be one of {EXTRA_ROTAMERS}.")
    excluded = [r.strip().upper() for r in _split(ctx.args.exclude_resnames)]
    refused = sorted(set(excluded) & POLYMER_RESNAMES)
    if refused:
        raise ValueError(
            f"--exclude-resnames names polymer residue {','.join(refused)}; the flag "
            f"is for non-polymer heteroatoms (a ligand, ion, lipid or cofactor). "
            f"Deleting an interior polymer residue exposes backbone its neighbours "
            f"occluded, so every remaining residue's SASA, CMS and SC would change for "
            f"reasons unrelated to the interface."
        )
    if BLANK_CHAIN in binder + target:
        raise ValueError(
            f"{BLANK_CHAIN!r} is polyv's rendering of a BLANK chain ID in the "
            f"chains_found column, not a selectable chain, and it is also the "
            f"separator in InterfaceAnalyzerMover's partner string. Give the "
            f"structures real chain IDs (see the chainsel tool) instead."
        )

    column = ctx.args.input_column
    if column == NO_DEFAULT_COLUMN or column not in ctx.df.columns:
        available = ", ".join(
            str(c) for c in ctx.df.columns if str(c).endswith("_path")
        )
        raise ValueError(
            f"polyv has no default input column: pass -i/--input-column with the "
            f"structure column holding the backbone to triage (got {column!r}, which "
            f"table '{ctx.args.table}' does not have). Structure columns available: "
            f"{available or '(none)'}."
        )

    ready = ctx.ready
    members: list[tuple[str, str]] = []
    for name in ready.index:
        name = cast(str, name)
        src = Path(str(ready.at[name, column]))
        if not src.exists():
            print(f"{name}: MISSING {src} (skipping)")
            continue
        # No CIF->PDB staging: PyRosetta reads .pdb and .cif, and nothing may be
        # written next to the input file.
        members.append((name, str(volume_path(src))))

    # One sub-manifest per task, and a top-level row per task pointing at it
    # (the cms pattern).
    tasks_dir = ctx.out_dir / "polyv_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    x = ctx.args.designs_per_task
    manifest_rows: list[tuple[str, ...]] = []
    for t, i in enumerate(range(0, len(members), x)):
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            for name, src in members[i : i + x]:
                f.write(f"{name}\t{src}\n")
        manifest_rows.append(
            (
                str(volume_path(task_file)),
                ",".join(binder),
                ",".join(target),
                "keep" if ctx.args.keep_gly_pro else "mutate",
                ctx.args.extra_rotamers,
                # May be empty, so never last.
                ",".join(excluded),
                # Kept last: always non-empty.
                str(ctx.args.distance_weight),
            )
        )

    return manifest_rows
