#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) that rebuilds each design's multi-chain complex
sequence.

ProteinMPNN redesigns only the binder chain, so ``proteinmpnn_sequence`` is the binder
monomer alone and a predictor would fold it in isolation. ``mkcomplex`` writes a NEW
sequence column that puts the fixed target chain(s) back around the design's own
sequence, in the ``/``-separated chainbreak form the predictors parse (see
``prosapia.utils.chains``). It is pure string work, so the whole per-design step lives
in ``mkcomplex.sh``.

Fixed chains come as literals (``--prepend-seqs``/``--append-seqs``) or from a FASTA on
the runs volume (``--prepend-fasta``/``--append-fasta``), one of the two per side.
``--repeat`` emits N copies of every fixed chain and ``--design-copies`` N copies of
THIS ROW's sequence. The layout is prepend chains, design copies, append chains, and
``design_chain_index`` is the index of the FIRST design copy. Flags and traps:
``.claude/skills/mkcomplex/SKILL.md``.

Usage:
    sapia run mkcomplex outputs/20260928_100140_7ojg_binder \
        --table table1 \
        --prepend-seqs MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ \
        --repeat 2

    sapia run mkcomplex outputs/20260928_100140_7ojg_binder \
        --table table1 -i trifuse_sequence \
        --prepend-seqs T,T,T --design-copies 3
"""

from argparse import ArgumentParser
from pathlib import Path
from typing import cast

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path


class MkcomplexArgs(CommonArgs):
    prepend_seqs: str
    prepend_fasta: str
    append_seqs: str
    append_fasta: str
    separator: str
    repeat: int
    design_copies: int


def add_run_mkcomplex_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--prepend-seqs",
        type=str,
        default="",
        help="Literal chain sequences placed BEFORE the design's own sequence, "
        "comma-joined and ordered, e.g. 'SEQA,SEQB'. One-letter codes, uppercase; "
        "a chain with any other character is reported as an error for that design "
        "rather than silently assembled. Mutually exclusive with --prepend-fasta.",
    )
    parser.add_argument(
        "--prepend-fasta",
        type=str,
        default="",
        help="FASTA file whose records are the fixed chains placed BEFORE the "
        "design's own sequence, in file order. Must live on the runs volume under "
        "--executor modal (the path is resolved into the task container's /runs). "
        "Mutually exclusive with --prepend-seqs.",
    )
    parser.add_argument(
        "--append-seqs",
        type=str,
        default="",
        help="Literal chain sequences placed AFTER the design's own sequence, "
        "comma-joined and ordered. Mutually exclusive with --append-fasta.",
    )
    parser.add_argument(
        "--append-fasta",
        type=str,
        default="",
        help="FASTA file whose records are the fixed chains placed AFTER the "
        "design's own sequence, in file order. Mutually exclusive with "
        "--append-seqs.",
    )
    parser.add_argument(
        "--separator",
        type=str,
        default="/",
        help="Chainbreak character joining the chains. Default '/', which is what "
        "boltz/alphafold3/colabfold parse into chains. An input sequence that "
        "already contains it is refused (it is already a multi-chain string).",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=1,
        help="How many consecutive copies of EACH fixed chain to emit (default 1 = "
        "one copy each). It applies to both sides and to both ways of supplying "
        "them, and NEVER to the design's own sequence (that is --design-copies, "
        "which defaults to 1): '--prepend-seqs A,B --repeat 2' gives "
        "'A/A/B/B/<design>', so "
        "'--prepend-seqs X --repeat 2' is exactly '--prepend-seqs X,X'. Chains "
        "identical in sequence are later collapsed into one entity with N chain "
        "ids by the predictor, so --repeat 2 on a homodimer target asks for a "
        "2-chain entity, not two entities.",
    )
    parser.add_argument(
        "--design-copies",
        type=int,
        default=1,
        help="How many consecutive copies of THE DESIGN'S OWN sequence to emit "
        "(default 1 = the single copy mkcomplex has always emitted). Unlike "
        "--prepend-seqs/--append-seqs, which are one literal for the whole run, "
        "this copies the per-row value of --input-column, so every row gets its "
        "own chain repeated: '--prepend-seqs T,T,T --design-copies 3' gives "
        "'T/T/T/F/F/F' with F this row's sequence. The layout stays prepend "
        "chains, then the design copies, then append chains, and "
        "mkcomplex_design_chain_index remains the index of the FIRST copy. "
        "Copies are identical by construction, so the predictor collapses them "
        "into ONE entity with N chain ids -- same caveat as --repeat.",
    )


def _resolve_fasta(flag: str, value: str) -> str:
    """Validate a FASTA flag at submit time and return its in-container path."""
    path = Path(value)
    if not path.is_file():
        raise FileNotFoundError(
            f"{flag} {path} does not exist. Under --executor modal it must be a "
            f"path on the runs volume (the workstation's /runs)."
        )
    # Must resolve inside a task container, not just in the submitting shell.
    return str(volume_path(path))


def _side(flag: str, seqs: str, fasta: str) -> tuple[str, str]:
    """One side's (literals, fasta_path) fields, refusing both sources at once."""
    if seqs and fasta:
        raise ValueError(
            f"--{flag}-seqs and --{flag}-fasta are mutually exclusive: give the "
            f"fixed chains either as literals or as a FASTA, not both."
        )
    return seqs, _resolve_fasta(f"--{flag}-fasta", fasta) if fasta else ""


def build_mkcomplex_manifest(ctx: ManifestCtx[MkcomplexArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool

    prepend_seqs, prepend_fasta = _side(
        "prepend", ctx.args.prepend_seqs, ctx.args.prepend_fasta
    )
    append_seqs, append_fasta = _side(
        "append", ctx.args.append_seqs, ctx.args.append_fasta
    )
    if not (prepend_seqs or prepend_fasta or append_seqs or append_fasta):
        raise ValueError(
            "no fixed chains given: pass --prepend-seqs/--prepend-fasta (chains "
            "before the design) and/or --append-seqs/--append-fasta (after it). "
            "Without them mkcomplex would only copy the input column."
        )

    if not ctx.args.separator:
        raise ValueError("--separator must not be empty (default '/').")
    if ctx.args.repeat < 1:
        raise ValueError(f"--repeat must be >= 1, got {ctx.args.repeat}.")
    if ctx.args.design_copies < 1:
        raise ValueError(
            f"--design-copies must be >= 1, got {ctx.args.design_copies}. "
            f"Zero copies would drop the design's own chain and assemble a "
            f"complex that contains no design at all."
        )

    ready = ctx.ready

    manifest_rows: list[tuple[str, ...]] = []
    for name in ready.index:
        name = cast(str, name)
        sequence = str(ready.at[name, ctx.args.input_column]).strip()
        manifest_rows.append(
            (
                name,
                sequence,
                prepend_seqs,
                prepend_fasta,
                append_seqs,
                append_fasta,
                # Last: always non-empty, so the line never ends on an empty field.
                ctx.args.separator,
                str(ctx.args.repeat),
                str(ctx.args.design_copies),
            )
        )

    return manifest_rows
