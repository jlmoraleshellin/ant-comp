#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) measuring each rfd3 design's motif CA RMSD.

Each task runs motifrmsd_worker.py on one design: it reads rfd3's own per-design
sidecar JSON (``diffused_index_map``), superposes the design's motif CA atoms onto
the reference's motif CA atoms with a Kabsch fit, and reports the RMSD of those
atoms after that fit -- the literal "motif RMSD after alignment on the motif"
gate. Use ``sapia collect motifrmsd`` to merge the metrics back into the table.

SCOPE. This tool only works on **rfdiffusion3 motif-scaffolding output**, because
the residue correspondence it measures over comes from rfd3's sidecar JSON and
from nowhere else. Pointed at any other structure column it errors per design
rather than falling back to a whole-structure alignment -- that comparison is
``usalign``'s, a different premise.

FINDING THE SIDECAR. The table's ``<leaf>_path`` column points at a **converted
PDB** under ``<run_dir>/.cif_to_pdb/``, not at the ``.cif.gz`` the sidecar sits
beside. prosapia's ``ensure_pdb`` writes a ``<cached>.src`` record next to every
converted file holding the absolute path of the source it was made from, and this
builder follows that record to the ``.cif.gz`` and then to its ``.json`` sibling.
That route is used rather than globbing ``results_shard_*/`` by file stem because
**the stem is not unique**: two rfd3 runs in one run_dir (two contigs, two
``--table-label``s) produce files with identical names in different output dirs,
distinguished in the cache only by the digest in the cached file's name. A
stem-glob would pick one of the two at random. (Verified on
``outputs/20261008_085450_5M2J_vhh2_motif_graft``: both
``...diff_0_model_0.cif.8ecf9c29.src`` and ``...diff_0_model_0.cif.ac344f9c.src``
exist, pointing into ``table0_seg2`` and ``table0_seg3`` respectively.)

Resolution happens here, at submit time, not in the task: the path rewriting a
Modal task container needs (``volume_path``) lives in prosapia, which the task's
slim image does not carry. A design whose sidecar cannot be resolved is still
submitted, with an empty sidecar field, so the task writes a per-design
``error:`` row rather than vanishing from the table.

Usage:
    sapia run motifrmsd outputs/20261008_085450_5M2J_vhh2_motif_graft \
        --table table0_seg2 \
        --input-column rfdiffusion3_seg2_path \
        --ref-structure inputs/5M2J_D_vhh.pdb
"""

import re
from argparse import ArgumentParser
from pathlib import Path
from typing import cast

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path


class MotifRmsdArgs(CommonArgs):
    ref_structure: str
    allow_partial: bool
    motif_chain: str


# A chain id as ``diffused_index_map``'s keys spell it: letters only, matching the
# worker's RES_KEY. Anything else (a digit, a comma-list, a residue range) can
# never match a key, so it is refused at submit time rather than per design.
MOTIF_CHAIN_RE = re.compile(r"^[A-Za-z]+$")


def add_run_motifrmsd_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--ref-structure",
        type=str,
        default="",
        help="Reference structure the motif was lifted from, in ITS OWN chain "
        "ids and numbering (the left-hand side of rfd3's diffused_index_map, "
        "e.g. chain D of 5M2J). Empty (the default) falls back per design to the "
        "sidecar's own specification.input -- which is what rfd3 actually "
        "diffused. When given, the flag WINS, and its content is compared "
        "(sha256) against the sidecar's input: see motifrmsd_ref_matches_spec. "
        "Must live on the runs volume under --executor modal.",
    )
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Measure the matched subset when some residues of "
        "diffused_index_map are absent from the design or the reference, with a "
        "'warn:' status (so the row is still excluded by == 'OK'). OFF by "
        "default: a count mismatch means the correspondence did not resolve, "
        "and an RMSD over a quietly smaller atom set is the failure this tool "
        "exists to avoid.",
    )
    parser.add_argument(
        "--motif-chain",
        type=str,
        default="",
        help="Restrict design_fixed_positions to the motif residues lifted from "
        "ONE reference chain, e.g. 'D'. A REFERENCE-frame chain letter -- the "
        "left/key side of diffused_index_map, the same frame --ref-structure and "
        "the contig use, NOT the design's output chain. Affects "
        "design_fixed_positions (and the design-side letter reported in "
        "motifrmsd_fixed_positions_chain) and NOTHING else: the RMSD, the "
        "segments and every other column keep measuring the WHOLE mapping. Needed "
        "for binder scaffolding against a held-fixed target, where rfd3 records "
        "the entire target in diffused_index_map too, so the motif spans several "
        "design chains and the unrestricted column is always NA. Empty (the "
        "default) is the unrestricted behaviour, unchanged.",
    )


def _json_beside(path: Path) -> Path:
    """The ``.json`` sibling of an rfd3 output, handling the ``.cif.gz`` double
    suffix (``x.cif.gz`` -> ``x.json``)."""
    name = path.name
    for suffix in (".cif.gz", ".pdb.gz"):
        if name.endswith(suffix):
            return path.with_name(name[: -len(suffix)] + ".json")
    return path.with_suffix(".json")


def resolve_sidecar(design: Path) -> tuple[Path | None, str]:
    """The rfd3 sidecar JSON for a design, or ``(None, reason)``.

    Two routes, in order, and no guessing beyond them:

    1. ``<design>.src`` -- the record prosapia's CIF->PDB cache writes beside
       every converted file, naming the exact source ``.cif.gz``. Authoritative,
       and the only route that is unambiguous when two runs share a file stem.
    2. the ``.json`` sibling of the design path itself, for a column that already
       points straight at an rfd3 output (a ``.cif.gz`` or a PDB in the results
       dir) rather than at the conversion cache.

    Never raises for one bad artifact: an empty, undecodable or path-less
    ``.src`` record is returned as ``(None, reason)``, so that design is still
    submitted and collects as a per-design ``error:`` row instead of aborting
    the whole submit. A present-but-broken ``.src`` does NOT fall through to
    route 2 -- the record exists, so the cache owns this design, and a sibling
    guess is exactly what route 1 exists to avoid.
    """
    src_record = design.with_suffix(".src")
    if src_record.is_file():
        try:
            text = src_record.read_text().strip()
        except (OSError, UnicodeDecodeError) as e:
            return None, f"cache record {src_record} unreadable: {e}"
        if not text:
            return None, f"cache record {src_record} is empty"
        source = Path(text)
        try:
            sidecar = _json_beside(source)
        except ValueError as e:
            # e.g. a record holding "/" or "." -- no file name to take a sibling of.
            return None, f"cache record {src_record} holds no file path ({text!r}): {e}"
        if sidecar.is_file():
            return sidecar, ""
        return None, f"no sidecar JSON beside the cached source {source}"

    sidecar = _json_beside(design)
    if sidecar.is_file():
        return sidecar, ""
    return (
        None,
        f"no CIF->PDB cache record ({src_record.name}) and no sidecar JSON "
        f"beside {design.name}",
    )


def build_motifrmsd_manifest(ctx: ManifestCtx[MotifRmsdArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool

    ref_path = ""
    if ctx.args.ref_structure:
        ref = Path(ctx.args.ref_structure)
        if not ref.exists():
            raise FileNotFoundError(
                f"--ref-structure {ref} does not exist. Under --executor modal it "
                f"must be a path on the runs volume (the workstation's /runs)."
            )
        # Must resolve inside a task container, not just in the submitting shell.
        ref_path = str(volume_path(ref))

    # A malformed letter can never match a map key, so it is a typo, not a
    # measurement failure: fail here, before anything is queued. A WELL-FORMED
    # letter that simply is not in a design's map is per-design data (the worker
    # NAs that design's column with a reason), because which chains a map names
    # is a property of each design, not of the flag.
    motif_chain = ctx.args.motif_chain.strip()
    if motif_chain and not MOTIF_CHAIN_RE.match(motif_chain):
        raise ValueError(
            f"--motif-chain {ctx.args.motif_chain!r} is not a chain id. Give one "
            f"REFERENCE-frame chain id as diffused_index_map spells it (letters "
            f"only, e.g. 'D') -- not a residue range, a list or a design chain."
        )

    ready = ctx.ready

    manifest_rows: list[tuple[str, ...]] = []
    n_unresolved = 0
    for name in ready.index:
        name = cast(str, name)
        design = Path(str(ready.at[name, ctx.args.input_column]))
        if not design.exists():
            print(f"{name}: MISSING {design} (skipping)")
            continue
        try:
            sidecar, reason = resolve_sidecar(design)
        except (OSError, ValueError) as e:
            # Belt and braces: resolve_sidecar reports every failure it knows of
            # as data, but one unforeseen bad artifact among thousands must still
            # cost one row, never the whole submit.
            sidecar, reason = None, f"sidecar resolution raised {type(e).__name__}: {e}"
        if sidecar is None:
            # Submitted anyway: the task records 'error: ...' as data, so the
            # design appears in the table as a failure rather than as a gap.
            n_unresolved += 1
            print(f"{name}: no rfd3 sidecar JSON -- {reason}")
        manifest_rows.append(
            (
                name,
                str(volume_path(design)),
                "" if sidecar is None else str(volume_path(sidecar)),
                ref_path,
                motif_chain,
                # Kept last: always non-empty, so the manifest line never ends
                # on an empty field.
                "yes" if ctx.args.allow_partial else "no",
            )
        )

    if n_unresolved:
        print(
            f"WARNING: {n_unresolved} of {len(manifest_rows)} design(s) have no "
            f"resolvable rfd3 sidecar JSON. They will collect as "
            f"'error: rfd3 sidecar JSON not found'. If that is all of them, the "
            f"input column is probably not rfd3 motif-scaffolding output."
        )

    return manifest_rows
