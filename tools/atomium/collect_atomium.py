"""Collect AtomiUM FASTA outputs into a child table, one row per designed sequence.

Lineage is ``<output_dir>/grp_<g>/seqs/<stem>.fa``: the submitter symlinks each input
as ``<design_name>.pdb``, so the FASTA stem IS the parent row (stamped ``parent_name``).

The FASTA is not ProteinMPNN's, so this does not share its parser::

    >_init_sequence, chains=D                  <- echo of the native input, skipped
    MKT...
    >gen_sample_0, temperature=0.1, seq_rec=0.43
    MRE...

The echo is recognised by the first comma-separated field (case-insensitive), so the
``chains=`` suffix and any further field are tolerated. There is no score column, and
``seq_rec`` is similarity to the INPUT sequence, not quality. ``sample`` recurs once
per ``--sampling-temp`` value, so rows are keyed by entry index (``<parent>_a<i>``)
and (``sample``, ``temperature``) identify a draw. Details:
``.claude/skills/atomium/SKILL.md``.

Usage:
    sapia collect atomium outputs/<run> --table <the table the run reserved>
"""

from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import pandas as pd

from prosapia.core import (
    Collected,
    CollectCtx,
    CollectEach,
    DesignCtx,
)

# The header AtomiUM gives the native sequence it echoes first in every FASTA.
INIT_HEADER = "_init_sequence"

# Trailing ``key=value`` fields on a generated header.
HEADER_FIELDS: Dict[str, type] = {
    "temperature": float,
    "seq_rec": float,
}


def parse_atomium_header(header: str) -> Dict[str, Any]:
    """Parse ``gen_sample_0, temperature=0.1, seq_rec=0.43``; bad fields become NA."""
    out: Dict[str, Any] = {"sample": pd.NA}
    out.update({k: pd.NA for k in HEADER_FIELDS})
    for part in header.split(","):
        part = part.strip()
        if part.startswith("gen_sample_"):
            try:
                out["sample"] = int(part[len("gen_sample_") :])
            except ValueError:
                pass
            continue
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        key, value = key.strip(), value.strip()
        if key in HEADER_FIELDS:
            try:
                out[key] = HEADER_FIELDS[key](value)
            except ValueError:
                out[key] = pd.NA
    return out


def parse_fasta(fasta_path: Path) -> List[Tuple[str, str]]:
    entries: List[Tuple[str, str]] = []
    header = ""
    seq_lines: List[str] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    entries.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header:
        entries.append((header, "".join(seq_lines)))
    return entries


def collect_atomium(ctx: CollectCtx) -> CollectEach:
    """Per-parent collector: one child row ``<parent>_a<i>`` per sampled sequence."""
    # The staged input is symlinked as <design>.pdb: the FASTA stem is the parent row.
    fasta_by_parent: Dict[str, Path] = {}
    for subdir in sorted(p for p in ctx.out_dir.iterdir() if p.is_dir()):
        if subdir.name in ("atomium_logs", "atomium_tasks"):
            continue
        seqs_dir = subdir / "seqs"
        if not seqs_dir.is_dir():
            continue
        for fasta_path in seqs_dir.glob("*.fa"):
            fasta_by_parent[fasta_path.stem] = fasta_path

    def one(d: DesignCtx) -> Iterable[Collected]:
        fasta_path = fasta_by_parent.get(d.name)
        if fasta_path is None:
            print(f"{d.name}: no FASTA found (skipping)")
            return

        entries = parse_fasta(fasta_path)
        n_entries = len(entries)
        n_designed = 0

        for i, (header, sequence) in enumerate(entries):
            # Skip AtomiUM's echo of the native input sequence.
            if header.split(",", 1)[0].strip().casefold() == INIT_HEADER.casefold():
                continue
            n_designed += 1

            # Bare column names; the driver leaf-prefixes them.
            data: Dict[str, Any] = {"iteration": i, "sequence": sequence}
            data.update(parse_atomium_header(header))
            yield Collected(
                name=f"{d.name}_a{i}",
                parent=d.name,
                path=fasta_path,
                data=data,
            )

        if n_designed == 0:
            print(
                f"{d.name}: FASTA has {n_entries} entr(y/ies), "
                "none designed (skipping)"
            )

    return one
