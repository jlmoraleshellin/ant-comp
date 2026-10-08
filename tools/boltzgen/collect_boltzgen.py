#!/usr/bin/env python3
"""
Collect BoltzGen nanobody designs into the (child) design table.

``run_boltzgen.py`` points each target's ``--output`` at ``<out_dir>/runs/<name>/``,
where BoltzGen writes (``README.md`` "Pipeline output", verified against the pinned
SHA):

    runs/<name>/intermediate_designs_inverse_folded/
        <id>.cif                     the inverse-folded design complex
        <id>.npz                     design_mask + source_context for that design
        refold_cif/<id>.cif          the REFOLDED complex (target + binder)
        aggregate_metrics_analyze.csv

``refold_cif/<id>.cif`` is what every analysis metric was computed on, so it is what
``boltzgen_path`` points at (converted to PDB; the mmCIF stays as
``boltzgen_cif_path``). **There is no binder-alone structure.** The ``design_folding``
step is disabled for nanobody binders and ``refold_design_cif/`` is never written
(upstream issue #148). A self-consistency check therefore needs ``chainsel`` to cut
the nanobody chain out first.

There is also no ``final_ranked_designs/``: ``filtering`` cannot run without ESMFold2
scores and cannot be configured out of needing them, so this tool never runs it.
Every design BoltzGen produced is collected, ranked by nothing. Rank with ``-f``.

The join is BY NAME, never by row order
---------------------------------------
``file_name == id + ".cif"`` (``analyze.py:586-588,628-629``) and ``<id>.npz`` sits
beside the design CIF, so one ``id`` addresses the CSV row, the design CIF, the npz
and the refolded complex. Nothing here indexes into the CSV positionally.

The mapping problem this collector exists to make visible
---------------------------------------------------------
Three different naming systems have to be reconciled before a single sequence can be
believed, and getting it wrong yields a plausible, wrong answer every time:

1. the CSV's sequence columns are keyed by **feature asym id**, not PDB chain label
   (``docs/sequence-analysis.md``);
2. ``source_context[...]["output_chain"]`` is BoltzGen's **internal** chain name;
3. the written mmCIF label is ``re.sub(r"\\d+", "", internal_name)`` and, **on a
   collision, an arbitrary letter from a reversed-alphabet pool**
   (``data/write/mmcif.py:63-80``). So the label in the file need not be either of
   the other two.

So the binder chain is resolved from physical evidence and then cross-checked three
ways, and the result is reported:

* the chain is the one whose atoms carry a **design-mask B-factor of 1** in the
  design CIF (``model/validation/design.py``/``writer.py``: ``structure.atoms
  ["bfactor"] = atom_design_mask``). Exactly one such chain, or the design is an
  error. Only that file is ever read for the mask: the refolded complex carries
  **pLDDT** in the same column, so ``read_chains`` takes ``design_mask`` explicitly
  and returns ``designed=None`` when it is false.
* its designed-residue count must agree with the npz ``design_mask`` **and** with the
  CSV's own ``num_design``. A count mismatch is an error, never a trimmed fit; so is
  an ABSENT ``num_design``, because the three-witness guarantee would silently become
  a two-witness one.
* the design chain and the refolded chain must be identical **residue for residue**,
  not merely the same length -- the designed positions are read off the former and
  reported against the latter, so a length-preserving shift would invert
  ``boltzgen_fixed_resnums`` (fixing the CDRs, freeing the framework).
* its sequence must equal the CSV's ``designed_chain_sequence``.
  ``boltzgen_mapping_seq_identity`` reports that as a number; anything below 1.0 is an
  error, the sequence comes back ``NA``, **and the measured identity is still written
  on the error row** so the number that condemned it is visible.

``boltzgen_scaffold`` (which of the VHH frameworks was drawn, at random, for this
design) is read from ``json.loads(npz["source_context"].item())``, taking the entry
positionally matched to the binder chain and reading its ``source``. It is in no CSV
and no log line. Because position is the only link, the name is validated twice before
it is written -- against the scaffold list the run recorded for THIS GROUP in
``runs/<group>/.submit.json``, and
against the other chains' sources -- and ``boltzgen_scaffold_src`` says which checks
ran. Any failure is ``NA`` with a verdict that names it, never a scaffold name.

Safe to re-run: rows are rebuilt from the CSV and the structures on disk.

Usage:
    sapia collect boltzgen outputs/RUN --table <the table the run reserved>
"""

import hashlib
import json
import re
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, Iterable, NamedTuple

import gemmi
import numpy as np
import pandas as pd
import yaml

from prosapia.core import (
    RUN_META_FILENAME,
    CollectArgs,
    CollectCtx,
    CollectEach,
    Collected,
    DesignCtx,
)
from prosapia.utils import ensure_pdb

# Layout run_boltzgen.py imposes on the out_dir -- keep the two in step.
RUNS_DIRNAME = "runs"
DESIGN_DIRNAME = "intermediate_designs_inverse_folded"
REFOLD_DIRNAME = "refold_cif"
METRICS_CSV = "aggregate_metrics_analyze.csv"

# CSV key columns. Their absence is a hard error: without them nothing can be joined.
ID_COL = "id"
FILE_COL = "file_name"
# The designed chain's FULL sequence (framework + CDRs). Under nanobody-anything
# exactly one chain has designed residues, so this carries no ':' -- analyze.py builds
# it only from chains with at least one designed token (``designed_chain_ids``).
CHAIN_SEQ_COL = "designed_chain_sequence"
# Only the designed residues, concatenated.
DESIGN_SEQ_COL = "designed_sequence"
# Upstream's own count of designed tokens (analyze.py:766) -- a third, independent
# witness for the designed-residue count.
NUM_DESIGN_COL = "num_design"

# An atom belongs to a designed residue when its B-factor is the design mask's 1.
# The mask is written as a float 0.0/1.0, so compare with a midpoint rather than ==.
DESIGN_BFACTOR_CUTOFF = 0.5

# Metrics collected by default, spelled exactly as upstream spells them.
#
# This list is for the step set this tool forces (design, inverse_folding, folding,
# analysis). Three groups are DELIBERATELY ABSENT and must not be added back without
# changing the step set, because upstream does not compute them when inverse folding
# runs (cli/boltzgen.py:1308-1316 passes allatom_fold_metrics / delta_sasa_original /
# noncovalents_original = skip_inverse_folding, which is False here):
#   * the all-atom RMSDs `rmsd`, `rmsd_design`, `rmsd_target` -- use the `bb_` ones;
#   * `delta_sasa_original`, `design_sasa_*_original`;
#   * `plip_hbonds`, `plip_saltbridge`, `plip_hydrophobic` -- use the `_refolded` ones.
# `esmfold2_*` are absent for the same structural reason (that step never runs), and
# `designfolding-*` because design_folding is disabled for nanobody binders.
CORE_METRICS = (
    # --- confidence, from the Boltz-2 refold (eval_keys_confidence) -------------
    "design_ptm",
    "design_iptm",
    "design_iiptm",
    "design_to_target_iptm",
    "design_residue_iptm",
    "target_ptm",
    "ptm",
    "iptm",
    "protein_iptm",
    "min_design_to_target_pae",
    "min_interaction_pae",
    "design_ipsae_min",
    "design_to_target_ipsae",
    "target_to_design_ipsae",
    "complex_plddt",
    "complex_iplddt",
    "complex_pde",
    "complex_ipde",
    # --- backbone RMSDs: designed structure vs its own refold ------------------
    "bb_rmsd",
    "bb_rmsd_design",
    "bb_rmsd_target",
    "bb_rmsd_design_target",
    "bb_target_aligned_rmsd_design",
    # --- interface burial and contacts, on the refolded complex ----------------
    "delta_sasa_refolded",
    "design_sasa_bound_refolded",
    "design_sasa_unbound_refolded",
    "plip_hbonds_refolded",
    "plip_saltbridge_refolded",
    "plip_hydrophobic_refolded",
    # --- sequence-level developability -----------------------------------------
    "design_hydrophobicity",
    "design_chain_hydrophobicity",
    "liability_score",
    "liability_num_violations",
    "liability_details",
    "free_cysteines",
    # --- composition / provenance ----------------------------------------------
    "num_design",
    "num_tokens",
    "num_prot_tokens",
    "helix",
    "sheet",
    "loop",
    "target_id",
)

# Upstream writes threshold booleans whose names contain '<' and '>'
# ("rmsd<2.5", "design_ptm>80"). They would make column names a --filter cannot
# address, and they are all recomputable from the number beside them, so --all-metrics
# drops them rather than writing columns nobody can select on.
_SAFE_COLUMN = re.compile(r"^[A-Za-z_][A-Za-z0-9_.-]*$")

# The mmCIF writer's chain renaming: digits stripped, then a collision pool.
_STRIP_DIGITS = re.compile(r"\d+")


class MappingError(ValueError):
    """A design whose binder chain this collector refuses to guess at. Becomes an
    ``error`` status on that design's row -- never a best-effort mapping.

    ``data`` carries whatever was genuinely MEASURED before the refusal (e.g. the
    sequence identity that failed), so the error row shows the number that condemned it
    instead of an empty cell. It never carries a value that would read as a result.
    """

    def __init__(self, message: str, data: dict[str, Any] | None = None) -> None:
        super().__init__(message)
        self.data: dict[str, Any] = data or {}


class ChainView(NamedTuple):
    """One chain of a structure, as the collector needs to see it."""

    label: str
    seq: str  # one-letter, in file order
    designed: list[bool]  # per residue, from the design-mask B-factor


class BinderMapping(NamedTuple):
    """How the nanobody chain was identified, and what was read off it."""

    chain: str  # the label in the structure boltzgen_path points at
    sequence: str  # the designed chain's full sequence (framework + CDRs)
    design_positions: list[int]  # 1-based ordinals within that chain
    binder_len: int
    scaffold: Any  # basename of the chosen scaffold .cif, or pd.NA
    chain_map: str  # CHAIN_MAP_* -- how the chain <-> source_context link was tied
    scaffold_src: str  # SCAFFOLD_* -- how far the scaffold name was validated
    seq_identity: float


# `chain_map`: how the binder chain was tied to its source_context entry. Only the two
# `ok*` values mean the link holds; it rests on chain POSITION, which is how
# source_context is built (data/source_context.py:198-206, zipped strict=True).
# `update_designed` writes `output_chain` UNCONDITIONALLY for every entry
# (``data/source_context.py:200-205``) and ``DesignWriter`` guards the call on the
# context existing at all, not on any chain (``writer.py:506-516``). So a per-chain
# absence CANNOT be produced: either the whole context is JSON null, or every entry has
# the field. CHAIN_MAP_MISSING_OUTPUT is therefore defensive -- if it ever fires,
# something is happening that nobody here understands, and it must not read as success.
CHAIN_MAP_CONFIRMED = "ok"  # output_chain, digits stripped, == the written label
CHAIN_MAP_RENAMED = "ok_renamed"  # output_chain present and different: writer re-pooled
CHAIN_MAP_MISSING_OUTPUT = "output_chain_missing"  # cannot happen -> raises
CHAIN_MAP_NONE = "no_source_context"  # whole context null: provenance LOST, not absent

# The chain_map verdicts under which the binder chain / sequence may be believed.
CHAIN_MAP_OK = (CHAIN_MAP_CONFIRMED, CHAIN_MAP_RENAMED)

# `scaffold_src`: how far boltzgen_scaffold itself was validated. Only the first two
# carry a scaffold name; the rest are NA.
SCAFFOLD_VERIFIED = "in_run_scaffold_list"  # matched a scaffold THIS RUN offered
SCAFFOLD_DISTINCT = "distinct_from_target"  # no allow-list (--spec), but not a target
SCAFFOLD_NOT_IN_RUN = "not_in_run_scaffold_list"  # -> NA
SCAFFOLD_IS_TARGET = "is_a_target_chain_source"  # -> NA
SCAFFOLD_NONE = "no_source_context"  # -> NA


class BoltzGenCollectArgs(CollectArgs):
    metrics: str
    all_metrics: bool


def add_collect_boltzgen_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--metrics",
        type=str,
        default="",
        help="Comma-separated extra columns from aggregate_metrics_analyze.csv to "
        "collect on top of the core set (e.g. 'seq_recovery,novelty_refolded'). Names "
        "are BoltzGen's own, as spelled in that CSV.",
    )
    parser.add_argument(
        "--all-metrics",
        action="store_true",
        help="Collect every column of aggregate_metrics_analyze.csv instead of the "
        "core set. Wide (~100 columns). Columns whose names contain '<' or '>' "
        "(upstream's threshold booleans, e.g. 'design_iptm>80') are still skipped: a "
        "--filter cannot address them and the number they threshold is collected.",
    )


# ------------------------------------------------------------------ structures


def one_letter(resname: str) -> str:
    """Three-letter residue name -> one-letter code, ``X`` for anything else."""
    info = gemmi.find_tabulated_residue(resname)
    code = info.one_letter_code if info is not None else ""
    return code.upper() if code and code.isalpha() else "X"


def read_chains(path: Path, *, design_mask: bool) -> list[ChainView]:
    """The first model's polymer chains, in file order.

    **Every non-water chain is returned, polymer or not**, because the index into this
    list is used against ``source_context["chains"]`` -- which also holds an entry per
    non-polymer chain (``data/source_context.py:14-30`` has no type filter;
    ``mmcif.py:1083-1115`` builds NONPOLYMER chains; waters are removed upstream at
    ``mmcif.py:1013``). ``select``/``update_designed`` zip the context against the
    structure's chains with ``strict=True``, so the correspondence is positional and
    **enforced**, not conventional -- but only if both sides count the same chains.
    Counting polymers only, as an earlier version did, offset the whole lookup on any
    target carrying a ligand, ion or glycan chain. A non-polymer chain has ``seq == ""``
    and no designed residues, so it can never be mistaken for the binder.

    ``design_mask`` says whether this file is one whose B-factor column IS the
    atom-level design mask -- true only for the design mmCIF written by
    ``model/validation/design.py`` / ``task/predict/writer.py:453``. When it is false,
    ``ChainView.designed`` is ``None`` rather than a list.

    That distinction is load-bearing and deliberately not a default. The REFOLDED
    complex carries per-atom **pLDDT** in the same column
    (``task/predict/writer.py:168-176``), so reading it with the same ``>= 0.5`` rule
    would report every residue of every chain as designed -- a perfectly plausible,
    entirely wrong answer. ``None`` makes that a ``TypeError`` instead.
    """
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_waters()
    out: list[ChainView] = []
    for chain in st[0]:
        polymer = chain.get_polymer()
        residues = list(polymer) if len(polymer) else []
        seq = "".join(one_letter(res.name) for res in residues)
        designed = (
            [
                any(atom.b_iso >= DESIGN_BFACTOR_CUTOFF for atom in res)
                for res in residues
            ]
            if design_mask
            else None
        )
        out.append(ChainView(chain.name, seq, designed))
    return out


def identity(a: str, b: str) -> float:
    """Fraction of identical residues between two sequences, 0.0 .. 1.0.

    Equal-length sequences are compared position by position. Unequal ones go through
    a Needleman-Wunsch global alignment (match 1, mismatch 0, gap -1) implemented here
    rather than pulling in biopython for one function; identity is then matches over
    the length of the longer sequence, so an alignment that only matches a fragment
    cannot read as 1.0.
    """
    if not a or not b:
        return 0.0
    if len(a) == len(b):
        return sum(x == y for x, y in zip(a, b)) / len(a)

    n, m = len(a), len(b)
    prev = [-j for j in range(m + 1)]
    prev_hits = [0] * (m + 1)
    for i in range(1, n + 1):
        cur = [-i] + [0] * m
        cur_hits = [0] * (m + 1)
        for j in range(1, m + 1):
            hit = a[i - 1] == b[j - 1]
            diag = prev[j - 1] + (1 if hit else 0)
            up, left = prev[j] - 1, cur[j - 1] - 1
            best = max(diag, up, left)
            cur[j] = best
            if best == diag:
                cur_hits[j] = prev_hits[j - 1] + (1 if hit else 0)
            elif best == up:
                cur_hits[j] = prev_hits[j]
            else:
                cur_hits[j] = cur_hits[j - 1]
        prev, prev_hits = cur, cur_hits
    return prev_hits[m] / max(n, m)


def render_positions(nums: list[int]) -> str:
    """``[26,27,28,51,52]`` -> ``'26:28,51:52'``.

    The prosapia position mini-language (``,`` fragments, ``start:end`` inclusive,
    1-indexed within the chain) -- the same string ``proteinmpnn --fixed-positions``
    parses. A single position stays a bare number.
    """
    if not nums:
        return ""
    runs: list[list[int]] = []
    for n in sorted(set(nums)):
        if runs and n == runs[-1][1] + 1:
            runs[-1][1] = n
        else:
            runs.append([n, n])
    return ",".join(str(lo) if lo == hi else f"{lo}:{hi}" for lo, hi in runs)


# -------------------------------------------------------------------- mapping


def binder_chain_index(chains: list[ChainView]) -> int:
    """The index of the single chain carrying designed residues.

    Zero designed chains means the design mask never reached the written B-factors;
    more than one means the spec designed something other than a lone nanobody. Both
    are errors: the sequence columns and the scaffold would otherwise be attached to a
    chain nobody verified.
    """
    designed = [i for i, c in enumerate(chains) if any(c.designed)]
    if len(designed) == 1:
        return designed[0]
    raise MappingError(
        f"expected exactly one designed chain in the design mmCIF, found "
        f"{len(designed)} ({[chains[i].label for i in designed]} of "
        f"{[c.label for c in chains]}); the design-mask B-factor cannot identify the "
        f"nanobody chain"
    )


def npz_design_count(npz_path: Path) -> int:
    """Number of designed tokens in the per-design ``.npz``.

    ``design_mask`` is per token over the WHOLE complex (``writer.py``'s metadata
    dict), so only its True count is used -- never its positions, which would assume a
    token-to-residue correspondence this collector does not need.
    """
    with np.load(npz_path, allow_pickle=False) as data:
        if "design_mask" not in data:
            raise MappingError(f"{npz_path.name} has no 'design_mask'")
        return int(np.asarray(data["design_mask"]).astype(bool).sum())


def npz_source_context(npz_path: Path) -> dict | None:
    """The ``source_context`` dict, or ``None`` when it is absent or unreadable."""
    try:
        with np.load(npz_path, allow_pickle=False) as data:
            if "source_context" not in data:
                return None
            context = json.loads(str(data["source_context"].item()))
    except Exception as e:  # noqa: BLE001
        # Broad on purpose: a damaged .npz raises whatever numpy/zipfile feel like
        # (upstream issue #223 is an intermittent zipfile.BadZipFile, which subclasses
        # Exception, not ValueError). An unreadable context costs the SCAFFOLD, not the
        # design -- the caller reports `no_source_context` and writes pd.NA. Printed,
        # because otherwise a typo'd key and a file BoltzGen never wrote look identical
        # from the table: a whole scaffold column of NA with nothing to explain it.
        print(f"  {npz_path.name}: source_context unreadable ({type(e).__name__}: {e})")
        return None
    return context if isinstance(context, dict) else None


def scaffold_of(
    context: dict | None,
    index: int,
    label: str,
    allowed_stems: frozenset[str] | None,
) -> tuple[Any, str, str]:
    """``(scaffold basename or NA, chain_map verdict, scaffold_src verdict)``.

    ``source_context["chains"]`` is built by zipping the structure's chains in order
    (``data/source_context.py:198-206``, ``strict=True``), so **chain position is the
    link** and ``output_chain`` is only a confirmation of it.

    Position alone is not enough to publish a scaffold name. If the ordering ever
    differed -- by input entity rather than by output chain, say -- this would return
    the TARGET file's basename for every design, under ``status == OK``, and the
    diversity analysis the skill recommends would report one framework at 100%. A
    partial permutation would be worse: every design attributed to the wrong one of
    five frameworks, and the framework decision made on noise. So the name is checked
    twice before it is written:

    * **against the run's own scaffold list.** ``run_boltzgen.py`` records the resolved
      scaffold YAML paths in that group's ``.submit.json``; a ``source`` whose stem is
      not among them
      is not a scaffold this run offered. (Every shipped scaffold's ``.yaml`` and its
      ``path:`` ``.cif`` share a stem -- checked for all seven -- so stems compare.)
      ``allowed_stems`` is ``None`` for a ``--spec`` run, where the scaffolds come from
      the user's YAML and no list exists; the verdict then says so.
    * **against the other chains' sources.** The binder's ``source`` must differ from
      every non-binder chain's. This catches the "positional lookup returned the
      target" failure with no allow-list at all, so it runs in both cases.

    Either check failing yields ``pd.NA`` and a verdict that names the failure --
    never a scaffold name, and never a verdict that reads as success.
    """
    chains = (context or {}).get("chains")
    if not isinstance(chains, list) or index >= len(chains):
        return pd.NA, CHAIN_MAP_NONE, SCAFFOLD_NONE
    entry = chains[index]
    if not isinstance(entry, dict) or not entry.get("source"):
        return pd.NA, CHAIN_MAP_NONE, SCAFFOLD_NONE

    source = str(entry["source"])
    raw_output_chain = entry.get("output_chain")
    if raw_output_chain is None:
        # Upstream cannot produce this (see CHAIN_MAP_MISSING_OUTPUT). Refuse rather
        # than invent a verdict for a state we cannot explain.
        raise MappingError(
            f"source_context entry {index} has no 'output_chain' "
            f"({CHAIN_MAP_MISSING_OUTPUT}). "
            f"update_designed writes it unconditionally for every chain "
            f"(source_context.py:200-205), so this npz was not written by the pipeline "
            f"this tool targets and nothing positional about it can be trusted"
        )
    if _STRIP_DIGITS.sub("", str(raw_output_chain)) == label:
        chain_map = CHAIN_MAP_CONFIRMED
    else:
        chain_map = CHAIN_MAP_RENAMED

    others = {
        str(c["source"])
        for i, c in enumerate(chains)
        if i != index and isinstance(c, dict) and c.get("source")
    }
    if source in others:
        # The entry at the binder's position names a file that a non-binder chain also
        # came from -- i.e. the target. The positional link is broken.
        return pd.NA, chain_map, SCAFFOLD_IS_TARGET

    stem = Path(source).stem
    if allowed_stems is None:
        return Path(source).name, chain_map, SCAFFOLD_DISTINCT
    if stem not in allowed_stems:
        return pd.NA, chain_map, SCAFFOLD_NOT_IN_RUN
    return Path(source).name, chain_map, SCAFFOLD_VERIFIED


def map_binder(
    design_cif: Path,
    npz_path: Path,
    refold: Path,
    csv_seq: str,
    csv_num_design: Any,
    allowed_stems: frozenset[str] | None,
) -> BinderMapping:
    """Identify the nanobody chain and everything read off it, or raise MappingError.

    Order matters: the chain is found from the design mask, all three independent
    witnesses of the designed-residue count are checked against it, the design chain is
    matched **residue for residue** to the refolded one, and only then is a sequence
    attached to a chain label.

    A raised ``MappingError`` carries whatever WAS measured in ``.data``, so an error
    row still shows the number that condemned it rather than an empty cell.
    """
    if ":" in csv_seq:
        raise MappingError(
            f"'{CHAIN_SEQ_COL}' holds {csv_seq.count(':') + 1} designed chains; this "
            f"tool's premise is a single designed nanobody chain, so there is no "
            f"unambiguous binder sequence"
        )

    design_chains = read_chains(design_cif, design_mask=True)
    if not design_chains:
        raise MappingError(f"no polymer chains in {design_cif.name}")
    index = binder_chain_index(design_chains)
    binder = design_chains[index]
    assert binder.designed is not None  # read with design_mask=True
    positions = [i for i, flag in enumerate(binder.designed, start=1) if flag]

    n_npz = npz_design_count(npz_path)
    if n_npz != len(positions):
        raise MappingError(
            f"designed-residue count disagrees: {len(positions)} from the design-mask "
            f"B-factor of chain {binder.label}, {n_npz} from {npz_path.name}'s "
            f"design_mask"
        )
    # The THIRD witness is not optional. The skill promises three, and `num_design` is
    # written unconditionally (analyze.py:766) -- its absence means a different
    # upstream version, which is exactly when the other two stop being trustworthy too.
    if csv_num_design is None or pd.isna(csv_num_design):
        raise MappingError(
            f"the CSV has no usable '{NUM_DESIGN_COL}', so the designed-residue count "
            f"has only two witnesses instead of three; this tool's guarantee does not "
            f"hold on this output"
        )
    n_csv = int(csv_num_design)
    if n_csv != len(positions):
        raise MappingError(
            f"designed-residue count disagrees: {len(positions)} from the "
            f"design-mask B-factor, {n_csv} from the CSV's '{NUM_DESIGN_COL}'"
        )

    refold_chains = read_chains(refold, design_mask=False)
    if index >= len(refold_chains):
        raise MappingError(
            f"the refolded complex has {len(refold_chains)} polymer chain(s) but the "
            f"binder is chain index {index} of the design mmCIF"
        )
    refold_binder = refold_chains[index]
    # Residue for residue, not merely the same length. The designed POSITIONS are read
    # off the design mmCIF and reported against the refolded complex, so a
    # length-preserving shift between the two would pass a length check, and
    # boltzgen_fixed_resnums would then fix the CDRs and free the framework -- the
    # exact inverse of the intended experiment, with nothing to show for it.
    if refold_binder.seq != binder.seq:
        raise MappingError(
            f"the binder chain differs between the design mmCIF (chain "
            f"{binder.label}, {len(binder.seq)} aa) and the refolded complex (chain "
            f"{refold_binder.label}, {len(refold_binder.seq)} aa, identity "
            f"{identity(binder.seq, refold_binder.seq):.3f}); the designed positions "
            f"cannot be carried across"
        )

    seq_identity = identity(csv_seq, refold_binder.seq)
    if seq_identity < 1.0:
        raise MappingError(
            f"'{CHAIN_SEQ_COL}' is not the sequence of chain {refold_binder.label} of "
            f"the refolded complex (identity {seq_identity:.3f}); the asym-id -> "
            f"chain-label mapping cannot be trusted for this design",
            data={"mapping_seq_identity": seq_identity},
        )

    scaffold, chain_map, scaffold_src = scaffold_of(
        npz_source_context(npz_path), index, binder.label, allowed_stems
    )
    return BinderMapping(
        chain=refold_binder.label,
        sequence=csv_seq,
        design_positions=positions,
        binder_len=len(refold_binder.seq),
        scaffold=scaffold,
        chain_map=chain_map,
        scaffold_src=scaffold_src,
        seq_identity=seq_identity,
    )


# -------------------------------------------------------------------- collect


SUBMIT_RECORD = ".submit.json"
SPECS_DIRNAME = "specs"


def read_submit_record(group_dir: Path) -> tuple[dict | None, str]:
    """``(record, state)`` for ``runs/<group>/.submit.json``.

    ``state`` is ``"present"``, ``"absent"`` or a description of the damage, so the
    caller can say WHICH it was instead of reporting every case as "missing".

    The record lives beside the designs it describes, not in the run-level
    ``.meta.json``, because the latter is keyed by group NAME, merged across submits,
    rewritten by the framework before the tool's own builder runs, and written without
    a lock. See ``run_boltzgen.read_submit_record`` for the three silent failures that
    caused.
    """
    path = group_dir / SUBMIT_RECORD
    if not path.is_file():
        return None, "absent"
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as e:
        return None, f"unreadable ({type(e).__name__})"
    if not isinstance(data, dict):
        return None, f"not a JSON object (got {type(data).__name__})"
    return data, "present"


def claims_from_spec(spec_path: Path) -> set[str] | None:
    """Which steering claims the SPEC asks for, or ``None`` if it cannot be read.

    The spec YAML is what BoltzGen actually consumed, so it is the independent witness
    that an epitope was requested at all -- the submit record carries the residue NAMES,
    which the spec does not.
    """
    try:
        data = yaml.safe_load(spec_path.read_text())
    except (OSError, yaml.YAMLError):
        return None
    if not isinstance(data, dict):
        return None
    claims: set[str] = set()
    for entity in data.get("entities") or []:
        block = (entity or {}).get("file") if isinstance(entity, dict) else None
        for item in (block or {}).get("binding_types") or []:
            chain = (item or {}).get("chain") if isinstance(item, dict) else None
            for key in ("binding", "not_binding"):
                if chain and key in chain:
                    claims.add(key)
    return claims


class MetaUnreadable(RuntimeError):
    """The run-level sidecar is absent or unparseable.

    Survivable. It carries only run-level facts -- ``runs_root`` (which defaults to
    ``<out_dir>/runs`` anyway), ``protocol``, ``steps`` and ``root_designs``. Nothing a
    design's trust rests on is in it: the epitope attestation, the residue and chain
    frames, the chain-id translation and the scaffold allow-list all live per group in
    ``runs/<group>/.submit.json``, precisely because the framework rewrites this file
    before the tool's own builder runs and silently replaces a corrupt one
    (``base_run.py:345`` vs ``:377``).

    The collector therefore warns and continues; a missing PER-GROUP record is the hard
    stop.
    """


def _run_meta(out_dir: Path) -> dict:
    """The run-level sidecar. Raises ``MetaUnreadable`` when absent or damaged."""
    path = out_dir / RUN_META_FILENAME
    if not path.is_file():
        raise MetaUnreadable(
            f"{path} does not exist, so runs_root falls back to <out_dir>/runs. The "
            f"usual cause is a -l/--dir-label or -t/--table resolving to a different "
            f"out_dir than the run used. Per-group attestation is unaffected -- it "
            f"lives in runs/<group>/{SUBMIT_RECORD}."
        )
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as e:
        raise MetaUnreadable(
            f"{path} cannot be read ({type(e).__name__}: {e}); runs_root falls back to "
            f"<out_dir>/runs. Per-group attestation is unaffected -- it lives in "
            f"runs/<group>/{SUBMIT_RECORD}."
        ) from e
    if not isinstance(data, dict):
        raise MetaUnreadable(f"{path} is not a JSON object (got {type(data).__name__})")
    return data


def _cell(row: pd.Series, col: str) -> Any:
    """A CSV cell, with an absent column or a NaN read as ``None``."""
    if col not in row.index:
        return None
    value = row[col]
    try:
        return None if value is None or pd.isna(value) else value
    except (TypeError, ValueError):  # array-like: not an NA scalar
        return value


def _metric_data(row: pd.Series, columns: list[str]) -> dict[str, Any]:
    """The declared metric columns, with an absent one left out entirely.

    A column this tool declares but the CSV does not carry is NOT written as 0 --
    prosapia leaves it NA, which is what "not measured" means here.
    """
    return {col: row[col] for col in columns if col in row.index}


def collect_boltzgen(ctx: CollectCtx[BoltzGenCollectArgs]) -> CollectEach:
    """Per-target BoltzGen collector. A create tool: the framework iterates the ready
    parents (table rows, or the root design groups the run recorded) and this mints one
    child row per design BoltzGen produced, carrying ``parent`` for lineage. A target
    with no designs yields no rows. The framework stamps status/path/parent_name."""
    # The run-level sidecar carries only run-level facts (runs_root, protocol, steps,
    # root_designs). Everything a design's trust rests on -- the epitope attestation,
    # the frames, the chain-id translation and the scaffold allow-list -- lives per
    # group in runs/<group>/.submit.json, for the reasons in read_submit_record. So an
    # unreadable .meta.json is survivable here; a missing per-group record is not.
    #
    # There is deliberately NO fallback to an older top-level layout. This tool has
    # never run on real data, so there are no out_dirs to support, and a fallback would
    # reintroduce exactly the failure this design closed: a group with no record of its
    # own, in an out_dir that happens to carry stale top-level keys, silently inheriting
    # an attestation from a different run.
    try:
        meta = _run_meta(ctx.out_dir)
    except MetaUnreadable as e:
        print(f"  {e}")
        meta = {}
    recorded_root = meta.get("runs_root")
    runs_root = Path(recorded_root) if recorded_root else ctx.out_dir / RUNS_DIRNAME
    specs_dir = ctx.out_dir / SPECS_DIRNAME

    run_dir = ctx.args.run_dir
    extra_metrics = [m.strip() for m in ctx.args.metrics.split(",") if m.strip()]
    print(f"Collecting BoltzGen designs from {runs_root}")
    warned_absent: set[str] = set()

    def group_record(name: str) -> tuple[dict, str | None]:
        """``(record, error)`` for one group.

        ``error`` is a message when the record is missing, damaged, or disagrees with
        the spec BoltzGen consumed; the caller turns it into an ``error:`` status on
        every design of that group and writes NONE of the attestation columns, rather
        than publishing a claim that does not describe these designs.
        """
        group_dir = runs_root / name
        spec_path = specs_dir / f"{name}.yaml"
        record, state = read_submit_record(group_dir)

        if record is None:
            return {}, (
                f"{SUBMIT_RECORD} beside {group_dir} is {state}. The usual cause is a "
                f"-l/--dir-label or -t/--table resolving to a different out_dir than "
                f"the run used -- check those first; nothing is necessarily wrong with "
                f"the outputs. Otherwise this group was not submitted from here and "
                f"nothing attests its epitope."
            )

        # Cross-check against the spec BoltzGen actually consumed. specs/<group>.yaml is
        # rewritten on every submit, so a disagreement means the record and the spec
        # came from different runs.
        if not spec_path.is_file():
            print(
                f"{name}: {spec_path} is missing, so the spec cross-check is SKIPPED "
                f"for this group -- its attestation rests on {SUBMIT_RECORD} alone."
            )
        else:
            sha = hashlib.sha256(spec_path.read_bytes()).hexdigest()
            if record.get("spec_sha256") and record["spec_sha256"] != sha:
                return record, (
                    f"{SUBMIT_RECORD} was written for a different spec than "
                    f"{spec_path.name} currently holds (recorded "
                    f"{record['spec_sha256'][:12]}, on disk {sha[:12]}); the epitope "
                    f"attestation does not describe these designs"
                )
            spec_claims = claims_from_spec(spec_path)
            recorded_claims = set(record.get("claims") or [])
            if spec_claims is None:
                print(
                    f"{name}: {spec_path.name} could not be parsed, so the claims "
                    f"cross-check is SKIPPED for this group."
                )
            elif spec_claims != recorded_claims:
                return record, (
                    f"the spec requests {sorted(spec_claims) or 'no'} steering but "
                    f"{SUBMIT_RECORD} records {sorted(recorded_claims) or 'none'}; the "
                    f"attestation and the spec disagree"
                )

        # Per CLAIM: a --not-binding-only run has no binding resnames and that is
        # correct, not a missing attestation.
        for claim, key in (
            ("binding", "binding_resnames"),
            ("not_binding", "not_binding_resnames"),
        ):
            if claim in (record.get("claims") or []) and not record.get(key):
                return record, (
                    f"the spec requests '{claim}' steering but {SUBMIT_RECORD} records "
                    f"no {key}; its residue frame and chain-id translation cannot be "
                    f"attested. A hand-written --spec YAML carrying binding_types "
                    f"always lands here: this tool only resolves residue NAMES for the "
                    f"epitope it built itself."
                )
        return record, None

    def one(d: DesignCtx) -> Iterable[Collected]:
        design_dir = runs_root / d.name / DESIGN_DIRNAME
        csv_path = design_dir / METRICS_CSV
        if not csv_path.is_file():
            print(f"{d.name}: no {METRICS_CSV} under {design_dir}, skipping")
            return
        try:
            df = pd.read_csv(csv_path)
        except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            print(f"{d.name}: unreadable {csv_path}: {e}")
            return
        missing_keys = [c for c in (ID_COL, FILE_COL, CHAIN_SEQ_COL) if c not in df.columns]
        if missing_keys:
            print(
                f"{d.name}: {METRICS_CSV} has no {missing_keys} column(s); nothing can "
                f"be joined by name. Not collecting this target."
            )
            return

        if ctx.args.all_metrics:
            columns = [
                c
                for c in df.columns
                if c not in (ID_COL, FILE_COL, CHAIN_SEQ_COL) and _SAFE_COLUMN.match(c)
            ]
        else:
            columns = [c for c in (*CORE_METRICS, *extra_metrics) if c in df.columns]
            absent = [
                c
                for c in (*CORE_METRICS, *extra_metrics)
                if c not in df.columns and c not in warned_absent
            ]
            if absent:
                warned_absent.update(absent)
                print(f"  not in this CSV (left NA): {', '.join(sorted(absent))}")

        record, record_error = group_record(d.name)
        # The absence of the trust record must not be readable as the absence of the
        # thing it attests. boltzgen_binding_resnames is NA when no steering was asked
        # for -- so a group whose record is missing, stale or at odds with its own spec
        # is an ERROR, not a design with an unrecorded epitope.
        #
        # And an error row carries NONE of it. The whole content of `record_error` is
        # "this attestation does not describe these designs"; shipping that same
        # attestation in the row's columns would contradict the message beside it, and
        # the skill promises an error row is all-NA.
        if record_error:
            print(f"{d.name}: {record_error}")
            epitope = anti_epitope = ""
        else:
            epitope = record.get("binding_resnames") or ""
            anti_epitope = record.get("not_binding_resnames") or ""
        # The scaffold allow-list is the one this GROUP was submitted with, not the
        # latest run's: a later submit with a different --scaffolds must not re-point an
        # earlier group's provenance check at a set it never used.
        group_scaffolds = list(record.get("scaffolds") or [])
        allowed_stems = (
            frozenset(Path(s).stem for s in group_scaffolds)
            if group_scaffolds
            else None
        )
        print(
            f"{d.name}: scaffold allow-list "
            + (
                f"({len(allowed_stems)}) {', '.join(sorted(allowed_stems))}"
                if allowed_stems
                else f"NONE recorded -- scaffold_src can only reach "
                f"'{SCAFFOLD_DISTINCT}', not '{SCAFFOLD_VERIFIED}'"
            )
        )

        # A duplicated `id` would collapse two designs into one row, last-wins, with
        # nothing to show it happened: the by-name join (file_name == id + '.cif') is
        # ambiguous for every artifact of that design.
        duplicated = sorted(set(df[ID_COL][df[ID_COL].duplicated()].astype(str)))
        if duplicated:
            print(
                f"{d.name}: {METRICS_CSV} repeats id(s) {duplicated}; the by-name join "
                f"is ambiguous for them and they are collected as errors."
            )

        n = 0
        emitted_ids: set[str] = set()
        scaffold_verdicts: dict[str, int] = {}
        for _, row in df.sort_values(ID_COL).iterrows():
            design_id = str(row[ID_COL])
            file_name = str(row[FILE_COL])
            row_name = (
                design_id if design_id.startswith(d.name) else f"{d.name}_{design_id}"
            )
            data: dict[str, Any] = _metric_data(row, columns)
            data["id"] = design_id
            if epitope:
                data["binding_resnames"] = epitope
            # A SEPARATE column: not_binding is the opposite claim in identical
            # `A91=LYS` notation, and sharing a column would let a reader checking
            # where the design was steered TO read anti-epitope residues as epitope.
            if anti_epitope:
                data["not_binding_resnames"] = anti_epitope
            if epitope or anti_epitope:
                # Unreachable on an error row: `epitope`/`anti_epitope` are forced
                # empty when record_error is set, so none of the attestation columns
                # is written there.
                # The identifier witnesses belong to BOTH claims: an anti-epitope is
                # written in the same `A91=LYS` notation and rests on exactly the same
                # residue frame and chain-id translation, so a --not-binding-only run
                # must not attest residues with no frame beside them.
                data["binding_frame"] = record.get("res_frame") or pd.NA
                data["binding_chain_frame"] = record.get("chain_frame") or pd.NA
                # Both sides of the author -> BoltzGen chain-id translation, so a wrong
                # one is visible in the table rather than inferred.
                data["chain_id_map"] = record.get("chain_id_map") or pd.NA

            if design_id in duplicated:
                if design_id in emitted_ids:
                    continue  # one error row per ambiguous id, not one per occurrence
                emitted_ids.add(design_id)
                yield Collected(
                    name=row_name,
                    parent=d.name,
                    path="",
                    status=(
                        f"error: '{design_id}' appears more than once in "
                        f"{METRICS_CSV}; the by-name join to its CIF/npz is ambiguous"
                    ),
                    data=data,
                )
                n += 1
                continue

            if record_error:
                yield Collected(
                    name=row_name,
                    parent=d.name,
                    path="",
                    status=f"error: {record_error}",
                    data=data,
                )
                n += 1
                continue

            design_cif = design_dir / file_name
            npz_path = design_dir / f"{design_id}.npz"
            refold_cif = design_dir / REFOLD_DIRNAME / f"{design_id}.cif"

            missing = [
                p.name for p in (design_cif, npz_path, refold_cif) if not p.is_file()
            ]
            if missing:
                print(f"{d.name}: {design_id} missing {missing}")
                # NOT the path that was just found absent: filter_ready keys on a
                # column being non-empty, not on <leaf>_status, so a downstream
                # `-i boltzgen_cif_path` would pick this row up and die on it.
                data["cif_path"] = ""
                yield Collected(
                    name=row_name,
                    parent=d.name,
                    path="",
                    status=f"error: missing output file(s) {missing}",
                    data=data,
                )
                n += 1
                continue

            data["cif_path"] = str(refold_cif)
            try:
                refold_pdb = ensure_pdb(refold_cif, run_dir)
            except Exception as e:  # noqa: BLE001 -- see the comment on map_binder below
                print(f"{d.name}: {design_id} unreadable refolded mmCIF: {e}")
                yield Collected(
                    name=row_name,
                    parent=d.name,
                    path="",
                    status=f"error: unreadable refolded mmCIF ({type(e).__name__})",
                    data=data,
                )
                n += 1
                continue

            try:
                mapping = map_binder(
                    design_cif,
                    npz_path,
                    refold_pdb,
                    str(row[CHAIN_SEQ_COL]),
                    _cell(row, NUM_DESIGN_COL),
                    allowed_stems,
                )
            except Exception as e:  # noqa: BLE001
                # Deliberately broad. A corrupt input file is an ERROR ROW, never an
                # aborted collect: the exceptions that reach here are whatever gemmi,
                # numpy and zipfile raise on a damaged file, and they do not share a
                # base class. gemmi raises RuntimeError on a parse failure, and
                # upstream issue #223 is specifically an intermittent
                # `zipfile.BadZipFile` on BoltzGen `.npz` in long runs -- which
                # subclasses Exception, not ValueError. Letting one of those out of
                # this loop propagates past the driver's save_frame, so ONE bad file in
                # two hundred throws away the whole collect of a 12-hour run, and
                # re-running collects nothing until someone deletes the file by hand.
                # KeyboardInterrupt/SystemExit are BaseException and still propagate.
                #
                # The metrics are still true of SOME chain pair, but nothing that names
                # a chain or a residue can be trusted, so none of it is written and the
                # path is left empty rather than inviting a downstream step.
                print(f"{d.name}: {design_id} {type(e).__name__}: {e}")
                data.update(getattr(e, "data", {}) or {})
                yield Collected(
                    name=row_name,
                    parent=d.name,
                    path="",
                    status=f"error: {type(e).__name__}: {e}",
                    data=data,
                )
                n += 1
                continue

            fixed = [
                i
                for i in range(1, mapping.binder_len + 1)
                if i not in set(mapping.design_positions)
            ]
            data.update(
                {
                    # Always one of CHAIN_MAP_OK here: every other verdict either
                    # raised above or is CHAIN_MAP_NONE, which is carried through so a
                    # caller can exclude it -- the documented trust filter does.
                    "chain_map": mapping.chain_map,
                    "scaffold_src": mapping.scaffold_src,
                    "mapping_seq_identity": mapping.seq_identity,
                    "binder_chain": mapping.chain,
                    "sequence": mapping.sequence,
                    "cdr_sequence": _cell(row, DESIGN_SEQ_COL),
                    "scaffold": mapping.scaffold,
                    "design_resnums": render_positions(mapping.design_positions),
                    "fixed_resnums": render_positions(fixed),
                    "n_design_res": len(mapping.design_positions),
                    "binder_len": mapping.binder_len,
                }
            )
            scaffold_verdicts[mapping.scaffold_src] = (
                scaffold_verdicts.get(mapping.scaffold_src, 0) + 1
            )
            yield Collected(
                name=row_name, parent=d.name, path=refold_pdb, data=data
            )
            n += 1
        print(f"{d.name}: OK ({n} design(s))")
        named = sum(
            c for v, c in scaffold_verdicts.items()
            if v in (SCAFFOLD_VERIFIED, SCAFFOLD_DISTINCT)
        )
        if scaffold_verdicts and not named:
            # Every scaffold NA is almost never a BoltzGen problem. The binder's index
            # is its position among POLYMER chains, and it indexes source_context's
            # chain list; a target carrying a ligand, glycan or ion chain can make the
            # two disagree and send every lookup to the wrong entry. Fail-safe (NA, not
            # a wrong name), but say where to look.
            print(
                f"{d.name}: NO design got a validated scaffold "
                f"({dict(sorted(scaffold_verdicts.items()))}). If the target carries "
                f"non-polymer chains (ligand, glycan, ion), the binder's polymer index "
                f"may not line up with source_context's chain list. The sequences and "
                f"chain labels are unaffected; only boltzgen_scaffold is."
            )

    return one
