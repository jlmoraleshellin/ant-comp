#!/usr/bin/env python3
"""
Collect BoltzGen binder designs into the (child) design table.

``run_boltzgen.py`` points each target's ``--output`` at ``<out_dir>/runs/<name>/``,
where BoltzGen writes (``README.md`` "Pipeline output", verified against the pinned
SHA):

    runs/<name>/intermediate_designs_inverse_folded/
        <id>.cif                     the inverse-folded design complex
        <id>.npz                     design_mask + source_context for that design
        refold_cif/<id>.cif          the REFOLDED complex (target + binder)
        refold_design_cif/<id>.cif   the BINDER ALONE -- protein-anything ONLY
        aggregate_metrics_analyze.csv

``refold_cif/<id>.cif`` is what every analysis metric was computed on, so it is what
``boltzgen_path`` points at (converted to PDB; the mmCIF stays as
``boltzgen_cif_path``).

THE PROTOCOL IS A COLUMN, and it decides what several others mean
------------------------------------------------------------------
``boltzgen_protocol`` is read back from the group's own ``runs/<group>/.submit.json``
-- written at submit, beside the designs it describes -- and **not** from a flag at
collect time or from the run-level ``.meta.json``, which is merged across submits. A
run_dir can hold nanobody and minibinder groups side by side; a run-level value would
label both with whichever ran last. Read it before anything else. Three columns change
meaning with it (``scaffold``/``scaffold_src``, ``cdr_sequence``, ``liability_*``) and
two exist only under ``protein-anything`` (``binder_refold_cif``,
``designfolding-*``).

**The binder-alone structure exists only under ``protein-anything``.** Upstream runs
``design_folding`` for ``protein-anything`` and ``protein-small_molecule`` only
(``cli/boltzgen.py:1237``), and that step is what writes ``refold_design_cif/``
(``task/predict/writer.py:94-104``). Under ``nanobody-anything`` it is never written
(upstream issue #148) and a self-consistency check still needs ``chainsel`` to cut the
nanobody chain out first; under ``protein-anything`` the answer is already in the
table as ``boltzgen_designfolding-bb_rmsd_design``.

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
# The BINDER FOLDED ALONE, written only by the `design_folding` step -- which only
# protein-anything runs (cli/boltzgen.py:1237). const.refold_design_cif_dirname;
# task/predict/writer.py:94-104 writes it beside fold_out_design_npz/<id>.npz and
# publishes the npz LAST, so an npz implies this cif. Absent under nanobody-anything
# (upstream issue #148) and that is not a failure there.
REFOLD_DESIGN_DIRNAME = "refold_design_cif"
METRICS_CSV = "aggregate_metrics_analyze.csv"

# The two protocols this tool's run module can submit, as recorded per group in
# .submit.json. Read back into boltzgen_protocol -- a TRUST column.
NANOBODY_PROTOCOL = "nanobody-anything"
PROTEIN_PROTOCOL = "protein-anything"
# Every protocol string this collector knows how to interpret. A record naming
# anything else errors its designs rather than being read under the nearest-looking
# protocol's rules -- see the fail-closed note in collect_boltzgen.
KNOWN_PROTOCOLS = frozenset({NANOBODY_PROTOCOL, PROTEIN_PROTOCOL})

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
#   * `plip_hbonds`, `plip_saltbridge` -- use the `_refolded` ones.
# `esmfold2_*` are absent for the same structural reason (that step never runs).
#
# THREE NAMES WERE REMOVED FROM THIS LIST ON 2026-10-09, after a real run showed them
# missing from both shards' tables while this module printed "left NA" about them.
# They are not missing because of a step set; upstream NEVER writes them:
#   * `plip_hydrophobic_refolded` -- `count_noncovalents` (analyze_utils.py:183-225)
#     returns exactly `plip_hbonds` and `plip_saltbridge`. There is no hydrophobic
#     contact metric anywhere in BoltzGen. The name was invented here.
#   * `target_id` -- added to the frame at `analyze.py:450`, which is AFTER
#     `df.to_csv(aggregate_metrics_analyze.csv)` at `:421`. It reaches
#     per_target_metrics_analyze.csv, never this one.
#   * `free_cysteines` -- guarded by `self.free_cys` (`analyze.py:920`), whose
#     constructor default is False (`:124`) and which NO protocol config and NO key in
#     analysis.yaml ever sets. It is unreachable through this tool. (`filter_cysteine`
#     in the nanobody/peptide filtering config is a different thing and `filtering`
#     never runs here.)
# Declaring a column upstream cannot produce is not conservative, it is wrong. Since
# 2026-10-09 a declared-but-absent name is written as pd.NA rather than dropped, so it
# no longer breaks a filter -- but it now produces a permanently all-NA column that
# looks exactly like a metric that was not measured. Every name in this tuple and in
# PROTEIN_PROTOCOL_METRICS is a claim that upstream writes it; see _metric_data.
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
    # --- sequence-level developability -----------------------------------------
    # liability_* uses the modality the PROTOCOL selects, and the two are not
    # comparable. nanobody-anything -> `antibody` rules (cli/boltzgen.py:99-105).
    # protein-anything -> the base default `peptide` (analysis.yaml), which on a
    # folded 80-225 aa protein is close to a lysine+arginine count. Read
    # boltzgen_protocol before comparing these across rows.
    "design_hydrophobicity",
    "design_chain_hydrophobicity",
    "liability_score",
    "liability_num_violations",
    "liability_details",
    # --- composition / provenance ----------------------------------------------
    "num_design",
    "num_tokens",
    "num_prot_tokens",
    "helix",
    "sheet",
    "loop",
)

# Collected ONLY under protein-anything, because only that protocol computes them.
#
# `designfolding-*` come from the `design_folding` step: the designed chain is folded
# ALONE, without the target, and compared back to the designed backbone
# (analyze.py:971-1039). `allatom_fold_metrics` is False here, so only the `bb_`
# family exists; `get_fold_metrics` does NOT prefix its confidence block
# (analyze_utils.py:178-180), which is why `designfolding-design_ptm` has no `bb_`.
#
# **The hyphen is upstream's, not a typo.** `analyze.py:1037-1039` builds the key as
# f"designfolding-{k}". Collected under that spelling so the column matches the raw
# CSV and --all-metrics does not produce a second name for the same number. A --filter
# addresses it as df["boltzgen_designfolding-bb_rmsd_design"]; it cannot be used in a
# df.query() string.
#
# `design_largest_hydrophobic_patch_refolded` is computed because
# `largest_hydrophobic_refolded` is true in analysis.yaml and protein-anything's
# protocol config is `{}`, so nothing turns it off. nanobody-anything DOES turn it off
# (cli/boltzgen.py:99-105), which is why it is not in CORE_METRICS.
PROTEIN_PROTOCOL_METRICS = (
    # Binder self-consistency: designed backbone vs the SAME SEQUENCE folded alone.
    # This is the number that nanobody-anything needs chainsel + usalign to get.
    "designfolding-bb_rmsd_design",
    "designfolding-bb_rmsd",
    # Confidence of the binder folded alone (a monomer prediction).
    "designfolding-design_ptm",
    "designfolding-complex_plddt",
    # Surface hydrophobicity of the refolded binder, A^2.
    "design_largest_hydrophobic_patch_refolded",
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
# protein-anything: there IS no framework, so there is nothing to validate. Its own
# verdict, NOT SCAFFOLD_DISTINCT.
#
# Without it this column lies by accident. A de-novo binder entity gets a
# source_context entry whose `source` is the literal string "inline"
# (schema.py:1607 -> source_context.from_structure, default source="inline"), with no
# scaffold allow-list recorded, so the existing path would return
# ("inline", distinct_from_target) -- a scaffold NAME and a verdict that reads as
# "checked and fine" for a protocol that has no scaffolds at all.
SCAFFOLD_NOT_APPLICABLE = "not_applicable_no_scaffold"  # -> NA


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
    more than one means the spec designed something other than a lone binder. Both
    are errors: the sequence columns and the scaffold would otherwise be attached to a
    chain nobody verified. Protocol-independent -- nanobody-anything and
    protein-anything both design exactly one chain.
    """
    designed = [i for i, c in enumerate(chains) if any(c.designed)]
    if len(designed) == 1:
        return designed[0]
    raise MappingError(
        f"expected exactly one designed chain in the design mmCIF, found "
        f"{len(designed)} ({[chains[i].label for i in designed]} of "
        f"{[c.label for c in chains]}); the design-mask B-factor cannot identify the "
        f"binder chain"
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
    *,
    has_scaffolds: bool = True,
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

    ``has_scaffolds=False`` (``protein-anything``) short-circuits the whole thing:
    the protocol has no framework, so there is no name to publish and no check to
    report the result of. The ``chain_map`` verdict is still computed, because it is
    about the chain LINK and holds for a de-novo entity exactly as for a scaffold one
    (``update_designed`` writes ``output_chain`` for every chain).
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

    if not has_scaffolds:
        # protein-anything. `source` here is "inline" -- the de-novo entity's
        # placeholder, not a file. Publishing it as boltzgen_scaffold would put a
        # scaffold name on a design that has no framework.
        return pd.NA, chain_map, SCAFFOLD_NOT_APPLICABLE

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
    *,
    has_scaffolds: bool = True,
) -> BinderMapping:
    """Identify the binder chain and everything read off it, or raise MappingError.

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
            f"tool's premise is a SINGLE designed binder chain under either protocol, "
            f"so there is no unambiguous binder sequence"
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
        npz_source_context(npz_path),
        index,
        binder.label,
        allowed_stems,
        has_scaffolds=has_scaffolds,
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


# EVERY non-metric column this tool documents, in table order. Seeded to pd.NA on
# EVERY row before anything else is written, so the column exists whatever happened to
# that design.
#
# WHY (added 2026-10-09, closing the hole the pd.NA metric change left open). The
# metric fill covered only the declared METRIC names. The trust and attestation
# columns were still written on their success paths only -- `chain_map` inside the
# post-mapping `data.update(...)`, the `binding_*` family inside `if epitope or
# anti_epitope`. So a collect in which EVERY design errors produced a table with no
# `boltzgen_chain_map` column at all, and the mandatory trust filter this tool's skill
# tells you to run FIRST raised KeyError.
#
# That is not a rare shape. The single most common cause of an all-error collect is a
# `-l/--dir-label` or `-t/--table` resolving to a different out_dir than the run used
# -- which this module's own error message names as the first thing to check. So the
# failure landed precisely on the run the filter was most needed for.
#
# Seeded UNCONDITIONALLY, including the protocol-specific `binder_refold_cif`, so the
# guarantee needs no qualifier: every structural, trust and attestation column exists
# on every row of every table this tool writes. (Metric columns are a separate,
# protocol-dependent set -- see `declared` in collect_boltzgen -- and are NA-filled
# against whichever set applies to that row's protocol.)
ROW_KEYS = (
    "id",
    "protocol",
    "chain_map",
    "scaffold_src",
    "mapping_seq_identity",
    "binding_resnames",
    "not_binding_resnames",
    "binding_frame",
    "binding_chain_frame",
    "chain_id_map",
    "cif_path",
    "binder_refold_cif",
    "binder_chain",
    "sequence",
    "cdr_sequence",
    "scaffold",
    "design_resnums",
    "fixed_resnums",
    "n_design_res",
    "binder_len",
)


def _metric_data(
    row: pd.Series, columns: list[str], absent: list[str]
) -> dict[str, Any]:
    """The declared metric columns, with an absent one written as ``pd.NA``.

    WHY pd.NA AND NOT A DROPPED KEY -- changed 2026-10-09
    -----------------------------------------------------
    It used to be ``{col: row[col] for col in columns if col in row.index}``, with
    ``columns`` already filtered to what the CSV carried. So a declared-but-absent
    name never reached ``Collected.data``, ``prosapia`` never minted the column, and
    the table's shape was a function of **what one CSV happened to contain** rather
    than of the tool version. Measured on the Volume: ``table0_a100s1.tsv`` and
    ``table0_h100.tsv`` of ``20261009_071201_boltzgen_pilot_1tnf_d143`` have 57
    columns each and carry no ``boltzgen_free_cysteines``,
    ``boltzgen_plip_hydrophobic_refolded`` or ``boltzgen_target_id`` under any
    spelling -- while this module printed "left NA" about all three.

    The consequence was not cosmetic: a ``-f`` filter naming such a column raised
    ``KeyError`` instead of filtering, which is worst exactly on the run whose metrics
    went missing. A filter has to be writable against the tool's documented column
    set, so that set must exist in every table the tool writes.

    ``pd.NA`` is also the honest value: it means "not measured", which is what
    happened. **NA is never 0.**

    THE SHARP EDGE THIS INTRODUCES, and what replaces the lost signal
    -----------------------------------------------------------------
    The ``KeyError`` was loud. An all-NA column is not: at a glance it is
    indistinguishable from a metric that was genuinely not produced. And five of the
    ``PROTEIN_PROTOCOL_METRICS`` names are read off ``analyze.py`` and have **never
    been observed in a real CSV**, so a misspelling would now present as a plausible
    all-NA column instead of an error.

    So the printed line in ``collect_boltzgen`` is now the only detector, and it
    **names every absent column** rather than counting them. Do not reduce it to a
    count. An entirely-NA column has two causes -- the metric was not produced, or
    the declared name here is wrong -- and only a ``--all-metrics`` collect, which
    takes the CSV's own column names, tells them apart.
    """
    # Seed FIRST, so a later real value overwrites the NA rather than the reverse.
    data: dict[str, Any] = {key: pd.NA for key in ROW_KEYS}
    data.update({col: row[col] for col in columns if col in row.index})
    # Declared and missing -> NA, so the column EXISTS and a filter naming it does
    # not raise. Written per row rather than once, because prosapia mints a column
    # from the keys of each row's Collected.data.
    for col in absent:
        data[col] = pd.NA
    return data


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

        record, record_error = group_record(d.name)

        # WHICH PROTOCOL THESE DESIGNS ARE, read from the group's own .submit.json --
        # which sits beside the designs it describes -- and NOT from a flag or from
        # .meta.json. A run_dir can hold nanobody and minibinder groups side by side,
        # and .meta.json is merged across submits, so a run-level value would label
        # both with whichever ran last. That is the exact failure this column exists
        # to prevent.
        #
        # Taken from the raw record even when `record_error` is set. The error says the
        # EPITOPE attestation does not describe these designs; the protocol is a
        # separate, run-level fact and withholding it would leave an error row unable
        # to say what class of binder it even was.
        #
        # NA when the record is missing/unreadable (an error row anyway), or when it
        # predates this column -- the 2026-10-08/09 pilot runs. NA is the absence of a
        # record, never a protocol.
        protocol = (record or {}).get("protocol") or None
        if protocol is None and record:
            print(
                f"{d.name}: {SUBMIT_RECORD} records no 'protocol'; it was written "
                f"before boltzgen_protocol existed. The column is NA for this group "
                f"-- it is NOT a claim that the protocol was nanobody-anything."
            )
        elif protocol is not None and protocol not in KNOWN_PROTOCOLS:
            # A protocol string this collector does not understand means its branches
            # -- the scaffold question, the designfolding columns, the cdr_sequence
            # suppression -- were written for a different experiment. Carried into the
            # table verbatim (the table should say what the record says) but every
            # design of the group ERRORS: there is no safe default reading.
            record_error = record_error or (
                f"{SUBMIT_RECORD} records protocol {protocol!r}, which this collector "
                f"does not know. Its scaffold, sequence and designfolding columns are "
                f"written per protocol and none of them can be interpreted for this "
                f"one. Collected as errors rather than under the nearest-looking "
                f"protocol's rules."
            )
            print(f"{d.name}: {record_error}")

        # FAIL CLOSED. `protocol == NANOBODY`, never `protocol != PROTEIN`.
        #
        # The inverted form treated an ABSENT or UNRECOGNISED protocol as nanobody,
        # and nanobody is the branch that PUBLISHES A PROVENANCE VERDICT. A de-novo
        # design landing there gets scaffold "inline" with scaffold_src
        # `distinct_from_target` -- a verdict that reads as *checked and passed* -- on
        # a protocol that has no scaffolds at all. The safe reading of "I cannot tell
        # which protocol this was" is MAKE NO SCAFFOLD CLAIM, never "assume nanobody":
        # a withheld name costs a diversity analysis, a wrong one corrupts it.
        #
        # So only an explicit, recognised nanobody record enables the scaffold check.
        # Everything else -- protein-anything, an absent protocol (a pre-2026-10-09
        # record), an unknown string, a missing record -- reports
        # `not_applicable_no_scaffold` and NA.
        has_scaffolds = protocol == NANOBODY_PROTOCOL
        if protocol is None:
            print(
                f"{d.name}: no protocol recorded, so the scaffold check is NOT run "
                f"(scaffold NA, scaffold_src '{SCAFFOLD_NOT_APPLICABLE}'). This is "
                f"fail-closed: assuming nanobody would publish a provenance verdict "
                f"nothing attests."
            )

        # Columns depend on the protocol: only protein-anything runs design_folding and
        # the largest-hydrophobic-patch metric, so only it can produce those columns.
        declared = (
            (*CORE_METRICS, *PROTEIN_PROTOCOL_METRICS, *extra_metrics)
            if protocol == PROTEIN_PROTOCOL
            else (*CORE_METRICS, *extra_metrics)
        )
        if ctx.args.all_metrics:
            # --all-metrics takes the CSV's OWN column names, so nothing can be
            # "declared but absent" and there is nothing to fill with NA. This is also
            # the diagnostic for an all-NA column below: it shows what is really there.
            columns = [
                c
                for c in df.columns
                if c not in (ID_COL, FILE_COL, CHAIN_SEQ_COL) and _SAFE_COLUMN.match(c)
            ]
            absent: list[str] = []
        else:
            columns = [c for c in declared if c in df.columns]
            # NOT filtered by warned_absent: every absent name must be written as NA on
            # every row of every group, or the column's presence would depend on which
            # group happened to be collected first. Only the PRINTING is deduped.
            absent = [c for c in declared if c not in df.columns]
            unwarned = [c for c in absent if c not in warned_absent]
            if unwarned:
                warned_absent.update(unwarned)
                # This line is the ONLY remaining detector for a misspelled declared
                # name -- since 2026-10-09 these columns are written as NA rather than
                # dropped, so they no longer raise KeyError in a filter. It names every
                # column deliberately; do not reduce it to a count.
                print(
                    f"  not in this CSV, so written as NA (the column EXISTS and a -f "
                    f"filter naming it will not raise; an ALL-NA column means either "
                    f"the metric was not produced OR the name declared here is wrong "
                    f"-- re-collect with --all-metrics to see the CSV's real names): "
                    f"{', '.join(sorted(unwarned))}"
                )

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
        if not has_scaffolds:
            print(
                f"{d.name}: protocol {protocol} has no framework; boltzgen_scaffold "
                f"is NA and scaffold_src is '{SCAFFOLD_NOT_APPLICABLE}' for every row"
            )
        else:
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
            data: dict[str, Any] = _metric_data(row, columns, absent)
            data["id"] = design_id
            # On EVERY row, including every error row. A row that cannot say what
            # class of binder it is is the problem this column was added to solve.
            data["protocol"] = protocol if protocol else pd.NA
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
            # The BINDER ALONE. Only protein-anything writes it (the design_folding
            # step), so it is required only there. Required rather than optional:
            # writer.py publishes the npz LAST, and analysis would have dropped the
            # design entirely without the npz (analyze.py:972-977), so a row in the
            # CSV means the cif should exist. Its absence is an anomaly, not a
            # protocol difference, and an anomaly is an error row.
            binder_refold_cif = (
                design_dir / REFOLD_DESIGN_DIRNAME / f"{design_id}.cif"
            )
            required = [design_cif, npz_path, refold_cif]
            if protocol == PROTEIN_PROTOCOL:
                required.append(binder_refold_cif)

            # The names are reported RELATIVE TO design_dir, not as bare basenames:
            # <id>.cif, refold_cif/<id>.cif and refold_design_cif/<id>.cif all share
            # a basename, so "missing ['d1.cif']" would not say which one.
            missing = [
                str(p.relative_to(design_dir)) for p in required if not p.is_file()
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
            if protocol == PROTEIN_PROTOCOL:
                # The one real new capability of protein-anything: a binder-alone
                # structure, so self-consistency does not need a chainsel + fold
                # round-trip. Left as mmCIF (usalign reads it); the complex is
                # converted to PDB because boltzgen_path is the column downstream
                # tools default to.
                data["binder_refold_cif"] = str(binder_refold_cif)
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
                    has_scaffolds=has_scaffolds,
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

            # THE FOURTH WITNESS, and it is free -- but only under protein-anything.
            #
            # That protocol designs the WHOLE chain (the spec's `protein:` entity is
            # all designed residues, schema.py:528-556 marks every one), so the
            # designed-residue count and the chain length are two independent
            # measurements of the same object: one from the design mask, one from the
            # refolded chain's residues. They must be equal.
            #
            # An inequality means the mask and the chain do not describe the same
            # thing, which displaces EVERY residue-indexed column (design_resnums,
            # fixed_resnums) and makes binder_len a length for some other chain. The
            # three existing witnesses all count designed residues, so they agree with
            # each other and miss this entirely. An error row, never a warning: a
            # count mismatch is an error.
            if (
                protocol == PROTEIN_PROTOCOL
                and len(mapping.design_positions) != mapping.binder_len
            ):
                msg = (
                    f"under {PROTEIN_PROTOCOL} the whole binder chain is designed, so "
                    f"the designed-residue count must equal the chain length -- but "
                    f"{len(mapping.design_positions)} residues are masked as designed "
                    f"in a chain of {mapping.binder_len}. The design mask and the "
                    f"chain do not describe the same object, so every residue-indexed "
                    f"column would be displaced"
                )
                print(f"{d.name}: {design_id} {msg}")
                # CARRY WHAT WAS GENUINELY MEASURED, exactly as MappingError.data does
                # for the failures raised inside map_binder. `map_binder` ran to
                # completion here -- this witness fires AFTER it -- so three results
                # exist and discarding them would make the one row that needs
                # diagnosing the least diagnosable row in the table. The first
                # question anyone asks of this error is "was the chain mapping sound
                # and only the count wrong, or was the whole mapping adrift?", and
                # mapping_seq_identity is the column that answers it.
                #
                #   * chain_map / mapping_seq_identity -- both concern the chain LINK
                #     and the CSV-vs-structure sequence identity, neither of which is
                #     a count. They were established before this check and are
                #     unaffected by it.
                #   * scaffold_src -- published because under THIS protocol it is not
                #     a positional lookup at all: scaffold_of short-circuits to
                #     `not_applicable_no_scaffold` before it ever indexes
                #     source_context["chains"] (there is no framework to validate), so
                #     it is a statement about the protocol, not about this design's
                #     chain index. `scaffold` is NA beside it, as always here.
                #
                # Everything residue-indexed or count-derived stays NA: binder_chain,
                # sequence, cdr_sequence, design_resnums, fixed_resnums, n_design_res,
                # binder_len. Those are precisely what the error is about, and naming
                # a chain or a residue on a row whose mapping is in doubt is the thing
                # this collector refuses to do anywhere else.
                data.update(
                    {
                        "chain_map": mapping.chain_map,
                        "mapping_seq_identity": mapping.seq_identity,
                        "scaffold_src": mapping.scaffold_src,
                        "scaffold": mapping.scaffold,
                    }
                )
                yield Collected(
                    name=row_name,
                    parent=d.name,
                    path="",
                    status=f"error: MappingError: {msg}",
                    data=data,
                )
                n += 1
                continue

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
                    # NANOBODY-ONLY BY MEANING, written as NA elsewhere rather than
                    # omitted, so the column exists under both protocols and a -f
                    # filter naming it cannot raise KeyError.
                    #
                    # Under protein-anything the whole chain is designed, so upstream's
                    # `designed_sequence` equals `designed_chain_sequence` and putting
                    # it in a column called cdr_sequence would assert a CDR/framework
                    # split that does not exist. boltzgen_sequence already carries it.
                    "cdr_sequence": (
                        pd.NA
                        if protocol == PROTEIN_PROTOCOL
                        else _cell(row, DESIGN_SEQ_COL)
                    ),
                    # NA under protein-anything via scaffold_of's has_scaffolds branch;
                    # scaffold_src says `not_applicable_no_scaffold` beside it.
                    "scaffold": mapping.scaffold,
                    # Honest under BOTH protocols, but they mean different things.
                    # nanobody: the CDRs, and fixed_resnums is the framework that
                    # `proteinmpnn --fixed-positions` consumes. protein-anything: the
                    # whole chain is designed, so design_resnums is `1:binder_len` and
                    # fixed_resnums is EMPTY -- which is correct, not a bug, and means
                    # there is nothing to hand proteinmpnn.
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
        if n == 0:
            # A header-only aggregate_metrics_analyze.csv is NOT a successful collect
            # of nothing, and must not print the word OK. It is the DOCUMENTED
            # signature of a protein-anything run whose design_folding step did not
            # happen: analysis is launched with designfolding_metrics=true from the
            # PROTOCOL (cli/boltzgen.py:1237,1310) while --steps only filters the step
            # list, so analyze.py:972-977 returns None for every design and writes a
            # CSV with columns and no rows -- after the whole GPU spend, with no
            # traceback anywhere.
            print(
                f"{d.name}: NO DESIGNS -- {METRICS_CSV} exists but has zero rows. "
                f"This is not an empty success. Under {PROTEIN_PROTOCOL} it is the "
                f"signature of a run whose design_folding step did not produce "
                f"fold_out_design_npz/: analysis then skips every design and writes a "
                f"header-only CSV. Check the task .out for repeated 'Folded path does "
                f"not exist. Skipping', and check that the submit's step list included "
                f"design_folding."
            )
        else:
            print(f"{d.name}: OK ({n} design(s))")
        named = sum(
            c for v, c in scaffold_verdicts.items()
            if v in (SCAFFOLD_VERIFIED, SCAFFOLD_DISTINCT)
        )
        if scaffold_verdicts and not named and has_scaffolds:
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
