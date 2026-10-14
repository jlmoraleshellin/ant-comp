#!/usr/bin/env python3
"""
Submit BoltzGen binder-design runs, one array task per target.

PREMISE -- read this before using any column this tool writes
------------------------------------------------------------
BoltzGen is campaign-driven, not design-driven: one ``boltzgen run <spec>.yaml``
process generates N designs against one target and only then are there rows. So the
unit of work is **one BoltzGen run per target**, the row count is known at collect
time, and this is a ``create`` tool.

TWO PROTOCOLS, chosen by ``--protocol``, and they are different experiments
---------------------------------------------------------------------------
``nanobody-anything`` (**the default, unchanged**) grafts new CDRs onto a fixed,
randomly drawn VHH framework. ``protein-anything`` designs a **de-novo minibinder**
chain of a length you choose, with no framework at all. Everything the two share --
the epitope resolution, the chain-id translation, the ``label_seq_id`` residue frame,
the single-designed-chain mapping guard, the ``boltzgen check`` pre-flight -- is
protocol-independent and runs identically. Everything that differs is listed here and
is marked per column in the skill.

``boltzgen_protocol`` is a COLLECTED COLUMN, read back from the group's own
``.submit.json``, so every row states which class of binder it is. That column exists
because a 60-design nanobody run was once handed to a caller who wanted minibinders
and nothing in the table said which they had got. **Read it before any other column.**

**Under ``--spec`` it attests the FLAG, not the YAML.** A hand-written spec decides
what is actually designed -- its entities are used verbatim -- while ``--protocol``
only selects BoltzGen's step list and config defaults. The two can disagree and
nothing here can detect it: the column records what was requested, which is the only
thing this tool knows. (The resubmit guard compares the recorded protocol explicitly
for exactly this reason: a ``--spec`` YAML hashes identically under either protocol,
so the sha alone would let a second submit relabel the first's designs.)

The premise, and the conditions under which the numbers are meaningless:

* **The binder is ONE designed chain**, under both protocols. Under
  ``nanobody-anything`` the framework sequence AND structure come from a scaffold
  ``.cif`` and only CDR1/2/3 are designed; under ``protein-anything`` the whole chain
  is designed from nothing. If a spec is supplied (``--spec``) whose designed region
  is not a single chain, the collector refuses the design rather than guessing which
  chain the sequence belongs to.
* **``protein-anything`` draws ONE binder length per BATCH, not per design.**
  Upstream's own ``--diffusion_batch_size`` help says so in as many words: "for design
  tasks that randomly sample the binder length ... all designs generated in the same
  batch will share the same length". Mechanically, ``data.cfg.multiplicity`` is
  ``ceil(num_designs / diffusion_batch_size)`` (``cli/boltzgen.py:1056-1059``), each
  dataset item re-parses the YAML (``data_from_yaml.py:247-266``) and the length is
  drawn once per parse (``schema.py:546``). So a wide ``--binder-length`` with the
  default batch size of 10 (which upstream picks for ``--num-designs >= 100``) samples
  only ``num_designs / 10`` distinct lengths. **Pass ``--diffusion-batch-size 1`` when
  the length range matters.** Same trap as rfd3's ``--length min-max``.
* **``protein-anything`` scores developability with PEPTIDE liability rules.**
  ``protocol_configs["protein-anything"]`` is ``{}`` (``cli/boltzgen.py:85``) -- the
  base config -- and ``analysis.yaml`` sets ``liability_modality: peptide``. Nothing
  carries the nanobody protocol's ``antibody`` rules across; the two protocols are
  separate dict entries with no inheritance. But the peptide motif set
  (``analyze_utils.py:921-941``) is ``ProtTryp`` on every K/R, an N-terminal
  cyclization motif and a DPP4 cleavage motif -- rules for a short synthetic peptide.
  On a folded 80-225 aa protein ``liability_score`` is essentially a count of lysines
  and arginines. **It is collected under upstream's spelling and it means nothing
  here.** Use ``--config-set 'analysis liability_analysis=false'`` if the noise is
  unwanted.
* **Residue indices in ``--binding``/``--not-binding`` are 1-based ``label_seq_id``,
  NOT ``auth_seq_id``** (upstream README, "IMPORTANT"). A crystallographic target
  numbered 18..155 is 1..138 to BoltzGen. This builder resolves each requested
  position against the target structure at SUBMIT time and records the residue NAME
  it found, so a silently wrong epitope shows up as the wrong amino acids in
  ``boltzgen_binding_resnames`` rather than as a plausible but misplaced design. A
  position that does not exist is an error, never a dropped hotspot.
* **Under ``nanobody-anything`` the scaffold is drawn at RANDOM per design**
  (``data/parse/schema.py:1947-1949`` does ``path = random.choice(path)``). That is
  wanted -- it is the diversity axis -- and there is no seed to pin it with (BoltzGen
  has no ``--seed`` flag at all). Which scaffold a given design actually used is
  recovered per design by the collector into ``boltzgen_scaffold``. **Designs are
  therefore not reproducible, deliberately.** Do not "fix" this. Under
  ``protein-anything`` there is NO scaffold: ``--scaffolds`` and ``--scaffold-dir``
  are REFUSED at submit time rather than ignored, ``boltzgen_scaffold`` is NA and
  ``boltzgen_scaffold_src`` reports ``not_applicable_no_scaffold`` -- never
  ``distinct_from_target``, which would read as a verdict that had been checked.
* **Target geometry is not guaranteed.** Open issue #226: on a multi-chain target a
  second protomer can MOVE during design, giving a complex whose binding geometry is
  wrong. Unanswered upstream. ``boltzgen_bb_rmsd_target`` is the column that catches
  it and it must be read before any interface number. Protocol-independent.

ESMFold2 is EXCLUDED, by construction and not by accident -- under BOTH protocols
---------------------------------------------------------------------------------
Every run gets ``--config analysis esmfold2_metrics=false`` and a step list that
omits ``esmfold2_scoring`` and ``filtering``. Three separate facts force it, verified
against the pinned SHA, and **all three were re-checked for ``protein-anything``**:

* generation never touches esmfold2 (``task/predict/`` has no reference to it);
* ``analysis`` is launched with ``esmfold2_metrics={not use_affinity}`` = true and
  ``analyze.py:638-652`` then unconditionally reads
  ``esmfold2_scores/<id>.input.json``, so it crashes unless that key is set false.
  ``use_affinity`` is ``protocol in ["protein-small_molecule"]``
  (``cli/boltzgen.py:1259``), so it is False for BOTH protocols here;
* ``filtering`` sets ``self.use_esmfold2 = not use_affinity`` at ``filter.py:218``
  -- **hard-wired to the protocol, reachable by no config key** -- and raises
  ``ValueError("Missing ESMFold2 scores; ...")``. It can only be omitted. Same
  ``use_affinity``, so this holds for ``protein-anything`` too.

THE STEP LIST IS PROTOCOL-DEPENDENT, and getting it wrong EMPTIES the CSV
-------------------------------------------------------------------------
``nanobody-anything``: ``design inverse_folding folding analysis``.
``protein-anything``: the same plus ``design_folding``.

That extra step is **not optional and not a nicety**. ``cli/boltzgen.py:1237`` sets
``do_design_folding = protocol in ["protein-anything", "protein-small_molecule"]``
and then passes ``designfolding_metrics={do_design_folding}`` to the analysis step at
``:1310``. ``filter_steps`` (``:1363``) only drops steps from the list -- it does not
rebuild their arguments -- so running ``protein-anything`` WITHOUT ``design_folding``
leaves ``designfolding_metrics=true`` on an analysis step with no
``fold_out_design_npz/`` to read. ``analyze.py:972-977`` then prints "Folded path does
not exist. Skipping" and ``return None`` for **every** design, so
``aggregate_metrics_analyze.csv`` comes back with zero rows after all the GPU work.
Not a crash -- an empty table. This is the same class of failure as the esmfold2 one
and the reason ``--steps`` is not a user flag.

What that step BUYS is a real new capability: it writes
``refold_design_cif/<id>.cif`` -- **the binder folded ALONE, without the target**
(``task/predict/writer.py:94-104``, ``const.refold_design_cif_dirname``) -- plus the
``designfolding-*`` metrics. So under ``protein-anything`` binder self-consistency is
in the table directly (``boltzgen_designfolding-bb_rmsd_design``) and does not need a
separate ``chainsel`` + fold round-trip. Under ``nanobody-anything`` it is still
absent (issue #148), and ``chainsel`` + ``usalign`` is still the route.

``--steps`` is therefore NOT a flag of this tool: anyone who could set it could only
set it wrong. ``--esmfold2_acceleration off`` does not help -- it selects a numerical
backend, it does not disable the model. ``--extra-flags`` and ``--config-set`` are
**guarded** (``BLOCKED_EXTRA_FLAGS``, ``check_config_set``): because the task script
interpolates them after the forced flags and argparse is last-wins, an unguarded
escape hatch would let ``--extra-flags '--protocol protein-anything'`` select a
protocol whose step list this tool never wrote, and the run would collect nothing
while ``boltzgen_protocol`` still said nanobody. ``--protocol`` is a FIRST-CLASS FLAG
precisely so that the step list, the spec shape and the recorded column move together.

The spec this tool writes is also validated before any GPU work: ``boltzgen.sh`` runs
``boltzgen check`` on it first and fails the task if it does not pass. That is the only
independent test that the YAML nesting below is what BoltzGen's schema expects -- a
parser that ignores an unknown key would otherwise turn a shape error into "no epitope
was requested" while ``boltzgen_binding_resnames`` still reported one.

Consequence: there is no ``final_ranked_designs/``, no ``results_overview.pdf`` and
no ``--budget``. Ranking happens in the table, with ``-f``, which is where it belongs.

Chain and residue specs
-----------------------
``--target-chains`` uses the prosapia chain mini-language (``,`` separates, ``:`` is
an inclusive letter range): ``A``, ``A,B,C``, ``A:C``. ``--binding`` /
``--not-binding`` use the workspace's chain-prefixed residue syntax, the same one
``epitope`` and ``bindcraft2`` take: ``A91,A128,B91-95``, with ``{expr}`` islands
resolved per design against the table lineage (``A{epitope_start}-{epitope_end}``).

``--binder-length`` (``protein-anything`` only) takes ``N`` or ``MIN-MAX`` or
BoltzGen's own ``MIN..MAX``, and ``{expr}`` islands resolve in it too. It becomes the
``sequence:`` of a ``protein:`` entity. It is **refused** under ``nanobody-anything``,
where the framework fixes the length and the CDR lengths are sampled.

Usage:
    # root run: one BoltzGen run against a target not in any table yet
    sapia run boltzgen outputs/RUN -e modal \\
        --target-structure inputs/trimer.cif --target-chains A,B,C \\
        --binding A91,A128,B91,B128 --num-designs 50

    # de-novo minibinders, 80-225 aa, against the same target. The batch size is
    # NOT decoration: without it a 200-design run samples only 20 lengths.
    sapia run boltzgen outputs/RUN -e modal \\
        --protocol protein-anything --binder-length 80-225 \\
        --diffusion-batch-size 1 \\
        --target-structure inputs/trimer.cif --target-chains A,B,C \\
        --binding A91,A128,B91,B128 --num-designs 200

    # child run: one BoltzGen run per target structure already in a table
    # -i is the column holding the TARGET structure. default_input_column is
    # 'pdb_path'; point it at whatever actually holds the target in your table
    # (chainsel_path for a cropped target, mkcomplex_path, ...).
    sapia run boltzgen outputs/RUN -t table0 -i chainsel_path -e modal \\
        --binding 'A{epitope_start}-{epitope_end}' --num-designs 200

    # a spec YAML you wrote yourself, used verbatim (still {expr}-resolved)
    sapia run boltzgen outputs/RUN -e modal --spec specs/mytarget.yaml
"""

import hashlib
import json
import os
import re
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, NamedTuple, cast

import gemmi
import yaml

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.data_manager import LookupFn
from prosapia.core.executors import volume_path
from prosapia.utils import expand_chain_spec, resolve_template

# Layout this tool imposes on its out_dir. collect_boltzgen.py reads the same
# convention -- keep the two in step.
RUNS_DIRNAME = "runs"
SPECS_DIRNAME = "specs"

# NOTE the BoltzGen output layout the collector reads (intermediate_designs_inverse_
# folded/, refold_cif/, aggregate_metrics_analyze.csv) is declared ONCE, in
# collect_boltzgen.py. It was duplicated here and the copy went unused; two copies of a
# layout contract is one copy too many.

# The two protocols this tool supports, and the step list each one forces.
#
# NANOBODY IS THE DEFAULT AND ITS LIST IS BYTE-FOR-BYTE WHAT IT ALWAYS WAS. The
# 2026-10-08 real-data verification and the 2026-10-09 60-design run were made with
# exactly `--protocol nanobody-anything --steps design inverse_folding folding
# analysis`, and that is what an invocation with no --protocol still produces.
#
# protein-anything adds `design_folding` and it is NOT optional: upstream derives
# `designfolding_metrics` from the PROTOCOL (cli/boltzgen.py:1237,1310) while
# `--steps` only filters the step list (`:1363`), so omitting it leaves analysis
# reading a fold_out_design_npz/ that was never written and `analyze.py:972-977`
# returns None for every design -- an EMPTY aggregate_metrics_analyze.csv after all
# the GPU work. See the module docstring.
NANOBODY_PROTOCOL = "nanobody-anything"
PROTEIN_PROTOCOL = "protein-anything"
DEFAULT_PROTOCOL = NANOBODY_PROTOCOL

_BASE_STEPS = ("design", "inverse_folding", "folding", "analysis")
PROTOCOL_STEPS: dict[str, tuple[str, ...]] = {
    NANOBODY_PROTOCOL: _BASE_STEPS,
    # design_folding sits between folding and analysis in upstream's own ordering;
    # `--steps` is a membership set so the order of the tokens does not matter, but
    # keep it readable.
    PROTEIN_PROTOCOL: (
        "design",
        "inverse_folding",
        "folding",
        "design_folding",
        "analysis",
    ),
}
PROTOCOLS = tuple(PROTOCOL_STEPS)

# Kept for anything that imported it before --protocol existed; it is the DEFAULT, not
# "the" protocol. Prefer args.protocol.
PROTOCOL = DEFAULT_PROTOCOL
STEPS = PROTOCOL_STEPS[DEFAULT_PROTOCOL]
FORCED_CONFIG = (("analysis", "esmfold2_metrics=false"),)

# The five VHH scaffolds example/nanobody/penguinpox.yaml uses, verbatim and in its
# order. example/nanobody_scaffolds/ also ships 8coh, 8z8v and gefurulimab; they are
# reachable by naming them in --scaffolds. A scaffold YAML references its own .cif
# RELATIVELY ("file references are relative to the location of the .yaml file"), so
# the YAMLs must be pointed at IN PLACE -- never copied into the run_dir.
DEFAULT_SCAFFOLDS = (
    "7eow",
    "7xl0",
    "gontivimab",
    "isecarosmab",
    "sonelokimab",
)
# Inside the Modal image (see modal_image.py). It does not exist where `sapia run`
# executes, so scaffold names are NOT checked for existence at submit time.
DEFAULT_SCAFFOLD_DIR = "/opt/boltzgen/example/nanobody_scaffolds"

# A {expr} placeholder island, resolved per-design up the table lineage. Meaningless
# in a root run (no table), so we reject it there.
_HAS_PLACEHOLDER = re.compile(r"\{[^}]*\}")

# Chain-prefixed residue syntax, as in epitope/bindcraft2: `A96` / `A96-99`.
_SINGLE = re.compile(r"^([A-Za-z])(\d+)$")
_RANGE = re.compile(r"^([A-Za-z])(\d+)-([A-Za-z]?)(\d+)$")

# --binder-length. `120`, `80-225` (workspace convention, as in rfd3's --length) and
# BoltzGen's own `80..225` are all accepted; nothing else is.
_LEN_SINGLE = re.compile(r"^(\d+)$")
_LEN_RANGE = re.compile(r"^(\d+)\s*(?:-|\.\.)\s*(\d+)$")

# Chain ids to draw a de-novo binder's chain from when --binder-chain is `auto`.
# Reversed so a one- or two-chain target gets 'Z', visibly not a target chain, rather
# than 'A'. Any letter already written by the target is skipped -- see
# pick_binder_chain for why that matters.
_BINDER_CHAIN_POOL = "ZYXWVUTSRQPONMLKJIHGFEDCBA"


# `boltzgen run` flags this tool owns. argparse is last-wins and --extra-flags is
# interpolated AFTER them, so without this list a user could quietly redefine what the
# tool is: `--extra-flags '--protocol protein-anything'` changes the protocol while
# every docstring, the skill, the column prefix and the single-designed-chain guard
# still say nanobody. Refused at submit time, by name.
BLOCKED_EXTRA_FLAGS = (
    "--protocol",
    "--steps",
    "--output",
    "--num_designs",
    "--config",
    "--cache",
    # --overwrite / --reuse would FALSIFY THE RESUBMIT GUARD'S OWN PREMISE. That guard
    # refuses a top-up into a group that already holds designs, and its stated reason
    # is "BoltzGen is never run with --overwrite, so those designs would survive". If
    # --extra-flags can supply --overwrite, the surviving designs the guard protects
    # are instead silently destroyed -- and --reuse changes which steps are skipped on
    # a re-run, so the outputs the collector reads need not come from this spec at all.
    "--overwrite",
    "--reuse",
    # There is a dedicated --diffusion-batch-size, and argparse is LAST-WINS, so
    # `--extra-flags '--diffusion_batch_size 10'` silently overrides it. The submit log
    # computes its per-batch length warning from the TOOL flag, so the caller is told
    # one thing and gets another -- on the exact parameter that decides whether a
    # length range is sampled per design or per batch of ten.
    "--diffusion_batch_size",
)
# The one config key the tool's contract depends on: analysis crashes when it is true
# (analyze.py:638-652 reads esmfold2_scores/<id>.input.json unconditionally).
PROTECTED_CONFIG_KEY = "esmfold2_metrics"


# How a chain's residue positions were established. Collected as
# boltzgen_binding_frame, per design group.
FRAME_DECLARED = "declared_sequence"
FRAME_SYNTHESIZED = "synthesized_from_model"
# Which identifier a YAML `chain: id:` is resolved against -- and it depends on the
# INPUT FORMAT. Collected as boltzgen_binding_chain_frame.
#   mmCIF: parse_mmcif keeps gemmi's subchain id, i.e. label_asym_id
#          (mmcif.py:1055,1071 -> chains["name"], resolved at schema.py:2301).
#   PDB:   parse_pdb goes through make_mmcif_block() -- which renames A to Axp -- and
#          then patch_chain_names() takes chain["name"][0], the FIRST CHARACTER, with a
#          collision pool (pdb_parser.py:102,107-126). So a PDB's chain ids are its own
#          one-letter names, as long as no two subchains share a first character.
CHAIN_FRAME_LABEL = "label_asym_id"   # mmCIF input: translate author -> label
CHAIN_FRAME_AUTHOR = "author_letter"  # PDB input: patch_chain_names gives it back


class Residue(NamedTuple):
    """One position of a chain's DECLARED polymer sequence.

    ``modelled`` is False for a position that the file declares (``_entity_poly_seq``
    / SEQRES) but gives no coordinates for -- BoltzGen parses it with
    ``is_present=False``, so it occupies a ``label_seq_id`` slot while being unable to
    contact anything.
    """

    name: str
    modelled: bool


class TargetFrames(NamedTuple):
    """Everything about a target's chains that the spec and the attestation need.

    ``fmt`` is gemmi's DETECTED format (``st.input_format``), not a suffix test: the
    whole chain-id translation pivots on it, and ``pdb1abc.ent``, ``target.pdb.gz`` and
    ``target.PDB1`` all look like mmCIF to a ``.suffix`` comparison.

    ``first_chars`` maps a first character to every WRITTEN chain name starting with it,
    over **all** subchains including non-polymers -- that is what upstream's
    ``patch_chain_names`` keys on for a PDB input, so a ligand subchain sharing the
    polymer's first letter has to be visible here.

    ``nonpolymer`` maps an author chain to its non-polymer subchains, so the builder can
    say which cofactors an explicit ``--target-chains`` will drop.
    """

    chains: dict[str, "ChainFrame"]
    fmt: str
    first_chars: dict[str, list[str]]
    nonpolymer: dict[str, list[str]]


class ChainFrame(NamedTuple):
    """One polymer chain of a target, in the frames BoltzGen actually uses.

    ``label_asym`` is the id a YAML ``chain: id:`` is resolved against; ``auth_asym``
    is the letter a structure viewer shows. They are routinely different -- on 1TUP the
    author chains ``E,F,A,B`` are label chains ``A,B,C,D``.
    """

    label_asym: str
    auth_asym: str
    residues: dict[int, Residue]
    res_frame: str


class SpecConfigError(ValueError):
    """A run-wide spec misconfiguration that applies to every target (e.g. --spec
    together with the spec-building flags). Unlike a per-row error it is NOT swallowed
    by the warn-and-skip loop -- it fails the whole submit up front."""


class BoltzGenArgs(CommonArgs):
    protocol: str
    target_structure: Path | None
    target_chains: str
    binding: str | None
    not_binding: str | None
    binder_length: str | None
    binder_chain: str | None
    # None means "not passed". The defaults are applied in resolved_scaffolds() rather
    # than by argparse, so a flag that is MEANINGLESS under the chosen protocol can be
    # refused for having been passed at all -- comparing a value against its default
    # cannot tell `--scaffolds <the default string>` from silence.
    scaffolds: str | None
    scaffold_dir: str | None
    num_designs: int
    diffusion_batch_size: int | None
    num_workers: int
    spec: str | None
    extra_spec: str | None
    config_set: list[str]
    extra_flags: str


def add_run_boltzgen_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--protocol",
        type=str,
        default=DEFAULT_PROTOCOL,
        choices=list(PROTOCOLS),
        help="Which BoltzGen protocol to run. Default %(default)s -- new CDRs on a "
        "fixed, randomly drawn VHH framework. 'protein-anything' designs a DE-NOVO "
        "minibinder chain instead, with no framework: it requires --binder-length and "
        "REFUSES --scaffolds/--scaffold-dir. The two are different experiments and "
        "every row records which one it was, in boltzgen_protocol.",
    )
    parser.add_argument(
        "--binder-length",
        type=str,
        default=None,
        help="protein-anything ONLY: the length of the de-novo binder chain, as 'N' "
        "(exactly N residues), 'MIN-MAX', or BoltzGen's own 'MIN..MAX'. Inclusive at "
        "both ends. {expr} islands resolve per design up the lineage. **A RANGE DRAWS "
        "ONE LENGTH PER BATCH, NOT PER DESIGN** -- upstream's --diffusion_batch_size "
        "help says so, and upstream picks a batch size of 10 for --num-designs >= "
        "100, so a 200-design run would sample only 20 lengths. Pass "
        "--diffusion-batch-size 1 when the range matters. REFUSED under "
        "nanobody-anything, where the framework fixes the length and the CDR lengths "
        "are sampled from the scaffold instead.",
    )
    parser.add_argument(
        "--binder-chain",
        type=str,
        default=None,
        help="protein-anything ONLY: the chain id to give the de-novo binder entity. "
        "Default: the first letter of 'ZYX...' that the target does not already use, "
        "so BoltzGen never has to rename a chain. Pass a letter to pin it. A letter "
        "the target already uses is REFUSED -- BoltzGen would silently rename one of "
        "them (schema.py:1633-1638) and the chain labels in the output would not be "
        "the ones written here.",
    )
    parser.add_argument(
        "--target-structure",
        type=Path,
        default=None,
        help="Target structure (.cif or .pdb) for a ROOT run (no --table). With "
        "--table the target comes from the table's --input-column instead. The "
        "design group is named `<stem>_bg`.",
    )
    parser.add_argument(
        "--target-chains",
        type=str,
        default="all",
        help="Which chains of the target structure to include, in the prosapia chain "
        "mini-language (',' separates, ':' is an inclusive letter range): 'A', "
        "'A,B,C', 'A:C'. Defaults to 'all' (every chain in the file, which is also "
        "BoltzGen's own behaviour when `include:` is omitted).",
    )
    parser.add_argument(
        "--binding",
        type=str,
        default=None,
        help="Epitope the nanobody should bind, chain-prefixed: 'A91,A128,B91-95'. "
        "May embed {expr} islands resolved per design up the lineage. **Residue "
        "numbers are 1-based positions in the target's DECLARED polymer sequence "
        "(mmCIF _entity_poly_seq / PDB SEQRES) -- i.e. label_seq_id, NOT auth_seq_id, "
        "and NOT the count of modelled residues.** The builder resolves each one the "
        "way BoltzGen's own parser does and records the residue name it found in "
        "boltzgen_binding_resnames. A position that is absent, or that is declared but "
        "has NO COORDINATES, is an error at submit time. "
        "Note upstream issue #215: `binding` is a soft preference, not a constraint.",
    )
    parser.add_argument(
        "--not-binding",
        type=str,
        default=None,
        help="Target regions the nanobody should AVOID, same syntax and same "
        "label_seq_id frame as --binding. Upstream issue #217 reports that `binding` "
        "and `not_binding` may not work together on the SAME chain; put them on "
        "different chains, or verify on a small run first.",
    )
    parser.add_argument(
        "--scaffolds",
        type=str,
        default=None,
        help="nanobody-anything ONLY: comma-separated VHH scaffold YAMLs BoltzGen "
        "draws from, AT RANDOM, one per design. Bare names resolve against "
        "--scaffold-dir; a name containing '/' is used as a path verbatim. Defaults "
        "to the five in example/nanobody/penguinpox.yaml ("
        f"{','.join(DEFAULT_SCAFFOLDS)}); 8coh, 8z8v and gefurulimab also ship. Pass "
        "a single name to pin the framework -- but then every design shares it, which "
        "is the opposite of the diversity this protocol is for. REFUSED under "
        "protein-anything, which has no framework at all.",
    )
    parser.add_argument(
        "--scaffold-dir",
        type=str,
        default=None,
        help="nanobody-anything ONLY: directory bare --scaffolds names resolve "
        "against, INSIDE the task container. Defaults to $BOLTZGEN_SCAFFOLD_DIR or "
        f"{DEFAULT_SCAFFOLD_DIR}. Not checked for existence at submit time -- it does "
        "not exist where `sapia run` runs. REFUSED under protein-anything.",
    )
    parser.add_argument(
        "--num-designs",
        type=int,
        default=50,
        help="Designs to generate per target (-> --num_designs). Defaults to "
        "%(default)s, which is upstream's own smoke-test size; upstream's DEFAULT is "
        "10000 and its advice for real work is 10,000-60,000. Deliberately small "
        "here: per-design runtime on this image is UNMEASURED. Note issue #139 -- "
        "values above 10000 are rejected.",
    )
    parser.add_argument(
        "--diffusion-batch-size",
        type=int,
        default=None,
        help="-> --diffusion_batch_size. Omitted by default (upstream picks 1 below "
        "100 designs, else 10). **Every design in one batch shares the same sampled "
        "CDR length**, so a large batch relative to --num-designs under-samples "
        "lengths -- the same trap as rfd3's --length min-max.",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=1,
        help="-> --num_workers (DataLoader workers). Defaults to %(default)s, "
        "matching upstream's own default.",
    )
    parser.add_argument(
        "--spec",
        type=str,
        default=None,
        help="Path to a BoltzGen design YAML to use VERBATIM instead of writing one "
        "from the flags above. Its {expr} islands are still resolved per design and "
        "the result is written into the run_dir. Mutually exclusive with "
        "--target-structure/--target-chains/--binding/--not-binding/--scaffolds. "
        "Paths inside it are resolved by BoltzGen relative to the ORIGINAL yaml's "
        "directory, which is not where the resolved copy lives -- so use absolute "
        "paths inside a --spec file.",
    )
    parser.add_argument(
        "--extra-spec",
        type=str,
        default=None,
        help="Path to a YAML/JSON file whose top-level mapping is merged into EVERY "
        "written spec (e.g. extra `constraints:` or `templates:` blocks). String "
        "values may embed {expr}. A key this tool also sets ('entities') is an error.",
    )
    parser.add_argument(
        "--config-set",
        action="append",
        default=[],
        metavar="'STEP KEY=VALUE'",
        help="Extra BoltzGen step config override, emitted as `--config STEP "
        "KEY=VALUE`. Repeatable. Appended after the forced "
        "`--config analysis esmfold2_metrics=false`; setting that key to anything but "
        "false is REFUSED at submit time, because the analysis step would then crash "
        "after all four GPU steps have already run.",
    )
    parser.add_argument(
        "--extra-flags",
        type=str,
        default="",
        help="Extra `boltzgen run` tokens appended verbatim (space-split in the task "
        "script), e.g. '--use_kernels false'. Escape hatch for flags without a "
        "dedicated option. "
        f"REFUSED at submit time: {', '.join(BLOCKED_EXTRA_FLAGS)} -- they are this "
        "tool's contract and argparse is last-wins, so passing them here would "
        "silently change what the tool IS while every column, docstring and guard "
        "still described the other behaviour. --overwrite/--reuse additionally "
        "falsify the resubmit guard's own premise, and --diffusion_batch_size has a "
        "dedicated flag whose value the submit log reports.",
    )


# ---------------------------------------------------------------- residue specs


def parse_residue_spec(spec: str, flag: str) -> dict[str, list[int]]:
    """``'A91,A128,B91-95'`` -> ``{'A': [91, 128], 'B': [91, 92, 93, 94, 95]}``.

    Order within a chain is preserved; duplicates and malformed tokens raise. The
    epitope IS the question being asked, and a silently dropped or doubled residue
    changes what was designed against.
    """
    out: dict[str, list[int]] = {}
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if (m := _SINGLE.match(token)) is not None:
            chain, nums = m.group(1), [int(m.group(2))]
        elif (m := _RANGE.match(token)) is not None:
            chain, chain2 = m.group(1), m.group(3)
            start, end = int(m.group(2)), int(m.group(4))
            if chain2 and chain2 != chain:
                raise ValueError(
                    f"{flag}: range {token!r} spans two chains ({chain} and {chain2}); "
                    f"write one range per chain"
                )
            if end < start:
                raise ValueError(f"{flag}: range {token!r} ends before it starts")
            nums = list(range(start, end + 1))
        else:
            raise ValueError(
                f"{flag}: malformed residue token {token!r} (expected 'A96' or "
                f"'A96-99', chain letter then 1-based label_seq_id)"
            )
        seen = out.setdefault(chain, [])
        for n in nums:
            if n in seen:
                raise ValueError(f"{flag}: residue {chain}{n} listed twice")
            seen.append(n)
    if not out:
        raise ValueError(f"{flag}: empty residue spec")
    return out


def label_seq_residues(structure_path: Path) -> TargetFrames:
    """``{label_asym_id: ChainFrame}`` for the first model's polymers.

    **This replicates BoltzGen's own parser; it does not assume anything.** Both of
    BoltzGen's input paths index residues by position in the entity's DECLARED polymer
    sequence, 1-based, including positions with no coordinates:

    * mmCIF -- ``data/parse/mmcif.py:686-706`` iterates ``_entity_poly_seq`` and emits
      one residue per entry, present or not, at ``res_idx = label_seq - 1``;
    * PDB -- ``data/parse/pdb_parser.py:34-80`` runs ``setup_entities``, synthesizes
      ``full_sequence`` from the first non-empty subchain **when SEQRES is absent**,
      aligns each subchain into it with ``gemmi.align_sequence_to_polymer`` and writes
      ``label_seq``, then goes through the same mmCIF path.

    ``parse_range`` (``schema.py:649-667``) then turns the user's 1-based number into a
    0-based offset into *that* array. So the frame is the declared sequence's
    numbering, i.e. ``label_seq_id`` -- **not** ``auth_seq_id``, and **not** the ordinal
    of modelled residues.

    The distinction is not academic. Measured on 1TUP (p53, fetched 2026-10-08): chain
    A declares 219 residues and models 196; chains B and C model residues whose
    ``label_seq`` starts at 3 and 2, so their first modelled residue is position 3 and
    2, not 1. Counting modelled residues -- which an earlier version of this function
    did -- would have reported a different, perfectly plausible residue for every
    position after the first gap.

    **The CHAIN identifier is label_asym_id too, and that is not the letter a viewer
    shows.** ``mmcif.py:1055,1071`` parses chains from ``structure[0].subchains()`` and
    passes ``chain_id=subchain_id`` -- gemmi's subchain id IS the ``label_asym_id`` --
    which becomes ``structure.chains["name"]``, the array ``schema.py:2301`` resolves a
    YAML ``chain: id:`` against. Measured on 1TUP: author chains ``E,F,A,B`` are
    ``label_asym_id`` ``A,B,C,D``. A user who writes ``A`` meaning the protein would
    hit the DNA strand. This function therefore keys its result by the **author letter**
    -- what a caller reads off a viewer and writes in a flag -- and records the
    ``label_asym_id`` beside it, so ``resolve_chain_ids`` can TRANSLATE one into the
    other rather than guess. Matching is case-sensitive.

    ``ChainFrame.res_frame`` is ``declared_sequence`` when that chain's entity carried a
    polymer sequence (mmCIF ``_entity_poly_seq``, or PDB ``SEQRES``) and
    ``synthesized_from_model`` when it did not and the sequence had to be synthesized
    from coordinates, exactly as ``parse_pdb`` does. Per chain, not per file. Collected
    as ``boltzgen_binding_frame``.
    """
    try:
        st = gemmi.read_structure(str(structure_path))
    except (RuntimeError, ValueError, OSError) as e:
        raise ValueError(
            f"{structure_path.name}: gemmi cannot read this file ({type(e).__name__}: "
            f"{e}). The chain-id translation pivots on the DETECTED format, so an "
            f"unclassifiable target is refused rather than guessed at."
        ) from e
    fmt = str(st.input_format).rsplit(".", 1)[-1]
    if fmt not in ("Pdb", "Mmcif"):
        raise ValueError(
            f"{structure_path.name}: gemmi detected format {fmt!r}, which is neither "
            f"Pdb nor Mmcif. BoltzGen resolves chain ids differently for those two and "
            f"this tool will not guess for a third."
        )
    st.setup_entities()
    st.assign_label_seq_id()

    # pdb_parser.py:50-54: a file with no declared sequence gets one synthesized from
    # the first non-empty subchain of each polymer entity. Do the same, then force a
    # re-assignment so the alignment (and any per-copy gaps) is gemmi's, not ours.
    # Tracked PER ENTITY, not per file: a target can mix a deposited polymer (which
    # declares its sequence) with a modelled one (which does not), and reporting one
    # file-wide frame would mislabel the honest chains.
    synthesized_subchains: set[str] = set()
    for entity in st.entities:
        if entity.entity_type.name != "Polymer" or len(entity.full_sequence) >= 1:
            continue
        for model in st:
            for subchain_id in entity.subchains:
                sub = model.get_subchain(subchain_id)
                if len(sub):
                    entity.full_sequence = [res.name for res in sub]
                    synthesized_subchains.update(entity.subchains)
                    break
            if len(entity.full_sequence):
                break
    if synthesized_subchains:
        st.assign_label_seq_id(True)

    # Every written subchain's first character -- what patch_chain_names keys on for a
    # PDB input. Non-polymers included: a ligand subchain sharing the polymer's first
    # letter is exactly what makes the collision pool fire.
    first_chars: dict[str, list[str]] = {}
    nonpolymer: dict[str, list[str]] = {}
    for chain in st[0]:
        polymer_subs = {res.subchain for res in chain.get_polymer()}
        for res in chain:
            sub = res.subchain
            if sub and sub not in first_chars.setdefault(sub[:1], []):
                first_chars[sub[:1]].append(sub)
            if sub and sub not in polymer_subs:
                others = nonpolymer.setdefault(chain.name, [])
                if sub not in others:
                    others.append(sub)

    out: dict[str, ChainFrame] = {}
    for chain in st[0]:
        polymer = chain.get_polymer()
        if len(polymer) == 0:
            continue
        # A gemmi polymer span can cover more than one subchain. The whole translation
        # assumes one label_asym_id per author chain, so say so instead of picking.
        subs = sorted({res.subchain for res in polymer})
        if len(subs) != 1 or not subs[0]:
            raise ValueError(
                f"{structure_path.name}: author chain {chain.name!r} covers "
                f"{len(subs)} polymer subchain(s) ({subs}). This tool translates one "
                f"author letter to one label_asym_id; with several there is no single "
                f"id to write, and substituting the author letter would put an "
                f"untranslated id in the spec while the table said label_asym_id. "
                f"Split the chains, or crop the target."
            )
        if chain.name in out:
            raise ValueError(
                f"{structure_path.name}: author chain name {chain.name!r} appears more "
                f"than once. An epitope on it would resolve against whichever came "
                f"last; refusing rather than picking."
            )
        # Keyed by the AUTHOR letter, because that is what a caller reads off a
        # viewer and writes in --target-chains/--binding. The label_asym_id is carried
        # alongside; resolve_chain_ids translates one into the other.
        label_asym = subs[0]
        entity = st.get_entity_of(polymer)
        raw_declared = [str(m) for m in entity.full_sequence] if entity else []

        # MICROHETEROGENEITY. gemmi collapses several _entity_poly_seq rows that share
        # a `num` into one comma-joined full_sequence entry, i.e. ONE slot. BoltzGen
        # does not: mmcif.py:988-999 appends one [num, mon_id] per ROW with no
        # de-duplication, so its residue array gains an EXTRA slot per alternate and
        # every later position is shifted. Its own parse_polymer then matches the
        # modelled residue to the first row and treats the second as absent
        # (mmcif.py:690-706), so even upstream's frame is self-inconsistent there.
        # Refuse rather than attest a frame that is off by one.
        alternates = [m for m in raw_declared if "," in m]
        if alternates:
            raise ValueError(
                f"{structure_path.name}: chain {label_asym} has microheterogeneity in "
                f"its declared sequence ({alternates[:3]}{'...' if len(alternates) > 3 else ''}): "
                f"two _entity_poly_seq rows share a residue number. BoltzGen appends "
                f"one slot per ROW while every reader collapses them, so every residue "
                f"after the first such site is at a DIFFERENT index in BoltzGen's frame "
                f"than in the file's. Clean the alternates out of the target first."
            )

        declared = [m.strip() for m in raw_declared]
        labels = [res.label_seq for res in polymer]

        if not declared or any(label is None for label in labels):
            # Nothing to index against but the coordinates. BoltzGen would be in the
            # same position, so the ordinal IS the frame here -- but say so, PER CHAIN:
            # a file can mix declared and undeclared polymers.
            out[chain.name] = ChainFrame(
                label_asym=label_asym,
                auth_asym=chain.name,
                residues={
                    i: Residue(res.name, True)
                    for i, res in enumerate(polymer, start=1)
                },
                res_frame=FRAME_SYNTHESIZED,
            )
            continue

        table = {i: Residue(name, False) for i, name in enumerate(declared, start=1)}
        for res in polymer:
            pos = int(res.label_seq)
            if pos not in table:
                raise ValueError(
                    f"{structure_path.name}: chain {label_asym} has a residue at "
                    f"label_seq_id {pos} but its entity declares only "
                    f"{len(declared)} residue(s); the file is internally inconsistent "
                    f"and its residue numbering cannot be trusted"
                )
            if table[pos].name != res.name:
                raise ValueError(
                    f"{structure_path.name}: chain {label_asym} position {pos} is "
                    f"{res.name} in the coordinates but {table[pos].name} in the "
                    f"declared sequence; the numbering frame is ambiguous"
                )
            table[pos] = Residue(res.name, True)
        out[chain.name] = ChainFrame(
            label_asym=label_asym,
            auth_asym=chain.name,
            residues=table,
            res_frame=(
                FRAME_SYNTHESIZED
                if label_asym in synthesized_subchains
                else FRAME_DECLARED
            ),
        )

    return TargetFrames(
        chains=out, fmt=fmt, first_chars=first_chars, nonpolymer=nonpolymer
    )


def resolve_chain_ids(
    wanted: set[str], frames: TargetFrames, structure_path: Path
) -> tuple[dict[str, str], str]:
    """``({author letter: the id BoltzGen resolves}, frame)``.

    **This translation is the whole point, and it is documented nowhere upstream.**
    A caller reads chain letters off a viewer, which shows ``auth_asym_id``. What
    BoltzGen matches a YAML ``chain: id:`` against (``schema.py:2299-2303``, against
    ``ParsedChain.name``) depends on the input format -- which is taken from gemmi's
    DETECTION, never from the filename:

    * **mmCIF** -- ``parse_mmcif`` builds chains from gemmi **subchains**
      (``mmcif.py:1051-1053,1073,1112``), so the id is ``label_asym_id``. Measured on
      1TUP: author chains ``E,F,A,B`` are label chains ``A,B,C,D``. A caller asking for
      "A" meaning the p53 protein would otherwise hit the DNA strand -- silently, with
      a perfectly plausible nucleotide name attested back. **The common case is milder
      and just as wrong: an author chain split into a polymer subchain plus ligand
      subchains shifts every later chain's label.**
    * **PDB** -- ``parse_pdb`` round-trips through ``make_mmcif_block()`` (which renames
      ``A`` to ``Axp``) and then ``patch_chain_names`` truncates each name to its FIRST
      CHARACTER with a collision pool (``pdb_parser.py:102,107-126``), i.e. back to the
      author letter.

    Matching is **case-sensitive**: a file carrying both ``A`` and ``a`` is legal, and
    quietly folding one onto the other would attest a residue from a chain the caller
    did not name. Refusal is otherwise reserved for the genuinely unpredictable case --
    a PDB where two WRITTEN subchains share a first character, where
    ``patch_chain_names`` re-letters one from a reversed-alphabet pool.
    """
    is_pdb = frames.fmt == "Pdb"
    frame = CHAIN_FRAME_AUTHOR if is_pdb else CHAIN_FRAME_LABEL

    mapping: dict[str, str] = {}
    missing: list[str] = []
    for letter in sorted(wanted):
        chain = frames.chains.get(letter)
        if chain is None:
            missing.append(letter)
            continue
        mapping[letter] = letter if is_pdb else chain.label_asym

    if missing:
        folded = sorted(
            k for k in frames.chains if k.lower() in {m.lower() for m in missing}
        )
        raise ValueError(
            f"{structure_path.name}: chain(s) {missing} are not author chains of this "
            f"file (it has {sorted(frames.chains) or 'none'})."
            + (f" Matching is case-sensitive; did you mean {folded}?" if folded else "")
            + f" Chain letters are the AUTHOR letters a viewer shows; this tool "
            f"translates them to what BoltzGen resolves ({frame})."
        )

    if is_pdb:
        # patch_chain_names keys on the first character of every WRITTEN chain name --
        # polymer and non-polymer alike -- and re-letters a collision from a pool. So
        # with a ligand subchain beside its polymer the id BoltzGen ends up with is not
        # predictable from here, and it may be the POLYMER that gets re-lettered.
        clashing = sorted(
            letter
            for letter in mapping
            if len(frames.first_chars.get(frames.chains[letter].label_asym[:1], [])) > 1
        )
        if clashing:
            raise ValueError(
                f"{structure_path.name}: chain(s) {clashing} share a first character "
                f"with another written subchain ({frames.first_chars}). For a .pdb "
                f"input BoltzGen truncates chain names to one character and re-letters "
                f"collisions from a pool -- and it may be the polymer that is "
                f"re-lettered, putting the epitope on a chain nobody named. Convert the "
                f"target to mmCIF, or crop it to the chains you want."
            )
    return mapping, frame


def resolve_binding_resnames(
    wanted: dict[str, list[int]],
    frames: TargetFrames,
    structure_path: Path,
    flag: str,
) -> tuple[list[str], str]:
    """``(['A91=LYS', 'A128=GLU', ...], res_frame)`` -- the residue at each listed
    position, in the order listed, plus the frame it was read in.

    Labels stay in the caller's own AUTHOR letters: the attestation is there to be
    eyeballed against what they asked for, and ``boltzgen_chain_id_map`` carries the
    translation separately.

    Raises when a chain is absent, when a position is outside the declared sequence,
    or when a position exists but is **unmodelled**. All three are the same mistake:
    the numbering frame is wrong, and every design made against it would be against a
    site nobody asked for. An unmodelled position raises rather than warns because
    BoltzGen accepts it (the residue is parsed with ``is_present=False``) and a residue
    with no coordinates cannot be an epitope -- it would steer nothing while the table
    asserted an epitope.

    Every message names the FRAME, because the frame decides what the user has to do
    about it: renumber against the declared sequence, or against the modelled ordinals.
    """
    labels: list[str] = []
    seen_frames: set[str] = set()
    for chain, nums in wanted.items():
        if chain not in frames.chains:
            raise ValueError(
                f"{flag}: chain {chain!r} has no polymer in {structure_path.name} "
                f"(its author chains are {sorted(frames.chains) or 'none'}). Chain "
                f"letters here are the AUTHOR letters a viewer shows, matched "
                f"case-sensitively; the tool translates them for BoltzGen."
            )
        frame = frames.chains[chain]
        by_seq = frame.residues
        seen_frames.add(frame.res_frame)
        where = (
            f"chain {chain} of {structure_path.name} (frame {frame.res_frame}, "
            f"positions 1..{len(by_seq)})"
        )
        missing = [n for n in nums if n not in by_seq]
        if missing:
            raise ValueError(
                f"{flag}: position(s) {missing} absent from {where}. Positions are "
                f"1-based indices into the "
                + (
                    "DECLARED polymer sequence (mmCIF _entity_poly_seq / PDB SEQRES)"
                    if frame.res_frame == FRAME_DECLARED
                    else "sequence synthesized from the coordinates, because this file "
                    "declares none"
                )
                + " -- NOT auth_seq_id. Note a HET or modified residue that the parser "
                f"excludes from the polymer is not in this frame at all and reports the "
                f"same way."
            )
        unmodelled = [n for n in nums if not by_seq[n].modelled]
        if unmodelled:
            raise ValueError(
                f"{flag}: position(s) {unmodelled} on {where} are declared but have NO "
                f"COORDINATES (a disordered loop, or a terminus the depositor could not "
                f"resolve). BoltzGen would accept them and steer nothing. Pick modelled "
                f"residues, or crop the target first."
            )
        labels.extend(f"{chain}{n}={by_seq[n].name}" for n in nums)
    return labels, "+".join(sorted(seen_frames)) if seen_frames else ""


def _ranges(nums: list[int]) -> str:
    """``[26,27,28,51,52]`` -> ``'26..28,51..52'`` (BoltzGen's res_index syntax)."""
    out: list[str] = []
    for n in sorted(set(nums)):
        if out and n == int(out[-1].split("..")[-1]) + 1:
            head = out[-1].split("..")[0]
            out[-1] = f"{head}..{n}"
        else:
            out.append(str(n))
    return ",".join(out)


# ---------------------------------------------------------------- spec building


def load_yaml_mapping(path_str: str | None, flag: str) -> dict[str, Any]:
    """Parse a YAML/JSON file into a mapping; ``{}`` when unset or empty."""
    if not path_str:
        return {}
    path = Path(path_str)
    if not path.is_file():
        raise FileNotFoundError(f"{flag} file not found: {path_str}")
    data = yaml.safe_load(path.read_text())
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"{flag} must be a mapping, got {type(data).__name__}")
    return data


def _tree_has_placeholder(obj: Any) -> bool:
    """True if any string key/value anywhere in ``obj`` embeds a ``{expr}`` island."""
    if isinstance(obj, dict):
        return any(
            _tree_has_placeholder(k) or _tree_has_placeholder(v) for k, v in obj.items()
        )
    if isinstance(obj, list):
        return any(_tree_has_placeholder(v) for v in obj)
    if isinstance(obj, str):
        return bool(_HAS_PLACEHOLDER.search(obj))
    return False


def _resolve_tree(obj: Any, lookup: LookupFn, name: str) -> Any:
    """Recursively resolve ``{expr}`` islands; non-string scalars pass through."""
    if isinstance(obj, dict):
        return {
            resolve_template(str(k), lookup, name): _resolve_tree(v, lookup, name)
            for k, v in obj.items()
        }
    if isinstance(obj, list):
        return [_resolve_tree(v, lookup, name) for v in obj]
    if isinstance(obj, str):
        return resolve_template(obj, lookup, name)
    return obj


def resolved_scaffolds(args: BoltzGenArgs) -> tuple[str, str]:
    """``(--scaffolds, --scaffold-dir)`` with their defaults applied.

    Both flags default to ``None`` in the parser rather than to their values, so that
    ``_reject_protocol_mismatched_flags`` can tell "passed" from "not passed". The
    defaults are applied here instead.
    """
    scaffolds = args.scaffolds or ",".join(DEFAULT_SCAFFOLDS)
    scaffold_dir = args.scaffold_dir or os.getenv(
        "BOLTZGEN_SCAFFOLD_DIR", DEFAULT_SCAFFOLD_DIR
    )
    return scaffolds, scaffold_dir


def parse_binder_length(raw: str) -> str:
    """``--binder-length`` -> the ``sequence:`` value a ``protein:`` entity takes.

    ``'120'`` -> ``'120'``; ``'80-225'`` and ``'80..225'`` -> ``'80..225'``.

    **The hyphen form is translated, not passed through, and that is the whole point
    of this function.** BoltzGen tokenises a ``sequence:`` with
    ``re.findall(r"\\d+\\.\\.\\d+|\\d+|[a-zA-Z]", part)`` (``schema.py:947``), so the
    string ``'80-225'`` is read as TWO separate design blocks of 80 and 225 residues
    and silently produces a **305-residue** binder. A caller writing the workspace's
    own ``MIN-MAX`` convention (rfd3's ``--length`` uses it) would get a plausible run
    that designs the wrong thing. So the hyphen is handled explicitly and anything
    richer than a number or a range is refused rather than forwarded.

    The richer syntax BoltzGen supports -- interleaved fixed residues
    (``3..5C6C3``), comma-separated blocks -- is deliberately NOT exposed. It
    describes a constrained peptide, not a de-novo minibinder, and a tool that
    accepted it could no longer say what ``boltzgen_binder_len`` means. Use ``--spec``
    for those.
    """
    raw = raw.strip()
    if not raw:
        raise ValueError("--binder-length is empty; pass 'N' or 'MIN-MAX'")

    single = _LEN_SINGLE.match(raw)
    if single:
        n = int(single.group(1))
        if n < 1:
            raise ValueError(f"--binder-length {raw!r}: length must be >= 1")
        return str(n)

    span = _LEN_RANGE.match(raw)
    if not span:
        raise ValueError(
            f"--binder-length {raw!r} is malformed. Accepted: 'N' (exactly N "
            f"residues), 'MIN-MAX', or BoltzGen's own 'MIN..MAX'. Nothing else is "
            f"forwarded: BoltzGen tokenises a sequence spec by regex, so a form this "
            f"function does not recognise would be re-read as several design blocks "
            f"and silently design a chain of a different length. For interleaved "
            f"fixed residues or several blocks, write the YAML yourself and use "
            f"--spec."
        )
    lo, hi = int(span.group(1)), int(span.group(2))
    if lo < 1:
        raise ValueError(f"--binder-length {raw!r}: minimum must be >= 1")
    if hi < lo:
        raise ValueError(
            f"--binder-length {raw!r}: maximum {hi} is below minimum {lo}"
        )
    return f"{lo}..{hi}"


def pick_binder_chain(requested: str | None, taken: set[str]) -> str:
    """The chain id for the de-novo binder entity.

    ``taken`` is every chain id anything in the target FILE could contribute, in every
    frame it could be contributed under -- the ``label_asym_id`` of each polymer chain
    (what an mmCIF resolves to), the author letters (what a PDB resolves to), the
    resolved ids for the chains this run names, **and every non-polymer subchain id**
    (ligands, metals, glycans), plus the first character of each, which is what
    ``patch_chain_names`` keys on for a PDB. It is over-inclusive on purpose; the
    caller assembles it.

    A collision is not benign: ``Structure.concatenate(..., return_renaming=True)``
    (``schema.py:1633-1638``) renames one of the two and only logs it, so the chain
    labels in the output would not be the ones this tool wrote into the spec, and the
    attested epitope would name a chain that no longer exists under that id.

    A request that collides is REFUSED rather than quietly re-pointed -- a caller who
    named a chain wants that chain.
    """
    if requested is not None:
        requested = requested.strip()
        if len(requested) != 1 or not requested.isalpha():
            raise ValueError(
                f"--binder-chain {requested!r} must be a single letter"
            )
        if requested in taken:
            raise ValueError(
                f"--binder-chain {requested!r} is already used by the target "
                f"({sorted(taken)}). BoltzGen would rename one of the two and only "
                f"log it, so the chain labels in the output would not be the ones "
                f"written here. Pick a free letter, or omit --binder-chain and let "
                f"the tool choose."
            )
        return requested
    for letter in _BINDER_CHAIN_POOL:
        if letter not in taken:
            return letter
    raise ValueError(
        f"the target uses all 26 chain letters ({sorted(taken)}); there is no free id "
        f"for the de-novo binder. Crop the target with --target-chains."
    )


def scaffold_paths(scaffolds: str, scaffold_dir: str) -> list[str]:
    """``'7eow,gontivimab'`` -> the YAML paths BoltzGen draws from, in order.

    A token containing '/' is a path verbatim; anything else is a bare name resolved
    against ``scaffold_dir`` with a ``.yaml`` suffix. Never checked for existence:
    these paths live inside the task container, not where `sapia run` executes.
    """
    names = [s.strip() for s in scaffolds.split(",") if s.strip()]
    if not names:
        raise ValueError("--scaffolds is empty; pass at least one scaffold")
    out: list[str] = []
    for name in names:
        if "/" in name:
            out.append(str(Path(name)))
        else:
            stem = name[:-5] if name.endswith(".yaml") else name
            out.append(f"{scaffold_dir.rstrip('/')}/{stem}.yaml")
    return out


def build_entities(
    target_path: Path,
    args: BoltzGenArgs,
    binding: dict[str, list[int]] | None,
    not_binding: dict[str, list[int]] | None,
    chain_ids: dict[str, str],
    binder_chain: str | None = None,
    binder_sequence: str | None = None,
) -> list[dict[str, Any]]:
    """The spec's ``entities:`` list: the target file entity, then the BINDER entity.

    The binder entity is what the two protocols disagree about, and it is the only
    thing they disagree about in the YAML:

    * ``nanobody-anything`` -- a ``file:`` whose ``path:`` is a LIST of scaffold
      YAMLs. That list is what makes BoltzGen pick a framework at random per design
      (``data/parse/schema.py:1947-1949``). Mirrors
      ``example/nanobody_against_penguinpox``.
    * ``protein-anything`` -- a ``protein:`` with an ``id:`` and a ``sequence:`` that
      is a length or a length range. Mirrors ``example/vanilla_protein/1g13prot.yaml``.

    The target entity is byte-identical either way, which is why the epitope
    attestation, the chain-id translation and the residue frame are untouched by the
    protocol.
    """
    target: dict[str, Any] = {"path": str(target_path)}

    if args.target_chains.strip().lower() != "all":
        # The caller's AUTHOR letters, translated to whatever BoltzGen resolves for
        # this input format. See resolve_chain_ids -- on an mmCIF these are usually
        # different letters, and getting it wrong relocates the whole target.
        letters = expand_chain_spec(args.target_chains)
        target["include"] = [
            {"chain": {"id": chain_ids[letter]}} for letter in letters
        ]

    blocks: dict[str, dict[str, Any]] = {}
    for spec, key in ((binding, "binding"), (not_binding, "not_binding")):
        for chain, nums in (spec or {}).items():
            resolved = chain_ids[chain]
            blocks.setdefault(resolved, {"id": resolved})[key] = _ranges(nums)
    if blocks:
        target["binding_types"] = [
            {"chain": blocks[c]} for c in sorted(blocks)
        ]

    if args.protocol == PROTEIN_PROTOCOL:
        assert binder_chain is not None and binder_sequence is not None
        binder: dict[str, Any] = {
            "protein": {"id": binder_chain, "sequence": binder_sequence}
        }
    else:
        scaffolds, scaffold_dir = resolved_scaffolds(args)
        binder = {"file": {"path": scaffold_paths(scaffolds, scaffold_dir)}}

    # Target FIRST. Entity order decides which side loses a chain-id collision
    # (schema.py:1633-1638 renames the LATER one), and pick_binder_chain has already
    # made sure there is no collision -- but if one ever slipped through, renaming the
    # binder is recoverable (the collector finds it from the design mask) while
    # renaming a target chain would move the attested epitope.
    return [{"file": target}, binder]


class SpecResult(NamedTuple):
    """One target's spec plus everything the table has to be told about it."""

    spec: dict[str, Any]
    binding_resnames: list[str]  # the epitope, in the caller's AUTHOR letters
    not_binding_resnames: list[str]  # the ANTI-epitope -- a different claim
    res_frame: str
    chain_frame: str
    chain_id_map: str  # 'A>C,B>D' -- author letter > the id BoltzGen resolves
    dropped_nonpolymer: list[str]  # subchains an explicit --target-chains leaves out
    # protein-anything only; "" under nanobody-anything. Recorded in .submit.json as
    # provenance. NOT collected as a column: the collector reads the binder chain back
    # off the written structure's design mask, which is a measurement rather than a
    # request, and two columns claiming the same thing invite the reader to trust the
    # weaker one.
    binder_chain: str
    binder_sequence: str


def build_spec(
    name: str,
    target_path: Path,
    args: BoltzGenArgs,
    lookup: LookupFn,
    extra_fields: dict[str, Any],
) -> SpecResult:
    """One target's BoltzGen design spec and its epitope provenance.

    ``binding`` and ``not_binding`` resnames are kept APART. They are opposite claims
    about the same target in the same ``A91=LYS`` notation, so concatenating them would
    let a reader checking "where was this steered TO" read anti-epitope residues as
    epitope residues -- with nothing in the string to tell them apart.

    Raises ``ValueError`` on a per-row problem (an unresolvable ``{expr}``, an epitope
    position absent from or unmodelled in this target, an epitope chain that
    ``--target-chains`` excludes, a target whose chain conventions disagree) and the
    run-wide ``SpecConfigError`` on an --extra-spec key collision.
    """
    frames = label_seq_residues(target_path)
    binding = not_binding = None
    found: dict[str, list[str]] = {"--binding": [], "--not-binding": []}
    res_frames: set[str] = set()

    # The chains the spec will actually include. An epitope on a chain outside this set
    # is refused: BoltzGen's binding_types lookup would either raise ("Specified chain
    # id X not in file", schema.py:2301) or -- once `include` has cropped the file --
    # steer nothing at all while the table asserted an epitope.
    included = (
        None
        if args.target_chains.strip().lower() == "all"
        else set(expand_chain_spec(args.target_chains))
    )

    specs = {
        flag: parse_residue_spec(resolve_template(raw, lookup, name), flag)
        for flag, raw in (
            ("--binding", args.binding),
            ("--not-binding", args.not_binding),
        )
        if raw is not None
    }

    # Every chain letter this run will put in the YAML, checked for the chain-identifier
    # trap BEFORE any residue is resolved against it.
    named = set(included or ()) | {c for parsed in specs.values() for c in parsed}
    if not named:
        # `--target-chains all` with no epitope: BoltzGen takes every chain and no id
        # is ever written, so there is nothing to translate and nothing to get wrong.
        named = set(frames.chains)
    chain_ids, chain_frame = resolve_chain_ids(named, frames, target_path)

    # An explicit --target-chains writes `include:` with one entry per POLYMER subchain,
    # so any cofactor, metal or glycan on those author chains is left out -- while the
    # default `all` keeps everything. The two design against materially different
    # targets; name what is being dropped rather than letting it be noticed later.
    dropped = (
        sorted(
            sub
            for letter in (included or ())
            for sub in frames.nonpolymer.get(letter, [])
        )
        if included is not None
        else []
    )

    for flag, parsed in specs.items():
        if included is not None:
            outside = sorted(set(parsed) - included)
            if outside:
                raise ValueError(
                    f"{flag}: chain(s) {outside} are not in --target-chains "
                    f"({args.target_chains} -> {sorted(included)}), so they are not in "
                    f"the spec this run builds. Add them to --target-chains, or move "
                    f"the epitope onto a chain that is included."
                )
        labels, res_frame = resolve_binding_resnames(parsed, frames, target_path, flag)
        found[flag] = labels
        if res_frame:
            res_frames.add(res_frame)
        if flag == "--binding":
            binding = parsed
        else:
            not_binding = parsed

    # The de-novo binder entity, for protein-anything only.
    binder_chain = binder_sequence = ""
    if args.protocol == PROTEIN_PROTOCOL:
        assert args.binder_length is not None  # guarded in the manifest builder
        binder_sequence = parse_binder_length(
            resolve_template(args.binder_length, lookup, name)
        )
        # Every chain id ANYTHING IN THE TARGET FILE could contribute, in every frame
        # it could be contributed under. Deliberately over-inclusive: the pool has 26
        # entries and the cost of a collision is not a clash but a SILENT RENAME
        # (Structure.concatenate, schema.py:1633-1638, logs and continues), which is
        # the case pick_binder_chain exists to refuse.
        #
        # Four sources, and the last two are the ones an earlier version missed:
        #   * the label_asym_id of each POLYMER chain -- the id an mmCIF target
        #     resolves to, and the only one the translation covers;
        #   * the author letters, which is what a PDB target resolves to;
        #   * the resolved ids for the chains this run actually names;
        #   * **every NON-POLYMER subchain id.** frames.chains holds polymer chains
        #     only, so a target carrying ligands, metals or glycans contributes
        #     subchain ids that were invisible here. On an mmCIF with many hetero
        #     subchains those run well into the alphabet, and BoltzGen parses them as
        #     chains like any other (mmcif.py builds NONPOLYMER chains; source_context
        #     has no type filter). first_chars is keyed by first character over ALL
        #     written subchains, polymer or not, so its VALUES are the full id list.
        taken = {frame.label_asym for frame in frames.chains.values()}
        taken |= set(frames.chains)  # author letters, for a PDB input
        taken |= set(chain_ids.values())
        taken |= {sub for subs in frames.first_chars.values() for sub in subs}
        taken |= {sub for subs in frames.nonpolymer.values() for sub in subs}
        # A PDB's ids are patch_chain_names' FIRST CHARACTERS, so a two-character
        # subchain id like 'Ba' makes 'B' unavailable too.
        taken |= {sub[:1] for sub in list(taken) if sub}
        binder_chain = pick_binder_chain(args.binder_chain, taken)

    spec: dict[str, Any] = {
        "entities": build_entities(
            target_path,
            args,
            binding,
            not_binding,
            chain_ids,
            binder_chain=binder_chain or None,
            binder_sequence=binder_sequence or None,
        )
    }

    resolved_extra = _resolve_tree(extra_fields, lookup, name)
    collisions = sorted(set(spec) & set(resolved_extra))
    if collisions:
        raise SpecConfigError(
            f"--extra-spec sets key(s) {collisions} that this tool also writes; "
            f"use --spec to supply the whole YAML instead"
        )
    spec.update(resolved_extra)
    return SpecResult(
        spec=spec,
        binding_resnames=found["--binding"],
        not_binding_resnames=found["--not-binding"],
        res_frame="+".join(sorted(res_frames)),
        chain_frame=chain_frame,
        chain_id_map=",".join(f"{a}>{b}" for a, b in sorted(chain_ids.items())),
        dropped_nonpolymer=dropped,
        binder_chain=binder_chain,
        binder_sequence=binder_sequence,
    )


def _reject_protocol_mismatched_flags(args: BoltzGenArgs) -> None:
    """Refuse a flag that has NO MEANING under the chosen protocol.

    Silently ignoring a flag the caller passed is the failure mode this workbench
    keeps getting bitten by: the run proceeds, the numbers look plausible, and nothing
    in the table says the request was dropped. Both directions are refused:

    * ``--scaffolds`` / ``--scaffold-dir`` under ``protein-anything`` -- there is no
      framework to draw from, so the flag would select nothing. A caller who passes it
      has misunderstood which experiment they are running, which is worth an error.
    * ``--binder-length`` / ``--binder-chain`` under ``nanobody-anything`` -- the VHH
      framework fixes the chain and its length, and the CDR lengths come from the
      scaffold's own sampling. There is no ``sequence:`` for the flag to reach, so it
      would be written nowhere.
    * ``--binder-length`` MISSING under ``protein-anything`` -- there is no sensible
      default length for a de-novo binder, and BoltzGen has no default either: the
      ``protein:`` entity's ``sequence:`` is mandatory.

    Checked for ALL of them, so one error names every mismatch rather than making the
    caller fix them one submit at a time.
    """
    protein = args.protocol == PROTEIN_PROTOCOL
    wrong = [
        flag
        for flag, value, only_under in (
            ("--scaffolds", args.scaffolds, NANOBODY_PROTOCOL),
            ("--scaffold-dir", args.scaffold_dir, NANOBODY_PROTOCOL),
            ("--binder-length", args.binder_length, PROTEIN_PROTOCOL),
            ("--binder-chain", args.binder_chain, PROTEIN_PROTOCOL),
        )
        if value is not None and only_under != args.protocol
    ]
    if wrong:
        raise SpecConfigError(
            f"--protocol {args.protocol} does not use {sorted(wrong)}. They would be "
            f"silently ignored, so they are refused instead: "
            f"--scaffolds/--scaffold-dir name a VHH framework and {PROTEIN_PROTOCOL} "
            f"has none; --binder-length/--binder-chain describe a de-novo chain and "
            f"{NANOBODY_PROTOCOL} designs CDRs on a fixed framework instead."
        )
    if protein and args.binder_length is None and args.spec is None:
        raise SpecConfigError(
            f"--protocol {PROTEIN_PROTOCOL} needs --binder-length: the de-novo binder "
            f"chain's `sequence:` is what sets its length and there is no default. "
            f"Pass 'N' or 'MIN-MAX' (and --diffusion-batch-size 1 if a range matters "
            f"-- one length is drawn per BATCH, not per design)."
        )


def _reject_spec_building_flags(args: BoltzGenArgs) -> None:
    """--spec owns the whole YAML, so the flags that write one must be absent."""
    conflicting = [
        flag
        for flag, value in (
            ("--target-structure", args.target_structure),
            ("--binding", args.binding),
            ("--not-binding", args.not_binding),
            ("--scaffolds", args.scaffolds),
            ("--scaffold-dir", args.scaffold_dir),
            ("--binder-length", args.binder_length),
            ("--binder-chain", args.binder_chain),
        )
        if value is not None
    ]
    if args.target_chains.strip().lower() != "all":
        conflicting.append("--target-chains")
    if conflicting:
        raise SpecConfigError(
            f"--spec supplies the whole design YAML, so {sorted(conflicting)} would "
            f"be silently ignored. Drop them, or drop --spec and let this tool write "
            f"the YAML."
        )


# ---------------------------------------------------------------- CLI assembly


def check_config_set(raw: str) -> str:
    """Validate one ``--config-set`` entry and return it.

    Refuses anything that would reinstate ``esmfold2_metrics=true``: the analysis step
    reads ``esmfold2_scores/<id>.input.json`` unconditionally when it is on
    (``analyze.py:638-652``), and this tool never runs ``esmfold2_scoring``, so turning
    it back on is a guaranteed crash four GPU steps in.
    """
    tokens = raw.split()
    if len(tokens) < 2 or "=" not in " ".join(tokens[1:]):
        raise ValueError(
            f"--config-set {raw!r} must be '<step> <key>=<value>', e.g. "
            f"'analysis liability_modality=antibody'"
        )
    for kv in tokens[1:]:
        key, _, value = kv.partition("=")
        if key.strip() == PROTECTED_CONFIG_KEY and value.strip().lower() != "false":
            raise ValueError(
                f"--config-set {raw!r} would set {PROTECTED_CONFIG_KEY}={value!r}. "
                f"This tool never runs the esmfold2_scoring step, and the analysis "
                f"step reads esmfold2_scores/<id>.input.json unconditionally when "
                f"{PROTECTED_CONFIG_KEY} is on -- the run would crash after the GPU "
                f"work. It is forced to false and must stay there."
            )
    return raw


def check_extra_flags(extra: str) -> str:
    """Validate ``--extra-flags`` and return it stripped.

    The tokens are appended AFTER the flags this tool owns and argparse is last-wins,
    so an unguarded ``--extra-flags`` can silently redefine the tool: ``--protocol
    protein-anything`` would change what is being designed while the docstrings, the
    skill, the ``boltzgen_`` column prefix and the single-designed-chain guard all
    still say nanobody, and the guard would still pass.
    """
    extra = extra.strip()
    if not extra:
        return ""
    offending = sorted(
        {
            flag
            for token in extra.split()
            for flag in BLOCKED_EXTRA_FLAGS
            if token == flag or token.startswith(flag + "=")
        }
    )
    if offending:
        raise ValueError(
            f"--extra-flags may not contain {offending}: these are set by the tool "
            f"itself and argparse is last-wins, so passing them here would silently "
            f"change what this tool IS while every column, docstring and guard still "
            f"said nanobody. Use --config-set for step config; for anything else, "
            f"the tool needs a dedicated flag."
        )
    return extra


def cli_flags(args: BoltzGenArgs) -> str:
    """The run-wide `boltzgen run` flags, identical for every target.

    ``--protocol``, ``--steps``, ``--output`` and ``--cache`` are in the task script,
    not here, and ``check_extra_flags`` refuses to let them back in: they are the
    tool's contract, not a per-run choice (see the module docstring).
    """
    flags: list[str] = [f"--num_designs {args.num_designs}"]
    if args.diffusion_batch_size is not None:
        flags.append(f"--diffusion_batch_size {args.diffusion_batch_size}")
    flags.append(f"--num_workers {args.num_workers}")
    for step, kv in FORCED_CONFIG:
        flags.append(f"--config {step} {kv}")
    for raw in args.config_set:
        flags.append(f"--config {check_config_set(raw)}")
    extra = check_extra_flags(args.extra_flags)
    if extra:
        flags.append(extra)
    return " ".join(flags)


def target_groups(ctx: ManifestCtx[BoltzGenArgs]) -> list[tuple[str, Path | None]]:
    """The (name, target_path) groups this run designs against.

    One per ready table row for a child run; a single group for a root run. With
    ``--spec`` the target lives inside the YAML, so the path is ``None``.
    """
    args = ctx.args

    if args.table is not None:
        if args.target_structure is not None:
            raise SpecConfigError(
                "--target-structure is only valid for a root run (no --table); with "
                "--table, targets come from the table's --input-column. Drop one."
            )
        groups: list[tuple[str, Path | None]] = []
        for name in ctx.ready.index:
            name = cast(str, name)
            if args.spec is not None:
                # The row parameterises the run (it decides readiness and supplies the
                # lineage {expr} resolves against); the target is the YAML's.
                groups.append((name, None))
                continue
            target_path = volume_path(str(ctx.ready.at[name, args.input_column]))
            if not target_path.exists():
                print(f"{name}: MISSING {target_path} (skipping)")
                continue
            groups.append((name, target_path))
        return groups

    # Root run: {expr} islands resolve up a table lineage this run does not have.
    for raw, flag in (
        (args.binding, "--binding"),
        (args.not_binding, "--not-binding"),
    ):
        if raw and _HAS_PLACEHOLDER.search(raw):
            raise SpecConfigError(
                f"{flag} contains a {{expr}} placeholder, but this is a root run (no "
                f"--table) with no table lineage to resolve it against. Use literal "
                f"values, or run with --table."
            )

    if args.spec is not None:
        return [(f"{Path(args.spec).stem}_bg", None)]
    if args.target_structure is None:
        raise SpecConfigError(
            "A root run (no --table) needs a target: pass --target-structure <file> "
            "or --spec <file>. With --table, targets come from --input-column."
        )
    target_path = volume_path(args.target_structure)
    if not target_path.exists():
        raise FileNotFoundError(f"--target-structure {target_path} does not exist.")
    return [(f"{target_path.stem}_bg", target_path)]


SUBMIT_RECORD = ".submit.json"
# Where BoltzGen writes its designs under a group's --output. Duplicated from
# collect_boltzgen.DESIGN_DIRNAME on purpose: the builder needs it to tell a group that
# has already produced outputs from a fresh one, and importing the collector into the
# run module would couple the two phases.
_DESIGN_DIRNAME = "intermediate_designs_inverse_folded"


def read_submit_record(group_dir: Path) -> tuple[dict | None, str]:
    """``(record, state)`` -- a group's own submit record, and why it is absent.

    ``state`` is ``"present"``, ``"absent"``, or a description of the damage. The caller
    MUST distinguish them: "no readable record" and "a record that disagrees" are the
    same situation for an output dir that already holds designs, and treating only the
    second as a problem is what let a deleted record re-enable a stale top-up.

    **Why this lives beside the group's outputs and not in the run-level sidecar.**
    The epitope attestation is a claim about the designs in ``runs/<group>/``, but it
    used to be stored in ``<out_dir>/.meta.json`` keyed by group NAME. Group names are
    deterministic (``<stem>_bg``, or the parent row name), so re-running a target into
    the same out_dir -- the ordinary thing to do -- merged a new submit's records over
    an old one's and left any group it did not resubmit untouched. Three consequences,
    all silent:

    * a group resubmitted WITHOUT ``--binding`` kept the previous epitope string, so
      unsteered designs collected as ``OK`` carrying a certified epitope;
    * ``prosapia``'s ``write_run_meta`` runs BEFORE the manifest builder
      (``base_run.py:345`` vs ``:377``) and swallows a corrupt read, so a damaged
      sidecar was replaced by a fresh object and every per-group record vanished before
      this module could notice;
    * the read-modify-write is unlocked (no ``flock`` anywhere in ``prosapia``), so two
      overlapping submits could drop one another's groups entirely.

    Per-group it dies with the outputs it describes, a resubmit replaces it, and there
    is no shared file to race on.
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


def claims_in_spec(spec: dict[str, Any]) -> list[str]:
    """Which steering claims a resolved spec dict asks for: ``binding`` / ``not_binding``.

    Derived from the spec itself, never from the flags, so the record always describes
    the YAML BoltzGen is handed -- including a ``--spec`` YAML the user wrote, whose
    ``binding_types`` the flags know nothing about. ``collect_boltzgen.claims_from_spec``
    is the same walk over the YAML on disk and cross-checks this; a test asserts the two
    agree on the same spec. **Keep the two in step.**
    """
    claims: set[str] = set()
    for entity in spec.get("entities") or []:
        block = (entity or {}).get("file") if isinstance(entity, dict) else None
        for item in (block or {}).get("binding_types") or []:
            chain = (item or {}).get("chain") if isinstance(item, dict) else None
            for key in ("binding", "not_binding"):
                if chain and key in chain:
                    claims.add(key)
    return sorted(claims)


def _spec_sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def build_boltzgen_manifest(ctx: ManifestCtx[BoltzGenArgs]) -> list[tuple[str, ...]]:
    args = ctx.args

    # Before anything else: a flag that means nothing under this protocol is an error,
    # never a silent no-op. Run-wide, so it fails the submit rather than every row.
    _reject_protocol_mismatched_flags(args)

    if args.spec is not None:
        _reject_spec_building_flags(args)
        verbatim = load_yaml_mapping(args.spec, "--spec")
        if not verbatim:
            raise SpecConfigError(f"--spec {args.spec} parsed to an empty mapping")
    else:
        verbatim = {}
    extra_fields = load_yaml_mapping(args.extra_spec, "--extra-spec")

    if args.table is None and (
        _tree_has_placeholder(extra_fields) or _tree_has_placeholder(verbatim)
    ):
        raise SpecConfigError(
            "--spec/--extra-spec contain a {expr} placeholder, but this is a root run "
            "(no --table) with no table lineage to resolve it against."
        )

    groups = target_groups(ctx)
    if not groups:
        return []

    specs_dir = ctx.out_dir / SPECS_DIRNAME
    specs_dir.mkdir(parents=True, exist_ok=True)
    runs_root = volume_path(ctx.out_dir) / RUNS_DIRNAME

    # The scaffold allow-list the collector validates `source_context` against. Under
    # --spec the scaffolds come from the user's YAML, not from --scaffolds, so there is
    # NO allow-list and the collector must say so rather than check the wrong one.
    # Under protein-anything there is no framework at all, so the list is likewise
    # empty -- but the collector distinguishes the two from `protocol`, not from the
    # empty list, so that "no allow-list to check" and "nothing to check" do not share
    # a verdict.
    # Recorded PER GROUP, because a later submit with a different --scaffolds would
    # otherwise re-point an earlier group's allow-list at a set it never used.
    if verbatim or args.protocol == PROTEIN_PROTOCOL:
        scaffolds: list[str] = []
    else:
        scaffolds = scaffold_paths(*resolved_scaffolds(args))

    steps = PROTOCOL_STEPS[args.protocol]
    flags = cli_flags(args)
    manifest_rows: list[tuple[str, ...]] = []
    submitted: list[str] = []
    refused: list[tuple[str, str]] = []

    for name, target_path in groups:
        group_dir = runs_root / name
        try:
            if verbatim:
                spec = _resolve_tree(verbatim, ctx.lookup, name)
                result = None
            else:
                assert target_path is not None  # target_groups guarantees it
                result = build_spec(name, target_path, args, ctx.lookup, extra_fields)
                spec = result.spec

            spec_text = yaml.safe_dump(spec, sort_keys=False)
            sha = _spec_sha256(spec_text)

            # An output dir that already holds designs may be topped up ONLY by a
            # submit that PROVES it is the same spec. BoltzGen is never given
            # --overwrite, so anything already there survives and would collect under
            # whatever attestation this run writes.
            #
            # The predicate is deliberately "require a readable match", not "forbid a
            # confirmed mismatch". A record that is absent, truncated, half-written or
            # deleted gives no evidence either way, and the earlier form -- gated on a
            # readable previous sha -- silently PASSED every one of those cases. The
            # worst of them was self-inflicted: the refusal message used to say "delete
            # {group_dir}", which reads as "delete the stale record", and doing that
            # turned the refusal off. Hence the advice below names the whole directory
            # and says the record alone is not enough.
            previous, state = read_submit_record(group_dir)
            if (group_dir / _DESIGN_DIRNAME).is_dir():
                recorded = (previous or {}).get("spec_sha256")
                # THE PROTOCOL IS PART OF WHAT A RESUBMIT MUST PROVE, and the spec
                # hash does NOT cover it.
                #
                # Under --spec the written YAML is the user's file verbatim, so it is
                # byte-identical whichever --protocol was passed. Without this check,
                #     sapia run boltzgen R --spec x.yaml
                #     sapia run boltzgen R --spec x.yaml --protocol protein-anything
                # passes the sha comparison, rewrites .submit.json, and stamps
                # boltzgen_protocol = protein-anything onto designs that actually drew
                # a random VHH framework -- while forcing their perfectly recoverable
                # scaffold provenance to NA. That is verbatim the failure
                # boltzgen_protocol was introduced to prevent, committed by the column
                # itself.
                #
                # Compared explicitly rather than folded into the hash on purpose: a
                # hash would also invalidate every record written before `protocol`
                # existed, and it could not NAME the mismatch in the refusal. `steps`
                # is checked beside it because it is what the task script actually
                # runs; today it is derived from the protocol, but a future step-list
                # change under an unchanged protocol must not slip through either.
                rec_protocol = (previous or {}).get("protocol")
                rec_steps = (previous or {}).get("steps")
                # Same doctrine as the sha: REQUIRE A READABLE MATCH, never merely
                # "forbid a confirmed mismatch". A record predating these keys gives
                # no evidence either way, so it is refused like any other record that
                # cannot prove what it describes.
                missing = (
                    "missing a spec_sha256"
                    if not recorded
                    else "missing a protocol"
                    if not rec_protocol
                    else "missing a steps list"
                    if not rec_steps
                    else None
                )
                if state != "present" or missing:
                    raise ValueError(
                        f"{group_dir} already holds designs but its "
                        f"{SUBMIT_RECORD} is {state if state != 'present' else missing}, "
                        f"so there is no evidence they were made from this run's spec "
                        f"and protocol. BoltzGen is never run with --overwrite and "
                        f"those designs would survive to collect under this run's "
                        f"epitope and protocol. Use -l/--dir-label for a fresh leaf, "
                        f"or delete the WHOLE directory {group_dir} -- deleting "
                        f"{SUBMIT_RECORD} on its own leaves the old designs in place "
                        f"and is exactly what this check exists to catch. (A record "
                        f"written before 2026-10-09 has no protocol/steps and lands "
                        f"here; that is intended -- it cannot attest what it never "
                        f"recorded.)"
                    )
                if recorded != sha:
                    raise ValueError(
                        f"{group_dir} already holds designs made from a DIFFERENT spec "
                        f"(recorded {recorded[:12]}, this run {sha[:12]}). BoltzGen is "
                        f"never run with --overwrite, so those designs would survive "
                        f"and collect under this run's epitope. Use -l/--dir-label for "
                        f"a fresh leaf, or delete the WHOLE directory {group_dir} -- "
                        f"deleting {SUBMIT_RECORD} on its own leaves the old designs "
                        f"in place and would turn this check off."
                    )
                # SHAPE FIRST, then compare. `list()` on a JSON number raises
                # TypeError, which the caller's `except (ValueError, OSError,
                # RuntimeError)` does not catch -- so a single hand-corrupted record
                # (`"steps": 5`) would kill the WHOLE submit with a traceback instead
                # of refusing one group. A string degrades differently and just as
                # badly: `list("design folding")` is a per-character comparison that
                # refuses with an incomprehensible message. Both are loud, neither is
                # graceful, and the record is a file a human can edit.
                #
                # Named rather than described: the refusal says what the record
                # actually contains, because "steps is malformed" sends the reader
                # back to the file to find out what is wrong with it.
                if not isinstance(rec_steps, list) or not all(
                    isinstance(s, str) for s in rec_steps
                ):
                    raise ValueError(
                        f"{group_dir}/{SUBMIT_RECORD} has a malformed 'steps': "
                        f"expected a list of strings, found {type(rec_steps).__name__} "
                        f"{rec_steps!r}. The record cannot attest which steps produced "
                        f"the designs already in that directory, so this run will not "
                        f"top them up. Use -l/--dir-label for a fresh leaf, or delete "
                        f"the WHOLE directory {group_dir}."
                    )
                if not isinstance(rec_protocol, str):
                    raise ValueError(
                        f"{group_dir}/{SUBMIT_RECORD} has a malformed 'protocol': "
                        f"expected a string, found {type(rec_protocol).__name__} "
                        f"{rec_protocol!r}. Use -l/--dir-label for a fresh leaf, or "
                        f"delete the WHOLE directory {group_dir}."
                    )
                if rec_protocol != args.protocol or list(rec_steps) != list(steps):
                    raise ValueError(
                        f"{group_dir} already holds designs made under a DIFFERENT "
                        f"PROTOCOL (recorded {rec_protocol!r} / {list(rec_steps)}, "
                        f"this run {args.protocol!r} / {list(steps)}). The spec hash "
                        f"matches, which under --spec it will for ANY protocol -- the "
                        f"YAML is yours and does not mention one. Topping up would "
                        f"stamp this run's protocol onto designs made by the other "
                        f"one, which is precisely what boltzgen_protocol exists to "
                        f"prevent. Use -l/--dir-label for a fresh leaf, or delete the "
                        f"WHOLE directory {group_dir}."
                    )
        except SpecConfigError:
            # A run-wide misconfiguration hits every row identically: fail fast
            # instead of silently skipping the entire table.
            raise
        except (ValueError, OSError, RuntimeError) as e:
            # One bad row (an unresolvable {expr}, an epitope absent from THIS target)
            # must not sink the whole array -- but a skipped group leaves NO status
            # anywhere, because a create tool cannot stamp a parent it never submitted.
            # So every skip is counted and reported in the summary below.
            print(f"{name}: REFUSED -- {e}")
            refused.append((name, str(e)))
            continue

        spec_yaml = specs_dir / f"{name}.yaml"
        spec_yaml.write_text(spec_text)

        record: dict[str, Any] = {
            "spec_sha256": sha,
            # WHICH PROTOCOL THESE DESIGNS ARE. Collected as boltzgen_protocol, from
            # here and not from a flag re-read at collect time: the flag describes
            # whatever the LAST submit asked for, while this file sits beside the
            # designs it describes and is replaced only when they are. A run_dir can
            # hold nanobody and minibinder groups side by side, and a per-run value
            # would label both with whichever ran last.
            "protocol": args.protocol,
            "steps": list(steps),
            "scaffolds": scaffolds,
            # Derived from the SPEC, not the flags: a --spec YAML carries its own
            # binding_types, and recording [] for it made every design error at collect
            # with a message blaming a stale record for a bookkeeping gap here.
            "claims": claims_in_spec(spec),
            "binding_resnames": ",".join(result.binding_resnames) if result else "",
            "not_binding_resnames": (
                ",".join(result.not_binding_resnames) if result else ""
            ),
            "res_frame": result.res_frame if result else "",
            "chain_frame": result.chain_frame if result else "",
            "chain_id_map": result.chain_id_map if result else "",
            # protein-anything provenance. Not collected -- the collector measures the
            # binder chain and its length off the written structure instead -- but
            # recorded so a disagreement between request and result is diagnosable.
            "binder_chain_requested": result.binder_chain if result else "",
            "binder_sequence_spec": result.binder_sequence if result else "",
        }
        group_dir.mkdir(parents=True, exist_ok=True)
        (group_dir / SUBMIT_RECORD).write_text(json.dumps(record, indent=2))

        # Fields 5 and 6 are new: the protocol and the step list used to be hard-coded
        # in boltzgen.sh, which cannot be right once there are two protocols. Fields
        # 1-4 are unchanged and the script still cuts them at the same indices, so a
        # nanobody submit produces exactly the same `boltzgen run` command line it did
        # before -- the .sh interpolates these two where the literals used to be.
        manifest_rows.append(
            (
                name,
                str(volume_path(spec_yaml)),
                str(group_dir),
                flags,
                args.protocol,
                " ".join(steps),
            )
        )
        submitted.append(name)

        if args.protocol == PROTEIN_PROTOCOL and result is not None:
            print(
                f"{name}: de-novo binder chain {result.binder_chain} "
                f"sequence {result.binder_sequence}"
                + (
                    "  [RANGE: one length per BATCH, not per design -- "
                    "--diffusion-batch-size 1 samples it per design]"
                    if ".." in result.binder_sequence
                    and (args.diffusion_batch_size or 0) != 1
                    else ""
                )
            )

        if result is not None:
            if result.binding_resnames:
                print(f"{name}: epitope [{result.res_frame}] "
                      f"{record['binding_resnames']}")
            if result.not_binding_resnames:
                print(f"{name}: ANTI-epitope [{result.res_frame}] "
                      f"{record['not_binding_resnames']}")
            if result.chain_frame == CHAIN_FRAME_LABEL and any(
                a != b
                for a, b in (
                    pair.split(">") for pair in result.chain_id_map.split(",") if pair
                )
            ):
                print(
                    f"{name}: chain ids TRANSLATED author>label_asym_id "
                    f"{result.chain_id_map} -- BoltzGen resolves label_asym_id for an "
                    f"mmCIF target, and they differ in this file"
                )
            if result.dropped_nonpolymer:
                print(
                    f"{name}: --target-chains DROPS non-polymer subchain(s) "
                    f"{result.dropped_nonpolymer} (cofactors, metals, glycans). The "
                    f"default --target-chains all would keep them."
                )

    total = len(groups)
    print(
        f"Submitting {len(submitted)} of {total} target(s)"
        + (f"; {len(refused)} REFUSED" if refused else "")
    )
    for name, reason in refused:
        print(f"  refused {name}: {reason}")
    if refused and not submitted:
        # A filter that fails shut reads exactly like a filter that found nothing. All
        # of them refused is a configuration error, not a result.
        raise SpecConfigError(
            f"all {total} target(s) were refused; nothing would be submitted. First "
            f"reason: {refused[0][1]}"
        )

    # Run-level facts only. `protocol` is here for a human reading the sidecar; the
    # COLUMN is read from each group's own .submit.json, because .meta.json is merged
    # across submits and a run_dir can hold both protocols.
    ctx.write_meta(
        runs_root=str(runs_root),
        protocol=args.protocol,
        steps=list(steps),
        num_designs=args.num_designs,
    )
    if args.table is None:
        # A root run has no parent table for collect to iterate; record the group names
        # so `sapia collect` can find their outputs. REBUILT FROM DISK rather than from
        # the old sidecar, because the framework rewrites .meta.json before this
        # function runs and a corrupt one is silently replaced -- every group that ever
        # left a submit record is still here.
        on_disk = {
            d.name
            for d in runs_root.glob("*")
            if d.is_dir() and (d / SUBMIT_RECORD).is_file()
        }
        ctx.write_meta(root_designs=sorted(on_disk | set(submitted)))

    return manifest_rows
