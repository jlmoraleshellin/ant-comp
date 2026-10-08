#!/usr/bin/env python3
"""
Submit BoltzGen nanobody-binder design runs, one array task per target.

PREMISE -- read this before using any column this tool writes
------------------------------------------------------------
BoltzGen is campaign-driven, not design-driven: one ``boltzgen run <spec>.yaml``
process generates N designs against one target and only then are there rows. So the
unit of work is **one BoltzGen run per target**, the row count is known at collect
time, and this is a ``create`` tool.

This tool wraps exactly ONE of BoltzGen's protocols: ``nanobody-anything`` --
grafting new CDRs onto a fixed VHH framework and docking it onto a (possibly
multi-chain) protein target. Its premise, and the conditions under which its numbers
are meaningless:

* **The binder is a nanobody made of ONE designed chain.** Framework sequence AND
  framework structure come from a scaffold ``.cif``; only CDR1/2/3 are designed, in
  sequence, structure and length. Nothing here designs a de-novo backbone -- for
  that use ``rfdiffusion3``/``bindcraft2``. If a spec is supplied (``--spec``) whose
  designed region is not a single chain, the collector refuses the design rather
  than guessing which chain the sequence belongs to.
* **Residue indices in ``--binding``/``--not-binding`` are 1-based ``label_seq_id``,
  NOT ``auth_seq_id``** (upstream README, "IMPORTANT"). A crystallographic target
  numbered 18..155 is 1..138 to BoltzGen. This builder resolves each requested
  position against the target structure at SUBMIT time and records the residue NAME
  it found, so a silently wrong epitope shows up as the wrong amino acids in
  ``boltzgen_binding_resnames`` rather than as a plausible but misplaced design. A
  position that does not exist is an error, never a dropped hotspot.
* **The scaffold is drawn at RANDOM per design** (``data/parse/schema.py:1947-1949``
  does ``path = random.choice(path)``). That is wanted here -- it is the diversity
  axis -- and there is no seed to pin it with (BoltzGen has no ``--seed`` flag at
  all). Which scaffold a given design actually used is recovered per design by the
  collector into ``boltzgen_scaffold``. **Designs are therefore not reproducible,
  deliberately.** Do not "fix" this.
* **Target geometry is not guaranteed.** Open issue #226: on a multi-chain target a
  second protomer can MOVE during design, giving a complex whose binding geometry is
  wrong. Unanswered upstream. ``boltzgen_bb_rmsd_target`` is the column that catches
  it and it must be read before any interface number.

ESMFold2 is EXCLUDED, by construction and not by accident
---------------------------------------------------------
Every run is ``--steps design inverse_folding folding analysis`` plus
``--config analysis esmfold2_metrics=false``. Three separate facts force that exact
list, verified against the pinned SHA:

* generation never touches esmfold2 (``task/predict/`` has no reference to it);
* ``analysis`` is launched with ``esmfold2_metrics={not use_affinity}`` = true and
  ``analyze.py:638-652`` then unconditionally reads
  ``esmfold2_scores/<id>.input.json``, so it crashes unless that key is set false;
* ``filtering`` sets ``self.use_esmfold2 = not use_affinity`` at ``filter.py:218``
  -- **hard-wired to the protocol, reachable by no config key** -- and raises
  ``ValueError("Missing ESMFold2 scores; ...")``. It can only be omitted.

``--steps`` is therefore NOT a flag of this tool: anyone who could set it could only
set it wrong. ``--esmfold2_acceleration off`` does not help -- it selects a numerical
backend, it does not disable the model. ``--extra-flags`` and ``--config-set`` are
**guarded** (``BLOCKED_EXTRA_FLAGS``, ``check_config_set``): because the task script
interpolates them after the forced flags and argparse is last-wins, an unguarded
escape hatch would let ``--extra-flags '--protocol protein-anything'`` change what the
tool designs while every column name, docstring and runtime guard still said nanobody.

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

Usage:
    # root run: one BoltzGen run against a target not in any table yet
    sapia run boltzgen outputs/RUN -e modal \\
        --target-structure inputs/trimer.cif --target-chains A,B,C \\
        --binding A91,A128,B91,B128 --num-designs 50

    # child run: one BoltzGen run per target structure already in a table
    sapia run boltzgen outputs/RUN -t table0 -i pdb_path -e modal \\
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

# The protocol and step list are NOT configurable; see the module docstring.
PROTOCOL = "nanobody-anything"
STEPS = ("design", "inverse_folding", "folding", "analysis")
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
    target_structure: Path | None
    target_chains: str
    binding: str | None
    not_binding: str | None
    scaffolds: str
    scaffold_dir: str
    num_designs: int
    diffusion_batch_size: int | None
    num_workers: int
    spec: str | None
    extra_spec: str | None
    config_set: list[str]
    extra_flags: str


def add_run_boltzgen_args(parser: ArgumentParser) -> None:
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
        default=",".join(DEFAULT_SCAFFOLDS),
        help="Comma-separated VHH scaffold YAMLs BoltzGen draws from, AT RANDOM, one "
        "per design. Bare names resolve against --scaffold-dir; a name containing '/' "
        "is used as a path verbatim. Defaults to the five in "
        "example/nanobody/penguinpox.yaml (%(default)s); 8coh, 8z8v and gefurulimab "
        "also ship. Pass a single name to pin the framework -- but then every design "
        "shares it, which is the opposite of the diversity this protocol is for.",
    )
    parser.add_argument(
        "--scaffold-dir",
        type=str,
        default=os.getenv("BOLTZGEN_SCAFFOLD_DIR", DEFAULT_SCAFFOLD_DIR),
        help="Directory bare --scaffolds names resolve against, INSIDE the task "
        "container. Defaults to $BOLTZGEN_SCAFFOLD_DIR or %(default)s. Not checked "
        "for existence at submit time -- it does not exist where `sapia run` runs.",
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
        "script), e.g. '--use_kernels false --overwrite'. Escape hatch for flags "
        "without a dedicated option. "
        f"REFUSED at submit time: {', '.join(BLOCKED_EXTRA_FLAGS)} -- they are this "
        "tool's contract and argparse is last-wins, so passing them here would "
        "silently change what the tool IS while every column, docstring and guard "
        "still said nanobody.",
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
) -> list[dict[str, Any]]:
    """The spec's ``entities:`` list: the target file entity, then the scaffold list.

    The scaffold entity is a ``file:`` whose ``path:`` is a LIST -- that list is what
    makes BoltzGen pick a framework at random per design
    (``data/parse/schema.py:1947-1949``).
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

    return [
        {"file": target},
        {"file": {"path": scaffold_paths(args.scaffolds, args.scaffold_dir)}},
    ]


class SpecResult(NamedTuple):
    """One target's spec plus everything the table has to be told about it."""

    spec: dict[str, Any]
    binding_resnames: list[str]  # the epitope, in the caller's AUTHOR letters
    not_binding_resnames: list[str]  # the ANTI-epitope -- a different claim
    res_frame: str
    chain_frame: str
    chain_id_map: str  # 'A>C,B>D' -- author letter > the id BoltzGen resolves
    dropped_nonpolymer: list[str]  # subchains an explicit --target-chains leaves out


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

    spec: dict[str, Any] = {
        "entities": build_entities(target_path, args, binding, not_binding, chain_ids)
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
    )


def _reject_spec_building_flags(args: BoltzGenArgs) -> None:
    """--spec owns the whole YAML, so the flags that write one must be absent."""
    conflicting = [
        flag
        for flag, value in (
            ("--target-structure", args.target_structure),
            ("--binding", args.binding),
            ("--not-binding", args.not_binding),
        )
        if value is not None
    ]
    if args.target_chains.strip().lower() != "all":
        conflicting.append("--target-chains")
    if args.scaffolds != ",".join(DEFAULT_SCAFFOLDS):
        conflicting.append("--scaffolds")
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
    # Recorded PER GROUP, because a later submit with a different --scaffolds would
    # otherwise re-point an earlier group's allow-list at a set it never used.
    scaffolds = [] if verbatim else scaffold_paths(args.scaffolds, args.scaffold_dir)

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
                if state != "present" or not recorded:
                    raise ValueError(
                        f"{group_dir} already holds designs but its "
                        f"{SUBMIT_RECORD} is {state if state != 'present' else 'missing a spec_sha256'}, "
                        f"so there is no evidence they were made from this run's spec. "
                        f"BoltzGen is never run with --overwrite and those designs "
                        f"would survive to collect under this run's epitope. Use "
                        f"-l/--dir-label for a fresh leaf, or delete the WHOLE "
                        f"directory {group_dir} -- deleting {SUBMIT_RECORD} on its own "
                        f"leaves the old designs in place and is exactly what this "
                        f"check exists to catch."
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
        }
        group_dir.mkdir(parents=True, exist_ok=True)
        (group_dir / SUBMIT_RECORD).write_text(json.dumps(record, indent=2))

        manifest_rows.append(
            (name, str(volume_path(spec_yaml)), str(group_dir), flags)
        )
        submitted.append(name)

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

    ctx.write_meta(
        runs_root=str(runs_root),
        protocol=PROTOCOL,
        steps=list(STEPS),
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
