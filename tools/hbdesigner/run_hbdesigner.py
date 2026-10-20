#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) that designs BURIED HYDROGEN-BOND NETWORKS onto
each input backbone with HBDesigner, spawning a new child table.

PREMISE -- read this before using the numbers it writes.

HBDesigner (RosettaCommons/HBDesigner, MIT; PyRosetta is free for academic use) is
a GNN + PyRosetta algorithm that picks a small set of positions on an EXISTING
backbone and assigns polar residues there such that they form one connected,
buried hydrogen-bond network. It designs ``--n-res`` positions (2-6) and NOTHING
else: every other position in the output structure is written as GLYCINE. The
output is therefore a *network stub on a backbone*, not a foldable sequence. It is
an input to sequence design, not a replacement for it -- the intended chain is

    rfdiffusion3 -> hbdesigner -> proteinmpnn (network positions fixed) -> boltz

and the ``hbdesigner_fixed_positions`` / ``hbdesigner_fix<i>`` columns exist to
make that hand-off exact.

Its numbers are meaningless if:

  * you read the output as a designed protein. It is poly-glycine outside the
    network (``to_pdb(unk_to_gly=True)``, verified in ``rank_and_save``). Anything
    that computes on the sequence -- a predictor, a composition filter, a
    hydrophobicity score -- will see poly-Gly and be wrong, not broken-looking.
  * the backbone has no buried volume. HBDesigner buries a network; on a backbone
    with no core (a short helix, an extended peptide) it will find nothing and the
    run collects as ``error: no networks`` rather than as a bad score.
  * you designed one side of an interface without grafting the partner back in
    (see ``--graft-chains``) -- the partner chain returns as poly-glycine and any
    downstream interface metric is computed against a chain that is not there.

action ``create``: each input backbone yields up to ``--top-k`` RANKED designs,
each a different network, i.e. new entities -- so it mints a child table with rows
``<parent>_hb1 .. <parent>_hb<k>`` (mirroring proteinmpnn's ``<parent>_f<i>``).
Rank 1 is the best by HBDesigner's own ordering (fewest buried unsats, then
highest saturation, then lowest HB_Score_full).

WHAT IT DOES NOT DO: it does not design the rest of the sequence, it does not
relax or validate the fold, it does not predict a structure, and it does not score
an interface. It has no opinion on whether the network is worth having -- the
filters it applies (``--max-buns``, ``--min-sat``, ...) are the only judgement in
the tool, and they are yours to set.

RESIDUE NUMBERING -- the one thing that silently goes wrong here.
HBDesigner speaks PDB chain+resnum (``--guide-res A12,B13``, ``--anchor-res B5``,
and the ``network`` string it reports, e.g. ``A12S:A16T``). ProteinMPNN speaks
1-based index WITHIN each parsed chain. They agree only when the input PDB is
numbered 1..L per chain -- which rfdiffusion3 output happens to be, and a trimmed
crystal target is not. The collector therefore never copies HBDesigner's resnums
into ``fixed_positions``: it maps them through the output PDB's own residue order
and reports the mapping as ``resnum_offset`` / ``resnum_shift_max`` so you can see
whether a shift happened.

TWO BACKENDS, TWO HBDESIGNERS -- the numbers are not automatically comparable.

  * Modal runs the RosettaCommons checkout pinned in ``modal_image.py``
    (``ed65fa05``), with gemmi, biotite and hydride in the image.
  * vib runs the lab's conda env ``hbdesigner_gpu`` (activation
    ``SAPIA_ACTIVATE_HBDESIGNER``), which has the KUHLMAN-LAB checkout at
    ``b21ae46c`` and NO gemmi, NO biotite and NO hydride. The worker therefore
    reads PDBs with Biopython there instead of gemmi (``pdb_parser`` says which).

    The vib inference path, traced through the installed source:
    ``cli()`` (inference_hbdesigner.py:1085) -> ``HBDesRunner.run`` (266) ->
    ``validate_inputs`` sets ``opts.packer = "hbpacker"`` UNCONDITIONALLY
    (263-264; there is no CLI flag for it) -> ``pack_with_hbpacker`` (644-720)
    -> ``Pool.map(minimize_and_score_network)`` (704-714) ->
    ``minimize_and_score_network`` (inference/packing.py:35-155), which calls
    ``score_network`` (packing.py:122, defined 404-448) for ``HB_Score_full``,
    ``HB_Score_hb`` and ``Avg_Burial``, and ``get_satisfaction``
    (packing.py:124, defined data/hbnet.py:67-198) for ``saturation``,
    ``buried_heavy_unsats`` and ``buried_unsat_Hpol`` -> ``filter_packs``
    (887-892 applies --max_BUNs / --max_BUPHs / --min_sat / --max_hb_score) ->
    ``rank_and_save`` (382-474) writes the stats CSV.

    biotite and hydride are NOT needed on vib: no file under ``inference/``
    imports either, ``biotite_hbond_detect`` (hbnet.py:281) has no caller
    anywhere, and ``run_hydride`` is reachable only through
    ``pack_with_rosetta(mode="hydride")``, whose only callers are the training
    scripts. ``hbnet.score_protein`` and ``hbnet.rosetta_hbond_detect`` are
    likewise training-only (scripts/train_hbpacker.py:466,
    scripts/train_hbdesigner.py:804) and are NOT the vib scoring path.

    Every flag this builder emits exists in that version's ``--help``; the stats
    CSV it writes carries the columns the collector reads. What is NOT
    established is whether the MODAL pin scores through these same functions, so
    treat rows from the two backends as separate populations until a run has
    compared them on the same backbones.

``--guide-res`` NEEDS AT LEAST TWO DISTINCT RESIDUES. The guide atom is placed at
the Cb CENTROID of the guide residues, so one residue (or the same residue twice)
is rejected inside the task (``Only 1 provided. You must provide at least 2 guide
res for centroid calculation``). The builder refuses it at submit time instead,
after ``{expr}`` resolution, and refuses ``--guide-radius`` on a design whose
``--guide-res`` resolved to nothing. Residue tokens are ``<chain letter><number>``
only: no ranges, no insertion codes, no negative numbers. For a single residue use
``--anchor-res``, which forces that residue INTO every network.

Usage:
    # monomer: 3-residue network on each rfd3 backbone, 5 ranked designs each
    sapia run hbdesigner outputs/<run> --table table0 \\
        --n-res 3 --n-samples 200 --top-k 5

    # one-sided interface: network must include target residue B5, chain B is not
    # designable, and chain B's real sequence is grafted back afterwards
    sapia run hbdesigner outputs/<run> --table table0 \\
        -i rfdiffusion3_path \\
        --n-res 3 --n-samples 200 --anchor-res B5 --omit-chains B
"""

import re
from argparse import ArgumentParser
from pathlib import Path
from typing import cast

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path
from prosapia.utils import ensure_pdb, expand_chain_spec, resolve_template

# The worker's PDB reader, reused at submit time so "which chains does this
# structure have" means exactly the same thing here and in the task. Its
# module-level imports are standard library only, so this costs nothing.
from .hbdesigner_worker import ResidueMap

# Mirrors RESOURCES["cpu"] in modal_image.py. Used only when neither --n-workers
# nor -c/--cpus-per-task is given, so the packing pool matches the container.
DEFAULT_CPUS = 16

# --graft-chains sentinel: graft exactly the chains --omit-chains made poly-Gly.
GRAFT_AUTO = "auto"
GRAFT_NONE = "none"

DESIGN_MODELS = ("design_002", "design_020")

# upstream's --guide_radius default; the value that means "no distance cap".
GUIDE_RADIUS_OFF = 1e6

# The ONLY residue token HBDesigner's validate_residues() accepts. It splits each
# comma-separated token with re.split("(\\d+)") and unpacks exactly two pieces,
# then looks the first up in PDB_CHAIN_IDS -- so one chain LETTER followed by one
# run of digits, and nothing else: no ranges ('A12-15' unpacks to four pieces),
# no negative numbers ('A-5' yields chain 'A-'), no insertion codes ('A12A'
# yields a trailing piece). Checked at submit time because the failure otherwise
# costs a GPU task.
RESIDUE_TOKEN_RE = re.compile(r"^[A-Za-z]\d+$")

# A design name becomes a directory the task CREATES AND DELETES, a file name,
# and the prefix of every child row. One path component only: letters, digits,
# '.', '-' and '_', starting with a letter or digit. The task script repeats this
# check before its `rm -rf`, but refusing here means it never reaches a GPU.
DESIGN_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

# The minimum for --guide-res. Upstream refuses fewer than two TOKENS
# ("Only 1 provided. You must provide at least 2 guide res for centroid
# calculation.", validate_residues(mode='guide')); it counts tokens, so
# 'A12,A12' passes upstream. Requiring two DISTINCT residues is this tool's
# POLICY, not upstream's rule: the guide atom is the Cb centroid of the guide
# residues, and the centroid of one residue listed twice is that residue, so a
# repeated token is a typo rather than an instruction. A design whose
# --guide-res collapses to one residue ABORTS THE SUBMIT (it is not skipped):
# a silently dropped design reads as a filter problem later.
MIN_GUIDE_RES = 2

# upstream's `restype_hb_sc`: the only residues HBDesigner will place in a
# network, and therefore the whole alphabet --guide-seq and --omit-aa may use.
HB_SIDECHAIN_AAS = frozenset("RNDQEHKSTWY")


class HBDesignerArgs(CommonArgs):
    design_model: str
    cpu: bool
    n_workers: int | None
    n_samples: int
    top_k: int
    n_res: int
    t_range: list[float]
    min_burial: float
    guide_res: str
    guide_radius: float
    guide_seq: str
    max_buns: int
    max_buphs: int
    min_sat: float
    max_hb_energy: float
    max_hb_score: float
    min_core_res: int
    anchor_res: str
    symm_chains: str
    symm_file: str
    sel_chains: str
    omit_chains: str
    omit_aa: str
    seed: int | None
    graft_chains: str
    set: list[str]


def add_run_hbdesigner_args(parser: ArgumentParser) -> None:
    # --- model / sampling -------------------------------------------------
    parser.add_argument(
        "--design-model",
        choices=DESIGN_MODELS,
        default="design_020",
        help="Which design checkpoint to sample from. 'design_020' (default, "
        "moderate noise) is the general-purpose model; 'design_002' (low noise) "
        "is more conservative. Both ship in the image.",
    )
    parser.add_argument(
        "--n-res",
        type=int,
        choices=range(2, 7),
        default=2,
        metavar="{2..6}",
        help="Size of the designed network, in residues (2-6, default 2). This is "
        "the ONLY thing hbdesigner designs -- every other position comes back as "
        "glycine. Bigger networks pack far less often: upstream's advice is "
        "n_res=2 -> --n-samples 100, 3 -> 200, 4/5 -> 500, 6 -> 1000.",
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=100,
        help="Networks sampled before packing/scoring (default 100). Only a "
        "fraction survive packing -- upstream reports ~10 of 200 as typical -- so "
        "this is the main knob when a design returns nothing.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=5,
        help="How many ranked designs to KEEP per input backbone (default 5). You "
        "may get FEWER: upstream takes min(top_k, n_surviving), so the child table "
        "is not n_parents * top_k rows. Never assume a fixed fan-out.",
    )
    parser.add_argument(
        "--t-range",
        type=float,
        nargs=2,
        default=[0.1, 1.0],
        metavar=("LO", "HI"),
        help="Sampling temperature range (default 0.1 1.0). Lower = more "
        "conservative, higher = more diverse. Passed through as --T_range.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed. Off by default. Note upstream's own caveat: with "
        "parallel packing and Rosetta scoring, runs are NOT bit-reproducible even "
        "with a seed.",
    )

    # --- where the network may go ----------------------------------------
    parser.add_argument(
        "--min-burial",
        type=float,
        default=0.0,
        help="Minimum burial (Rosetta sidechain-neighbour count) for a position to "
        "be designable (default 0.0 = no constraint). Raise it to force the "
        "network into the core.",
    )
    parser.add_argument(
        "--min-core-res",
        type=int,
        default=0,
        help="Minimum number of CORE residues a network must contain (default 0).",
    )
    parser.add_argument(
        "--guide-res",
        default="",
        help="Residues whose Cb centroid a virtual guide atom is placed at, so the "
        "network forms near them. PDB chain+resnum, comma-joined: 'A3,A26'. "
        "AT LEAST TWO are required (a centroid of one is undefined; upstream "
        "refuses it and so does this builder, at submit time). Tokens are "
        "<chain letter><number> only: no ranges, no insertion codes, no negative "
        "numbers. {expr} islands are resolved per design up the lineage (e.g. "
        "'A{motif_end}') and the count is checked AFTER that. NOTE: PDB "
        "numbering, NOT the 1..L-per-chain indexing proteinmpnn's "
        "--fixed-positions uses. Using it makes HBDesigner append a marker atom "
        "(HETATM V1/ORI) to every output PDB; the task strips it and reports "
        "hbdesigner_guide_atoms_removed.",
    )
    parser.add_argument(
        "--guide-radius",
        type=float,
        default=1e6,
        help="Hard cap (Angstrom) on the Cb distance from the guide atom for a "
        "designable position. Default 1e6 = off. Meaningless without --guide-res "
        "(there is no guide atom to measure from), so passing it alone is a "
        "submit-time error rather than a run that quietly ignores it.",
    )
    parser.add_argument(
        "--anchor-res",
        default="",
        help="Residues that EVERY returned network must contain, PDB chain+resnum, "
        "comma-joined: 'B5'. One residue is fine here (unlike --guide-res), and "
        "tokens are <chain letter><number> only. {expr} islands resolved per "
        "design. This is what makes "
        "one-sided interface design work: anchor on the target residue you want "
        "hydrogen-bonded. Anchors must be polar "
        "(upstream rejects hydrophobic anchors). Anchors keep their input identity, "
        "so they appear in the reported network alongside the designed positions.",
    )

    # --- what the network may be -----------------------------------------
    parser.add_argument(
        "--guide-seq",
        default="",
        help="Per-position amino-acid conditioning, comma-joined, one entry per "
        "--n-res position: 'S,N,T' (exactly those), 'X,T,X' (at least one THR), "
        "'S,N|Q,T' (S, then N or Q, then T). Default = all X (unconstrained). "
        "Upstream's advice: SER/THR networks pack far more often, so 'S,X,X' is a "
        "cheap success-rate win when the identities do not matter.",
    )
    parser.add_argument(
        "--omit-aa",
        default="",
        help="Amino acids the network may NOT use, comma-joined: 'R,K'. Passed "
        "through as --omit_AA. Must not contradict --guide-seq.",
    )

    # --- chains -----------------------------------------------------------
    parser.add_argument(
        "--omit-chains",
        default="",
        help="Chains HBDesigner SEES but may not take network residues from "
        "(except via --anchor-res). Chain mini-language: ':' is an inclusive "
        "letter range, ',' separates ('B', 'A:C', 'A,C'). This is the one-sided "
        "interface-design flag. The omitted chain still comes back POLY-GLYCINE -- "
        "see --graft-chains, which by default repairs exactly these chains.",
    )
    parser.add_argument(
        "--sel-chains",
        default="",
        help="Chains to run HBDesigner on; every other chain is REMOVED from the "
        "input before design and put back once design has finished, in the "
        "input's own chain order (measured). Chain "
        "mini-language. Different from --omit-chains: here the model never sees "
        "the other chains at all, so it cannot avoid clashing with them. "
        "CAUTION, measured on the vib checkout: the unused chains come back "
        "POLY-GLYCINE, not with their own sequence, and this tool offers no way "
        "to restore them (graft_seq.py is refused here). It is also unsafe with "
        "--anchor-res, whose index is resolved before the extraction. Prefer "
        "cutting the input with chainsel and using --omit-chains.",
    )
    parser.add_argument(
        "--symm-chains",
        default="",
        help="Symmetrize the designed network across chains after design: 'A,B' "
        "ties A with B; 'A,B;C,D' ties two pairs independently. Off by default. "
        "Experimental upstream, and it can drop every network (see the zero-output "
        "trap in the skill).",
    )
    parser.add_argument(
        "--symm-file",
        default="",
        help="Rosetta .symm file for STRICT symmetry (made with "
        "make_symmdef_file.pl). Requires --symm-chains. The path must be readable "
        "from the task container -- put it inside the run_dir on the Volume.",
    )

    # --- acceptance filters (the only judgement in the tool) ---------------
    parser.add_argument(
        "--max-buns",
        type=int,
        default=0,
        help="Max buried unsatisfied HEAVY atoms allowed in a network (default 0, "
        "i.e. none). Passed as --max_BUNs. Raising it is more permissive.",
    )
    parser.add_argument(
        "--max-buphs",
        type=int,
        default=5,
        help="Max buried unsatisfied POLAR HYDROGENS allowed (default 5). Passed "
        "as --max_BUPHs.",
    )
    parser.add_argument(
        "--min-sat",
        type=float,
        default=0.5,
        help="Minimum network saturation, in [0, 2] (default 0.5). Higher is "
        "STRICTER -- the opposite direction from the two --max-* filters.",
    )
    parser.add_argument(
        "--max-hb-energy",
        type=float,
        default=0.0,
        help="Rosetta score:hb_max_energy -- the per-hydrogen-bond energy ceiling "
        "used while scoring (default 0.0).",
    )
    parser.add_argument(
        "--max-hb-score",
        type=float,
        default=0.0,
        help="Max HB_Score a returned network may have (default 0.0). More "
        "negative is stricter.",
    )

    # --- resources / post-processing / escape hatch ------------------------
    parser.add_argument(
        "--n-workers",
        type=int,
        default=None,
        help="Worker processes for CPU packing (DataLoader workers, and the pool "
        "used for lazy symmetrization). Default: match the task's CPU count "
        f"(-c/--cpus-per-task, else {DEFAULT_CPUS}). MUST be >= 1, on both "
        "backends: the Modal pin builds its DataLoader with "
        "persistent_workers=True, which raises on num_workers=0, and the vib "
        "version asserts '--n_workers must be >=1' itself (its --seed help text "
        "suggesting n_workers=0 contradicts that assert -- do not follow it).",
    )
    parser.add_argument(
        "--cpu",
        action="store_true",
        help="Run inference on CPU (loads the *_cpu.yaml configs) and request no "
        "GPU. Much slower; for debugging or a GPU-less backend.",
    )
    parser.add_argument(
        "--graft-chains",
        default=GRAFT_AUTO,
        help="After design, restore the original sequence AND sidechains of these "
        f"chains from the tool's own input PDB (upstream's graft_seq.py). "
        f"'{GRAFT_AUTO}' (default) grafts exactly the chains named by "
        f"--omit-chains, and does nothing otherwise; '{GRAFT_NONE}' never grafts; "
        "or give an explicit chain list. WITHOUT this, a one-sided interface "
        "design hands the next tool a target chain that is poly-glycine. Grafting "
        "only overwrites positions that came back as glycine, so designed and "
        "anchor residues are preserved. REFUSED together with --sel-chains: "
        "graft_seq.py copies BY POSITION between two files it asserts have the "
        "same residue count, and this tool does not establish that a "
        "--sel-chains output still lines up residue for residue with the input. "
        "Note that leaves the unused chains POLY-GLYCINE with no way to restore "
        "them here -- cut the input with chainsel and use --omit-chains instead.",
    )
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="TOKEN",
        help="Forward a raw run_hbdesigner flag verbatim, repeatable: "
        "--set --design_model_ckpt --set /runs/<path>.pt. The escape hatch for "
        "upstream flags this tool does not name (the two custom-checkpoint flags). "
        "Tokens may not contain whitespace. NOTE --n_workers is appended AFTER "
        "these tokens, so argparse takes the tool's value and a '--set "
        "--n_workers' of your own is silently overridden; use -c/--cpus-per-task "
        "or --n-workers instead.",
    )


def _chain_list(spec: str, flag: str) -> str:
    """Chain mini-language -> the comma-joined list HBDesigner wants ('A,C')."""
    spec = spec.strip()
    if not spec:
        return ""
    try:
        chains = expand_chain_spec(spec)
    except ValueError as e:
        raise ValueError(f"{flag}: {e}") from None
    if len(set(chains)) != len(chains):
        raise ValueError(f"{flag} names a chain twice: {chains}")
    return ",".join(chains)


def _residue_tokens(spec: str, flag: str, name: str) -> list[str]:
    """Split a resolved residue list and check every token HBDesigner will parse."""
    tokens = [t.strip() for t in spec.split(",") if t.strip()]
    for token in tokens:
        if RESIDUE_TOKEN_RE.match(token):
            continue
        hint = ""
        if "-" in token or ":" in token or ".." in token:
            hint = (
                " RANGES ARE NOT SUPPORTED: HBDesigner parses each token with "
                "re.split('(\\d+)') and takes exactly <chain letter><number>, so "
                "'A12-15' and 'A-5' both fail inside the task. List every residue "
                "separately: 'A12,A13,A14,A15'."
            )
        raise ValueError(
            f"{name}: {flag} token {token!r} is not <chain letter><residue number> "
            f"(e.g. 'A12').{hint}"
        )
    return tokens


def _check_residues_exist(
    name: str, rmap: ResidueMap, tokens: list[str], flag: str
) -> None:
    """Every resolved residue token must actually be in the input structure.

    ``validate_residues`` looks each token up in the scaffold and raises
    ``Parsing error! 0 matches found for res <token>.`` when it is absent -- but
    inside the task, after the GPU is allocated and the models are loaded. The
    same lookup costs nothing here.
    """
    missing = [t for t in tokens if (t[0], int(t[1:])) not in rmap.index]
    if not missing:
        return
    detail = []
    for t in missing:
        chain, num = t[0], int(t[1:])
        if chain not in rmap.seq:
            detail.append(f"{t} (no chain {chain}; chains: {','.join(rmap.chain_order)})")
        else:
            nums = rmap.resnums[chain]
            detail.append(
                f"{t} (chain {chain} has {len(nums)} residues numbered "
                f"{nums[0]}-{nums[-1]})"
            )
    raise ValueError(
        f"{name}: {flag} names residue(s) that are not in the input structure: "
        f"{'; '.join(detail)}. HBDesigner would raise 'Parsing error! 0 matches "
        f"found' for these inside the task, after the GPU had been allocated. "
        f"Note the numbering is the INPUT PDB's own, which for an rfd3 backbone "
        f"restarts at 1 per chain and need not match a crystal numbering."
    )


def _check_chain_span(
    name: str,
    input_pdb: Path,
    omit_chains: str,
    sel_chains: str,
    anchor_tokens: list[str],
    guide_tokens: list[str] | None = None,
) -> None:
    """Refuse a run whose network could never satisfy upstream's chain-span rule.

    ``filter_packs`` keeps a packed network only if
    ``r["n_chains"] == total_chains``, where ``total_chains`` is the number of
    chains in the scaffold AS HBDESIGNER SEES IT (after ``--sel_chains``
    extraction). A network must therefore touch EVERY chain of that scaffold.

    Residues can only be taken from chains that are designable (everything not
    in ``--omit-chains``) or named by ``--anchor-res``. So if some chain of the
    scaffold is neither, no sampled network can ever pass the filter and the run
    returns nothing -- after the GPU task has loaded the models and sampled.
    Nine runs were lost to this before the check existed.

    Monomer and two-sided design (no ``--omit-chains``) are unaffected: every
    chain is designable, so the set always matches.
    """
    rmap = ResidueMap(input_pdb)
    file_chains = list(rmap.chain_order)
    if not file_chains:
        raise ValueError(
            f"{name}: {input_pdb.name} has no amino-acid chains that the tool's "
            f"PDB reader recognises."
        )

    # The ANCHOR is resolved against the FULL file: upstream calls
    # validate_residues(mode='anchor') BEFORE extract_chains, so an anchor on a
    # chain that --sel-chains later removes still resolves at that point (it is
    # the index surviving the extraction that goes wrong, and the --sel-chains
    # guard further down covers that).
    _check_residues_exist(name, rmap, anchor_tokens, "--anchor-res")

    selected = [c.strip() for c in sel_chains.split(",") if c.strip()]
    if selected:
        missing = [c for c in selected if c not in file_chains]
        if missing:
            raise ValueError(
                f"{name}: --sel-chains names chain(s) {','.join(missing)}, which "
                f"are not in {input_pdb.name} (chains: {','.join(file_chains)})."
            )
        # extract_chains masks by chain index, so file order is preserved.
        scaffold = [c for c in file_chains if c in selected]
    else:
        scaffold = file_chains

    # The GUIDE residues are resolved AFTER extract_chains (run() line ~302
    # follows the --sel_chains block at ~291), so they must exist in the
    # SELECTED scaffold, not merely in the file. The asymmetry with --anchor-res
    # above is upstream's, not ours: a guide residue on a chain that
    # --sel-chains removed raises 'Parsing error! 0 matches found' inside the
    # task even though the residue is plainly in the input file.
    if guide_tokens:
        off_scaffold = [t for t in guide_tokens if t[0] not in scaffold]
        if off_scaffold:
            raise ValueError(
                f"{name}: --guide-res names residue(s) {','.join(off_scaffold)} on "
                f"chain(s) --sel-chains removed. HBDesigner resolves the guide "
                f"residues AFTER the extraction, against the selected scaffold "
                f"({','.join(scaffold)}), so these would raise 'Parsing error! 0 "
                f"matches found' inside the task even though they are in the "
                f"input file. Either select their chain or guide on a selected one."
            )
        # Chain membership is settled above; this checks the residue numbers
        # themselves, which --sel-chains does not change (extract_chains masks
        # rows, it does not renumber).
        _check_residues_exist(name, rmap, guide_tokens, "--guide-res")

    omitted = {c.strip() for c in omit_chains.split(",") if c.strip()}
    anchor_chains = {t[0] for t in anchor_tokens}
    reachable = {c for c in scaffold if c not in omitted} | (
        anchor_chains & set(scaffold)
    )
    unreachable = [c for c in scaffold if c not in reachable]
    if unreachable:
        raise ValueError(
            f"{name}: no network could ever pass HBDesigner's own filter with "
            f"these flags. Upstream's filter_packs keeps a network only when it "
            f"touches EVERY chain of the scaffold (n_chains == total_chains). "
            f"The scaffold here is {','.join(scaffold)}"
            + (f" (--sel-chains cut it from {','.join(file_chains)})" if selected else "")
            + f", but chain(s) {','.join(unreachable)} are omitted from design "
            f"and carry no --anchor-res, so no residue can be placed on them and "
            f"every sampled network is discarded. Three ways out: (1) cut the "
            f"input to just the chains you want, e.g. "
            f"'sapia run chainsel ... --chains {','.join(sorted(reachable)) or 'X,D'}' "
            f"and run hbdesigner on chainsel_path; (2) put an --anchor-res on "
            f"every omitted chain; or (3) do not omit {','.join(unreachable)}."
        )

    # --sel-chains + --anchor-res: the anchor index is resolved against the FULL
    # scaffold (run() line ~276) and used after extract_chains has shrunk it
    # (featurize_inference: `p.aatype[anchor_res] = p_copy.aatype[anchor_res]`),
    # so an index past the new length raises IndexError inside the task. Only an
    # anchor in the file's FIRST chain is guaranteed to keep its index.
    # The general invariant is: every chain preceding the anchor's chain in the
    # FILE is selected, and the anchor's chain is itself selected. Only then does
    # the absolute index survive the extraction unchanged. The safe subset
    # implemented here is the anchor's chain being the file's FIRST chain and
    # selected; anything else is refused rather than reasoned about, because the
    # failure is silent in one direction: an anchor on an unselected chain does
    # NOT raise, it pins whatever residue happens to sit at that index in the
    # extracted scaffold (e.g. --sel-chains B,D --anchor-res A50 would anchor
    # B50), which is a plausible wrong answer rather than an error.
    if selected and anchor_tokens:
        bad = [
            t
            for t in anchor_tokens
            if t[0] != file_chains[0] or t[0] not in selected
        ]
        if bad:
            raise ValueError(
                f"{name}: --sel-chains cannot be combined with --anchor-res "
                f"{','.join(bad)}. HBDesigner resolves the anchor to an absolute "
                f"residue index on the FULL scaffold ({','.join(file_chains)}) "
                f"and applies it after --sel_chains has cut it down to "
                f"{','.join(scaffold)}. If the index lands past the end the task "
                f"dies with 'IndexError: index N is out of bounds for axis 0 "
                f"with size M' (observed); if it lands inside, it silently pins "
                f"the WRONG residue -- an anchor on an unselected chain is not an "
                f"error upstream, it just anchors whatever sits at that index. "
                f"Only an anchor in the file's first chain ({file_chains[0]}), "
                f"with that chain selected, keeps its index. Instead, cut the "
                f"input first with chainsel so the chains you want are the whole "
                f"file, then use --omit-chains, not --sel-chains."
            )


def _check_guide_seq(args: "HBDesignerArgs") -> None:
    """The --guide-seq / --omit-aa rules upstream asserts, checked at submit time.

    All three come straight from the installed ``validate_inputs``: the guide
    sequence is split on ',' and must have exactly ``--n-res`` entries (its
    default is ``["X"] * n_res``, NOT the 'X,X,X' its --help claims); every
    alternative in an entry must be in ``restype_hb_sc`` or 'X'; ``--omit-aa``
    must be a subset of ``restype_hb_sc``; and no guide-sequence letter may also
    be omitted. Each of these otherwise fails inside the task, after the GPU has
    been allocated and the models loaded.
    """
    # Checked on the RAW token: upstream compares against `restype_hb_sc` with no
    # case folding, so a lowercase 's' is rejected there, not silently accepted.
    omit = [a.strip() for a in args.omit_aa.split(",") if a.strip()]
    bad_omit = [a for a in omit if a not in HB_SIDECHAIN_AAS]
    if bad_omit:
        raise ValueError(
            f"--omit-aa: {','.join(bad_omit)} is not a hydrogen-bonding amino "
            f"acid. HBDesigner only ever places {','.join(sorted(HB_SIDECHAIN_AAS))}, "
            f"so only those can be omitted (upper case: upstream does not fold "
            f"case, so 'r' is as invalid there as 'Z')."
        )

    spec = args.guide_seq.strip()
    if not spec:
        return
    entries = spec.split(",")
    if len(entries) != args.n_res:
        raise ValueError(
            f"--guide-seq has {len(entries)} entr{'y' if len(entries) == 1 else 'ies'} "
            f"({spec!r}) but --n-res is {args.n_res}: upstream requires exactly one "
            f"entry per designed position. Use 'X' for a position you do not want "
            f"to constrain, e.g. {','.join(['X'] * args.n_res)}."
        )
    for entry in entries:
        for alt in entry.split("|"):
            alt = alt.strip()
            if alt != "X" and alt not in HB_SIDECHAIN_AAS:
                raise ValueError(
                    f"--guide-seq entry {entry!r} contains {alt!r}, which is "
                    f"neither 'X' nor one of {','.join(sorted(HB_SIDECHAIN_AAS))} "
                    f"(the only residues HBDesigner places in a network). Upper "
                    f"case only -- upstream does not fold case."
                )
            if alt in omit:
                raise ValueError(
                    f"--guide-seq asks for {alt!r} at entry {entry!r} while "
                    f"--omit-aa forbids it. Drop it from one of the two."
                )


def _symm_chain_groups(spec: str) -> str:
    """'A,B;C,D' -> each ';' group expanded through the chain mini-language."""
    spec = spec.strip()
    if not spec:
        return ""
    groups = [g for g in (part.strip() for part in spec.split(";")) if g]
    if not groups:
        raise ValueError("--symm-chains is empty after parsing")
    return ";".join(_chain_list(group, "--symm-chains") for group in groups)


def _resolve_graft_chains(args: HBDesignerArgs, omit_chains: str) -> str:
    """Which chains the task script grafts back, as a comma-joined list ('' = none)."""
    spec = args.graft_chains.strip()
    if spec.lower() == GRAFT_NONE:
        return ""
    chains = omit_chains if spec.lower() == GRAFT_AUTO else _chain_list(spec, "--graft-chains")
    if chains and args.sel_chains.strip():
        raise ValueError(
            "--graft-chains cannot be combined with --sel-chains. --sel-chains "
            "removes the unused chains before design and puts them back once "
            "design has finished, and upstream's graft_seq.py copies BY POSITION between "
            "two files it asserts have the same residue count -- so a graft here "
            "is only correct if the two still line up residue for residue, which "
            "this tool does not establish. Pass --graft-chains none. Note that "
            "this leaves the unused chains POLY-GLYCINE (measured on the vib "
            "checkout: a --sel_chains D run returned chains A,B,C,D with "
            "non-glycine residues only at the three network positions), and this "
            "tool then has no way to restore them. If you need those chains with "
            "their sequence, cut the input with chainsel and use --omit-chains, "
            "whose graft-back is supported."
        )
    return chains


def _flag_tokens(args: HBDesignerArgs, guide_res: str, anchor_res: str) -> list[str]:
    """The run_hbdesigner argv (minus --pdb/--out_dir) for one design.

    Every value is passed explicitly, defaults included, so the submitted command
    is fully readable in the manifest and an upstream default change cannot move
    a run's meaning underneath it.
    """
    tokens: list[str] = [
        "--design_model", args.design_model,
        "--n_res", str(args.n_res),
        "--n_samples", str(args.n_samples),
        "--top_k", str(args.top_k),
        "--T_range", str(args.t_range[0]), str(args.t_range[1]),
        "--min_burial", str(args.min_burial),
        "--min_core_res", str(args.min_core_res),
        "--guide_radius", str(args.guide_radius),
        "--max_BUNs", str(args.max_buns),
        "--max_BUPHs", str(args.max_buphs),
        "--min_sat", str(args.min_sat),
        "--max_hb_energy", str(args.max_hb_energy),
        "--max_hb_score", str(args.max_hb_score),
    ]
    optional: list[tuple[str, str]] = [
        ("--guide_res", guide_res),
        ("--guide_seq", args.guide_seq.strip()),
        ("--anchor_res", anchor_res),
        ("--omit_AA", args.omit_aa.strip()),
        ("--sel_chains", _chain_list(args.sel_chains, "--sel-chains")),
        ("--omit_chains", _chain_list(args.omit_chains, "--omit-chains")),
        ("--symm_chains", _symm_chain_groups(args.symm_chains)),
    ]
    for flag, value in optional:
        if value:
            tokens += [flag, value]
    if args.symm_file.strip():
        tokens += ["--symm_file", str(volume_path(Path(args.symm_file.strip())))]
    if args.seed is not None:
        tokens += ["--seed", str(args.seed)]
    if args.cpu:
        tokens.append("--cpu")
    tokens += [str(tok) for tok in args.set]
    return tokens


def build_hbdesigner_manifest(ctx: ManifestCtx[HBDesignerArgs]) -> list[tuple[str, ...]]:
    args = ctx.args

    if args.cpu:
        args.gpus_per_task = 0

    # Fail loudly on a column that matches nothing, instead of repeating the known
    # proteinmpnn bug where a wrong default silently submits zero tasks.
    column = args.input_column
    if column not in ctx.df.columns:
        available = ", ".join(str(c) for c in ctx.df.columns if str(c).endswith("_path"))
        raise ValueError(
            f"hbdesigner: table {args.table!r} has no column {column!r}. Pass "
            f"-i/--input-column with the backbone column to design onto. Structure "
            f"columns available: {available or '(none)'}."
        )

    if args.n_workers is not None and args.n_workers < 1:
        raise ValueError(
            f"--n-workers must be at least 1 (got {args.n_workers}): upstream builds "
            f"its packing DataLoader with persistent_workers=True, which raises on "
            f"num_workers=0. Use --n-workers 1 for 'no parallelism'."
        )
    n_workers = args.n_workers or args.cpus_per_task or DEFAULT_CPUS

    if args.symm_file.strip() and not args.symm_chains.strip():
        raise ValueError(
            "--symm-file is only used for strict symmetry and needs --symm-chains "
            "to say which chains it ties."
        )
    if args.top_k < 1:
        raise ValueError(f"--top-k must be at least 1 (got {args.top_k})")

    if args.guide_radius <= 0:
        raise ValueError(
            f"--guide-radius must be greater than 0 (got {args.guide_radius}); "
            f"upstream asserts 0.0 < guide_radius."
        )

    _check_guide_seq(args)

    omit_chains = _chain_list(args.omit_chains, "--omit-chains")
    # Expanded through the chain mini-language once, here, so every consumer
    # below sees 'A,B,C' and never the raw 'A:C'.
    sel_chains = _chain_list(args.sel_chains, "--sel-chains")
    graft_chains = _resolve_graft_chains(args, omit_chains)

    manifest_rows: list[tuple[str, ...]] = []
    for name in ctx.ready.index:
        name = cast(str, name)
        if not DESIGN_NAME_RE.match(name):
            raise ValueError(
                f"design name {name!r} is not usable by this tool: a name must "
                f"start with a letter or digit and hold only letters, digits, "
                f"'.', '-' and '_'. It becomes a directory that the task CREATES "
                f"AND DELETES under the output dir, a staged file name, and the "
                f"prefix of every child row, so a '/', a space or a '..' in it "
                f"would reach outside this design's own files. Rename the row, or "
                f"exclude it with -f."
            )
        src = Path(str(ctx.ready.at[name, column]))
        if not src.exists():
            print(f"{name}: MISSING {src} (skipping)")
            continue
        # HBDesigner reads PDB only (Protein.from_pdb_file).
        input_pdb = ensure_pdb(src, args.run_dir)

        # Per-design {expr} resolution happens HERE, so a bad column name kills the
        # submit instead of 100 tasks. Both flags take HBDesigner's own
        # chain+resnum syntax verbatim; only the {...} islands are ours.
        guide_res = resolve_template(args.guide_res.strip(), ctx.lookup, name)
        anchor_res = resolve_template(args.anchor_res.strip(), ctx.lookup, name)

        # Checked AFTER {expr} resolution, because that is what decides what is
        # actually emitted: '{hotspot_a},{hotspot_b}' may collapse to one residue,
        # or to the SAME residue twice, on some rows and not others.
        guide_tokens = _residue_tokens(guide_res, "--guide-res", name)
        distinct = len(set(guide_tokens))
        if guide_tokens and distinct < MIN_GUIDE_RES:
            repeated = "" if distinct == len(guide_tokens) else (
                f" They are not distinct: {','.join(guide_tokens)} names "
                f"{distinct} residue. Two copies of one residue have the same Cb "
                f"centroid as one copy."
            )
            raise ValueError(
                f"{name}: --guide-res resolved to {distinct} distinct residue "
                f"({','.join(sorted(set(guide_tokens)))}), but HBDesigner needs at "
                f"least {MIN_GUIDE_RES}: the guide atom is placed at the Cb "
                f"CENTROID of the guide residues, which is undefined for a single "
                f"one (upstream: 'You must provide at least 2 guide res for "
                f"centroid calculation').{repeated} Add a second residue near the "
                f"one you want the network steered to, or drop --guide-res and use "
                f"--anchor-res, which takes one residue and forces it INTO every "
                f"network."
            )
        # The radius is only applied when a guide atom exists, and whether one
        # exists is decided per design by what --guide-res RESOLVED to -- so this
        # is checked here rather than once against the unresolved template.
        if args.guide_radius != GUIDE_RADIUS_OFF and not guide_tokens:
            raise ValueError(
                f"{name}: --guide-radius {args.guide_radius} was given but "
                f"--guide-res resolved to nothing for this design "
                f"({args.guide_res.strip()!r}). The radius is a hard cap on the Cb "
                f"distance FROM THE GUIDE ATOM, and there is no guide atom without "
                f"--guide-res, so upstream would ignore it and the run would "
                f"silently not be the one you asked for. Pass --guide-res with at "
                f"least {MIN_GUIDE_RES} residues, or drop --guide-radius."
            )
        anchor_tokens = _residue_tokens(anchor_res, "--anchor-res", name)

        # Reads the input once per design: the chains are a property of the FILE,
        # and a table can mix inputs with different chain sets. The chain specs
        # are passed EXPANDED ('A:C' -> 'A,B,C'), never raw, or a valid range
        # spec would be read as a single chain named 'A:C'.
        _check_chain_span(
            name, input_pdb, omit_chains, sel_chains, anchor_tokens, guide_tokens
        )

        tokens = _flag_tokens(args, guide_res, anchor_res) + ["--n_workers", str(n_workers)]
        bad = [tok for tok in tokens if any(ch.isspace() for ch in tok)]
        if bad:
            raise ValueError(
                f"{name}: flag token(s) {bad} contain whitespace; the task script "
                f"word-splits the manifest's argument field, so they cannot be "
                f"passed. Remove the spaces (e.g. --guide-seq 'S,N|Q,T', not 'S, N')."
            )

        manifest_rows.append(
            (
                name,
                str(volume_path(input_pdb)),
                # May be empty, so never last.
                graft_chains,
                omit_chains,
                # Always non-empty: the manifest line never ends on an empty field.
                " ".join(tokens),
            )
        )

    return manifest_rows
