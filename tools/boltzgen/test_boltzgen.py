#!/usr/bin/env python3
"""
Self-test for the boltzgen tool.

WHAT THESE CHECKS PROVE
  * the residue/chain spec parsers and the two position renderers (BoltzGen's
    ``26..33`` and prosapia's ``26:33``) are arithmetically right, including the
    complement that feeds ``proteinmpnn --fixed-positions``;
  * **the epitope frame is the DECLARED polymer sequence, not the modelled ordinal.**
    A SEQRES fixture with an unmodelled residue proves position 3 resolves to the
    declared ASP and not to the third modelled residue (GLU), that an unmodelled
    epitope position ERRORS, and that a position past the declared sequence errors;
  * an epitope on a chain ``--target-chains`` excludes errors at submit time;
  * the written design YAML has the shape BoltzGen's schema wants -- in particular a
    scaffold entity whose ``path`` is a LIST, which is what makes the framework a
    random draw per design;
  * ``map_binder`` finds the nanobody chain from the design-mask B-factor, carries it
    into a refolded complex whose chain has been RENAMED, and refuses -- with an
    error, never a plausible answer -- a designed-residue count mismatch, an ABSENT
    third witness, two designed chains, a **length-preserving substitution** between
    design and refold, and a sequence that is not the chain's (with the measured
    identity carried onto the error row);
  * **the scaffold name is never a positional guess**: a source outside the run's
    recorded scaffold list, or one that is also a target chain's source, yields NA
    with a verdict that is not a success word -- including with no allow-list at all;
  * a missing ``output_chain`` is its own verdict, not disguised as a rename;
  * the refolded complex's pLDDT column is never read as a design mask
    (``designed`` is ``None`` there, not "everything is designed");
  * the forced step/config contract never leaks into a user-settable flag, and
    ``--extra-flags`` / ``--config-set`` REFUSE every token that would redefine it;
  * the collector turns a corrupt ``.npz``, an unparseable mmCIF, missing files, a
    duplicated ``id`` and a lost epitope record into ERROR ROWS -- the whole collect
    still completes, and every error row has an empty path and an empty ``cif_path``
    when the file is absent.

WHAT THEY CANNOT PROVE -- these are unknown until a real run
  * that the pinned image BUILDS, that torch 2.10.0 + cuequivariance-ops-cu12 0.11.0 +
    numpy 2.4.6 + numba 0.65.1 + pandas 3.0.6 + gemmi 0.7.5 import together, or that
    boltzgen runs on an H100;
  * that a REAL BoltzGen design CIF carries the design mask in its B-factor column
    the way ``writer.py`` says it does, and that the refolded complex really has the
    same chain ORDER as the design CIF. Both are assumptions this collector is built
    on and both are asserted here only against structures this file wrote;
  * that ``source_context["chains"]`` is ordered by OUTPUT CHAIN rather than by input
    entity. The guards above make a wrong order visible (NA + a verdict) instead of
    silent, which is the most a test can do without real output;
  * that ``boltzgen check`` accepts the YAML this tool writes -- the spec nesting is
    unverified against upstream's schema, which is exactly why the task script runs
    ``check`` before spending GPU time;
  * that ``aggregate_metrics_analyze.csv`` really carries every CORE_METRICS name
    under this step set (the absent ones are predicted from the code, not observed),
    nor any observed RANGE for any of them -- no threshold anywhere is calibrated;
  * that ``designed_chain_sequence`` is colon-free on a real nanobody run;
  * that the five default scaffold YAMLs resolve inside the image.

Run:  uv run python tools/boltzgen/test_boltzgen.py
"""

import importlib.util
import json
import sys
import tempfile
from argparse import ArgumentParser
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
import yaml

from prosapia.core import CollectCtx, DesignCtx

HERE = Path(__file__).resolve().parent

# The scaffold allow-list the run would have recorded in .meta.json.
STEMS = frozenset({"gontivimab", "7eow", "7xl0"})

# A minimal mmCIF whose AUTHOR chain letter ('H') differs from its label_asym_id ('A') --
# the 1TUP situation, self-contained so the test owns it. BoltzGen resolves a YAML
# `chain: id:` as label_asym_id for an mmCIF, so a caller writing 'H' must reach 'A'.
DIVERGING_CIF = """data_test
_entry.id test
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer test
loop_
_entity_poly.entity_id
_entity_poly.type
1 polypeptide(L)
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1 ALA n
1 2 CYS n
1 3 ASP n
1 4 GLU n
1 5 PHE n
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 ? 0.000 0.000 0.000 1.00 0.00 101 H 1
ATOM 2 C CA . CYS A 1 2 ? 3.800 0.000 0.000 1.00 0.00 102 H 1
ATOM 3 C CA . ASP A 1 3 ? 7.600 0.000 0.000 1.00 0.00 103 H 1
ATOM 4 C CA . GLU A 1 4 ? 11.400 0.000 0.000 1.00 0.00 104 H 1
ATOM 5 C CA . PHE A 1 5 ? 15.200 0.000 0.000 1.00 0.00 105 H 1
#
"""

THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}

PASS = FAIL = 0


def _load(mod: str):
    spec = importlib.util.spec_from_file_location(f"_bg_{mod}", HERE / f"{mod}.py")
    assert spec and spec.loader
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


_RUN = None  # set in main(); TF() is a convenience for the many resolve_* calls


def TF(path: Path):
    """label_seq_residues' TargetFrames for a path (tests call it repeatedly)."""
    assert _RUN is not None
    return _RUN.label_seq_residues(path)


def check(label: str, got, want) -> None:
    global PASS, FAIL
    if got == want:
        PASS += 1
        print(f"  ok    {label}")
    else:
        FAIL += 1
        print(f"  FAIL  {label}\n          got  {got!r}\n          want {want!r}")


def raises(label: str, fn, substring: str) -> None:
    global PASS, FAIL
    try:
        fn()
    except Exception as e:  # noqa: BLE001 -- the point is what message comes out
        if substring.lower() in str(e).lower():
            PASS += 1
            print(f"  ok    {label}")
        else:
            FAIL += 1
            print(f"  FAIL  {label}\n          raised {e!r}\n          wanted {substring!r}")
        return
    FAIL += 1
    print(f"  FAIL  {label} -- did not raise")


def write_pdb(path: Path, chains: dict[str, tuple[str, list[float]]]) -> Path:
    """One CA atom per residue; ``chains`` maps label -> (sequence, per-residue bfac)."""
    lines, serial = [], 1
    for label, (seq, bfacs) in chains.items():
        for i, (aa, bfac) in enumerate(zip(seq, bfacs), start=1):
            lines.append(
                f"ATOM  {serial:5d}  CA  {THREE[aa]} {label}{i:4d}    "
                f"{i * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{bfac:6.2f}           C"
            )
            serial += 1
        lines.append("TER")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return path


def write_cif(path: Path, chains: dict[str, tuple[str, list[float]]]) -> Path:
    """The same structure as write_pdb, re-emitted as real mmCIF (ensure_pdb wants one)."""
    pdb = write_pdb(path.with_name(path.stem + "_tmp.pdb"), chains)
    st = gemmi.read_structure(str(pdb))
    st.setup_entities()
    st.make_mmcif_document().write_file(str(path))
    pdb.unlink()
    return path


def write_npz(path: Path, n_designed: int, n_other: int, sources: list[str],
              output_chains: list[str]) -> Path:
    mask = np.array([1.0] * n_designed + [0.0] * n_other, dtype=np.float32)
    context = {
        "version": 1,
        "chains": [
            {"source": s, "output_chain": c} for s, c in zip(sources, output_chains)
        ],
    }
    np.savez(path, design_mask=mask, source_context=np.asarray(json.dumps(context)))
    return path


def main() -> int:
    global _RUN
    run = _RUN = _load("run_boltzgen")
    col = _load("collect_boltzgen")
    img = _load("modal_image")

    print("\n-- residue / chain spec parsing ------------------------------------")
    check(
        "singles and a range expand in order",
        run.parse_residue_spec("A91,A128,B5-8", "--binding"),
        {"A": [91, 128], "B": [5, 6, 7, 8]},
    )
    check(
        "A5-8 and B5-8 stay on their own chains",
        run.parse_residue_spec("A5-6,B7-8", "--binding"),
        {"A": [5, 6], "B": [7, 8]},
    )
    raises("a cross-chain range is refused",
           lambda: run.parse_residue_spec("A5-B9", "--binding"), "spans two chains")
    raises("a reversed range is refused",
           lambda: run.parse_residue_spec("A9-5", "--binding"), "ends before")
    raises("a duplicate residue is refused",
           lambda: run.parse_residue_spec("A5,A5", "--binding"), "listed twice")
    raises("a bare number is refused",
           lambda: run.parse_residue_spec("91", "--binding"), "malformed")
    raises("an empty spec is refused",
           lambda: run.parse_residue_spec(" , ", "--binding"), "empty")

    print("\n-- position rendering ---------------------------------------------")
    check("BoltzGen res_index syntax collapses runs",
          run._ranges([26, 27, 28, 51, 52, 97]), "26..28,51..52,97")
    check("prosapia position syntax collapses runs",
          col.render_positions([26, 27, 28, 51, 52, 97]), "26:28,51:52,97")
    check("a single position stays bare", col.render_positions([7]), "7")
    check("an empty list renders empty", col.render_positions([]), "")
    # The complement is what `proteinmpnn --fixed-positions` consumes.
    designed = [3, 4, 5, 9]
    fixed = [i for i in range(1, 11) if i not in set(designed)]
    check("designed + fixed partition the chain",
          (col.render_positions(designed), col.render_positions(fixed)),
          ("3:5,9", "1:2,6:8,10"))

    print("\n-- sequence identity ----------------------------------------------")
    check("identical sequences", col.identity("ACDEF", "ACDEF"), 1.0)
    check("one substitution in five", col.identity("ACDEF", "ACDEG"), 0.8)
    check("an empty sequence is 0", col.identity("", "ACDEF"), 0.0)
    # Unequal lengths go through the alignment; a fragment must NOT read as 1.0.
    check("a 5-mer inside a 10-mer is 0.5, not 1.0",
          round(col.identity("ACDEF", "ACDEFGHIKL"), 6), 0.5)
    check("a gapped alignment recovers the matches",
          round(col.identity("ACDEF", "ACDXEF"), 6), round(5 / 6, 6))

    print("\n-- scaffold paths -------------------------------------------------")
    check("bare names resolve against --scaffold-dir",
          run.scaffold_paths("7eow,gontivimab", "/opt/bg/sc"),
          ["/opt/bg/sc/7eow.yaml", "/opt/bg/sc/gontivimab.yaml"])
    check("a .yaml suffix is not doubled",
          run.scaffold_paths("7eow.yaml", "/opt/bg/sc"), ["/opt/bg/sc/7eow.yaml"])
    check("a path is taken verbatim",
          run.scaffold_paths("/my/own/vhh.yaml", "/opt/bg/sc"), ["/my/own/vhh.yaml"])
    raises("an empty scaffold list is refused",
           lambda: run.scaffold_paths(" ", "/opt/bg/sc"), "at least one")
    check("the default list is penguinpox's five",
          run.DEFAULT_SCAFFOLDS,
          ("7eow", "7xl0", "gontivimab", "isecarosmab", "sonelokimab"))

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)

        print("\n-- label_seq_id resolution ----------------------------------------")
        # A target whose AUTHOR numbering would start at 1 but which we always read as
        # a 1-based ordinal within the chain's polymer.
        target = write_pdb(
            tmpdir / "target.pdb",
            {"A": ("ACDEFGHIKL", [0.0] * 10), "B": ("MNPQRSTVWY", [0.0] * 10)},
        )
        tframes = run.label_seq_residues(target).chains
        check("label_seq_id 1..N per chain", tframes["A"].residues[1].name, "ALA")
        check("second chain is numbered from 1 again",
              tframes["B"].residues[1].name, "MET")
        check("the resolved names are reported in the order listed",
              run.resolve_binding_resnames(
                  {"A": [1, 3], "B": [2]}, TF(target), target, "--binding")[0],
              ["A1=ALA", "A3=ASP", "B2=ASN"])
        check("a file with no declared sequence says so, PER CHAIN",
              run.resolve_binding_resnames(
                  {"A": [1]}, TF(target), target, "--binding")[1],
              "synthesized_from_model")
        raises("a position past the end of the chain is an error",
               lambda: run.resolve_binding_resnames(
                   {"A": [11]}, TF(target), target, "--binding"),
               "absent from chain A")
        raises("and the error names the FRAME, not just the position",
               lambda: run.resolve_binding_resnames(
                   {"A": [11]}, TF(target), target, "--binding"),
               "frame synthesized_from_model")
        raises("and tells the reader a HET residue reports the same way",
               lambda: run.resolve_binding_resnames(
                   {"A": [11]}, TF(target), target, "--binding"),
               "HET or modified residue")
        raises("an absent chain is an error",
               lambda: run.resolve_binding_resnames(
                   {"Z": [1]}, TF(target), target, "--binding"),
               "no polymer")
        raises("and says chain ids ARE the author letters a viewer shows",
               lambda: run.resolve_binding_resnames(
                   {"Z": [1]}, TF(target), target, "--binding"),
               "AUTHOR letters a viewer shows")

        print("\n-- the written design spec ----------------------------------------")
        class Args:
            target_chains = "A,B"
            scaffolds = ",".join(run.DEFAULT_SCAFFOLDS)
            scaffold_dir = "/opt/bg/sc"

        IDS = {"A": "A", "B": "B"}
        entities = run.build_entities(
            Path("/runs/t.cif"), Args(), {"A": [91, 92, 128]}, {"B": [5]}, IDS
        )
        check("two entities: target then scaffolds", len(entities), 2)
        check("the target include list follows --target-chains",
              entities[0]["file"]["include"],
              [{"chain": {"id": "A"}}, {"chain": {"id": "B"}}])
        check("binding and not_binding land on their own chains",
              entities[0]["file"]["binding_types"],
              [{"chain": {"id": "A", "binding": "91..92,128"}},
               {"chain": {"id": "B", "not_binding": "5"}}])
        check("the scaffold entity's path is a LIST (this is the random draw)",
              isinstance(entities[1]["file"]["path"], list), True)
        check("all five scaffolds are offered",
              len(entities[1]["file"]["path"]), 5)
        plain = run.build_entities(Path("/runs/t.cif"), Args(), None, None, IDS)
        check("no --binding means no binding_types block",
              "binding_types" in plain[0]["file"], False)
        class AllChains(Args):
            target_chains = "all"
        check("--target-chains all means no include block",
              "include" in run.build_entities(
                  Path("/runs/t.cif"), AllChains(), None, None, IDS)[0]["file"],
              False)
        # It has to survive a YAML round trip, since that is how it reaches BoltzGen.
        check("the spec round-trips through YAML",
              yaml.safe_load(yaml.safe_dump({"entities": entities}, sort_keys=False)),
              {"entities": entities})

        print("\n-- binder-chain mapping -------------------------------------------")
        # A 2-chain target (A, B) plus a 12-residue nanobody whose residues 4..6 and
        # 10 are designed. The design CIF carries the design mask in its B-factors.
        nb_seq = "QVQLVESGGGL"  # 11 residues
        designed_flags = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        design = write_pdb(
            tmpdir / "design.pdb",
            {
                "A": ("ACDEFGHIKL", [0.0] * 10),
                "B": ("MNPQRSTVWY", [0.0] * 10),
                "C": (nb_seq, designed_flags),
            },
        )
        # The refolded complex: SAME chain order, but the writer re-pooled the label
        # from C to Z, and its B-factors are pLDDT, not the design mask.
        refold = write_pdb(
            tmpdir / "refold.pdb",
            {
                "A": ("ACDEFGHIKL", [88.0] * 10),
                "B": ("MNPQRSTVWY", [91.0] * 10),
                "Z": (nb_seq, [77.0] * 11),
            },
        )
        npz = write_npz(
            tmpdir / "design.npz",
            n_designed=4,
            n_other=27,
            sources=["/t/t.cif", "/t/t.cif", "/opt/bg/sc/gontivimab.cif"],
            output_chains=["A0", "B0", "C0"],
        )
        m = col.map_binder(design, npz, refold, nb_seq, 4, STEMS)
        check("the binder chain is the REFOLD's label, not the design's", m.chain, "Z")
        check("the designed positions are 1-based within the binder chain",
              m.design_positions, [4, 5, 6, 10])
        check("the scaffold is the basename of source", m.scaffold, "gontivimab.cif")
        check("source_context's output_chain confirms the positional link",
              m.chain_map, "ok")
        check("sequence identity is exactly 1.0", m.seq_identity, 1.0)
        check("binder length", m.binder_len, 11)
        check("the prosapia-syntax designed positions",
              col.render_positions(m.design_positions), "4:6,10")

        # The mmCIF writer re-pooled the label (internal "C1" -> written "C" would
        # agree; here the written label is NOT recoverable from output_chain). The
        # mapping still holds -- it rests on chain POSITION, which is how
        # source_context is built -- but the reader is told it was not confirmed.
        npz_renamed = write_npz(
            tmpdir / "design_renamed.npz", 4, 27,
            ["/t/t.cif", "/t/t.cif", "/opt/bg/sc/7eow.cif"], ["A", "B", "Q0"],
        )
        m_renamed = col.map_binder(design, npz_renamed, refold, nb_seq, 4, STEMS)
        check("an unconfirmable output_chain is reported, not hidden",
              m_renamed.chain_map, "ok_renamed")
        check("and the scaffold is still read positionally",
              m_renamed.scaffold, "7eow.cif")
        check("and it was checked against the run's scaffold list",
              m_renamed.scaffold_src, "in_run_scaffold_list")

        # No source_context -> the scaffold is NA, and that is SAID, not guessed.
        bare = tmpdir / "bare.npz"
        np.savez(bare, design_mask=np.array([1.0] * 4 + [0.0] * 27, dtype=np.float32))
        m_bare = col.map_binder(design, bare, refold, nb_seq, 4, STEMS)
        check("a missing source_context is named", m_bare.chain_map, "no_source_context")
        check("and so is the scaffold's provenance",
              m_bare.scaffold_src, "no_source_context")
        check("and the scaffold itself is NA, not a guess",
              pd.isna(m_bare.scaffold), True)

        print("\n-- the mapping failures that must NEVER return a number ------------")
        wrong_npz = write_npz(
            tmpdir / "wrong.npz", 7, 24, ["/t/t.cif"] * 3, ["A", "B", "C"]
        )
        raises("npz design_mask count != B-factor count",
               lambda: col.map_binder(design, wrong_npz, refold, nb_seq, 7, STEMS),
               "designed-residue count disagrees")
        raises("CSV num_design != B-factor count",
               lambda: col.map_binder(design, npz, refold, nb_seq, 9, STEMS),
               "designed-residue count disagrees")
        raises("a sequence that is not the chain's",
               lambda: col.map_binder(design, npz, refold, "QVQLVESGGGW", 4, STEMS),
               "identity")
        raises("a colon in designed_chain_sequence (more than one designed chain)",
               lambda: col.map_binder(design, npz, refold, f"{nb_seq}:AAAA", 4, STEMS),
               "single designed nanobody chain")

        two_designed = write_pdb(
            tmpdir / "two.pdb",
            {
                "A": ("ACDEFGHIKL", [0.0] * 9 + [1.0]),
                "C": (nb_seq, designed_flags),
            },
        )
        raises("two designed chains",
               lambda: col.map_binder(two_designed, npz, refold, nb_seq, 4, STEMS),
               "exactly one designed chain")

        none_designed = write_pdb(
            tmpdir / "none.pdb", {"C": (nb_seq, [0.0] * 11)}
        )
        raises("no designed chain at all",
               lambda: col.map_binder(none_designed, npz, refold, nb_seq, 4, STEMS),
               "exactly one designed chain")

        short_refold = write_pdb(
            tmpdir / "short.pdb",
            {
                "A": ("ACDEFGHIKL", [88.0] * 10),
                "B": ("MNPQRSTVWY", [91.0] * 10),
                "Z": (nb_seq[:8], [77.0] * 8),
            },
        )
        raises("the refold's binder chain is a different length",
               lambda: col.map_binder(design, npz, short_refold, nb_seq, 4, STEMS),
               "cannot be carried across")

        fewer_chains = write_pdb(
            tmpdir / "fewer.pdb", {"A": ("ACDEFGHIKL", [88.0] * 10)}
        )
        raises("the refold has fewer chains than the design",
               lambda: col.map_binder(design, npz, fewer_chains, nb_seq, 4, STEMS),
               "chain index")


        print("\n-- BLOCKER 1: the scaffold name is never a positional guess ---------")
        # The entry at the binder's position names a file that is NOT one of the
        # scaffolds this run offered -> NA, with a verdict that is not a success word.
        npz_foreign = write_npz(
            tmpdir / "foreign.npz", 4, 27,
            ["/t.cif", "/t.cif", "/somewhere/else/8z8v.cif"], ["A", "B", "D"],
        )
        m_foreign = col.map_binder(design, npz_foreign, refold, nb_seq, 4, STEMS)
        check("a scaffold outside the run's list is refused",
              m_foreign.scaffold_src, "not_in_run_scaffold_list")
        check("and its name is NA, never emitted",
              pd.isna(m_foreign.scaffold), True)

        # The catastrophic ordering failure: source_context ordered by input entity, so
        # the binder's position resolves to the TARGET file. Caught with no allow-list
        # at all, because the binder's source would equal a non-binder chain's.
        npz_target = write_npz(
            tmpdir / "istarget.npz", 4, 27,
            ["/t.cif", "/t.cif", "/t.cif"], ["A", "B", "D"],
        )
        m_target = col.map_binder(design, npz_target, refold, nb_seq, 4, STEMS)
        check("the target's own file is never reported as a scaffold",
              m_target.scaffold_src, "is_a_target_chain_source")
        check("and it is NA", pd.isna(m_target.scaffold), True)
        check("the same check fires with NO allow-list (a --spec run)",
              col.map_binder(design, npz_target, refold, nb_seq, 4, None).scaffold_src,
              "is_a_target_chain_source")
        check("--spec says the allow-list did not run, not that it passed",
              col.map_binder(design, npz_foreign, refold, nb_seq, 4, None).scaffold_src,
              "distinct_from_target")

        # "output_chain missing" must not masquerade as "output_chain renamed".
        no_oc = tmpdir / "nooc.npz"
        np.savez(
            no_oc,
            design_mask=np.array([1.0] * 4 + [0.0] * 27, dtype=np.float32),
            source_context=np.asarray(json.dumps({
                "version": 1,
                "chains": [{"source": s} for s in
                           ["/t.cif", "/t.cif", "/opt/bg/sc/gontivimab.cif"]],
            })),
        )
        raises("a missing output_chain CANNOT happen upstream, so it is refused",
               lambda: col.map_binder(design, no_oc, refold, nb_seq, 4, STEMS),
               "no 'output_chain'")
        check("and only the two ok* verdicts are trusted",
              col.CHAIN_MAP_OK, ("ok", "ok_renamed"))

        print("\n-- item 3: design -> refold transfer is residue-for-residue --------")
        # Same LENGTH, one substitution. A length-only check passes this; the designed
        # positions would then be carried onto a chain that is not the same chain.
        shifted = write_pdb(
            tmpdir / "shifted.pdb",
            {
                "A": ("ACDEFGHIKL", [88.0] * 10),
                "B": ("MNPQRSTVWY", [91.0] * 10),
                "Z": ("QVQLVESGGGW", [77.0] * 11),
            },
        )
        raises("a length-preserving substitution is refused",
               lambda: col.map_binder(design, npz, shifted, nb_seq, 4, STEMS),
               "cannot be carried across")

        print("\n-- item 8: the third witness is not optional -----------------------")
        raises("an absent num_design is an error, not a two-witness pass",
               lambda: col.map_binder(design, npz, refold, nb_seq, None, STEMS),
               "two witnesses instead of three")
        raises("a NaN num_design is too",
               lambda: col.map_binder(design, npz, refold, nb_seq, float("nan"), STEMS),
               "two witnesses instead of three")

        print("\n-- item 11: the pLDDT column is never read as a design mask --------")
        check("a design mmCIF carries the mask",
              col.read_chains(design, design_mask=True)[2].designed is None, False)
        check("the refolded complex does NOT, and says None rather than 'all designed'",
              col.read_chains(refold, design_mask=False)[2].designed, None)

        print("\n-- item 13: the measured identity survives onto the error row ------")
        try:
            col.map_binder(design, npz, refold, "QVQLVESGGGW", 4, STEMS)
        except col.MappingError as e:
            check("the identity that condemned the row is carried in .data",
                  round(e.data["mapping_seq_identity"], 4), round(10 / 11, 4))
        else:
            check("the identity that condemned the row is carried in .data", None, "raise")

        print("\n-- BLOCKER 2: the epitope frame is the DECLARED sequence ------------")
        # SEQRES declares 5 residues; the 3rd has no coordinates. The 3rd MODELLED
        # residue is GLU at label_seq 4 -- counting modelled residues would call
        # position 3 'GLU' when BoltzGen binds the (unmodelled) ASP.
        gapped = tmpdir / "gapped.pdb"
        gapped.write_text(
            "SEQRES   1 A    5  ALA CYS ASP GLU PHE\n"
            "ATOM      1  CA  ALA A  11       0.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      2  CA  CYS A  12       3.800   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  CA  GLU A  14       7.600   0.000   0.000  1.00  0.00           C\n"
            "ATOM      4  CA  PHE A  15      11.400   0.000   0.000  1.00  0.00           C\n"
            "TER\nEND\n"
        )
        gframes = run.label_seq_residues(gapped).chains
        table = gframes["A"].residues
        check("the frame is the file's declared sequence",
              gframes["A"].res_frame, "declared_sequence")
        check("the declared sequence keeps its 5 positions", len(table), 5)
        check("position 3 is the DECLARED residue, not the 3rd modelled one",
              table[3].name, "ASP")
        check("and it is marked unmodelled", table[3].modelled, False)
        check("position 4 is the 3rd modelled residue", table[4].name, "GLU")
        check("an ordinal count would have said GLU at 3 -- it does not",
              table[3].name == "GLU", False)
        check("modelled positions resolve normally",
              run.resolve_binding_resnames(
                  {"A": [4, 5]}, TF(gapped), gapped, "--binding")[0],
              ["A4=GLU", "A5=PHE"])
        raises("a DECLARED but UNMODELLED epitope residue is an error",
               lambda: run.resolve_binding_resnames(
                   {"A": [3]}, TF(gapped), gapped, "--binding"),
               "NO COORDINATES")
        raises("a position past the declared sequence is an error",
               lambda: run.resolve_binding_resnames(
                   {"A": [6]}, TF(gapped), gapped, "--binding"),
               "positions 1..5")

        print("\n-- item 6: an epitope chain must be in --target-chains --------------")
        class TwoChain:
            target_chains = "A"
            binding = "B3"
            not_binding = None
            scaffolds = ",".join(run.DEFAULT_SCAFFOLDS)
            scaffold_dir = "/opt/bg/sc"
        raises("an epitope on an excluded chain is refused at submit time",
               lambda: run.build_spec("d", target, TwoChain(), lambda n, c: None, {}),
               "not in --target-chains")
        raises("and a chain that is not in the file at all is refused too",
               lambda: run.build_spec(
                   "d", target, type("X", (TwoChain,), {"binding": "Z3"})(),
                   lambda n, c: None, {}),
               "not author chains of this file")

        class InChain(TwoChain):
            target_chains = "A,B"
            binding = "A3"
        res_ok = run.build_spec("d", target, InChain(), lambda n, c: None, {})
        check("an epitope on an included chain is accepted",
              res_ok.binding_resnames, ["A3=ASP"])
        check("and the frame travels with it",
              res_ok.res_frame, "synthesized_from_model")
        check("and so does the chain frame (.pdb -> the author letter)",
              res_ok.chain_frame, "author_letter")
        check("with both sides of the id translation recorded",
              res_ok.chain_id_map, "A>A,B>B")

    print("\n-- the forced step/config contract ---------------------------------")
    check("the step list omits esmfold2_scoring AND filtering",
          run.STEPS, ("design", "inverse_folding", "folding", "analysis"))
    check("analysis is forced to skip esmfold2 metrics",
          run.FORCED_CONFIG, (("analysis", "esmfold2_metrics=false"),))
    check("the protocol is nanobody-anything", run.PROTOCOL, "nanobody-anything")

    class Flags:
        num_designs = 50
        diffusion_batch_size = None
        num_workers = 1
        config_set: list = []
        extra_flags = ""

    flags = run.cli_flags(Flags())
    check("--num_designs is passed", "--num_designs 50" in flags, True)
    check("the esmfold2 opt-out is always emitted",
          "--config analysis esmfold2_metrics=false" in flags, True)
    check("--steps never comes from the flags (it is in the .sh)",
          "--steps" in flags, False)
    check("--protocol never comes from the flags either",
          "--protocol" in flags, False)
    check("--budget is never emitted (filtering does not run)",
          "--budget" in flags, False)
    raises("a malformed --config-set is refused at submit time",
           lambda: run.cli_flags(type("F", (Flags,), {"config_set": ["analysis"]})()),
           "<step> <key>=<value>")

    sh = (HERE / "boltzgen.sh").read_text()
    check("the .sh pins the protocol", "--protocol nanobody-anything" in sh, True)
    check("the .sh pins the step list",
          "--steps design inverse_folding folding analysis" in sh, True)
    check("the .sh never mentions filtering",
          " filtering" in sh.split("--steps")[1].split("\n")[0], False)


    print("\n-- item 5: the escape hatches cannot redefine the tool --------------")
    for bad, why in (
        ("--protocol protein-anything", "changes WHAT is designed"),
        ("--steps design inverse_folding folding analysis filtering", "reinstates filtering"),
        ("--output /tmp/elsewhere", "moves the output the collector reads"),
        ("--num_designs 10000", "shadows --num-designs"),
        ("--cache /tmp/nowhere", "moves the weight cache off the Volume"),
        ("--config analysis esmfold2_metrics=true", "re-enables esmfold2"),
        ("--use_kernels false --protocol antibody-anything", "hidden among good tokens"),
        ("--protocol=protein-anything", "the = form"),
    ):
        raises(f"--extra-flags refuses {bad.split()[0]} ({why})",
               lambda b=bad: run.check_extra_flags(b), "may not contain")
    check("a harmless --extra-flags still passes through",
          run.check_extra_flags(" --use_kernels false --overwrite "),
          "--use_kernels false --overwrite")
    check("an empty --extra-flags is empty", run.check_extra_flags("   "), "")
    raises("--config-set may not reinstate esmfold2_metrics",
           lambda: run.check_config_set("analysis esmfold2_metrics=true"),
           "must stay there")
    raises("not even spelled oddly",
           lambda: run.check_config_set("analysis esmfold2_metrics=TRUE"),
           "must stay there")
    check("but it may re-state the forced value",
          run.check_config_set("analysis esmfold2_metrics=false"),
          "analysis esmfold2_metrics=false")
    check("and any other key is fine",
          run.check_config_set("analysis liability_modality=antibody"),
          "analysis liability_modality=antibody")

    print("\n-- the collector's error rows (a bad file never sinks a collect) ----")
    with tempfile.TemporaryDirectory() as tmp2:
        tmpdir2 = Path(tmp2)
        run_dir = tmpdir2 / "run"
        out_dir = run_dir / "boltzgen"
        out_dir.mkdir(parents=True)
        group = "g"
        design_dir = out_dir / "runs" / group / col.DESIGN_DIRNAME
        (design_dir / col.REFOLD_DIRNAME).mkdir(parents=True)

        NB2 = "QVQLVESGGGL"
        FL2 = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        chains_design = {"A": ("ACDEFGHIKL", [0.0] * 10), "D": (NB2, FL2)}
        chains_refold = {"A": ("ACDEFGHIKL", [88.0] * 10), "D": (NB2, [70.0] * 11)}

        ids = ["g_0", "g_1", "g_2", "g_3"]
        for did in ids:
            write_cif(design_dir / f"{did}.cif", chains_design)
            write_cif(design_dir / col.REFOLD_DIRNAME / f"{did}.cif", chains_refold)
            write_npz(design_dir / f"{did}.npz", 4, 17,
                      ["/t.cif", "/opt/bg/sc/7eow.cif"], ["A", "D"])
        # g_1: a CORRUPT .npz with a real ZIP magic number, so numpy gets as far as
        # zipfile and raises zipfile.BadZipFile -- upstream issue #223 exactly, and a
        # subclass of Exception, not ValueError. If this escapes the per-design try,
        # the WHOLE collect of a 12-hour run is lost.
        (design_dir / "g_1.npz").write_bytes(b"PK\x03\x04" + b"\x00" * 64)
        # g_2: a refolded mmCIF gemmi cannot parse -> RuntimeError, not ValueError.
        (design_dir / col.REFOLD_DIRNAME / "g_2.cif").write_text("garbage\n")
        # g_3: no files at all.
        for suffix in ("cif", "npz"):
            (design_dir / f"g_3.{suffix}").unlink()
        (design_dir / col.REFOLD_DIRNAME / "g_3.cif").unlink()

        rows_csv = [
            {"id": did, "file_name": f"{did}.cif", "designed_chain_sequence": NB2,
             "designed_sequence": "VESG", "num_design": 4, "design_ptm": 0.8}
            for did in ids
        ]
        # a duplicated id: two designs collapsing into one row, last-wins
        rows_csv.append(dict(rows_csv[0], design_ptm=0.1))
        pd.DataFrame(rows_csv).to_csv(design_dir / col.METRICS_CSV, index=False)

        def collect(record, meta: dict | None = None, spec: dict | None = None) -> dict:
            specs = out_dir / "specs"
            specs.mkdir(exist_ok=True)
            spec_path = specs / f"{group}.yaml"
            if spec is None:
                spec_path.unlink(missing_ok=True)
            else:
                spec_path.write_text(yaml.safe_dump(spec, sort_keys=False))
            (out_dir / ".meta.json").write_text(json.dumps(
                meta if meta is not None else {"runs_root": str(out_dir / "runs")}))
            rec_path = design_dir.parent / ".submit.json"
            if record is None:
                rec_path.unlink(missing_ok=True)
            elif isinstance(record, str):
                rec_path.write_text(record)  # deliberately damaged
            else:
                rec_path.write_text(json.dumps(record))
            cp = ArgumentParser()
            cp.add_argument("--force", action="store_true")
            col.add_collect_boltzgen_args(cp)
            cargs = cp.parse_args([])
            cargs.run_dir = run_dir
            cargs.table = "t0"
            cctx = CollectCtx(
                df=pd.DataFrame(), args=cargs, table_name="t0", out_dir=out_dir,
                status_col="boltzgen_status", path_col="boltzgen_path",
                parent_table=None, parent_df=pd.DataFrame(),
                lookup=lambda n, c: None, creates_table=True,
                default_input_column="pdb_path",
            )
            got = list(col.collect_boltzgen(cctx)(
                DesignCtx(name=group, out_dir=out_dir, lookup=cctx.lookup)))
            return {c.name: c for c in got}

        base_record = {
            "spec_sha256": "", "scaffolds": ["/opt/bg/sc/7eow.yaml"],
            "claims": [], "binding_resnames": "", "not_binding_resnames": "",
            "res_frame": "", "chain_frame": "", "chain_id_map": "",
        }
        got = collect(base_record)
        check("every design produced a row, none lost to a bad file", len(got), 4)
        check("the corrupt .npz is an error row, not an exception",
              got["g_1"].status.startswith("error:"), True)
        check("and it names BadZipFile specifically (issue #223)",
              "BadZipFile" in got["g_1"].status, True)
        check("the unparseable refold mmCIF is an error row too",
              got["g_2"].status.startswith("error:"), True)
        check("the missing-files row is an error row",
              got["g_3"].status.startswith("error:"), True)
        check("a missing refold leaves cif_path EMPTY, not a path that is not there",
              got["g_3"].data["cif_path"], "")
        check("the duplicated id is an error, never a silent last-wins collapse",
              "appears more than once" in got["g_0"].status, True)
        check("a duplicated id yields ONE row, not one per occurrence",
              sum(1 for k in got if k == "g_0"), 1)
        check("every error row has an empty path",
              [k for k, v in got.items() if v.status != "OK" and str(v.path)], [])

        # BLOCKER 1/4: a claim with no resnames beside it is unattested.
        got2 = collect(dict(base_record, claims=["binding"]))
        check("a 'binding' claim with no resnames makes every design an error",
              all(v.status.startswith("error:") for v in got2.values()), True)
        check("and says which claim", "'binding'" in got2["g_2"].status, True)

        # BLOCKER 4: a --not-binding-ONLY run is legitimate and must collect.
        got2b = collect(dict(base_record, claims=["not_binding"],
                             not_binding_resnames="B5=SER",
                             res_frame="declared_sequence",
                             chain_frame="label_asym_id", chain_id_map="B>D"))
        check("a --not-binding-only run is NOT an error",
              got2b["g_0"].status.startswith("error: '"), True)  # dup-id only
        check("...and its other designs collect",
              got2b["g_2"].status.startswith("error: unreadable"), True)
        check("the anti-epitope lands in its OWN column",
              got2b["g_2"].data["not_binding_resnames"], "B5=SER")
        # item 5: the identifier witnesses belong to BOTH claims.
        check("with the residue frame beside it",
              got2b["g_2"].data["binding_frame"], "declared_sequence")
        check("and the chain frame",
              got2b["g_2"].data["binding_chain_frame"], "label_asym_id")
        check("and the chain-id translation",
              got2b["g_2"].data["chain_id_map"], "B>D")
        check("and NO binding_resnames, because none was asked for",
              "binding_resnames" in got2b["g_2"].data, False)

        # BLOCKER 1: a MISSING per-group record is an error, not a silent NA.
        got2c = collect(None, meta={"runs_root": str(out_dir / "runs")})
        check("a missing .submit.json makes every design an error",
              all(v.status.startswith("error:") for v in got2c.values()), True)
        check("and leads with the routine cause, not a 12-hour re-run",
              "-l/--dir-label" in got2c["g_2"].status, True)
        check("and says WHICH state the record was in", "absent" in got2c["g_2"].status, True)
        # A damaged record is distinguished from an absent one.
        got2d = collect("{truncated")
        check("a truncated record errors too",
              all(v.status.startswith("error:") for v in got2d.values()), True)
        check("and says it was unreadable, not absent",
              "unreadable" in got2d["g_2"].status, True)

        # An error row carries NONE of the attestation it is complaining about.
        attest = ("binding_resnames", "not_binding_resnames", "binding_frame",
                  "binding_chain_frame", "chain_id_map")
        # A record whose spec_sha256 does not match specs/<group>.yaml: the realistic
        # "this attestation came from a different run" case, WITH resnames populated.
        stale = collect(
            dict(base_record, claims=["binding"], binding_resnames="A3=ASP",
                 res_frame="declared_sequence", chain_frame="label_asym_id",
                 chain_id_map="A>C", spec_sha256="deadbeef" * 8),
            spec={"entities": [{"file": {"path": "/t.cif", "binding_types": [
                {"chain": {"id": "A", "binding": "3"}}]}}]},
        )
        check("a stale record makes every design an error",
              all(v.status.startswith("error:") for v in stale.values()), True)
        check("and says the attestation does not describe these designs",
              "does not describe these designs" in stale["g_2"].status, True)
        check("and NOT ONE attestation column rides along on the error row",
              [k for k in attest if k in stale["g_2"].data], [])

        got3 = collect(dict(base_record, claims=["binding"],
                            binding_resnames="A3=ASP",
                            res_frame="declared_sequence", chain_id_map="A>C"))
        check("with the record present, the good designs collect",
              got3["g_0"].data["binding_resnames"], "A3=ASP")
        check("and the frame travels into the table",
              got3["g_0"].data["binding_frame"], "declared_sequence")
        check("and both sides of the chain-id translation do too",
              got3["g_0"].data["chain_id_map"], "A>C")


    print("\n-- BLOCKERS 1-3: the attestation lives WITH the outputs -------------")
    with tempfile.TemporaryDirectory() as tmp4:
        out5 = Path(tmp4) / "boltzgen"
        out5.mkdir(parents=True)
        t_a = write_pdb(Path(tmp4) / "ta.pdb", {"A": ("ACDEFGHIKL", [0.0] * 10)})
        t_b = write_pdb(Path(tmp4) / "tb.pdb", {"A": ("MNPQRSTVWY", [0.0] * 10)})

        from functools import partial
        from prosapia.core import ManifestCtx
        from prosapia.core.base_run import write_run_meta

        def submit(target: Path, extra: list[str]) -> list:
            ap = ArgumentParser()
            ap.add_argument("--table", default=None)
            ap.add_argument("--input-column", default="pdb_path")
            ap.add_argument("--force", action="store_true")
            ap.add_argument("--executor", default="modal")
            run.add_run_boltzgen_args(ap)
            a = ap.parse_args(["--target-structure", str(target)] + extra)
            a.run_dir = Path(tmp4)
            return run.build_boltzgen_manifest(ManifestCtx(
                df=pd.DataFrame(), args=a, out_dir=out5,
                lookup=lambda n, c: None,
                write_meta=partial(write_run_meta, out5),
            ))

        submit(t_a, ["--binding", "A3"])
        submit(t_b, ["--binding", "A5"])
        runs = out5 / "runs"
        rec_a, st_a = run.read_submit_record(runs / "ta_bg")
        rec_b, st_b = run.read_submit_record(runs / "tb_bg")
        check("both records are present", (st_a, st_b), ("present", "present"))
        check("each group gets its OWN record, beside its own outputs",
              (rec_a["binding_resnames"], rec_b["binding_resnames"]),
              ("A3=ASP", "A5=ARG"))
        check("so a second submit cannot overwrite the first's attestation",
              run.read_submit_record(runs / "ta_bg")[0]["binding_resnames"], "A3=ASP")
        check("and each carries its own scaffold allow-list",
              bool(rec_a["scaffolds"]) and bool(rec_b["scaffolds"]), True)
        check("root_designs is rebuilt FROM DISK, not from the clobbered sidecar",
              sorted(json.loads((out5 / ".meta.json").read_text())["root_designs"]),
              ["ta_bg", "tb_bg"])

        # BLOCKER 1, the core case: resubmit a group WITHOUT --binding. The old shared
        # sidecar kept the previous epitope and certified unsteered designs with it.
        submit(t_a, [])
        check("a resubmit without --binding CLEARS the previous epitope",
              run.read_submit_record(runs / "ta_bg")[0]["binding_resnames"], "")
        check("and records no claim",
              run.read_submit_record(runs / "ta_bg")[0]["claims"], [])

        # BLOCKER 1, mirror case: a changed spec over EXISTING designs is refused,
        # because BoltzGen is never given --overwrite and the old designs would survive.
        (runs / "tb_bg" / "intermediate_designs_inverse_folded").mkdir(parents=True)
        raises("a changed spec over EXISTING designs is refused, not topped up",
               lambda: submit(t_b, ["--binding", "A7"]),
               "already holds designs made from a DIFFERENT spec")
        check("and the old record is left intact, not half-replaced",
              run.read_submit_record(runs / "tb_bg")[0]["binding_resnames"], "A5=ARG")
        # The same change into a group with no designs yet is fine.
        rows2 = submit(t_a, ["--binding", "A7"])
        check("while a group with no designs yet accepts the new spec", len(rows2), 1)
        check("and its attestation is the NEW one",
              run.read_submit_record(runs / "ta_bg")[0]["binding_resnames"], "A7=HIS")

        # BLOCKER 5: the refusal must require a readable MATCH, not merely forbid a
        # confirmed mismatch. Deleting the record -- which the old message's "delete
        # {group_dir}" advice invited -- used to turn the check off entirely.
        (runs / "tb_bg" / ".submit.json").unlink()
        raises("deleting the record does NOT re-enable a stale top-up",
               lambda: submit(t_b, ["--binding", "A7"]),
               "there is no evidence they were made from this run's spec")
        raises("and the advice names the WHOLE directory, not the record",
               lambda: submit(t_b, ["--binding", "A7"]),
               "delete the WHOLE directory")
        raises("and says deleting the record alone is the bug",
               lambda: submit(t_b, ["--binding", "A7"]),
               "leaves the old designs in place")
        (runs / "tb_bg" / ".submit.json").write_text("{truncated")
        raises("a truncated record is refused the same way",
               lambda: submit(t_b, ["--binding", "A7"]), "unreadable")
        (runs / "tb_bg" / ".submit.json").write_text(json.dumps({"claims": []}))
        raises("so is a record with no spec_sha256",
               lambda: submit(t_b, ["--binding", "A7"]), "missing a spec_sha256")
        check("and read_submit_record names which state it was",
              run.read_submit_record(runs / "tb_bg")[1], "present")
        check("absent is its own state",
              run.read_submit_record(runs / "nope_bg"), (None, "absent"))

        # BLOCKER 2: the framework rewrites .meta.json BEFORE the builder runs and
        # swallows a corrupt read, so a damaged sidecar cannot take the attestation
        # with it any more.
        (out5 / ".meta.json").write_text("{corrupt")
        submit(t_a, ["--binding", "A7"])
        check("a corrupt run sidecar does not touch any group's record",
              run.read_submit_record(runs / "ta_bg")[0]["binding_resnames"], "A7=HIS")
        check("nor the group it did not resubmit",
              run.read_submit_record(runs / "ta_bg")[1], "present")
        check("and root_designs is still complete, rebuilt from disk",
              sorted(json.loads((out5 / ".meta.json").read_text())["root_designs"]),
              ["ta_bg", "tb_bg"])

    print("\n-- item 7: chain ids are TRANSLATED, not assumed --------------------")
    with tempfile.TemporaryDirectory() as tmp5:
        d5 = Path(tmp5)
        # A .pdb: patch_chain_names gives the author letter back, so the map is identity.
        pdb5 = write_pdb(d5 / "p.pdb", {"A": ("ACDEFGHIKL", [0.0] * 10)})
        f5 = run.label_seq_residues(pdb5)
        ids, frame = run.resolve_chain_ids({"A"}, f5, pdb5)
        check("a .pdb resolves to the author letter",
              (ids, frame), ({"A": "A"}, "author_letter"))
        raises("and an author letter that is not in the file is refused",
               lambda: run.resolve_chain_ids({"Q"}, f5, pdb5), "not author chains")

        # An mmCIF whose author letter and label_asym_id DIVERGE -- the 1TUP case,
        # self-contained: author 'H', label_asym_id 'A'.
        cif5 = d5 / "diverge.cif"
        cif5.write_text(DIVERGING_CIF)
        f6 = run.label_seq_residues(cif5)
        check("the table is keyed by the AUTHOR letter a viewer shows",
              sorted(f6.chains), ["H"])
        check("and carries the label_asym_id beside it", f6.chains["H"].label_asym, "A")
        ids6, frame6 = run.resolve_chain_ids({"H"}, f6, cif5)
        check("an mmCIF translates author -> label_asym_id",
              (ids6, frame6), ({"H": "A"}, "label_asym_id"))
        # The whole point: 'A' is what the FILE calls the label chain, but the caller
        # must write 'H'. Writing 'A' is refused rather than silently resolved.
        raises("a label_asym_id passed as if it were an author letter is refused",
               lambda: run.resolve_chain_ids({"A"}, f6, cif5), "not author chains")
        check("residue names resolve against the author-keyed chain",
              run.resolve_binding_resnames({"H": [1, 3]}, f6, cif5, "--binding")[0],
              ["H1=ALA", "H3=ASP"])

        # End to end: the YAML must carry the TRANSLATED id, not the author letter.
        class DivArgs:
            target_chains = "H"
            binding = "H3"
            not_binding = None
            scaffolds = ",".join(run.DEFAULT_SCAFFOLDS)
            scaffold_dir = "/opt/bg/sc"
        res_div = run.build_spec("d", cif5, DivArgs(), lambda n, c: None, {})
        target_entity = res_div.spec["entities"][0]["file"]
        check("the include block carries the label_asym_id",
              target_entity["include"], [{"chain": {"id": "A"}}])
        check("and so does binding_types",
              target_entity["binding_types"], [{"chain": {"id": "A", "binding": "3"}}])
        check("while the ATTESTATION stays in the caller's own letters",
              res_div.binding_resnames, ["H3=ASP"])
        check("and both sides of the translation are recorded",
              res_div.chain_id_map, "H>A")
        check("with the frame that produced it",
              res_div.chain_frame, "label_asym_id")

    print("\n-- claims come from the SPEC, and both walkers agree ----------------")
    spec_both = {"entities": [{"file": {"path": "/t.cif", "binding_types": [
        {"chain": {"id": "A", "binding": "3,5"}},
        {"chain": {"id": "B", "not_binding": "7"}},
    ]}}, {"file": {"path": ["/s/7eow.yaml"]}}]}
    spec_none = {"entities": [{"file": {"path": "/t.cif"}}]}
    check("a spec with both steering kinds reports both",
          run.claims_in_spec(spec_both), ["binding", "not_binding"])
    check("a spec with neither reports none", run.claims_in_spec(spec_none), [])
    with tempfile.TemporaryDirectory() as tmp8:
        d8 = Path(tmp8)
        for name, spec in (("both", spec_both), ("none", spec_none)):
            y = d8 / f"{name}.yaml"
            y.write_text(yaml.safe_dump(spec, sort_keys=False))
            # The builder walks a dict, the collector walks the YAML on disk. They are
            # separate implementations by design (the two phases do not import each
            # other); this is the guard that keeps them in step.
            check(f"the run and collect walkers agree on '{name}'",
                  set(run.claims_in_spec(spec)), col.claims_from_spec(y))

    # A --spec run's claims come from the YAML the user wrote, not from flags that
    # never saw it. Recording [] there made every design error at collect, blaming a
    # stale record for a bookkeeping gap in the builder.
    with tempfile.TemporaryDirectory() as tmp9:
        d9 = Path(tmp9)
        out9 = d9 / "boltzgen"
        out9.mkdir(parents=True)
        user_spec = d9 / "mine.yaml"
        user_spec.write_text(yaml.safe_dump(spec_both, sort_keys=False))

        from functools import partial
        from prosapia.core import ManifestCtx
        from prosapia.core.base_run import write_run_meta

        ap = ArgumentParser()
        ap.add_argument("--table", default=None)
        ap.add_argument("--input-column", default="pdb_path")
        ap.add_argument("--force", action="store_true")
        ap.add_argument("--executor", default="modal")
        run.add_run_boltzgen_args(ap)
        a = ap.parse_args(["--spec", str(user_spec)])
        a.run_dir = d9
        run.build_boltzgen_manifest(ManifestCtx(
            df=pd.DataFrame(), args=a, out_dir=out9,
            lookup=lambda n, c: None, write_meta=partial(write_run_meta, out9)))
        rec, _ = run.read_submit_record(out9 / "runs" / "mine_bg")
        check("a --spec run records the claims ITS OWN YAML makes",
              rec["claims"], ["binding", "not_binding"])
        check("and no resnames, because the tool did not build that epitope",
              (rec["binding_resnames"], rec["not_binding_resnames"]), ("", ""))
        check("so the spec and the record AGREE -- no false disagreement",
              col.claims_from_spec(out9 / "specs" / "mine_bg.yaml"),
              set(rec["claims"]))

    print("\n-- items 8/9/15/6/7: the frame guards --------------------------------")
    with tempfile.TemporaryDirectory() as tmp7:
        d7 = Path(tmp7)
        # item 8: the format is gemmi's DETECTION, not a suffix test. A PDB named .ent
        # must still be classified Pdb, or the whole chain translation flips.
        ent = d7 / "pdb1abc.ent"
        ent.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      2  CA  CYS A   2       3.800   0.000   0.000  1.00  0.00           C\n"
            "TER\nEND\n"
        )
        f_ent = run.label_seq_residues(ent)
        check("a .ent file is detected as Pdb, not assumed mmCIF", f_ent.fmt, "Pdb")
        check("so it resolves to the author letter",
              run.resolve_chain_ids({"A"}, f_ent, ent)[1], "author_letter")
        # A proper .cif is detected and translated.
        cif7 = d7 / "t.cif"
        cif7.write_text(DIVERGING_CIF)
        f_cif = run.label_seq_residues(cif7)
        check("an mmCIF is detected as Mmcif", f_cif.fmt, "Mmcif")
        check("so it translates author -> label_asym_id",
              run.resolve_chain_ids({"H"}, f_cif, cif7), ({"H": "A"}, "label_asym_id"))
        # An extension neither gemmi nor we can classify is REFUSED. The old suffix
        # test silently classified `target.PDB1` as mmCIF and flipped the whole
        # translation; refusing is the safe outcome, not detection.
        odd = d7 / "target.PDB1"
        odd.write_text(DIVERGING_CIF)
        raises("an extension that cannot be classified is refused, not assumed mmCIF",
               lambda: run.label_seq_residues(odd), "cannot read this file")
        garbage = d7 / "x.txt"
        garbage.write_text("this is not a structure\n")
        raises("an unclassifiable file is refused, not guessed at",
               lambda: run.label_seq_residues(garbage), "cannot read this file")

        # item 9: chain ids are matched CASE-SENSITIVELY. A file with both A and a is
        # legal; folding one onto the other attests a residue from a chain nobody named.
        mixed = write_pdb(d7 / "mixed.pdb",
                          {"A": ("ACDEF", [0.0] * 5), "B": ("MNPQR", [0.0] * 5)})
        f_mixed = run.label_seq_residues(mixed)
        raises("a lowercase id is NOT silently folded to uppercase",
               lambda: run.resolve_chain_ids({"a"}, f_mixed, mixed),
               "not author chains")
        raises("and the message offers the case-sensitive near miss",
               lambda: run.resolve_chain_ids({"a"}, f_mixed, mixed),
               "case-sensitive")
        check("the residue parser no longer uppercases either",
              run.parse_residue_spec("a3,B5", "--binding"), {"a": [3], "B": [5]})

        # item 6: patch_chain_names keys on the first character of EVERY written
        # subchain, so a ligand sharing the polymer's first letter is unpredictable --
        # and it may be the POLYMER that gets re-lettered.
        lig = d7 / "lig.pdb"
        lig.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      2  CA  CYS A   2       3.800   0.000   0.000  1.00  0.00           C\n"
            "TER\n"
            "HETATM    3 ZN    ZN A 101       9.000   0.000   0.000  1.00  0.00          ZN\n"
            "END\n"
        )
        f_lig = run.label_seq_residues(lig)
        clashes = [k for k, v in f_lig.first_chars.items() if len(v) > 1]
        check("a ligand subchain IS visible in first_chars", bool(clashes), True)
        raises("and a .pdb whose subchains share a first character is refused",
               lambda: run.resolve_chain_ids({"A"}, f_lig, lig),
               "share a first character")
        # item 10: and the non-polymer subchain an explicit --target-chains drops is named
        check("the non-polymer subchain is recorded so the submit can name it",
              bool(f_lig.nonpolymer.get("A")), True)

    print("\n-- item 6: microheterogeneity is refused at submit time -------------")
    with tempfile.TemporaryDirectory() as tmp6:
        micro = Path(tmp6) / "micro.cif"
        # Two _entity_poly_seq rows sharing a num. Upstream APPENDS one array element
        # per row (mmcif.py:998) and then trips its own
        # `assert len(set(res_indices)) == len(res_indices)` (mmcif.py:803-806), so this
        # file is a hard parse failure inside a GPU task. Refusing at submit time costs
        # a submit instead -- and if that assert is ever relaxed, the silent off-by-N
        # returns and this guard is what stands between it and the table.
        micro.write_text(DIVERGING_CIF.replace("1 2 CYS n", "1 2 CYS n\n1 2 SER y"))
        raises("a target with microheterogeneity is refused before the GPU spend",
               lambda: run.label_seq_residues(micro), "microheterogeneity")
        raises("and the message says why the frame would shift",
               lambda: run.label_seq_residues(micro), "DIFFERENT index")

    print("\n-- the pin ---------------------------------------------------------")
    check("BOLTZGEN_REF is a full 40-char SHA", len(img.BOLTZGEN_REF), 40)
    check("and it is the SHA the recon pinned",
          img.BOLTZGEN_REF.startswith("a07fea1"), True)
    reqs = (HERE / "boltzgen_requirements.txt").read_text()
    check("torch is pinned to the last cu12 release", "torch==2.10.0" in reqs, True)
    pins = [ln for ln in reqs.splitlines() if ln.strip() and not ln.startswith("#")]
    check("no CUDA-13 wheel leaked into the closure",
          any("cu13" in ln for ln in pins), False)
    check("cuequivariance stays on cu12",
          "cuequivariance-ops-cu12==" in reqs, True)
    check("the closure is exact (every line is ==)",
          all("==" in ln for ln in pins), True)
    check("boltzgen itself is NOT in the requirements (installed --no-deps)",
          any(ln.startswith("boltzgen") for ln in reqs.splitlines()), False)

    print(f"\n{PASS} passed, {FAIL} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
