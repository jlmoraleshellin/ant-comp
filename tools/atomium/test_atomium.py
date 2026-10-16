#!/usr/bin/env python3
"""Self-test for the atomium tool: the _init_sequence skip, the trust columns,
the submit-time refusals and staging isolation.

WHAT THIS PROVES, AND WHAT IT CANNOT.

It proves the *contracts*, against frames this file writes:

  * the ``_init_sequence`` echo is skipped when its header is the REAL one,
    ``">_init_sequence, chains=D"`` -- the regression that matters most, because
    the old exact-equality test never fired and the echo was collected as a
    design row (``<parent>_a0``, status OK, carrying the native input sequence);
  * draws keep their documented names ``_a1``, ``_a2``, ... (the enumerate index
    over ALL entries, echo included), because row names are lineage keys;
  * each trust column reads 1 on a good row, reads 0 ONLY on a row whose status is
    an error, and reads NA (never 0, never a number) when the check could not be
    made at all;
  * every refusal raises at manifest-build time, before anything is queued: an
    out-of-range position, a group count that does not match the designed chains,
    both-flag mutual exclusion, a missing column, a missing per-row value, and a
    chain letter the design's own structure does not have;
  * two submits to one leaf stage into separate subgroup dirs with separate,
    content-addressed task files, and a changed signature for an already-staged
    design is refused rather than leaving two FASTAs for one design.

It CANNOT prove the premise or the mapping. Every FASTA and every jsonl here is
written by this file, so it agrees with this file by construction. These remain
UNKNOWN until the tool is re-run on real data:

  * whether a MULTI-CHAIN designed FASTA really joins its chains with ``/`` and in
    which order -- the probe that settled the single-chain case had one designed
    chain, so ``len_matches_parsed`` is deliberately order-independent and
    ``fixed_held`` is NA with more than one designed chain;
  * whether ``assigned_pdbs.jsonl`` / ``fixed_pdbs.jsonl`` carry exactly the shapes
    assumed here (``{design: [[designed],[fixed]]}`` and
    ``{design: {chain: [positions]}}``) in the real container, for the pinned
    AtomiUM commit -- they are ProteinMPNN's helper formats, read from a checkout
    this file never sees;
  * whether gemmi's polymer chain lengths (used at submit time for the bounds
    check) equal ``parse_multiple_chains.py``'s residue counts on a real input
    with hetero residues, altlocs or insertion codes. A disagreement makes the
    bounds check refuse a valid spec (loud) -- it cannot make a bad one pass
    silently, because the collector re-checks lengths against the task's own file;
  * whether ``--symmetry`` broadcast and ``--tied-positions`` behave as assumed in
    the real helper scripts.

Run:  uv run python tools/atomium/test_atomium.py
"""

import json
import sys
import tempfile
from argparse import ArgumentParser, Namespace
from pathlib import Path

import pandas as pd

from prosapia.core import CollectCtx, DesignCtx, ManifestCtx
from prosapia.core.tool_registry import _load as load_tool_spec

HERE = Path(__file__).resolve().parent

# Load the tool exactly as `sapia` discovers it, through the registry's own spec
# loader: the modules under test are the ones a real run would use.
TOOL = load_tool_spec(HERE / "spec.py")
run_mod = sys.modules[TOOL.build_manifest_fn.__module__]
collect_mod = sys.modules[TOOL.collect_fn.__module__]

FAILURES: list[str] = []
CHECKS = 0


def check(label: str, got, want) -> None:
    global CHECKS
    CHECKS += 1
    if got != want:
        FAILURES.append(f"{label}: got {got!r}, want {want!r}")


def check_true(label: str, value: bool) -> None:
    global CHECKS
    CHECKS += 1
    if not value:
        FAILURES.append(f"{label}: expected True")


def check_na(label: str, value) -> None:
    global CHECKS
    CHECKS += 1
    if not pd.isna(value):
        FAILURES.append(f"{label}: expected NA, got {value!r}")


def check_raises(label: str, fn, *needles: str) -> None:
    """The call must raise, and the message must contain every needle."""
    global CHECKS
    CHECKS += 1
    try:
        fn()
    except Exception as exc:  # noqa: BLE001 - we assert on the message
        message = str(exc)
        for needle in needles:
            if needle not in message:
                FAILURES.append(f"{label}: message {message!r} lacks {needle!r}")
        return
    FAILURES.append(f"{label}: expected a raise, got none")


# ---------------------------------------------------------------------------
# Fixtures: a fake out_dir shaped exactly like a real grp_ dir
# ---------------------------------------------------------------------------

# A 20-residue "chain D" and a 12-residue "chain A" (the target, never designed).
NATIVE_D = "QVQLVESGGGLVQAGGS"  # 17
NATIVE_A = "MKTAYIAKQRQI"  # 12
# The designed draw: same length, with positions 3 and 11 (1-based) preserved.
DESIGNED_1 = "QVQAAAAAAALAAAAAA"
DESIGNED_2 = "QVQEEEEEEELEEEEEE"
FIXED_POS = [3, 11]


def _pdb(path: Path, chains=(("A", 12), ("D", 17))) -> None:
    """A minimal multi-chain PDB; gemmi reads the chain names and residue counts."""
    lines = []
    serial = 1
    for chain, n in chains:
        for resnum in range(1, n + 1):
            lines.append(
                f"ATOM  {serial:5d}  CA  GLY {chain}{resnum:4d}    "
                f"{resnum * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C"
            )
            serial += 1
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def write_group(
    out_dir: Path,
    design: str = "d1",
    *,
    grp: str = "grp_aaaa",
    init: str | None = None,
    draws: tuple[str, ...] = (DESIGNED_1, DESIGNED_2),
    designed_chains: list[str] | None = ["D"],
    fixed: dict[str, list[int]] | None = None,
    parsed_chains: dict[str, str] | None = None,
    parsed_name: str | None = None,
    init_header: str = "_init_sequence, chains=D",
    staged_pdb: Path | None = None,
) -> Path:
    """Write one grp_ dir the way atomium.sh leaves it: the three helper jsonls
    plus seqs/<design>.fa. Every knob here exists to break exactly one assumption."""
    grp_dir = out_dir / grp
    (grp_dir / "seqs").mkdir(parents=True, exist_ok=True)
    (grp_dir / "inputs").mkdir(parents=True, exist_ok=True)
    if staged_pdb is not None:
        link = grp_dir / "inputs" / f"{design}.pdb"
        if not link.exists():
            link.symlink_to(staged_pdb)

    chains = parsed_chains if parsed_chains is not None else {"A": NATIVE_A, "D": NATIVE_D}
    entry = {"name": parsed_name or design}
    for c, seq in chains.items():
        entry[f"seq_chain_{c}"] = seq
    (grp_dir / "parsed_pdbs.jsonl").write_text(json.dumps(entry) + "\n")

    if designed_chains is not None:
        fixed_chains = [c for c in chains if c not in designed_chains]
        (grp_dir / "assigned_pdbs.jsonl").write_text(
            json.dumps({design: [designed_chains, fixed_chains]}) + "\n"
        )
    if fixed is not None:
        full = {c: fixed.get(c, []) for c in chains}
        (grp_dir / "fixed_pdbs.jsonl").write_text(json.dumps({design: full}) + "\n")

    if init is None:
        init = "/".join(chains[c] for c in (designed_chains or list(chains)))
    lines = [f">{init_header}", init]
    for k, seq in enumerate(draws):
        lines.append(f">gen_sample_{k}, temperature=0.1, seq_rec=0.43")
        lines.append(seq)
    (grp_dir / "seqs" / f"{design}.fa").write_text("\n".join(lines) + "\n")
    return grp_dir


def collect(out_dir: Path, design: str = "d1") -> list:
    """Run the real collector over a fake out_dir and return the Collected rows."""
    ctx = CollectCtx(
        df=pd.DataFrame(),
        args=Namespace(run_dir=out_dir.parent, table="table1", dir_label="", force=False),
        table_name="table1",
        out_dir=out_dir,
        status_col="atomium_status",
        path_col="atomium_path",
        parent_table="table0",
        parent_df=pd.DataFrame(index=pd.Index([design])),
        lookup=lambda name, column: None,
        creates_table=True,
        default_input_column="rfdiffusion3_path",
    )
    one = collect_mod.collect_atomium(ctx)
    return list(one(DesignCtx(design, out_dir, ctx.lookup)))


# ---------------------------------------------------------------------------
# Fixtures: manifest builder
# ---------------------------------------------------------------------------


def make_frame(tmp: Path, names=("d1", "d2"), **columns) -> pd.DataFrame:
    inputs = tmp / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    rows = []
    for name in names:
        pdb = inputs / f"{name}.pdb"
        _pdb(pdb)
        row: dict = {"name": name, "rfdiffusion3_path": str(pdb)}
        for col, values in columns.items():
            row[col] = values[name]
        rows.append(row)
    return pd.DataFrame(rows).set_index("name")


def make_args(run_dir: Path, **over) -> Namespace:
    base = dict(
        run_dir=run_dir,
        input_column="rfdiffusion3_path",
        force=False,
        designs_per_task=10,
        chains_to_design="",
        chains_to_design_column="",
        fixed_positions="",
        fixed_positions_column="",
        tied_positions="",
        symmetry=None,
        bias_aa="",
        set=[],
        model_noise="n05",
        num_seq_per_target=2,
        sampling_temp="0.1",
        seed=37,
    )
    base.update(over)
    return Namespace(**base)


def make_ctx(df: pd.DataFrame, args: Namespace, out_dir: Path) -> ManifestCtx:
    out_dir.mkdir(parents=True, exist_ok=True)
    return ManifestCtx(
        df=df,
        args=args,
        out_dir=out_dir,
        lookup=lambda name, column: None,
        write_meta=lambda **kw: None,
    )


def subgroups(out_dir: Path, rows: list[tuple[str, ...]]) -> list[dict]:
    """Every subgroup this submit staged, read the way atomium.sh reads it."""
    out: list[dict] = []
    for row in rows:
        for line in Path(row[0]).read_text().splitlines():
            if not line:
                continue
            grp, chains, fixedp, tie, tied, temp, extra = line.split("\t")
            grp_dir = Path(grp)
            out.append(
                {
                    "task": Path(row[0]).name,
                    "grp": grp_dir,
                    "chains": chains,
                    "fixed": fixedp,
                    "members": {
                        link.stem: str(link.readlink())
                        for link in sorted((grp_dir / "inputs").glob("*.pdb"))
                    },
                }
            )
    return out


# ---------------------------------------------------------------------------
# Collector tests
# ---------------------------------------------------------------------------


def test_init_sequence_is_skipped_with_the_real_header(tmp: Path) -> None:
    """THE regression: the real header carries ', chains=D' after the name."""
    out = tmp / "init_real"
    write_group(out, fixed={"D": FIXED_POS})
    rows = collect(out)
    check("real header: row count", len(rows), 2)
    check("real header: names", [r.name for r in rows], ["d1_a1", "d1_a2"])
    check_true(
        "real header: no row carries the native echo",
        all(r.data["sequence"] != NATIVE_D for r in rows),
    )
    check_true("real header: all OK", all(r.status == "OK" for r in rows))
    # And the predicate itself, directly, on every shape seen or plausible.
    check_true("is_init: real", collect_mod.is_init_header("_init_sequence, chains=D"))
    check_true("is_init: bare", collect_mod.is_init_header("_init_sequence"))
    check_true(
        "is_init: more trailing fields",
        collect_mod.is_init_header("_init_sequence, chains=D, foo=1"),
    )
    check_true(
        "is_init: padded",
        collect_mod.is_init_header(" _init_sequence , chains=D"),
    )
    check_true(
        "is_init: a generated header never matches",
        not collect_mod.is_init_header("gen_sample_0, temperature=0.1, seq_rec=0.43"),
    )
    check_true(
        "is_init: a design named like the echo never matches",
        not collect_mod.is_init_header("gen_sample_0, name=_init_sequence"),
    )


def test_bare_init_header_still_skipped(tmp: Path) -> None:
    """The old header shape keeps working, so an older out_dir re-collects the same."""
    out = tmp / "init_bare"
    write_group(out, init_header="_init_sequence", fixed={"D": FIXED_POS})
    rows = collect(out)
    check("bare header: row count", len(rows), 2)
    check("bare header: names", [r.name for r in rows], ["d1_a1", "d1_a2"])


def test_row_names_are_a1_a2(tmp: Path) -> None:
    """Draw naming is the enumerate index over ALL entries, echo included.

    Deliberate: row names are lineage keys and every existing atomium table (and
    the skill) documents _a1.. . Renumbering from _a0 would re-key them.
    """
    out = tmp / "names"
    write_group(out, draws=(DESIGNED_1, DESIGNED_2, DESIGNED_1))
    rows = collect(out)
    check("names: three draws", [r.name for r in rows], ["d1_a1", "d1_a2", "d1_a3"])
    check("names: iteration column", [r.data["iteration"] for r in rows], [1, 2, 3])


def test_trust_columns_on_a_good_row(tmp: Path) -> None:
    out = tmp / "good"
    write_group(out, fixed={"D": FIXED_POS})
    row = collect(out)[0]
    check("good: status", row.status, "OK")
    check("good: init_matches_parsed", row.data["init_matches_parsed"], 1)
    check("good: len_matches_parsed", row.data["len_matches_parsed"], 1)
    check("good: fixed_held", row.data["fixed_held"], 1)
    check("good: n_fixed", row.data["n_fixed"], 2)
    check("good: n_chains_designed", row.data["n_chains_designed"], 1)
    check("good: parsed_len", row.data["parsed_len"], len(NATIVE_D))
    check("good: seq_len", row.data["seq_len"], len(DESIGNED_1))
    check("good: fixed_positions", row.data["fixed_positions"], "D:3,11")
    check(
        "good: fixed_residues are chain-prefixed",
        row.data["fixed_residues"],
        "D:" + "".join(NATIVE_D[p - 1] for p in FIXED_POS),
    )
    check("good: init_order_matches", row.data["init_order_matches"], 1)
    check("good: sequence is the full chain", row.data["sequence"], DESIGNED_1)
    check("good: header fields", row.data["temperature"], 0.1)
    check("good: sample", row.data["sample"], 0)


def test_no_fixed_positions(tmp: Path) -> None:
    """Nothing pinned: n_fixed is a measured 0, fixed_held is NA (never a vacuous 1)."""
    out = tmp / "nofixed"
    write_group(out, fixed=None)
    row = collect(out)[0]
    check("nofixed: status", row.status, "OK")
    check("nofixed: n_fixed", row.data["n_fixed"], 0)
    check_na("nofixed: fixed_held", row.data["fixed_held"])
    check_na("nofixed: fixed_positions is NA, not ''", row.data["fixed_positions"])
    check_na("nofixed: fixed_residues is NA, not ''", row.data["fixed_residues"])
    check("nofixed: init_matches_parsed", row.data["init_matches_parsed"], 1)


def test_all_chains_designed_without_assigned_file(tmp: Path) -> None:
    """No --chains-to-design: every parsed chain is designed, so the FASTA is both."""
    out = tmp / "allchains"
    both = NATIVE_A + "/" + NATIVE_D
    write_group(
        out,
        designed_chains=None,
        init=both,
        draws=("A" * len(NATIVE_A) + "/" + "C" * len(NATIVE_D),),
        init_header="_init_sequence, chains=A,D",
    )
    row = collect(out)[0]
    check("allchains: status", row.status, "OK")
    check("allchains: n_chains_designed", row.data["n_chains_designed"], 2)
    check("allchains: parsed_len", row.data["parsed_len"], len(NATIVE_A) + len(NATIVE_D))
    check("allchains: len_matches_parsed", row.data["len_matches_parsed"], 1)
    check_na("allchains: fixed_held is NA with >1 designed chain", row.data["fixed_held"])


def test_fixed_position_not_held_is_an_error(tmp: Path) -> None:
    """A pinned residue that moved: error row, fixed_held 0 -- never a silent sequence."""
    out = tmp / "notheld"
    broken = list(DESIGNED_1)
    broken[FIXED_POS[0] - 1] = "W"  # was Q
    write_group(out, draws=("".join(broken),), fixed={"D": FIXED_POS})
    row = collect(out)[0]
    check_true("notheld: status is an error", row.status.startswith("error:"))
    check("notheld: fixed_held", row.data["fixed_held"], 0)
    check("notheld: len_matches_parsed still 1", row.data["len_matches_parsed"], 1)
    check_true("notheld: the message names the position", "3:" in row.status)
    check_true(
        "notheld: the sequence AtomiUM wrote is still recorded",
        row.data["sequence"] == "".join(broken),
    )


def test_truncated_sequence_is_an_error(tmp: Path) -> None:
    """The failure the skill used to claim was normal: a sequence short of the chain."""
    out = tmp / "short"
    short = "".join(c for i, c in enumerate(DESIGNED_1) if (i + 1) not in FIXED_POS)
    write_group(out, draws=(short,), fixed={"D": FIXED_POS})
    row = collect(out)[0]
    check_true("short: status is an error", row.status.startswith("error:"))
    check("short: len_matches_parsed", row.data["len_matches_parsed"], 0)
    check_na("short: fixed_held not guessed on a bad frame", row.data["fixed_held"])
    check("short: init_matches_parsed unaffected", row.data["init_matches_parsed"], 1)


def test_init_mismatch_is_an_error(tmp: Path) -> None:
    """If the echo is not the parsed chain, the row's frame is unproven -> error."""
    out = tmp / "initbad"
    write_group(out, init=NATIVE_D[:-1], fixed={"D": FIXED_POS})
    rows = collect(out)
    check_true("initbad: status is an error", rows[0].status.startswith("error:"))
    check("initbad: init_matches_parsed", rows[0].data["init_matches_parsed"], 0)
    check_true("initbad: rows still minted", len(rows) == 2)


def test_missing_init_entry_is_an_error(tmp: Path) -> None:
    out = tmp / "noinit"
    grp = write_group(out, fixed={"D": FIXED_POS})
    fa = grp / "seqs" / "d1.fa"
    kept = [ln for ln in fa.read_text().splitlines()][2:]
    fa.write_text("\n".join(kept) + "\n")
    row = collect(out)[0]
    check_true("noinit: status is an error", row.status.startswith("error:"))
    check_na("noinit: init_matches_parsed is NA, not 0", row.data["init_matches_parsed"])


def test_missing_parsed_entry_is_an_error_with_na(tmp: Path) -> None:
    """A check that could not be made is NA, never 0 and never a number."""
    out = tmp / "noparsed"
    write_group(out, parsed_name="somebody_else", fixed={"D": FIXED_POS})
    row = collect(out)[0]
    check_true("noparsed: status is an error", row.status.startswith("error:"))
    for col in (
        "init_matches_parsed",
        "len_matches_parsed",
        "fixed_held",
        "n_fixed",
        "parsed_len",
        "n_chains_designed",
    ):
        check_na(f"noparsed: {col} is NA", row.data[col])


def test_missing_parsed_file_is_an_error_with_na(tmp: Path) -> None:
    out = tmp / "nofile"
    grp = write_group(out, fixed={"D": FIXED_POS})
    (grp / "parsed_pdbs.jsonl").unlink()
    row = collect(out)[0]
    check_true("nofile: status is an error", row.status.startswith("error:"))
    check_na("nofile: n_fixed is NA", row.data["n_fixed"])
    check_na("nofile: len_matches_parsed is NA", row.data["len_matches_parsed"])


def test_wrong_designed_chain_is_an_error(tmp: Path) -> None:
    """A chain letter the structure does not have: the task exits 0 and designs
    nothing, so the collector has to be the one that notices."""
    out = tmp / "wrongchain"
    write_group(out, designed_chains=["Z"], init=NATIVE_D, fixed=None)
    row = collect(out)[0]
    check_true("wrongchain: status is an error", row.status.startswith("error:"))
    check_true("wrongchain: names the chain", "Z" in row.status)
    check_na("wrongchain: n_fixed is NA", row.data["n_fixed"])


def test_fixed_position_out_of_range_is_an_error(tmp: Path) -> None:
    out = tmp / "oob"
    write_group(out, fixed={"D": [3, 99]})
    row = collect(out)[0]
    check_true("oob: status is an error", row.status.startswith("error:"))
    check("oob: n_fixed still counted", row.data["n_fixed"], 2)
    check_na("oob: fixed_held is NA", row.data["fixed_held"])


def test_the_probe_shape(tmp: Path) -> None:
    """The 2-design probe of 2026-10-14, reproduced at its measured dimensions.

    chain D length 144, generated length 144, 17 pinned positions whose residues
    read TSNYWYENTNLIKRSPF in the native chain AND in the generated sequence. This
    is the shape the repaired collector must call OK -- and the shape the old skill
    wrongly described as "the FASTA omits the fixed residues", which is the second
    half of this test: a sequence 17 residues short is an ERROR row, not a binder.
    """
    motif = "TSNYWYENTNLIKRSPF"
    positions = [26, 27, 28, 31, 33, 50, 52, 53, 54, 57, 58, 99, 100, 101, 102, 103, 104]
    native = list("G" * 144)
    designed = list("A" * 144)
    for pos, aa in zip(positions, motif):
        native[pos - 1] = aa
        designed[pos - 1] = aa
    native_s, designed_s = "".join(native), "".join(designed)

    out = tmp / "probe"
    write_group(
        out,
        parsed_chains={"A": NATIVE_A, "D": native_s},
        init=native_s,
        draws=(designed_s,),
        fixed={"D": positions},
    )
    row = collect(out)[0]
    check("probe: status", row.status, "OK")
    check("probe: parsed_len", row.data["parsed_len"], 144)
    check("probe: seq_len", row.data["seq_len"], 144)
    check("probe: n_fixed", row.data["n_fixed"], 17)
    check(
        "probe: fixed_residues spell the motif, chain-prefixed",
        row.data["fixed_residues"],
        "D:" + motif,
    )
    check("probe: init_matches_parsed", row.data["init_matches_parsed"], 1)
    check("probe: len_matches_parsed", row.data["len_matches_parsed"], 1)
    check("probe: fixed_held", row.data["fixed_held"], 1)

    # The old (false) premise, had it been true: 127 residues, motif absent.
    truncated = "".join(c for i, c in enumerate(designed_s) if (i + 1) not in positions)
    out2 = tmp / "probe_truncated"
    write_group(
        out2,
        parsed_chains={"A": NATIVE_A, "D": native_s},
        init=native_s,
        draws=(truncated,),
        fixed={"D": positions},
    )
    row2 = collect(out2)[0]
    check("probe: truncated seq_len", row2.data["seq_len"], 127)
    check_true("probe: truncated is an error", row2.status.startswith("error:"))
    check("probe: truncated len_matches_parsed", row2.data["len_matches_parsed"], 0)

    # fixed_residues follows the LISTED order, which the helper does not sort.
    out3 = tmp / "probe_order"
    write_group(
        out3,
        parsed_chains={"A": NATIVE_A, "D": native_s},
        init=native_s,
        draws=(designed_s,),
        fixed={"D": list(reversed(positions))},
    )
    row3 = collect(out3)[0]
    check(
        "probe: residues follow the listed order",
        row3.data["fixed_residues"],
        "D:" + motif[::-1],
    )
    check("probe: reversed order still held", row3.data["fixed_held"], 1)


def test_input_pdb_records_the_staged_structure(tmp: Path) -> None:
    """F2: which structure was this sequence designed on? No column said so."""
    out = tmp / "inputpdb"
    src = tmp / "inputpdb_src"
    src.mkdir(parents=True, exist_ok=True)
    pdb = src / "backbone.pdb"
    _pdb(pdb)
    write_group(out, fixed={"D": FIXED_POS}, staged_pdb=pdb)
    row = collect(out)[0]
    check("input_pdb: records the symlink target", row.data["input_pdb"], str(pdb))
    # And it survives on an error row, because it does not depend on the jsonls.
    out2 = tmp / "inputpdb_err"
    grp = write_group(out2, fixed={"D": FIXED_POS}, staged_pdb=pdb)
    (grp / "parsed_pdbs.jsonl").unlink()
    row2 = collect(out2)[0]
    check_true("input_pdb: present on an error row", row2.status.startswith("error:"))
    check("input_pdb: still recorded", row2.data["input_pdb"], str(pdb))
    # Absent link -> NA, not "" and not an error.
    out3 = tmp / "inputpdb_none"
    write_group(out3, fixed={"D": FIXED_POS})
    row3 = collect(out3)[0]
    check("input_pdb: row still OK without a link", row3.status, "OK")
    check_na("input_pdb: NA when the link is gone", row3.data["input_pdb"])


def test_fixed_jsonl_gaps_are_errors_not_zero(tmp: Path) -> None:
    """F4: a present fixed_pdbs.jsonl means pins WERE requested."""
    out = tmp / "fixedgap"
    grp = write_group(out, fixed={"D": FIXED_POS})
    (grp / "fixed_pdbs.jsonl").write_text('{"somebody_else": {"D": [3, 11]}}\n')
    row = collect(out)[0]
    check_true("fixedgap: design missing -> error", row.status.startswith("error:"))
    check_na("fixedgap: n_fixed is NA, not 0", row.data["n_fixed"])
    check_true("fixedgap: names the file", "fixed_pdbs.jsonl" in row.status)

    out2 = tmp / "fixedgap_chain"
    grp2 = write_group(out2, fixed={"D": FIXED_POS})
    (grp2 / "fixed_pdbs.jsonl").write_text('{"d1": {"A": []}}\n')
    row2 = collect(out2)[0]
    check_true("fixedgap: chain missing -> error", row2.status.startswith("error:"))
    check_na("fixedgap: n_fixed NA for the chain gap", row2.data["n_fixed"])

    out3 = tmp / "assignedgap"
    grp3 = write_group(out3, fixed={"D": FIXED_POS})
    (grp3 / "assigned_pdbs.jsonl").write_text('{"somebody_else": [["D"], ["A"]]}\n')
    row3 = collect(out3)[0]
    check_true("assignedgap: design missing -> error", row3.status.startswith("error:"))
    check_na("assignedgap: n_chains_designed is NA", row3.data["n_chains_designed"])


def test_init_order_matches_is_informational(tmp: Path) -> None:
    """F3: the element-wise companion settles the multi-chain order question
    without ever failing a row on an unverified assumption."""
    out = tmp / "order"
    both_reversed = NATIVE_D + "/" + NATIVE_A  # designed order is A, D
    write_group(
        out,
        designed_chains=["A", "D"],
        init=both_reversed,
        draws=("C" * len(NATIVE_D) + "/" + "C" * len(NATIVE_A),),
        init_header="_init_sequence, chains=A,D",
    )
    row = collect(out)[0]
    check("order: status is still OK", row.status, "OK")
    check("order: multiset check passes", row.data["init_matches_parsed"], 1)
    check("order: element-wise check reports the difference", row.data["init_order_matches"], 0)
    check("order: length check passes too", row.data["len_matches_parsed"], 1)


def test_no_column_value_is_a_pandas_na_token(tmp: Path) -> None:
    """F8: 'NA' is in read_csv's default na_values, so a bare two-residue pin of
    Asn+Ala would come back from the TSV as a null."""
    out = tmp / "natoken"
    native = list("G" * 20)
    native[2], native[10] = "N", "A"
    designed = list("C" * 20)
    designed[2], designed[10] = "N", "A"
    write_group(
        out,
        parsed_chains={"D": "".join(native)},
        init="".join(native),
        draws=("".join(designed),),
        fixed={"D": [3, 11]},
    )
    row = collect(out)[0]
    check("natoken: status", row.status, "OK")
    check("natoken: prefixed, so never a bare NA", row.data["fixed_residues"], "D:NA")
    na_tokens = {"", "NA", "N/A", "NULL", "NaN", "nan", "None", "null", "<NA>"}
    for col in ("fixed_positions", "fixed_residues", "input_pdb", "status"):
        value = row.data[col] if col != "status" else row.status
        check_true(
            f"natoken: {col} is not a pandas NA token",
            pd.isna(value) or str(value) not in na_tokens,
        )


def test_no_fasta_yields_nothing(tmp: Path) -> None:
    out = tmp / "nofasta"
    write_group(out, design="d1")
    check("nofasta: an unrelated design mints no rows", len(collect(out, "d9")), 0)


# ---------------------------------------------------------------------------
# Manifest-builder tests
# ---------------------------------------------------------------------------


def test_positions_are_bounds_checked(tmp: Path) -> None:
    """Chain D is 17 residues here; position 18 must never reach a container."""
    root = tmp / "bounds"
    df = make_frame(root)
    args = make_args(root, chains_to_design="D", fixed_positions="3,18")
    check_raises(
        "bounds: out-of-range position refused",
        lambda: run_mod.build_atomium_manifest(make_ctx(df, args, root / "out")),
        "out of bounds",
        "17",
    )
    ok = make_args(root, chains_to_design="D", fixed_positions="3,17")
    rows = run_mod.build_atomium_manifest(make_ctx(df, ok, root / "out_ok"))
    check("bounds: in-range accepted", len(rows), 1)
    check(
        "bounds: position list",
        subgroups(root / "out_ok", rows)[0]["fixed"],
        "3 17",
    )


def test_bounds_checked_per_design_from_a_column(tmp: Path) -> None:
    root = tmp / "boundscol"
    df = make_frame(root, col={"d1": "3,17", "d2": "3,40"})
    args = make_args(root, chains_to_design="D", fixed_positions_column="col")
    check_raises(
        "bounds col: the offending design is named",
        lambda: run_mod.build_atomium_manifest(make_ctx(df, args, root / "out")),
        "'d2'",
        "out of bounds",
    )


def test_group_count_must_match_the_designed_chains(tmp: Path) -> None:
    """'3/11' with one designed chain silently pinned only 3 before this."""
    root = tmp / "groups"
    df = make_frame(root)
    args = make_args(root, chains_to_design="D", fixed_positions="3/11")
    check_raises(
        "groups: extra group refused",
        lambda: run_mod.build_atomium_manifest(make_ctx(df, args, root / "out")),
        "2 '/'-separated chain group(s)",
        "1 chain(s) are designed",
    )
    args2 = make_args(root, chains_to_design="A,D", fixed_positions="3")
    check_raises(
        "groups: missing group refused",
        lambda: run_mod.build_atomium_manifest(make_ctx(df, args2, root / "out2")),
        "2 chain(s) are designed",
    )
    both = make_args(root, chains_to_design="A,D", fixed_positions="3/11")
    rows = run_mod.build_atomium_manifest(make_ctx(df, both, root / "out3"))
    check(
        "groups: one group per chain accepted",
        subgroups(root / "out3", rows)[0]["fixed"],
        "3, 11",
    )


def test_mutual_exclusion(tmp: Path) -> None:
    root = tmp / "mutex"
    df = make_frame(root, col={"d1": "3", "d2": "3"}, cc={"d1": "D", "d2": "D"})
    check_raises(
        "mutex: fixed-positions + column",
        lambda: run_mod.build_atomium_manifest(
            make_ctx(
                df,
                make_args(
                    root,
                    chains_to_design="D",
                    fixed_positions="3",
                    fixed_positions_column="col",
                ),
                root / "o1",
            )
        ),
        "mutually exclusive",
        "--fixed-positions-column",
    )
    check_raises(
        "mutex: chains-to-design + column",
        lambda: run_mod.build_atomium_manifest(
            make_ctx(
                df,
                make_args(root, chains_to_design="D", chains_to_design_column="cc"),
                root / "o2",
            )
        ),
        "mutually exclusive",
        "--chains-to-design-column",
    )


def test_column_refusals(tmp: Path) -> None:
    root = tmp / "cols"
    df = make_frame(root, col={"d1": "3", "d2": float("nan")})
    check_raises(
        "cols: a missing value refuses the submit",
        lambda: run_mod.build_atomium_manifest(
            make_ctx(
                df,
                make_args(root, chains_to_design="D", fixed_positions_column="col"),
                root / "o1",
            )
        ),
        "d2",
        "missing or empty",
    )
    check_raises(
        "cols: an absent column names the available ones",
        lambda: run_mod.build_atomium_manifest(
            make_ctx(
                df,
                make_args(root, chains_to_design="D", fixed_positions_column="nope"),
                root / "o2",
            )
        ),
        "not a column of the source table",
        "rfdiffusion3_path",
    )


def test_unknown_chain_letter_refuses(tmp: Path) -> None:
    root = tmp / "chains"
    df = make_frame(root)
    check_raises(
        "chains: a run-wide letter the structure lacks",
        lambda: run_mod.build_atomium_manifest(
            make_ctx(df, make_args(root, chains_to_design="B"), root / "o1")
        ),
        "not polymer chains",
        "present: A, D",
    )
    df2 = make_frame(root, cc={"d1": "D", "d2": "Q"})
    check_raises(
        "chains: a per-design letter the structure lacks",
        lambda: run_mod.build_atomium_manifest(
            make_ctx(
                df2, make_args(root, chains_to_design_column="cc"), root / "o2"
            )
        ),
        "'d2'",
        "Q",
    )


def test_two_submits_are_isolated(tmp: Path) -> None:
    """The shard pattern: submit, collect, submit the next shard into one leaf."""
    root = tmp / "shards"
    out = root / "out"
    df = make_frame(root, names=("d1", "d2", "d3", "d4"))
    args = make_args(root, chains_to_design="D")

    rows1 = run_mod.build_atomium_manifest(
        make_ctx(df.loc[["d1", "d2"]], args, out)
    )
    rows2 = run_mod.build_atomium_manifest(
        make_ctx(df.loc[["d3", "d4"]], args, out)
    )
    sg1, sg2 = subgroups(out, rows1), subgroups(out, rows2)
    check("shards: submit 1 members", sorted(sg1[0]["members"]), ["d1", "d2"])
    check("shards: submit 2 members", sorted(sg2[0]["members"]), ["d3", "d4"])
    check_true("shards: separate grp dirs", sg1[0]["grp"] != sg2[0]["grp"])
    check_true("shards: separate task files", sg1[0]["task"] != sg2[0]["task"])
    check(
        "shards: no counter-named task file survives",
        sorted(p.name for p in (out / "atomium_tasks").glob("task_0.tsv")),
        [],
    )
    check(
        "shards: each grp dir holds only its own inputs",
        [len(list((sg["grp"] / "inputs").glob("*.pdb"))) for sg in (sg1[0], sg2[0])],
        [2, 2],
    )


def test_identical_resubmit_is_idempotent(tmp: Path) -> None:
    root = tmp / "idem"
    out = root / "out"
    df = make_frame(root)
    args = make_args(root, chains_to_design="D")
    rows1 = run_mod.build_atomium_manifest(make_ctx(df, args, out))
    rows2 = run_mod.build_atomium_manifest(make_ctx(df, args, out))
    check("idem: same task file", rows1, rows2)
    check(
        "idem: one grp dir only",
        len([p for p in out.glob("grp_*")]),
        1,
    )
    check(
        "idem: one task file only",
        len(list((out / "atomium_tasks").glob("*.tsv"))),
        1,
    )


def test_changed_signature_is_refused_then_forced(tmp: Path) -> None:
    root = tmp / "resig"
    out = root / "out"
    df = make_frame(root, names=("d1",))
    first = make_args(root, chains_to_design="D", fixed_positions="3")
    run_mod.build_atomium_manifest(make_ctx(df, first, out))
    second = make_args(root, chains_to_design="D", fixed_positions="11")
    check_raises(
        "resig: refused without --force",
        lambda: run_mod.build_atomium_manifest(make_ctx(df, second, out)),
        "DIFFERENT subgroup",
        "--force",
    )
    forced = make_args(root, chains_to_design="D", fixed_positions="11", force=True)
    rows = run_mod.build_atomium_manifest(make_ctx(df, forced, out))
    check("resig: --force leaves one grp dir", len(list(out.glob("grp_*"))), 1)
    check("resig: and it is the new one", subgroups(out, rows)[0]["fixed"], "11")


def test_moved_input_path_is_repointed(tmp: Path) -> None:
    root = tmp / "moved"
    out = root / "out"
    df = make_frame(root, names=("d1",))
    args = make_args(root, chains_to_design="D")
    run_mod.build_atomium_manifest(make_ctx(df, args, out))
    moved = root / "elsewhere"
    moved.mkdir()
    _pdb(moved / "d1.pdb")
    df2 = df.copy()
    df2.at["d1", "rfdiffusion3_path"] = str(moved / "d1.pdb")
    rows = run_mod.build_atomium_manifest(make_ctx(df2, args, out))
    link = sorted((subgroups(out, rows)[0]["grp"] / "inputs").glob("*.pdb"))[0]
    check("moved: symlink repointed", str(link.readlink()), str(moved / "d1.pdb"))


def test_stale_staged_input_is_removed(tmp: Path) -> None:
    """_symlink_member only ever appended; staging now makes the dir MATCH."""
    root = tmp / "stale"
    out = root / "out"
    df = make_frame(root, names=("d1",))
    args = make_args(root, chains_to_design="D")
    rows = run_mod.build_atomium_manifest(make_ctx(df, args, out))
    grp = subgroups(out, rows)[0]["grp"]
    intruder = grp / "inputs" / "intruder.pdb"
    _pdb(intruder)
    rows = run_mod.build_atomium_manifest(make_ctx(df, args, out))
    check(
        "stale: intruder removed",
        sorted(p.name for p in (grp / "inputs").glob("*.pdb")),
        ["d1.pdb"],
    )


def test_chain_spec_expanding_to_nothing_is_refused(tmp: Path) -> None:
    """F1, the blocking one. '[]' is non-empty, expands to no chains, and an empty
    --chain_list means 'design EVERY chain' downstream -- pinning nothing at OK."""
    root = tmp / "emptychains"
    for i, spec in enumerate(("[]", ",", "[ ]")):
        df = make_frame(root, cc={"d1": "D", "d2": spec})
        check_raises(
            f"emptychains: per-design {spec!r} refused",
            lambda d=df, n=i: run_mod.build_atomium_manifest(
                make_ctx(d, make_args(root, chains_to_design_column="cc"), root / f"o{n}")
            ),
            "'d2'",
            "expands to NO chains",
        )
    check_raises(
        "emptychains: the run-wide flag is refused once, before the loop",
        lambda: run_mod.build_atomium_manifest(
            make_ctx(
                make_frame(root), make_args(root, chains_to_design="[]"), root / "flag"
            )
        ),
        "--chains-to-design",
        "expands to NO chains",
    )
    # With positions in play the refusal must still come BEFORE any staging.
    out = root / "withpos"
    df3 = make_frame(root, cc={"d1": "D", "d2": "[]"})
    check_raises(
        "emptychains: refused with a fixed-positions spec too",
        lambda: run_mod.build_atomium_manifest(
            make_ctx(
                df3,
                make_args(root, chains_to_design_column="cc", fixed_positions="3"),
                out,
            )
        ),
        "expands to NO chains",
    )
    check("emptychains: nothing staged", sorted(p.name for p in out.glob("grp_*")), [])


def test_every_bad_design_is_named_in_one_raise(tmp: Path) -> None:
    """F6: one raise naming every offender, matching _column_specs' shape."""
    root = tmp / "batch"
    df = make_frame(root, names=("d1", "d2", "d3"), cc={"d1": "D", "d2": "Q", "d3": "Z"})
    check_raises(
        "batch: both offenders named once",
        lambda: run_mod.build_atomium_manifest(
            make_ctx(df, make_args(root, chains_to_design_column="cc"), root / "o1")
        ),
        "2 of 3 ready design(s) failed validation",
        "'d2'",
        "'d3'",
        "-f/--filter",
    )
    df2 = make_frame(root, names=("d1", "d2"), col={"d1": "99", "d2": "40"})
    check_raises(
        "batch: position failures are collected too",
        lambda: run_mod.build_atomium_manifest(
            make_ctx(
                df2,
                make_args(root, chains_to_design="D", fixed_positions_column="col"),
                root / "o2",
            )
        ),
        "2 of 2 ready design(s) failed validation",
        "'d1'",
        "'d2'",
    )


def test_pandas_na_cell_is_reported_as_missing(tmp: Path) -> None:
    """F7: pd.NA stringified to '<NA>' and died later with a misleading message."""
    root = tmp / "nacell"
    df = make_frame(root, col={"d1": "3", "d2": pd.NA})
    check_raises(
        "nacell: pd.NA is missing, not a spec",
        lambda: run_mod.build_atomium_manifest(
            make_ctx(
                df,
                make_args(root, chains_to_design="D", fixed_positions_column="col"),
                root / "o1",
            )
        ),
        "missing or empty",
        "d2",
    )


def test_bias_jsonl_is_content_addressed(tmp: Path) -> None:
    """F5: a fixed name is overwritten by the next arm into the same leaf."""
    root = tmp / "bias"
    out = root / "out"
    df = make_frame(root, names=("d1",))
    rows1 = run_mod.build_atomium_manifest(
        make_ctx(df, make_args(root, chains_to_design="D", bias_aa="H:1.39"), out)
    )
    rows2 = run_mod.build_atomium_manifest(
        make_ctx(df, make_args(root, chains_to_design="D", bias_aa="C:-1.0"), out)
    )
    files = sorted(p.name for p in out.glob("bias_AA*.jsonl"))
    check("bias: two arms, two files", len(files), 2)
    extras = [Path(r[0]).read_text().split("\t")[-1].strip() for r in (rows1[0:1] + rows2[0:1])]
    check_true("bias: the two argv tails name different files", extras[0] != extras[1])
    for extra in extras:
        path = Path(extra.split("--bias_AA_jsonl ")[1].split()[0])
        check_true(f"bias: {path.name} still exists", path.is_file())
    # Same bias twice -> same file, no accumulation.
    run_mod.build_atomium_manifest(
        make_ctx(df, make_args(root, chains_to_design="D", bias_aa="H:1.39"), out)
    )
    check("bias: idempotent", len(list(out.glob("bias_AA*.jsonl"))), 2)


def test_cpu_flags_and_spec() -> None:
    check("spec: name", TOOL.name, "atomium")
    check("spec: action", TOOL.action, "create")
    check("spec: default input column", TOOL.default_input_column, "rfdiffusion3_path")
    parser = ArgumentParser()
    TOOL.add_run_args_fn(parser)
    defaults = {a.dest: a.default for a in parser._actions}
    check("flags: fixed_positions_column default", defaults["fixed_positions_column"], "")
    check("flags: chains_to_design_column default", defaults["chains_to_design_column"], "")
    check("flags: model_noise default", defaults["model_noise"], "n05")
    check("flags: sampling_temp default", defaults["sampling_temp"], "0.1")
    check("flags: num_seq_per_target default", defaults["num_seq_per_target"], 2)


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        test_init_sequence_is_skipped_with_the_real_header(tmp)
        test_bare_init_header_still_skipped(tmp)
        test_row_names_are_a1_a2(tmp)
        test_trust_columns_on_a_good_row(tmp)
        test_no_fixed_positions(tmp)
        test_all_chains_designed_without_assigned_file(tmp)
        test_fixed_position_not_held_is_an_error(tmp)
        test_truncated_sequence_is_an_error(tmp)
        test_init_mismatch_is_an_error(tmp)
        test_missing_init_entry_is_an_error(tmp)
        test_missing_parsed_entry_is_an_error_with_na(tmp)
        test_missing_parsed_file_is_an_error_with_na(tmp)
        test_wrong_designed_chain_is_an_error(tmp)
        test_fixed_position_out_of_range_is_an_error(tmp)
        test_the_probe_shape(tmp)
        test_input_pdb_records_the_staged_structure(tmp)
        test_fixed_jsonl_gaps_are_errors_not_zero(tmp)
        test_init_order_matches_is_informational(tmp)
        test_no_column_value_is_a_pandas_na_token(tmp)
        test_no_fasta_yields_nothing(tmp)
        test_positions_are_bounds_checked(tmp)
        test_bounds_checked_per_design_from_a_column(tmp)
        test_group_count_must_match_the_designed_chains(tmp)
        test_mutual_exclusion(tmp)
        test_column_refusals(tmp)
        test_unknown_chain_letter_refuses(tmp)
        test_two_submits_are_isolated(tmp)
        test_identical_resubmit_is_idempotent(tmp)
        test_changed_signature_is_refused_then_forced(tmp)
        test_moved_input_path_is_repointed(tmp)
        test_stale_staged_input_is_removed(tmp)
        test_chain_spec_expanding_to_nothing_is_refused(tmp)
        test_every_bad_design_is_named_in_one_raise(tmp)
        test_pandas_na_cell_is_reported_as_missing(tmp)
        test_bias_jsonl_is_content_addressed(tmp)
        test_cpu_flags_and_spec()

    if FAILURES:
        print(f"FAILED {len(FAILURES)} of {CHECKS} checks:")
        for f in FAILURES:
            print(f"  - {f}")
        return 1
    print(f"OK: {CHECKS} checks passed.")
    print(
        "NOTE: this proves the refusal and trust-column contracts against frames "
        "this file wrote. It cannot prove the multi-chain FASTA layout, the real "
        "helper jsonl shapes, or that gemmi's chain lengths equal "
        "parse_multiple_chains.py's -- see the module docstring."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
