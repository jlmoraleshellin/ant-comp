#!/usr/bin/env python3
"""Self-test for the proteinmpnn fork: per-row column flags, and the staging fix.

WHAT THIS PROVES, AND WHAT IT CANNOT.

It proves the *wiring and the refusal contracts*: that with neither new flag the
fork's manifest, sub-manifests and staged layout mean exactly what the bundled
builder's do (same designs, same six fields, same batching -- only the generated
``grp_``/``task_`` NAMES differ, deliberately); that each design gets ITS OWN
spec rather than the first or last row's broadcast to all; that the colon/comma
expansion is upstream's ``_parse_positions`` object and not a reimplementation;
and that every refusal path raises before anything is queued, naming the design
and the column.

It also proves the **staging-collision fix** at the filesystem level:
``test_upstream_collides_across_submits`` reproduces the bug against UPSTREAM's
builder (two submits to one out_dir both write ``grp_0``, which ends up holding
both designs under the second submit's spec), and the fork is then shown to keep
them apart in BOTH shapes -- the observed crash case (different chain lengths)
and, separately, the SAME-LENGTH case, which under the old numbering would have
produced a wrong answer with no crash and no error status. Plus: ten shard
submits to one leaf, a ``--force`` re-run of the identical design (no duplicate
subgroup), a changed spec for an already-staged design (refused; replaced under
``--force``), a partial-overlap ``--force`` (still refused -- never deletes a
result the caller is not redoing), and a moved input path (symlink repointed in
place rather than orphaning a stale FASTA).

What those filesystem tests still CANNOT establish: that ProteinMPNN itself,
in the real container, then behaves as expected on the isolated directories
(only a real run shows that), and that no OTHER shared-per-leaf artifact
collides across submits. Two that demonstrably still do, both in prosapia and
not fixable from here, are the top-level manifest
``.manifests/<table>_<leaf>_manifest.txt`` and the task logs
``<script>_<taskid>.out/.err/.exit`` + ``<script>_modal.json`` -- task ids
restart at 1 every submit. They are safe only under the existing rule: submit,
wait, collect, then submit the next shard.

It CANNOT prove the premise or the mapping. Every frame here is written by this
file, so it agrees with this file by construction. These remain UNKNOWN until
the tool is run on real data:

  * whether the values in a real ``motifrmsd_design_fixed_positions`` column are
    in ProteinMPNN's numbering (1..L per chain, file order) for the designs in
    this campaign's table -- motifrmsd claims to guarantee it, this file does
    not check it, and nothing at run time can;
  * whether ``motifrmsd_fixed_positions_chain`` really names the chain
    ProteinMPNN's ``--chain_list`` indexes (an author chain letter that
    ``parse_multiple_chains.py`` does not see would fix nothing, loudly or
    otherwise);
  * whether ``make_fixed_positions_dict.py`` accepts a single-group
    ``--position_list`` against a single-chain ``--chain_list`` in the real
    container -- the group-count contract is ProteinMPNN's, exercised nowhere
    here;
  * what 338 single-design subgroups cost in wall time (one model load each).

Run:  uv run python tools/proteinmpnn/test_proteinmpnn.py
"""

import sys
import tempfile
from argparse import ArgumentParser, Namespace
from pathlib import Path

import pandas as pd

from prosapia.core import ManifestCtx
from prosapia.core.tool_registry import _load as load_tool_spec, get_builtin
from prosapia.tools.proteinmpnn import run_proteinmpnn as builtin_run
from prosapia.utils import parse_positions
from prosapia.utils.positions import parse_chain_positions

HERE = Path(__file__).resolve().parent

# Load the tool exactly as `sapia` discovers it: through the registry's own
# spec loader, which builds the synthetic package that makes spec.py's relative
# import of run_proteinmpnn resolve. The module under test is therefore the one
# a real run would use, not a second copy loaded by hand -- which is what makes
# the identity checks below meaningful.
TOOL = load_tool_spec(HERE / "spec.py")
fork = sys.modules[TOOL.build_manifest_fn.__module__]


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
                FAILURES.append(
                    f"{label}: message {message!r} lacks {needle!r}"
                )
        return
    FAILURES.append(f"{label}: expected a raise, got none")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Three designs, each with a DISTINCT three-fragment spec, modelled on the real
# campaign's motifrmsd_design_fixed_positions values.
SPECS = {
    "d1": "11:18,45:54,87:92",
    "d2": "22:29,68:77,119:124",
    "d3": "25:32,58:67,102:107",
}


def _pdb(path: Path) -> None:
    """A minimal two-chain PDB. Only --symmetry reads a structure, so the content
    matters only in that the file must exist to be symlinked."""
    lines = []
    serial = 1
    for chain in ("A", "D"):
        for resnum in (1, 2, 3):
            lines.append(
                f"ATOM  {serial:5d}  CA  GLY {chain}{resnum:4d}    "
                f"{resnum * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C"
            )
            serial += 1
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def make_frame(tmp: Path, fixed_col: dict[str, object] | None = None,
               chain_col: dict[str, object] | None = None) -> pd.DataFrame:
    inputs = tmp / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    rows = []
    for name in SPECS:
        pdb = inputs / f"{name}.pdb"
        _pdb(pdb)
        row: dict[str, object] = {"name": name, "rfdiffusion3_path": str(pdb)}
        if fixed_col is not None:
            row["motifrmsd_design_fixed_positions"] = fixed_col[name]
        if chain_col is not None:
            row["motifrmsd_fixed_positions_chain"] = chain_col[name]
        rows.append(row)
    return pd.DataFrame(rows).set_index("name")


def make_args(run_dir: Path, **over) -> Namespace:
    base = dict(
        run_dir=run_dir,
        input_column="rfdiffusion3_path",
        force=False,
        designs_per_task=10,
        chains_to_design="",
        fixed_positions="",
        tied_positions="",
        symmetry=None,
        bias_aa="",
        set=[],
        num_seq_per_target=2,
        sampling_temp="0.2",
        seed=37,
        batch_size=1,
        fixed_positions_column="",
        chains_to_design_column="",
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


def _redact(obj, out_dir: Path):
    """Replace out_dir (and anything under it) in a nested structure with <OUT>."""
    token = str(out_dir)
    if isinstance(obj, str):
        return obj.replace(token, "<OUT>")
    if isinstance(obj, Path):
        return str(obj).replace(token, "<OUT>")
    if isinstance(obj, (list, tuple)):
        return type(obj)(_redact(x, out_dir) for x in obj)
    if isinstance(obj, dict):
        return {_redact(k, out_dir): _redact(v, out_dir) for k, v in obj.items()}
    return obj


def snapshot(out_dir: Path) -> list[str]:
    """Everything the builder wrote under out_dir, with out_dir itself redacted.

    Sub-manifest text, grp_ dir names, the symlinks staged into each grp_/inputs
    and where they point. Used for IDEMPOTENCE checks (the same builder, the
    same out_dir, twice), where the names must be literally identical.
    """
    token = str(out_dir)
    out: list[str] = []
    for path in sorted(out_dir.rglob("*")):
        rel = str(path.relative_to(out_dir))
        if path.is_symlink():
            out.append(f"link {rel} -> {path.readlink()}")
        elif path.is_dir():
            out.append(f"dir  {rel}")
        else:
            out.append(f"file {rel}\n{path.read_text().replace(token, '<OUT>')}")
    return out


def subgroups(out_dir: Path, rows: list[tuple[str, ...]] | None = None) -> list[dict]:
    """Every subgroup this submit staged, read the way proteinmpnn.sh reads it.

    Follows the manifest rows (field 1 = a task sub-manifest) to each
    sub-manifest line (grp_dir, chains, fixed, tie_mode, tied, extra) and then
    to the symlinks actually staged in that grp dir -- i.e. the association the
    task script will act on. ``rows=None`` falls back to globbing every task
    file in the dir, which is only right when the dir holds ONE submit.
    """
    if rows is None:
        task_files = sorted((out_dir / "proteinmpnn_tasks").glob("task_*.tsv"))
    else:
        task_files = [Path(r[0]) for r in rows]
    out: list[dict] = []
    for task_file in task_files:
        for line in task_file.read_text().splitlines():
            if not line:
                continue
            grp, chains, fixed, tie, tied, extra = line.split("\t")
            grp_dir = Path(grp)
            out.append(
                {
                    "task": task_file.name,
                    "grp": grp_dir,
                    "chains": chains,
                    "fixed": fixed,
                    "tie": tie,
                    "tied": tied,
                    "extra": extra,
                    "members": {
                        link.stem: str(link.readlink())
                        for link in sorted((grp_dir / "inputs").glob("*.pdb"))
                    },
                }
            )
    return out


def semantics(out_dir: Path, rows: list[tuple[str, ...]]) -> list:
    """The naming-independent content of a submit: who got which params, batched how.

    Everything a subgroup means, with the generated grp_/task_ NAMES removed --
    the names are an implementation detail (upstream counts, the fork digests),
    and nothing downstream reads them. Tasks are compared as a set of subgroup
    sets so bin-packing is still checked.
    """
    by_task: dict[str, list] = {}
    for sg in subgroups(out_dir, rows):
        by_task.setdefault(sg["task"], []).append(
            (
                tuple(sorted(sg["members"].items())),
                sg["chains"],
                sg["fixed"],
                sg["tie"],
                sg["tied"],
                sg["extra"],
            )
        )
    return sorted(tuple(sorted(v)) for v in by_task.values())


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_reuse_not_reimplementation() -> None:
    """The fork calls upstream's parsing objects, it does not copy them."""
    check_true(
        "_parse_positions is upstream's object",
        fork._parse_positions is builtin_run._parse_positions,
    )
    check_true(
        "_parse_chains is upstream's object",
        fork._parse_chains is builtin_run._parse_chains,
    )
    check_true(
        "_resolve_chains is upstream's object",
        fork._resolve_chains is builtin_run._resolve_chains,
    )
    check_true(
        "_write_bias_jsonl is upstream's object",
        fork._write_bias_jsonl is builtin_run._write_bias_jsonl,
    )
    check_true(
        "_pack_subgroups is upstream's object",
        fork._pack_subgroups is builtin_run._pack_subgroups,
    )
    # and upstream's _parse_positions is itself the thin formatter over the
    # shared parse_positions -> parse_chain_positions, which the fork therefore
    # inherits rather than re-deriving.
    check_true(
        "upstream _parse_positions delegates to parse_positions",
        builtin_run.parse_positions is parse_positions,
    )
    check(
        "shared parser expands a 3-fragment spec",
        parse_chain_positions("11:18,45:54,87:92"),
        list(range(11, 19)) + list(range(45, 55)) + list(range(87, 93)),
    )


def test_no_flag_regression(tmp: Path) -> None:
    """With neither new flag, the fork's output is the built-in's, byte for byte."""
    root = tmp / "noflag"
    df = make_frame(root)

    args_a = make_args(root)
    out_a = root / "builtin_out"
    rows_a = builtin_run.build_proteinmpnn_manifest(make_ctx(df, args_a, out_a))

    args_b = make_args(root)
    out_b = root / "fork_out"
    rows_b = fork.build_proteinmpnn_manifest(make_ctx(df, args_b, out_b))

    # The generated grp_/task_ NAMES differ by design (upstream counts from 0,
    # the fork digests content -- see the staging-collision fix). Everything the
    # task script and the collector actually read is compared instead: which
    # designs are staged where, pointing at what, with which six fields, batched
    # into which tasks.
    check("no-flag: one task either way", len(rows_b), len(rows_a))
    check(
        "no-flag: semantics identical",
        _redact(semantics(out_b, rows_b), out_b),
        _redact(semantics(out_a, rows_a), out_a),
    )

    # Same again with the flag that already existed, exercised the old way.
    root2 = tmp / "noflag_fixed"
    df2 = make_frame(root2)
    out_c = root2 / "builtin_out"
    out_d = root2 / "fork_out"
    kw = dict(fixed_positions="5:12,30:37", chains_to_design="D")
    rows_c = builtin_run.build_proteinmpnn_manifest(
        make_ctx(df2, make_args(root2, **kw), out_c)
    )
    rows_d = fork.build_proteinmpnn_manifest(
        make_ctx(df2, make_args(root2, **kw), out_d)
    )
    check(
        "no-flag: --fixed-positions semantics identical",
        _redact(semantics(out_d, rows_d), out_d),
        _redact(semantics(out_c, rows_c), out_c),
    )
    # ...and it really did fix something (guards against both sides being empty).
    sg = subgroups(out_d, rows_d)
    check_true(
        "no-flag: --fixed-positions expanded",
        sg[0]["fixed"] == "5 6 7 8 9 10 11 12 30 31 32 33 34 35 36 37",
    )
    check("no-flag: one shared subgroup", len(sg), 1)


def _spec_by_design(out_dir: Path, rows=None) -> dict[str, str]:
    """Map each design name to the fixed-position list its subgroup carries."""
    return {
        name: sg["fixed"] for sg in subgroups(out_dir, rows) for name in sg["members"]
    }


def _chains_by_design(out_dir: Path, rows=None) -> dict[str, str]:
    return {
        name: sg["chains"] for sg in subgroups(out_dir, rows) for name in sg["members"]
    }


def test_per_design_specs(tmp: Path) -> None:
    """Each design gets ITS OWN spec -- not the first or last row broadcast."""
    root = tmp / "percol"
    df = make_frame(root, fixed_col=dict(SPECS))
    out = root / "out"
    rows = fork.build_proteinmpnn_manifest(
        make_ctx(
            df,
            make_args(
                root,
                chains_to_design="D",
                fixed_positions_column="motifrmsd_design_fixed_positions",
            ),
            out,
        )
    )
    got = _spec_by_design(out)

    def expand(spec: str) -> str:
        return " ".join(str(p) for p in parse_chain_positions(spec))

    for name, spec in SPECS.items():
        check(f"per-design spec for {name}", got.get(name), expand(spec))

    # The three specs must be three DIFFERENT strings -- a broadcast bug would
    # make them equal and still pass a per-design equality check against a
    # constant fixture, so assert distinctness explicitly.
    check("per-design: three distinct specs", len(set(got.values())), 3)
    # Distinct signatures => one subgroup each, bin-packed into one task.
    check("per-design: one task", len(rows), 1)
    check("per-design: three subgroups", len(subgroups(out, rows)), 3)
    check("per-design: one grp dir each", len(sorted(out.glob("grp_*"))), 3)
    # The chain list is still the run-wide flag.
    check("per-design: chains from the flag", set(_chains_by_design(out).values()), {"D"})


def test_chains_to_design_column(tmp: Path) -> None:
    """The symmetric per-row chain flag, and the heterogeneous case it exists for."""
    root = tmp / "chaincol"
    df = make_frame(
        root,
        fixed_col=dict(SPECS),
        chain_col={"d1": "D", "d2": "D", "d3": "A"},
    )
    out = root / "out"
    fork.build_proteinmpnn_manifest(
        make_ctx(
            df,
            make_args(
                root,
                fixed_positions_column="motifrmsd_design_fixed_positions",
                chains_to_design_column="motifrmsd_fixed_positions_chain",
            ),
            out,
        )
    )
    check("chain column: per-design chains", _chains_by_design(out),
          {"d1": "D", "d2": "D", "d3": "A"})
    # A per-row chain list satisfies the "positions need a chain list" guard, so
    # --chains-to-design was not required above; that it got here proves it.
    check("chain column: specs still per-design",
          len(set(_spec_by_design(out).values())), 3)

    # The chain mini-language is upstream's, applied per row.
    root2 = tmp / "chaincol_range"
    df2 = make_frame(root2, chain_col={"d1": "A:C", "d2": "A:C", "d3": "A:C"})
    out2 = root2 / "out"
    fork.build_proteinmpnn_manifest(
        make_ctx(
            df2,
            make_args(root2, chains_to_design_column="motifrmsd_fixed_positions_chain"),
            out2,
        )
    )
    check("chain column: range expands", set(_chains_by_design(out2).values()), {"A B C"})


def test_both_flags_refused(tmp: Path) -> None:
    """Both sources for one concept => a clear error, never silent precedence."""
    root = tmp / "both"
    df = make_frame(root, fixed_col=dict(SPECS))
    check_raises(
        "both fixed flags refused",
        lambda: fork.build_proteinmpnn_manifest(
            make_ctx(
                df,
                make_args(
                    root,
                    chains_to_design="D",
                    fixed_positions="5:12",
                    fixed_positions_column="motifrmsd_design_fixed_positions",
                ),
                root / "out1",
            )
        ),
        "--fixed-positions",
        "--fixed-positions-column",
        "mutually exclusive",
    )
    df2 = make_frame(root, chain_col={k: "D" for k in SPECS})
    check_raises(
        "both chain flags refused",
        lambda: fork.build_proteinmpnn_manifest(
            make_ctx(
                df2,
                make_args(
                    root,
                    chains_to_design="D",
                    chains_to_design_column="motifrmsd_fixed_positions_chain",
                ),
                root / "out2",
            )
        ),
        "--chains-to-design-column",
        "mutually exclusive",
    )
    # Refused BEFORE anything was staged: no grp dirs, no sub-manifests.
    check("both flags: nothing staged", sorted((root / "out1").glob("grp_*")), [])
    check(
        "both flags: no sub-manifests",
        (root / "out1" / "proteinmpnn_tasks").exists(),
        False,
    )


def test_neither_flag_is_full_redesign(tmp: Path) -> None:
    """Documented resolution of the ambiguity: absent means full redesign, as today.

    There is no way to detect "the caller meant to fix something and forgot", so
    the fork does not guess. Neither flag => an empty fixed-position list, no
    --fixed_positions_jsonl, every position redesigned -- the built-in's
    behaviour, unchanged.
    """
    root = tmp / "neither"
    df = make_frame(root, fixed_col=dict(SPECS))
    out = root / "out"
    fork.build_proteinmpnn_manifest(make_ctx(df, make_args(root), out))
    check("neither flag: no fixed positions", set(_spec_by_design(out).values()), {""})
    check("neither flag: no chain list", set(_chains_by_design(out).values()), {""})
    # All three share one signature, so they batch into a single subgroup.
    check("neither flag: one subgroup", len(sorted(out.glob("grp_*"))), 1)


def test_missing_value_refuses(tmp: Path) -> None:
    """A design with no value is refused by NAME -- never silently full-redesigned."""
    root = tmp / "missing"
    for label, bad in (
        ("NaN", float("nan")),
        ("None", None),
        ("empty string", ""),
        ("whitespace", "   "),
    ):
        specs = dict(SPECS)
        specs["d2"] = bad
        df = make_frame(root, fixed_col=specs)
        check_raises(
            f"missing value ({label}) refuses the submit",
            lambda df=df: fork.build_proteinmpnn_manifest(
                make_ctx(
                    df,
                    make_args(
                        root,
                        chains_to_design="D",
                        fixed_positions_column="motifrmsd_design_fixed_positions",
                    ),
                    root / f"out_{label.replace(' ', '_')}",
                )
            ),
            "d2",                                      # names the design
            "motifrmsd_design_fixed_positions",        # names the column
            "1 of 3",                                  # names the count
            "-f/--filter",                             # names the remedy
        )

    # Every offender is named in ONE raise, so the table is fixed in one pass.
    specs = dict(SPECS)
    specs["d1"] = float("nan")
    specs["d3"] = ""
    df = make_frame(root, fixed_col=specs)
    check_raises(
        "missing values: all offenders named at once",
        lambda: fork.build_proteinmpnn_manifest(
            make_ctx(
                df,
                make_args(
                    root,
                    chains_to_design="D",
                    fixed_positions_column="motifrmsd_design_fixed_positions",
                ),
                root / "out_multi",
            )
        ),
        "2 of 3",
        "d1, d3",
    )
    # Nothing was staged for the good rows either: the submit is all-or-nothing.
    check(
        "missing values: no partial staging",
        sorted((root / "out_multi").glob("grp_*")),
        [],
    )

    # A chain column is validated the same way.
    df2 = make_frame(root, chain_col={"d1": "D", "d2": float("nan"), "d3": "D"})
    check_raises(
        "missing chain value refuses the submit",
        lambda: fork.build_proteinmpnn_manifest(
            make_ctx(
                df2,
                make_args(
                    root,
                    chains_to_design_column="motifrmsd_fixed_positions_chain",
                ),
                root / "out_chain",
            )
        ),
        "d2",
        "motifrmsd_fixed_positions_chain",
    )


def test_missing_column_refuses(tmp: Path) -> None:
    """A typo'd column name is an error naming the available columns, not a no-op."""
    root = tmp / "nocol"
    df = make_frame(root, fixed_col=dict(SPECS))
    check_raises(
        "unknown column refused",
        lambda: fork.build_proteinmpnn_manifest(
            make_ctx(
                df,
                make_args(
                    root,
                    chains_to_design="D",
                    fixed_positions_column="motifrmsd_fixed_postions",  # typo
                ),
                root / "out",
            )
        ),
        "motifrmsd_fixed_postions",
        "not a column",
        "motifrmsd_design_fixed_positions",  # lists what IS available
    )


def test_requires_a_chain_list(tmp: Path) -> None:
    """The new flag inherits the guardrail: positions need chains to map onto.

    Without it, ProteinMPNN's make_fixed_positions_dict matches no chain, the
    fixed dict comes back empty for every chain, and the whole structure is
    redesigned at exit 0.
    """
    root = tmp / "nochains"
    df = make_frame(root, fixed_col=dict(SPECS))
    check_raises(
        "fixed-positions-column without chains refused",
        lambda: fork.build_proteinmpnn_manifest(
            make_ctx(
                df,
                make_args(
                    root,
                    fixed_positions_column="motifrmsd_design_fixed_positions",
                ),
                root / "out",
            )
        ),
        "--fixed-positions-column",
        "--chains-to-design",
    )


def test_tied_positions_still_guarded(tmp: Path) -> None:
    """The fork did not widen the existing guards while extending them."""
    root = tmp / "tied"
    df = make_frame(root)
    check_raises(
        "tied-positions without chains still refused",
        lambda: fork.build_proteinmpnn_manifest(
            make_ctx(df, make_args(root, tied_positions="1:8/1:8"), root / "out1")
        ),
        "--tied-positions",
    )
    check_raises(
        "symmetry + tied-positions still refused",
        lambda: fork.build_proteinmpnn_manifest(
            make_ctx(
                df,
                make_args(root, symmetry="2", tied_positions="1:8/1:8"),
                root / "out2",
            )
        ),
        "mutually exclusive",
    )
    check_raises(
        "bad symmetry still refused",
        lambda: fork.build_proteinmpnn_manifest(
            make_ctx(df, make_args(root, symmetry="C12"), root / "out3")
        ),
        "--symmetry",
    )


def test_symmetry_broadcast_applies_to_the_column(tmp: Path) -> None:
    """A per-row single-unit spec is broadcast across tied chains, as the flag is."""
    root = tmp / "symm"
    one = {k: "1:3" for k in SPECS}
    df = make_frame(root, fixed_col=one)
    out = root / "out"
    fork.build_proteinmpnn_manifest(
        make_ctx(
            df,
            make_args(
                root,
                chains_to_design="A,D",
                symmetry="2",
                fixed_positions_column="motifrmsd_design_fixed_positions",
            ),
            out,
        )
    )
    check(
        "symmetry: lone group broadcast to 2 chains",
        set(_spec_by_design(out).values()),
        {"1 2 3, 1 2 3"},
    )


# ---------------------------------------------------------------------------
# The staging-collision fix (separate submits to the same table + same leaf)
# ---------------------------------------------------------------------------


def _pdb_n(path: Path, nres: int, chain: str = "D") -> None:
    """A one-chain PDB of a chosen length, so "same length" can be made explicit."""
    lines = []
    for i in range(1, nres + 1):
        lines.append(
            f"ATOM  {i:5d}  CA  GLY {chain}{i:4d}    "
            f"{i * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def one_design(tmp: Path, name: str, spec: str, nres: int) -> pd.DataFrame:
    """A one-row source frame -- what ONE shard's submit sees as ctx.ready."""
    inputs = tmp / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    pdb = inputs / f"{name}.pdb"
    _pdb_n(pdb, nres)
    return pd.DataFrame(
        [
            {
                "name": name,
                "rfdiffusion3_path": str(pdb),
                "motifrmsd_design_fixed_positions": spec,
            }
        ]
    ).set_index("name")


def _submit(builder, df, root: Path, out: Path, **over):
    """One `sapia run` against a given out_dir (= one table + one -l label)."""
    args = make_args(
        root,
        chains_to_design="D",
        fixed_positions_column="motifrmsd_design_fixed_positions",
        **over,
    )
    return builder(make_ctx(df, args, out))


def test_upstream_collides_across_submits(tmp: Path) -> None:
    """The bug itself, reproduced against UPSTREAM's builder, so the fix has a baseline.

    Upstream numbers subgroups with a counter that restarts at 0 every submit,
    and _symlink_member appends without clearing. Two submits to one out_dir
    therefore share grp_0, and the second submit's single protein_mpnn_run call
    would process BOTH PDBs under the second submit's fixed-positions spec.
    """
    root = tmp / "upstream_collision"
    out = root / "out"
    a = one_design(root, "dA", "1:5", nres=105)
    b = one_design(root, "dB", "100:105", nres=200)

    # Upstream has no column flags: give each submit its own run-wide spec, which
    # is the same thing from the staging code's point of view.
    rows1 = builtin_run.build_proteinmpnn_manifest(
        make_ctx(a, make_args(root, chains_to_design="D", fixed_positions="1:5"), out)
    )
    rows2 = builtin_run.build_proteinmpnn_manifest(
        make_ctx(
            b, make_args(root, chains_to_design="D", fixed_positions="100:105"), out
        )
    )

    check("upstream: both submits used grp_0", sorted(q.name for q in out.glob("grp_*")),
          ["grp_0"])
    check(
        "upstream: grp_0 ends up holding BOTH designs",
        sorted(q.stem for q in (out / "grp_0" / "inputs").glob("*.pdb")),
        ["dA", "dB"],
    )
    # ...and submit 2's spec is the one that subgroup carries, so dA would be
    # redesigned under dB's spec. With these lengths that is the measured
    # IndexError; at equal lengths it is a silent wrong answer.
    check(
        "upstream: the shared subgroup carries submit 2's spec",
        subgroups(out, rows2)[0]["fixed"],
        "100 101 102 103 104 105",
    )
    check("upstream: submit 2 overwrote submit 1's task file", rows2, rows1)


def _isolation_checks(label: str, out: Path, rows1, rows2, want: dict[str, str]) -> None:
    """Both submits' subgroups are physically separate and carry their own spec."""
    sg1, sg2 = subgroups(out, rows1), subgroups(out, rows2)
    check(f"{label}: submit 1 staged one subgroup", len(sg1), 1)
    check(f"{label}: submit 2 staged one subgroup", len(sg2), 1)
    check_true(f"{label}: the two subgroups are different dirs", sg1[0]["grp"] != sg2[0]["grp"])
    check(f"{label}: exactly two subgroup dirs exist", len(sorted(out.glob("grp_*"))), 2)
    # No directory holds both designs -- the whole point.
    for sg in sg1 + sg2:
        check(f"{label}: {sg['grp'].name} holds exactly one design", len(sg["members"]), 1)
    # Each design's subgroup carries that design's OWN spec.
    got = {}
    for sg in sg1 + sg2:
        for name in sg["members"]:
            got[name] = sg["fixed"]
    check(f"{label}: each design keeps its own spec", got, want)
    # The task sub-manifests did not overwrite each other either.
    check(f"{label}: two distinct task files", len({rows1[0][0], rows2[0][0]}), 2)
    check_true(f"{label}: submit 1's task file still exists", Path(rows1[0][0]).is_file())


def test_second_submit_isolated_crash_case(tmp: Path) -> None:
    """The exact observed failure: two submits, different chain lengths."""
    root = tmp / "fix_crash"
    out = root / "out"
    a = one_design(root, "dA", "1:5", nres=105)
    b = one_design(root, "dB", "100:105", nres=200)
    rows1 = _submit(fork.build_proteinmpnn_manifest, a, root, out)
    rows2 = _submit(fork.build_proteinmpnn_manifest, b, root, out)
    _isolation_checks(
        "crash case",
        out,
        rows1,
        rows2,
        {"dA": "1 2 3 4 5", "dB": "100 101 102 103 104 105"},
    )


def test_second_submit_isolated_same_length(tmp: Path) -> None:
    """The case that would NOT have crashed -- the silent one.

    Two designs of the SAME chain length, each with an in-range spec. Under the
    old grp_0 collision ProteinMPNN would have run both PDBs in one call under
    submit 2's spec, exited 0 and written a status OK row: design dA redesigned
    with dB's motif held fixed, with nothing anywhere to show for it. The crash
    test alone cannot prove this case, which is why it is tested separately.
    """
    root = tmp / "fix_silent"
    out = root / "out"
    a = one_design(root, "dA", "1:5", nres=105)
    b = one_design(root, "dB", "6:10", nres=105)
    rows1 = _submit(fork.build_proteinmpnn_manifest, a, root, out)
    rows2 = _submit(fork.build_proteinmpnn_manifest, b, root, out)
    _isolation_checks(
        "same length",
        out,
        rows1,
        rows2,
        {"dA": "1 2 3 4 5", "dB": "6 7 8 9 10"},
    )
    # And the sanity check the crash case cannot make: both specs are within the
    # 105-residue chain, so neither submit would have raised IndexError.
    for sg in subgroups(out, rows1) + subgroups(out, rows2):
        check_true(
            "same length: spec is in range (so the old bug would have been silent)",
            max(int(x) for x in sg["fixed"].split()) <= 105,
        )


def test_ten_shards_never_collide(tmp: Path) -> None:
    """The production shape: ten separate submits, one shared leaf."""
    root = tmp / "shards"
    out = root / "out"
    all_rows = []
    for i in range(10):
        df = one_design(root, f"s{i}", f"{i + 1}:{i + 5}", nres=120)
        all_rows.append(_submit(fork.build_proteinmpnn_manifest, df, root, out))
    check("10 shards: ten subgroup dirs", len(sorted(out.glob("grp_*"))), 10)
    check("10 shards: ten task files",
          len(sorted((out / "proteinmpnn_tasks").glob("task_*.tsv"))), 10)
    placement = {}
    for rows in all_rows:
        for sg in subgroups(out, rows):
            check(f"10 shards: {sg['grp'].name} holds one design", len(sg["members"]), 1)
            for name in sg["members"]:
                placement[name] = sg["fixed"]
    check(
        "10 shards: every design kept its own spec",
        placement,
        {f"s{i}": " ".join(str(x) for x in range(i + 1, i + 6)) for i in range(10)},
    )


def test_force_rerun_is_idempotent(tmp: Path) -> None:
    """A genuine re-run of the SAME design does not accumulate a second subgroup."""
    root = tmp / "rerun"
    out = root / "out"
    df = one_design(root, "dA", "1:5", nres=105)
    rows1 = _submit(fork.build_proteinmpnn_manifest, df, root, out)
    before = snapshot(out)
    rows2 = _submit(fork.build_proteinmpnn_manifest, df, root, out, force=True)
    check("force rerun: same manifest rows", rows2, rows1)
    check("force rerun: still one subgroup", len(sorted(out.glob("grp_*"))), 1)
    check("force rerun: still one task file",
          len(sorted((out / "proteinmpnn_tasks").glob("task_*.tsv"))), 1)
    check("force rerun: tree unchanged", snapshot(out), before)
    check("force rerun: one staged input",
          len(subgroups(out, rows2)[0]["members"]), 1)

    # Same again WITHOUT --force (a design whose first task failed is still
    # ready, so this is a real path, not a hypothetical).
    rows3 = _submit(fork.build_proteinmpnn_manifest, df, root, out)
    check("plain rerun: same manifest rows", rows3, rows1)
    check("plain rerun: still one subgroup", len(sorted(out.glob("grp_*"))), 1)


def test_changed_signature_for_one_design_is_refused(tmp: Path) -> None:
    """Same design, same leaf, a DIFFERENT spec: refused, because the collector
    would otherwise find two FASTAs for it and keep an arbitrary one."""
    root = tmp / "respec"
    out = root / "out"
    _submit(fork.build_proteinmpnn_manifest, one_design(root, "dA", "1:5", 105), root, out)
    changed = one_design(root, "dA", "20:25", 105)
    check_raises(
        "changed signature refused",
        lambda: _submit(fork.build_proteinmpnn_manifest, changed, root, out),
        "dA",
        "already staged",
        "-l/--dir-label",
        "--force",
    )
    check("changed signature: nothing new staged", len(sorted(out.glob("grp_*"))), 1)

    # With --force, and because that subgroup holds nothing this submit is not
    # redoing, the superseded subgroup is replaced rather than duplicated.
    rows = _submit(fork.build_proteinmpnn_manifest, changed, root, out, force=True)
    check("changed signature + force: one subgroup again",
          len(sorted(out.glob("grp_*"))), 1)
    check("changed signature + force: new spec in effect",
          subgroups(out, rows)[0]["fixed"], "20 21 22 23 24 25")
    check("changed signature + force: one staged input",
          sorted(subgroups(out, rows)[0]["members"]), ["dA"])


def test_force_superset_rerun_replaces_the_old_subgroup(tmp: Path) -> None:
    """Re-running a SUPERSET of an earlier submit (one shared signature)."""
    root = tmp / "superset"
    out = root / "out"
    df2 = make_frame(root).iloc[:2]          # d1, d2 -- no column flags: one signature
    df3 = make_frame(root)                   # d1, d2, d3

    fork.build_proteinmpnn_manifest(make_ctx(df2, make_args(root), out))
    check("superset: first submit one subgroup", len(sorted(out.glob("grp_*"))), 1)

    check_raises(
        "superset without --force refused",
        lambda: fork.build_proteinmpnn_manifest(make_ctx(df3, make_args(root), out)),
        "already staged",
        "--force",
    )
    rows = fork.build_proteinmpnn_manifest(
        make_ctx(df3, make_args(root, force=True), out)
    )
    check("superset + force: the old subgroup was replaced",
          len(sorted(out.glob("grp_*"))), 1)
    check("superset + force: all three staged together",
          sorted(subgroups(out, rows)[0]["members"]), ["d1", "d2", "d3"])


def test_partial_overlap_is_never_silently_deleted(tmp: Path) -> None:
    """--force may replace a stale subgroup only if it holds nothing else."""
    root = tmp / "partial"
    out = root / "out"
    df3 = make_frame(root)                   # d1, d2, d3 in ONE subgroup
    fork.build_proteinmpnn_manifest(make_ctx(df3, make_args(root), out))
    # Now re-submit only d1 and d2 with --force: the stale subgroup also holds
    # d3, whose result this submit is not redoing, so it must NOT be removed.
    check_raises(
        "partial overlap refused even with --force",
        lambda: fork.build_proteinmpnn_manifest(
            make_ctx(df3.iloc[:2], make_args(root, force=True), out)
        ),
        "already staged",
        "not redoing",
    )
    check("partial overlap: the original subgroup survives",
          len(sorted(out.glob("grp_*"))), 1)
    check(
        "partial overlap: with all three designs still staged",
        sorted(q.stem for q in next(out.glob("grp_*")).glob("inputs/*.pdb")),
        ["d1", "d2", "d3"],
    )


def test_moved_input_path_is_repointed(tmp: Path) -> None:
    """Same design, same signature, a new input file: the symlink follows it.

    Member PATHS are out of the directory digest on purpose, so this reuses the
    directory (and overwrites its outputs) instead of orphaning the old one with
    a stale FASTA in it. Upstream's _symlink_member would have left the old
    link in place and silently redesigned the OLD structure.
    """
    root = tmp / "moved"
    out = root / "out"
    df = one_design(root, "dA", "1:5", nres=105)
    rows1 = _submit(fork.build_proteinmpnn_manifest, df, root, out)
    old_target = subgroups(out, rows1)[0]["members"]["dA"]

    moved = root / "inputs2"
    moved.mkdir(parents=True, exist_ok=True)
    _pdb_n(moved / "dA.pdb", 105)
    df2 = df.copy()
    df2.at["dA", "rfdiffusion3_path"] = str(moved / "dA.pdb")
    rows2 = _submit(fork.build_proteinmpnn_manifest, df2, root, out, force=True)

    check("moved input: same subgroup dir", rows2, rows1)
    check("moved input: still one subgroup", len(sorted(out.glob("grp_*"))), 1)
    sg = subgroups(out, rows2)[0]
    check("moved input: one staged input", len(sg["members"]), 1)
    check("moved input: symlink repointed", sg["members"]["dA"], str(moved / "dA.pdb"))
    check_true("moved input: and not the old target", sg["members"]["dA"] != old_target)


def test_names_are_stable_and_distinct(tmp: Path) -> None:
    """The digest depends on signature + member names, and on nothing else."""
    sig_a = ("D", "1 2 3", "", "")
    sig_b = ("D", "1 2 4", "", "")
    out = tmp / "digest"
    m1 = [("d1", Path("/x/d1.pdb"))]
    m1_moved = [("d1", Path("/elsewhere/d1.pdb"))]
    m2 = [("d2", Path("/x/d2.pdb"))]
    check_true(
        "digest: stable across calls",
        fork._subgroup_dir(out, sig_a, m1) == fork._subgroup_dir(out, sig_a, m1),
    )
    check_true(
        "digest: ignores the member's path",
        fork._subgroup_dir(out, sig_a, m1) == fork._subgroup_dir(out, sig_a, m1_moved),
    )
    check_true(
        "digest: differs on the member",
        fork._subgroup_dir(out, sig_a, m1) != fork._subgroup_dir(out, sig_a, m2),
    )
    check_true(
        "digest: differs on the signature",
        fork._subgroup_dir(out, sig_a, m1) != fork._subgroup_dir(out, sig_b, m1),
    )
    check_true(
        "digest: member order does not matter",
        fork._subgroup_dir(out, sig_a, m1 + m2)
        == fork._subgroup_dir(out, sig_a, m2 + m1),
    )
    check_true(
        "digest: name starts with grp_",
        fork._subgroup_dir(out, sig_a, m1).name.startswith("grp_"),
    )


def test_spec_shadows_the_builtin() -> None:
    """spec.py swaps two hooks and nothing else; the tool still shadows by name."""
    tool = TOOL
    builtin = get_builtin("proteinmpnn")
    check("spec: name unchanged (shadows)", tool.name, "proteinmpnn")
    check("spec: action unchanged", tool.action, builtin.action)
    check("spec: default_input_column unchanged",
          tool.default_input_column, builtin.default_input_column)
    check_true("spec: collect_fn is upstream's", tool.collect_fn is builtin.collect_fn)
    check("spec: default_script is upstream's",
          tool.default_script, builtin.default_script)
    check_true(
        "spec: bundled proteinmpnn.sh exists",
        Path(tool.default_script).is_file(),
    )
    check_true(
        "spec: modal_image.py sits beside that script",
        (Path(tool.default_script).parent / "modal_image.py").is_file(),
    )
    check_true("spec: build_manifest_fn replaced",
               tool.build_manifest_fn is fork.build_proteinmpnn_manifest)
    check_true("spec: add_run_args_fn replaced",
               tool.add_run_args_fn is fork.add_run_proteinmpnn_args)


def test_parser_keeps_every_upstream_flag() -> None:
    """The args wrapper ADDS; it does not replace upstream's flags."""
    before = ArgumentParser(add_help=False)
    builtin_run.add_run_proteinmpnn_args(before)
    after = ArgumentParser(add_help=False)
    fork.add_run_proteinmpnn_args(after)

    def options(p: ArgumentParser) -> set[str]:
        return {o for a in p._actions for o in a.option_strings}

    old, new = options(before), options(after)
    check("parser: no upstream flag lost", sorted(old - new), [])
    check(
        "parser: exactly the two new flags",
        sorted(new - old),
        ["--chains-to-design-column", "--fixed-positions-column"],
    )
    defaults = {a.dest: a.default for a in after._actions}
    check("parser: fixed_positions_column default", defaults["fixed_positions_column"], "")
    check("parser: chains_to_design_column default", defaults["chains_to_design_column"], "")


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        test_reuse_not_reimplementation()
        test_no_flag_regression(tmp)
        test_per_design_specs(tmp)
        test_chains_to_design_column(tmp)
        test_both_flags_refused(tmp)
        test_neither_flag_is_full_redesign(tmp)
        test_missing_value_refuses(tmp)
        test_missing_column_refuses(tmp)
        test_requires_a_chain_list(tmp)
        test_tied_positions_still_guarded(tmp)
        test_symmetry_broadcast_applies_to_the_column(tmp)
        test_upstream_collides_across_submits(tmp)
        test_second_submit_isolated_crash_case(tmp)
        test_second_submit_isolated_same_length(tmp)
        test_ten_shards_never_collide(tmp)
        test_force_rerun_is_idempotent(tmp)
        test_changed_signature_for_one_design_is_refused(tmp)
        test_force_superset_rerun_replaces_the_old_subgroup(tmp)
        test_partial_overlap_is_never_silently_deleted(tmp)
        test_moved_input_path_is_repointed(tmp)
        test_names_are_stable_and_distinct(tmp)
        test_spec_shadows_the_builtin()
        test_parser_keeps_every_upstream_flag()

    if FAILURES:
        print(f"FAILED {len(FAILURES)} of {CHECKS} checks:")
        for f in FAILURES:
            print(f"  - {f}")
        return 1
    print(f"OK: {CHECKS} checks passed.")
    print(
        "NOTE: this proves the wiring and the refusal contracts only. It cannot "
        "prove that a real motifrmsd column is in ProteinMPNN numbering, nor "
        "that its chain letter is the one parse_multiple_chains indexes -- see "
        "the module docstring."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
