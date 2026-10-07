#!/usr/bin/env python3
"""
Self-test for the two fork deviations of rfdiffusion3 (ant-comp, 2026-10-07).

1. $RFD3_CKPT is ignored when the run's EXECUTOR (-e/--executor, which itself defaults
   to $SAPIA_EXECUTOR) is modal, and honoured under every other executor. The gate
   keys on the parsed executor, never on the env var directly, so a stray
   SAPIA_EXECUTOR must not be able to disagree with an explicit -e.
   The probe builds the REAL `sapia run rfdiffusion3` parser (prosapia's
   build_run_parser + this fork's flags, loaded through tool discovery), so -e and its
   $SAPIA_EXECUTOR default are prosapia's own, not a mimic. RFD3_CKPT is read at module
   import (after load_dotenv(), which picks up this repo's .env -- where RFD3_CKPT is a
   vib cluster path), so every case runs in a fresh subprocess with its own
   environment. load_dotenv() never overrides a variable that is already set, so the
   values passed here win over .env.
2. modal_image.py pins rc-foundry and atomworks together.

What this does NOT establish: that the pinned image builds, that rfd3 imports and runs
inside it, or that rfd3 auto-discovers the checkpoint on the /checkpoints Volume when
ckpt_path is omitted. Those need a real Modal run.

Run:  uv run python tools/rfdiffusion3/test_rfdiffusion3.py
"""

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
VIB_PATH = "/data/groups/csb/anastassia.vorobieva/softwares/foundry/rfd3_latest.ckpt"
CLUSTER_CKPT = "/some/cluster/rfd3.ckpt"


def _load(mod: str):
    spec = importlib.util.spec_from_file_location(f"_rfd3_fork_{mod}", HERE / f"{mod}.py")
    assert spec and spec.loader
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _probe(argv: list[str]) -> None:
    """Child-process mode: parse argv with the real run parser and report the overrides."""
    from prosapia.core.base_run import build_run_parser
    from prosapia.core.tool_registry import BUILTIN_TOOLS_DIR, discover

    tool = discover(BUILTIN_TOOLS_DIR, HERE.parent)["rfdiffusion3"]
    run = sys.modules[tool.build_manifest_fn.__module__]
    parser = build_run_parser(
        tool.metadata, tool.default_script, tool.default_input_column, tool.add_run_args_fn
    )
    args = parser.parse_args(["dummy_run_dir", *argv])
    print(json.dumps({
        "module_file": run.__file__,
        "RFD3_CKPT": run.RFD3_CKPT,
        "executor": args.executor,
        "ckpt_path": args.ckpt_path,
        "overrides": run._cli_overrides(args),
    }))


def _run(env_over: dict[str, str | None], argv: list[str] = ()) -> dict:
    env = dict(os.environ)
    for k, v in env_over.items():
        if v is None:
            env.pop(k, None)
        else:
            env[k] = v
    out = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "--probe", *argv],
        env=env, capture_output=True, text=True,
    )
    if out.returncode != 0:
        raise RuntimeError(f"probe failed for env={env_over} argv={argv}:\n{out.stderr}")
    return json.loads(out.stdout.strip().splitlines()[-1])


FAILS: list[str] = []


def check(label: str, cond: bool, detail: object = "") -> None:
    print(f"{'PASS' if cond else 'FAIL'}  {label}" + (f"   [{detail}]" if not cond else ""))
    if not cond:
        FAILS.append(label)


def main() -> None:
    r = _run({})
    check("probe loads the FORK's run module, not the built-in",
          Path(r["module_file"]).resolve() == HERE / "run_rfdiffusion3.py", r["module_file"])

    # --- the executor actually used decides; SAPIA_EXECUTOR is irrelevant when -e given ---
    for env_ex in ("modal", "slurm", None):
        r = _run({"SAPIA_EXECUTOR": env_ex, "RFD3_CKPT": VIB_PATH}, ["-e", "modal"])
        check(f"-e modal, SAPIA_EXECUTOR={env_ex!r}: RFD3_CKPT suppressed",
              "ckpt_path=" not in r["overrides"], r)
        r = _run({"SAPIA_EXECUTOR": env_ex, "RFD3_CKPT": CLUSTER_CKPT}, ["-e", "slurm"])
        check(f"-e slurm, SAPIA_EXECUTOR={env_ex!r}: RFD3_CKPT honoured",
              f"ckpt_path={CLUSTER_CKPT}" in r["overrides"], r)

    # --- no -e: executor falls back to $SAPIA_EXECUTOR (prosapia's own default) ---
    r = _run({"SAPIA_EXECUTOR": "modal", "RFD3_CKPT": VIB_PATH})
    check("no -e, SAPIA_EXECUTOR=modal: executor=modal, RFD3_CKPT suppressed",
          r["executor"] == "modal" and "ckpt_path=" not in r["overrides"], r)
    for env_ex in ("slurm", None):
        r = _run({"SAPIA_EXECUTOR": env_ex, "RFD3_CKPT": CLUSTER_CKPT})
        check(f"no -e, SAPIA_EXECUTOR={env_ex!r}: RFD3_CKPT honoured",
              r["executor"] == "slurm" and f"ckpt_path={CLUSTER_CKPT}" in r["overrides"], r)
    r = _run({"SAPIA_EXECUTOR": "slurm", "RFD3_CKPT": ""}, ["-e", "slurm"])
    check("-e slurm + empty RFD3_CKPT: no override", "ckpt_path=" not in r["overrides"], r)

    # --- RFD3_CKPT only from the repo .env (not in the process env) ---
    # Control first: off Modal the same .env value IS picked up -- proves load_dotenv()
    # really reads a non-empty RFD3_CKPT here, so the Modal case below is not vacuous.
    # (Skipped if this checkout's .env has no RFD3_CKPT.)
    r_ctl = _run({"SAPIA_EXECUTOR": None, "RFD3_CKPT": None}, ["-e", "slurm"])
    if r_ctl["RFD3_CKPT"]:
        check("control: -e slurm + .env-only RFD3_CKPT is used", "ckpt_path=" in r_ctl["overrides"], r_ctl)
        r = _run({"SAPIA_EXECUTOR": None, "RFD3_CKPT": None}, ["-e", "modal"])
        check("-e modal + .env-only RFD3_CKPT: still ignored", "ckpt_path=" not in r["overrides"], r)
    else:
        print("SKIP  .env-only cases: no RFD3_CKPT in this checkout's .env")

    # --- an explicit --ckpt-path wins on every executor ---
    for ex in ("modal", "slurm"):
        r = _run({"RFD3_CKPT": VIB_PATH}, ["-e", ex, "--ckpt-path", "/checkpoints/x.ckpt"])
        check(f"-e {ex} + explicit --ckpt-path honoured",
              "ckpt_path=/checkpoints/x.ckpt" in r["overrides"]
              and r["overrides"].count("ckpt_path=") == 1, r)

    # --- modal image pins ---
    img = _load("modal_image")
    check("rc-foundry pinned to 0.2.0", img.RC_FOUNDRY_VERSION == "0.2.0", img.RC_FOUNDRY_VERSION)
    check("atomworks pinned to 2.2.1", img.ATOMWORKS_VERSION == "2.2.1", img.ATOMWORKS_VERSION)
    src = (HERE / "modal_image.py").read_text()
    check("image() pins both in one pip_install",
          'f"rc-foundry[all]=={RC_FOUNDRY_VERSION}",\n            f"atomworks[ml]=={ATOMWORKS_VERSION}"' in src)
    check("no unpinned rc-foundry install remains", '"rc-foundry[all]")' not in src)
    img.image()  # builds the lazy Image spec; no network, no Modal auth needed
    check("image() constructs", True)

    print(f"\n{'ALL PASS' if not FAILS else f'{len(FAILS)} FAILED'}")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--probe":
        _probe(sys.argv[2:])
    else:
        main()
