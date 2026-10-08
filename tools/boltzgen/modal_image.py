"""Modal image for boltzgen (used by ``--executor modal``).

THE PIN, AND WHY IT IS A GIT SHA AND NOT PyPI
---------------------------------------------
``BOLTZGEN_REF`` is ``a07fea1153111ef453aaea7a15697ee05de9291f`` (main, 2026-10-06,
"Merge PR #286 ... boltzgen-numpy-compat"). **Not PyPI.** PyPI's latest is 0.3.2
uploaded 2026-05-01; ``main``'s ``pyproject.toml`` *also* says ``version = "0.3.2"``
while the content differs by five months. The published wheel contains no ``esmfold2``
module and no ``esmfold2_scoring`` step, so the README in the repo does not describe
``pip install boltzgen``, and the ``liability_modality=antibody`` fix for
``nanobody-anything`` (upstream issue #228) is in git and not in the release. A
version string cannot distinguish them; only a SHA can.

THE TRANSITIVE CLOSURE IS PINNED TOO, in ``boltzgen_requirements.txt``
---------------------------------------------------------------------
This is the same failure class recorded in ``tools/rfdiffusion3/modal_image.py``,
where an unbounded ``atomworks[ml]>=2.1.1`` let a freshly built image pick up a new
major release that had removed a function rc-foundry imports, and every task died at
import on an image that had worked for weeks. BoltzGen's ``pyproject.toml`` leaves
~20 dependencies with no upper bound, ``torch>=2.4.1`` among them.

So the image installs an exact, fully resolved requirement set and only then the
package itself with ``--no-deps``. Generated 2026-10-08 with ``uv pip compile
--exclude-newer 2026-10-08``; the file carries its own provenance header.

**torch is forced to 2.10.0.** BoltzGen depends on ``cuequivariance_ops_cu12`` and
``cuequivariance_ops_torch_cu12`` -- CUDA-**12** builds -- and its own ``Dockerfile``
is ``nvidia/cuda:12.2.2-cudnn8`` with a cu121 index. An unconstrained resolve on
2026-10-08 gives torch 2.14.1, whose PyPI wheel carries the **CUDA 13** runtime
(``nvidia-cudnn-cu13``, ``nvidia-nccl-cu13``), i.e. a cu13 torch and cu12
cuequivariance kernels in one process. Probing every release: 2.9.1 and 2.10.0 pull
``nvidia-cudnn-cu12``, 2.11.0 and later pull cu13. 2.10.0 is therefore the newest
version that keeps the whole stack on CUDA 12, and it satisfies ``torch>=2.4.1``.
Bump the SHA and the requirements file together, never one alone.

Other resolved versions worth knowing before reading a traceback: numpy 2.4.6 (the
top of BoltzGen's own ``<2.5``), numba 0.65.1, pandas **3.0.6** (a major release,
2026-01-21; BoltzGen's ``pandas`` dependency has no upper bound at all, so this is
what any current install gets), gemmi 0.7.5, biotite 1.7.1.

**None of this has been executed.** The pin is reproducible, not verified. The first
real run on it is the verification run.

WEIGHTS
-------
~6.4 GB of checkpoints (``boltzgen1_diverse``, ``boltzgen1_adherence``,
``boltzgen1_ifold``, ``boltz2_conf_final``, SolubleMPNN, ``mols.zip``) are NOT baked
in: they download on first use into the Volume mounted at ``/boltzgen_cache``
(``SAPIA_MODAL_VOLUME_BOLTZGEN_CACHE``, default ``boltzgen-cache``), following the
``boltz-cache`` pattern. The Volume is mounted **read-write** deliberately -- upstream
documents no offline/read-only mode, exposes no ``HF_HUB_OFFLINE`` flag, and its own
Dockerfile chowns the cache to the runtime user, so write access is what it expects.

Warm it with a single small run before fanning out, or several first containers each
pull the same 6.4 GB. Upstream issue #216 reports ``mols.zip`` hanging at 0% behind a
TLS-intercepting proxy -- the same failure mode this network has with Modal blob
storage. If that happens the fix is to stage the Volume once from a container, not to
weaken TLS.

THE MOUNT POINT MUST NOT EXIST AT BUILD TIME, and the cache path is carried as
``BOLTZGEN_CACHE`` rather than ``HF_HOME``/``XDG_CACHE_HOME`` on purpose: a cache
variable set in ``.env()`` is seen by the image builder's own tooling, which then
writes into the mount point and makes Modal refuse to mount the Volume -- the
container never starts and no ``.out``/``.err``/``.exit`` is written at all. That cost
three submits on bindcraft2 (see ``tools/bindcraft2/modal_image.py``). ``boltzgen.sh``
turns ``BOLTZGEN_CACHE`` into the ``--cache`` flag at runtime, and the build ends with
a guard that fails the BUILD if the path exists.
"""

from pathlib import Path

import modal

from prosapia.core.executors.modal import get_named_volume

BOLTZGEN_REPO = "https://github.com/HannesStark/boltzgen.git"
# main @ 2026-10-06. Bump together with boltzgen_requirements.txt, never alone.
BOLTZGEN_REF = "a07fea1153111ef453aaea7a15697ee05de9291f"

# The checkout stays on disk after install: the five VHH scaffold YAMLs under
# example/nanobody_scaffolds/ reference their own .cif RELATIVELY, so they have to be
# pointed at in place. run_boltzgen.py's --scaffold-dir default is SCAFFOLD_DIR below.
SRC_DIR = "/opt/boltzgen"
SCAFFOLD_DIR = f"{SRC_DIR}/example/nanobody_scaffolds"
REQUIREMENTS = "/opt/boltzgen_requirements.txt"

CACHE_DIR = "/boltzgen_cache"

# A nanobody run is long by construction: --num-designs designs through four GPU
# steps. The 12 h timeout is a guard, not an estimate -- per-design runtime is
# UNDOCUMENTED upstream (the only figure is an A100 bar chart image). Size it with
# -T/--time per run once the verification run has measured it. H100 per the campaign
# decision; override per run with the base flag --modal-gpu.
RESOURCES = {"gpu": "H100", "cpu": 8, "memory": "64G", "timeout": "12:00:00"}


def image() -> modal.Image:
    return (
        modal.Image.debian_slim(python_version="3.12")
        # build-essential: several of the pinned deps (biotraj via biotite,
        # pdbeccdutils) do not publish a wheel for every platform.
        .apt_install("git", "build-essential")
        # copy=True so the following run_commands can see the file at build time.
        .add_local_file(
            str(Path(__file__).parent / "boltzgen_requirements.txt"),
            REQUIREMENTS,
            copy=True,
        )
        .run_commands(
            f"git clone --filter=blob:none {BOLTZGEN_REPO} {SRC_DIR}",
            f"git -C {SRC_DIR} checkout --detach {BOLTZGEN_REF}",
            # Deps first, exactly pinned; then the package with --no-deps so pip
            # cannot quietly re-resolve anything around them.
            f"pip install --no-cache-dir -r {REQUIREMENTS}",
            f"pip install --no-cache-dir --no-deps {SRC_DIR}",
            # Fail the BUILD, not someone's 12-hour run, if the pins do not line up.
            "python -c \"import torch, boltzgen, cuequivariance_torch; "
            "print('torch', torch.__version__)\"",
            "boltzgen --help > /dev/null",
            # The scaffolds the default --scaffolds names must actually be here.
            f"test -f {SCAFFOLD_DIR}/7eow.yaml && test -f {SCAFFOLD_DIR}/sonelokimab.yaml",
        )
        .env(
            {
                # Read by boltzgen.sh and turned into --cache at RUNTIME. Deliberately
                # NOT named HF_HOME or XDG_CACHE_HOME: see the module docstring.
                "BOLTZGEN_CACHE": CACHE_DIR,
                # run_boltzgen.py's --scaffold-dir default reads this when it happens
                # to be set; inside the container it documents where they live.
                "BOLTZGEN_SCAFFOLD_DIR": SCAFFOLD_DIR,
                "NVIDIA_VISIBLE_DEVICES": "all",
                "NVIDIA_DRIVER_CAPABILITIES": "compute,utility",
            }
        )
        .run_commands(
            # The Volume mounts at CACHE_DIR, so the build must leave that path
            # ABSENT -- not merely emptied. Emptying is what failed twice on
            # bindcraft2: if a layer created it, a later `rm` may only write overlay
            # whiteouts and a layer-wise check still sees content.
            f"echo '--- mount point {CACHE_DIR} must not exist at build time ---'; "
            f"if [ -e {CACHE_DIR} ]; then ls -la {CACHE_DIR}; "
            f"echo 'ERROR: {CACHE_DIR} exists at build time; the Volume cannot mount "
            f"there and the container will never start'; exit 1; "
            f"else echo '(absent, good)'; fi"
        )
    )


def volumes() -> dict[str, modal.Volume]:
    return {
        CACHE_DIR: get_named_volume(
            "SAPIA_MODAL_VOLUME_BOLTZGEN_CACHE", "boltzgen-cache"
        )
    }
