"""Modal image for rfdiffusion3 (used by ``--executor modal``).

Checkpoints live on a Volume (``SAPIA_MODAL_VOLUME_RFD3_CKPT``, default
``sapia-rfd3-checkpoints``) mounted at ``/checkpoints``; populate it once with
``foundry install rfd3 --checkpoint-dir /checkpoints`` from a container of this
image. Under Modal the ``RFD3_CKPT`` env var is ignored (see run_rfdiffusion3.py);
to force a specific checkpoint pass ``--ckpt-path /checkpoints/<file>`` explicitly,
since the value is read at submit time and used inside the container.

FORK (ant-comp, 2026-10-07) -- deviation from the prosapia built-in: the install is
PINNED. Upstream installs bare ``rc-foundry[all]``, whose only atomworks constraint is
``atomworks[ml]>=2.1.1`` (no upper bound). atomworks 3.0.0 (PyPI, 2026-10-05) removed
``atomworks.io.utils.bonds.get_inferred_polymer_bonds``, which rc-foundry 0.2.0's
``rfd3/inference/input_parsing.py:14`` imports, so every task on a freshly built image
died at import (diagnostic run ``outputs/20261007_174543_rfd3_modal_check``).

How the pair was chosen (searched PyPI release history 2025-08-15 .. 2026-10-07):
  * atomworks releases: ... 2.2.0 (2025-12-19), 2.2.1 (2026-06-08), 3.0.0 (2026-10-05).
  * rc-foundry releases: ... 0.1.12 (2026-03-09), 0.2.0 (2026-06-05) = latest.
  * Resolving bare ``rc-foundry[all]`` for py3.12/linux with uv ``--exclude-newer
    2026-09-28`` (the day after the last known-good run,
    ``outputs/20260927_211428_rfd3_denovo``) gives exactly rc-foundry 0.2.0 +
    atomworks 2.2.1. Without the cut-off it gives atomworks 3.0.0.
  * The atomworks 2.2.1 wheel's ``atomworks/io/utils/bonds.py`` defines
    ``get_inferred_polymer_bonds`` (line 57, listed in ``__all__``); 3.0.0's does not.
  * A static AST check of every ``from atomworks... import X`` in rc-foundry 0.2.0's
    rfd3/foundry/foundry_cli packages: 190/190 resolve against 2.2.1; against 3.0.0,
    two are missing (that one, and ``atomworks.ml.conditions.base.Level``).
Both pins go in ONE pip_install so the resolver sees them together; rc-foundry's
``>=2.1.1`` admits 2.2.1, so the explicit pin is not overridden (checked with
``uv pip compile``). Other transitive deps remain unpinned (as of 2026-10-07 they
differ from the 09-27 resolution only by patch releases, e.g. torch 2.14.0 -> 2.14.1).
Revisit when an rc-foundry release supports atomworks 3.x.
"""

import modal

from prosapia.core.executors.modal import get_named_volume

CHECKPOINT_DIR = "/checkpoints"
RESOURCES = {"gpu": "A10", "cpu": 8, "memory": "32G", "timeout": "04:00:00"}

# Pinned together: see the module docstring. Bump both at once, never one alone.
RC_FOUNDRY_VERSION = "0.2.0"
ATOMWORKS_VERSION = "2.2.1"


def image() -> modal.Image:
    return (
        modal.Image.debian_slim(python_version="3.12")
        .apt_install("git")
        .pip_install(
            f"rc-foundry[all]=={RC_FOUNDRY_VERSION}",
            f"atomworks[ml]=={ATOMWORKS_VERSION}",
        )
        .env({"FOUNDRY_CHECKPOINT_DIRS": CHECKPOINT_DIR})
    )


def volumes() -> dict[str, modal.Volume]:
    return {CHECKPOINT_DIR: get_named_volume("SAPIA_MODAL_VOLUME_RFD3_CKPT", "sapia-rfd3-checkpoints")}
