"""Modal image for polyv (used by ``--executor modal``).

Deliberately IDENTICAL, layer for layer, to the bundled ``pyrosetta`` tool's image
(``prosapia/tools/pyrosetta/modal_image.py``): same base, same ``pip_install``, same
``run_commands`` string. Modal keys its build cache on those commands, so polyv reuses
the ~1.5 GB PyRosetta wheel this workspace already built for ``pyrosetta`` instead of
downloading it again. **Keep the two in step** -- changing a character here forks the
cache and costs another multi-minute build.

PyRosetta is installed with ``pyrosetta-installer``, which downloads the release wheel
from RosettaCommons at image-build time. Rosetta is single-threaded and this tool is
CPU-only (the manifest builder forces ``gpus_per_task = 0``), but a task holds
``--designs-per-task`` designs, so the timeout is generous.

PyRosetta is free for non-commercial use; commercial use needs a Rosetta license.
"""

import modal

RESOURCES = {"cpu": 2, "memory": "8G", "timeout": "04:00:00"}


def image() -> modal.Image:
    return (
        modal.Image.debian_slim(python_version="3.12")
        .pip_install("pyrosetta-installer")
        .run_commands(
            "python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'"
        )
    )
