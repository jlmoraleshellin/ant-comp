"""Modal image for chainsel (used by ``--executor modal``).

The worker is pure Python -- gemmi to read PDB/mmCIF, select chains and write the
subset back out -- so the image is a slim Debian with that one (pinned) wheel and
nothing else: no numpy, no GPU, and it builds in seconds. The manifest builder sets
``gpus_per_task = 0``.
"""

import modal

RESOURCES = {"cpu": 1, "memory": "4G", "timeout": "00:20:00"}


def image() -> modal.Image:
    # Pinned: the worker's split mode and its tests were measured on gemmi 0.7.5's
    # reader (it merges/deletes residues whose numbers repeat), and the golden
    # byte-identical test is tied to 0.7.5's writer. PyPI has a cp312
    # manylinux_2_27/2_28 x86_64 wheel for 0.7.5. Upgrade deliberately (see the skill).
    return modal.Image.debian_slim(python_version="3.12").pip_install("gemmi==0.7.5")
