"""Modal image for motifrmsd (used by ``--executor modal``).

The worker is pure Python -- gemmi to read PDB/mmCIF, numpy for the SVD/Kabsch
fit and the RMSDs -- so the image is a slim Debian with those two wheels and
nothing else: no GPU, no torch, no scipy, no biopython, and it builds in seconds.
The manifest builder sets ``gpus_per_task = 0``.
"""

import modal

RESOURCES = {"cpu": 2, "memory": "4G", "timeout": "00:20:00"}


def image() -> modal.Image:
    return modal.Image.debian_slim(python_version="3.12").pip_install("gemmi", "numpy")
