"""Modal image for hisbond (used by ``--executor modal``).

The worker is pure Python (gemmi + numpy): slim Debian, no GPU, no scipy or biopython.
gemmi is pinned to 0.7.5, whose structure reader MERGES residues repeating a number in
a chain; the worker therefore reads atom records itself so centre indices count
residues in file order. The manifest builder sets ``gpus_per_task = 0``.
"""

import modal

RESOURCES = {"cpu": 1, "memory": "2G", "timeout": "00:20:00"}


def image() -> modal.Image:
    return modal.Image.debian_slim(python_version="3.12").pip_install("gemmi==0.7.5", "numpy")
