"""Modal image for propka (used by ``--executor modal``).

PROPKA3 is pure Python with no runtime dependencies and the worker parses PDB text
itself, so the image is a slim Debian with one wheel: no GPU, numpy, gemmi or biopython.

The version is PINNED EXACTLY and echoed into propka_version on every row: PROPKA's
parameters ship with the package, so a bump can move every pKa and leaves built on
different versions must not be pooled. 3.5.1 is the release the worker's reading of the
API was written against (``conformations['AVR']``, ``group.pka_value`` / ``buried`` /
``energy_volume`` / ``energy_local`` / ``determinants``; ``.pdb`` input only). The
manifest builder sets ``gpus_per_task = 0``; the timeout is generous because each design
is up to two PROPKA runs over a whole complex.
"""

import modal

RESOURCES = {"cpu": 1, "memory": "4G", "timeout": "02:00:00"}


def image() -> modal.Image:
    return modal.Image.debian_slim(python_version="3.12").pip_install("propka==3.5.1")
