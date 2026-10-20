"""Modal image for boltzconf (used by ``--executor modal``).

The worker is pure Python (gemmi for the structure's atom records, numpy for boltz's npz
arrays and the distance matrix): slim Debian, no GPU, no torch, no boltz, since it reads
files boltz already wrote. Pinned to the versions the tests were measured on (gemmi
0.7.5, numpy 2.5.3) so an upstream change cannot move a column that already gates a
selection. The manifest builder sets ``gpus_per_task = 0``.
"""

import modal

RESOURCES = {"cpu": 2, "memory": "4G", "timeout": "00:30:00"}


def image() -> modal.Image:
    return modal.Image.debian_slim(python_version="3.12").pip_install(
        "gemmi==0.7.5", "numpy==2.5.3"
    )
