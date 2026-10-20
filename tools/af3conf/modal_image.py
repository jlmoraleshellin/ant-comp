"""Modal image for af3conf (used by ``--executor modal``).

The worker is pure Python (gemmi for AF3's CIF, numpy for the PAE matrix and ipSAE, JSON
via the standard library): slim Debian, no GPU, no jax, no AlphaFold3, since it reads
files AF3 already wrote. Pinned to boltzconf's versions (gemmi 0.7.5, numpy 2.5.3) so
the shared ipSAE arithmetic runs on one numpy. The manifest builder sets
``gpus_per_task = 0``.
"""

import modal

RESOURCES = {"cpu": 2, "memory": "4G", "timeout": "00:30:00"}


def image() -> modal.Image:
    return modal.Image.debian_slim(python_version="3.12").pip_install(
        "gemmi==0.7.5", "numpy==2.5.3"
    )
