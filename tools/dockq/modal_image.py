"""Modal image for dockq (used by ``--executor modal``).

The worker imports the published DockQ v2 package and nothing else of its own. DockQ
2.1.3 ships on PyPI as an SDIST ONLY: pip compiles its Cython extension
(``DockQ.operations``) at install time, so the image needs a C compiler
(``build-essential``; ``python3-dev`` is cheap insurance for the compile). The package's fallback for a missing extension imports
``operations_nocy`` as a top-level module, which fails when DockQ is installed as a
package, so a broken build shows up as an import error on every row, never as
different numbers. The package and its runtime dependencies are pinned to the versions
the tool was measured with (numpy must stay < 2.0: DockQ's own pin). The manifest
builder sets ``gpus_per_task = 0``, so callers never need ``-g 0``.
"""

import modal

RESOURCES = {"cpu": 1, "memory": "2G", "timeout": "00:20:00"}


def image() -> modal.Image:
    return (
        modal.Image.debian_slim(python_version="3.12")
        .apt_install("build-essential", "python3-dev")
        .pip_install(
            "numpy==1.26.4",
            "biopython==1.88",
            "networkx==3.7",
            "parallelbar==2.5",
            "tqdm==4.70.1",
            "colorama==0.4.6",
            "DockQ==2.1.3",
        )
    )
