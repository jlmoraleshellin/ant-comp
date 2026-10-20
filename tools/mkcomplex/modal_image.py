"""Modal image for mkcomplex (used by ``--executor modal``).

The task script is bash plus one awk one-liner, so the image is bare slim Debian:
nothing installed, no GPU (``gpus_per_task = 0``), built in seconds.
"""

import modal

RESOURCES = {"cpu": 1, "memory": "1G", "timeout": "00:10:00"}


def image() -> modal.Image:
    return modal.Image.debian_slim(python_version="3.12")
