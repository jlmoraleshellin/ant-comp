"""Modal image for seqcluster (used by ``--executor modal``).

A slim Debian with ONE static MMseqs2 binary on PATH and nothing else: the worker is
stdlib-only Python, so there is no pip layer at all.

The binary is the official static release ``18-8cc5c``, **sse4.1 build**
(``mmseqs-linux-sse41.tar.gz``), PINNED by URL and sha256 so a rebuild can never
silently swap the clustering engine under a column a campaign has already filtered on.
sse4.1 rather than avx2 deliberately: Modal does not promise a CPU generation, and an
avx2 binary on an older host dies with an illegal instruction. The speed cost is
negligible at the scale this tool is used for (a few thousand short sequences).

The archive is fetched with Python's urllib (debian_slim ships no curl/wget) and
checked with ``sha256sum -c``, which fails the BUILD on any mismatch. The build then
runs ``mmseqs version`` and a two-sequence ``easy-cluster`` as a smoke check, so a
binary that cannot run on the build host fails here, not in a campaign.

sha256 of the asset, verified twice (2026-10-09 by the first build attempt, and
2026-10-10 by re-downloading the release asset into an empty directory):
0ece6a6af8f5d198bccabc98cb5fddfc9212ce04bb332e486bfdcc468c0d7b08.
``mmseqs version`` of that binary prints 8cc5ce367b5638c4306c2d7cfc652dd099a4643f.

No GPU. The manifest builder sets ``gpus_per_task = 0``.
"""

import modal

MMSEQS_RELEASE = "18-8cc5c"
MMSEQS_URL = (
    f"https://github.com/soedinglab/MMseqs2/releases/download/{MMSEQS_RELEASE}/"
    "mmseqs-linux-sse41.tar.gz"
)
MMSEQS_SHA256 = "0ece6a6af8f5d198bccabc98cb5fddfc9212ce04bb332e486bfdcc468c0d7b08"

RESOURCES = {"cpu": 4, "memory": "8G", "timeout": "01:00:00"}


def image() -> modal.Image:
    return (
        modal.Image.debian_slim(python_version="3.12")
        .run_commands(
            "python -c \"import urllib.request; urllib.request.urlretrieve("
            f"'{MMSEQS_URL}', '/tmp/mmseqs.tar.gz')\"",
            f"echo '{MMSEQS_SHA256}  /tmp/mmseqs.tar.gz' | sha256sum -c -",
            "tar -xzf /tmp/mmseqs.tar.gz -C /opt",
            "ln -s /opt/mmseqs/bin/mmseqs /usr/local/bin/mmseqs",
            "rm /tmp/mmseqs.tar.gz",
            # Smoke check: the binary runs, and easy-cluster works end to end.
            "mmseqs version",
            "printf '>s1\\nMKTAYIAKQRQISFVKSHFSRQ\\n>s2\\nMKTAYIAKQRQISFVKSHFSRA\\n' > /tmp/smoke.fa"
            " && mmseqs easy-cluster /tmp/smoke.fa /tmp/smoke /tmp/smoke_tmp --threads 1"
            " >/dev/null && cat /tmp/smoke_cluster.tsv"
            " && rm -rf /tmp/smoke*",
        )
        .env({"MMSEQS": "/usr/local/bin/mmseqs"})
    )
