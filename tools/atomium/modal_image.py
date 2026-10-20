"""Modal image for atomium (used by ``--executor modal``).

The private repo is cloned at BUILD time (``github-token`` / ``github-username``
Secrets; the Modal executor has no per-tool secrets hook) and baked in, weights
included, so tasks neither re-clone nor mount a Volume. Every PyG binary has a
prebuilt cp312 wheel, hence no gcc.
"""

import modal

# Pinned against rebuild drift; bump deliberately (pure_wo_jit HEAD, 2026-09-25).
ATOMIUM_REPO = "github.com/AndreiSokolovskii/develop_atomium.git"
ATOMIUM_BRANCH = "pure_wo_jit"
ATOMIUM_COMMIT = "def1d6a802552a0954e08104029cc3abefcbb73d"

TORCH_VERSION = "2.11.0"
PYG_WHEELS = f"https://data.pyg.org/whl/torch-{TORCH_VERSION}+cu128.html"

RESOURCES = {"gpu": "L4", "cpu": 8, "memory": "16G", "timeout": "01:00:00"}


def image() -> modal.Image:
    return (
        modal.Image.micromamba(python_version="3.12")
        .apt_install("git")
        .pip_install(
            f"torch=={TORCH_VERSION}",
            index_url="https://download.pytorch.org/whl/cu128",
        )
        .pip_install(
            "torch_geometric",
            "pyg_lib",
            "torch_scatter",
            "torch_sparse",
            # Required: model_lib.py:747 uses torch_geometric.nn.radius_graph.
            "torch_cluster",
            find_links=PYG_WHEELS,
        )
        .run_commands(
            f"git clone -b {ATOMIUM_BRANCH} "
            f"https://$GITHUB_USERNAME:$GITHUB_TOKEN@{ATOMIUM_REPO} /opt/atomium",
            f"git -C /opt/atomium checkout {ATOMIUM_COMMIT}",
            "git -C /opt/atomium rev-parse HEAD > /opt/atomium/COMMIT",
            # Drops the credentialed remote from .git/config.
            "rm -rf /opt/atomium/.git",
            "test -f '/opt/atomium/model_weights/model_noised_=_n05.pt'",
            secrets=[
                modal.Secret.from_name("github-token"),
                modal.Secret.from_name("github-username"),
            ],
        )
        .env({"ATOMIUM": "/opt/atomium"})
    )
