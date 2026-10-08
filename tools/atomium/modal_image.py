"""Modal image for atomium (used by ``--executor modal``).

AtomiUM is a private ProteinMPNN-like sequence designer (a PyG re-implementation).
The repo is cloned at BUILD time with the ``github-token`` / ``github-username``
Modal Secrets, which matters twice over: the prosapia Modal executor hardcodes the
task secret to the run's ``.env`` (``executors/modal.py``) and offers no per-tool
``secrets()`` hook, so a named Secret is only reachable during the build -- and
baking the checkout in means tasks don't re-clone a private repo on every container.

Weights ship IN the repo (``model_weights/model_noised_=_n0X.pt``, ~52 MB x 9, real
files, not LFS pointers), so there is no weights Volume. ``atomium.py`` resolves them
relative to its own ``__file__``, so the task script needs no particular cwd.

All four PyG binaries have prebuilt cp312 linux wheels at the ``find_links`` index
for torch 2.11.0+cu128, so nothing compiles from source (hence no gcc/g++ here).
"""

import modal

# Pinned so a rebuild can't silently pick up new work on an actively developed
# branch. Bump deliberately: HEAD of `pure_wo_jit`, 2026-09-25 "Update add_utils.py".
ATOMIUM_REPO = "github.com/AndreiSokolovskii/develop_atomium.git"
ATOMIUM_BRANCH = "pure_wo_jit"
ATOMIUM_COMMIT = "dee5a4d67ae02293cee59f86c86085c8fa11879d"

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
