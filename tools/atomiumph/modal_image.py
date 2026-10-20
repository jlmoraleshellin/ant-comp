"""Modal image for atomiumph (used by ``--executor modal``).

A SEPARATE image from ``tools/atomium``'s with its own commit pin, so bumping
``design_ph.py`` cannot move the validated inverse-folding tool. The private repo is
cloned at BUILD time (``github-token`` / ``github-username`` Secrets; the Modal
executor has no per-tool secrets hook) and baked in.

The build asserts what would otherwise fail 15 minutes into a GPU task: ``design_ph.py``
exists, all four ``ft_weights_for_pH/ft_model_epoch10_{n03,n05,n03_v2,n05_v2}.pt``
(the checkpoints it loads, not ``model_weights/``) and ``v2_protonation_vocab.py`` are
present, and ``design_ph.py`` knows ``protect_titratable`` (the tool always passes it).
Upstream's ``--help`` is dumped to ``/opt/atomium_ph/DESIGN_PH_HELP.txt`` (``|| true``,
it imports torch) and echoed into every task's ``.out`` by ``atomiumph.sh``.
"""

import modal

# Pinned against drift, independently of atomium's pin; see atomiumph_image_commit.
ATOMIUM_REPO = "github.com/AndreiSokolovskii/develop_atomium.git"
ATOMIUM_BRANCH = "pure_wo_jit"
ATOMIUM_COMMIT = "7ffb756b0a706266a96ee24bf6237cc551dcb92f"

TORCH_VERSION = "2.11.0"
PYG_WHEELS = f"https://data.pyg.org/whl/torch-{TORCH_VERSION}+cu128.html"

# 2 h (atomium: 1 h): the coupling scan is costly; a timeout loses the whole chunk.
RESOURCES = {"gpu": "L4", "cpu": 8, "memory": "16G", "timeout": "02:00:00"}


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
            # Required, easy to drop: model_lib.py needs torch-cluster (radius_graph).
            "torch_cluster",
            find_links=PYG_WHEELS,
        )
        .run_commands(
            f"git clone -b {ATOMIUM_BRANCH} "
            f"https://$GITHUB_USERNAME:$GITHUB_TOKEN@{ATOMIUM_REPO} /opt/atomium_ph",
            f"git -C /opt/atomium_ph checkout {ATOMIUM_COMMIT}",
            "git -C /opt/atomium_ph rev-parse HEAD > /opt/atomium_ph/COMMIT",
            # Drops the credentialed remote from .git/config.
            "rm -rf /opt/atomium_ph/.git",
            # Fails the build in seconds if the pin predates the pH entry point.
            "test -f /opt/atomium_ph/design_ph.py",
            # What design_ph.py ACTUALLY loads (not the 21-class model_weights/).
            "test -f /opt/atomium_ph/ft_weights_for_pH/ft_model_epoch10_n03.pt",
            "test -f /opt/atomium_ph/ft_weights_for_pH/ft_model_epoch10_n05.pt",
            "test -f /opt/atomium_ph/ft_weights_for_pH/ft_model_epoch10_n03_v2.pt",
            "test -f /opt/atomium_ph/ft_weights_for_pH/ft_model_epoch10_n05_v2.pt",
            "test -f /opt/atomium_ph/v2_protonation_vocab.py",
            # Always passed: fail here, not in argparse inside every GPU container.
            "grep -q protect_titratable /opt/atomium_ph/design_ph.py "
            "|| { echo 'atomiumph: the pinned commit of design_ph.py does not know "
            "--protect_titratable, which this tool always passes. Bump ATOMIUM_COMMIT "
            "in tools/atomiumph/modal_image.py, or drop the flag from "
            "_design_ph_extra in run_atomiumph.py.' >&2; exit 1; }",
            # Upstream --help dump; `|| true`, since rendering it imports torch.
            "(cd /opt/atomium_ph && python design_ph.py --help "
            "> /opt/atomium_ph/DESIGN_PH_HELP.txt 2>&1) || true",
            "touch /opt/atomium_ph/DESIGN_PH_HELP.txt",
            secrets=[
                modal.Secret.from_name("github-token"),
                modal.Secret.from_name("github-username"),
            ],
        )
        .env({"ATOMIUMPH": "/opt/atomium_ph"})
    )
