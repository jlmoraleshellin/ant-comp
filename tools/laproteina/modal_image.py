"""Modal image for laproteina (used by ``--executor modal``).

Built the way LPC's own ``env/build_uv_env.sh --minimal`` builds it, in the same
order and with the same pins, because the ordering is load-bearing: the PyTorch
Geometric wheels are published per torch build, so ``torch==2.7.0+cu126`` has to be
in place before ``-f https://data.pyg.org/whl/torch-2.7.0+cu126.html`` resolves.

``--minimal`` is the right build here. The dependencies it skips -- JAX, ColabFold
and tmol -- are LPC's **reward** models, used for inference-time search and for its
own AF2 refolding. This tool stops at backbones and the campaign validates with
``boltz``, so none of them is ever reached. Skipping them avoids ~5 GB of AlphaFold
parameters and a RoseTTAFold3 checkpoint per container.

**Two things are NOT baked into the image:**

* **The checkpoints** (``complexa.ckpt`` 2.7 GB + ``complexa_ae.ckpt`` 3.8 GB) live
  on a Volume (``SAPIA_MODAL_VOLUME_LPC_CKPT``, default ``lpc-checkpoints``) mounted
  at ``/ckpts``. They are public on NGC and need no credentials, but 6.5 GB per
  container is not something to re-download per task. Populate the Volume once; see
  the laproteina skill.

* **The target registry.** LPC resolves ``--task-name`` against
  ``configs/targets/targets_dict.yaml`` *inside its own checkout*, so a target added
  locally does not exist in a freshly cloned image. The assets Volume
  (``SAPIA_MODAL_VOLUME_LPC_ASSETS``, default ``lpc-assets``) is mounted at
  ``/lpc_assets`` and the task script overlays it onto the checkout before running.
  Put the campaign's ``targets_dict.yaml`` and its ``data/target_data/`` tree there.
"""

import modal

from prosapia.core.executors.modal import get_named_volume

CKPT_DIR = "/ckpts"
ASSETS_DIR = "/lpc_assets"
LPC_HOME = "/opt/Proteina-Complexa"

LPC_REPO = "https://github.com/NVIDIA-BioNeMo/Proteina-Complexa.git"
# The repo's default branch is `dev` and it moves; pin a commit so a rebuild cannot
# silently pick up new work. Bump deliberately.
LPC_REF = "dev"

TORCH_VERSION = "2.7.0+cu126"
TORCH_INDEX = "https://download.pytorch.org/whl/cu126"
PYG_FIND_LINKS = "https://data.pyg.org/whl/torch-2.7.0+cu126.html"

# One generation job of a few binders against a ~500-residue trimeric target. The
# 4 h timeout is a guard, not an estimate -- size it with -T/--time per run once the
# first real run has been measured.
RESOURCES = {"gpu": "A100", "cpu": 8, "memory": "32G", "timeout": "04:00:00"}


def image() -> modal.Image:
    return (
        modal.Image.debian_slim(python_version="3.12")
        .apt_install("git", "build-essential", "wget")
        .run_commands(
            f"git clone --depth 1 --branch {LPC_REF} {LPC_REPO} {LPC_HOME}",
            # Mirrors build_uv_env.sh steps 3-6, minus the --minimal skips.
            f"pip install torch=={TORCH_VERSION} torchvision torchaudio "
            f"--index-url {TORCH_INDEX}",
            f"pip install -e {LPC_HOME}",
            f"pip install torch_geometric torch_scatter torch_sparse torch_cluster "
            f"-f {PYG_FIND_LINKS}",
            "pip install graphein==1.7.7 --no-deps",
            "pip install 'atomworks[ml,openbabel,dev]'",
            "pip install biotite==1.6.0",
            # Fail the build, not someone's run. torch_cluster in particular is easy
            # to lose: nothing imports it by name, but the model's radius_graph call
            # is a torch-cluster binding, so without it every task dies at model
            # construction with a confusing error.
            "python -c \"import torch, torch_cluster, torch_geometric; "
            "print('torch', torch.__version__)\"",
            "complexa --help > /dev/null",
        )
        .env(
            {
                "COMPLEXA_HOME": LPC_HOME,
                "LPC_ASSETS": ASSETS_DIR,
                # LPC's configs reference these through ${oc.env:...}; Hydra fails
                # at config-resolution time if one is missing, so set them even
                # where this tool never reaches the code that reads them.
                "LOCAL_CODE_PATH": LPC_HOME,
                "DATA_PATH": f"{LPC_HOME}/data",
                "COMMUNITY_MODELS_PATH": f"{LPC_HOME}/community_models",
                "NVIDIA_VISIBLE_DEVICES": "all",
                "NVIDIA_DRIVER_CAPABILITIES": "compute,utility",
            }
        )
    )


def volumes() -> dict[str, modal.Volume]:
    return {
        CKPT_DIR: get_named_volume("SAPIA_MODAL_VOLUME_LPC_CKPT", "lpc-checkpoints"),
        ASSETS_DIR: get_named_volume("SAPIA_MODAL_VOLUME_LPC_ASSETS", "lpc-assets"),
    }
