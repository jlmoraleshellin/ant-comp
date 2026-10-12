"""Modal image for atomiumph (used by ``--executor modal``).

Deliberately a SEPARATE image from ``tools/atomium``'s, with its own commit pin, even
though both clone the same private repo. Bumping the pin to pick up work on
``design_ph.py`` must not silently move the validated plain-inverse-folding tool onto
newer upstream code; two files, two pins, two decisions.

The clone happens at BUILD time with the ``github-token`` / ``github-username`` Modal
Secrets. That is not a preference: prosapia's Modal executor hardcodes the task secret
to the run's ``.env`` (``executors/modal.py``) and offers no per-tool ``secrets()``
hook, so a named Secret is reachable only during the build -- and baking the checkout in
means tasks don't re-clone a private repo on every container.

Build-time assertions, both cheap and both there to turn a 15-minute GPU failure into a
15-second build failure:

  * ``design_ph.py`` exists. If the pinned commit predates the pH feature the build says
    so immediately, instead of every task dying at launch with "No such file".
  * **every fine-tuned pH checkpoint this tool accepts** --
    ``ft_weights_for_pH/ft_model_epoch10_{n03,n05,n03_v2,n05_v2}.pt`` -- and
    ``v2_protonation_vocab.py``. All four are asserted because all four are selectable
    with ``--model-noise``; asserting only the pair the first run used would let a
    ``_v2`` request die inside a GPU container.
  * ``design_ph.py`` mentions ``protect_titratable``. The tool always passes that flag,
    and a pin that does not know it would make every task die in argparse -- after the
    image pull and the checkpoint load. A grep costs a second.

The build also dumps upstream's own ``--help`` to ``/opt/atomium_ph/DESIGN_PH_HELP.txt``
and ``atomiumph.sh`` echoes it into every task's ``.out``. That text is the only
authority on what the flags this tool forwards actually DO -- several are documented
here as unverified. The dump is ``|| true``: it imports torch to print a help string,
and a tool whose help cannot be rendered at build time must still be runnable.

WHICH CHECKPOINT (READ FROM SOURCE at 7ffb756, design_ph.py:2041-2042, not yet observed
in a run)::

    model_folder_path = os.path.join(os.path.dirname(file_path), 'ft_weights_for_pH')
    checkpoint_path = os.path.join(
        model_folder_path, f'ft_model_epoch10_{args.checkpoint_model_noise}.pt')

So ``design_ph.py`` loads ``ft_weights_for_pH/ft_model_epoch10_<noise>.pt`` -- the
fine-tuned 30-class (protonation-aware) head -- and **never** touches
``model_weights/model_noised_=_*.pt``, which is the 21-class base that only
``atomium.py`` reads. The earlier build asserted the base weights, i.e. the wrong
directory: it would have passed on an image with no pH checkpoints at all. The task
script now prints the resolved path and its byte size, so the ``.out`` finally records
which file was loaded instead of listing files nobody opened.

``n03_v2`` / ``n05_v2`` are DIFFERENT weights, not newer builds of ``n03``/``n05``.
All four are accepted (``design_ph.py:2132-2134`` lists exactly these), and
``atomiumph_model_noise`` is what keeps rows from two of them out of one ranking.
"""

import modal

# Pinned so a rebuild can't silently pick up new work on an actively developed branch.
# This pin and tools/atomium's are DIFFERENT commits of the same branch and are
# independent decisions: bump this constant deliberately, and only this one.
# Every run records the commit it actually built from as `atomiumph_image_commit`, so
# a table can always say which code produced a row -- and two different values in one
# table mean the rows are not comparable.
ATOMIUM_REPO = "github.com/AndreiSokolovskii/develop_atomium.git"
ATOMIUM_BRANCH = "pure_wo_jit"
ATOMIUM_COMMIT = "7ffb756b0a706266a96ee24bf6237cc551dcb92f"

TORCH_VERSION = "2.11.0"
PYG_WHEELS = f"https://data.pyg.org/whl/torch-{TORCH_VERSION}+cu128.html"

# 2 h, not atomium's 1 h: --redesign-scope coupling costs n_candidates x n_centers
# forward passes per structure. A timeout loses the WHOLE chunk, including structures
# the task had already finished, because outputs only return when the call completes.
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
            # Required and easy to drop: nothing imports it directly, but
            # model_lib.py uses torch_geometric.nn.radius_graph, a torch-cluster
            # binding. Without it every task dies at model construction.
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
            # The checkpoints design_ph.py ACTUALLY loads (see the module docstring).
            # Asserting model_weights/ instead, as this build used to, would pass on an
            # image carrying only the 21-class base weights -- and every task would then
            # die on a missing file after pulling the image.
            "test -f /opt/atomium_ph/ft_weights_for_pH/ft_model_epoch10_n03.pt",
            "test -f /opt/atomium_ph/ft_weights_for_pH/ft_model_epoch10_n05.pt",
            "test -f /opt/atomium_ph/ft_weights_for_pH/ft_model_epoch10_n03_v2.pt",
            "test -f /opt/atomium_ph/ft_weights_for_pH/ft_model_epoch10_n05_v2.pt",
            "test -f /opt/atomium_ph/v2_protonation_vocab.py",
            # The tool ALWAYS passes --protect_titratable. A pin that does not know it
            # would fail in argparse inside every GPU container instead of here.
            "grep -q protect_titratable /opt/atomium_ph/design_ph.py "
            "|| { echo 'atomiumph: the pinned commit of design_ph.py does not know "
            "--protect_titratable, which this tool always passes. Bump ATOMIUM_COMMIT "
            "in tools/atomiumph/modal_image.py, or drop the flag from "
            "_design_ph_extra in run_atomiumph.py.' >&2; exit 1; }",
            # The only authority on what the forwarded flags mean. `|| true` on the
            # dump ALONE: rendering help imports torch, and that must not be able to
            # fail a build whose tasks would have run fine.
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
