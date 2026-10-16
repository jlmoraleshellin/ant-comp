"""Modal image for bindcraft2 (used by ``--executor modal``).

Built the way BindCraft2's own ``containers/Dockerfile`` builds it, for the same two
reasons that file gives:

* **Editable install.** ``settings/`` and ``scaffolds/`` live at the repository root,
  not inside the package, and the runtime resolves them relative to it. A plain
  ``pip install`` copies only the package and orphans them, so the checkout has to
  stay in place and on the import path.
* **ldconfig.** The jax cuda13 wheels keep their CUDA and cuDNN libraries under
  ``nvidia/*/lib`` inside site-packages, where the loader does not look. Without the
  ``ld.so.conf.d`` entry jax silently falls back to the CPU -- it warns once, runs a
  hundred times slower, and the campaign still starts. The ``ldconfig -p | grep
  libcupti`` at the end fails the *build* instead of someone's campaign.

The AlphaFold parameters (~5.3 GB) are NOT baked in: they download on first use into
``$XDG_CACHE_HOME/bindcraft``, which is a Volume
(``SAPIA_MODAL_VOLUME_BINDCRAFT_CACHE``, default ``bindcraft-cache``) so later tasks
reuse them. Warm it with a single one-target campaign before fanning out, or several
first containers will each pull the same 5.3 GB. ``BINDCRAFT_AF2_PARAMS`` in ``.env``
overrides the location entirely.
"""

import modal

from prosapia.core.executors.modal import get_named_volume

# Mount point for the AlphaFold-parameter Volume.
#
# ROOT CAUSE, observed in a build log on 2026-09-30 -- do not re-litigate this.
# Three submits died at container start with
# `cannot mount volume on non-empty path` -- app detached, 0 tasks, no
# container, and therefore NO .out/.err/.exit written at all, because the
# container never starts. The build guard finally printed the directory:
#
#     --- mount point /bc2_cache must not exist at build time ---
#     drwxr-xr-x 4 root root 83  uv
#
# **Setting `XDG_CACHE_HOME` to the mount point in `.env()` is what did it.**
# The builder's own `uv` writes `$XDG_CACHE_HOME/uv` into the image, so the
# mount point is populated before any task runs. Renaming the mount point does
# not help: the directory follows the variable. Two earlier hypotheses were
# tested and are WRONG -- it is not the legacy path being populated (the log
# shows `/bindcraft_cache` `(absent)`), and it is not `rm` leaving overlay
# whiteouts (nothing was ever deleted successfully).
#
# Also wrong, and worth naming because it is the reasoning that produced the
# first bad fix: "`.env()` is the last image step, so no build command sees the
# variable". The builder's own tooling does see it.
#
# So the image carries the path as **BC2_CACHE_DIR**, a name `uv` ignores, and
# `bindcraft2.sh` turns it into `XDG_CACHE_HOME` at RUNTIME. Never put
# `XDG_CACHE_HOME` in this image's `.env()` again.
CACHE_DIR = "/bc2_cache"
LEGACY_CACHE_DIR = "/bindcraft_cache"

BINDCRAFT_REPO = "https://github.com/PacesaLab/BindCraft2.git"
# Pinned to a release so a rebuild can't silently pick up new work on main. Bump
# deliberately: v1.0.3, ce3150f8d1132a0c900d66ad8f6d85ec08c9ad17.
BINDCRAFT_REF = "v1.0.3"

# A campaign is long by construction: it keeps spending trajectories until
# number_of_final_designs are accepted or max_trajectories are used up. The 12 h
# timeout is a guard, not an estimate -- size it with -T/--time per run.
RESOURCES = {"gpu": "A100", "cpu": 8, "memory": "32G", "timeout": "12:00:00"}

# The Dockerfile's ldconfig step, verbatim in spirit: point the loader at the CUDA
# libraries the jax wheels carry inside site-packages.
_LDCONFIG = (
    "python -c \"import nvidia, pathlib; print('\\n'.join(sorted(str(path) "
    "for root in nvidia.__path__ for path in pathlib.Path(root).glob('*/lib'))))\" "
    "> /etc/ld.so.conf.d/bindcraft-cuda.conf "
    "&& ldconfig && ldconfig -p | grep -q libcupti"
)


def image() -> modal.Image:
    return (
        modal.Image.debian_slim(python_version="3.12")
        # build-essential is here for the same reason the Dockerfile keeps it:
        # biotraj (pulled in by biotite) does not publish a wheel everywhere.
        .apt_install("git", "build-essential")
        .run_commands(
            f"git clone --depth 1 --branch {BINDCRAFT_REF} {BINDCRAFT_REPO} /opt/bindcraft",
            "pip install -e '/opt/bindcraft[cuda13]'",
            _LDCONFIG,
            # Fail the build, not the campaign, if the wheels don't line up.
            "python -c \"import jax, bindcraft.proteinmpnn; print('jax', jax.__version__)\"",
            "bindcraft --help > /dev/null",
        )
        .env(
            {
                # Where `bindcraft` looks for the AlphaFold parameters it downloads,
                # carried under a name `uv` does NOT key off. `bindcraft2.sh` turns
                # it into XDG_CACHE_HOME at runtime.
                #
                # DO NOT set XDG_CACHE_HOME here. Doing so is what made three
                # submits die at container start: the builder's own `uv` writes
                # $XDG_CACHE_HOME/uv into the image, leaving the mount point
                # non-empty. See the comment on CACHE_DIR.
                "BC2_CACHE_DIR": CACHE_DIR,
                "NVIDIA_VISIBLE_DEVICES": "all",
                "NVIDIA_DRIVER_CAPABILITIES": "compute,utility",
            }
        )
        # Modal refuses to mount a Volume onto a non-empty path, and a campaign
        # submitted 2026-09-30 died at container start with exactly that:
        # `cannot mount volume on non-empty path: "/bindcraft_cache"`, retried 7
        # times, no .out/.err/.exit because the container never started.
        #
        # Note the `.env()` above is the LAST env step, so none of the build
        # commands ran with XDG_CACHE_HOME set -- they cached into /root/.cache,
        # not here. Something else populates this path and it is NOT yet known
        # what, so the listing below is deliberate: it puts the diagnosis in the
        # build log rather than costing a separate probe container.
        #
        # Clearing is safe under every candidate cause: the AlphaFold parameters
        # are deliberately NOT baked in (see the module docstring) -- they
        # download on first use into the Volume that mounts here. Anything found
        # at this path at build time is therefore litter by construction.
        .run_commands(
            # Diagnostic ONLY: name what actually populates the old mount point,
            # since that was never observed. Nothing mounts here any more, so a
            # non-empty listing is now information rather than a failure.
            f"echo '--- LEGACY {LEGACY_CACHE_DIR} at end of build (diagnostic only) ---'; "
            f"ls -la {LEGACY_CACHE_DIR} 2>/dev/null || echo '(absent)'; "
            f"du -sh {LEGACY_CACHE_DIR} 2>/dev/null || true; "
            # The real guard. The Volume mounts at CACHE_DIR, so the build must
            # leave that path ABSENT -- not merely emptied. Emptying is what
            # failed twice: if a layer created it, a later `rm` may only write
            # overlay whiteouts and a layer-wise check still sees content.
            # Fail the BUILD, not someone's campaign: at container start this
            # error costs a submit, a queue wait and silent retries, with no
            # .out/.err/.exit written at all because the container never starts.
            f"echo '--- mount point {CACHE_DIR} must not exist at build time ---'; "
            f"if [ -e {CACHE_DIR} ]; then "
            f"ls -la {CACHE_DIR}; "
            f"echo 'ERROR: {CACHE_DIR} exists at build time; the Volume cannot mount there'; "
            f"exit 1; "
            f"else echo '(absent, good)'; fi"
        )
    )


def volumes() -> dict[str, modal.Volume]:
    return {
        CACHE_DIR: get_named_volume(
            "SAPIA_MODAL_VOLUME_BINDCRAFT_CACHE", "bindcraft-cache"
        )
    }
