"""Modal image for ssprofile (used by ``--executor modal``).

gemmi to read PDB/mmCIF, numpy for the radius of gyration, and mdtraj for the
secondary-structure assignment. mdtraj is here, rather than a hand-written DSSP,
because it IS the SS engine behind rpxdock's frac_helix/frac_sheet/frac_loop:
willutil's ``wu.dssp`` (what rpxdock calls) is a thin wrapper around
``mdtraj.compute_dssp`` -- see ssprofile_worker.py. Using the same engine keeps the
two tools' fractions on the same basis. willutil itself is NOT installed: it is
git-only and drags in numba, xarray, netcdf4, seaborn, Bio, cppimport and more for a
15-line wrapper.

mdtraj requires scipy, pyparsing and packaging; all ship as wheels, so the image still
builds in well under a minute. mdtraj is PINNED so a rebuild cannot silently change
DSSP behaviour under a column a campaign has already filtered on. NOTE: rpxdock's image
installs willutil from git with mdtraj UNPINNED, so the exact mdtraj version behind
rpxdock's fractions is not recorded anywhere -- parity between the two tools is "same
code path", not "same version", until someone checks.

No GPU. The manifest builder sets ``gpus_per_task = 0``.
"""

import modal

MDTRAJ_PIN = "mdtraj==1.11.0"

RESOURCES = {"cpu": 2, "memory": "4G", "timeout": "00:30:00"}


def image() -> modal.Image:
    return (
        modal.Image.debian_slim(python_version="3.12")
        .pip_install("gemmi", "numpy", MDTRAJ_PIN)
        # Fail the BUILD, not a campaign, if the compiled DSSP kernel is not there.
        .run_commands(
            "python -c 'import gemmi, numpy, mdtraj; "
            "from mdtraj.geometry import _geometry; print(mdtraj.__version__)'"
        )
    )
