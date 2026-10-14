from pathlib import Path

from prosapia.core import Tool

from .collect_ssprofile import collect_ssprofile
from .run_ssprofile import (
    NO_DEFAULT_COLUMN,
    add_run_ssprofile_args,
    build_ssprofile_manifest,
)

TOOL = Tool(
    name="ssprofile",
    action="update",
    description=(
        "CA radius of gyration (mass-unweighted) and 3-state DSSP helix/strand/loop "
        "fractions of a chain selection, measured on that selection alone."
    ),
    default_script=str(Path(__file__).parent / "ssprofile.sh"),
    # Sentinel, as in cms/chainsel/epitope: no structure column is a defensible
    # default (rfdiffusion3_path, boltz_path, bindcraft2_path, chainsel_path are all
    # normal inputs), so the manifest builder raises unless -i/--input-column is given.
    default_input_column=NO_DEFAULT_COLUMN,
    build_manifest_fn=build_ssprofile_manifest,
    add_run_args_fn=add_run_ssprofile_args,
    collect_fn=collect_ssprofile,
)
