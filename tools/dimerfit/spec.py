from pathlib import Path

from prosapia.core import Tool

from .collect_dimerfit import collect_dimerfit
from .run_dimerfit import (
    NO_DEFAULT_COLUMN,
    add_run_dimerfit_args,
    build_dimerfit_manifest,
)

TOOL = Tool(
    name="dimerfit",
    action="update",
    description="Place a C2 dock back into the binder/target frame and measure "
    "whether the partner protomer occludes the target binding site and whether the "
    "two protomers can be linked.",
    default_script=str(Path(__file__).parent / "dimerfit.sh"),
    # No honest default: the dock column differs per campaign, and a wrong one fails
    # silently. The builder refuses the run unless -i names a real column.
    default_input_column=NO_DEFAULT_COLUMN,
    build_manifest_fn=build_dimerfit_manifest,
    add_run_args_fn=add_run_dimerfit_args,
    collect_fn=collect_dimerfit,
)
