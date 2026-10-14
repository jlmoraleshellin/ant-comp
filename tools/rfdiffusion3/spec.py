from pathlib import Path

from prosapia.core import Tool

from .collect_rfdiffusion3 import collect_rfd3
from .run_rfdiffusion3 import add_run_rfd3_args, build_rfd3_manifest

TOOL = Tool(
    name="rfdiffusion3",
    action="create",
    description="Run RFdiffusion3.",
    default_script=str(Path(__file__).parent / "rfdiffusion3.sh"),
    default_input_column="pdb_path",
    build_manifest_fn=build_rfd3_manifest,
    add_run_args_fn=add_run_rfd3_args,
    collect_fn=collect_rfd3,
)
