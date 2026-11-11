from pathlib import Path

from prosapia.core import Tool

from .collect_bindcraft2 import add_collect_bindcraft2_args, collect_bindcraft2
from .run_bindcraft2 import add_run_bindcraft2_args, build_bindcraft2_manifest

TOOL = Tool(
    name="bindcraft2",
    action="create",
    description="Run BindCraft2 binder-design campaigns (one campaign per target).",
    default_script=str(Path(__file__).parent / "bindcraft2.sh"),
    default_input_column="pdb_path",
    build_manifest_fn=build_bindcraft2_manifest,
    add_run_args_fn=add_run_bindcraft2_args,
    collect_fn=collect_bindcraft2,
    add_collect_args_fn=add_collect_bindcraft2_args,
)
