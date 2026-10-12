from pathlib import Path

from prosapia.core import Tool

from .collect_dockq import collect_dockq
from .run_dockq import NO_DEFAULT_COLUMN, add_run_dockq_args, build_dockq_manifest

TOOL = Tool(
    name="dockq",
    action="update",
    description=(
        "Did the predicted binder keep the REFERENCE binding mode: DockQ v2 (the "
        "published package) per native binder-target interface -- Fnat, iRMSD, LRMSD, "
        "DockQ and CAPRI class -- against a reference complex column, with the "
        "homo-oligomer chain mapping searched and recorded."
    ),
    default_script=str(Path(__file__).parent / "dockq.sh"),
    # Sentinel, as in hisbond: the model column must be named (-i).
    default_input_column=NO_DEFAULT_COLUMN,
    build_manifest_fn=build_dockq_manifest,
    add_run_args_fn=add_run_dockq_args,
    collect_fn=collect_dockq,
)
