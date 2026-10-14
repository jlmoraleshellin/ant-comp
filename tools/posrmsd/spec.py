from pathlib import Path

from prosapia.core import Tool

from .collect_posrmsd import add_collect_posrmsd_args, collect_posrmsd
from .run_posrmsd import add_run_posrmsd_args, build_posrmsd_manifest

TOOL = Tool(
    name="posrmsd",
    action="update",
    description="CA RMSD of structure A fitted onto structure B over an explicit "
    "list of 1-indexed file-order chain positions (same residue ordering assumed, "
    "no alignment), with per-segment RMSDs and a residue-identity trust check.",
    default_script=str(Path(__file__).parent / "posrmsd.sh"),
    # Rows are selected by --col-a, as for usalign; this sentinel matches no
    # column, so filter_ready leaves the frame unchanged.
    default_input_column="not applicable",
    build_manifest_fn=build_posrmsd_manifest,
    add_run_args_fn=add_run_posrmsd_args,
    collect_fn=collect_posrmsd,
    add_collect_args_fn=add_collect_posrmsd_args,
)
