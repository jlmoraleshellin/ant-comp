from pathlib import Path

from prosapia.core import Tool

from .collect_dockrmsd import add_collect_dockrmsd_args, collect_dockrmsd
from .run_dockrmsd import add_run_dockrmsd_args, build_dockrmsd_manifest

TOOL = Tool(
    name="dockrmsd",
    action="update",
    description="Docking RMSD: superpose the model's TARGET chains onto the "
    "reference's, apply that one transform to the model BINDER and report the "
    "binder RMSD with no further fitting, minimised over every permutation of the "
    "target chains. Roles are forced by the chain flags, never inferred from "
    "length; binder residues correspond by file-order index (equal length "
    "required), target residues by author numbering.",
    default_script=str(Path(__file__).parent / "dockrmsd.sh"),
    # Rows are selected by --col-a, as for usalign/posrmsd; this sentinel matches
    # no column, so filter_ready leaves the frame unchanged.
    default_input_column="not applicable",
    build_manifest_fn=build_dockrmsd_manifest,
    add_run_args_fn=add_run_dockrmsd_args,
    collect_fn=collect_dockrmsd,
    add_collect_args_fn=add_collect_dockrmsd_args,
)
