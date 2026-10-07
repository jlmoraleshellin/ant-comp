from pathlib import Path

from prosapia.core import Tool

from .collect_polyv import collect_polyv
from .run_polyv import NO_DEFAULT_COLUMN, add_run_polyv_args, build_polyv_manifest

TOOL = Tool(
    name="polyv",
    action="update",
    description=(
        "Sequence-blind backbone triage: mutate the binder to poly-valine, repack "
        "only its side chains against a frozen target, and measure interface "
        "geometry (CMS, dSASA, SC) with PyRosetta. Records scores as a floor for "
        "discarding backbones that barely touch the target, not a ranker for picking "
        "good ones. Not an interface score for a designed sequence, and not "
        "comparable with cms or pyrosetta's if_dSASA."
    ),
    default_script=str(Path(__file__).parent / "polyv.sh"),
    # Sentinel, as in cms/usalign/chainsel: no structure column is a defensible
    # default, so the manifest builder raises unless -i/--input-column is given.
    default_input_column=NO_DEFAULT_COLUMN,
    build_manifest_fn=build_polyv_manifest,
    add_run_args_fn=add_run_polyv_args,
    collect_fn=collect_polyv,
)
