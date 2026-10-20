from pathlib import Path

from prosapia.core import Tool

from .collect_atomium import collect_atomium
from .run_atomium import add_run_atomium_args, build_atomium_manifest

TOOL = Tool(
    name="atomium",
    action="create",
    description="Run AtomiUM (ProteinMPNN-like sequence design) on a set of backbones.",
    default_script=str(Path(__file__).parent / "atomium.sh"),
    # rfd3's real column (the bundled proteinmpnn's older default submits nothing).
    default_input_column="rfdiffusion3_path",
    build_manifest_fn=build_atomium_manifest,
    add_run_args_fn=add_run_atomium_args,
    collect_fn=collect_atomium,
)
