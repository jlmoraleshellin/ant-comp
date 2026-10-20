from pathlib import Path

from prosapia.core import Tool

from .collect_chainsel import collect_chainsel
from .run_chainsel import (
    NO_DEFAULT_COLUMN,
    add_run_chainsel_args,
    build_chainsel_manifest,
)

TOOL = Tool(
    name="chainsel",
    action="update",
    description=(
        "Extract a named subset of chains into a new structure file, optionally "
        "merging several source chains into one output chain, or cutting one "
        "source chain into several by residue position (--split-groups)."
    ),
    default_script=str(Path(__file__).parent / "chainsel.sh"),
    # Sentinel, as in usalign: no column is a defensible default here (the tool is
    # equally at home on rfdiffusion3_path, boltz_path or alphafold3_path), so the
    # manifest builder raises unless -i/--input-column is given.
    default_input_column=NO_DEFAULT_COLUMN,
    build_manifest_fn=build_chainsel_manifest,
    add_run_args_fn=add_run_chainsel_args,
    collect_fn=collect_chainsel,
)
