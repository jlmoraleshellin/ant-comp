from pathlib import Path

from prosapia.core import Tool

from .collect_atomiumph import collect_atomiumph
from .run_atomiumph import add_run_atomiumph_args, build_atomiumph_manifest

TOOL = Tool(
    name="atomiumph",
    action="create",
    description=(
        "Design pH-SELECTIVE binder sequences with AtomiUM's design_ph.py: pin "
        "titratable binder residues into a protonation state and redesign around them "
        "so the state is favoured in the complex more than in the free binder. Input "
        "must be a complex (binder chains + target chains); one child row per "
        "(lambda, sample) draw, carrying dE/ddE."
    ),
    default_script=str(Path(__file__).parent / "atomiumph.sh"),
    default_input_column="rfdiffusion3_path",
    build_manifest_fn=build_atomiumph_manifest,
    add_run_args_fn=add_run_atomiumph_args,
    collect_fn=collect_atomiumph,
)
