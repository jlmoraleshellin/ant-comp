from pathlib import Path

from prosapia.core import Tool

from .collect_boltzgen import add_collect_boltzgen_args, collect_boltzgen
from .run_boltzgen import add_run_boltzgen_args, build_boltzgen_manifest

TOOL = Tool(
    name="boltzgen",
    action="create",
    description=(
        "Design binders against a protein target with BoltzGen. --protocol "
        "nanobody-anything (default) grafts new CDRs onto a fixed, randomly drawn "
        "VHH framework; --protocol protein-anything designs a de-novo minibinder "
        "chain of a length you choose (--binder-length). One BoltzGen run per "
        "target; the designs are the child rows, and every row records which "
        "protocol made it in boltzgen_protocol."
    ),
    default_script=str(Path(__file__).parent / "boltzgen.sh"),
    # The target a run designs against is a structure, so the default is the generic
    # structure column. Point -i at whatever column actually holds the target in your
    # table (e.g. chainsel_path, mkcomplex_path).
    default_input_column="pdb_path",
    build_manifest_fn=build_boltzgen_manifest,
    add_run_args_fn=add_run_boltzgen_args,
    collect_fn=collect_boltzgen,
    add_collect_args_fn=add_collect_boltzgen_args,
)
