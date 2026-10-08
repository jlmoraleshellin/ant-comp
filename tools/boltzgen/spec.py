from pathlib import Path

from prosapia.core import Tool

from .collect_boltzgen import add_collect_boltzgen_args, collect_boltzgen
from .run_boltzgen import add_run_boltzgen_args, build_boltzgen_manifest

TOOL = Tool(
    name="boltzgen",
    action="create",
    description=(
        "Design nanobody (VHH) binders against a protein target with BoltzGen's "
        "nanobody-anything protocol: new CDRs on a fixed, randomly drawn VHH "
        "framework. One BoltzGen run per target; the designs are the child rows."
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
