from pathlib import Path

from prosapia.core import Tool

from .collect_laproteina import collect_laproteina
from .run_laproteina import add_run_laproteina_args, build_laproteina_manifest

TOOL = Tool(
    name="laproteina",
    action="create",
    description=(
        "Generate binder backbones with La-Proteina Complexa (one job per task)."
    ),
    default_script=str(Path(__file__).parent / "laproteina.sh"),
    # Root-only: LPC resolves its target by name out of its own registry, not from
    # a table column, so no input column is ever read. The builder rejects -t.
    default_input_column="not applicable",
    build_manifest_fn=build_laproteina_manifest,
    add_run_args_fn=add_run_laproteina_args,
    collect_fn=collect_laproteina,
)
