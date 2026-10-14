from pathlib import Path

from prosapia.core import Tool

from .collect_motifrmsd import collect_motifrmsd
from .run_motifrmsd import add_run_motifrmsd_args, build_motifrmsd_manifest

TOOL = Tool(
    name="motifrmsd",
    action="update",
    description="Motif CA RMSD of an rfd3 motif-scaffolding design against its "
    "reference, over rfd3's own diffused_index_map correspondence.",
    default_script=str(Path(__file__).parent / "motifrmsd.sh"),
    default_input_column="rfdiffusion3_path",
    build_manifest_fn=build_motifrmsd_manifest,
    add_run_args_fn=add_run_motifrmsd_args,
    collect_fn=collect_motifrmsd,
)
