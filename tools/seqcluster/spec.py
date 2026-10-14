from pathlib import Path

from prosapia.core import Tool

from .collect_seqcluster import collect_seqcluster
from .run_seqcluster import (
    NO_DEFAULT_COLUMN,
    add_run_seqcluster_args,
    build_seqcluster_manifest,
)

TOOL = Tool(
    name="seqcluster",
    action="update",
    description=(
        "Cluster a table's single-chain design sequences by sequence identity with "
        "MMseqs2 easy-cluster, in ONE task over the whole set, and label every design "
        "with its cluster, representative and cluster size."
    ),
    default_script=str(Path(__file__).parent / "seqcluster.sh"),
    # Sentinel, as in chainsel/ssprofile/epitope: no sequence column is a defensible
    # default (proteinmpnn_sequence, proteinmpnn_redesign_full_sequence, a bindcraft2
    # or boltzgen sequence are all normal inputs), so the manifest builder raises
    # unless -i/--input-column is given.
    default_input_column=NO_DEFAULT_COLUMN,
    build_manifest_fn=build_seqcluster_manifest,
    add_run_args_fn=add_run_seqcluster_args,
    collect_fn=collect_seqcluster,
)
