from pathlib import Path

from prosapia.core import Tool

from .collect_hisbond import collect_hisbond
from .run_hisbond import (
    NO_DEFAULT_COLUMN,
    add_run_hisbond_args,
    build_hisbond_manifest,
)

TOOL = Tool(
    name="hisbond",
    action="update",
    description=(
        "Do the designed His centres hydrogen-bond the TARGET in a predicted "
        "complex, and does the His ACCEPT (pH-sensitive) or DONATE (persistent): "
        "ND1/NE2 heavy-atom distance plus a ring-bisector orientation proxy, with "
        "the role read from the tautomer when the file has hydrogens and from the "
        "partner's chemistry when it does not."
    ),
    default_script=str(Path(__file__).parent / "hisbond.sh"),
    # Sentinel, as in epitope/cms/chainsel: the builder raises unless -i is given.
    default_input_column=NO_DEFAULT_COLUMN,
    build_manifest_fn=build_hisbond_manifest,
    add_run_args_fn=add_run_hisbond_args,
    collect_fn=collect_hisbond,
)
