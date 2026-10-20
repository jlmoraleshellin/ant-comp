from pathlib import Path

from prosapia.core import Tool

from .collect_propka import collect_propka
from .run_propka import (
    NO_DEFAULT_COLUMN,
    add_run_propka_args,
    build_propka_manifest,
)

TOOL = Tool(
    name="propka",
    action="update",
    description=(
        "What pKa does PROPKA3 predict for each designed histidine of the binder in the "
        "COMPLEX, what for the same residue in the FREE BINDER (its chain(s) alone, "
        "same coordinates, rigid), and how big is the shift (dpka = complex - free)? "
        "Centres: a spec (with {expr} islands) or every histidine. --free-only runs "
        "the free state alone."
    ),
    default_script=str(Path(__file__).parent / "propka.sh"),
    # Sentinel, as in hisbond/epitope/cms/chainsel: any structure column is valid, so -i is required.
    default_input_column=NO_DEFAULT_COLUMN,
    build_manifest_fn=build_propka_manifest,
    add_run_args_fn=add_run_propka_args,
    collect_fn=collect_propka,
)
