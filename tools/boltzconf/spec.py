from pathlib import Path

from prosapia.core import Tool

from .collect_boltzconf import collect_boltzconf
from .run_boltzconf import (
    NO_DEFAULT_COLUMN,
    add_run_boltzconf_args,
    build_boltzconf_manifest,
)

TOOL = Tool(
    name="boltzconf",
    action="update",
    description=(
        "Binder-specific confidence of a Boltz prediction: binder pLDDT, the "
        "binder's chains_ptm, the symmetrised binder-target pair_chains_iptm "
        "(overall and per target chain), the PAE restricted to the residue "
        "pairs that touch, and ipSAE (Dunbrack's interface-restricted pTM, "
        "per chain pair and per direction) -- the numbers boltz_iptm dilutes."
    ),
    default_script=str(Path(__file__).parent / "boltzconf.sh"),
    # Sentinel, as in epitope/cms/chainsel: several boltz leaves, none a default.
    default_input_column=NO_DEFAULT_COLUMN,
    build_manifest_fn=build_boltzconf_manifest,
    add_run_args_fn=add_run_boltzconf_args,
    collect_fn=collect_boltzconf,
)
