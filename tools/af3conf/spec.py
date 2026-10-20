from pathlib import Path

from prosapia.core import Tool

from .collect_af3conf import collect_af3conf
from .run_af3conf import (
    NO_DEFAULT_COLUMN,
    add_run_af3conf_args,
    build_af3conf_manifest,
)

TOOL = Tool(
    name="af3conf",
    action="update",
    description=(
        "Binder-specific confidence of an AlphaFold3 prediction, read from the "
        "files AF3 wrote beside the structure: binder pLDDT, the binder's "
        "chain_ptm, the symmetrised binder-target chain_pair_iptm (overall and "
        "per target chain) and ipSAE (Dunbrack's interface-restricted pTM, per "
        "chain pair and per direction) -- the numbers alphafold3_iptm dilutes. "
        "The AF3 counterpart of boltzconf; CPU only, no re-prediction."
    ),
    default_script=str(Path(__file__).parent / "af3conf.sh"),
    # Sentinel, as in boltzconf: several alphafold3 leaves exist, none is a default.
    default_input_column=NO_DEFAULT_COLUMN,
    build_manifest_fn=build_af3conf_manifest,
    add_run_args_fn=add_run_af3conf_args,
    collect_fn=collect_af3conf,
)
