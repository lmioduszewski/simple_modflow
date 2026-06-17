"""Compatibility imports for the package-level canonical model builder."""

from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    MasterExampleConfig,
    build_canonical_model,
    build_transient_model,
    representative_cells,
)

__all__ = [
    "CanonicalModelConfig",
    "MasterExampleConfig",
    "build_canonical_model",
    "build_transient_model",
    "representative_cells",
]
