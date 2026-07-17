"""Warned facade: `InterpolatedSurface` moved to `mf6/grid/interpolated_surface`.

Implementation plan 4.3 consolidated the three `surfaces` modules; this legacy
path resolves with a DeprecationWarning (hidden from completion per D12 —
see `docs/deprecation_policy.md`). `rasterize_points` had no importers and was
retired to `attic/rasterize_points.py`.
"""

from __future__ import annotations

from myflopy._deprecation import deprecated_module_getattr

_DEPRECATED = {
    "InterpolatedSurface": (
        "myflopy.modflow.mf6.grid.interpolated_surface:InterpolatedSurface",
        "myflopy.modflow.mf6.grid.interpolated_surface.InterpolatedSurface",
        "0.2.0",
    ),
}
__getattr__, __dir__ = deprecated_module_getattr(_DEPRECATED, __name__)
