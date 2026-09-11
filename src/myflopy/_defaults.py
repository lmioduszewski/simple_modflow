"""Defaults more than one layer needs, each defined once.

EPSG:2927 (NAD83(HARN) / Washington South, ftUS) is the default CRS everywhere.
Location-agnostic models are placed at its origin on purpose -- they plot in the
Pacific, which gives a map a clean blue background -- so a defaulted CRS is never
refused or warned about. It was spelled 19 ways at 19 sites -- ``"EPSG:2927"``,
``2927``, ``"2927"`` -- until ledger 225; ``tests/test_default_crs.py`` keeps it here.
"""

from __future__ import annotations

DEFAULT_EPSG = 2927
DEFAULT_CRS = f"EPSG:{DEFAULT_EPSG}"
