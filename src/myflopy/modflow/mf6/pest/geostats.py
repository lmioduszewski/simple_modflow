"""Geostatistical helpers for the first ``myflopy`` PEST slice."""

from __future__ import annotations

from myflopy.modflow.mf6.pest.specs import ExpGeoStruct


def _import_pyemu():
    """Import ``pyemu`` lazily with a clear error message."""

    from myflopy._optional import require

    return require("pyemu", feature="myflopy PEST workflows")


def build_geostruct(spec: ExpGeoStruct):
    """Convert an :class:`ExpGeoStruct` spec into a pyEMU ``GeoStruct``."""

    pyemu = _import_pyemu()
    variogram = pyemu.geostats.ExpVario(
        contribution=spec.contribution,
        a=spec.range,
        anisotropy=spec.anisotropy,
        bearing=spec.bearing,
    )
    return pyemu.geostats.GeoStruct(
        variograms=variogram,
        nugget=spec.nugget,
        transform=spec.transform,
    )
