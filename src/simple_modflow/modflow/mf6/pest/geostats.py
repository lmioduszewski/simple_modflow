"""Geostatistical helpers for the first ``simple_modflow`` PEST slice."""

from __future__ import annotations

from simple_modflow.modflow.mf6.pest.specs import ExpGeoStruct


def _import_pyemu():
    """Import ``pyemu`` lazily with a clear error message."""

    try:
        import pyemu
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "pyemu is required for simple_modflow PEST workflows. Install pyemu "
            "in the active Python environment before building a PestProject."
        ) from exc
    return pyemu


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
