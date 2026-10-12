"""Geostatistical helpers for the first ``myflopy`` PEST slice."""

from __future__ import annotations

import math

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


def checked_nugget(nugget) -> float:
    """Return ``nugget`` as a float, refusing a negative or non-finite one.

    Parameters
    ----------
    nugget : float
        A nugget variance, anything ``float()`` accepts.

    Returns
    -------
    float
        ``float(nugget)``.

    Raises
    ------
    ValueError
        When it is negative, NaN or infinite: a negative nugget would make two
        distinct points more correlated than one point with itself, and the prior's
        covariance could not be drawn (ledger 489).
    """

    value = float(nugget)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(
            f"nugget={nugget!r}: a nugget is a variance and must be a finite number of "
            "at least 0.0 (0.0, the default, is no nugget)."
        )
    return value


def unit_sill(spec: ExpGeoStruct) -> ExpGeoStruct:
    """The same structure with its contribution and nugget scaled to sum to ``1.0``.

    pyEMU's covariance draws (``helpers.geostatistical_draws``, which ``PstFrom.draw``
    and so :meth:`~myflopy.modflow.mf6.pest.project.PestProject.run_ies`/``prior``
    reach, and ``geostatistical_prior_builder``) build the structure's covariance and
    then SET each diagonal entry to the variance the parameter's bounds imply,
    scaling the off-diagonal by the same standard deviations. A structure of sill
    ``contribution + nugget`` therefore lost its nugget there: the diagonal was
    replaced and the correlation of two points stayed ``exp(-h / range)`` whatever
    the nugget (ledger 489; measured 0.718 at ``h = range / 3`` for nuggets 0,
    0.5 and 1). Scaled to a sill of 1 the bound variance is the total and the
    correlation is ``contribution * exp(-h / range) / (contribution + nugget)`` --
    the rule pyEMU's own spectral simulation applies (``SpecSim2d`` rescales
    contribution and nugget to unity).

    Parameters
    ----------
    spec : ExpGeoStruct
        The structure as declared: ``contribution`` the variogram's sill, ``nugget``
        the variance added on top of it at zero separation. ``contribution`` must be
        positive and ``nugget`` finite and not negative.

    Returns
    -------
    ExpGeoStruct
        A copy whose ``contribution`` is ``contribution / (contribution + nugget)``
        and ``nugget`` is ``nugget / (contribution + nugget)``; range, anisotropy,
        bearing and transform unchanged.

    Raises
    ------
    ValueError
        When ``nugget`` is negative or not finite, or ``contribution`` is not
        positive and finite: no correlation between distinct points can then exceed
        the one at zero separation.
    """

    contribution = float(spec.contribution)
    nugget = checked_nugget(spec.nugget)
    if not math.isfinite(contribution) or contribution <= 0.0:
        raise ValueError(
            f"contribution={spec.contribution!r}: the variogram's sill must be a finite "
            "number above 0.0."
        )
    total = contribution + nugget
    return ExpGeoStruct(
        range=spec.range,
        contribution=contribution / total,
        anisotropy=spec.anisotropy,
        bearing=spec.bearing,
        nugget=nugget / total,
        transform=spec.transform,
    )
