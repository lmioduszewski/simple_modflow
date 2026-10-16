"""How PEST++ scores a simulated observation: the ONE bound rule and the ONE phi rule.

:class:`~myflopy.modflow.mf6.pest.ies.IesResults` (the review) and
:meth:`PestProject.run_template <myflopy.modflow.mf6.pest.project.PestProject.run_template>`
(one run of the template, ledger 605) both score through here, so a bound and a
weight mean the same thing in a review, in a template run and in PESTPP-IES's own
``<case>.phi.group.csv``:

* :func:`constraint_sides` -- which observations are one-sided bounds, by pyEMU's own
  tag rule (ledger 564);
* :func:`pest_residuals` -- simulated minus measured as PEST++ scores it: a bound
  contributes 0 where it holds;
* :func:`observation_phi` -- each observation's weighted squared residual, beside its
  group, value, weight and both residuals.

pyEMU is imported inside the function that needs it (myflopy imports it only where a
PEST step does).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def constraint_sides(pst) -> pd.Series:
    """Each observation's side as PESTPP-IES scores it: ``"at_most"``, ``"at_least"`` or ``None``.

    The ONE tag rule (ledger 564), and pyEMU's own
    (``pyemu.pst.pst_handler.get_constraint_tags``, which
    ``Pst.less_than_obs_constraints``/``greater_than_obs_constraints`` read): a
    NONZERO-weight observation whose group starts ``l_``, ``less`` (or ``>@``) is
    simulated at most its value; ``g_``, ``greater`` (or ``<@``) at least it.

    Parameters
    ----------
    pst : pyemu.Pst
        The run's control file -- anything carrying its ``observation_data``
        (``obgnme``, ``weight``).

    Returns
    -------
    pandas.Series
        Indexed by observation name, every observation of the control file.
    """

    from pyemu.pst.pst_handler import get_constraint_tags

    obs = pst.observation_data
    group = obs["obgnme"].astype(str).str.lower()
    weighted = pd.to_numeric(obs["weight"], errors="coerce").fillna(0.0) != 0.0
    at_most = (group.str.startswith(get_constraint_tags("lt")) & weighted).to_numpy()
    at_least = (group.str.startswith(get_constraint_tags("gt")) & weighted).to_numpy()
    return pd.Series(
        np.where(at_most, "at_most", np.where(at_least, "at_least", None)),
        index=obs.index,
        dtype=object,
    )


def pest_residuals(pst, simulated: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Simulated minus measured as PEST++ scores it, and raw (ledger 564).

    The ONE place a residual is scored: a bound observation (:func:`constraint_sides`)
    contributes 0 where it holds -- ``"at_most"``: simulated <= measured;
    ``"at_least"``: simulated >= measured -- and its whole residual where it is
    violated; any other is simulated minus measured. Without it a satisfied bound read
    as a misfit: ``phi_contributions`` gave an ``l_`` group satisfied by every
    realization 2,423,185 where PESTPP-IES's own file says 0
    (``docs/reproducers/pest_inequality_targets.py``).

    Parameters
    ----------
    pst : pyemu.Pst
        The run's control file (``obsval``, ``obgnme``, ``weight``).
    simulated : pandas.Series
        Simulated values indexed by observation name (one realization).

    Returns
    -------
    tuple of pandas.Series
        ``(residual, raw)``, both indexed as ``simulated``: PEST++'s residual and
        simulated minus measured.
    """

    measured = pst.observation_data["obsval"].astype(float).reindex(simulated.index)
    raw = simulated.astype(float) - measured
    sides = constraint_sides(pst).reindex(simulated.index)
    residual = raw.copy()
    at_most = (sides == "at_most").to_numpy()
    at_least = (sides == "at_least").to_numpy()
    residual[at_most] = raw[at_most].clip(lower=0.0)
    residual[at_least] = raw[at_least].clip(upper=0.0)
    return residual, raw


def observation_phi(pst, simulated: pd.Series) -> pd.DataFrame:
    """Every observation ``simulated`` holds, scored as PEST++ scores it.

    The ONE phi rule: ``phi = (weight * residual) ** 2`` with :func:`pest_residuals`'
    residual, so a bound holding contributes 0 and a zero-weight row (a forecast, an
    unmeasured period) contributes 0. Summed by ``group`` it is PESTPP-IES's phi per
    group at the control file's weights -- what its ``<case>.phi.group.csv`` writes
    (``docs/reproducers/run_template_is_pestpp_once.py``).

    Parameters
    ----------
    pst : pyemu.Pst
        The control file (``obsval``, ``obgnme``, ``weight``), indexed by observation
        name as ``simulated`` is.
    simulated : pandas.Series
        Simulated values indexed by observation name.

    Returns
    -------
    pandas.DataFrame
        Indexed as ``simulated``: ``group`` (``obgnme``), ``obsval``, ``weight``,
        ``simulated``, ``residual`` (PEST++'s), ``raw_residual`` (simulated minus
        ``obsval``) and ``phi``.
    """

    obs = pst.observation_data
    residual, raw = pest_residuals(pst, simulated)
    weight = pd.to_numeric(obs["weight"], errors="coerce").fillna(0.0).reindex(simulated.index)
    return pd.DataFrame({
        "group": obs["obgnme"].reindex(simulated.index).astype(str),
        "obsval": obs["obsval"].astype(float).reindex(simulated.index),
        "weight": weight.astype(float),
        "simulated": simulated.astype(float),
        "residual": residual.astype(float),
        "raw_residual": raw.astype(float),
        "phi": (weight.astype(float) * residual.astype(float)) ** 2,
    }, index=simulated.index)


__all__ = ["constraint_sides", "observation_phi", "pest_residuals"]
