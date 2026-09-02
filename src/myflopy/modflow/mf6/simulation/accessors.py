"""Shared accessor helpers used by ``SimulationBase`` and loaded MF6 runs."""

from __future__ import annotations

import pandas as pd

from myflopy.modflow.mf6.budget import Budget
from myflopy.modflow.mf6.headsplus import ConcResults, TempResults
from myflopy.modflow.mf6.headsplus import HeadsPlus as Hp
from myflopy.modflow.mf6.package_explorer import ModelPackages
from myflopy.modflow.mf6.package_results import ModelBudgetNamespace
from myflopy.modflow.utils.inputs import Inputs
from myflopy.modflow.utils.outputs import LakOutputData, SFROutputData, UzfOutputData


class ModelOutputs:
    """Namespace for package-specific output helper objects.

    This keeps ``SimulationBase`` package access separate from output-helper
    access. For example, ``model.uzf`` can refer to the MF6 UZF package while
    ``model.outputs.uzf`` exposes convenience output helpers such as
    ``ifno_to_cellid``.
    """

    def __init__(self, model):
        """Bind the package-output accessor namespace (``.lak`` / ``.sfr`` / ``.uzf``) to ``model``."""

        self.model = model

    @property
    def lak(self) -> LakOutputData:
        """Return the lake output helper."""

        return get_lak_output(self.model)

    @property
    def uzf(self) -> UzfOutputData:
        """Return the UZF output helper."""

        return get_uzf_output(self.model)

    @property
    def sfr(self) -> SFROutputData:
        """Return the SFR output helper."""

        return get_sfr_output(self.model)


def get_hds(model):
    """Build and cache the ``HeadsPlus`` helper for ``model``."""
    model._hds = Hp(model=model, vor=model.vor)
    return model._hds


def get_conc(model, unit: str | None = None):
    """Build and cache the ``ConcResults`` (GWT concentration) helper for ``model``."""
    model._conc = ConcResults(model=model, vor=model.vor, unit=unit)
    return model._conc


def get_temp(model, unit: str | None = None):
    """Build and cache the ``TempResults`` (GWE temperature) helper for ``model``."""
    model._temp = TempResults(model=model, vor=model.vor, unit=unit)
    return model._temp


def get_all_heads(model):
    """Return the full heads table from the cached ``HeadsPlus`` helper."""
    return get_hds(model).all_heads


def get_all_conc(model):
    """Return the full concentration table from the cached ``ConcResults`` helper."""
    return get_conc(model).all_conc


def get_all_temp(model):
    """Return the full temperature table from the cached ``TempResults`` helper."""
    return get_temp(model).all_temp


def get_inputs(model):
    """Return an ``Inputs`` helper for inspecting model inputs."""
    return Inputs(model)


def get_lak_output(model):
    """Return the lake-package output helper."""
    return LakOutputData(model)


def get_uzf_output(model):
    """Return the UZF-package output helper."""
    return UzfOutputData(model)


def get_sfr_output(model):
    """Return the SFR-package output helper."""
    return SFROutputData(model)


def get_outputs(model):
    """Return the grouped output-helper namespace for a model."""

    return ModelOutputs(model)


def get_packages(model):
    """Return the preferred package-exploration namespace for a model."""

    return ModelPackages(model)


def get_budget(model, package: str = None):
    """Return a budget helper, optionally scoped to one package."""
    if package is None:
        return Budget(model)
    return Budget(model, package)


def get_model_budget(model):
    """Return the model-budget term namespace behind ``model.budget.<term>``."""
    return ModelBudgetNamespace(model)


def get_budget_cumulative(model):
    """Return the cumulative listing budget as a DataFrame."""
    return pd.DataFrame(model.gwf.output.list().get_cumulative())


def get_budget_incremental(model):
    """Return the incremental listing budget as a DataFrame."""
    return pd.DataFrame(model.gwf.output.list().get_incremental())


def field_reader(model):
    """The dependent-variable reader for this model's kind (heads/conc/temp).

    Lets kind-agnostic machinery (``kstpkper``, choropleth time axis) read the
    output field without hardcoding heads: GWT -> concentration, GWE ->
    temperature, everything else -> heads.
    """

    model_type = getattr(model, "model_type", "gwf6")
    if model_type == "gwt6":
        return get_conc(model)
    if model_type == "gwe6":
        return get_temp(model)
    return get_hds(model)


def get_kstpkper(model):
    """Return and cache available ``(kstp, kper)`` combinations."""
    if model._kstpkper is None:
        model._kstpkper = field_reader(model).kstpkper
    return model._kstpkper
