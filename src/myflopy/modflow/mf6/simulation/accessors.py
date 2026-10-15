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


def _current_reader(model, cls, slot: str, unit: str | None = None):
    """``model``'s reader of one field, kept in ``model.<slot>`` while its file is unchanged.

    Every access re-resolves the file the OC declares (``cls._declared_path``) and
    stats it; the kept reader comes back only when it is that file, unchanged since it
    was indexed (device, inode, size, mtime), on the model's current grid and with the
    ``unit`` asked for -- otherwise a new one replaces it. A re-run into a reused
    directory rewrites the file IN PLACE (same inode, new size and mtime, measured), so
    the stat is what notices. A new reader was built on every access before, and never
    closed: one ``model.hds.map()`` built five, a ResourceWarning each (ledger 598).
    The replaced reader is not closed -- a caller may still hold it; its file closes
    when it is collected, and it refuses to read a rewritten file.
    """

    path = cls._declared_path(model)
    vor = model.vor
    cached = getattr(model, slot, None)
    if (
        isinstance(cached, cls)
        and getattr(cached, "unit", None) == (unit or getattr(cls, "_default_unit", None))
        and cached._is_current(path, vor=vor)
    ):
        return cached
    options = {} if cls is Hp else {"unit": unit}
    reader = cls(path, model=model, vor=vor, **options)
    setattr(model, slot, reader)
    return reader


def get_hds(model):
    """``model``'s ``HeadsPlus``, kept as ``model._hds`` while its head file is unchanged.

    The file is the one the OC declares (``HEAD FILEOUT``); see
    :func:`_current_reader` for when a new reader replaces the kept one (ledger 598).
    Head TARGETS do not come through here -- they read the OC-declared file's
    period-end records through ``target_values.read_target_values`` (ledger 557).
    """

    return _current_reader(model, Hp, "_hds")


def get_conc(model, unit: str | None = None):
    """``model``'s ``ConcResults`` (GWT concentration), kept as ``model._conc``.

    As :func:`get_hds`, for the file ``CONCENTRATION FILEOUT`` declares.
    """

    return _current_reader(model, ConcResults, "_conc", unit)


def get_temp(model, unit: str | None = None):
    """``model``'s ``TempResults`` (GWE temperature), kept as ``model._temp``.

    As :func:`get_hds`, for the file ``TEMPERATURE FILEOUT`` declares.
    """

    return _current_reader(model, TempResults, "_temp", unit)


def get_all_heads(model):
    """Return the full heads table from the model's ``HeadsPlus`` (:func:`get_hds`)."""
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
    """Every saved ``(kstp, kper)`` of the model's own field file, zero-based.

    Read through :func:`field_reader` on every call, so it follows the file: it was
    cached on the model as ``_kstpkper`` and never invalidated (ledger 301), so a clock
    read while MF6 was still writing stayed short for the object's life (ledger 598).
    """

    return field_reader(model).get_kstpkper()
