"""Normalized configuration view of a model (``model.config``).

Extracts the *settings* of a model -- the scalar option/solver values that
control how a run behaves, as opposed to its cell data -- into one comparable
table. This is the source the :class:`~myflopy.project.model_diff.ModelDiff`
config tier diffs, and it is independently useful for inspecting one model.

What is captured:

* **tdis** -- ``options`` (time units, start datetime) + ``dimensions`` (nper) +
  one row per stress period for ``perlen`` / ``nstp`` / ``tsmult``.
* **ims** -- the full solver block set (``options`` + ``nonlinear`` + ``linear``:
  dvclose/maximum, under-relaxation, backtracking, linear acceleration, etc.).
* **oc** -- ``options`` + ``period`` save/print records.
* every other package -- its ``options`` block.

Output-file records (``*_filerecord``) are skipped: they hold model-specific
paths that would always differ and mean nothing for a configuration comparison.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

_SETTINGS_COLUMNS = ["section", "setting", "value"]


def _normalize(value):
    """Coerce a flopy-returned value into a plain, comparable Python object."""

    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        if value.dtype.names:  # structured / rec array (e.g. saverecord)
            return tuple(tuple(_normalize(x) for x in row) for row in value.tolist())
        return tuple(_normalize(x) for x in value.tolist())
    if isinstance(value, (list, tuple)):
        return tuple(_normalize(x) for x in value)
    if isinstance(value, str):
        return value
    return value


def _options_rows(section: str, pkg, blocks) -> list[tuple]:
    """Yield ``(section, setting, value)`` rows from a package's chosen blocks."""

    rows: list[tuple] = []
    pkg_blocks = getattr(pkg, "blocks", None)
    if not pkg_blocks:
        return rows
    for block_name in blocks:
        block = pkg_blocks.get(block_name)
        if block is None:
            continue
        for name, dataset in getattr(block, "datasets", {}).items():
            if "filerecord" in name:  # output paths -- model-specific noise
                continue
            try:
                value = dataset.get_data()
            except Exception:
                continue
            if value is None:
                continue
            rows.append((section, name, _normalize(value)))
    return rows


def _tdis_period_rows(tdis) -> list[tuple]:
    """Expand TDIS perioddata into per-period perlen/nstp/tsmult rows."""

    try:
        data = tdis.perioddata.get_data()
    except Exception:
        return []
    if data is None:
        return []
    rows: list[tuple] = []
    names = getattr(getattr(data, "dtype", None), "names", None)
    for index, record in enumerate(data):
        if names:
            values = {field: record[field] for field in names}
            perlen = values.get("perlen")
            nstp = values.get("nstp")
            tsmult = values.get("tsmult")
        else:  # positional fallback: (perlen, nstp, tsmult)
            perlen, nstp, tsmult = record[0], record[1], record[2]
        rows.append(("tdis", f"period[{index}].perlen", _normalize(perlen)))
        rows.append(("tdis", f"period[{index}].nstp", _normalize(nstp)))
        rows.append(("tdis", f"period[{index}].tsmult", _normalize(tsmult)))
    return rows


def _ims_package(model):
    """Return the model's IMS solver package, or ``None``."""

    ims = getattr(model, "ims", None)
    if ims is not None:
        return ims
    sim = getattr(model, "sim", None)
    if sim is None:
        return None
    for pkg in getattr(sim, "sim_package_list", None) or []:
        if str(getattr(pkg, "package_type", "")).lower() == "ims":
            return pkg
    return getattr(sim, "ims", None)


def _extract_settings(model) -> pd.DataFrame:
    """Build the ``(section, setting, value)`` settings table for a model."""

    rows: list[tuple] = []
    sim = getattr(model, "sim", None)
    gwf = getattr(model, "gwf", None)

    tdis = getattr(sim, "tdis", None) if sim is not None else None
    if tdis is not None:
        rows += _options_rows("tdis", tdis, ["options", "dimensions"])
        rows += _tdis_period_rows(tdis)

    ims = _ims_package(model)
    if ims is not None:
        rows += _options_rows("ims", ims, ["options", "nonlinear", "linear"])

    if gwf is not None:
        for pkg in getattr(gwf, "packagelist", None) or []:
            ptype = str(getattr(pkg, "package_type", "")).lower()
            section = str(getattr(pkg, "package_name", ptype) or ptype).lower()
            blocks = ["options", "period"] if ptype == "oc" else ["options"]
            rows += _options_rows(section, pkg, blocks)

    return pd.DataFrame(rows, columns=_SETTINGS_COLUMNS)


class ModelConfig:
    """A normalized, comparable view of one model's configuration settings."""

    def __init__(self, settings: pd.DataFrame):
        """Wrap a normalized ``(section, setting, value)`` settings table."""

        self._settings = settings

    @classmethod
    def from_model(cls, model) -> ModelConfig:
        """Extract the configuration settings from a live or loaded model."""

        return cls(_extract_settings(model))

    def settings(self, *, section: str | None = None) -> pd.DataFrame:
        """Return the ``(section, setting, value)`` table, optionally filtered."""

        frame = self._settings.copy()
        if section is not None:
            frame = frame[frame["section"] == str(section).lower()].reset_index(drop=True)
        return frame

    @property
    def sections(self) -> list[str]:
        """Return the config sections present (``tdis``, ``ims``, ``oc``, ...)."""

        return list(dict.fromkeys(self._settings["section"].tolist()))

    def section(self, name: str) -> dict:
        """Return one section's settings as a ``{setting: value}`` mapping."""

        frame = self.settings(section=name)
        return dict(zip(frame["setting"], frame["value"], strict=False))

    @property
    def tdis(self) -> dict:
        """The TDIS (time-discretization) settings as a ``{setting: value}`` mapping."""

        return self.section("tdis")

    @property
    def ims(self) -> dict:
        """The IMS (solver) settings as a ``{setting: value}`` mapping."""

        return self.section("ims")

    @property
    def oc(self) -> dict:
        """The OC (output-control) settings as a ``{setting: value}`` mapping."""

        return self.section("oc")

    def __repr__(self) -> str:
        """Show the config's sections and total setting count."""

        return f"ModelConfig(sections={self.sections!r}, settings={len(self._settings)})"
