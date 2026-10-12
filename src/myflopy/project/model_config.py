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
* every other package -- its ``options`` block, plus a LAK's ``implicit`` (MF6's
  ``IMPLICIT``, which the package's class writes and no FloPy dataset holds; ledger 326).

Output-file records (``*_filerecord``) are skipped: they hold model-specific
paths that would always differ and mean nothing for a configuration comparison.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from flopy.mf6.mfbase import FlopyException, MFDataException

from myflopy._logging import get_logger

logger = get_logger(__name__)

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
            except (MFDataException, FlopyException):
                # flopy raises these for a dataset that is declared in the dfn
                # but carries no data in this model. Nothing to compare, so the
                # row is simply absent from the config.
                logger.debug(
                    "no data for %s/%s; leaving it out of the config",
                    block_name, name,
                )
                continue
            if name == "newtonoptions":
                # The one reading of NEWTON (ledger 450): an EMPTY payload writes
                # no option and reads as absent, and a reloaded `[(True,)]` reads as
                # the same "NEWTON" an in-memory `[("NEWTON",)]` does.
                value = _newton_kwarg(value)
            if value is None:
                continue
            rows.append((section, name, _normalize(value)))
    return rows


def _tdis_period_rows(tdis) -> list[tuple]:
    """Expand TDIS perioddata into per-period perlen/nstp/tsmult rows."""

    try:
        data = tdis.perioddata.get_data()
    except (AttributeError, MFDataException):
        # AttributeError: a duck-typed or partially built TDIS with no
        # perioddata at all. MFDataException: flopy has the attribute but no
        # stored data behind it.
        logger.debug("no TDIS perioddata; the config carries no period rows")
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


#: The MODEL-level option names ``mf.gwf`` accepts -- MF6's own GWF name-file options.
#: One-to-one with ``package_api._model_options``, so a reader can carry every one of
#: them and there is no name this can silently fail to notice.
MODEL_OPTION_NAMES = (
    "list",
    "print_input",
    "print_flows",
    "save_flows",
    "newtonoptions",
    "nc_mesh2d_filerecord",
    "nc_structured_filerecord",
    "nc_filerecord",
)


def _newton_kwarg(value):
    """FloPy's ``newtonoptions`` payload as the ``str | list[str]`` ``mf.gwf`` wants.

    FloPy hands this back in more than one shape depending on how the model was
    produced -- ``[(True,)]`` and ``[(True, 'UNDER_RELAXATION')]`` from a file reload,
    ``(('NEWTON', 'UNDER_RELAXATION'),)`` elsewhere -- all nested one level, with the
    leading element being the NEWTON marker itself (``True`` or the word). Flattening
    and dropping that marker leaves the modifiers, which is what ``mf.gwf`` takes.

    ``"NEWTON"`` is returned for the bare case rather than ``[]`` or ``""``: **both of
    those write NOTHING at all** -- measured on all four spellings -- so a reader that
    normalized to the empty form would drop the very option it exists to carry.

    And for the same reason a payload holding NO record -- what ``newtonoptions=""``
    or ``[]`` leaves, ``get_data()`` answering ``[]`` while ``has_data()`` says True --
    returns ``None``: FloPy writes no NEWTON line for it, so the model runs without
    NEWTON, and reporting ``"NEWTON"`` turned it ON in a rebuild (ledger 450).
    One record, even an empty-looking ``[(True,)]`` from a file reload, is NEWTON.
    """

    if value is None or (isinstance(value, (list, tuple)) and len(value) == 0):
        return None
    if isinstance(value, str) and not value.strip():
        return None
    flat: list[str] = []
    stack = [value]
    while stack:
        item = stack.pop(0)
        if isinstance(item, (list, tuple)):
            stack = list(item) + stack
        elif isinstance(item, str):
            flat.append(item)
        elif item is not None and not isinstance(item, bool):
            flat.append(str(item))
    modifiers = [word.lower() for word in flat if word.upper() != "NEWTON"]
    return modifiers or "NEWTON"


def _model_option_kwargs(gwf) -> dict:
    """Every model-level option ``gwf`` SETS, as ``mf.gwf`` keyword arguments.

    Reads the name file's options block, which is not in ``gwf.packagelist`` and so is
    invisible to :func:`_extract_settings`'s package loop (ledger 286). Unlike
    :func:`_options_rows` this keeps the ``*_filerecord`` entries: they are noise for a
    DIFF and content for a REBUILD, and dropping them here would be the same silent
    loss in miniature.
    """

    name_file = getattr(gwf, "name_file", None)
    if name_file is None:
        return {}
    options: dict = {}
    for name in MODEL_OPTION_NAMES:
        dataset = getattr(name_file, name, None)
        if dataset is None:
            continue
        try:
            if hasattr(dataset, "has_data") and not dataset.has_data():
                continue
            value = dataset.get_data() if hasattr(dataset, "get_data") else dataset
        except (MFDataException, FlopyException):
            logger.debug("no data for model option %s; leaving it out", name)
            continue
        if name == "newtonoptions":
            value = _newton_kwarg(value)
        if value is None:
            continue
        options[name] = value
    return options


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
        # The MODEL's own options block -- `NEWTON` lives here, and nothing read it.
        # `gwf.packagelist` holds the model's PACKAGES, and the name file is not one of
        # them, so the loop below is structurally blind to it: two runs differing only
        # in `NEWTON` compared EQUAL through `ref.diff(new)`'s config tier. That is the
        # option whose absence dried 90,408 cells and diverged a real model, and the
        # tool that should have caught it could not see it (ledger 286).
        name_file = getattr(gwf, "name_file", None)
        if name_file is not None:
            rows += _options_rows("model", name_file, ["options"])

        for pkg in getattr(gwf, "packagelist", None) or []:
            ptype = str(getattr(pkg, "package_type", "")).lower()
            section = str(getattr(pkg, "package_name", ptype) or ptype).lower()
            blocks = ["options", "period"] if ptype == "oc" else ["options"]
            rows += _options_rows(section, pkg, blocks)
            # MF6's LAK IMPLICIT is no FloPy dataset: the package's CLASS writes it
            # (ledger 326), so the block walk above cannot see it and two runs differing
            # only in `mf.lak(implicit=True)` compared EQUAL -- ledger 286's NEWTON trap,
            # reopened. Read through the class marker (a LOADED package whose file holds
            # it is re-classed by `keep_implicit`), and like a FloPy keyword option, only
            # when set.
            if getattr(pkg, "writes_lak_implicit", False):
                rows.append((section, "implicit", True))

    return pd.DataFrame(rows, columns=_SETTINGS_COLUMNS)


class ModelConfig:
    """A normalized, comparable view of one model's configuration settings."""

    def __init__(self, settings: pd.DataFrame):
        """Wrap a normalized ``(section, setting, value)`` settings table."""

        self._settings = settings

    @classmethod
    def from_model(cls, model) -> ModelConfig:
        """Extract the configuration settings from a live or loaded model.

        Parameters
        ----------
        model : SimulationBase
            A myflopy model -- ``model.config`` passes itself. Read duck-typed:
            ``model.sim.tdis`` gives the ``tdis`` section, ``model.ims`` (else the
            first IMS in ``model.sim.sim_package_list``, else ``model.sim.ims``)
            the ``ims`` section, and ``model.gwf`` the ``model`` section (its name
            file) plus one section per package, named by the package's NAME in
            lower case. Each of the three is skipped when its attribute is
            missing, so an object exposing none of them -- a bare FloPy
            ``ModflowGwf``, which has no ``.sim``/``.gwf`` -- gives an EMPTY config
            with no error.
        """

        return cls(_extract_settings(model))

    def settings(self, *, section: str | None = None) -> pd.DataFrame:
        """Return the ``(section, setting, value)`` table, optionally filtered.

        Parameters
        ----------
        section : str, optional
            Keep only this section's rows, re-indexed from 0. Matched exactly after
            ``str(section).lower()``, so ``"IMS"`` reads ``ims``. The sections are
            ``tdis``, ``ims``, ``model`` and each package's NAME in lower case
            (``npf``, ``oc``, a second RIV's ``crisp``) -- see :attr:`sections`. A
            name no section carries returns an empty table, not an error. Default
            ``None``: every row. Always a copy; editing it leaves the config as it
            was.
        """

        frame = self._settings.copy()
        if section is not None:
            frame = frame[frame["section"] == str(section).lower()].reset_index(drop=True)
        return frame

    @property
    def sections(self) -> list[str]:
        """Return the config sections present (``tdis``, ``ims``, ``oc``, ...)."""

        return list(dict.fromkeys(self._settings["section"].tolist()))

    def section(self, name: str) -> dict:
        """Return one section's settings as a ``{setting: value}`` mapping.

        Parameters
        ----------
        name : str
            The section, matched as :meth:`settings` matches ``section=``: exactly,
            after ``str(name).lower()`` -- ``tdis``, ``ims``, ``model`` or a
            package's NAME in lower case (:attr:`sections` lists them). A name no
            section carries returns ``{}``, not an error, so a misspelt section
            reads as one with no settings.
        """

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

    @classmethod
    def model_options(cls, gwf) -> dict:
        """Every model-level option ``gwf`` SETS, as ``mf.gwf`` keyword arguments.

        A classmethod rather than a second module-level function so
        :class:`~myflopy.modflow.mf6.simulation.base.SimulationBase` reaches it through
        the import it ALREADY defers. Hoisting that import to module level looks safe --
        this module imports nothing from myflopy but ``_logging`` -- and is not: a bare
        ``import myflopy`` succeeds while the test suite's import order hits a circular
        import through ``myflopy.project.__init__``. Measured, not assumed.

        Parameters
        ----------
        gwf : flopy.mf6.ModflowGwf
            The FloPy flow model -- ``model.gwf``, or any object with a
            ``name_file``. Its name file's options block is read for the eight
            names in ``MODEL_OPTION_NAMES`` (``list``, ``print_input``,
            ``print_flows``, ``save_flows``, ``newtonoptions`` and the three
            NetCDF ``*_filerecord`` entries, which are kept here because a
            rebuild needs them). An option the model does not set is left out,
            while one set to ``False`` comes back as ``False``;
            ``newtonoptions`` comes back as ``"NEWTON"`` or its modifiers
            (``["under_relaxation"]``), the forms ``mf.gwf`` takes, and is left
            out for a model holding an EMPTY ``newtonoptions`` (``""``/``[]``),
            which FloPy writes as no option at all, so NEWTON is off there (ledger
            450). An object with no ``name_file`` gives ``{}``.
        """

        return _model_option_kwargs(gwf)

    @property
    def model(self) -> dict:
        """The MODEL-level settings -- the name file's own options block.

        ``NEWTON`` lives here, not in any package, which is why it went unseen: the
        extractor walked ``gwf.packagelist`` and the name file is not in it (ledger 286).
        """

        return self.section("model")

    def __repr__(self) -> str:
        """Show the config's sections and total setting count."""

        return f"ModelConfig(sections={self.sections!r}, settings={len(self._settings)})"
