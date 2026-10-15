"""Dependent-variable output readers built on FloPy's ``HeadFile``.

MODFLOW 6 writes heads (GWF), concentration (GWT), and temperature (GWE) as the
*same* binary layout -- FloPy's ``HeadFile`` reads all three, differing only by
its ``text=`` tag. So the file-agnostic machinery (tidy frames, the ``(kstpkper,
layer, cell)`` table, the unified grammar verbs, choropleth maps) lives once in
:class:`DependentVariableFile`; :class:`HeadsPlus`, :class:`ConcResults`, and
:class:`TempResults` are thin subclasses that set the file suffix / ``text`` tag /
value label and (for heads) add the head-only extras (observations, mounding,
legacy choropleth helpers).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:

    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor

from pathlib import Path

import flopy.utils.binaryfile as bf
import numpy as np
import pandas as pd

from myflopy._numbers import is_no_value
from myflopy._optional import require
from myflopy.modflow.mf6.heads_observations import (
    get_obs_cells as _get_obs_cells,
)
from myflopy.modflow.mf6.heads_observations import (
    get_obs_heads as _get_obs_heads,
)
from myflopy.modflow.mf6.heads_observations import (
    sort_dict_by_keys as _sort_dict_by_keys,
)
from myflopy.modflow.mf6.hydraulics import unit_labels
from myflopy.modflow.mf6.package_explorer_utils import check_view_layers
from myflopy.modflow.mf6.package_plotting import (
    LAYER_FIELD_MAP_PARAMS,
    LineProfileView,
    _apply_backend,
    refuse_noun_parameters,
    resolve_noun_hover,
)
from myflopy.modflow.mf6.period_end import (
    at_period_end,
    last_saved_step,
    period_end_record,
    period_ends,
    period_lengths,
    saved_steps,
)
from myflopy.run_record import nonconverged_steps
from myflopy.modflow.utils.datatypes.datalists import convert_nested_to_int
from myflopy.modflow.utils.datatypes.hover import (
    HoverSpec,
    check_hover_layers,
    conc_hover,
    head_hover,
    temp_hover,
)

idxx = pd.IndexSlice  # for easy index slicing in a MultiIndex DataFrame


#: True where a value is one of MODFLOW's no-value markers (``|value| >= 1e29``). The
#: rule is :func:`myflopy._numbers.is_no_value`; this name is kept for the modules
#: that import it from here.
_is_sentinel = is_no_value


def _as_layer_cell_heads(data, *, nlay: int, ncpl: int):
    """Return dependent-variable data as a consistent ``(nlay, ncpl)`` array."""

    values = np.asarray(data)
    values = np.squeeze(values)

    if values.ndim == 1:
        if nlay == 1 and values.size == ncpl:
            return values.reshape(1, ncpl)
        if values.size == nlay * ncpl:
            return values.reshape(nlay, ncpl)

    if values.ndim == 2:
        if values.shape == (nlay, ncpl):
            return values
        if values.shape == (ncpl, nlay):
            return values.T
        if nlay == 1 and values.size == ncpl:
            return values.reshape(1, ncpl)
        if values.size == nlay * ncpl:
            return values.reshape(nlay, ncpl)

    if values.ndim >= 3:
        if values.shape[0] == nlay and int(np.prod(values.shape[1:])) == ncpl:
            return values.reshape(nlay, ncpl)
        if values.shape[-1] == nlay and int(np.prod(values.shape[:-1])) == ncpl:
            return np.moveaxis(values, -1, 0).reshape(nlay, ncpl)

    if values.size == nlay * ncpl:
        return values.reshape(nlay, ncpl)

    raise ValueError(
        "Could not reshape values to layer/cell form. "
        f"Got shape={values.shape}, expected nlay={nlay}, ncpl={ncpl}."
    )


def _exact_kstpkper(kstpkper, saved) -> tuple[int, int]:
    """``kstpkper`` as a saved ``(kstp, kper)``, refused when it is not two integers or not saved."""

    import numbers

    if isinstance(kstpkper, (str, bytes)) or not hasattr(kstpkper, "__len__") or len(kstpkper) != 2:
        raise TypeError(
            f"kstpkper must be two zero-based integers, (kstp, kper); got {kstpkper!r}."
        )
    for value in kstpkper:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Integral):
            raise TypeError(
                f"kstpkper must be two zero-based integers, (kstp, kper); got {kstpkper!r}."
            )
    key = (int(kstpkper[0]), int(kstpkper[1]))
    available = [(int(item[0]), int(item[1])) for item in saved]
    if key not in available:
        raise ValueError(f"kstpkper {key} is unavailable. Available: {available}")
    return key


def _run_directory(output_path) -> Path | None:
    """The simulation directory holding ``mfsim.lst`` for a field file, or ``None``.

    The file sits in its model's output folder, which is the simulation's own directory
    or, for a model of a coupled run, a subdirectory of it.
    """

    folder = Path(output_path).resolve().parent
    for candidate in (folder, folder.parent):
        if (candidate / "mfsim.lst").is_file():
            return candidate
    return None


def read_saved_field(
    reader,
    *,
    per: int = -1,
    kstpkper: tuple[int, int] | None = None,
    perlen=None,
    ncpl: int,
    label: str = "head",
    run_dir=None,
) -> tuple[np.ndarray, tuple[int, int], str]:
    """The ``(nlay, ncpl)`` numbers one record of a field file holds, and which record.

    The read ``model.hds.remap`` makes, and the lazy recipe after it: the record that ENDS
    stress period ``per`` (:func:`~myflopy.modflow.mf6.period_end.period_end_record`, that
    period alone), or an exact ``kstpkper``. The file's own numbers come back, MODFLOW's
    no-value markers included.

    Parameters
    ----------
    reader : flopy.utils.HeadFile
        An open reader of the field file (a :class:`DependentVariableFile` is one).
    per : int, default -1
        The stress period, zero-based; a negative one counts back from the last.
        Ignored when ``kstpkper`` is given.
    kstpkper : tuple of (int, int), optional
        An exact zero-based saved time instead, unchecked against the period's end.
        Default ``None``.
    perlen : sequence of float, optional
        Each period's TDIS length. Default ``None``: no TDIS, so the period's last saved
        step is taken unchecked, and the description says so.
    ncpl : int
        Cells per layer, for the reshape.
    label : str, default "head"
        The file's field (``"head"``, ``"concentration"``, ``"temperature"``): names the
        file and the output-control word in a refusal.
    run_dir : str or pathlib.Path, optional
        The simulation directory, whose ``mfsim.lst`` names every time step the solution
        did NOT converge on (:func:`~myflopy.run_record.nonconverged_steps`): MF6 saves
        such a step's numbers -- a failed run its last step, one told to ``CONTINUE``
        steps mid-run -- and reading one is refused. Default ``None``: not checked.

    Returns
    -------
    tuple of (numpy.ndarray, tuple of (int, int), str)
        The ``(nlay, ncpl)`` floats, the zero-based ``(kstp, kper)`` read, and a
        description of that time for a message.

    Raises
    ------
    TypeError
        ``per`` or ``kstpkper`` is not made of integers.
    ValueError
        As :func:`~myflopy.modflow.mf6.period_end.period_end_record`; a ``kstpkper`` the
        file did not save; or the record is a step the run did not converge on.
    """

    saved = saved_steps(reader)
    if kstpkper is not None:
        key = _exact_kstpkper(kstpkper, saved)
        when = f"time step {key[0]} of stress period {key[1]} (zero-based)"
    else:
        key = period_end_record(
            saved, per, perlen=perlen, file=f"the {label} file", saves=label.upper()
        )
        when = (
            f"the end of stress period {key[1]}"
            if perlen is not None
            else f"the last saved step of stress period {key[1]} (no TDIS to check that it "
            "ends the period)"
        )
    if run_dir is not None and key in nonconverged_steps(run_dir):
        raise ValueError(
            f"Time step {key[0]} of stress period {key[1]} (zero-based) is a step the "
            f"run did NOT converge on: {Path(run_dir) / 'mfsim.lst'} says so, and MF6 "
            f"still saves such a step's {label}s, at the time it would have ended -- a "
            "failed run writes its last step that way, and one told to CONTINUE writes "
            "them mid-run. Read a period the run solved (per=...), or rerun the model."
        )
    raw = reader.get_data(kstpkper=key)
    return (
        np.asarray(_as_layer_cell_heads(raw, nlay=int(reader.nlay), ncpl=int(ncpl)), dtype=float),
        key,
        when,
    )


class DependentVariableFile(LineProfileView, bf.HeadFile):
    """A per-cell dependent-variable output reader (heads / concentration / temperature).

    Also a full :class:`LineProfileView` leaf in the unified grammar --
    ``.get()/.summary()`` (tables), ``.map()/.plot()/.section()`` (panels),
    ``.mosaic()/.animate()`` (composers). ``LineProfileView`` is first in the MRO
    deliberately so the grammar's ``plot`` shadows flopy's legacy
    ``LayerFile.plot``. Subclasses set the class attributes below.
    """

    #: public value label + tidy-frame column used by the unified grammar
    value_name = "value"
    #: internal column name in the ``(kstpkper, layer, cell)`` table
    store_column = "value"
    #: output-file suffix appended to ``<model name>`` when no path is given
    _output_suffix = ".bin"
    #: FloPy ``HeadFile`` binary ``text=`` selector for this field
    _binary_text = "head"
    #: label used by :meth:`summary`
    _summary_label = "value"
    #: ``Choro`` type key so the map reads this field from the model
    _choro_type = "hds"
    #: Whether :meth:`remap` gives a DRY cell (-1e30) the value of the first wet layer
    #: below it -- the water table. A head's; a dry cell holds no concentration or
    #: temperature at all.
    _dry_takes_water_table = False

    def __init__(
            self,
            path: Path = None,
            model=None,
            vor: Vor = None,
    ):
        """Parameters
        ----------
        path
            Path to a binary MF6 output file. If omitted, ``model`` is used with
            this class's ``_output_suffix``.
        model
            Optional parent model used to infer file paths and geometry helpers.
        vor
            Optional Voronoi/grid helper, defaulting to ``model.vor``.
        """
        from myflopy.modflow.mf6.simulation.base import SimulationBase
        if path is None:
            if model is None:
                raise ValueError("Must provide an output file or model")
            assert isinstance(model, SimulationBase), 'no valid model provided'
            self.output_path = model.model_output_folder_path / f'{model.name}{self._output_suffix}'
        else:
            self.output_path = path

        super().__init__(filename=self.output_path, text=self._binary_text)

        if model is None:
            self.model = None
        elif isinstance(model, SimulationBase):
            self.model = model
        else:
            raise ValueError("model must be an instance of SimulationBase")

        self.kstpkper = convert_nested_to_int(self.get_kstpkper())
        self.vor = self.model.vor if vor is None else vor
        self._value_table = None
        self.nper = pd.DataFrame(self.get_kstpkper()).iloc[:, 1].max() + 1
        self.numstp = pd.DataFrame(self.get_kstpkper()).iloc[:, 0].max() + 1
        self.vor_list = self.vor.gdf_vorPolys.geometry.to_list()
        self.cell_list = [i for i in range(len(self.vor_list))]
        self.area_list = [cell.area for cell in self.vor_list]
        self.x_list = [cell.centroid.xy[0][0] for cell in self.vor_list]
        self.y_list = [cell.centroid.xy[1][0] for cell in self.vor_list]

    @property
    def crs(self):
        """The CRS of the cells this file is read on: the grid's ``crs``, read live.

        It was copied from ``vor.crs`` when the reader was built, so relabelling the
        grid afterwards (``grid.crs = ...``, ledger 368) left the reader on the old
        label (ledger 431). Setting it here is refused: the CRS is the grid's,
        so relabel the grid.
        """

        return self.vor.crs

    @crs.setter
    def crs(self, value):
        raise AttributeError(
            f"{type(self).__name__}.crs = {value!r}: a results reader has no CRS of its "
            "own -- it reads its grid's. Relabel the grid instead: "
            f"grid.crs = {value!r} (this reader's grid is its .vor)."
        )

    def _check_layers(self, layer):
        """``layer`` as a list of ints, refused when the views do not number it.

        ``check_view_layers`` (ledger 496) against the model's views when
        the reader has a model, else against the file's own ``nlay`` -- the same
        count, since MF6 writes a DISU's nodes as one layer.
        """

        gwf = getattr(self.model, "gwf", None) if self.model is not None else None
        return check_view_layers(gwf if gwf is not None else int(self.nlay), layer)

    @property
    def all_values(self):
        """The full ``kstpkper/layer/cell`` value table, built and cached on first access.

        The file's own numbers, MODFLOW's no-value sentinels (+1e30 inactive, -1e30
        dry) included: every view of it -- :meth:`get`, :meth:`long`, :meth:`wide`,
        :meth:`array`, :meth:`get_obs_heads`, ``model.region_heads``, the maps,
        the sections and the group and diff tiers -- reads them as ``NaN``
        (ledger 315; :meth:`get_obs_heads` and ``model.region_heads`` since
        ledger 495). A calibration
        TARGET is not a view: what a head target hands PEST keeps the number
        (ledgers 351, 355).
        """

        if self._value_table is None:
            self._value_table = self._build_value_table()
        return self._value_table

    def _default_hover(self, *, layers: str, surfaces: bool):
        """The default hover spec for this field's :meth:`map` (overridden per kind)."""

        return HoverSpec(
            primary=self.value_name,
            title=self.value_name,
            layers=layers,
            surfaces=surfaces,
            footer=("period", "date"),
        )

    # -- unified grammar: data verbs ---------------------------------------
    def get(
        self,
        *,
        per: int | list[int] | None = None,
        layer: int | list[int] | None = None,
        cells: int | list[int] | None = None,
    ) -> pd.DataFrame:
        """Return period-end values as a tidy ``per``/``layer``/``cell``/``<value>`` frame.

        Each stress period is its END-of-period value: the record that ends it on
        this file's own clock, by the rule every target reads
        (:func:`~myflopy.modflow.mf6.period_end.period_ends`, ledger 336). A period
        whose end output control did not save (``SAVE HEAD FIRST``, ``STEPS`` short
        of it, ``FREQUENCY``) is refused rather than read at its last saved step,
        which it used to be -- step 1 of 3 labelled as the period. A reader opened
        from a path with no model has no TDIS to check against and takes each
        period's last saved step. MODFLOW's no-value sentinels (+1e30 inactive,
        -1e30 dry; ``|value| >= 1e29``) become ``NaN``, by the same rule as
        :meth:`array` and the maps (ledger 315). The three selectors combine with
        AND; floats inside a sequence are truncated to int.

        Parameters
        ----------
        per : int or sequence of int, optional
            Stress period(s) to keep, zero-based. An int (a numpy integer too) or
            any iterable of them: list, tuple, range, ndarray. ``None`` (default)
            keeps every period. A period the file does not hold selects nothing, so
            the frame comes back empty with no error. ``-1`` is NOT the last period
            here, unlike :meth:`array`; it selects nothing.
        layer : int or sequence of int, optional
            Layer(s) to keep, zero-based (0 is the top), in the same forms. ``None``
            (default) keeps every layer. A layer the model's views do not number
            raises ``ValueError`` naming the layer count -- ONE layer on a DISU,
            whatever FloPy infers from ``IHC`` -- where it selected nothing
            (ledger 496).
        cells : int or sequence of int, optional
            Cell index(es) to keep, zero-based within a layer (``0`` to
            ``ncpl - 1``), in the same forms. ``None`` (default) keeps every cell.

        Returns
        -------
        pandas.DataFrame
            Columns ``per``, ``layer`` and ``cell`` (int), and :attr:`value_name`
            (``head``, ``conc`` or ``temp``; float, ``NaN`` where the file holds no
            value). One row per selected period, layer and cell.

        Raises
        ------
        TypeError
            If a selector is a bare float, or any other value that is neither an
            integer nor iterable.
        ValueError
            If ``layer`` names a layer the views do not number; or if a period
            the file saved lacks its final step -- naming the period,
            the steps saved and the ``SAVE HEAD LAST`` (``CONCENTRATION``,
            ``TEMPERATURE``) setting that fixes it. Checked over the whole file,
            whatever ``per`` selects.

        Examples
        --------
        >>> model.hds.get(per=0, layer=0)            # one period, the top layer
        >>> model.hds.get(cells=[5, 7])["head"]      # two cells, every period
        """

        value = self.value_name
        self._check_layers(layer)
        frame = self.all_values.reset_index()
        ends = period_ends(
            self.model,
            reader=lambda: self,
            file=f"the {self._binary_text} file",
            saves=self._binary_text.upper(),
            fallback=frame["kstpkper"],
        )
        frame = at_period_end(frame, ends).copy()
        frame["per"] = [int(key[1]) for key in frame["kstpkper"]]
        frame = frame.rename(columns={self.store_column: value})
        frame[value] = pd.to_numeric(frame[value], errors="coerce")
        frame.loc[_is_sentinel(frame[value]), value] = np.nan
        for column, selector in (("per", per), ("layer", layer), ("cells", cells)):
            if selector is None:
                continue
            values = (
                [int(selector)]
                if isinstance(selector, (int, np.integer))
                else [int(item) for item in selector]
            )
            frame = frame[frame["cell" if column == "cells" else column].isin(values)]
        return frame[["per", "layer", "cell", value]].reset_index(drop=True)

    def summary(self) -> pd.DataFrame:
        """Return a one-row summary of the saved field."""

        value = self.value_name
        frame = self.get()
        data = frame[value].to_numpy(dtype=float)
        finite = data[np.isfinite(data)]
        return pd.DataFrame(
            [
                {
                    "label": self._summary_label,
                    "records": int(len(frame)),
                    "periods": int(frame["per"].nunique()),
                    "layers": int(frame["layer"].nunique()),
                    "cells": int(frame["cell"].nunique()),
                    "min": float(finite.min()) if finite.size else float("nan"),
                    "max": float(finite.max()) if finite.size else float("nan"),
                    "mean": float(finite.mean()) if finite.size else float("nan"),
                }
            ]
        )

    # -- unified grammar: dimension + series hooks --------------------------
    def _spatial_layers(self) -> list[int]:
        """Every grid layer -- the field is saved for all layers."""

        return list(range(int(self.nlay)))

    def _spatial_periods(self) -> list[int]:
        """Stress periods present in the saved output."""

        return sorted({int(key[1]) for key in self.kstpkper})

    def _series_default_agg(self) -> str:
        """Collapse cells within a plotted line by mean (averaging the field, not summing)."""

        return "mean"  # heads/conc/temp average over cells; summing them is meaningless

    def _noun_path(self) -> str:
        """``model.hds`` / ``model.conc`` / ``model.temp``, for messages."""

        return f"model.{self._choro_type}"

    def _sections(
        self,
        model=None,
        *,
        line=None,
        cells: int | list[int] | None = None,
        per: int | None = None,
        layer: int | list[int] = 0,
        **kwargs,
    ):
        """Return this model's :class:`XSection` for the grammar's ``xs`` verbs."""

        del model  # single-model surface; SpatialView validates the selector
        if self.model is None:
            raise ValueError(f"{type(self).__name__}.section() requires a parent model.")
        from myflopy.modflow.utils.datatypes.xsections import XSection

        return {
            self.model.name: XSection(
                model=self.model, line=line, cells=cells, per=per, layer=layer, **kwargs
            )
        }

    def long(self, *, values: str | None = None) -> pd.Series:
        """Return the field as a long series indexed by ``kstpkper/layer/cell``.

        Every SAVED time step, not one value per period (:meth:`get` is the
        period-end table). MODFLOW's no-value sentinels (+1e30 inactive, -1e30
        dry; ``|value| >= 1e29``) are ``NaN``, by the rule :meth:`get`,
        :meth:`array` and the maps follow (ledger 315); they were kept as numbers
        until 2026-10-11. The file's own numbers are :attr:`all_values`.

        Parameters
        ----------
        values : str, optional
            The value-table column to return. The table holds ONE column,
            ``store_column`` -- ``"elev"`` for heads, ``"conc"``, ``"temp"`` -- so
            that is the only name accepted; any other raises ``KeyError``,
            ``"head"`` (the name :meth:`get` uses) included. Default ``None``:
            that column. The series is named after it.
        """

        values = values or self.store_column
        if values not in self.all_values.columns:
            raise KeyError(f"Value column {values!r} was not found.")
        series = pd.to_numeric(self.all_values[values], errors="coerce")
        series = series.mask(_is_sentinel(series))
        series.name = values
        return series

    def wide(
        self,
        *,
        index: list[str] | tuple[str, ...] = ("layer", "cell"),
        values: str | None = None,
        agg: str = "first",
    ) -> pd.DataFrame:
        """Pivot the field to one row per layer/cell and one column per ``kstpkper``.

        Columns are named ``kstpkper_<kstp>_<kper>``, one per SAVED time step,
        after the ``index`` columns. MODFLOW's no-value sentinels (+1e30
        inactive, -1e30 dry; ``|value| >= 1e29``) are ``NaN`` before ``agg``
        runs, as in :meth:`long` (ledger 315), so an ``agg`` combining layers
        skips a dry cell -- it folded its -1e30 into the number until
        2026-10-11. A row and step with no value left is ``NaN``.

        Parameters
        ----------
        index : sequence of str, default ('layer', 'cell')
            The row keys, from the value table's own levels: ``"layer"`` and/or
            ``"cell"`` (zero-based). Naming one alone collapses the other through
            ``agg`` -- ``("cell",)`` is one row per cell. A name the table lacks
            raises ``KeyError`` listing the missing ones, ``"per"`` included (this
            table is keyed by ``kstpkper``); ``"kstpkper"`` and the value column
            are in it but are not row keys -- alone, pandas raises ``ValueError``.
        values : str, optional
            The value column, as in :meth:`long`: only ``store_column``
            (``"elev"`` for heads) exists, and any other name raises ``KeyError``.
            Default ``None``: that column.
        agg : str or callable, default 'first'
            The pandas ``pivot_table`` ``aggfunc`` reducing the values that share
            one row and one time step -- a name (``"mean"``, ``"max"``, ``"sum"``)
            or a callable. Under the default ``index`` every ``(layer, cell,
            kstpkper)`` holds one value, so it acts only when ``index`` drops a
            level.
        """

        values = values or self.store_column
        if values not in self.all_values.columns:
            raise KeyError(f"Value column {values!r} was not found.")
        frame = self.all_values.reset_index()
        missing_index = [column for column in index if column not in frame.columns]
        if missing_index:
            raise KeyError(f"Wide index columns were not found: {missing_index}")
        frame[values] = pd.to_numeric(frame[values], errors="coerce")
        frame[values] = frame[values].mask(_is_sentinel(frame[values]))
        wide = frame.pivot_table(
            index=list(index),
            columns="kstpkper",
            values=values,
            aggfunc=agg,
            dropna=False,
        )
        wide.columns = [f"kstpkper_{kstp}_{kper}" for kstp, kper in wide.columns]
        return wide.reset_index()

    def array(
        self,
        *,
        layer: int = 0,
        per: int | None = -1,
        kstpkper: tuple[int, int] | None = None,
        masked: bool = True,
    ) -> np.ndarray:
        """Return one layer's field at one time as a 1-D ``(ncpl,)`` array.

        The quick spatial accessor behind the choropleth maps: it picks a saved
        time, picks a ``layer``, reshapes to one value per Voronoi cell, and (by
        default) converts MODFLOW's no-value sentinels (+1e30 inactive, -1e30 dry)
        to ``NaN`` so the result drops straight into a plot or a per-cell calculation.

        Parameters
        ----------
        layer : int, default 0
            Zero-based model layer (0 == layer 1). A negative layer counts from the
            bottom (``-1`` is the lowest). One at or past the layer count the
            model's views number -- ONE on a DISU, whatever FloPy infers from
            ``IHC`` -- raises ``ValueError`` naming that count, as below ``-nlay``
            does (it raised numpy's bare ``IndexError``, ledger 496).
        per : int or None, default -1
            Zero-based stress period; its LAST SAVED time step is returned
            (:func:`~myflopy.modflow.mf6.period_end.last_saved_step`, ledger 337),
            not a check that the period's end was saved. ``-1`` (default) or
            ``None`` returns the last saved time overall; any other negative
            value, or a period the file does not hold, raises ``ValueError``
            naming the periods available. A numpy integer is accepted. Ignored
            when ``kstpkper`` is given.
        kstpkper : tuple of (int, int), optional
            Exact zero-based ``(kstp, kper)`` to read, taking precedence over
            ``per``; a time the file did not save raises ``ValueError`` listing
            the saved ones. Default ``None``: ``per`` decides.
        masked : bool, default True
            Replace MODFLOW's no-value sentinels (``|value| >= 1e29``: +1e30
            inactive, -1e30 dry) with ``NaN``. ``False`` returns the file's raw
            numbers.

        Returns
        -------
        numpy.ndarray
            A ``(ncpl,)`` float array, one value per Voronoi cell.
        """

        # Use the normalized kstpkper list the class trusts (the value table reads
        # with these); flopy's raw get_kstpkper() is offset for some files.
        keys = [tuple(int(v) for v in key) for key in self.kstpkper]
        if kstpkper is not None:
            key = tuple(int(v) for v in kstpkper)
            if key not in keys:
                raise ValueError(f"kstpkper {key} is unavailable. Available: {keys}")
        elif per in (None, -1):
            key = keys[-1]
        else:
            key = last_saved_step(keys, per)

        if isinstance(layer, (int, np.integer)) and not isinstance(layer, bool) and layer < 0:
            layer = int(layer) + int(self.nlay)
        (layer,) = self._check_layers(layer)
        ncpl = int(self.vor.ncpl)
        field = _as_layer_cell_heads(self.get_data(kstpkper=key), nlay=self.nlay, ncpl=ncpl)
        values = np.asarray(field[layer], dtype=float)
        if masked:
            values[_is_sentinel(values)] = np.nan
        return values

    def remap(
        self,
        grid,
        *,
        per: int = -1,
        layers: int | list[int] | None = None,
        method: str = "linear",
        uncovered: str = "nearest",
        kstpkper: tuple[int, int] | None = None,
    ) -> np.ndarray:
        """Carry this run's solved field onto ANOTHER grid -- starting values for a model there.

        ``mf.ic(strt=old.hds.remap(new_vor))`` starts a model on ``new_vor`` from the heads
        this run ENDED with (``model.conc.remap`` / ``model.temp.remap`` the same for a
        transport or energy model's starting concentration or temperature). It reads the
        field at the END of stress period ``per``, carries each requested layer onto
        ``grid`` and returns an array with a number in every cell.

        **What the file holds no number for.** An inactive cell (+1e30) is a gap. A DRY head
        (-1e30) takes the head of the first wet layer below it in its column -- the water
        table, normally below the dry cell's bottom, so the new model starts that cell dry
        too (MF6's IC documentation: a head below the cell bottom starts the cell dry);
        without rewetting MF6 can keep a cell dry after the water below has risen past its
        bottom, and that cell then starts wet, where the water now is. A dry concentration
        or temperature is a gap. Under NEWTON MF6 writes a dry cell's head as an ordinary
        number below its bottom, and that number is carried as it is.

        **A step that did not converge** is saved too -- a failed run's last step, at the
        time it would have ended, and under IMS ``CONTINUE`` steps mid-run; ``mfsim.lst``
        names each, and reading one is refused.

        **How a value crosses.** ``"linear"`` (default) interpolates between the old cells'
        centres -- the points MF6 solved the field at -- onto each new cell's centre; where
        a new centre lies outside the old centres' hull, or in a triangle reaching past
        their neighbourhood (across an inactive strip or a notch), it takes the value of the
        old cell covering most of it, moved along that cell's own gradient -- never across
        the gap, and never outside the range the old layer holds, so a concentration cannot
        come back negative.
        ``"area"`` takes the area-weighted mean of the old cells under each new cell: the
        rule ``with_grid(remap=True)`` carries an input ``ic.strt`` by. A new cell with NO old
        cell holding a value under it -- outside the old model, or over cells with none --
        is ``uncovered``.

        The array is EAGER: a simulation holding it builds and runs, but ``project.save()``
        skips it (an array cannot be serialized, ledger 193) and the array names its cells
        by position, so it does not follow a later ``with_grid``. To save one, register it
        as a library package (``project.add_package(key, mf.ic(strt=...))`` +
        ``mf.ref(key)``) -- and register a built grid too (``project.add_grid`` +
        ``grid=mf.grid_ref(key)``), which a saved simulation cannot hold either.

        Parameters
        ----------
        grid : VoronoiGridPlus
            The BUILT grid to carry the field onto, in this model's CRS: what
            ``GridSpec.resolve()`` returns, ``model.vor`` of another model, or
            ``mf.load_mf6_run(ws).vor``. A ``GridSpec``, a ``grid_ref``, a model, a FloPy
            grid or a GeoDataFrame raises ``TypeError`` naming the spelling that gives a
            built grid; a grid in another CRS, or sharing no area with this one, raises
            ``ValueError``.
        per : int, default -1
            The stress period, zero-based, read at its END -- the step whose time into the
            period is the period's TDIS length, by the rule :meth:`get` reads (ledger 336),
            but checked for THIS period only. A negative period counts back from the last
            the simulation DECLARES, so ``-1`` (default) is the last stress period, not
            the last one the file reached. A period whose end the file does not hold is
            refused naming ``SAVE HEAD LAST`` (``CONCENTRATION``, ``TEMPERATURE``) -- and,
            when it is the last period the file reached, a run that stopped inside it. A
            period the file holds nothing for is refused naming the file's last record. A
            reader opened from a path with no model takes the period's last saved step,
            unchecked, and the warning says so. Ignored when ``kstpkper`` is given.
        layers : int or sequence of int, optional
            The OLD layer each output row reads, zero-based (a negative counts up from the
            bottom). Repeats are allowed, so a sequence is the new->old layer map when the
            new model has another layer count: ``[0, 0, 1, 2]`` fills four layers from three.
            An int returns one ``(ncpl,)`` row. Default ``None``: every layer, in order --
            right only when the new model has the same number of layers, which this cannot
            check (``grid`` has no layers).
        method : {"linear", "area"}, default "linear"
            ``"linear"`` interpolates between cell centres; ``"area"`` takes the
            area-weighted mean of the old cells under each new cell. Away from a pumping
            well, linear measured 6.4x closer than the area mean to a finer mesh's own
            solved heads (refining; 2.6x coarsening), and 57x closer on a field with a known
            answer; at the well cell neither wins
            (``docs/reproducers/heads_remap_linear_vs_area.py``). The area mean keeps a
            field's area average over each new cell, which can matter for a concentration
            when coarsening.
        uncovered : {"nearest", "raise"}, default "nearest"
            A new cell no old cell holding a value reaches: ``"nearest"`` takes the value of
            the nearest such old cell in its layer and warns ONCE per call, naming how many
            cells lie outside the old model, how many sat over cells with no value, and the
            farthest fill distance; ``"raise"`` refuses with ``RemapRefused`` (a
            ``ValueError``) naming the counts. The same word and default as
            ``with_grid(remap=True)``.
        kstpkper : tuple of (int, int), optional
            An exact zero-based ``(kstp, kper)`` to read instead, unchecked against the
            period's end -- for a run whose output control did not save it. Two integers
            (a ``bool``, a float or a string is refused with ``TypeError``); a time the
            file did not save raises ``ValueError`` listing the saved ones. Default
            ``None``: ``per`` decides.

        Returns
        -------
        numpy.ndarray
            ``(len(layers), grid.ncpl)`` floats -- every layer when ``layers`` is ``None``
            -- or ``(grid.ncpl,)`` for an int ``layers``, in this field's units. Finite in
            every cell: ``mf.ic(strt=...)`` takes it as it is.

        Raises
        ------
        TypeError
            ``grid`` is not a built grid; ``per``, ``kstpkper`` or a layer is not made of
            integers.
        ValueError
            ``method``/``uncovered`` is not one of its choices (checked first); the model
            is DISU (MF6 writes its nodes as one layer, so a layer has no plan-view cells
            of its own) -- refused through the model, or, for a reader with none, because
            the grid's cells overlap; the period or layer is out of range, or the period's
            end is missing (above); the record is a step the run did not converge on; the
            grids' CRSs differ or they share no area; ``RemapRefused`` under
            ``uncovered="raise"`` or when a requested layer holds no value anywhere.

        Examples
        --------
        >>> strt = old.hds.remap(new_vor)                       # every layer, last period
        >>> strt = old.hds.remap(new_vor, per=3)                # end of period 3
        >>> strt = old.hds.remap(new_vor, layers=[0, 0, 1, 2])  # 4 new layers from 3 old
        >>> ic = mf.ic(strt=strt)
        """

        # Deferred, not at module level: importing any `grid` module runs the package's
        # __init__, which imports the plotting stack and through it this module -- a
        # cycle (measured: ImportError on `_is_sentinel`). The literal defaults above are
        # the engine's DEFAULT_FIELD_METHOD / DEFAULT_UNCOVERED, held equal by a test.
        from myflopy.modflow.mf6.grid.remap import check_field_choices, remap_field

        verb = f"{self._noun_path()}.remap"
        check_field_choices(method, uncovered)
        if self.model is not None and getattr(self.model, "grid_type", None) == "disu":
            raise ValueError(
                f"{verb}: this model is DISU, and MF6 writes a DISU's nodes as one layer, "
                "so its cells have no plan-view layer to carry. Carry a DIS or DISV model's "
                "field."
            )
        raw, _, when = read_saved_field(
            self,
            per=per,
            kstpkper=kstpkper,
            perlen=period_lengths(self.model),
            ncpl=int(self.vor.ncpl),
            label=self._binary_text,
            run_dir=_run_directory(self.output_path),
        )
        field, _ = remap_field(
            raw,
            self.vor,
            grid,
            layers=layers,
            method=method,
            uncovered=uncovered,
            water_table=self._dry_takes_water_table,
            label=self._binary_text,
            verb=verb,
            when=when,
        )
        return field

    def to_xugrid(
        self,
        *,
        layers=None,
        times=None,
        name: str | None = None,
        masked: bool = True,
        surfaces: bool = False,
    ):
        """Export the simulated field across layers and time as an xugrid object.

        Stacks the saved field into a single ``(time, layer, cell)``
        :class:`xugrid.UgridDataArray` on this model's Voronoi mesh -- the
        convenience behind ``model.to_xugrid()``. Unlike :meth:`array` (one
        layer at one time), this returns the whole history at once, ready for
        xarray-style slicing (``uda.isel(time=-1, layer=0)``), native unstructured
        plotting (``.ugrid.plot()``), and UGRID-NetCDF sharing
        (``.ugrid.to_netcdf(...)``). Uses the same normalized ``kstpkper`` and
        dry-cell masking as :meth:`array`, and the grid topology from
        :meth:`~myflopy.modflow.mf6.grid.voronoi.VoronoiGridPlus.ugrid2d`.

        Parameters
        ----------
        layers : sequence of int, optional
            Zero-based layers to include, in the order given. Default ``None``:
            all ``nlay`` layers. Each is a numpy index, so ``-1`` reads the lowest
            layer (and the ``layer`` coordinate keeps ``-1``), and one at or past
            ``nlay`` raises ``IndexError``.
        times : sequence of tuple of (int, int), optional
            Zero-based ``(kstp, kper)`` keys to include (default: every saved
            time step -- not one per period). Order is preserved; an unavailable
            key raises ``ValueError`` naming the saved ones.
        name : str, optional
            Variable name for the field. Default ``None``: this reader's
            ``value_name`` (``"head"``, ``"conc"`` or ``"temp"``).
        masked : bool, default True
            Replace MODFLOW's no-value sentinels (``|value| >= 1e29``: +1e30
            inactive, -1e30 dry) with ``NaN``; ``False`` keeps the raw numbers.
        surfaces : bool, default False
            Also attach the grid's layer frame -- ``top`` and ``botm_1`` ...
            ``botm_n``, one FLAT variable per layer
            (:meth:`~myflopy.modflow.mf6.grid.voronoi.VoronoiGridPlus.layer_surface_fields`)
            -- and return a ``UgridDataset`` rather than a ``UgridDataArray``.
            For QGIS, reduce the field to one value per cell in the same breath;
            see the Notes.

        Returns
        -------
        xugrid.UgridDataArray or xugrid.UgridDataset
            Dims ``("time", "layer", <face_dim>)``, with ``kstp``/``kper`` as
            coordinates on the ``time`` axis and the zero-based layer index on
            ``layer``. With ``surfaces=True`` a ``UgridDataset`` carrying that
            field plus one flat elevation variable per layer.

        Raises
        ------
        ImportError
            If the optional ``xugrid`` / ``xarray`` packages are not installed.
        ValueError
            If a requested ``times`` key is not a saved output time.

        Notes
        -----
        **Analysis vs. plotting shape.** The result keeps its full
        ``(time, layer, cell)`` shape -- that is what you want for differencing,
        aggregation, and NetCDF export. But ``.ugrid.plot()`` draws *one value per
        cell*, so reduce to a single field first with ``.isel(time=-1, layer=0)``
        (or ``.sel(...)``). xugrid does **not** facet over extra dimensions, so
        ``.ugrid.plot(col="layer")`` will not work; plot several layers by looping
        and passing ``ax=`` (see Examples). The same applies to NetCDF if you only
        want one slice -- though the full 3-D array also writes fine.

        **Writing a file QGIS can read.** A mesh dataset there is one value per
        face, so QGIS DROPS a variable that still carries ``time`` or ``layer``,
        and drops it silently -- measured against QGIS 4.2.2, which listed no
        elevation datasets for a layered array and all of them once flattened.
        The elevations ``surfaces=True`` attaches are already flat; the field is
        not, so reduce it on the way out::

            model.to_xugrid(surfaces=True).isel(time=-1, layer=0).ugrid.to_netcdf("model.nc")

        ``isel`` leaves the flat elevation variables untouched -- they carry
        neither dimension -- so that one call flattens the field and keeps the
        whole layer frame.

        Examples
        --------
        >>> uda = model.to_xugrid()                  # (time, layer, cell)
        >>> uda.isel(time=-1, layer=0).ugrid.plot()  # one map: last time, top layer
        >>> uda.isel(time=-1).ugrid.to_netcdf("heads.nc")   # share all layers

        Heads and the layer frame, flattened for QGIS:

        >>> uds = model.to_xugrid(surfaces=True)            # head + top + botm_n
        >>> uds.isel(time=-1, layer=0).ugrid.to_netcdf("model.nc")

        Plot every layer as its own subplot:

        >>> import matplotlib.pyplot as plt
        >>> final = uda.isel(time=-1)                       # (layer, cell)
        >>> fig, axes = plt.subplots(1, final.sizes["layer"])
        >>> for k, ax in zip(final["layer"].values, axes):
        ...     final.sel(layer=k).ugrid.plot(ax=ax)
        """

        xu = require("xugrid", feature="to_xugrid() unstructured-grid export")
        import xarray as xr  # guaranteed present: xarray is a xugrid dependency

        name = name or self.value_name
        ncpl = int(self.vor.ncpl)
        available = [tuple(int(v) for v in key) for key in self.kstpkper]
        if times is None:
            keys = available
        else:
            keys = [tuple(int(v) for v in key) for key in times]
            missing = [key for key in keys if key not in available]
            if missing:
                raise ValueError(
                    f"to_xugrid(): kstpkper {missing} unavailable. "
                    f"Available: {available}"
                )
        lays = list(range(self.nlay)) if layers is None else [int(layer) for layer in layers]

        data = np.empty((len(keys), len(lays), ncpl), dtype=float)
        for ti, key in enumerate(keys):
            field = _as_layer_cell_heads(self.get_data(kstpkper=key), nlay=self.nlay, ncpl=ncpl)
            for li, layer in enumerate(lays):
                data[ti, li, :] = np.asarray(field[layer], dtype=float)
        if masked:
            data[_is_sentinel(data)] = np.nan

        grid = self.vor.ugrid2d()
        face_dim = grid.face_dimension
        array = xr.DataArray(
            data,
            dims=("time", "layer", face_dim),
            coords={
                "time": np.arange(len(keys)),
                "kstp": ("time", [key[0] for key in keys]),
                "kper": ("time", [key[1] for key in keys]),
                "layer": lays,
            },
            name=name,
        )
        if not surfaces:
            return xu.UgridDataArray(array, grid)
        dataset = xr.Dataset({name: array})
        for field, values in self.vor.layer_surface_fields().items():
            dataset[field] = xr.DataArray(
                np.asarray(values, dtype=float), dims=(face_dim,)
            )
        return xu.UgridDataset(dataset, grids=[grid])

    def map(
        self,
        *,
        # -- which numbers ---------------------------------------------------
        per: int | None = None,
        kstpkper: tuple[int, int] | None = None,
        per_timestep: str | int = "last",
        layer: int = 0,
        # -- colour ----------------------------------------------------------
        zmin: float | None = None,
        zmax: float | None = None,
        colorscale: str | list | tuple | None = None,
        logscale: bool = False,
        # -- contours --------------------------------------------------------
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str | None = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        # -- highlighting ----------------------------------------------------
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        # -- overlays --------------------------------------------------------
        features=None,
        hillshade_path=None,
        bgs: bool = False,
        # -- framing ---------------------------------------------------------
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        extent=None,
        # -- hover -----------------------------------------------------------
        hover=None,
        hover_layers: str | None = None,
        hover_surfaces: bool = False,
        show_layer_elevs: bool | None = None,
        show_mounding: bool = False,
        # -- renderer --------------------------------------------------------
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """This field as a choropleth, optionally with contours and overlays.

        Every parameter is named rather than forwarded through ``**kwargs``:
        PyCharm and Pylance read the ``def`` line and never run the module, so a
        parameter that arrives through a tail is one an editor can never offer
        (plan 8.8). ``show_mounding`` is the reason it matters here -- it does not
        decorate the map, it rewrites every value on it.

        ``hover`` accepts a :class:`~myflopy.modflow.utils.datatypes.hover.HoverSpec`
        for full control of the sectioned, styled hover; otherwise this field's
        default spec is used. ``hover_layers`` (``"active"`` | ``"active+strip"`` |
        ``"all"`` | ``"none"``; anything else raises) sets how the vertical
        profile shows, and ``hover_surfaces=True`` adds the model-top/layer-bottom
        column (merged on ``"all"``) -- on a ``hover=`` spec too, as at
        ``model.plot.map``.

        Parameters
        ----------
        backend : {'plotly', 'mpl'}, default 'plotly'
            Renderer. ``'mpl'`` returns a Matplotlib figure instead of a Picture.
        **trace_kwargs
            The open tail, forwarded to ``model.plot.map``. ``values=``,
            ``type=`` and the legacy ``custom_hover``/``hover_heads``/``hover_ks``
            raise ``TypeError`` (see Raises), and ``hover_spec=`` is read as a
            legacy spelling of ``hover`` (``hover`` wins when both are given). A
            key naming a :func:`myflopy.plot.map` parameter this method does not
            name -- ``zoom`` (``hover_fields`` raises ``TypeError``,
            ledger 498), ``animation_kstpkpers``, or
            ``rch_scale``, which acts on a recharge map only -- reaches it as that
            parameter. Anything else rides to the ``go.Choroplethmap`` trace
            (``zmid``, ``colorbar``, ``reversescale``), where Plotly refuses an
            unknown name with ``ValueError`` when the figure is built.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        per : int, optional
            Stress period to read (0-based). Mutually exclusive with ``kstpkper``;
            with neither, the model's first output time is used.
        kstpkper : tuple of (int, int), optional
            Exact ``(timestep, period)`` to read, as MODFLOW reports it. Use
            ``model.kstpkper`` to list what is available.
        per_timestep : {'last', 'first'} or int, default 'last'
            Which timestep WITHIN ``per`` to read, when a period has several. Ignored
            when ``kstpkper`` names the timestep outright.
        layer : int, default 0
            Zero-based layer. Layer 0 is the top. It must be one of the layers the
            model's views number -- the grid's, and ONE on a DISU model, whose nodes
            MF6 writes as one record whatever layering FloPy infers -- or the call
            raises ``ValueError`` naming the count (ledger 448); ``-1`` is
            accepted only beside ``show_mounding=True``. A bare grid checks nothing.
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            Any Plotly colorscale name -- ``'earth'`` (the default), ``'Viridis'``,
            ``'Turbo'``, ``'RdBu'``, ...; append ``_r`` to reverse one -- or explicit
            ``[[position, color], ...]`` stops. The static backend builds its
            colormap from the same stops, so a name renders the same way on both.
            An unknown name warns and falls back to ``'earth'``.
        logscale : bool, default False
            Color on a log10 scale, for a field that spans orders of magnitude
            (hydraulic conductivity, recharge). The colorbar is labelled in real
            units at whole decades on both backends and the hover keeps the real
            values. ``zmin``/``zmax`` (and a ``zmid``) are then in DATA units --
            ``zmin=1e-4, zmax=10`` -- and a non-positive limit raises. Non-positive
            values draw as gaps.
        contours : bool or str, default False
            Overlay contour lines. ``True`` contours the mapped values; a string
            names a different field to contour instead.
        contour_values : sequence of float, optional
            Contour a supplied array rather than the mapped one.
        contour_levels : int or float or list of float, default 10
            A count of levels, a fixed interval, or explicit level values.
        contour_color : str, default 'black'
            Line colour for the contour trace.
        contour_width : float, default 1.5
            Line width for the contour trace.
        contour_name : str, optional
            Legend name for the contour trace.
        contour_clip : bool, default True
            Clip contours to the active domain instead of the full grid extent.
        contour_resolution : int, default 150
            Grid size used to build the contours, per axis. Higher is smoother and
            slower. **Only acts under ``contour_method="cubic"``**, which is the one
            that interpolates onto a square grid before contouring; the default
            linear method triangulates the cell centres directly, so there is no
            grid for this to size and passing it changes nothing (measured: 338
            contour points at 40, 150 and 400 alike).
        contour_method : {'linear', 'cubic'}, default 'linear'
            How the scattered cell values become contours. ``'linear'``
            triangulates the cell centres and contours the triangulation --
            fast, exact at the centres, and faceted. ``'cubic'`` interpolates onto
            a ``contour_resolution``-square grid first (Clough-Tocher) and contours
            that -- smoother, slower, and able to overshoot between cells.
            ``'tri'``/``'tricontour'`` and ``'clough'``/``'clough_tocher'``/
            ``'cloughtocher'`` are accepted as aliases. Anything else raises naming
            both; ``'nearest'`` in particular was documented here for a while and
            has never been implemented.
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight. Cell indices, a boolean mask of length ``ncpl``, the
            name of a registered model region (``"all_streams"``), a path to a
            vector file, or a shapely/GeoPandas geometry, by coverage. Highlights
            nothing (with a warning) when the selection is empty.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            How the highlight is drawn. ``'outline'`` traces the dissolved boundary
            of the selection and leaves every cell at full opacity. ``'dim'`` fades
            the *unselected* cells to 20% instead, which suits a bare grid where the
            selection is the subject, but on a field it costs a measured 6.9x of
            readable contrast everywhere you did not select -- and one box/lasso
            gesture in the browser overwrites it. ``'both'`` draws each.
        select_color : str, optional
            Highlight colour, defaulting to ``viz.PALETTE.highlight``. Worth setting
            when the default red collides with a red-blue diverging colorscale.
        features : path, GeoDataFrame, geometry, Feature, or list/dict of those, optional
            Custom geometry drawn over the cells -- wells, pond footprints, a
            boundary. A path to any vector file GeoPandas reads (name a GeoPackage
            layer with ``Feature(path, layer=...)``), a ``GeoDataFrame`` or
            ``GeoSeries``, a shapely geometry, or a :class:`myflopy.plot.Feature`
            carrying its own ``color``/``width``/``size``/``fill``/``label``/
            ``hover``. A list draws several layers; a dict names them. Points
            become markers, lines polylines, polygons outlines (or fills). Each
            layer gets a legend entry, a stable house color and a hover listing its
            attribute columns. A source without a CRS is taken to be in the grid's.
            Drawn by both backends and kept by ``mosaic`` and ``animate``.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        bgs : bool, default False
            Map depth below ground: each cell's head minus the model top, negative
            below it. Blank where the layer has no head, dry or inactive (ledger 315).
            ``show_mounding`` takes precedence.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``. The
            fit is for the figure's own size, or plotly's default 700x450 when it
            names none -- a larger display only adds margin around the grid -- and
            ``.save()`` to an image refits for the export's size.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        extent : tuple, points, geometry, GeoDataFrame, area or str, optional
            The window the map OPENS on, on both backends -- the one way to frame a
            static (matplotlib) map, which cannot be zoomed. Four numbers
            ``(xmin, ymin, xmax, ymax)`` in the grid's CRS (shapely's order, NOT
            FloPy's ``(xmin, xmax, ymin, ymax)``); a run of ``(x, y)`` points (their
            bounding box); any shapely geometry -- a polygon, a transect line -- by its
            bounds; a GeoDataFrame, reprojected from its own CRS; or the name of an
            area of interest (``model.aoi``). It sets the VIEW only: no cell is
            cropped and the colour scale still spans the whole field, so pass
            ``zmin``/``zmax`` to tighten it for a close-up. ``bounds_padding`` pads it
            as it pads the grid. A box off the grid raises, printing both boxes and
            both CRSs. In a ``mosaic``, one extent frames every panel and a LIST of
            them frames each panel on its own.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        hover_layers : {'active', 'active+strip', 'all', 'none'}, optional
            How the hover reads a field with one value per layer (heads,
            concentration, temperature, transmissivity): ``'active'`` the mapped
            layer's value alone; ``'active+strip'`` adds one line with the other
            layers' values; ``'all'`` a table with one row per layer, the mapped one
            accented; ``'none'`` no layer readout -- the value line goes too, leaving
            the title, ``hover_fields`` and the footer (with ``hover_surfaces``, the
            elevations alone). Default ``None`` keeps the map's own default
            (``'active+strip'`` on those fields). Any other value raises
            ``ValueError`` at the call, naming the four. No effect on a field with
            one value per cell, nor on the legacy ``custom_hover``.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.
        show_mounding : bool, default False
            Map mounding instead of head: each cell's head above its layer's bottom
            (above the model top with ``layer=-1``), clipped at 0, and add it to the
            hover. A DRY cell reads 0 and an inactive one is blank (ledger 315).

        Returns
        -------
        Choro or matplotlib.figure.Figure
            With ``backend='plotly'`` (the default), a
            :class:`~myflopy.viz.Picture`: it renders itself in Jupyter, and answers
            ``.fig``, ``.show()``, ``.save(path)`` and ``.html(path)``. It also
            carries ``.plot_mpl()`` for a static rendering and ``.ani`` for the
            animation over periods.

            With ``backend='mpl'``, a bare Matplotlib ``Figure`` -- not a Picture, so
            use ``.savefig(path)`` and ``.axes[0]`` rather than the picture verbs.

            .. end generated

        Raises
        ------
        TypeError
            If given ``values=`` or ``type=`` -- this noun draws its own field, so
            an override would repaint the cells while the hover, title and
            colorscale went on describing the real one -- or one of the legacy
            hover arguments, which a noun's own ``hover_spec`` supersedes.

        See Also
        --------
        section : the same field as a vertical slice.
        array : the values behind the picture, as a plain array.
        myflopy.plot.map : the free verb, which DOES take ``values=``.

        Examples
        --------
        >>> model.hds.map()                                  # last time, layer 0
        >>> model.hds.map(per=3, layer=1)
        >>> model.hds.map(kstpkper=(0, 5), contours=True, contour_levels=10)
        >>> model.hds.map(zmin=100, zmax=125)                # a scale held steady
        >>> model.hds.map(show_mounding=True)                # head above initial
        >>> model.hds.map(select="all_streams", select_style="both")
        >>> model.hds.map(backend="mpl").savefig("heads.png")
        """

        if self.model is None:
            raise ValueError(f"{type(self).__name__}.map() requires a parent model.")
        refuse_noun_parameters(
            f"model.{self._choro_type}", self._choro_type, trace_kwargs,
            honours=LAYER_FIELD_MAP_PARAMS,
        )
        # Checked at the call, beside `hover=` too (ledger 400).
        layers = (
            "active+strip" if hover_layers is None
            else check_hover_layers(hover_layers, name="hover_layers")
        )
        hover = resolve_noun_hover(hover, trace_kwargs)
        if hover is None:
            hover = self._default_hover(layers=layers, surfaces=hover_surfaces)
        choro = self.model.plot.map(
            per=per,
            kstpkper=kstpkper,
            per_timestep=per_timestep,
            layer=layer,
            type=self._choro_type,
            zmin=zmin,
            zmax=zmax,
            colorscale=colorscale,
            logscale=logscale,
            contours=contours,
            contour_values=contour_values,
            contour_levels=contour_levels,
            contour_color=contour_color,
            contour_width=contour_width,
            contour_name=contour_name,
            contour_clip=contour_clip,
            contour_resolution=contour_resolution,
            contour_method=contour_method,
            select=select,
            select_style=select_style,
            select_color=select_color,
            features=features,
            hillshade_path=hillshade_path,
            bgs=bgs,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
            show_layer_elevs=show_layer_elevs,
            show_mounding=show_mounding,
            hover_spec=hover,
            # Applied on top of a caller's `hover=` as the free verb applies them:
            # the noun dropped both beside one, so `hover_layers="none"` drew the
            # layer readout here and not at `model.plot.map` (ledger 402).
            # `hover_surfaces` is a bool here, so only True can be told from unset.
            hover_layers=hover_layers,
            hover_surfaces=True if hover_surfaces else None,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)

    def _build_value_table(self):
        """Build the full ``(kstpkper, layer, cell)`` value table for this model."""

        vor_cell_list = list(self.vor.gdf_vorPolys.index)

        """set generic MultiIndex for all stress periods and all cells"""
        mdx = pd.MultiIndex.from_product(
            iterables=[
                self.kstpkper,
                list(range(self.nlay)),
                vor_cell_list
            ],
            names=['kstpkper', 'layer', 'cell']
        )
        """Set up a MultiIndex DataFrame to hold the values
            for all cells and stress periods"""
        df_values = pd.DataFrame(
            index=mdx,
            columns=[self.store_column]
        )
        """get data for each stress period"""
        for kstpkper in self.kstpkper:
            sp_values = _as_layer_cell_heads(
                self.get_data(kstpkper=kstpkper),
                nlay=self.nlay,
                ncpl=len(vor_cell_list),
            )
            """copy and paste this stress period data to the MultiIndex DataFrame"""
            for layer in range(self.nlay):
                assert self.vor.ncpl == sp_values.shape[1], (
                    'Are you using the wrong voronoi grid??? \n'
                    f'The provided vor grid has {self.vor.ncpl} cells, but there are {sp_values.shape[1]} values '
                    f'in the model'
                )
                df_values.loc[idxx[kstpkper, layer, :]] = sp_values[layer].reshape(-1, 1)
        return df_values


class HeadsPlus(DependentVariableFile):
    """Extended heads-file reader (GWF) with head-only observation + choropleth extras.

    ``model.hds``. Adds the head observation helpers, the legacy ``choropleth``
    builders, and the head-below-bottom "dry" affordances on top of the shared
    :class:`DependentVariableFile` grammar.
    """

    value_name = "head"
    store_column = "elev"
    _output_suffix = ".hds"
    _binary_text = "head"
    _summary_label = "hds"
    _choro_type = "hds"
    _dry_takes_water_table = True

    def __init__(
            self,
            hds_path: Path = None,
            model=None,
            vor: Vor = None,
            obs_path: Path = None
    ):
        """Parameters
        ----------
        hds_path
            Path to a binary MF6 heads file. If omitted, ``model`` is used.
        model
            Optional parent model used to infer file paths and geometry helpers.
        vor
            Optional Voronoi/grid helper, defaulting to ``model.vor``.
        obs_path
            Optional observation-point file used by observation helper methods.
        """

        super().__init__(path=hds_path, model=model, vor=vor)
        # Historical aliases: ``.hds_path`` and a ``.hds`` self-reference (no
        # second binary reader is opened for the same file).
        self.hds_path = self.output_path
        self.hds = self
        self.obs_heads_df = None
        self._obs = {}
        self._obs_heads = None
        self.obs_path = obs_path

    def _default_hover(self, *, layers: str, surfaces: bool):
        """The default heads hover spec, in the model's own length unit."""

        return head_hover(
            layers=layers, surfaces=surfaces, unit=unit_labels(self.model)[0]
        )

    @property
    def all_heads(self):
        """The full ``kstpkper/layer/cell`` head table (alias of ``all_values``)."""

        return self.all_values

    def get_all_heads(self):
        """Build the full head table (alias of the generic value-table builder)."""

        return self._build_value_table()

    @property
    def obs(self):
        """Mapping of registered head-observation location sets."""

        return self._obs

    @property
    def obs_path(self):
        """Filesystem path to the head observation output file, if configured."""

        return self._obs_path

    @obs_path.setter
    def obs_path(self, val):
        """Set the path to the head observation output file."""

        self._obs_path = val

    @property
    def obs_heads(self):
        """Observed/simulated heads at the observation locations (loaded and cached)."""

        if self._obs_heads is None:
            self._obs_heads = self.get_obs_heads()
        return self._obs_heads

    @staticmethod
    def sort_dict_by_keys(
            dict_to_sort: dict = None
    ):
        """Sorts the given dict by its keys and returns the sorted dict.

        Parameters
        ----------
        dict_to_sort : dict, optional
            The mapping to order. Its keys go through ``sorted()``, so they must
            compare with one another (``str`` with ``int`` raises ``TypeError``).
            Not modified: a NEW ``dict`` is returned, in ascending key order.
            Default ``None`` returns ``{}``.
        """
        return _sort_dict_by_keys(dict_to_sort)

    def get_obs_cells(self, locs: Path, crs: str = None, loc_name_field='ExploName'):
        """Compatibility wrapper for the dedicated observation-cell helper.

        Parameters
        ----------
        locs : pathlib.Path, str or None
            Required -- no default is declared. A vector file GeoPandas reads
            (its FIRST layer) holding the observation points. A point gets the ONE
            cell holding it, by the grid's own selection rule. A feature reaching
            several cells (a line or polygon) RAISES ``ValueError`` naming it --
            it got the first of them in silence until ledger 476. A feature
            reaching no cell is left out, named in a logged warning. ``None``
            reads :attr:`obs_path`, and raises ``ValueError`` when that is unset
            too.
        crs : str, int or pyproj.CRS, optional
            The CRS of a file that carries NONE (ledger 476; it changed no
            cell before). A file carrying its own CRS is reprojected from it to
            the grid's, and a ``crs`` naming another one raises ``ValueError``.
            Default ``None``: a file with no CRS is taken to be in the grid's.
        loc_name_field : str, default 'ExploName'
            The attribute column naming each point; the names are the result's
            keys. A file without that column is keyed by row position
            ``0..n-1`` instead, with no error.

        Returns
        -------
        dict
            ``{name: cell}``, zero-based cell ids within a layer, in ascending
            name order.
        """

        return _get_obs_cells(self, locs=locs, crs=crs, loc_name_field=loc_name_field)

    def get_obs_heads(
            self, locs: Path | list[int] = None,
            crs: str = None,
            loc_name_field='ExploName',
            long_format=False
    ):
        """Compatibility wrapper for the dedicated observation-head helper.

        Every SAVED time step in every layer: no period-end selection
        (:meth:`get` makes one). MODFLOW's no-value sentinels (+1e30 inactive,
        -1e30 dry) read ``NaN``, as in :meth:`long` (ledger 315); they were the
        file's numbers until ledger 495, so a residual against a dry cell
        read 1e30.

        Parameters
        ----------
        locs : str, pathlib.Path or list of int, optional
            The locations. A path (``str`` or ``pathlib.Path``) is stored as
            :attr:`obs_path` and resolved through :meth:`get_obs_cells` NOW, so
            the columns are the point names -- a second path on one reader reads
            that path (a ``str`` was ignored for the constructor's file, and a
            second path raised, until ledger 476). A ``list`` of Python
            ``int`` is read as zero-based cell ids and names the columns itself; a
            numpy integer, or an id that is not a cell, raises ``AssertionError``.
            Default ``None``: the locations this reader resolved last, else
            :attr:`obs_path`, else ``ValueError``. Anything else raises
            ``TypeError``.
        crs : str, int or pyproj.CRS, optional
            Forwarded to :meth:`get_obs_cells`: the CRS of a points file that
            carries none. Unused for a list. Default ``None``.
        loc_name_field : str, default 'ExploName'
            The point-name column, forwarded to :meth:`get_obs_cells`. Unused for
            a list.
        long_format : bool, default False
            ``False``: index ``(layer, kstpkper)``, one column per location in
            sorted order. ``True``: melted to index ``(locs, layer, kstpkper)``
            with one ``elev`` column.

        Returns
        -------
        pandas.DataFrame
            The heads, in the model's length units; ``NaN`` where the cell holds
            no head.
        """

        return _get_obs_heads(
            self,
            locs=locs,
            crs=crs,
            loc_name_field=loc_name_field,
            long_format=long_format,
        )



class ConcResults(DependentVariableFile):
    """GWT concentration output reader (``model.conc``).

    The transport twin of :class:`HeadsPlus`: same grammar (``get/summary/map/
    xs/plot/mosaic/animate``), reading the ``.ucn`` concentration binary. ``unit``
    is model-dependent (mass/volume) -- the ``"mg/L"`` default is a convention,
    thread the real unit through when known.
    """

    value_name = "conc"
    store_column = "conc"
    _output_suffix = ".ucn"
    _binary_text = "concentration"
    _summary_label = "conc"
    _choro_type = "conc"
    _default_unit = "mg/L"

    def __init__(self, path: Path = None, model=None, vor: Vor = None, unit: str | None = None):
        super().__init__(path=path, model=model, vor=vor)
        self.unit = unit or self._default_unit

    def _default_hover(self, *, layers: str, surfaces: bool):
        """The default concentration hover spec (unit threaded from the reader)."""

        return conc_hover(unit=self.unit, layers=layers, surfaces=surfaces)

    @property
    def all_conc(self):
        """The full ``kstpkper/layer/cell`` concentration table (alias of ``all_values``)."""

        return self.all_values


class TempResults(DependentVariableFile):
    """GWE temperature output reader (``model.temp``).

    The energy-transport twin of :class:`HeadsPlus`: same grammar, reading the
    temperature binary. ``unit`` defaults to ``"°C"`` by convention.
    """

    value_name = "temp"
    store_column = "temp"
    _output_suffix = ".ucn"
    _binary_text = "temperature"
    _summary_label = "temp"
    _choro_type = "temp"
    _default_unit = "°C"

    def __init__(self, path: Path = None, model=None, vor: Vor = None, unit: str | None = None):
        super().__init__(path=path, model=model, vor=vor)
        self.unit = unit or self._default_unit

    def _default_hover(self, *, layers: str, surfaces: bool):
        """The default temperature hover spec (unit threaded from the reader)."""

        return temp_hover(unit=self.unit, layers=layers, surfaces=surfaces)

    @property
    def all_temp(self):
        """The full ``kstpkper/layer/cell`` temperature table (alias of ``all_values``)."""

        return self.all_values
