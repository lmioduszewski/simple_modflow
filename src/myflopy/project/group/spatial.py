"""Grouped spatial views: the two bases + `GroupHeads`.

``_GroupSpatialView`` is the *faceting* base (group-vs-single);
``_GroupFieldView`` is the *field* base (heads-vs-conc-vs-temp), which
``GroupHeads`` here and ``GroupConc``/``GroupTemp`` in their own modules
configure with five class attributes each.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from myflopy.modflow.mf6.hydraulics import unit_labels
from myflopy.modflow.mf6.package_explorer import (
    LineProfileView,
    SpatialView,
    build_group_input_compare_map_payload,
    get_default_group_compare_colorscale,
)
from myflopy.modflow.utils.datatypes.hover import compare_hover
from myflopy.project.group._shared import (
    _coerce_kstpkper,
    _ensure_group_map_compatible,
    _reduce_to_period_end,
    _resolve_group_compare_target,
)

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


class _GroupSpatialView(SpatialView):
    """:class:`SpatialView` wired for a :class:`ModelGroup`.

    Faceting defaults to one panel per model (``by="model"``) and the model
    axis is the group's members; the host implements ``_spatial_map(per, layer,
    model)`` to render one model's field as a ``Choro``. ``model=None`` means the
    reference model.
    """

    def _spatial_models(self):
        """The group's member model names -- the model axis for faceting."""

        return list(self.group.models)

    def _spatial_reference_model(self):
        """The reference model whose grid/layers frame the panels."""

        return self.group.models[self.group.reference]

    def _spatial_default_facet(self):
        """Group mosaics default to one panel per model."""

        return "model"

    def _spatial_periods(self):
        """Stress periods present in the data (``[0]`` when there is no period axis)."""

        frame = self.get()
        columns = getattr(frame, "columns", [])
        if "per" in columns and not frame.empty:
            return sorted({int(value) for value in frame["per"].dropna().tolist()})
        return [0]

    def _group_target(self, model):
        """Resolve a model selector to a concrete model name (default reference)."""

        if model is None:
            return self.group.reference
        name = str(model)
        if name not in self.group.models:
            raise KeyError(f"Model {name!r} is not in the group.")
        return name


class _GroupFieldView(LineProfileView, _GroupSpatialView):
    """A grouped DEPENDENT-VARIABLE accessor, for one field kind.

    The group-side counterpart of :class:`~myflopy.modflow.mf6.headsplus.\
DependentVariableFile`, and built for the same reason: heads, concentration and
    temperature differ only in *which reader they read and what the value column
    is called*, so the difference belongs in class attributes rather than in
    three copies of the same hundred lines.

    ``_GroupSpatialView`` above is the *faceting* base (group-vs-single); this is
    the *field* base (heads-vs-conc-vs-temp). Subclasses set the five attributes
    below and inherit ``get``/``compare``/``compare_map`` plus the whole
    ``SpatialView`` verb set unchanged.
    """

    #: attribute on a member model holding this field's reader; kind-GATED on
    #: purpose, so a mixed-kind group fails with the reader's own clear error
    #: rather than silently mapping the wrong physics.
    reader_attribute: str = "hds"
    #: attribute on a member model holding the full value table
    table_attribute: str = "all_heads"
    #: the value column inside that table (``DependentVariableFile.store_column``)
    value_column: str = "elev"
    #: short label for the mapped quantity, used in hovers and axis titles
    value_label: str = "head"
    #: unit shown in hovers; the per-kind conventions the hover factories use.
    #: ``None`` (a head) reads the drawn model's own length unit instead -- it was
    #: ``"ft"`` on a model in metres too (ledger 434).
    value_unit: str | None = None

    def _value_unit(self, model) -> str | None:
        """The hover unit of this field on ``model``: :attr:`value_unit`, else its length unit."""

        return self.value_unit if self.value_unit is not None else unit_labels(model)[0]

    def _noun_path(self) -> str:
        """``group.hds`` / ``group.conc`` / ``group.temp``, for messages."""

        return f"group.{self.reader_attribute}"

    def __init__(self, group: ModelGroup):
        """Bind the grouped field accessor to its :class:`ModelGroup`."""

        self.group = group

    # -- unified grammar hooks ----------------------------------------------
    def _spatial_map(self, *, per: int = 0, layer: int = 0, model=None, **kwargs):
        """Render one member's choropleth of the raw field (not deltas)."""

        target = self.group.models[self._group_target(model)]
        reader = getattr(target, self.reader_attribute)
        return reader.map(per=int(per), layer=int(layer), **kwargs)

    def _spatial_periods(self) -> list[int]:
        """Stress periods present in the reference model's saved output."""

        reference = self.group.models[self.group.reference]
        return sorted({int(key[1]) for key in reference.kstpkper})

    def _spatial_value_label(self) -> str:
        """The mapped quantity's label."""

        return self.value_label

    def _series_table(self) -> pd.DataFrame:
        """The period-end table backing ``plot()`` (one value per model/per/layer/cell)."""

        # one value per (model, per, layer, cell): reduce to period-end saves
        return _reduce_to_period_end(self.get())

    def _series_value_column(self, frame) -> str:
        """The column ``plot()`` draws."""

        return self.value_column

    def _series_default_agg(self) -> str:
        """Collapse cells within a plotted line by mean (averaging the field)."""

        return "mean"

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
        """Return member :class:`XSection` objects for the ``xs`` verbs.

        ``model=None`` overlays every member (labeled by model name);
        ``model="F9b"`` sections just that member.
        """

        from myflopy.modflow.utils.datatypes.xsections import XSection

        return {
            name: XSection(
                model=self.group.models[name],
                section_name=name,
                line=line,
                cells=cells,
                per=per,
                layer=layer,
                **kwargs,
            )
            for name in self._resolve_models(model)
        }

    def get(
        self,
        *,
        per: int | None = None,
        kstpkper: tuple[int, int] | None = None,
        layer: int | list[int] | None = None,
        cells: list[int] | None = None,
    ) -> pd.DataFrame:
        """Return this field, aligned across all models in the group.

        Each member's full value table (``all_heads``/``all_conc``/``all_temp``)
        with a ``model`` column: one row per saved step, layer and cell, the
        values as the file holds them -- a DRY or INACTIVE cell keeps MF6's
        ``-1e30``/``+1e30`` here; nothing masks it (ledger 315's rule is not
        applied at this tier, ledger 439).

        Parameters
        ----------
        per : int or None, optional
            Keep one stress period, zero-based, read at its LAST saved step on
            each member's own clock (``period_end.last_saved_step``, ledger 337),
            so members of different ``nstp`` each give their own period end. A
            member that saved nothing in it raises ``ValueError`` naming the
            periods it has. Ignored when ``kstpkper`` is given. Default ``None``:
            every saved step.
        kstpkper : tuple of (int, int) or None, optional
            Keep one exact saved step, zero-based ``(kstp, kper)``, on every
            member -- which, across members of different ``nstp``, need not be
            the same time in the period; a member that did not save it adds no
            rows. Wins over ``per`` when both are given. Default ``None``.
        layer : int or list of int or None, optional
            Keep these zero-based layers -- an int or a LIST (a tuple raises
            ``TypeError``). Default ``None``: every layer.
        cells : list of int or None, optional
            Keep these zero-based cell numbers within a layer -- any iterable of
            ints (a bare int raises ``TypeError``). Default ``None``: every cell.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``kstpkper``, ``layer``, ``cell`` and the field's value
            column (``elev`` for heads, in model length; ``conc``; ``temp``).
        """

        frames = []
        layer_values = None
        if isinstance(layer, list):
            layer_values = [int(value) for value in layer]
        elif layer is not None:
            layer_values = [int(layer)]

        for model_name, model in self.group.models.items():
            frame = getattr(model, self.table_attribute).reset_index().copy()
            frame["model"] = model_name
            selected = _coerce_kstpkper(model, per=per, kstpkper=kstpkper)
            if selected is not None:
                frame = frame[frame["kstpkper"] == selected]
            if layer_values is not None:
                frame = frame[frame["layer"].isin(layer_values)]
            if cells is not None:
                frame = frame[frame["cell"].isin([int(value) for value in cells])]
            frames.append(frame)

        if not frames:
            return pd.DataFrame(
                columns=["kstpkper", "layer", "cell", self.value_column, "model"]
            )

        combined = pd.concat(frames, ignore_index=True)
        return combined[["model", "kstpkper", "layer", "cell", self.value_column]]

    def compare(
        self,
        *,
        per: int | None = None,
        kstpkper: tuple[int, int] | None = None,
        layer: int | list[int] | None = None,
        cells: list[int] | None = None,
    ) -> pd.DataFrame:
        """Compare this field for all models against the group's reference model.

        The rows :meth:`get` selects are reduced to each member's own LAST saved
        step of every period (``_reduce_to_period_end``) and paired with the
        reference's on ``per``, ``layer`` and ``cell``, so members of different
        ``nstp`` align by stress period. The values are differenced as the file
        holds them: a cell DRY or INACTIVE in one member and not the other
        differences against MF6's ``1e30`` sentinel (a ``diff`` of order
        ``1e30``), and one dry in both reads 0.

        Parameters
        ----------
        per : int or None, optional
            Compare one stress period, zero-based -- each member's own last saved
            step of it (:meth:`get`). Default ``None``: every period.
        kstpkper : tuple of (int, int) or None, optional
            Compare one exact saved step, zero-based ``(kstp, kper)``, on every
            member; it then stands for period ``kper`` whatever its ``kstp``.
            Wins over ``per`` when both are given. Default ``None``.
        layer : int or list of int or None, optional
            Compare these zero-based layers -- an int or a LIST (a tuple raises
            ``TypeError``). Default ``None``: every layer.
        cells : list of int or None, optional
            Compare these zero-based cell numbers within a layer (an iterable of
            ints; a bare int raises ``TypeError``). Default ``None``: every cell.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, ``kstpkper`` (the compared member's
            step), ``per``, ``layer``, ``cell``, the value column,
            ``reference_<value>`` and ``diff`` -- always member minus reference.
            Empty, with those columns, when :meth:`get` is.
        """

        data = self.get(per=per, kstpkper=kstpkper, layer=layer, cells=cells)
        reference_column = f"reference_{self.value_column}"
        columns = [
            "model", "reference_model", "kstpkper", "per",
            "layer", "cell", self.value_column, reference_column, "diff",
        ]
        if data.empty:
            return pd.DataFrame(columns=columns)

        # Align on the stress PERIOD, not the full ``(kstp, kper)`` tuple: models
        # with the same physics but different time discretization save the
        # period-end value at a different ``kstp`` (e.g. kstp=13 vs kstp=9), so an
        # exact tuple match finds nothing. Reduce each model to its period-end
        # value and merge periods.
        data = _reduce_to_period_end(data)

        reference = self.group.reference
        ref = (
            data[data["model"] == reference]
            .rename(columns={self.value_column: reference_column})
            .drop(columns=["model", "kstpkper"])
        )
        comp = data[data["model"] != reference].merge(
            ref,
            on=["per", "layer", "cell"],
            how="inner",
        )
        comp["reference_model"] = reference
        # ALWAYS model - reference, never the other way round: the diff map's
        # diverging scale puts red at negative and blue at positive, so flipping
        # this silently inverts every difference map while the shared colorscale
        # helper still tests green.
        comp["diff"] = (
            comp[self.value_column].astype(float) - comp[reference_column].astype(float)
        )
        return comp[columns]

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        colorscale: str | None = None,
        **kwargs,
    ):
        """Diverging choropleth of this field's differences (model - reference).

        Colors every Voronoi cell by ``diff`` (:meth:`compare`) for one stress
        period and layer, using the same diff-map payload + rendering as the
        package diff maps. A cell dry or inactive in one member only draws a
        difference of order ``1e30`` (either sign), and the default colour range
        then spans plus/minus ``1e30``,
        which draws every real difference white -- pass ``zmin``/``zmax``.

        Internal plumbing: the public route is the one ``diff()`` verb --
        ``group.diff().hds.map(...)`` / ``.conc`` / ``.temp``.

        Parameters
        ----------
        model_name : str or None, optional
            The non-reference member to compare. Default ``None`` takes the
            group's one non-reference member and raises ``ValueError`` when there
            are several; the reference itself raises ``ValueError`` and a name not
            in the group ``KeyError``. The member must have the reference's
            discretization FILE content, read from each workspace (a changed
            ``botm`` or ``idomain`` refuses too), grid type and ``ncpl``, or
            ``ValueError``.
        per : int, default 0
            The stress period drawn, zero-based -- each member's own last saved
            step of it; a member that saved nothing in it raises ``ValueError``.
        layer : int, default 0
            The zero-based layer drawn.
        multiplier : float, default 1.0
            Scales every drawn difference -- a unit conversion or a sign flip --
            and the filled cells with them, except when the selection is empty
            (every cell is then ``fill_value`` as given).
        fill_value : float, default 0.0
            The value of a cell with no comparison row; every cell has one in a
            member's value table, so it acts only when a selection is empty.
            ``float("nan")`` leaves such cells blank.
        colorscale : str or list or None, optional
            Any Plotly colorscale name or explicit stops. Default ``None``: the
            group diff scale -- red at the negative end (lower than the
            reference), white at zero, blue at the positive end.
        **kwargs
            Forwarded to the reference model's ``model.plot.map(...)``: its named
            drawing options (``contours``, ``select``, ``features``, ``extent``,
            ``backend``, ...) and its Plotly trace tail. ``zmin``/``zmax`` default
            to plus/minus the largest absolute value drawn (a filled cell's
            included), ``zmid`` to 0 and ``hover_spec`` to the comparison hover
            (``value_unit`` from the class, else the reference model's own length
            unit for heads -- ledger 434 -- and none when it declares none);
            passing them overrides those. The names this method fills
            itself beyond its own parameters (``type``, ``custom_zs``,
            ``custom_hover``, ``hover_heads``, ``hover_ks``) raise ``TypeError``;
            ``values=`` REPLACES the drawn differences while the hover keeps them.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            The map picture, or a bare Matplotlib figure with ``backend="mpl"``.
        """

        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        selected = self.compare(per=per, layer=layer)
        if not selected.empty:
            selected = selected[selected["model"] == target_name]
        values, hover, absmax = build_group_input_compare_map_payload(
            selected,
            ncpl=reference_model.vor.ncpl,
            value_column=self.value_column,
            diff_column="diff",
            model_name=target_name,
            reference_model=self.group.reference,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg="first",  # one value per cell per (per, layer) -- do not sum
        )
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        reference_column = f"reference_{self.value_column}"
        unit = self._value_unit(reference_model)
        kwargs.setdefault(
            "hover_spec",
            compare_hover(
                self.value_column,
                "diff",
                title=f"Δ {self.value_label} vs reference",
                units={} if unit is None else {
                    self.value_column: unit,
                    reference_column: unit,
                    "diff": unit,
                },
                labels={"diff": f"Δ {self.value_label}"},
            ),
        )
        return reference_model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_group_compare_colorscale(),
            **kwargs,
        )


class GroupHeads(_GroupFieldView):
    """Heads accessor for :class:`ModelGroup` (``group.hds``).

    Mirrors the single-model ``model.hds`` leaf: aligned multi-model tables via
    :meth:`get`/:meth:`compare`, plus the unified grammar -- ``map``/``plot``/
    ``xs`` panels and ``mosaic``/``animate`` composers, faceting over the
    group's members (reference by default).

    Every method lives on :class:`_GroupFieldView`; this class is the heads
    *configuration* of it. Its transport twins are
    :class:`~myflopy.project.group.conc.GroupConc` and
    :class:`~myflopy.project.group.temp.GroupTemp`.
    """

    reader_attribute = "hds"
    table_attribute = "all_heads"
    value_column = "elev"
    value_label = "head"
    value_unit = None  # the model's own length unit
