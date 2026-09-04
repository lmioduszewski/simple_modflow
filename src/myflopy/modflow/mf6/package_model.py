"""Top-level model package explorer wiring."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

import geopandas as gpd
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy._logging import get_logger
from myflopy.modflow.mf6.package_inputs import (
    CellPackageInputsExplorer,
    StaticArrayFieldExplorer,
    UzfInputsNamespace,
)
from myflopy.modflow.mf6.package_plotting import (
    _apply_backend,
    as_mpl_figure,
    normalize_backend,
    refuse_noun_parameters,
    resolve_noun_hover,
)
from myflopy.modflow.mf6.package_registry import (
    get_package_explorer_spec,
)
from myflopy.modflow.mf6.package_results import (
    CellPackageResultsNamespace,
    UzfResultsNamespace,
)
from myflopy.modflow.mf6.package_surface_water import (
    LakPackageExplorer,
    SfrPackageExplorer,
    SurfaceWaterPackageExplorer,
)
from myflopy.viz import mosaic as _mosaic

logger = get_logger(__name__)


class PackageExplorer:
    """Namespace for one package's preferred exploration helpers."""

    def __init__(self, model: SimulationBase, package_name: str):
        """Bind a generic package explorer (``.inputs`` / ``.results``) to ``model``."""

        self.model = model
        self.package_name = str(package_name).lower()

    @property
    def inputs(self) -> CellPackageInputsExplorer:
        """Return normalized input helpers for this package."""

        return CellPackageInputsExplorer(self.model, self.package_name)

    @property
    def results(self) -> CellPackageResultsNamespace:
        """Return normalized result helpers for this package."""

        return CellPackageResultsNamespace(self.model, self.package_name)


class StaticArrayPackageExplorer:
    """Namespace for static array fields in one package."""

    def __init__(
        self,
        model: SimulationBase,
        package_name: str,
        fields: Mapping[str, dict[str, str]],
    ):
        """Bind a static-array explorer to ``model`` for a package's declared array ``fields``."""

        self.model = model
        self.package_name = str(package_name).lower()
        self._fields = dict(fields)

    @property
    def declared_fields(self) -> list[str]:
        """The array fields this package DECLARES, without opening the package.

        :attr:`fields` reports what is actually present, which needs the package
        loaded; this is the cheap answer, for callers that only need to know
        what could be drawn (``model.packages.summary()``).
        """

        return list(self._fields)

    def _available_field_items(self) -> list[tuple[str, dict[str, str]]]:
        """The declared ``(field_name, metadata)`` pairs that actually exist on the package."""

        package = self.model.package(self.package_name)
        available = []
        for field_name, metadata in self._fields.items():
            data = getattr(package, field_name, None)
            if data is None:
                continue
            available.append((field_name, metadata))
        return available

    @property
    def fields(self) -> pd.DataFrame:
        """Return supported static array fields available on this package."""

        return pd.DataFrame(
            [
                {
                    "field": field_name,
                    "label": metadata.get("label"),
                    "colorscale": metadata.get("colorscale", "Viridis"),
                }
                for field_name, metadata in self._available_field_items()
            ]
        ).reindex(columns=["field", "label", "colorscale"])

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported array field."""

        frames = [
            self._field(field_name).summary()
            for field_name in self.fields["field"].tolist()
        ]
        if not frames:
            return pd.DataFrame(
                columns=[
                    "label",
                    "records",
                    "periods",
                    "layers",
                    "cells",
                    "value_columns",
                ]
            )
        summary = pd.concat(frames, ignore_index=True)
        summary["field_name"] = summary["value_columns"].apply(
            lambda values: values[0] if values else None
        )
        return summary.merge(
            self.fields.rename(columns={"field": "field_name"}),
            on="field_name",
            how="left",
        )

    def _field(self, field_name: str) -> StaticArrayFieldExplorer:
        """A field-pinned explorer for one static array (raises if unknown or absent on the package)."""

        metadata = self._fields.get(str(field_name))
        if metadata is None:
            raise AttributeError(
                f"{type(self).__name__!s} has no array field {field_name!r}"
            )
        package = self.model.package(self.package_name)
        if getattr(package, str(field_name), None) is None:
            raise AttributeError(
                f"Package {self.package_name!r} has no available array field {field_name!r}"
            )
        return StaticArrayFieldExplorer(
            self.model,
            self.package_name,
            str(field_name),
            label=metadata.get("label"),
            colorscale=metadata.get("colorscale", "Viridis"),
        )

    def __getattr__(self, field_name: str) -> StaticArrayFieldExplorer:
        """Return a supported static array field explorer."""

        return self._field(field_name)


class UzfPackageExplorer:
    """Top-level UZF package explorer namespace."""

    def __init__(self, model: SimulationBase):
        """Bind the top-level UZF explorer (``.inputs`` / ``.results``) to ``model``."""

        self.model = model

    @property
    def inputs(self) -> UzfInputsNamespace:
        """Return the UZF input exploration namespace."""

        return UzfInputsNamespace(self.model)

    @property
    def results(self) -> UzfResultsNamespace:
        """Return the UZF result exploration namespace."""

        return UzfResultsNamespace(self.model)


class HfbResultsExplorer:
    """Simulated flow across each horizontal flow barrier.

    There is no HFB budget record to read -- MODFLOW 6 writes none, because a
    barrier modifies conductance rather than adding flux. The flow is recovered
    instead from ``FLOW-JA-FACE``, which carries one value per cell-to-cell
    connection: a barrier is a connection, so its flow is simply that entry.

    Mapping connection to position needs the model's own ``IA``/``JA``, read from
    the binary grid file MODFLOW 6 writes beside its output. Reconstructing them
    from the grid would be wrong wherever ``idomain`` removes cells, because
    MODFLOW 6 renumbers.
    """

    def __init__(self, model: SimulationBase):
        """Bind the HFB results explorer to ``model``."""

        self.model = model
        self.package_name = "hfb"

    def _connectivity(self):
        """The model's own ``(ia, ja)``, from the binary grid file."""

        import glob

        from flopy.mf6.utils import MfGrdFile

        workspace = str(self.model.sim.sim_path)
        found = glob.glob(workspace + "/*.grb")
        if not found:
            raise FileNotFoundError(
                "no binary grid file (*.grb) beside the model output, so barrier "
                "connections cannot be located in FLOW-JA-FACE. MODFLOW 6 writes one "
                "next to its other output; re-run the model if it is missing."
            )
        grid_file = MfGrdFile(found[0])
        return np.asarray(grid_file.ia), np.asarray(grid_file.ja)

    @property
    def q(self) -> HfbResultsExplorer:
        """The flow noun. Present so the grammar reads like every other package."""

        return self

    def get(self, *, kstpkper=None, per: int | None = None) -> pd.DataFrame:
        """Return the simulated flow across every barrier, per time step.

        Parameters
        ----------
        kstpkper : tuple, optional
            One ``(kstp, kper)``. All saved steps by default.
        per : int, optional
            Restrict to one stress period.

        Returns
        -------
        pandas.DataFrame
            ``kstp``, ``per``, ``layer``, ``cell1``, ``cell2``, ``hydchr`` and
            ``q_cell1_to_cell2`` -- MODFLOW 6's own sign, in the frame the column
            name states: positive is flow FROM ``cell1`` INTO ``cell2``.
        """

        barriers = HfbPackageExplorer(self.model).get()
        if barriers.empty:
            return pd.DataFrame(
                columns=["kstp", "per", "layer", "cell1", "cell2", "hydchr", "q_cell1_to_cell2"]
            )

        ia, ja = self._connectivity()
        offset = 1 if ja.min() == 1 else 0
        reader = self.model._get_budget_reader()
        steps = [kstpkper] if kstpkper is not None else list(self.model._get_budget_kstpkper())
        ncpl = int(self.model.vor.ncpl)

        rows = []
        for step in steps:
            if per is not None and int(step[1]) != int(per):
                continue
            records = reader.get_data(text="FLOW-JA-FACE", kstpkper=step)
            if not records:
                continue
            flows = np.ravel(records[0])
            for barrier in barriers.itertuples():
                node_a = barrier.layer * ncpl + barrier.cell1
                node_b = barrier.layer * ncpl + barrier.cell2
                lo, hi = int(ia[node_a]) - offset, int(ia[node_a + 1]) - offset
                position = next(
                    (k for k in range(lo, hi) if int(ja[k]) - offset == node_b), None
                )
                rows.append(
                    {
                        "kstp": int(step[0]),
                        "per": int(step[1]),
                        "layer": barrier.layer,
                        "cell1": barrier.cell1,
                        "cell2": barrier.cell2,
                        "hydchr": barrier.hydchr,
                        "q_cell1_to_cell2": float(flows[position]) if position is not None else np.nan,
                    }
                )
        return pd.DataFrame(rows)

    def summary(self, **kwargs) -> pd.DataFrame:
        """One row per time step: how much water crossed the barriers, and the extremes."""

        table = self.get(**kwargs)
        if table.empty:
            return pd.DataFrame(columns=["per", "kstp", "barriers", "q_net", "q_abs_total", "q_max_abs"])
        rows = []
        for (period, step), block in table.groupby(["per", "kstp"]):
            flows = block["q_cell1_to_cell2"]
            rows.append(
                {
                    "per": int(period),
                    "kstp": int(step),
                    "barriers": len(block),
                    "q_net": float(flows.sum()),
                    "q_abs_total": float(flows.abs().sum()),
                    "q_max_abs": float(flows.abs().max()),
                }
            )
        return pd.DataFrame(rows)

    def map(
        self,
        *,
        # -- which records ---------------------------------------------------
        kstpkper=None,
        per: int | None = None,
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
        # -- overlays and framing --------------------------------------------
        locs=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        # -- hover -----------------------------------------------------------
        hover=None,
        # -- renderer --------------------------------------------------------
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """Draw the barriers coloured by the flow crossing them.

        Each barrier is drawn on the face it occupies, thickness fixed and colour
        carrying |q| -- a barrier is a line, so it is drawn as one. The cells
        beneath are the plain grid; the drawing parameters below style THAT
        backdrop, which is what makes a barrier legible against it.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        HFB is edge-indexed and MODFLOW 6 writes no HFB budget record, so this
        noun reads ``FLOW-JA-FACE`` through the model's own ``IA``/``JA`` (ledger
        162). There is no ``values=``: the flow across the barrier IS the subject,
        and an override would repaint it while the hover went on reporting q.

        Parameters
        ----------
        kstpkper : tuple of (int, int), optional
            Exact ``(timestep, period)`` to read, as MODFLOW reports it. With
            neither this nor ``per``, the last saved time is drawn -- the LAST
            ``(per, kstp)`` pair, not the largest ``kstp``, because a model whose
            periods each hold one time step has ``kstp == 0`` throughout.
        per : int, optional
            Stress period to read, zero-based. Mutually exclusive with
            ``kstpkper``.
        layer : int, default 0
            Zero-based layer. A barrier sits in one layer, so this selects which
            set of faces is drawn; barriers in other layers are not shown.
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            A Plotly colorscale name, or explicit stops. **Pass stops for a diverging
            scale** -- names round-trip through a plotly-to-matplotlib table that
            maps ``'rdbu'`` to the REVERSED colormap, so a named diverging scale
            renders mirrored between backends (ledger 69/70).
        logscale : bool, default False
            Color on a log scale. Non-positive values are masked.
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
            vector file, or a shapely/GeoPandas geometry to intersect. Highlights
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
        locs : Path or GeoDataFrame, optional
            Point locations to mark -- wells, observations, samples. A path is read
            as a vector file.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map. Two differences are worth knowing before you rely on
            one. The Matplotlib branch draws in **model coordinates** with no
            basemap, so ``bgs`` and ``zoom`` have nothing to act on there; and a
            NAMED diverging colorscale renders mirrored between the two, because the
            name round-trips through a plotly-to-matplotlib table that maps
            ``'rdbu'`` to the reversed colormap -- pass explicit stops when the
            direction carries meaning (ledger 69/70).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            A :class:`~myflopy.viz.Picture` under the default backend, or a bare
            Matplotlib figure with ``backend="mpl"``. The Matplotlib branch draws
            the barrier lines itself rather than deferring to ``Choro.plot_mpl``,
            which reads cell values only and would drop them.

        See Also
        --------
        get : the per-barrier flows behind the picture, as a DataFrame.
        myflopy.plot.map : the same grid with no barriers on it.

        Examples
        --------
        >>> model.packages.hfb.results.q.map()
        >>> model.packages.hfb.results.q.map(per=5, layer=1)
        >>> model.packages.hfb.results.q.map(kstpkper=(0, 5))
        >>> model.packages.hfb.results.q.map(contours=True, contour_levels=6)
        >>> model.packages.hfb.results.q.map(select="all_streams")
        >>> model.packages.hfb.results.q.map(backend="mpl").savefig("hfb_q.png")
        """

        import plotly.graph_objects as go

        refuse_noun_parameters("hfb.results.q", "barrier flow", trace_kwargs)
        hover = resolve_noun_hover(hover, trace_kwargs)
        grid = self.model.vor
        table = self.get(kstpkper=kstpkper, per=per)
        table = table[table["layer"] == int(layer)]
        picture = grid.plot.map(
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
            locs=locs,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            hover=hover,
            **trace_kwargs,
        )
        if table.empty:
            logger.warning("no barrier flow to draw for layer %d", layer)
            return _apply_backend(picture, backend)

        # The LAST (per, kstp) pair, not the largest kstp: a model whose periods
        # each hold one time step has kstp == 0 throughout, so filtering on kstp
        # alone keeps every period and draws each barrier once per period.
        last = table.sort_values(["per", "kstp"]).iloc[-1]
        latest = table[(table["per"] == last["per"]) & (table["kstp"] == last["kstp"])]
        if normalize_backend(backend) == "mpl":
            faces = [
                face for face in
                (grid.shared_face(r.cell1, r.cell2) for r in latest.itertuples())
                if face is not None
            ]
            return _barriers_on_mpl(picture, faces, color="black", width=4,
                                    name="hfb flow")
        for record in latest.itertuples():
            segment = grid.shared_face(record.cell1, record.cell2)
            if segment is None:
                continue
            frame = gpd.GeoSeries([segment], crs=getattr(grid.gdf_vorPolys, "crs", None))
            if frame.crs is not None and str(frame.crs).upper() != "EPSG:4326":
                frame = frame.to_crs("EPSG:4326")
            x, y = frame.iloc[0].xy
            picture.add_overlay(
                go.Scattermap(
                    mode="lines",
                    lon=list(x),
                    lat=list(y),
                    line={"width": 4},
                    name="hfb flow",
                    showlegend=False,
                    hovertemplate=(
                        f"cell {record.cell1} &#8594; {record.cell2}<br>"
                        f"q {record.q_cell1_to_cell2:.4g}<extra></extra>"
                    ),
                )
            )
        return picture


def _barriers_on_mpl(picture, geometries, *, color, width, name=None):
    """Render an HFB map on Matplotlib with the barrier lines actually on it.

    Not a plain :func:`_apply_backend` call, and the difference matters:
    ``Choro.plot_mpl`` reads the cell values and NOTHING else -- overlays are a
    Plotly-side concept and it has never drawn them. Handing an HFB map to it
    unchanged returns a picture of the cells with the barriers silently missing,
    which is the one thing the picture is of.

    Both backends end up in the same frame, which is what makes this a second
    renderer rather than a second drawing: ``plot_mpl`` draws in MODEL
    coordinates (measured -- x limits -105..2205 against grid bounds 0..2100 on
    EPSG:2927), so the segments go on unprojected, where the Plotly overlay needs
    them reprojected to EPSG:4326 like every other ``Scattermap`` trace.
    """

    figure = as_mpl_figure(picture.plot_mpl())
    axes = figure.axes[0]
    for geometry in geometries:
        xs, ys = geometry.xy
        axes.plot(list(xs), list(ys), color=color, linewidth=width, zorder=5,
                  label=name)
    return figure


class HfbPackageExplorer:
    """Exploration namespace for the horizontal-flow-barrier package.

    Bespoke rather than registry-backed, because HFB is not shaped like a
    boundary condition. Its records are cell *pairs* and its geometry is the
    shared EDGE between them, so the generic cell-keyed explorer -- which knows
    exactly one cellid per record -- cannot read it.

    It also has **no results tier, and cannot have one**: MODFLOW 6 writes no HFB
    budget record at all. A barrier reduces the aquifer's own intercell
    conductance rather than adding flux, so its only trace in the output is a
    correction to ``FLOW-JA-FACE``. Asking for ``.results`` says so rather than
    returning an empty frame.
    """

    def __init__(self, model: SimulationBase):
        """Bind the HFB explorer to ``model``."""

        self.model = model
        self.package_name = "hfb"

    @property
    def results(self) -> HfbResultsExplorer:
        """Simulated flow across each barrier.

        MODFLOW 6 writes no HFB *budget record* -- a barrier reduces the aquifer's
        own intercell conductance rather than adding flux. But the flow across a
        barrier is not lost: the barrier sits ON a cell-to-cell connection, and
        ``FLOW-JA-FACE`` carries the flow on every connection. Looking each
        barrier's pair up in that array is what this tier does.
        """

        return HfbResultsExplorer(self.model)

    @property
    def inputs(self) -> HfbPackageExplorer:
        """The barriers. Present so ``.inputs.<verb>()`` reads the same as elsewhere."""

        return self

    def get(self) -> pd.DataFrame:
        """Return one row per barrier: period, layer, the two cells, and ``hydchr``.

        Returns
        -------
        pandas.DataFrame
            ``per``, ``layer``, ``cell1``, ``cell2``, ``hydchr`` -- plus
            ``length``, the shared face's own length, which is what MODFLOW 6
            multiplies ``hydchr`` by to get the barrier's conductance.
        """

        package = self.model.package(self.package_name)
        data = package.stress_period_data.get_data()
        grid = getattr(self.model, "vor", None)

        rows = []
        for period, records in sorted(data.items()):
            if records is None:
                continue
            for record in records:
                first, second, hydchr = record[0], record[1], record[2]
                layer_a, cell_a = int(first[0]), int(first[-1])
                layer_b, cell_b = int(second[0]), int(second[-1])
                length = None
                if grid is not None and layer_a == layer_b:
                    segment = grid.shared_face(cell_a, cell_b)
                    length = None if segment is None else float(segment.length)
                rows.append(
                    {
                        "per": int(period),
                        "layer": layer_a,
                        "cell1": cell_a,
                        "cell2": cell_b,
                        "hydchr": float(hydchr),
                        "length": length,
                        "vertical": layer_a != layer_b,
                    }
                )
        return pd.DataFrame(
            rows,
            columns=["per", "layer", "cell1", "cell2", "hydchr", "length", "vertical"],
        )

    def summary(self) -> pd.DataFrame:
        """Return one compact row per period: how many barriers, where, how tight."""

        table = self.get()
        if table.empty:
            return pd.DataFrame(
                columns=["per", "barriers", "layers", "hydchr_min", "hydchr_max", "total_length"]
            )
        rows = []
        for period, block in table.groupby("per"):
            rows.append(
                {
                    "per": int(period),
                    "barriers": len(block),
                    "layers": ", ".join(str(int(v) + 1) for v in sorted(block["layer"].unique())),
                    "hydchr_min": float(block["hydchr"].min()),
                    "hydchr_max": float(block["hydchr"].max()),
                    "total_length": float(block["length"].sum(skipna=True)),
                }
            )
        return pd.DataFrame(rows)

    def segments(self, *, per: int | None = None):
        """Return the barriers as a GeoDataFrame of the faces they sit on.

        The shape worth keeping: a barrier's geometry is the shared edge, and an
        edge survives a change of grid in a way a cell-pair index does not.
        """

        import geopandas as gpd

        grid = self.model.vor
        table = self.get()
        if per is not None:
            table = table[table["per"] == int(per)]
        rows, geoms = [], []
        for record in table.itertuples():
            if record.vertical:
                continue
            segment = grid.shared_face(record.cell1, record.cell2)
            if segment is None:
                continue
            rows.append(
                {
                    "per": record.per,
                    "layer": record.layer,
                    "cell1": record.cell1,
                    "cell2": record.cell2,
                    "hydchr": record.hydchr,
                }
            )
            geoms.append(segment)
        return gpd.GeoDataFrame(rows, geometry=geoms, crs=getattr(grid.gdf_vorPolys, "crs", None))

    def map(
        self,
        values=None,
        *,
        # -- which barriers --------------------------------------------------
        per: int | None = None,
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
        # -- overlays and framing --------------------------------------------
        locs=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        # -- hover -----------------------------------------------------------
        hover=None,
        # -- renderer --------------------------------------------------------
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """Draw the barriers on the grid, as lines over whatever the cells show.

        A barrier is an edge, so it is drawn as an edge -- over the cell field
        rather than instead of it. ``values`` colours the cells as any other map
        would; without one the grid is drawn plain beneath the barriers.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        ``values`` survives here where the other nouns REFUSE it, and the
        difference is real rather than an oversight: those nouns fix a field, so
        an override would repaint the cells while the hover and colorbar went on
        describing the real one. This noun's subject is the BARRIERS. The cells
        are a backdrop, and choosing what the backdrop shows contradicts nothing
        -- exactly as on ``vor.plot.map``, which this delegates to.

        Parameters
        ----------
        values : array-like, optional
            Per-cell values to colour the cells by -- heads, K, a zone id, any
            array of length ``vor.ncpl``. With none, the grid is drawn plain and
            the barriers are the only thing carrying colour.
        per : int, optional
            Which period's barriers to draw, zero-based. The first defined period
            by default. Barriers are usually time-invariant, so this matters only
            for a model that redefines them.
        layer : int, default 0
            Which layer's barriers to draw, zero-based. A barrier sits on the
            faces of one layer; those in other layers are not shown.
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            A Plotly colorscale name, or explicit stops. **Pass stops for a diverging
            scale** -- names round-trip through a plotly-to-matplotlib table that
            maps ``'rdbu'`` to the REVERSED colormap, so a named diverging scale
            renders mirrored between backends (ledger 69/70).
        logscale : bool, default False
            Color on a log scale. Non-positive values are masked.
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
            vector file, or a shapely/GeoPandas geometry to intersect. Highlights
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
        locs : Path or GeoDataFrame, optional
            Point locations to mark -- wells, observations, samples. A path is read
            as a vector file.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map. Two differences are worth knowing before you rely on
            one. The Matplotlib branch draws in **model coordinates** with no
            basemap, so ``bgs`` and ``zoom`` have nothing to act on there; and a
            NAMED diverging colorscale renders mirrored between the two, because the
            name round-trips through a plotly-to-matplotlib table that maps
            ``'rdbu'`` to the reversed colormap -- pass explicit stops when the
            direction carries meaning (ledger 69/70).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            A :class:`~myflopy.viz.Picture` under the default backend, or a bare
            Matplotlib figure with ``backend="mpl"``. The Matplotlib branch draws
            the barrier lines itself rather than deferring to ``Choro.plot_mpl``,
            which reads cell values only and would drop them.

        See Also
        --------
        segments : the barrier geometries behind the picture, as a GeoDataFrame.
        myflopy.plot.map : the same grid with no barriers on it.

        Examples
        --------
        >>> model.packages.hfb.inputs.map()
        >>> model.packages.hfb.inputs.map(values=model.hds.array(layer=0))
        >>> model.packages.hfb.inputs.map(layer=1, per=0)
        >>> model.packages.hfb.inputs.map(values=k, logscale=True)
        >>> model.packages.hfb.inputs.map(select="all_streams", select_style="both")
        >>> model.packages.hfb.inputs.map(backend="mpl").savefig("barriers.png")
        """

        import plotly.graph_objects as go

        from myflopy.viz import HIGHLIGHT_WIDTH, PALETTE

        refuse_noun_parameters("hfb.inputs", "the barriers", trace_kwargs)
        hover = resolve_noun_hover(hover, trace_kwargs)
        grid = self.model.vor
        picture = grid.plot.map(
            values,
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
            locs=locs,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            hover=hover,
            **trace_kwargs,
        )

        frame = self.segments(per=per)
        if len(frame):
            frame = frame[frame["layer"] == int(layer)]
        if not len(frame):
            return _apply_backend(picture, backend)
        if normalize_backend(backend) == "mpl":
            # `PALETTE.mpl_highlight`, not `.highlight`: the palette carries a
            # matplotlib twin of each colour precisely because Plotly's
            # `rgb(214,39,40)` spelling is not one matplotlib accepts, and
            # passing it raises rather than drawing the wrong colour.
            return _barriers_on_mpl(
                picture, list(frame.geometry),
                color=PALETTE.mpl_highlight, width=HIGHLIGHT_WIDTH,
                name=f"hfb (layer {int(layer) + 1})",
            )
        if frame.crs is not None and str(frame.crs).upper() != "EPSG:4326":
            frame = frame.to_crs("EPSG:4326")
        lon, lat = [], []
        for geom in frame.geometry:
            x, y = geom.xy
            lon.extend([*x, None])
            lat.extend([*y, None])
        picture.add_overlay(
            go.Scattermap(
                mode="lines",
                lon=lon,
                lat=lat,
                line={"color": PALETTE.highlight, "width": HIGHLIGHT_WIDTH},
                name=f"hfb (layer {int(layer) + 1})",
                hoverinfo="skip",
                showlegend=True,
            )
        )
        return picture

    def __repr__(self) -> str:
        """Short identity: how many barriers, in which layers."""

        table = self.get()
        layers = sorted(table["layer"].unique()) if len(table) else []
        return (
            f"HfbPackageExplorer({len(table)} barrier(s), "
            f"layer(s) {', '.join(str(int(v) + 1) for v in layers) or 'none'})"
        )


class ModelPackages:
    """Preferred package exploration namespace for one model/run.

    Examples
    --------
    ``model.packages.rch.inputs.get()``
        Normalized recharge input table.
    ``model.packages.rch.inputs.map(per=0)``
        Recharge choropleth using the shared choropleth styling.
    ``model.packages.uzf.inputs.finf.summary()``
        Compact summary of UZF infiltration inputs.
    """

    def __init__(self, model: SimulationBase):
        """Bind the preferred ``model.packages`` exploration namespace to ``model``."""

        self.model = model

    #: Packages whose fields are static ARRAYS rather than stress-period
    #: records. They answer `model.packages.npf.k.map()` -- one level shallower
    #: than the record packages' `.inputs.<field>.map()` -- so `summary` has to
    #: ask them for their fields differently.
    _STATIC_ARRAY_PACKAGES = ("npf", "ic", "sto")

    def _explorer_or_none(self, package_name: str):
        """The explorer for ``package_name``, or ``None`` if it has none."""

        try:
            return getattr(self, package_name)
        except AttributeError:
            return None

    def summary(self, *, detail: str = "fields") -> pd.DataFrame:
        """One row per package: what it is, and what can be drawn from it.

        The map from "what is in this model" to "what can I look at", which
        otherwise means knowing that record packages answer
        ``.inputs.<field>.map()`` while array packages answer ``.<field>.map()``,
        and that several packages have no explorer at all.

        Parameters
        ----------
        detail : {'fields', 'data'}, default 'fields'
            ``'fields'`` is CHEAP: the registry and the package list only, no
            package data read, so it stays usable on a lazily loaded run.
            ``'data'`` additionally reports ``records``, ``periods`` and
            ``layers`` per package, which requires reading each package's
            records -- and is what tells you that ``wel`` has nothing in layer 0
            before you draw an empty map of it.

        Returns
        -------
        pandas.DataFrame
            Columns ``package``, ``kind``, ``mappable``, ``fields``, ``results``
            (plus ``records``, ``periods``, ``layers`` when ``detail='data'``).

        Raises
        ------
        ValueError
            If ``detail`` is neither value.
        """

        if detail not in ("fields", "data"):
            raise ValueError(
                f"detail must be 'fields' or 'data', not {detail!r}."
            )

        rows = []
        for raw in self.model.package_names:
            name = str(raw).lower()
            spec = get_package_explorer_spec(name)
            explorer = self._explorer_or_none(name)

            if spec is not None:
                kind = spec.kind
                fields = list(spec.inputs)
                results = list(spec.results)
            elif name in self._STATIC_ARRAY_PACKAGES and explorer is not None:
                kind = "static_array"
                fields = explorer.declared_fields
                results = []
            else:
                kind = ""
                fields = []
                results = []

            rows.append({
                "package": str(raw).upper(),
                "kind": kind,
                "mappable": self._is_mappable(explorer, kind, fields),
                "fields": ", ".join(fields),
                "results": ", ".join(results),
            })

        frame = pd.DataFrame(rows).reindex(
            columns=["package", "kind", "mappable", "fields", "results"]
        )
        if detail == "fields":
            return frame
        counts = pd.DataFrame([self._package_data_counts(row) for row in rows])
        # Nullable Int64: a package with no readable records must stay blank
        # rather than turning the whole column into floats and printing "504.0".
        counts["records"] = counts["records"].astype("Int64")
        return frame.join(counts)

    @staticmethod
    def _is_mappable(explorer, kind: str, fields: list) -> bool:
        """Whether this package can actually be drawn.

        PROBED, not inferred from the registry. HFB is deliberately absent from
        `package_registry` -- it is face-indexed, so it has no cellid, and MF6
        writes it no budget record (ledger 162) -- yet
        `model.packages.hfb.inputs.map()` exists and works. Keying off registry
        fields alone would report the one package whose entry is entirely
        hand-written as undrawable, which is exactly backwards.
        """

        if explorer is None:
            return False
        if callable(getattr(getattr(explorer, "inputs", None), "map", None)):
            return True
        if kind == "static_array" and fields:
            return True
        return callable(getattr(getattr(explorer, "results", None), "map", None))

    def _package_data_counts(self, row: dict) -> dict:
        """``records``/``periods``/``layers`` for one package row (reads data)."""

        blank = {"records": None, "periods": "", "layers": ""}
        if not row["mappable"] or row["kind"] == "static_array":
            return blank
        explorer = self._explorer_or_none(str(row["package"]).lower())
        inputs = getattr(explorer, "inputs", None)
        if inputs is None:
            return blank
        # UZF (and any other per-field namespace) has no whole-package `get()`;
        # its records hang off each field. Counting the first declared field is
        # right for a summary: every UZF field carries one row per cell-period.
        source = inputs
        if not hasattr(source, "get"):
            first = (row["fields"].split(", ") or [""])[0]
            source = getattr(inputs, first, None)
            if source is None or not hasattr(source, "get"):
                return blank
        try:
            frame = source.get()
        except (KeyError, ValueError, AttributeError, OSError):
            # A package whose records cannot be read is reported as blank rather
            # than failing the whole table: `summary` is the thing you reach for
            # WHEN something is off, so it must survive one bad package.
            logger.debug("summary: %r records unreadable", row["package"])
            return blank

        def _uniq(column):
            if column not in frame.columns:
                return ""
            return ", ".join(str(v) for v in sorted(frame[column].dropna().unique()))

        return {
            "records": int(len(frame)),
            "periods": _uniq("per"),
            "layers": _uniq("layer"),
        }

    def mosaic(self, *, packages=None, ncols: int = 3, title: str | None = None,
               sync_views: bool = True, **kwargs):
        """Every mappable package, drawn as one grid of small multiples.

        The combinator over :meth:`summary` -- it draws what that table says is
        mappable. Record packages contribute their default input field; array
        packages contribute their first declared array.

        Parameters
        ----------
        packages : list of str, optional
            Restrict to these packages. Default: everything ``summary()`` marks
            mappable.
        ncols : int, default 3
            Grid width.
        title : str, optional
            Overall figure title.
        sync_views : bool, default True
            Pan and zoom the panels together.
        **kwargs
            Forwarded to each package's ``map()`` (e.g. ``per=``, ``layer=``).

        Returns
        -------
        viz.Fig
            One figure. A package that cannot be drawn is skipped with a DEBUG
            log rather than failing the grid.
        """

        table = self.summary()
        wanted = (
            [str(p).lower() for p in packages]
            if packages is not None
            else [str(p).lower() for p in table.loc[table["mappable"], "package"]]
        )

        panels = []
        for name in wanted:
            explorer = self._explorer_or_none(name)
            if explorer is None:
                logger.debug("mosaic: no explorer for %r, skipping", name)
                continue
            drawer = getattr(getattr(explorer, "inputs", None), "map", None)
            if drawer is None and name in self._STATIC_ARRAY_PACKAGES:
                declared = explorer.declared_fields
                drawer = getattr(getattr(explorer, declared[0], None), "map", None) \
                    if declared else None
            if drawer is None:
                logger.debug("mosaic: %r has no map(), skipping", name)
                continue
            try:
                panels.append((name, drawer(**kwargs)))
            except (KeyError, ValueError, AttributeError, OSError) as error:
                # One unmappable package must not lose the other nineteen.
                logger.debug("mosaic: %r could not be drawn (%s)", name, error)

        if not panels:
            raise ValueError(
                "no package could be drawn; `summary()` shows what is mappable."
            )
        return _mosaic(panels, ncols=ncols, title=title, sync_views=sync_views)

    def __getattr__(self, package_name: str) -> PackageExplorer:
        """Return a registry-backed generic package explorer."""

        spec = get_package_explorer_spec(package_name)
        if spec is None or spec.kind != "cell_stress":
            raise AttributeError(
                f"{type(self).__name__!s} has no package {package_name!r}"
            )
        return PackageExplorer(self.model, spec.name)

    @property
    def rch(self) -> PackageExplorer:
        """Recharge package exploration helpers."""

        return PackageExplorer(self.model, "rch")

    @property
    def chd(self) -> PackageExplorer:
        """Constant-head package exploration helpers."""

        return PackageExplorer(self.model, "chd")

    @property
    def drn(self) -> PackageExplorer:
        """Drain package exploration helpers."""

        return PackageExplorer(self.model, "drn")

    @property
    def ghb(self) -> PackageExplorer:
        """General-head boundary package exploration helpers."""

        return PackageExplorer(self.model, "ghb")

    @property
    def riv(self) -> PackageExplorer:
        """River package exploration helpers."""

        return PackageExplorer(self.model, "riv")

    @property
    def wel(self) -> PackageExplorer:
        """Well package exploration helpers."""

        return PackageExplorer(self.model, "wel")

    @property
    def evt(self) -> PackageExplorer:
        """Evapotranspiration package exploration helpers."""

        return PackageExplorer(self.model, "evt")

    @property
    def uzf(self) -> UzfPackageExplorer:
        """UZF package exploration helpers."""

        return UzfPackageExplorer(self.model)

    @property
    def hfb(self) -> HfbPackageExplorer:
        """Horizontal flow barrier exploration helpers.

        Not registry-backed: HFB records are cell PAIRS and MODFLOW 6 writes no
        HFB budget record, so neither the generic cell-keyed explorer nor a
        results tier applies.
        """

        return HfbPackageExplorer(self.model)

    @property
    def ic(self) -> StaticArrayPackageExplorer:
        """Initial conditions array exploration helpers."""

        return StaticArrayPackageExplorer(
            self.model,
            "ic",
            {
                "strt": {"label": "Starting head", "colorscale": "Viridis"},
            },
        )

    @property
    def npf(self) -> StaticArrayPackageExplorer:
        """Node property flow array exploration helpers."""

        return StaticArrayPackageExplorer(
            self.model,
            "npf",
            {
                "k": {
                    "label": "Horizontal hydraulic conductivity",
                    "colorscale": "Viridis",
                },
                "k22": {
                    "label": "Horizontal hydraulic conductivity K22",
                    "colorscale": "Viridis",
                },
                "k33": {
                    "label": "Vertical hydraulic conductivity",
                    "colorscale": "Viridis",
                },
            },
        )

    @property
    def sto(self) -> StaticArrayPackageExplorer:
        """Storage package array exploration helpers."""

        return StaticArrayPackageExplorer(
            self.model,
            "sto",
            {
                "ss": {"label": "Specific storage", "colorscale": "Viridis"},
                "sy": {"label": "Specific yield", "colorscale": "Viridis"},
            },
        )

    @property
    def lak(self) -> LakPackageExplorer:
        """LAK package exploration helpers."""

        return LakPackageExplorer(self.model)

    @property
    def sfr(self) -> SfrPackageExplorer:
        """SFR package exploration helpers."""

        return SfrPackageExplorer(self.model)

    @property
    def surface_water(self) -> SurfaceWaterPackageExplorer:
        """Combined SFR/LAK exploration helpers."""

        return SurfaceWaterPackageExplorer(self.model)


__all__ = [
    "PackageExplorer",
    "HfbPackageExplorer",
    "HfbResultsExplorer",
    "StaticArrayPackageExplorer",
    "UzfPackageExplorer",
    "ModelPackages",
]
