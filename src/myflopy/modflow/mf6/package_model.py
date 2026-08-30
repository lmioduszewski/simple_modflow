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

    def map(self, *, kstpkper=None, per: int | None = None, layer: int = 0, **kwargs):
        """Draw the barriers coloured by the flow crossing them.

        Each barrier is drawn on the face it occupies, thickness fixed and colour
        carrying |q| -- a barrier is a line, so it is drawn as one.

        Returns
        -------
        Choro
        """

        import plotly.graph_objects as go

        grid = self.model.vor
        table = self.get(kstpkper=kstpkper, per=per)
        table = table[table["layer"] == int(layer)]
        picture = grid.plot.map(**kwargs)
        if table.empty:
            logger.warning("no barrier flow to draw for layer %d", layer)
            return picture

        # The LAST (per, kstp) pair, not the largest kstp: a model whose periods
        # each hold one time step has kstp == 0 throughout, so filtering on kstp
        # alone keeps every period and draws each barrier once per period.
        last = table.sort_values(["per", "kstp"]).iloc[-1]
        latest = table[(table["per"] == last["per"]) & (table["kstp"] == last["kstp"])]
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

    def map(self, values=None, *, per: int | None = None, layer: int = 0, **kwargs):
        """Draw the barriers on the grid, as lines over whatever the cells show.

        A barrier is an edge, so it is drawn as an edge -- over the cell field
        rather than instead of it. ``values`` colours the cells as any other map
        would; without one the grid is drawn plain beneath the barriers.

        Parameters
        ----------
        values : array-like, optional
            Per-cell values to colour the cells by.
        per : int, optional
            Which period's barriers to draw. The first defined period by default.
        layer : int, default 0
            Which layer's barriers to draw.
        **kwargs
            Forwarded to ``vor.plot.map``.

        Returns
        -------
        Choro
        """

        import plotly.graph_objects as go

        from myflopy.viz import HIGHLIGHT_WIDTH, PALETTE

        grid = self.model.vor
        picture = grid.plot.map(values, **kwargs)

        frame = self.segments(per=per)
        if len(frame):
            frame = frame[frame["layer"] == int(layer)]
        if len(frame):
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
