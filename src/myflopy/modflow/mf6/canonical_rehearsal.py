"""A recalibration rehearsal: a Ten Trails-shaped truth calibrated back through the API.

The acceptance test of the PEST preparation (ledger 589). A small
model with the SHAPE of the user's Ten Trails model -- DISV, five layers, NEWTON,
an array RCH on an average-year cycle, two RIV packages (``riv`` and ``crisp``), a
DRN, a GHB, a CHD, layered storage and two excavated lakes under one LAK package --
is built twice in one :class:`~myflopy.workspace.Project`: a START model and a
TRUTH model that differs from it by KNOWN multipliers. Synthetic data are read off
the truth run, and the start is calibrated back to the truth with the same API the
user will copy for Ten Trails (``examples/mf6/notebooks/canonical_08_recalibration_rehearsal.ipynb``).

**The truth is declared through the package-first API, never through PEST** (ledger
589): layer 2's K is the start's raster times ``mf.Polygons(value=)``
factors, recharge is a second ``year.by_month`` table, the other variants are forks
of the start's library packages. A parameter routed to the wrong package, lake,
month or layer therefore cannot cancel between the truth and the calibration, and
the truth lies INSIDE the parameter space the notebook declares, so a forward run
of the template at the truth's multipliers reproduces the noise-free data exactly
(:func:`run_template`).

What is here:

* :func:`build_rehearsal` -- the project, the two simulations, the wells, lakes,
  zones and the forecast, as a :class:`RehearsalCase`;
* :func:`synthetic_hydrographs` / :func:`synthetic_series` -- the data, dated as a
  field program reads them, off the truth run;
* :func:`truth_parameter_values` and :func:`run_template` -- the truth round trip:
  the template run at the truth's multipliers;
* :func:`recovery` and :func:`forecast_recovery` -- one run's scalars and forecast
  against the truth, by one rule; :func:`pooled_recovery` -- the verdict on several
  converged runs at once, the user's criterion for whether the truth came back (ledger
  593), judged against what an exact posterior's finite ensemble gives
  (:func:`z_tail`, :func:`band_coverage`; ledger 594).

It is demo infrastructure like :mod:`~myflopy.modflow.mf6.canonical_calibration`
and is not exported top-level. It pins its own 21x21 grid and ignores
``SIMPLE_MODFLOW_CANONICAL_PROFILE``, so the slow lane keeps it under budget.
"""

from __future__ import annotations

import ast
import copy
import json
import math
import os
import shutil
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from statistics import NormalDist
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from scipy import stats
from shapely.geometry import LineString, Point

from myflopy._logging import get_logger
from myflopy._pythonpath import pythonpath_for_children
from myflopy.cycles import AverageYear
from myflopy.layers import Basin, LayerStack, Raster
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.observations import HeadTargets, LakeStageTargets, RivFlowTargets
from myflopy.modflow.mf6.pest.calibrated import (
    CalibratedInputs,
    _complete_values,
    _fill_templates,
    _is_model_command,
)
from myflopy.modflow.mf6.pest.ies import _pest_residuals
from myflopy.modflow.mf6.pest.project import _simulation_input_files
from myflopy.modflow.mf6.pest.zones import ZoneSpec
from myflopy.package_api import (
    chd,
    disv,
    drn,
    ghb,
    gwf,
    ic,
    ims,
    lak,
    npf,
    oc,
    rch,
    riv,
    simulation,
    sto,
)
from myflopy.sources import Polygons, RasterValues
from myflopy.specs import ModelContext, ref
from myflopy.workspace import Project, Run

logger = get_logger(__name__)

#: The flow model's name in both simulations.
MODEL_NAME = "flow"
#: The two lakes, by boundname: one LAK package holds both.
LAKE_NAMES = ("north", "south")
#: The two RIV packages, as Ten Trails carries them: a main stem and a tributary.
RIVER_PACKAGES = ("riv", "crisp")
#: Calendar months whose recharge the truth raises (``"wet"``) and lowers (``"dry"``).
WET_MONTHS = (11, 12, 1, 2, 3, 4)
DRY_MONTHS = (5, 6, 7, 8, 9, 10)
#: The start's recharge, ft/d by calendar month (January = 1).
MONTHLY_RECHARGE = {
    1: 0.0045, 2: 0.0040, 3: 0.0035, 4: 0.0025, 5: 0.0012, 6: 0.0006,
    7: 0.0003, 8: 0.0003, 9: 0.0006, 10: 0.0015, 11: 0.0035, 12: 0.0045,
}
#: Rainfall on, and maximum evaporation from, the lakes, ft/d by calendar month.
#: The same in the start and the truth: the net flux through the bed is what makes a
#: lake's stage depend on its bed leakance.
LAKE_RAINFALL = {month: 3.0 * rate for month, rate in MONTHLY_RECHARGE.items()}
LAKE_EVAPORATION = {
    1: 0.0005, 2: 0.0010, 3: 0.0020, 4: 0.0030, 5: 0.0040, 6: 0.0050,
    7: 0.0060, 8: 0.0060, 9: 0.0040, 10: 0.0020, 11: 0.0010, 12: 0.0005,
}
#: The start's layered properties, top layer first (sand, silt, gravel, clay, bedrock).
UNITS = ("sand", "silt", "gravel", "clay", "bedrock")
THICKNESS = (15.0, 20.0, 45.0, 10.0, 60.0)
K = (15.0, 2.0, None, 0.2, 1.0)  # None: the gravel's K is the raster k_gravel.tif
K33 = (1.5, 0.2, 4.0, 0.02, 0.1)
SY = (0.20, 0.12, 0.25, 0.05, 0.02)
SS = 1e-5
#: The zero-based layer whose K the truth changes by zone.
K_LAYER = 2
#: Each river package's conductance (ft2/d per cell) and how far its stage is below
#: the cell's top (ft); its bottom is 1.5 ft below the stage.
RIVER_SETTINGS = {"riv": (400.0, 2.0), "crisp": (100.0, 3.0)}
#: The start's lake-bed leakance (1/d), the same for both lakes.
BED_LEAKANCE = 0.2

#: What makes the truth: multipliers on the start, keyed by the TARGET spelling
#: ``cal.parameterize`` takes. ``k`` is by zone label of :attr:`RehearsalCase.zones`
#: on layer :data:`K_LAYER`; ``recharge`` by month group; ``riv.cond`` and
#: ``lak.bedleak`` by package and by lake.
TRUTH_MULTIPLIERS: dict[str, Any] = {
    "k": {1: 2.0, 2: 0.5},
    "recharge": {"wet": 1.3, "dry": 0.7},
    "sto.sy": 0.6,
    "riv.cond": {"riv": 1.0, "crisp": 3.0},
    "lak.bedleak": {"north": 0.5, "south": 1.0},
}

#: Observation wells: name, x, y (ft, in the grid's frame) and zero-based layer;
#: ``None`` is a well SCREENED across layers 1 and 2 (its screen is set from the
#: built elevations at its cell: the middle of layer 1 to the middle of layer 2).
WELLS = (
    ("W01", 300.0, 1000.0, 0),
    ("W02", 520.0, 1200.0, 2),
    ("W03", 880.0, 1480.0, 0),
    ("W04", 850.0, 720.0, 2),
    ("W05", 1300.0, 1600.0, 0),
    ("W06", 1420.0, 780.0, 2),
    ("W07", 1720.0, 1420.0, 0),
    ("W08", 1780.0, 560.0, 2),
    ("W09", 450.0, 1700.0, None),
    ("W10", 1600.0, 320.0, None),
    ("W11", 760.0, 300.0, 0),
    ("W12", 1250.0, 1280.0, 2),
    ("W13", 820.0, 1600.0, 0),
)
#: Wells a field program missed for a while: name -> calendar months with no reading.
GAPS = {"W05": (6, 7, 8), "W11": (12,)}
#: The forecast: a head nobody measures, in layer 2, at the run's last period.
FORECAST_WELL = ("F01", 1150.0, 430.0, 2)
#: The lakes' footprints: name, centre x, y, radius (ft) and the assumed depth (ft).
LAKES = (("north", 600.0, 1550.0, 170.0, 12.0), ("south", 650.0, 520.0, 170.0, 11.0))
#: Zone boxes on layer 2: label -> (xmin, ymin, xmax, ymax); a cell is in a zone when
#: its centre is, and each zone's polygon is the union of its cells, so the polygon
#: multiplier and the pilot-point zones pick the SAME cells.
ZONE_BOXES = {1: (250.0, 900.0, 1000.0, 1850.0), 2: (1200.0, 200.0, 1900.0, 900.0)}


def ground(x: Any, y: Any) -> np.ndarray:
    """The rehearsal's land surface (ft): a valley falling east, with two incised creeks.

    Parameters
    ----------
    x, y : float or numpy.ndarray
        Coordinates in the grid's frame (ft), broadcast against each other.

    Returns
    -------
    numpy.ndarray
        The elevation at each point: 128 ft at the west edge falling to 104 ft at
        the east, less a 4 ft trough along the main stem (``y = 1050``, east of
        ``x = 900``) and a 3 ft one along the tributary (``x = 1050``, north of
        ``y = 1000``).
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return (
        128.0
        - 24.0 * x / 2100.0
        - 4.0 * np.exp(-(((y - 1050.0) / 250.0) ** 2)) * (x > 900.0)
        - 3.0 * np.exp(-(((x - 1050.0) / 200.0) ** 2)) * (y > 1000.0)
    )


@dataclass(frozen=True)
class RehearsalCase:
    """A rehearsal's project, its start and truth, and where its data come from.

    Built by :func:`build_rehearsal`; read by :func:`synthetic_hydrographs`,
    :func:`synthetic_series` and :func:`truth_parameter_values`.

    Attributes
    ----------
    workspace : pathlib.Path
        The directory everything is written under: ``gis/`` (the rasters and
        vector files the packages read), ``project/`` (the
        :class:`~myflopy.workspace.Project`) and ``data/`` (the synthetic data).
    project : Project
        The project. Its package library holds every package the two simulations
        share or differ by -- ``npf``, ``sto``, ``rch``, ``crisp`` and ``lak`` for the
        start, each with a ``<key>/truth`` variant -- and its simulation library the
        two simulations.
    start_simulation : str
        The start simulation's name, ``"start"``: what the notebook calibrates.
    truth_simulation : str
        The truth simulation's name, ``"truth"``: the start with each
        ``<key>/truth`` package swapped in (``project.simulations.fork``).
    model_name : str
        The flow model's name in both, ``"flow"``.
    year : AverageYear
        The cycle both run: one steady period, then ``cycles`` x 12 months.
    truth_multipliers : dict
        :data:`TRUTH_MULTIPLIERS`: what the truth is, by parameter target.
    wells : geopandas.GeoDataFrame
        The observation wells: ``name``, ``layer`` (``<NA>`` for a screened well),
        ``screen_top``/``screen_bottom`` (ft; ``NaN`` for a layer well), ``weight``
        (1 / the head noise the notebook draws, 0.05 ft) and the point geometry.
    forecast : HeadTargets
        The prediction of interest: the head at ``F01`` in layer 2 at the last
        period, with no value (a forecast is never fitted).
    lakes : tuple of str
        The lakes' boundnames, :data:`LAKE_NAMES`.
    crisp_cells : tuple of int
        The zero-based cells (layer 0) the ``crisp`` RIV package holds -- the gauged
        reach of its baseflow target.
    zones : pathlib.Path
        The vector file of the layer-2 K zones: layer ``zones`` (column ``zone``, 1
        and 2) and one layer per zone, ``zone_1``/``zone_2``, each the union of the
        zone's cells.
    lake_rims : dict of str to float
        Each lake's rim (ft): the water surface its assumed bed was cut from. The
        notebook's one-sided target holds each lake's stage at most at its rim.
    """

    workspace: Path
    project: Project
    start_simulation: str
    truth_simulation: str
    model_name: str
    year: AverageYear
    truth_multipliers: Mapping[str, Any]
    wells: gpd.GeoDataFrame
    forecast: HeadTargets
    lakes: tuple[str, ...]
    crisp_cells: tuple[int, ...]
    zones: Path
    lake_rims: Mapping[str, float]
    _runs: dict = field(default_factory=dict, repr=False, compare=False)

    def truth_run(self) -> Run:
        """The truth simulation's run, prepared and executed on the first call.

        The run is ``project.prepare_run("truth", "truth", overwrite=True)`` -- a
        truth run an earlier build left in the project is rebuilt -- and is kept, so
        the data and the round trip read one run.

        Returns
        -------
        Run
            The executed run.

        Raises
        ------
        RuntimeError
            When MODFLOW 6 does not terminate normally, with the end of its report.
        """

        run = self._runs.get("truth")
        if run is None:
            run = self.project.prepare_run(
                self.truth_simulation, self.truth_simulation, overwrite=True
            )
            run.execute()
            if not run.success:
                raise RuntimeError(
                    "The rehearsal's truth run failed:\n" + "\n".join(run.report[-20:])
                )
            self._runs["truth"] = run
        return run

    def truth_model(self):
        """The truth run's flow model (``truth_run().model("flow")``).

        Returns
        -------
        SimulationBase
            The run's model view.
        """

        return self.truth_run().model(self.model_name)

    def month_groups(self) -> dict[str, list[int]]:
        """The transient periods of the WET and DRY months, as ``parameterize(periods=)`` takes them.

        Returns
        -------
        dict of str to list of int
            ``{"wet": [...], "dry": [...]}``: the zero-based periods of every cycle
            whose calendar month is in :data:`WET_MONTHS` / :data:`DRY_MONTHS` -- the
            groups the truth scales recharge by. The steady period is in neither (a
            month-grouped multiplier does not reach it; ledger 549).
        """

        months = self.year.periods.dropna(subset=["month"])
        return {
            label: [int(row.per) for row in months.itertuples() if int(row.month) in chosen]
            for label, chosen in (("wet", WET_MONTHS), ("dry", DRY_MONTHS))
        }

    def zone_spec(self, model):
        """The layer-2 K zones as a ``mf.ZoneSpec`` on ``model``'s grid (labels 1 and 2, else 0).

        Parameters
        ----------
        model : SimulationBase
            The model to label the cells of -- the start run's flow model.

        Returns
        -------
        ZoneSpec
            ``mf.ZoneSpec.from_polygons(<zones layer of case.zones>, model,
            column="zone")``. A zone's polygon is the union of its cells, so the
            centroid rule picks exactly the cells the truth's ``mf.Polygons``
            multiplier takes (it covers each whole, ledger 552).
        """

        return ZoneSpec.from_polygons(gpd.read_file(self.zones, layer="zones"), model, column="zone")


def _write_raster(path: Path, values: np.ndarray, crs: Any, *, origin: float, size: float) -> None:
    with rasterio.open(
        path, "w", driver="GTiff", height=values.shape[0], width=values.shape[1], count=1,
        dtype="float64", crs=crs, transform=from_origin(origin, origin + values.shape[0] * size, size, size),
    ) as raster:
        raster.write(values, 1)


def _gis(directory: Path, grid) -> dict[str, Any]:
    """Write the rasters and vector files the packages read; what they hold."""

    directory.mkdir(parents=True, exist_ok=True)
    size, origin = 10.0, -100.0
    centres = origin + size / 2 + size * np.arange(int(2300 / size))
    x, y = np.meshgrid(centres, centres[::-1])
    _write_raster(directory / "top.tif", ground(x, y), grid.crs, origin=origin, size=size)
    gravel = 30.0 + 15.0 * np.sin(x / 600.0) * np.cos(y / 800.0)
    _write_raster(directory / "k_gravel.tif", gravel, grid.crs, origin=origin, size=size)
    lakes = gpd.GeoDataFrame(
        {
            "name": [name for name, *_ in LAKES],
            "depth": [depth for *_, depth in LAKES],
            "water_surface": [float(ground(x0, y0)) - 0.5 for _, x0, y0, _, _ in LAKES],
        },
        geometry=[Point(x0, y0).buffer(radius, 128) for _, x0, y0, radius, _ in LAKES],
        crs=grid.crs,
    )
    lakes.to_file(directory / "lakes.gpkg")
    centres_xy = np.array([(cell[1], cell[2]) for cell in grid.get_disv_gridprops()["cell2d"]])
    polygons = grid.gdf_vorPolys.geometry
    zones = directory / "k_zones.gpkg"
    rows = []
    for label, (x0, y0, x1, y1) in ZONE_BOXES.items():
        inside = np.flatnonzero(
            (centres_xy[:, 0] >= x0) & (centres_xy[:, 0] <= x1)
            & (centres_xy[:, 1] >= y0) & (centres_xy[:, 1] <= y1)
        )
        shape = polygons.iloc[inside].union_all()
        rows.append({"zone": label, "geometry": shape})
        gpd.GeoDataFrame([{"zone": label}], geometry=[shape], crs=grid.crs).to_file(
            zones, layer=f"zone_{label}"
        )
    gpd.GeoDataFrame(rows, geometry="geometry", crs=grid.crs).to_file(zones, layer="zones")
    return {"lakes": directory / "lakes.gpkg", "rims": dict(zip(lakes["name"], lakes["water_surface"], strict=True)),
            "top": directory / "top.tif", "gravel": directory / "k_gravel.tif", "zones": zones,
            "xy": centres_xy}


def _cells_on(grid, line: LineString) -> list[int]:
    polygons = grid.gdf_vorPolys
    return [int(cell) for cell in polygons.index[polygons.intersects(line)]]


def _nearest(xy: np.ndarray, x: float, y: float) -> int:
    return int(np.argmin(np.hypot(xy[:, 0] - x, xy[:, 1] - y)))


def build_rehearsal(
    workspace: str | Path,
    *,
    cycles: int = 3,
    start_date: str = "2016-10-01",
    nlay: int = 5,
    seed: int = 2026,
) -> RehearsalCase:
    """Build the rehearsal's project: a start model, a truth model and their data sources.

    The model is the canonical valley's 21x21 irregular Voronoi grid (2,100 ft
    square) carrying a Ten Trails-shaped stack of its own:

    * ``nlay`` layers from a land-surface raster (:func:`ground`), each a unit of
      :data:`UNITS` with its :data:`THICKNESS`; two lakes (:data:`LAKES`) cut into it
      by ``stack.excavate`` on an assumed ``mf.Basin`` bowl bed, under ONE LAK
      package written ``IMPLICIT``, with rainfall and evaporation;
    * NPF ``k``/``k33`` per layer (layer 2's ``k`` the raster ``k_gravel.tif``, an
      ``mf.RasterValues`` recipe), convertible everywhere, under NEWTON and an IMS
      ``COMPLEX``/``BICGSTAB`` solver;
    * STO with a ``sy`` per layer; an ARRAY RCH (``mf.rch.array``) taking each
      calendar month's rate (:data:`MONTHLY_RECHARGE`) through
      ``AverageYear.by_month``;
    * a CHD on the west edge (120 ft), a GHB on the east edge (101 ft), one DRN (a
      spring), and TWO RIV packages, ``riv`` along the valley axis and ``crisp`` down
      the tributary, both saving their flows;
    * one steady period and ``cycles`` x 12 months (``mf.AverageYear``).

    The TRUTH is the start with :data:`TRUTH_MULTIPLIERS` applied, every one
    DECLARED through a package helper (ledger 589): layer 2's K is
    ``mf.RasterValues("k_gravel.tif") * mf.Polygons(zone_1, value=2.0, default=1.0)
    * mf.Polygons(zone_2, value=0.5, default=1.0)``, recharge a second
    ``by_month`` table whose STEADY period keeps the start's annual mean (a
    month-grouped multiplier does not reach it), Sy, crisp's conductance and the
    north lake's bed leakance scaled. Each variant is registered as ``<key>/truth``
    beside the start's ``<key>``, and the truth simulation is
    ``project.simulations.fork("start", "truth", swap=...)``.

    Nothing is run here; :meth:`RehearsalCase.truth_run` runs the truth on demand.

    Parameters
    ----------
    workspace : str or pathlib.Path
        Directory to build in; created if missing. ``gis/``, ``project/`` and
        ``data/`` are written under it, and an existing ``project/`` is reused as
        the project root (its libraries are replaced: every key is registered with
        ``overwrite=True``).
    cycles : int, default 3
        Average-year cycles after the steady period, ``>= 2`` (``ValueError``
        below: the equilibrium check compares the last two). Three by default
        because two are not enough here: the steady period is not the periodic
        state, and the truth's first two cycles differ by up to 0.59 ft (a lake's
        stage; 0.54 ft a head) while its second and third differ by 0.006 ft
        (measured; ledger 589).
    start_date : str, default "2016-10-01"
        The first day of the first cycle (``mf.AverageYear(start=)``): October
        makes it a water year, as Ten Trails runs.
    nlay : int, default 5
        Layers, 3 to 5 (``ValueError`` otherwise): the first ``nlay`` units of
        :data:`UNITS`. Layer 2 carries the K zones and the screened wells reach
        layers 1-2, so fewer than 3 cannot hold them; five is Ten Trails'.
    seed : int, default 2026
        Seed of the grid's jitter (``irregular_voronoi_grid(seed=)``; 2026 is the
        canonical mesh). Every feature is placed by coordinates and resolved on
        the mesh, so another seed is another mesh under the same model.

    Returns
    -------
    RehearsalCase
        The project and what the data and the round trip need.

    Raises
    ------
    ValueError
        For ``cycles < 2`` or ``nlay`` outside 3-5.
    """

    if int(cycles) != cycles or cycles < 2:
        raise ValueError(
            f"build_rehearsal(cycles={cycles!r}): at least 2 cycles, so the equilibrium "
            "check has two to compare."
        )
    if int(nlay) != nlay or not 3 <= nlay <= len(UNITS):
        raise ValueError(
            f"build_rehearsal(nlay={nlay!r}) takes 3 to {len(UNITS)} layers: layer 2 "
            "holds the K zones and the screened wells reach layers 1-2."
        )
    workspace = Path(workspace)
    grid = irregular_voronoi_grid(nrow=21, ncol=21, cell_size=100.0, seed=seed)
    gis = _gis(workspace / "gis", grid)
    xy = gis["xy"]
    year = AverageYear(int(cycles), start=start_date, nstp=1)

    stack = LayerStack(top=Raster(gis["top"]))
    for unit, thickness in zip(UNITS[:nlay], THICKNESS[:nlay], strict=True):
        stack = stack.add(unit, thickness=thickness)
    bed = Basin(gis["lakes"], "depth", water_surface="water_surface")
    stack = stack.excavate(bottom=bed, within=gis["lakes"], name="lakes")
    result = stack.build(grid, attach=True)
    context = ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    top = np.asarray(result.top, dtype=float)
    botm = np.asarray(result.botm, dtype=float)
    active = np.asarray(result.idomain) > 0

    west = np.flatnonzero(xy[:, 0] < 110.0)
    east = np.flatnonzero(xy[:, 0] > 2000.0)
    chd_records = [[(k, int(c)), 120.0] for c in west for k in range(nlay) if active[k, c]]
    ghb_records = [[(k, int(c)), 101.0, 800.0] for c in east for k in range(min(3, nlay))
                   if active[k, c]]
    crisp_cells = _cells_on(grid, LineString([(1050.0, 2080.0), (1050.0, 1150.0)]))
    riv_cells = [cell for cell in _cells_on(grid, LineString([(1100.0, 1050.0), (2050.0, 1050.0)]))
                 if cell not in crisp_cells]
    spring = _nearest(xy, 1500.0, 1350.0)

    def river(cells, package, factor=1.0):
        cond, depth = RIVER_SETTINGS[package]
        records = [[(0, c), float(top[c] - depth), cond * factor, float(top[c] - depth - 1.5)]
                   for c in cells]
        return riv(stress_period_data={0: records}, name=package, save_flows=True)

    def conductivity(zoned: bool) -> list:
        k = []
        for layer in range(nlay):
            if K[layer] is not None:
                k.append(K[layer])
                continue
            gravel = RasterValues(gis["gravel"])
            if zoned:
                factors = TRUTH_MULTIPLIERS["k"]
                for label in sorted(factors):
                    gravel = gravel * Polygons(gis["zones"], layer=f"zone_{label}",
                                               value=factors[label], default=1.0)
            k.append(gravel)
        return k

    def recharge(truth: bool):
        if not truth:
            return rch.array(recharge=year.by_month(MONTHLY_RECHARGE))
        groups = TRUTH_MULTIPLIERS["recharge"]
        monthly = {
            month: rate * (groups["wet"] if month in WET_MONTHS else groups["dry"])
            for month, rate in MONTHLY_RECHARGE.items()
        }
        start_mean = year.by_month(MONTHLY_RECHARGE)[0]
        return rch.array(recharge=year.by_month(monthly, steady=start_mean))

    def lakes(truth: bool):
        factors = TRUTH_MULTIPLIERS["lak.bedleak"] if truth else {}
        return lak(
            context=context, nper=year.nper, lakes=gis["lakes"], lake_id_field="name",
            starting_stage={name: gis["rims"][name] - 2.0 for name in LAKE_NAMES},
            bed_leakance={name: BED_LEAKANCE * factors.get(name, 1.0) for name in LAKE_NAMES},
            rainfall=year.by_month(LAKE_RAINFALL), evaporation=year.by_month(LAKE_EVAPORATION),
            implicit=True, unexcavated="raise",
        )

    def storage(truth: bool):
        factor = TRUTH_MULTIPLIERS["sto.sy"] if truth else 1.0
        return sto(ss=SS, sy=[value * factor for value in SY[:nlay]], iconvert=1,
                   steady_state=year.steady_state, transient=year.transient)

    project = Project(workspace / "project", name="rehearsal")
    library = {
        "npf": npf(k=conductivity(False), k33=list(K33[:nlay]), icelltype=1),
        "npf/truth": npf(k=conductivity(True), k33=list(K33[:nlay]), icelltype=1),
        "sto": storage(False), "sto/truth": storage(True),
        "rch": recharge(False), "rch/truth": recharge(True),
        "crisp": river(crisp_cells, "crisp"),
        "crisp/truth": river(crisp_cells, "crisp", TRUTH_MULTIPLIERS["riv.cond"]["crisp"]),
        "lak": lakes(False), "lak/truth": lakes(True),
    }
    for key, spec in library.items():
        project.add_package(key, spec, overwrite=True)
    model = gwf(MODEL_NAME, context=context, newtonoptions="NEWTON", packages=[
        disv(result),
        ref("npf"),
        ic(strt=112.0),
        ref("sto"),
        chd(stress_period_data={0: chd_records}),
        ghb(stress_period_data={0: ghb_records}),
        drn(stress_period_data={0: [[(0, spring), float(top[spring] - 4.0), 300.0]]}),
        river(riv_cells, "riv"),
        ref("crisp"),
        ref("rch"),
        ref("lak"),
        oc(head_filerecord=f"{MODEL_NAME}.hds", budget_filerecord=f"{MODEL_NAME}.cbc",
           saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    start = simulation(model, name="start", packages=[
        year.tdis(), ims(complexity="COMPLEX", linear_acceleration="BICGSTAB")])
    project.add_simulation(start, overwrite=True)
    project.simulations.fork(
        "start", "truth", overwrite=True,
        swap={key: f"{key}/truth" for key in ("npf", "sto", "rch", "crisp", "lak")},
    )

    middles = 0.5 * (np.vstack([top[None, :], botm[:-1]]) + botm)  # each layer's mid-depth
    wells = []
    for name, x, y, layer in WELLS:
        cell = _nearest(xy, x, y)
        row = {"name": name, "layer": pd.NA if layer is None else layer,
               "screen_top": np.nan, "screen_bottom": np.nan, "weight": 1.0 / 0.05,
               "geometry": Point(x, y)}
        if layer is None:
            row["screen_top"], row["screen_bottom"] = middles[1, cell], middles[2, cell]
        wells.append(row)
    wells = gpd.GeoDataFrame(wells, geometry="geometry", crs=grid.crs)
    wells["layer"] = wells["layer"].astype("Int64")
    name, x, y, layer = FORECAST_WELL
    forecast = HeadTargets(
        locations=gpd.GeoDataFrame({"name": [name], "layer": [layer]},
                                   geometry=[Point(x, y)], crs=grid.crs),
        values=pd.DataFrame({"time": [year.nper - 1], name: [np.nan]}),
    )
    return RehearsalCase(
        workspace=workspace, project=project, start_simulation="start",
        truth_simulation="truth", model_name=MODEL_NAME, year=year,
        truth_multipliers=copy.deepcopy(TRUTH_MULTIPLIERS), wells=wells, forecast=forecast,
        lakes=LAKE_NAMES, crisp_cells=tuple(crisp_cells), zones=gis["zones"],
        lake_rims=dict(gis["rims"]),
    )


# -- the data ----------------------------------------------------------------------------

#: The two binning modes, spelled as ``HeadTargets.from_hydrographs(mode=)`` spells them.
MODES = ("average_year", "history")
#: The offset of a history-mode distractor reading (ft): a rule reading it -- a
#: period's median or mean of its two readings -- misses by half of it.
DISTRACTOR_OFFSET = 3.0


def _truth_heads(case: RehearsalCase) -> pd.DataFrame:
    """The truth run's period-end head at every well, one row per period (``per``)."""

    names = list(case.wells["name"])
    placeholder = pd.DataFrame({"time": range(case.year.nper), **{name: np.nan for name in names}})
    sampler = HeadTargets(locations=case.wells, values=placeholder)
    return sampler.simulated_heads(case.truth_model()).set_index("per")[names]


def _truth_stages(case: RehearsalCase) -> pd.DataFrame:
    """The truth run's period-end stage of every lake, one row per period (``per``)."""

    sampler = LakeStageTargets(locations={name: name for name in case.lakes})
    frame = sampler.simulated_series(case.truth_model())
    return frame.rename(columns={"time": "per"}).set_index("per")[list(case.lakes)]


def _truth_crisp_flow(case: RehearsalCase) -> pd.Series:
    """The truth run's period-end ``crisp`` exchange (MF6's sign), indexed by period."""

    sampler = RivFlowTargets(locations={"crisp": list(case.crisp_cells)}, package="crisp")
    frame = sampler.simulated_series(case.truth_model())
    return frame.set_index("time")["crisp"]


def _average_year_readings(
    case: RehearsalCase, truth: pd.DataFrame, *, years: int, readings: int, noise_sd: float,
    rng: np.random.Generator, gaps: Mapping[str, tuple[int, ...]],
) -> pd.DataFrame:
    """Dated readings whose average-year binning is the truth's last cycle, plus noise.

    Each (well, calendar month) gets ONE noise draw, so its target is the truth plus
    that draw exactly: every year's readings of the month are that value plus
    offsets symmetric about it (an odd count, so the middle one is it), so each
    year's median is the value, and so is the median over the years.
    """

    months = case.year.months_of(-1)
    first = case.year.start.year - years
    columns: dict[str, dict[pd.Timestamp, float]] = {}
    for name in truth.columns:
        readings_of: dict[pd.Timestamp, float] = {}
        for month, per in months.items():
            if month in gaps.get(name, ()):
                continue
            centre = float(truth.loc[per, name]) + (rng.normal(0.0, noise_sd) if noise_sd else 0.0)
            for year in range(first, first + years):
                days = pd.Timestamp(year=year, month=month, day=1).days_in_month
                picked = np.sort(rng.choice(np.arange(1, days + 1), size=readings, replace=False))
                spread = rng.uniform(0.2, 1.0)
                offsets = spread * np.linspace(-1.0, 1.0, readings)
                rng.shuffle(offsets)
                for day, offset in zip(picked, offsets, strict=True):
                    stamp = pd.Timestamp(year=year, month=month, day=int(day), hour=9)
                    readings_of[stamp] = centre + float(offset)
        columns[name] = readings_of
    frame = pd.DataFrame(columns).sort_index()
    frame.index.name = "date"
    return frame


def _history_readings(
    case: RehearsalCase, truth: pd.DataFrame, *, noise_sd: float, rng: np.random.Generator,
    gaps: Mapping[str, tuple[int, ...]],
) -> pd.DataFrame:
    """Readings at every period END equal to the truth there (plus noise), with distractors.

    A distractor is a reading in the MIDDLE of a transient period, the truth there
    plus :data:`DISTRACTOR_OFFSET`: the nearest-end rule never reads it, a
    period's median or mean does.
    """

    periods = case.year.periods
    columns: dict[str, dict[pd.Timestamp, float]] = {}
    for name in truth.columns:
        readings_of: dict[pd.Timestamp, float] = {}
        for row in periods.itertuples():
            month = None if pd.isna(row.month) else int(row.month)
            if month is not None and month in gaps.get(name, ()):
                continue
            value = float(truth.loc[row.per, name])
            readings_of[pd.Timestamp(row.end)] = value + (rng.normal(0.0, noise_sd) if noise_sd else 0.0)
            if not row.steady:
                middle = pd.Timestamp(row.start) + (pd.Timestamp(row.end) - pd.Timestamp(row.start)) / 2
                readings_of[middle] = value + DISTRACTOR_OFFSET
        columns[name] = readings_of
    frame = pd.DataFrame(columns).sort_index()
    frame.index.name = "date"
    return frame


def _check_mode(mode: str, years, readings_per_month) -> tuple[int, int]:
    if mode not in MODES:
        raise ValueError(f"mode={mode!r}: use one of {MODES!r} (from_hydrographs' spellings).")
    if mode == "history":
        given = [name for name, value in (("years", years), ("readings_per_month", readings_per_month))
                 if value is not None]
        if given:
            raise ValueError(
                f"{', '.join(given)} given with mode='history': history readings sit at "
                "the model's own period ends, so neither can act; leave them out."
            )
        return 0, 0
    years = 5 if years is None else years
    readings = 3 if readings_per_month is None else readings_per_month
    if int(years) != years or years < 1:
        raise ValueError(f"years={years!r}: a whole number of years, 1 or more.")
    if int(readings) != readings or readings < 1 or readings % 2 == 0:
        raise ValueError(
            f"readings_per_month={readings!r}: an ODD count, so each month's median is "
            "exactly the truth (an even count averages the middle two)."
        )
    return int(years), int(readings)


def synthetic_hydrographs(
    case: RehearsalCase,
    *,
    mode: str = "average_year",
    years: int | None = None,
    readings_per_month: int | None = None,
    noise_sd: float = 0.0,
    seed: int = 2026,
) -> Path:
    """Write the wells' dated head readings, as a field program exports them.

    A WIDE CSV -- a ``date`` column and one column per well of
    :attr:`RehearsalCase.wells` -- under ``<workspace>/data/``, read off the TRUTH
    run (:meth:`RehearsalCase.truth_run`, run here if it has not been) at each well
    as ``HeadTargets.simulated_heads`` reads it (a screened well's T-weighted
    head). :data:`GAPS` leave some wells unread in some months, every year: those
    targets are unmeasured.

    Parameters
    ----------
    case : RehearsalCase
        From :func:`build_rehearsal`.
    mode : {"average_year", "history"}, default "average_year"
        The binning the readings are written for, spelled as
        ``HeadTargets.from_hydrographs(mode=)`` spells it (another word raises
        ``ValueError``). ``"average_year"``: readings on random days of
        ``years`` past years before the model's start, ``readings_per_month`` a
        month, offset SYMMETRICALLY about one value per (well, calendar month) --
        the truth's head at that month's period of the LAST cycle, plus noise --
        so the average-year rule (each year's monthly median, then the median over
        years) gives that value EXACTLY. ``"history"``: one reading at every
        period's END (the steady period's included), equal to the truth there plus
        noise, and one DISTRACTOR in the middle of each transient period, the truth
        plus :data:`DISTRACTOR_OFFSET` -- the nearest-end rule never reads it, a
        median or mean of the period does.
    years : int, optional
        Average year only: how many calendar years of readings, ``>= 1``; default
        ``None`` means 5. Given with ``mode="history"`` it raises ``ValueError``.
    readings_per_month : int, optional
        Average year only: readings per well, month and year, an ODD number
        ``>= 1`` (``ValueError`` otherwise: an even count's median averages the
        middle two, which is not the value); default ``None`` means 3. Given with
        ``mode="history"`` it raises ``ValueError``.
    noise_sd : float, default 0.0
        Standard deviation (ft) of the normal noise added to each target value --
        one draw per (well, month) in average-year mode, per (well, period) in
        history mode; the readings around it carry it whole. ``0.0`` writes the
        noise-free truth, which the template round trip compares against. Must be
        finite and ``>= 0`` (``ValueError``).
    seed : int, default 2026
        Seed of the generator that draws the days, offsets and noise; the same
        seed writes the same file.

    Returns
    -------
    pathlib.Path
        The CSV: ``<workspace>/data/heads_<mode>_sd<noise_sd>_seed<seed>.csv``.
    """

    years, readings = _check_mode(mode, years, readings_per_month)
    if not np.isfinite(noise_sd) or noise_sd < 0:
        raise ValueError(f"noise_sd={noise_sd!r}: a finite standard deviation >= 0.")
    rng = np.random.default_rng(seed)
    truth = _truth_heads(case)
    if mode == "average_year":
        frame = _average_year_readings(case, truth, years=years, readings=readings,
                                       noise_sd=noise_sd, rng=rng, gaps=GAPS)
    else:
        frame = _history_readings(case, truth, noise_sd=noise_sd, rng=rng, gaps=GAPS)
    path = case.workspace / "data" / f"heads_{mode}_sd{noise_sd:g}_seed{seed}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.reset_index().to_csv(path, index=False)
    return path


def synthetic_series(
    case: RehearsalCase,
    *,
    stage_sd: float = 0.0,
    flow_cv: float = 0.0,
    years: int = 5,
    readings_per_month: int = 3,
    seed: int = 2026,
) -> dict[str, Path]:
    """Write the lakes' dated stage readings and crisp's gauged baseflow, off the truth run.

    Parameters
    ----------
    case : RehearsalCase
        From :func:`build_rehearsal`.
    stage_sd : float, default 0.0
        Standard deviation (ft) of the noise on each lake's monthly stage, one draw
        per (lake, calendar month); finite and ``>= 0`` (``ValueError``).
    flow_cv : float, default 0.0
        Coefficient of variation of the noise on each monthly baseflow: each value
        is the truth times ``1 + N(0, flow_cv)``; finite and ``>= 0``
        (``ValueError``).
    years : int, default 5
        Calendar years of lake-stage readings (as :func:`synthetic_hydrographs`'
        average-year mode), ``>= 1``.
    readings_per_month : int, default 3
        Lake-stage readings per lake, month and year, ODD (as
        :func:`synthetic_hydrographs`).
    seed : int, default 2026
        Seed of the generator; the same seed writes the same files.

    Returns
    -------
    dict of str to pathlib.Path
        ``"lake_stage"``: a WIDE dated CSV (``date``, ``north``, ``south``) for
        ``LakeStageTargets.from_hydrographs(mode="average_year")``, binned exactly
        as the heads are. ``"crisp_flow"``: a period-keyed CSV (``time``,
        ``crisp``) of the ``crisp`` package's exchange at the end of each period of
        the LAST cycle, in MF6's sign -- NEGATIVE, the aquifer feeding the creek
        (ledger 350) -- for ``mf.RivFlowTargets(..., package="crisp")``; a gauge
        reads a period, not a date, so it is not binned.
    """

    _, readings = _check_mode("average_year", years, readings_per_month)
    for label, value in (("stage_sd", stage_sd), ("flow_cv", flow_cv)):
        if not np.isfinite(value) or value < 0:
            raise ValueError(f"{label}={value!r}: a finite number >= 0.")
    rng = np.random.default_rng(seed)
    stages = _average_year_readings(case, _truth_stages(case), years=int(years), readings=readings,
                                    noise_sd=stage_sd, rng=rng, gaps={})
    flow = _truth_crisp_flow(case)
    last = sorted(case.year.months_of(-1).values())
    noise = rng.normal(0.0, flow_cv, size=len(last)) if flow_cv else np.zeros(len(last))
    baseflow = pd.DataFrame({"time": last, "crisp": flow.loc[last].to_numpy() * (1.0 + noise)})
    directory = case.workspace / "data"
    directory.mkdir(parents=True, exist_ok=True)
    paths = {
        "lake_stage": directory / f"lake_stage_sd{stage_sd:g}_seed{seed}.csv",
        "crisp_flow": directory / f"crisp_flow_cv{flow_cv:g}_seed{seed}.csv",
    }
    stages.reset_index().to_csv(paths["lake_stage"], index=False)
    baseflow.to_csv(paths["crisp_flow"], index=False)
    return paths


# -- the round trips ---------------------------------------------------------------------


def _metadata(cal) -> dict:
    path = Path(cal.template_workspace) / "myflopy_pest_metadata.json"
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing: build() the calibration first.")
    return json.loads(path.read_text(encoding="utf-8"))


def _truth_of_group(record: dict, group: dict | None, case: RehearsalCase) -> float:
    """The truth's value of one parameter group of one ``parameterize`` record."""

    target, truth = record["target"], case.truth_multipliers
    neutral = 0.0 if record.get("additive") else 1.0
    if target == "k":
        if group is not None and "layer" in group:  # one pilot-point set: one layer, one zone
            if group["layer"] != K_LAYER:
                return neutral
            zone = group.get("zone")
            if zone not in truth["k"]:
                raise ValueError(
                    f"Pilot-point set {group['name']!r} is zone {zone!r} of layer {K_LAYER}, "
                    f"which the truth does not know (its zones are {sorted(truth['k'])}): "
                    "zone the layer with case.zone_spec(model)."
                )
            return float(truth["k"][zone])
        layers = record.get("layers")
        if layers is not None and K_LAYER not in layers:
            return neutral
        raise ValueError(
            f"Parameter {record['name']!r} scales 'k' on layer {K_LAYER} without the truth's "
            "zones, so no value of it is the truth: zoned pilot points on that layer, "
            "zones=case.zone_spec(model)."
        )
    if target == "recharge":
        label = None if group is None else group.get("label")
        if label not in truth["recharge"]:
            raise ValueError(
                f"Parameter {record['name']!r} scales recharge over periods labelled "
                f"{label!r}; the truth scales {sorted(truth['recharge'])} -- declare "
                "periods=case.month_groups()."
            )
        return float(truth["recharge"][label])
    if target in ("sto.sy", "sy"):
        return float(truth["sto.sy"])
    if target in ("riv.cond", "riv"):
        factors = {truth["riv.cond"][package] for package in record.get("packages", [])}
        if len(factors) != 1:
            raise ValueError(
                f"Parameter {record['name']!r} shares one multiplier over packages "
                f"{record.get('packages')} the truth scales differently: name one with "
                "package=."
            )
        return float(factors.pop())
    if target in ("lak.bedleak", "bedleak"):
        label = None if group is None else group.get("label")
        if label not in truth["lak.bedleak"]:
            raise ValueError(
                f"Parameter {record['name']!r} scales bed leakance as {label!r}; the truth "
                "scales each lake on its own -- use style='feature'."
            )
        return float(truth["lak.bedleak"][label])
    return neutral


def truth_parameter_values(case: RehearsalCase, cal) -> dict[str, float]:
    """The value of every parameter of a built calibration at which it IS the truth.

    Read from the build's own record of each ``parameterize`` call (the
    ``parameters`` of ``myflopy_pest_metadata.json``: its target, packages,
    pilot-point zones and layers, period labels and lake labels), and from
    :attr:`RehearsalCase.truth_multipliers` -- never from file names, so a
    parameter the library routed to the wrong file still gets the truth's value,
    and the round trip (:func:`run_template`) shows the misrouting as a residual.

    Parameters
    ----------
    case : RehearsalCase
        The case the calibration's model is the start of.
    cal : PestProject
        The calibration, after ``build()``.

    Returns
    -------
    dict of str to float
        Every parameter of ``cal.pst`` (lower-cased name) to its truth value: a
        zone's pilot points at the zone's ``k`` factor, a month group's recharge
        multiplier at its group's factor, ``sto.sy`` at 0.6, ``riv.cond`` at its
        package's factor, a lake's ``bedleak`` at its lake's -- and every parameter
        of a target the truth does not change at 1.0 (a multiplier) or 0.0 (an
        offset).

    Raises
    ------
    ValueError
        When a parameter cannot BE the truth, naming it and the spelling that can:
        ``k`` on layer 2 without the truth's zones, recharge not grouped wet/dry,
        one ``riv.cond`` shared by ``riv`` and ``crisp``, a ``bedleak`` over both
        lakes together. Also when a parameter group of the control file is in no
        build record.
    FileNotFoundError
        When the calibration has not been built.
    """

    records = _metadata(cal)["parameters"]
    values: dict[str, float] = {}
    by_group: dict[str, float] = {}
    for record in records:
        sets = record.get("pilot_points") or []
        if sets:
            for group in sets:
                by_group[group["pargp"]] = _truth_of_group(record, group, case)
            continue
        feature_sets = record.get("feature_sets") or []
        if feature_sets:
            for group in feature_sets:
                by_group[group["pargp"]] = _truth_of_group(record, group, case)
            continue
        for group in record.get("period_groups") or [{"label": None, "pargps": record["pargps"]}]:
            value = _truth_of_group(record, group if group.get("label") else None, case)
            for pargp in group["pargps"]:
                by_group[pargp] = value
    data = cal.pst.parameter_data
    missing = sorted(set(data["pargp"]) - set(by_group))
    if missing:
        raise ValueError(f"Parameter group(s) {missing} are in no build record of {cal.template_workspace}.")
    for name, pargp in zip(data["parnme"], data["pargp"], strict=True):
        values[str(name).lower()] = by_group[pargp]
    return values


def _forward_run_without_model(forward_run: Path) -> str:
    """``forward_run.py`` with the model command removed from ``main()``, every other line kept."""

    tree = ast.parse(forward_run.read_text(encoding="utf-8"))
    main = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "main")
    kept = [statement for statement in main.body if not _is_model_command(statement)]
    if len(kept) == len(main.body):
        raise RuntimeError(f"{forward_run}'s main() never runs the model; nothing to cut.")
    main.body = kept
    return ast.unparse(tree) + "\n"


def run_template(cal, values, workspace: str | Path, *, outputs_from: Any = None) -> pd.Series:
    """Run a built calibration's template at chosen parameter values; its simulated observations.

    The truth round trip of the rehearsal, and the check of a write-back. A copy of
    the template is made in ``workspace``, every template file filled at full
    precision (``repr`` of each value, as the write-back's replay fills them --
    pyEMU's own ``write_input_files`` writes a pilot point at 4 significant
    digits), and the forward run executed there as PESTPP-IES would; each
    instruction file then reads its output.

    Parameters
    ----------
    cal : PestProject
        The calibration, after ``build()``.
    values : mapping of str to float, or CalibratedInputs
        Each parameter's value, by name (any case): :func:`truth_parameter_values`,
        or a realization's ``ies.calibrated(...)`` (its values as replayed). A FIXED
        parameter left out is at its ``parval1``; an adjustable one left out raises
        ``ValueError``.
    workspace : str or pathlib.Path
        Directory for the copy; must not exist yet or be empty (``FileExistsError``
        otherwise: nothing is deleted here).
    outputs_from : Run, str, pathlib.Path or None, optional
        ``None`` (default): MODFLOW 6 runs in the copy. A run (or its directory):
        MODFLOW 6 does NOT run; that run's OUTPUT files -- every top-level file of
        its directory that is not one of its simulation's inputs -- are copied in
        and the forward run's other steps read them. This is how a model REBUILT
        from a write-back is scored by the calibration's own readers: the inputs
        are the realization's (filled and applied here), the heads, budgets and
        stages are the rebuilt run's.

    Returns
    -------
    pandas.Series
        The simulated value of every observation of ``cal.pst``, by its lower-cased
        name.

    Raises
    ------
    RuntimeError
        When the forward run fails, with the end of its error output.
    """

    template = Path(cal.template_workspace)
    target = Path(workspace)
    if target.exists() and any(target.iterdir()):
        raise FileExistsError(f"run_template(workspace={str(target)!r}) is not empty; pass a new folder.")
    if isinstance(values, CalibratedInputs):
        values = values._values
    pst = cal.pst
    complete = _complete_values(pst, pd.Series({str(k).lower(): float(v) for k, v in values.items()}))
    shutil.copytree(template, target, dirs_exist_ok=True)
    _fill_templates(pst, complete, template, target)
    script = "forward_run.py"
    if outputs_from is not None:
        source = Path(getattr(outputs_from, "workspace", outputs_from))
        inputs = {str(path) for path in _simulation_input_files(source)}
        for item in source.iterdir():
            if item.is_file() and item.name not in inputs and not item.name.startswith("."):
                shutil.copy2(item, target / item.name)
        (target / "forward_run_outputs.py").write_text(
            _forward_run_without_model(template / "forward_run.py"), encoding="utf-8"
        )
        script = "forward_run_outputs.py"
    environment = dict(os.environ)
    pythonpath = pythonpath_for_children()
    if pythonpath is not None:
        environment["PYTHONPATH"] = pythonpath
    done = subprocess.run([sys.executable, script], cwd=target, env=environment,
                          capture_output=True, text=True, check=False)
    if done.returncode != 0:
        raise RuntimeError(
            f"The forward run in {target} failed (exit {done.returncode}):\n{done.stderr[-4000:]}"
        )
    import pyemu  # optional: myflopy imports pyEMU only where a PEST step needs it

    simulated = []
    for instruction, output in zip(pst.instruction_files, pst.output_files, strict=True):
        frame = pyemu.pst_utils.InstructionFile(str(target / instruction)).read_output_file(str(target / output))
        simulated.append(frame["obsval"].astype(float))
    result = pd.concat(simulated)
    result.index = result.index.str.lower()
    return result


#: How closely a round trip must agree (ledger 589): an observation passes
#: when ``|simulated - expected| <= ROUND_TRIP_TOLERANCE + ROUND_TRIP_RTOL * |expected|``
#: -- 1e-3 in the model's units for a head, a stage or a difference, and a part in a
#: million of a flow's magnitude. Measured on the rehearsal: the truth round trip agrees
#: to 4e-8 (flows) and 2e-10 ft; a model rebuilt from a write-back to 8e-7 ft and
#: 2.7e-7 of a flow, which is FloPy writing list values at 9 significant digits and
#: the solver's closure, not a mapping (a misrouted parameter misses by feet).
ROUND_TRIP_TOLERANCE = 1e-3
ROUND_TRIP_RTOL = 1e-6


def round_trip(cal, simulated: pd.Series, expected: pd.Series | None = None) -> pd.DataFrame:
    """How far a set of simulated observations is from what it should reproduce.

    Two uses. With ``expected=None`` -- the TRUTH round trip -- every MEASURED
    observation of ``cal`` (non-zero weight) against its ``obsval``, scored as PEST++
    scores it through the review's one rule (``ies._pest_residuals``, ledger 564): a
    one-sided bound is 0 wherever it holds. With ``expected`` -- a WRITE-BACK check --
    every observation ``simulated`` holds against ``expected`` (the base realization's
    own values, :func:`base_observations`), bounds and unmeasured rows included.

    Parameters
    ----------
    cal : PestProject
        The calibration, after ``build()``.
    simulated : pandas.Series
        Simulated values by lower-cased observation name -- :func:`run_template`'s
        result.
    expected : pandas.Series or None, default None
        What ``simulated`` should equal, by lower-cased observation name. ``None``:
        the control file's measured values. An observation ``simulated`` holds and
        ``expected`` lacks raises ``KeyError``.

    Returns
    -------
    pandas.DataFrame
        One row per observation compared, indexed by its lower-cased name: ``group``
        (``obgnme``), ``expected``, ``simulated``, ``error`` (its absolute value as
        PEST++ would score it), ``allowed`` (``ROUND_TRIP_TOLERANCE +
        ROUND_TRIP_RTOL * |expected|``) and ``ratio`` (``error / allowed``: the
        round trip holds where every ratio is at most 1).

    Raises
    ------
    KeyError
        When a measured observation (``expected=None``) has no simulated value, or
        an observation of ``simulated`` has no expected one.
    """

    pst = cal.pst
    data = pst.observation_data.copy()
    data.index = data.index.str.lower()
    if expected is None:
        chosen = data.index[data["weight"].astype(float) > 0]
        missing = sorted(set(chosen) - set(simulated.index))
        if missing:
            raise KeyError(f"No simulated value for measured observation(s) {missing[:5]}.")
        values = simulated.loc[chosen].astype(float)
        residual, _ = _pest_residuals(pst, values)  # pyEMU names are lower-case
        reference = data.loc[chosen, "obsval"].astype(float)
    else:
        chosen = simulated.index
        missing = sorted(set(chosen) - set(expected.index))
        if missing:
            raise KeyError(f"No expected value for observation(s) {missing[:5]}.")
        values = simulated.astype(float)
        reference = expected.loc[chosen].astype(float)
        residual = values - reference
    allowed = ROUND_TRIP_TOLERANCE + ROUND_TRIP_RTOL * reference.abs()
    error = residual.abs()
    return pd.DataFrame({
        "group": data.loc[chosen, "obgnme"], "expected": reference, "simulated": values,
        "error": error, "allowed": allowed, "ratio": error / allowed,
    })


def base_observations(ies) -> pd.Series:
    """The posterior ``base`` realization's simulated observations, at full precision.

    Read from PESTPP-IES's residual file of the posterior iteration,
    ``<case>.<iteration>.base.rei`` (its ``modelled`` column, 16 digits; the
    observation ensemble holds 6).

    Parameters
    ----------
    ies : IesResults
        The run (``cal.run_ies(...)``).

    Returns
    -------
    pandas.Series
        Every observation's simulated value by lower-cased name.

    Raises
    ------
    FileNotFoundError
        When the run wrote no base residual file (``ies_include_base`` off).
    """

    import pyemu  # optional: myflopy imports pyEMU only where a PEST step needs it

    path = Path(ies.workspace) / f"{ies.case}.{ies.posterior_iteration}.base.rei"
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing: the run has no 'base' realization.")
    frame = pyemu.pst_utils.read_resfile(str(path))
    values = frame["modelled"].astype(float)
    values.index = values.index.str.lower()
    return values


#: The rehearsal's per-run stochastic pass criteria (ledger 589; the user's Q9).
PHI_RATIO_MAX = 0.2
GAP_CLOSED_MIN = 0.5

#: The two-sided 90% point of a standard normal, ``NormalDist().inv_cdf(0.95)``: a
#: scalar's 5-95% band is its posterior log10 mean +- this many log10 sds, and its
#: truth is inside exactly when ``|z|`` is at most this (ledger 593).
Z_BAND = NormalDist().inv_cdf(0.95)
#: How often each pooled criterion may fail an exact, converged posterior, by the
#: independence model :func:`pooled_recovery` computes it with (ledger 593).
#: The coverage entry is set low because a run's z values are NOT independent: with
#: the rehearsal's posterior correlations (rch_wet with Sy 0.68, with rch_dry 0.61,
#: measured over six converged runs) the same threshold fails an exact posterior more
#: often than the independent count says (ledger 594 has the Monte Carlo).
POOLED_FALSE_FAIL = {"coverage": 0.01, "largest_z": 0.01, "forecast": 0.015}


def _z_scale(members: int) -> float:
    """How much wider than a standard normal an exact posterior's ``z`` is.

    ``z`` measures a truth against the MEAN and SD of ``members`` realizations, not the
    posterior's own: for an exact posterior the truth is one more draw, so ``z`` is
    ``sqrt(1 + 1/members)`` times a Student t with ``members - 1`` degrees of freedom
    (ledger 594). This is that factor; the t is in :func:`z_tail`.
    """

    return math.sqrt(1.0 + 1.0 / members)


def z_tail(z: float, members: int) -> float:
    """How often an EXACT posterior's ``z`` is at least ``|z|`` away from 0.

    Parameters
    ----------
    z : float
        A truth's ``z`` (:func:`recovery`), the distance from the posterior
        realizations' log10 mean in their log10 sds. ``NaN`` reads as infinitely far.
    members : int
        How many realizations the mean and sd were taken over, at least 2.

    Returns
    -------
    float
        The two-sided tail probability ``P(|Z| >= |z|)`` for ``Z = sqrt(1 + 1/members)
        * T(members - 1)`` -- the ``z`` of one more draw from the posterior the
        realizations sample (ledger 594). At 63 realizations ``|z| = 1.645``
        is 0.108, not the normal's 0.100.

    Raises
    ------
    ValueError
        When ``members`` is below 2 (a sd needs two realizations).
    """

    if int(members) < 2:
        raise ValueError(f"a z needs at least 2 realizations, not {members}.")
    magnitude = math.inf if math.isnan(z) else abs(float(z))
    return float(2.0 * stats.t.sf(magnitude / _z_scale(int(members)), int(members) - 1))


def band_coverage(members: int) -> float:
    """The probability an EXACT posterior's 5-95% band, fitted to ``members``
    realizations, holds its truth: ``1 - z_tail(Z_BAND, members)``.

    0.9 only in the limit: 0.892 at 63 realizations, 0.866 at 15 (ledger
    594), because the band's mean and sd are themselves estimated.

    Parameters
    ----------
    members : int
        How many posterior realizations the band was fitted to, at least 2.

    Returns
    -------
    float
        The band's coverage of a truth drawn from the posterior the realizations sample.
    """

    return 1.0 - z_tail(Z_BAND, members)


def _realizations(ensemble) -> pd.DataFrame:
    """An IES ensemble's realizations, column names lower-cased: never its ``base`` row."""

    frame = ensemble._df.rename(columns=str.lower)
    return frame.drop(index=[label for label in frame.index if str(label).lower() == "base"])


def recovery(case: RehearsalCase, cal, ies) -> pd.DataFrame:
    """How the SCALAR parameters came back: truth, prior and posterior 5-95% bands.

    A scalar is a parameter group of ONE parameter -- every group but the pilot
    points -- read in its log10 space (each is a log-transformed multiplier, which
    PESTPP-IES draws and updates as a Gaussian in log space). **One rule** (ledger
    593): an ensemble's band is the normal fitted to its REALIZATIONS'
    log10 values -- the ``base`` realization is not a draw and is left out of every
    statistic -- so a band is its log10 mean +- :data:`Z_BAND` sds and ``inside`` is
    ``|z| <= Z_BAND``: the two cannot disagree. The fitted band rather than the
    ensemble's own 5th/95th percentiles, because 64 members put each percentile
    between two members (its sampling error is several times the mean's), and because
    the pooled criterion's largest ``|z|`` reaches past the ensemble's range, where a
    percentile has nothing to say.

    Parameters
    ----------
    case : RehearsalCase
        The case.
    cal : PestProject
        The calibration that was run.
    ies : IesResults
        Its run.

    Returns
    -------
    pandas.DataFrame
        One row per scalar parameter group (``pargp``): ``truth``, ``prior_p05``,
        ``prior_p95``, ``posterior_p05``, ``posterior_p95`` (the bands, as above, of
        the prior and posterior ensembles), ``posterior_mean`` (the GEOMETRIC mean of
        the posterior realizations), ``base`` (the posterior base realization's value,
        ``NaN`` for a run without one), ``inside`` (``|z| <= Z_BAND``: the truth within
        the posterior band), ``narrowed`` (the posterior's log10 sd below the prior's),
        ``gap_closed`` -- how much of the log distance from the start (1.0) to the
        truth the posterior mean covered, ``NaN`` where the truth is the start -- and
        ``z``, the truth's distance from the posterior log10 mean in posterior log10
        sds (``+-inf`` when every realization holds one value and the truth another,
        ``0.0`` when it is that value), and ``members``, how many posterior
        realizations the band and ``z`` were taken over -- what :func:`pooled_recovery`
        needs to know how far an exact posterior's ``z`` strays (:func:`z_tail`).
    """

    truth = truth_parameter_values(case, cal)
    data = cal.pst.parameter_data
    groups = data.groupby("pargp")["parnme"].apply(list)
    prior = _realizations(ies.prior_parameters)
    whole = ies.posterior_parameters._df.rename(columns=str.lower)
    posterior = _realizations(ies.posterior_parameters)
    bases = [label for label in whole.index if str(label).lower() == "base"]
    rows = []
    for group, names in groups.items():
        if len(names) != 1:
            continue
        name = str(names[0]).lower()
        log_post = np.log10(posterior[name].astype(float))
        log_prior = np.log10(prior[name].astype(float))
        mean, spread = float(log_post.mean()), float(log_post.std(ddof=1))
        prior_mean, prior_spread = float(log_prior.mean()), float(log_prior.std(ddof=1))
        value = truth[name]
        gap = float(np.log10(value))
        if spread > 0:
            z = (gap - mean) / spread
        else:
            z = 0.0 if gap == mean else float(np.copysign(np.inf, gap - mean))
        rows.append({
            "pargp": group, "truth": value,
            "prior_p05": 10 ** (prior_mean - Z_BAND * prior_spread),
            "prior_p95": 10 ** (prior_mean + Z_BAND * prior_spread),
            "posterior_p05": 10 ** (mean - Z_BAND * spread),
            "posterior_p95": 10 ** (mean + Z_BAND * spread),
            "posterior_mean": 10 ** mean,
            "base": float(whole.loc[bases[0], name]) if bases else np.nan,
            "inside": bool(abs(z) <= Z_BAND), "narrowed": bool(spread < prior_spread),
            "gap_closed": float(mean / gap) if gap != 0 else np.nan,
            "z": float(z), "members": int(log_post.size),
        })
    return pd.DataFrame(rows).set_index("pargp")


def forecast_recovery(ies, simulated: pd.Series, name: str | int = 0) -> dict:
    """Whether a forecast's truth is inside its posterior ensemble's range.

    The range is the posterior REALIZATIONS' minimum to maximum -- the ``base``
    realization left out, as in :func:`recovery` -- so for an exact posterior of
    ``members`` exchangeable realizations the truth falls outside it with probability
    ``2 / (members + 1)``, which is what :func:`pooled_recovery` allows for.

    Parameters
    ----------
    ies : IesResults
        The run (``cal.run_ies(...)``).
    simulated : pandas.Series
        The forecast's TRUE value, read by the observation's full name: the template
        run at the truth's multipliers (:func:`run_template` with
        :func:`truth_parameter_values`), noise-free.
    name : str or int, default 0
        Which forecast, as :meth:`IesResults.forecast` takes it: a name (exact, else
        a unique case-insensitive substring) or a position in its ``forecast_names``.

    Returns
    -------
    dict
        ``name`` (the full observation name), ``truth``, ``posterior_low`` /
        ``posterior_high`` and ``prior_low`` / ``prior_high`` (each ensemble's
        realizations' range), ``members`` (the number of posterior realizations) and
        ``covered`` (``posterior_low <= truth <= posterior_high``).
    """

    forecast = ies.forecast(name).name
    key = forecast.lower()
    posterior = _realizations(ies.posterior)[key].astype(float)
    prior = _realizations(ies.prior)[key].astype(float)
    truth = float(simulated[forecast])
    low, high = float(posterior.min()), float(posterior.max())
    return {"name": forecast, "truth": truth, "posterior_low": low, "posterior_high": high,
            "prior_low": float(prior.min()), "prior_high": float(prior.max()),
            "members": int(posterior.size), "covered": bool(low <= truth <= high)}


def _fewer_than(misses: Sequence[float]) -> list[float]:
    """``P(covered < m)`` for ``m = 0..len(misses)``: independent trials (a truth, a
    run's forecast), each missing with its own probability -- a Poisson-binomial, by
    convolution."""

    distribution = [1.0]
    for miss in misses:
        covered = 1.0 - miss
        nxt = [0.0] * (len(distribution) + 1)
        for count, probability in enumerate(distribution):
            nxt[count] += probability * miss
            nxt[count + 1] += probability * covered
        distribution = nxt
    return [float(sum(distribution[:m])) for m in range(len(distribution))]


@dataclass(frozen=True)
class PooledRecovery:
    """The pooled verdict on several rehearsal runs' recovery (:func:`pooled_recovery`).

    Attributes
    ----------
    runs : int
        How many runs were pooled.
    values : int
        How many scalar truths they hold in all (runs x scalars a run).
    inside : int
        How many of those sit inside their posterior 5-95% band (``|z| <= Z_BAND``).
    inside_min : int
        The fewest that pass: the largest count an exact posterior falls BELOW with
        probability at most ``POOLED_FALSE_FAIL["coverage"]`` -- each truth inside
        with its band's own coverage, :func:`band_coverage` of its run's realizations
        (0.892 at 63, not 0.9: the band's mean and sd are estimated; ledger
        594).
    largest_z : float
        The ``|z|`` furthest out for an exact posterior: the one with the smallest
        :func:`z_tail` -- the largest ``|z|`` when every run has as many realizations.
    largest_z_at : tuple of (int, str)
        Where it is: the run's position in the reports passed and the parameter group.
    z_max : float
        The largest ``|z|`` that passes at ``largest_z``'s run: the bound each of the
        ``values`` truths' ``z`` stays inside with probability ``(1 -
        POOLED_FALSE_FAIL["largest_z"]) ** (1 / values)`` (Sidak), read off the
        distribution an exact posterior's ``z`` has at that run's realizations
        (:func:`z_tail`): 3.88 at 63, where a standard normal's would be 3.63. Every
        ``|z|`` is judged against its own run's bound, and ``largest_z`` is the one
        nearest failing it.
    forecasts_covered : int
        How many runs' forecast truth is inside their posterior range.
    forecasts_min : int
        The fewest that pass: the largest count an exact posterior falls below with
        probability at most ``POOLED_FALSE_FAIL["forecast"]``, each run missing with
        probability ``2 / (members + 1)``.
    false_fail : dict of str to float
        How often an exact, converged posterior fails each criterion at these
        thresholds (``coverage``, ``largest_z``, ``forecast``) and, as their sum, at
        most how often it fails the verdict (``total``). The coverage rate assumes the
        scalars' z are independent; within one run they are correlated through the
        posterior, which raises it -- measured by the reproducer's
        ``correlated_false_fail`` (ledger 594) -- and is why its budget is
        the smallest.
    passes : dict of str to bool
        Each criterion's verdict: ``coverage`` (``inside >= inside_min``),
        ``largest_z`` (``largest_z <= z_max``), ``forecast``
        (``forecasts_covered >= forecasts_min``).
    """

    runs: int
    values: int
    inside: int
    inside_min: int
    largest_z: float
    largest_z_at: tuple[int, str]
    z_max: float
    forecasts_covered: int
    forecasts_min: int
    false_fail: Mapping[str, float]
    passes: Mapping[str, bool]

    @property
    def fraction(self) -> float:
        """The fraction of scalar truths inside their posterior 5-95% bands."""

        return self.inside / self.values

    @property
    def passed(self) -> bool:
        """Every pooled criterion passes."""

        return all(self.passes.values())

    def summary(self) -> pd.DataFrame:
        """One row per criterion: what was measured, the threshold, the verdict and
        how often an exact posterior fails it."""

        return pd.DataFrame({
            "measured": [f"{self.inside} of {self.values} ({self.fraction:.0%})",
                         f"{self.largest_z:.2f} (run {self.largest_z_at[0]}, "
                         f"{self.largest_z_at[1]})",
                         f"{self.forecasts_covered} of {self.runs}"],
            "threshold": [f">= {self.inside_min}", f"<= {self.z_max:.2f}",
                          f">= {self.forecasts_min}"],
            "passes": [self.passes["coverage"], self.passes["largest_z"],
                       self.passes["forecast"]],
            "false_fail": [self.false_fail["coverage"], self.false_fail["largest_z"],
                           self.false_fail["forecast"]],
        }, index=pd.Index(["truths inside their 5-95% bands", "largest |z|",
                           "forecast covered"], name="criterion"))


def pooled_recovery(reports: Sequence[Mapping]) -> PooledRecovery:
    """Judge the rehearsal's recovery over several CONVERGED runs at once.

    The user's criterion (ledger 593), in place of "every scalar truth
    inside its 5-95% band in every run", which an exact posterior meets only
    ``0.9 ** 6`` (53%) of the time with six scalars. Pooled over every run: the
    fraction of scalar truths inside their bands, the largest ``|z|``, and how many
    runs' forecast range holds its truth -- each against the threshold an EXACT,
    converged posterior fails with probability at most :data:`POOLED_FALSE_FAIL`'s
    entry (computed for the counts passed, so it holds for any number of runs and
    realizations). An exact posterior's ``z`` is measured against ESTIMATED moments --
    the truth is one more draw beside the realizations -- so it is ``sqrt(1 + 1/n)``
    times a Student t with ``n - 1`` degrees of freedom, not a standard normal
    (:func:`z_tail`; ledger 594): a normal null fails an exact posterior
    4.6% of the time at the thresholds it gives, not the 3.1% it claims. At six runs of
    six scalars at 63 realizations the thresholds are at least 27 of 36 inside, every
    ``|z|`` at most 3.88 and the forecast covered in 5 of 6: the three fail an exact
    posterior 2.7% of the time if a run's z were independent, and 3.4% with the
    correlation the rehearsal's posteriors carry. Run it on runs whose ensembles have
    converged: an unconverged one is over-confident (ledger 590) and fails it as it
    should.

    Parameters
    ----------
    reports : sequence of mapping
        One ``rehearsal_report`` per run, as
        ``canonical_08_recalibration_rehearsal.ipynb`` builds it, from independent
        seeds: each needs ``"scalars"`` (parameter group -> a mapping holding
        ``"z"`` and ``"members"``, from :func:`recovery`) and ``"forecast"`` (the mapping
        :func:`forecast_recovery` returns: ``"covered"`` and ``"members"``). A
        report's ``inside`` is re-derived from its ``z`` (``|z| <= Z_BAND``), the one
        rule. A ``z`` that is ``NaN`` is outside.

    Returns
    -------
    PooledRecovery
        The counts, the thresholds, each criterion's verdict and how often an exact
        posterior fails each; ``.passed`` is the verdict, ``.summary()`` a table.

    Raises
    ------
    ValueError
        When ``reports`` is empty, a report holds no scalars, or a scalar's
        ``"members"`` is below 2.
    KeyError
        When a report lacks ``"scalars"``, ``"forecast"``, a scalar's ``"z"`` or
        ``"members"`` or the forecast's ``"covered"``/``"members"`` -- a report from
        before ledger 593 carries ``forecast_truth``/``forecast_posterior``
        instead, and one from before ledger 594 no scalar ``"members"``.
    """

    reports = list(reports)
    if not reports:
        raise ValueError("pooled_recovery needs at least one run's report.")
    zs: list[tuple[int, str, float, int]] = []
    misses: list[float] = []
    covered = 0
    for index, report in enumerate(reports):
        for key in ("scalars", "forecast"):
            if key not in report:
                raise KeyError(f"report {index} has no {key!r}: build it as "
                               "canonical_08_recalibration_rehearsal.ipynb does "
                               "(recovery() and forecast_recovery()).")
        if not report["scalars"]:
            raise ValueError(f"report {index} holds no scalar parameters.")
        for group, row in report["scalars"].items():
            if "members" not in row:
                raise KeyError(f"report {index}'s scalar {group!r} has no 'members': take "
                               "the scalars' columns from recovery() with 'members' "
                               "(ledger 594).")
            zs.append((index, str(group), float(row["z"]), int(row["members"])))
        forecast = report["forecast"]
        covered += bool(forecast["covered"])
        misses.append(2.0 / (int(forecast["members"]) + 1))

    values = len(zs)
    magnitudes = [abs(z) if not math.isnan(z) else math.inf for _, _, z, _ in zs]
    tails = [z_tail(z, members) for _, _, z, members in zs]
    inside = sum(magnitude <= Z_BAND for magnitude in magnitudes)
    budget = POOLED_FALSE_FAIL
    # Coverage: each truth inside with its own band's coverage (a Poisson-binomial).
    fewer_inside = _fewer_than([1.0 - band_coverage(members) for *_, members in zs])
    inside_min = max(k for k, rate in enumerate(fewer_inside) if rate <= budget["coverage"])
    # Largest |z|: each truth's own tail against one per-value budget (Sidak).
    per_value = 1.0 - (1.0 - budget["largest_z"]) ** (1.0 / values)
    worst = min(range(values), key=lambda i: tails[i])
    members = zs[worst][3]
    z_max = _z_scale(members) * float(stats.t.isf(per_value / 2.0, members - 1))
    fewer = _fewer_than(misses)
    forecasts_min = max(m for m, rate in enumerate(fewer) if rate <= budget["forecast"])
    false_fail = {"coverage": fewer_inside[inside_min], "largest_z": budget["largest_z"],
                  "forecast": fewer[forecasts_min]}
    false_fail["total"] = sum(false_fail.values())
    return PooledRecovery(
        runs=len(reports), values=values, inside=inside, inside_min=inside_min,
        largest_z=magnitudes[worst], largest_z_at=(zs[worst][0], zs[worst][1]), z_max=z_max,
        forecasts_covered=covered, forecasts_min=forecasts_min, false_fail=false_fail,
        passes={"coverage": inside >= inside_min, "largest_z": tails[worst] >= per_value,
                "forecast": covered >= forecasts_min},
    )
