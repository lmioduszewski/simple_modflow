"""Canonical MODFLOW 6 PRT workflow.

PRT is a separate MF6 simulation that consumes a completed groundwater-flow
model's head and budget outputs through FMI. The resulting track CSV is exposed
as view nouns on :class:`PRTRunResults` (``pathlines``/``travel_time``/
``endpoints``/``capture``) and feeds the same plotting and 3D visualization layer
as other particle-tracking engines.
"""

from __future__ import annotations

import os
import re
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any

import flopy
import pandas as pd

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase


def _mfdata_value(value, default=None):
    """Unwrap a FloPy MFData value (calling ``get_data()`` if present), or return ``default``."""

    if value is None:
        return default
    getter = getattr(value, "get_data", None)
    return getter() if callable(getter) else value


def _default_prt_name(flow_name: str) -> str:
    """A default PRT model name derived from the flow model name (sanitized, ``_prt`` suffix)."""

    stem = re.sub(r"[^A-Za-z0-9_]", "_", str(flow_name))
    return f"{stem[:12]}_prt"


def _output_filename(model: SimulationBase, kind: str) -> str:
    """The GWF head/budget output filename from its OC record (falling back to ``<name>.hds/.cbc``)."""

    oc = getattr(model.gwf, "oc", None)
    record = getattr(oc, f"{kind}_filerecord", None) if oc is not None else None
    data = _mfdata_value(record)
    if data is not None and len(data):
        value = data[0][0]
        return value.decode() if isinstance(value, bytes) else str(value)
    suffix = "hds" if kind == "head" else "cbc"
    return f"{model.name}.{suffix}"


def _grid_filename(model: SimulationBase) -> str:
    """Return the GWF binary-grid filename required by PRT local-z releases."""

    grid = getattr(model.gwf, "dis", None) or getattr(model.gwf, "disv", None) or getattr(model.gwf, "disu", None)
    record = getattr(grid, "grb_filerecord", None) if grid is not None else None
    data = _mfdata_value(record)
    if data is not None and len(data):
        value = data[0][0]
        return value.decode() if isinstance(value, bytes) else str(value)
    package_type = getattr(grid, "package_type", "disv")
    return f"{model.name}.{package_type}.grb"


def _tdis_perioddata(model: SimulationBase) -> list[tuple[float, int, float]]:
    """The flow model's TDIS ``(perlen, nstp, tsmult)`` rows (defaulting to one steady period)."""

    tdis = getattr(model.sim, "tdis", None)
    data = _mfdata_value(getattr(tdis, "perioddata", None))
    if data is None:
        return [(1.0, 1, 1.0)]
    return [
        (float(row["perlen"]), int(row["nstp"]), float(row["tsmult"]))
        for row in data
    ]


def _tdis_time_units(model: SimulationBase) -> str:
    """The flow model's TDIS time units (``"DAYS"`` if unset)."""

    tdis = getattr(model.sim, "tdis", None)
    return str(_mfdata_value(getattr(tdis, "time_units", None), "DAYS"))


def _copy_grid_to_prt(flow_model: SimulationBase, prt_model):
    """Replicate the GWF DIS/DISV discretization onto a PRT model (raises for unsupported grids)."""

    gwf = flow_model.gwf
    if getattr(gwf, "disv", None) is not None:
        disv = gwf.disv
        return flopy.mf6.ModflowPrtdisv(
            prt_model,
            length_units=_mfdata_value(disv.length_units),
            xorigin=_mfdata_value(disv.xorigin),
            yorigin=_mfdata_value(disv.yorigin),
            angrot=_mfdata_value(disv.angrot),
            nlay=int(_mfdata_value(disv.nlay)),
            ncpl=int(_mfdata_value(disv.ncpl)),
            nvert=int(_mfdata_value(disv.nvert)),
            top=_mfdata_value(disv.top),
            botm=_mfdata_value(disv.botm),
            idomain=_mfdata_value(disv.idomain),
            vertices=_mfdata_value(disv.vertices),
            cell2d=_mfdata_value(disv.cell2d),
        )
    if getattr(gwf, "dis", None) is not None:
        dis = gwf.dis
        return flopy.mf6.ModflowPrtdis(
            prt_model,
            length_units=_mfdata_value(dis.length_units),
            xorigin=_mfdata_value(dis.xorigin),
            yorigin=_mfdata_value(dis.yorigin),
            angrot=_mfdata_value(dis.angrot),
            nlay=int(_mfdata_value(dis.nlay)),
            nrow=int(_mfdata_value(dis.nrow)),
            ncol=int(_mfdata_value(dis.ncol)),
            delr=_mfdata_value(dis.delr),
            delc=_mfdata_value(dis.delc),
            top=_mfdata_value(dis.top),
            botm=_mfdata_value(dis.botm),
            idomain=_mfdata_value(dis.idomain),
        )
    raise NotImplementedError("Canonical PRT currently supports GWF DIS and DISV grids.")


@dataclass(frozen=True)
class PRTReleasePoints:
    """The set of particle starting locations for a PRT (particle-tracking) run.

    Holds the PRP (particle release point) ``packagedata`` -- one record per
    released particle, giving its id and starting cell/coordinates -- in the form
    MF6's PRT model expects. Build it with the :meth:`from_cells` constructor
    (release particles at given cells) rather than assembling records by hand, then
    hand it to :class:`PRTProject` / ``ParticleTracking.prt(...)``.

    Attributes
    ----------
    packagedata
        Normalized PRP release records, one tuple per particle. Records are
        5-tuples ``(irpt, cellid, x, y, local_z)``, or 6-tuples ending in a
        release-group boundname when ``group=`` was given.
    """

    packagedata: tuple[tuple[Any, ...], ...]

    @property
    def has_groups(self) -> bool:
        """Whether any record carries a release-group boundname (a 6-tuple)."""

        return any(len(row) >= 6 for row in self.packagedata)

    @property
    def groups(self) -> tuple[str, ...]:
        """The distinct release-group labels, in first-seen order (empty when ungrouped).

        These are the labels MF6 echoes -- **uppercased** -- into the track CSV's
        ``name`` column, which is what ``results.capture`` groups by.
        """

        seen: list[str] = []
        for row in self.packagedata:
            if len(row) >= 6 and row[5] not in seen:
                seen.append(str(row[5]))
        return tuple(seen)

    @staticmethod
    def _resolve_groups(group, count: int) -> list[str] | None:
        """Broadcast ``group`` to one label per release point (``None`` stays ungrouped)."""

        if group is None:
            return None
        if isinstance(group, str):
            return [group] * count
        labels = [str(value) for value in group]
        if len(labels) != count:
            raise ValueError(
                f"group must be one label or one per release point; got {len(labels)} "
                f"labels for {count} points."
            )
        return labels

    @classmethod
    def merge(cls, *sets: PRTReleasePoints) -> PRTReleasePoints:
        """Combine release-point sets into one PRP, renumbering ``irpt`` sequentially.

        The practical way to build a multi-group release: make one grouped set per
        capture zone with :meth:`from_cells`, then merge them into the single PRP
        package MF6 wants. Group labels ride along untouched.

        Examples
        --------
        >>> west = mf.PRTReleasePoints.from_cells(model, [10, 11], group="west_wells")
        >>> east = mf.PRTReleasePoints.from_cells(model, [80, 81], group="east_wells")
        >>> rp = mf.PRTReleasePoints.merge(west, east)
        """

        rows: list[tuple[Any, ...]] = []
        for release_set in sets:
            for row in release_set.packagedata:
                rows.append((len(rows), *tuple(row)[1:]))
        return cls(tuple(rows))

    @classmethod
    def from_cells(
        cls,
        model: SimulationBase,
        cells: Sequence[int | tuple[int, ...]],
        *,
        layer: int = 0,
        local_z: float = 0.5,
        group: str | Sequence[str] | None = None,
    ) -> PRTReleasePoints:
        """Create particle release points at the centers of given model cells.

        Parameters
        ----------
        model : SimulationBase
            The flow model whose grid supplies the cell centers.
        cells : Sequence[int or tuple]
            Cell ids to release from -- zero-based DISV cell numbers, or structured
            ``(row, col)`` / ``(layer, row, col)`` tuples.
        layer : int, default 0
            Layer to place the particles in when ``cells`` are bare cell/``(row, col)``.
        local_z : float, default 0.5
            Vertical position within the cell (0 = bottom, 1 = top).
        group : str or Sequence[str], optional
            Release-group label(s) written as PRP boundnames -- one label for every
            point, or one per cell. MF6 echoes them (uppercased) into the track CSV,
            which is what ``results.capture.map()`` groups particles by.

        Returns
        -------
        PRTReleasePoints
            Pass it to ``model.particle_tracking.prt(release_points=...)``.

        Examples
        --------
        >>> rp = mf.PRTReleasePoints.from_cells(model, cells=[100, 120, 140],
        ...                                     layer=0, local_z=0.5)
        >>> west = mf.PRTReleasePoints.from_cells(model, cells=[100, 120],
        ...                                       group="west_wells")
        """

        cells = list(cells)
        groups = cls._resolve_groups(group, len(cells))
        grid = model.gwf.modelgrid
        rows = []
        for irpt, cell in enumerate(cells):
            if isinstance(cell, tuple):
                if len(cell) == 2:
                    cellid = (int(layer), int(cell[0]), int(cell[1]))
                elif len(cell) == 3:
                    cellid = tuple(int(value) for value in cell)
                else:
                    raise ValueError("Structured PRT cell tuples must be (row, col) or (layer, row, col).")
                row, col = cellid[-2:]
                x = grid.xcellcenters[row, col]
                y = grid.ycellcenters[row, col]
            else:
                cell = int(cell)
                cellid = (int(layer), cell)
                x = grid.xcellcenters[cell]
                y = grid.ycellcenters[cell]
            record: tuple[Any, ...] = (
                irpt,
                cellid,
                float(x),
                float(y),
                float(local_z),
            )
            rows.append(record if groups is None else (*record, groups[irpt]))
        return cls(tuple(rows))

    @classmethod
    def from_points(
        cls,
        model: SimulationBase,
        points,
        *,
        layer: int = 0,
        local_z: float = 0.5,
        group: str | Sequence[str] | None = None,
    ) -> PRTReleasePoints:
        """Create release points from point geometries or an iterable of ``(x, y)``.

        ``group`` labels the points as one release group (or one label per point) --
        see :meth:`from_cells`.
        """

        geometries = list(getattr(points, "geometry", points))
        groups = cls._resolve_groups(group, len(geometries))
        rows = []
        for irpt, point in enumerate(geometries):
            x, y = (point.x, point.y) if hasattr(point, "x") else point[:2]
            cell = model.gwf.modelgrid.intersect(float(x), float(y))
            if isinstance(cell, tuple):
                cellid = (int(layer), *(int(value) for value in cell)) if len(cell) == 2 else tuple(int(value) for value in cell)
            else:
                cellid = (int(layer), int(cell))
            record: tuple[Any, ...] = (irpt, cellid, float(x), float(y), float(local_z))
            rows.append(record if groups is None else (*record, groups[irpt]))
        return cls(tuple(rows))


@dataclass
class PRTRunResults:
    """Read access to the pathlines/output of a finished PRT simulation.

    A file-backed handle to one completed PRT run: it points at the track CSV (and
    optional budget) in ``workspace`` and exposes the simulated particle pathlines
    for plotting and analysis, joined back to the originating ``flow_model``'s grid.
    Created for you by :class:`PRTProject` after a run, or reopened from disk with
    :func:`open_prt_run` without rebuilding the simulation.

    Results are read through view nouns, each answering the usual
    ``get``/``summary``/``plot``/``map``/``mosaic`` verbs
    (``myflopy.modflow.mf6.prt_maps``)::

        results.pathlines.map()                  # the tracks over the head map
        results.pathlines.get()                  # normalized track records
        results.travel_time.map(stat="median")   # time-of-travel choropleth
        results.endpoints.get()                  # counts per terminating cell
        results.capture.map()                    # one panel per release group

    The three per-cell maps are time-integrated over the whole run, so they take
    no ``per``. :attr:`pathlines` is a **view**, not a frame -- call ``.get()``
    for the records (or :attr:`track_records` for the untouched CSV).

    Attributes
    ----------
    flow_model
        The flow model the particles were tracked through.
    workspace
        Directory holding the PRT output files.
    name
        PRT model name.
    track_csv_path, budget_path
        Paths to the particle-track CSV and (optional) budget output.
    """

    flow_model: SimulationBase
    workspace: Path
    name: str
    track_csv_path: Path
    budget_path: Path | None = None
    success: bool | None = None
    report: tuple[str, ...] = ()

    @property
    def engine(self) -> str:
        """The particle-tracking engine that produced these results (``"mf6-prt"``)."""

        return "mf6-prt"

    @cached_property
    def track_records(self) -> pd.DataFrame:
        """Load and cache the PRT track CSV exactly as MF6 wrote it.

        The raw table: one-based ``icell``, no derived columns. FloPy and PyVista
        consume it directly; for analysis prefer ``pathlines.get()``, which adds
        zero-based ``cell``/``layer``, ``travel_time``, and the release group.
        """

        if not self.track_csv_path.exists():
            raise FileNotFoundError(f"PRT track CSV not found: {self.track_csv_path}")
        return pd.read_csv(self.track_csv_path)

    def refresh(self) -> PRTRunResults:
        """Clear cached file-backed results so subsequent access rereads disk."""

        self.__dict__.pop("track_records", None)
        return self

    @property
    def terminal_points(self) -> pd.DataFrame:
        """The terminating pathline records (``ireason == 3``), one per particle endpoint."""

        data = self.track_records
        if "ireason" not in data.columns:
            return data.iloc[0:0].copy()
        return data.loc[data["ireason"] == 3].copy()

    @staticmethod
    def _maps():
        """The derived-map views module (imported lazily -- it pulls in the plot layer)."""

        from myflopy.modflow.mf6 import prt_maps

        return prt_maps

    @property
    def pathlines(self):
        """Pathline view: the tracks themselves (``get``/``summary``/``plot``/``map``/``mosaic``).

        A **view**, not a DataFrame -- ``pathlines.get()`` returns the normalized
        records and ``pathlines.map()`` draws them over the model map.
        """

        return self._maps().PRTPathlineView(self)

    @property
    def travel_time(self):
        """Travel-time view: ``get``/``summary``/``plot``/``map``/``mosaic`` (time of travel per cell)."""

        return self._maps().PRTTravelTimeView(self)

    @property
    def endpoints(self):
        """Endpoint view: particle-termination counts per cell, with the same verbs."""

        return self._maps().PRTEndpointsView(self)

    @property
    def capture(self):
        """Capture-zone view: endpoints split by release group (``map()`` gives one panel per group)."""

        return self._maps().PRTCaptureView(self)

    def scene(self, **kwargs):
        """The 3-D pathline scene, as a :class:`~myflopy.viz.VtkScene` Picture.

        The same picture as ``flow_model.plot.grid(pathlines=..., backend="vtk")``
        -- this is the spelling you reach for when you already have the run
        results in hand. Display it, ``.show()`` it, or ``.html(path)`` it.
        """

        from myflopy.modflow.mf6.interactive_plotting import build_particle_tracking_scene

        return build_particle_tracking_scene(self.flow_model, self.track_records, **kwargs)

    def plot_map(self, **kwargs):
        """Plot the particle pathlines on a plan-view map of the flow model."""

        from myflopy.modflow.mf6.interactive_plotting import plot_particle_pathlines

        return plot_particle_pathlines(self.flow_model, self.track_records, **kwargs)

    def export_3d_html(self, output_path: str | Path, **kwargs) -> Path:
        """Write the 3-D pathline scene to a standalone interactive HTML file.

        Kept as a named convenience because "give me a file I can send" is a
        distinct intent; it is exactly ``self.scene(**kwargs).html(output_path)``,
        and closes the plotter afterwards.
        """

        scene = self.scene(**kwargs)
        try:
            return scene.html(output_path)
        finally:
            scene.scene.close()


class PRTProject:
    """Build, run, and collect an MF6 particle-tracking run off an existing flow model.

    The high-level driver for native MF6 PRT: given a completed/loaded flow
    ``model`` and a set of ``release_points``, it assembles the PRT model (MIP/PRP/OC
    with the supplied porosity, tracking-time, and termination options), couples it
    to the flow model, writes and runs it in ``workspace``, and returns a
    :class:`PRTRunResults` for the pathlines. Construct it through
    ``model.particle_tracking.prt(...)`` (the :class:`ParticleTracking` front door)
    or directly. For the declarative spec API, see ``mf.prt(...)``.

    Parameters
    ----------
    model
        The flow model to track particles through.
    workspace
        Directory to write/run the PRT simulation in.
    release_points
        A :class:`PRTReleasePoints` (or raw PRP records) defining where particles start.
    porosity
        Aquifer porosity used to convert flow to seepage velocity.
    perioddata, time_units, release_perioddata, extend_tracking, stoptime, stoptraveltime, track_times
        Tracking time-discretization and termination controls.
    drape, local_z, stop_at_weak_sink
        Particle placement/termination options.
    name, exe_name
        PRT model name (<=16 chars) and the MODFLOW executable.
    """

    def __init__(
        self,
        model: SimulationBase,
        *,
        workspace: str | Path,
        name: str | None = None,
        release_points: PRTReleasePoints | Sequence[tuple[Any, ...]],
        porosity: float | Sequence[float] = 0.2,
        perioddata: Sequence[tuple[float, int, float]] | None = None,
        time_units: str | None = None,
        release_perioddata: dict | None = None,
        extend_tracking: bool = False,
        stoptime: float | None = None,
        stoptraveltime: float | None = None,
        drape: bool = True,
        local_z: bool = True,
        stop_at_weak_sink: bool = False,
        track_times: Sequence[float] | None = None,
        exe_name: str = "mf6",
    ):
        """Configure a PRT run off ``model``: validate the name/release points and set output paths.

        See the class docstring for the tracking-time, placement, and termination
        parameters. Raises if the name exceeds 16 characters or no release points
        are given.
        """

        self.flow_model = model
        self.workspace = Path(workspace)
        self.name = _default_prt_name(model.name) if name is None else str(name)
        if len(self.name) > 16:
            raise ValueError("MF6 PRT model names cannot exceed 16 characters.")
        self.release_points = (
            release_points
            if isinstance(release_points, PRTReleasePoints)
            else PRTReleasePoints(tuple(tuple(row) for row in release_points))
        )
        if not self.release_points.packagedata:
            raise ValueError("At least one PRT release point is required.")
        self.workspace.mkdir(parents=True, exist_ok=True)
        self.track_csv_path = self.workspace / f"{self.name}.trk.csv"
        self.budget_path = self.workspace / f"{self.name}.bud"

        self.sim = flopy.mf6.MFSimulation(
            sim_name=self.name,
            version="mf6",
            exe_name=exe_name,
            sim_ws=self.workspace,
        )
        self.tdis = flopy.mf6.ModflowTdis(
            self.sim,
            time_units=time_units or _tdis_time_units(model),
            nper=len(perioddata or _tdis_perioddata(model)),
            perioddata=list(perioddata or _tdis_perioddata(model)),
        )
        self.prt = flopy.mf6.ModflowPrt(
            self.sim,
            modelname=self.name,
            model_nam_file=f"{self.name}.nam",
        )
        self.discretization = _copy_grid_to_prt(model, self.prt)
        self.mip = flopy.mf6.ModflowPrtmip(self.prt, pname="mip", porosity=porosity)
        self.prp = flopy.mf6.ModflowPrtprp(
            self.prt,
            pname="prp",
            # Release-group labels are PRP boundnames; without BOUNDNAMES the sixth
            # field is rejected, so the flag follows the record width.
            boundnames=self.release_points.has_groups or None,
            nreleasepts=len(self.release_points.packagedata),
            packagedata=list(self.release_points.packagedata),
            perioddata={0: ["FIRST"]} if release_perioddata is None else release_perioddata,
            extend_tracking=extend_tracking,
            stoptime=stoptime,
            stoptraveltime=stoptraveltime,
            drape=drape,
            local_z=local_z,
            stop_at_weak_sink=stop_at_weak_sink,
            coordinate_check_method=None,
        )
        oc_kwargs = {}
        if track_times is not None:
            oc_kwargs.update(
                ntracktimes=len(track_times),
                tracktimes=[(float(value),) for value in track_times],
            )
        self.oc = flopy.mf6.ModflowPrtoc(
            self.prt,
            pname="oc",
            budget_filerecord=[self.budget_path.name],
            trackcsv_filerecord=[self.track_csv_path.name],
            saverecord=[("BUDGET", "ALL")],
            **oc_kwargs,
        )
        flow_workspace = Path(model.model_output_folder_path)
        self.fmi = flopy.mf6.ModflowPrtfmi(
            self.prt,
            pname="fmi",
            packagedata=[
                ("GWFHEAD", Path(os.path.relpath(flow_workspace / _output_filename(model, "head"), self.workspace)).as_posix()),
                (
                    "GWFBUDGET",
                    Path(os.path.relpath(flow_workspace / _output_filename(model, "budget"), self.workspace)).as_posix(),
                ),
                (
                    "GWFGRID",
                    Path(os.path.relpath(flow_workspace / _grid_filename(model), self.workspace)).as_posix(),
                ),
            ],
        )
        self.ems = flopy.mf6.ModflowEms(self.sim, pname="ems", filename=f"{self.name}.ems")

    def write(self, *, silent: bool = True) -> Path:
        """Write the PRT simulation files to the workspace and return its path."""

        self.sim.write_simulation(silent=silent)
        return self.workspace

    def results(self, *, success: bool | None = None, report: Sequence[str] = ()) -> PRTRunResults:
        """Build a :class:`PRTRunResults` handle over this project's output files."""

        return PRTRunResults(
            flow_model=self.flow_model,
            workspace=self.workspace,
            name=self.name,
            track_csv_path=self.track_csv_path,
            budget_path=self.budget_path,
            success=success,
            report=tuple(report),
        )

    def run(self, *, write: bool = True, silent: bool = True, report: bool = True) -> PRTRunResults:
        """Write (optionally) and run the PRT simulation, returning its :class:`PRTRunResults`.

        Requires the GWF head/budget/grid outputs to already exist (raises with the
        missing paths otherwise); raises ``RuntimeError`` if the MF6 run fails.
        """

        missing = [
            path
            for path in (
                Path(self.flow_model.model_output_folder_path) / _output_filename(self.flow_model, "head"),
                Path(self.flow_model.model_output_folder_path) / _output_filename(self.flow_model, "budget"),
                Path(self.flow_model.model_output_folder_path) / _grid_filename(self.flow_model),
            )
            if not path.exists()
        ]
        if missing:
            formatted = "\n".join(f"- {path}" for path in missing)
            raise FileNotFoundError(
                "PRT requires completed GWF head and budget outputs before it can run:\n"
                f"{formatted}"
            )
        if write:
            self.write(silent=silent)
        success, output = self.sim.run_simulation(silent=silent, report=report)
        if not success:
            tail = "\n".join(output[-20:]) if output else "No MF6 report was returned."
            raise RuntimeError(f"PRT simulation failed in {self.workspace}\n{tail}")
        return self.results(success=success, report=output)


def open_prt_run(
    model: SimulationBase,
    workspace: str | Path,
    *,
    name: str | None = None,
) -> PRTRunResults:
    """Reopen a finished PRT run from disk as :class:`PRTRunResults`.

    The review counterpart to :class:`PRTProject`: point it at a ``workspace`` that
    already contains PRT output and it locates the track/budget files and returns a
    results handle -- without re-assembling or re-running the PRT simulation. Use it
    to revisit pathlines from an earlier run.

    Parameters
    ----------
    model
        The flow model the particles were tracked through.
    workspace
        Directory containing the completed PRT output.
    name
        PRT model name; inferred from ``workspace`` when ``None``.

    Returns
    -------
    PRTRunResults
        A handle to the existing run's pathline output.
    """

    workspace = Path(workspace)
    if name is None:
        csv_candidates = sorted(workspace.glob("*.trk.csv"))
        if len(csv_candidates) != 1:
            raise FileNotFoundError(
                f"Expected one *.trk.csv file in {workspace}, found {len(csv_candidates)}."
            )
        track_path = csv_candidates[0]
        name = track_path.name.removesuffix(".trk.csv")
    else:
        track_path = workspace / f"{name}.trk.csv"
    return PRTRunResults(
        flow_model=model,
        workspace=workspace,
        name=name,
        track_csv_path=track_path,
        budget_path=workspace / f"{name}.bud",
    )


class ParticleTracking:
    """A model's particle-tracking front door (native MF6 PRT and legacy MP3DU).

    The small accessor exposed as ``model.particle_tracking``: it bundles the
    particle-tracking workflows available for a flow ``model`` so you do not have to
    import the project classes yourself. :meth:`prt` builds and runs a native MF6
    PRT simulation (returning a :class:`PRTProject`), while the MP3DU helpers cover
    the older MODPATH-style workflow.

    Parameters
    ----------
    model
        The flow model these particle-tracking workflows operate on.
    """

    def __init__(self, model: SimulationBase):
        """Bind the particle-tracking front door to a flow ``model``."""

        self.model = model

    def prt(self, *, workspace: str | Path, release_points, **kwargs) -> PRTProject:
        """Build a native MF6 PRT project for this model (see :class:`PRTProject`)."""

        return PRTProject(
            self.model,
            workspace=workspace,
            release_points=release_points,
            **kwargs,
        )

    def open_prt(self, workspace: str | Path, *, name: str | None = None) -> PRTRunResults:
        """Reopen a completed PRT run in ``workspace`` (see :func:`open_prt_run`)."""

        return open_prt_run(self.model, workspace, name=name)

    def mp3du(self, particles, **kwargs):
        """Prepare a legacy MP3DU (MODPATH-style) particle-tracking run for this model."""

        from myflopy.modflow.mp3du import prepare_particle_tracking

        return prepare_particle_tracking(model=self.model, particles=particles, **kwargs)
