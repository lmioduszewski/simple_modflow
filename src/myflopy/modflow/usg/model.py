"""The in-memory form of an imported MODFLOW-USG model."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

from myflopy._logging import get_logger
from myflopy.modflow.usg._io import NameFile
from myflopy.modflow.usg.cln import ClnData
from myflopy.modflow.usg.export import ExportManifest, export_gis
from myflopy.modflow.usg.packages import (
    BasData,
    DisuData,
    EtsData,
    LpfData,
    PeriodRecords,
    RchData,
    SmsData,
)

if TYPE_CHECKING:
    from myflopy.specs import SimulationSpec

logger = get_logger(__name__)

__all__ = ["UsgModel"]


@dataclass(slots=True)
class UsgModel:
    """A MODFLOW-USG model read into memory, ready to convert or to re-grid.

    Two things live here. The first is the model as MODFLOW-USG wrote it --
    :attr:`disu`, :attr:`lpf` and friends, in node numbering. The second, and
    the more useful one if you are moving to a new grid, is the same content
    reshaped to ``(nlay, ncpl)`` and paired with real coordinates: :attr:`k`,
    :attr:`botm`, :attr:`idomain`, the boundary-condition frames, and the CLN
    features as polygons. Those survive a change of grid; node numbers do not.

    Build one with :func:`~myflopy.modflow.usg.reader.read_usg` rather than
    directly.
    """

    name_file: NameFile
    disu: DisuData
    bas: BasData
    lpf: LpfData
    sms: SmsData | None = None
    grid: Any = None
    boundaries: dict[str, PeriodRecords] = field(default_factory=dict)
    rch: RchData | None = None
    ets: EtsData | None = None
    hfb: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))
    cln: ClnData | None = None
    oc_words: tuple[str, ...] = ()
    length_units: str = "unknown"
    time_units: str = "unknown"
    start_date_time: str | None = None

    # -- shape -------------------------------------------------------------

    @property
    def nlay(self) -> int:
        """Number of layers."""

        return self.disu.nlay

    @property
    def ncpl(self) -> int:
        """Cells per layer."""

        return self.disu.ncpl

    @property
    def nper(self) -> int:
        """Stress periods the model runs."""

        return self.disu.nper

    @property
    def nodes(self) -> int:
        """Total groundwater nodes."""

        return self.disu.nodes

    def _layered(self, values: np.ndarray) -> np.ndarray:
        """Reshape a per-node array to ``(nlay, ncpl)``."""

        return np.asarray(values, dtype=float).reshape(self.nlay, self.ncpl)

    # -- geometry ----------------------------------------------------------

    @property
    def top(self) -> np.ndarray:
        """Top of layer 1, shape ``(ncpl,)`` -- MODFLOW 6 DISV's ``top``."""

        return self._layered(self.disu.top)[0]

    @property
    def botm(self) -> np.ndarray:
        """Layer bottoms, shape ``(nlay, ncpl)``."""

        return self._layered(self.disu.bot)

    @property
    def idomain(self) -> np.ndarray:
        """Active-cell mask from IBOUND, shape ``(nlay, ncpl)``.

        MODFLOW-USG's negative IBOUND means *constant head*, which MODFLOW 6
        expresses through the CHD package rather than through IDOMAIN, so a
        negative value maps to active (1) here and is reported separately.
        """

        return (self._layered(self.bas.ibound) != 0).astype(int)

    @property
    def strt(self) -> np.ndarray:
        """Starting heads, shape ``(nlay, ncpl)``."""

        return self._layered(self.bas.strt)

    @property
    def thickness(self) -> np.ndarray:
        """Layer thickness, shape ``(nlay, ncpl)``."""

        tops = np.vstack([self.top[None, :], self.botm[:-1]])
        return tops - self.botm

    # -- properties --------------------------------------------------------

    @property
    def k(self) -> np.ndarray:
        """Horizontal hydraulic conductivity, shape ``(nlay, ncpl)``."""

        return self._layered(self.lpf.hk)

    @property
    def k33(self) -> np.ndarray:
        """Vertical hydraulic conductivity, shape ``(nlay, ncpl)``.

        ``LAYVKA = 0`` means the LPF array already *is* a conductivity; any
        other value makes it a ratio of horizontal to vertical, which is
        resolved here so the result is always a conductivity.
        """

        vka = self._layered(self.lpf.vka)
        ratio = np.asarray(self.lpf.layvka) != 0
        if not ratio.any():
            return vka
        out = vka.copy()
        out[ratio] = self.k[ratio] / np.where(vka[ratio] == 0, np.nan, vka[ratio])
        return out

    @property
    def ss(self) -> np.ndarray:
        """Specific storage, shape ``(nlay, ncpl)``."""

        return self._layered(self.lpf.ss)

    @property
    def sy(self) -> np.ndarray:
        """Specific yield, shape ``(nlay, ncpl)``."""

        return self._layered(self.lpf.sy)

    @property
    def icelltype(self) -> np.ndarray:
        """MODFLOW 6 cell type per layer: 0 confined, 1 convertible."""

        return np.where(np.asarray(self.lpf.laytyp) != 0, 1, 0)

    @property
    def uppermost_active(self) -> np.ndarray:
        """Layer index of the highest active cell in each column, shape ``(ncpl,)``.

        MODFLOW-USG's ``RCH``/``ETS`` can address "the highest active node"
        (``NRCHOP``/``NETSOP`` = 3); MODFLOW 6 has no such option and needs an
        explicit cell. IBOUND is static for the whole run, so resolving it once
        here is exact rather than approximate. A wholly inactive column returns
        layer 0 and is excluded by the converter.
        """

        active = self.idomain != 0
        first = np.argmax(active, axis=0)
        return np.where(active.any(axis=0), first, 0)

    @property
    def has_active_column(self) -> np.ndarray:
        """Boolean mask of columns with at least one active cell, shape ``(ncpl,)``."""

        return (self.idomain != 0).any(axis=0)

    # -- node numbering ----------------------------------------------------

    def to_cellid(self, node: np.ndarray | int) -> np.ndarray:
        """Convert 1-based USG node numbers to zero-based ``(layer, cell)`` pairs.

        Nodes above :attr:`nodes` belong to the CLN grid and have no MODFLOW 6
        cell; they come back as ``(-1, -1)`` and must be filtered by the caller.
        """

        node = np.asarray(node, dtype=np.int64)
        layer = (node - 1) // self.ncpl
        cell = (node - 1) % self.ncpl
        outside = node > self.nodes
        layer = np.where(outside, -1, layer)
        cell = np.where(outside, -1, cell)
        return np.column_stack([layer, cell])

    # -- grid-independent views -------------------------------------------

    def surfaces(self, *, interpolate: bool = True, method: str = "linear") -> list[Any]:
        """Return the layer contacts as myflopy :class:`~myflopy.surfaces.Surface` objects.

        ``nlay + 1`` surfaces, top first: the model top, then each layer bottom.

        Parameters
        ----------
        interpolate
            Build each surface from the cell centres as scattered points, so it
            can be evaluated on **any** grid. That is what makes an imported
            model a starting point for a model on a different mesh -- which is
            usually the reason to import one. ``False`` builds from the raw
            per-cell arrays instead, which is cheaper and exact but only valid
            on *this* grid: ``Surface.from_array`` is a same-grid constructor.
        method
            Interpolation method passed to :meth:`Surface.from_points`
            (``"linear"``, ``"nearest"``, ``"cubic"``). Ignored when
            ``interpolate`` is ``False``.

        Returns
        -------
        list[Surface]
        """

        from myflopy.surfaces import Surface

        if self.grid is None:
            raise ValueError("surfaces() needs the grid; pass gsf= to read_usg()")
        stack = [np.asarray(a, dtype=float) for a in (self.top, *self.botm)]
        if not interpolate:
            return [Surface.from_array(a) for a in stack]

        centres = self.grid.gdf_vorPolys.geometry.centroid
        xs = centres.x.to_numpy()
        ys = centres.y.to_numpy()
        return [Surface.from_points(xs, ys, a, method=method) for a in stack]

    def attach_layers_to_grid(self):
        """Publish this model's layer elevations onto ``grid.gdf_topbtm``.

        A grid built from a ``.gsf`` carries geometry and nothing else, but the
        model it came from has top and bottom for every cell -- it could not have
        run otherwise. Publishing them is what lets the layer-elevation hover,
        the mounding colorscale and the surface-aware SFR/LAK builders work on an
        imported grid, all of which read this one frame.

        The layout is the project's: ``geometry``, then integer column ``0`` for
        the model top and ``1..nlay`` for the layer bottoms -- the same shape
        :meth:`myflopy.layers.LayerStack.attach` writes.

        Called automatically by :func:`~myflopy.modflow.usg.reader.read_usg`;
        exposed because a grid rebuilt or re-gridded afterwards needs it again.
        """

        import geopandas as gpd

        if self.grid is None:
            raise ValueError("attach_layers_to_grid() needs the grid; pass gsf= to read_usg()")

        columns = {0: np.asarray(self.top, dtype=float)}
        for layer in range(self.nlay):
            columns[layer + 1] = np.asarray(self.botm[layer], dtype=float)
        self.grid.gdf_topbtm = gpd.GeoDataFrame(
            {"geometry": self.grid.gdf_vorPolys.geometry, **columns},
            geometry="geometry",
            crs=self.grid.crs,
        )
        logger.debug("published %d layer surfaces onto the grid", self.nlay + 1)
        return self.grid

    def boundary_frame(self, ftype: str):
        """Return one boundary condition as a GeoDataFrame with real geometry.

        The geometry is the cell polygon, so the frame can be re-mapped onto any
        other grid by intersection -- unlike the node numbers it came from.
        Only the periods that define their own records appear; a reused period
        carries the previous period's values by definition.
        """

        import geopandas as gpd
        import pandas as pd

        if self.grid is None:
            raise ValueError("boundary_frame() needs the grid; pass gsf= to read_usg()")
        records = self.boundaries.get(ftype.upper())
        if records is None:
            raise KeyError(f"no {ftype.upper()} package in this model")

        columns = _BC_COLUMNS.get(records.ftype, ("value1", "value2"))
        cells = self.grid.gdf_vorPolys
        frames = []
        for period, block in sorted(records.periods.items()):
            if len(block) == 0:
                continue
            cellid = self.to_cellid(block[:, 0].astype(np.int64))
            data = {"period": period, "layer": cellid[:, 0], "cell": cellid[:, 1]}
            for i, column in enumerate(columns, start=1):
                if block.shape[1] > i:
                    data[column] = block[:, i]
            frame = pd.DataFrame(data)
            frame["geometry"] = cells.geometry.to_numpy()[frame["cell"].to_numpy()]
            frames.append(frame)
        if not frames:
            return gpd.GeoDataFrame(
                {"period": [], "layer": [], "cell": [], "geometry": []},
                crs=getattr(cells, "crs", None),
            )
        return gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs=getattr(cells, "crs", None))

    def export_gis(
        self,
        directory,
        *,
        crs: str | None = None,
        stream_lines: dict | None = None,
        clip_to=None,
        resolution: float | None = None,
        bed_thickness: float = 1.0,
        overwrite: bool = False,
    ) -> ExportManifest:
        """Write this model's inputs as grid-independent vector and raster files.

        A converted model is bound to the grid it was converted on; this is the
        form that survives changing it. See
        :func:`~myflopy.modflow.usg.export.export_gis` for the parameters and for
        what each written file means.

        Returns
        -------
        ExportManifest

        Examples
        --------
        >>> usg = mf.read_usg("flow.nam", gsf="flow.gsf", crs="EPSG:2926")
        >>> written = usg.export_gis("from_usg", crs="EPSG:2927")
        >>> print(written.describe())
        """

        return export_gis(
            self,
            directory,
            crs=crs,
            stream_lines=stream_lines,
            clip_to=clip_to,
            resolution=resolution,
            bed_thickness=bed_thickness,
            overwrite=overwrite,
        )

    def cln_polygons(self):
        """Return the CLN features dissolved to polygons on the groundwater grid."""

        from myflopy.modflow.usg.cln import cln_polygons

        if self.cln is None:
            raise ValueError("this model has no CLN package")
        if self.grid is None:
            raise ValueError("cln_polygons() needs the grid; pass gsf= to read_usg()")
        return cln_polygons(self.cln, self.grid)

    # -- reporting ---------------------------------------------------------

    def report(self) -> str:
        """Return a plain-text account of what converts, what does not, and why.

        The deferrals are the point. A conversion that silently dropped a
        package would leave a model that runs and is wrong; this makes each
        omission countable, and says which cells it touched.
        """

        from myflopy.modflow.usg.convert import conversion_notes

        return conversion_notes(self)

    def to_mf6(
        self,
        name: str = "usg",
        *,
        crs: str | None = None,
        newton: bool | None = None,
        under_relaxation: bool = True,
        start_date_time: str | None = None,
        save_budget: bool = True,
        complexity: str = "COMPLEX",
        include: tuple[str, ...] | None = None,
        fix_for_mf6: bool = False,
        local_origin: bool = True,
    ) -> SimulationSpec:
        """Convert to a MODFLOW 6 :class:`~myflopy.specs.SimulationSpec`.

        Parameters
        ----------
        name
            Model and simulation name.
        crs
            Coordinate reference system to stamp on the grid, e.g. ``"EPSG:2927"``.
        newton
            Force the Newton formulation on or off. ``None`` follows the SMS
            package: a non-Picard ``NONLINMETH`` becomes ``NEWTON``.
        under_relaxation
            Add ``UNDER_RELAXATION`` to the Newton options.
        start_date_time
            ISO datetime for TDIS; defaults to what the DISU period labels say.
        save_budget
            Write a cell-by-cell budget file as well as heads.
        complexity
            IMS complexity preset. Defaults to ``"COMPLEX"``: a converted USG
            model is a hard nonlinear problem, and ``"MODERATE"`` makes MODFLOW 6
            die with SIGFPE on this class of model rather than converge slowly.
        include
            Restrict the converted boundary packages, e.g. ``("chd", "drn")``.
            ``None`` converts everything that mapped.
        fix_for_mf6
            Apply the smallest edits that make MODFLOW 6 accept the model --
            today, raising a GHB/RIV head that sits below its cell bottom up to
            that bottom. MODFLOW-USG tolerates such a boundary and MODFLOW 6
            refuses to run at all. Off by default, because a boundary head is
            not something to change without being told; :meth:`validate` names
            every case first, and each edit is logged when it is made.
        local_origin
            Write the mesh relative to its own lower-left corner and declare
            that corner as the DISV ``xorigin``/``yorigin``. Leave this on:
            MODFLOW 6 derives conductances from the raw vertex coordinates, and
            on a State Plane grid (~1.34 million feet here) that arithmetic
            loses enough precision to produce a NaN budget while still reporting
            "Normal termination". The model stays georeferenced either way.

        Returns
        -------
        SimulationSpec
        """

        from myflopy.modflow.usg.convert import to_mf6

        return to_mf6(
            self,
            name=name,
            crs=crs,
            newton=newton,
            under_relaxation=under_relaxation,
            start_date_time=start_date_time,
            save_budget=save_budget,
            complexity=complexity,
            include=include,
            fix_for_mf6=fix_for_mf6,
            local_origin=local_origin,
        )

    def validate(self) -> list[Any]:
        """Return what MODFLOW 6 will object to, without converting or running.

        Each entry says how many records are affected, why MODFLOW 6 minds, and
        the one-line fix. Entries whose ``blocks_run`` is true stop the run.
        """

        from myflopy.modflow.usg.convert import mf6_findings

        return mf6_findings(self)

    def __repr__(self) -> str:
        """Short identity: shape, timing, and the packages that were read."""

        return (
            f"UsgModel({self.name_file.path.name!r}, "
            f"{self.nlay} layers x {self.ncpl} cells, {self.nper} periods, "
            f"packages={list(self.boundaries)})"
        )


#: Field names after the node number, per list boundary condition.
_BC_COLUMNS = {
    "CHD": ("shead", "ehead"),
    "DRN": ("elev", "cond"),
    "GHB": ("bhead", "cond"),
    "RIV": ("stage", "cond", "rbot"),
    "WEL": ("q",),
}
