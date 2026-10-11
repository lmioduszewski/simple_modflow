"""The in-memory form of an imported MODFLOW-USG model."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

from myflopy._logging import get_logger
from myflopy.modflow.mf6.grid.geometry import top_active_layer
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

    Attributes
    ----------
    name_file : NameFile
        The parsed name file: its path, its directory (every relative path
        resolves against it), and its package rows and unit-number map.
    disu : DisuData
        ``DISU`` as read: node and layer counts, per-node ``top``/``bot``/
        ``area``, the connectivity (``iac``/``ja``), the stress-period timing
        (``perlen``/``nstp``/``tsmult``/``steady``) and the unit codes.
    bas : BasData
        ``BAS6`` as read: per-node ``ibound`` and ``strt``, and ``hnoflo``.
    lpf : LpfData
        ``LPF`` as read: per-layer flags (``laytyp``, ``layvka``, ...) and
        per-node ``hk``/``vka``/``ss``/``sy``.
    sms : SmsData or None
        ``SMS`` solver settings; ``None`` when the name file has no SMS.
    grid : VoronoiGridPlus or None
        The plan grid built from the ``.gsf``, with the layer elevations
        published onto ``grid.gdf_topbtm``. ``None`` when no ``.gsf`` was found
        or given; everything needing geometry (:meth:`surfaces`,
        :meth:`boundary_frame`, :meth:`export_gis`, :meth:`to_mf6`, ...) then
        raises ``ValueError``.
    boundaries : dict[str, PeriodRecords]
        The list boundary conditions read, keyed by upper-case file type --
        ``"CHD"``, ``"DRN"``, ``"GHB"``, ``"RIV"``, ``"WEL"`` -- with records
        in 1-based node numbering, GWF and CLN nodes kept apart. Empty when read
        with ``read_boundaries=False``.
    rch : RchData or None
        ``RCH`` as read: one recharge array per period that defines one.
        ``None`` when absent or not read.
    ets : EtsData or None
        ``ETS`` (segmented evapotranspiration) as read. ``None`` when absent or
        not read.
    hfb : numpy.ndarray
        ``HFB6`` barriers, shape ``(n, 3)``: 1-based ``node1``, ``node2`` and
        ``hydchr``, as read, as floats. Empty, shape ``(0, 3)``, when the name
        file lists no ``HFB6``/``HFB`` package or boundaries were not read (an
        HFB file declaring zero barriers gives an empty array too).
    cln : ClnData or None
        The CLN (connected linear network) as read; it has no MODFLOW 6
        counterpart and is not converted. ``None`` when the model has no CLN.
    oc_words : tuple of str
        The distinct output-control phrases the ``OC`` file uses, e.g.
        ``"SAVE HEAD"``. Empty without an OC file.
    length_units : str
        DISU ``LENUNI`` as MODFLOW 6's word: ``"feet"``, ``"meters"``,
        ``"centimeters"``, or ``"unknown"``.
    time_units : str
        DISU ``ITMUNI`` as MODFLOW 6's word: ``"seconds"``, ``"minutes"``,
        ``"hours"``, ``"days"``, ``"years"``, or ``"unknown"``.
    start_date_time : str or None
        ``"YYYY-MM-DD HH:MM:SS"`` inferred from the first ``NPER`` ``m/d/yyyy``
        dates in the DISU file (written on its stress-period lines): the first
        is read as the end of period 1, so the start is that date less the
        first period's length. Trusted only when there are at least two and
        they strictly increase; ``None`` otherwise. :meth:`to_mf6` uses it when
        not given one.
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

        The rule is the shared one, ``idomain > 0``
        (:func:`~myflopy.modflow.mf6.grid.geometry.top_active_layer`). A NEGATIVE
        USG IBOUND is a constant head -- an ACTIVE cell, unlike MF6's pass-through
        -1 -- which is why :attr:`idomain` maps it to 1 before this reads it; the
        ``!= 0`` this used before answered identically on that 0/1 array (ledger
        319).
        """

        first = top_active_layer(self.idomain)
        return np.where(first < 0, 0, first)

    @property
    def has_active_column(self) -> np.ndarray:
        """Boolean mask of columns with at least one active cell, shape ``(ncpl,)``."""

        return (self.idomain != 0).any(axis=0)

    # -- node numbering ----------------------------------------------------

    def to_cellid(self, node: np.ndarray | int) -> np.ndarray:
        """Convert 1-based USG node numbers to zero-based ``(layer, cell)`` pairs.

        Nodes above :attr:`nodes` belong to the CLN grid and have no MODFLOW 6
        cell; they come back as ``(-1, -1)`` and must be filtered by the caller.

        Parameters
        ----------
        node : int or array_like of int
            One or more 1-based USG node numbers, converted to ``int64``. The
            layered numbering is assumed: layer ``(node - 1) // ncpl``, cell
            ``(node - 1) % ncpl``. Values below 1 are not checked.

        Returns
        -------
        numpy.ndarray of int64, shape (n, 2)
            One ``(layer, cell)`` row per node -- a 2-D array even for a single
            node.
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
        interpolate : bool, optional
            Default ``True``: build each surface from the cell-polygon centroids
            as scattered points, so it can be evaluated on **any** grid. That is
            what makes an imported model a starting point for a model on a
            different mesh -- which is usually the reason to import one.
            ``False`` builds from the raw per-cell arrays instead, which is
            cheaper and exact but only valid on *this* grid:
            ``Surface.from_array`` is a same-grid constructor.
        method : str, optional
            Interpolation method passed to :meth:`Surface.from_points`
            (``"linear"``, the default, ``"nearest"`` or ``"cubic"``). Ignored
            when ``interpolate`` is ``False``.

        Returns
        -------
        list[Surface]

        Raises
        ------
        ValueError
            When the model has no grid (it was read without a ``.gsf``).
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

        Parameters
        ----------
        ftype : str
            The package's file type, case-insensitive: one of the keys of
            :attr:`boundaries` -- ``"CHD"``, ``"DRN"``, ``"GHB"``, ``"RIV"`` or
            ``"WEL"``, whichever the model has.

        Returns
        -------
        geopandas.GeoDataFrame
            One row per groundwater record: ``period`` (zero-based), ``layer``
            and ``cell`` (zero-based), the package's value columns in USG's
            names (``shead``/``ehead``, ``elev``/``cond``, ``bhead``/``cond``,
            ``stage``/``cond``/``rbot``, or ``q``) and the cell polygon as
            ``geometry``, in the CRS of ``grid.gdf_vorPolys``. Records on CLN
            nodes are not in it.
            A package with no groundwater records gives an empty frame with the
            ``period``/``layer``/``cell``/``geometry`` columns only.

        Raises
        ------
        ValueError
            When the model has no grid (it was read without a ``.gsf``).
        KeyError
            When the model has no package of type ``ftype``.
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
        form that survives changing it. Delegates to
        :func:`~myflopy.modflow.usg.export.export_gis`, whose module docstring
        explains what each written file means.

        Parameters
        ----------
        directory : str or pathlib.Path
            Where to write; created, with an ``arrays/`` subdirectory, if absent.
        crs : str, optional
            Coordinate reference system every vector and raster is written in;
            the data are reprojected into it from the grid's. Default ``None``
            uses the CRS of ``grid.gdf_vorPolys``. A ``.gsf`` carries no CRS, so
            check the grid's before relying on it -- the wrong one is silent and
            puts the export nowhere near the model it is meant to rebuild.
        stream_lines : dict of str to (str or pathlib.Path), optional
            ``{CLN stream name: vector file}``, mapping a CLN stream onto a
            centerline you would rather use. Each file is read, reprojected into
            the grid's CRS (assumed to be in it when it declares none), merged,
            and its longest piece kept; each stream node's ``station`` is then
            measured along *your* line, so interpolating onto your own reaches
            needs no reprojection. A name matching no CLN stream is ignored. The
            CLN's own centerline and ``station_cln`` are written either way.
            Default ``None`` measures along the CLN's own line.
        clip_to : str or pathlib.Path, optional
            A polygon file that each ``stream_lines`` file is clipped to before
            use (reprojected into the grid's CRS, or assumed to be in it when it
            declares none). A line empty after clipping is ignored with a
            warning. Ignored when ``stream_lines`` is not given.
        resolution : float, optional
            Pixel size of the per-cell array rasters under ``arrays/``, in the
            units of the written ``crs`` (the cells are reprojected before the
            raster is laid out). Default ``None`` takes a quarter of the median
            cell width, computed in the grid's own units, which resolves a
            typical cell and absorbs the smallest ones. The lakebed raster
            ``arrays/lakebed_bottom.tif`` is always written at that default,
            whatever this says.
        bed_thickness : float, optional
            Lakebed thickness, in model length units, assumed when deriving each
            lake's ``bedleak = fskin / bed_thickness`` in ``lakes.gpkg`` (both
            the per-lake and per-node layers); default ``1.0``. CLN has no bed
            thickness, so this is a choice, and ``fskin`` is written unmodified
            beside it. It does not reach ``streams.gpkg``, whose ``rbth`` column
            is always ``1.0``.
        overwrite : bool, optional
            Default ``False``: a ``directory`` that exists and is not empty is
            refused. ``True`` writes into it anyway; files of the same name are
            written over and nothing else there is removed.

        Returns
        -------
        ExportManifest
            What was written, file by file; the same content is saved as
            ``manifest.json`` beside a ``README.md``.

        Raises
        ------
        ValueError
            When the model has no grid, or when ``crs`` is not given and
            ``grid.gdf_vorPolys`` has no CRS.
        FileExistsError
            When ``directory`` is not empty and ``overwrite`` is ``False``.

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
        name : str, optional
            Name of the GWF model and of the simulation (default ``"usg"``). It
            also names the head file ``<name>.hds``, the budget file
            ``<name>.cbc`` and the model the IMS solution claims.
        crs : str, optional
            Coordinate reference system to stamp on the grid, e.g.
            ``"EPSG:2927"`` -- a RELABEL, saying the ``.gsf`` coordinates already
            are in that CRS; nothing is reprojected. It is assigned to
            ``self.grid.crs`` (the same grid object the returned spec's
            ``ModelContext`` holds, so this model's grid changes too), and the
            grid's setter carries it to every frame the grid caches:
            ``grid.gdf_vorPolys`` and ``grid.gdf_topbtm`` take the new label, and
            the lat/lon frames are rebuilt from it, so the grid reports ONE CRS
            (ledger 368). Default ``None`` leaves the grid's CRS as
            :func:`read_usg` set it.
        newton : bool, optional
            ``True``/``False`` forces the Newton formulation on or off in the
            model's ``NEWTONOPTIONS``. Default ``None`` follows the SMS package:
            a non-Picard ``NONLINMETH`` means Newton, and a model with no SMS
            gets none. The IMS under-relaxation and backtracking carried over
            from SMS follow SMS's own ``NONLINMETH``, not this.
        under_relaxation : bool, optional
            Default ``True`` writes ``NEWTON UNDER_RELAXATION``; ``False`` writes
            plain ``NEWTON``. Ignored when Newton is off.
        start_date_time : str, optional
            TDIS ``START_DATE_TIME``, e.g. ``"2020-01-01 00:00:00"``. Default
            ``None`` (or an empty string) uses :attr:`start_date_time`, which
            :func:`read_usg` infers from the DISU period labels and leaves
            ``None`` when they do not increase.
        save_budget : bool, optional
            Default ``True`` saves the cell-by-cell budget every time step to
            ``<name>.cbc`` as well as the heads to ``<name>.hds``. ``False``
            saves heads only.
        complexity : str, optional
            IMS complexity preset, ``"SIMPLE"``/``"MODERATE"``/``"COMPLEX"``,
            passed to :func:`myflopy.ims`. Defaults to ``"COMPLEX"``: a
            converted USG model is a hard nonlinear problem, and ``"MODERATE"``
            makes MODFLOW 6 die with SIGFPE on this class of model rather than
            converge slowly. SMS's outer closure, iteration caps and linear
            acceleration (and, for a Newton SMS that sets them, its
            delta-bar-delta under-relaxation and backtracking) are passed
            beside the preset and override its values.
        include : iterable of str, optional
            Restrict the optional packages to these names, case-insensitive:
            any of ``"chd"``, ``"drn"``, ``"ghb"``, ``"riv"``, ``"wel"``,
            ``"rch"``, ``"evt"`` and ``"hfb"``, e.g. ``("chd", "drn")``. A name
            outside that list is ignored without a message, and so a bare string
            (``"chd"``, read letter by letter) keeps none of them. DISV, NPF, IC, STO
            and OC are always written. Default ``None`` converts everything that
            mapped.
        fix_for_mf6 : bool, optional
            Apply the smallest edits that make MODFLOW 6 accept the model. Per
            period, a GHB head, a DRN invert or a RIV **bed** (``rbot``) below a
            floor just above its cell bottom -- the bottom plus
            ``max(|bottom| * 1e-6, 1e-6)``, so the written text cannot round
            below the bottom -- is raised to that floor. A river is lifted as a PAIR:
            ``rbot`` goes to the floor, and ``stage`` rises by the same delta
            only where the stage was below the floor too, because MODFLOW 6
            checks ``rbot`` against the cell bottom and ``stage >= rbot``
            separately, while a stage already inside its cell is left where it
            was drawn (ledger 291). A CHD record whose head (``ehead``) is not
            above that floor is omitted for that period, in confined and
            convertible layers alike; a period left with no CHD record is
            written EMPTY, so no constant head acts in it, rather than left out,
            which MODFLOW 6 would read as reusing the previous period's heads
            (ledger 367). Each kind of edit is counted in one logged
            warning per package.

            MODFLOW-USG tolerates such a boundary. Whether MODFLOW 6 refuses one
            depends on the layer -- see :meth:`validate` -- but the edits are
            applied wherever the elevation is outside its cell, not only where
            MODFLOW 6 would object, so the option means the same thing on every
            model. Default ``False``, because a boundary elevation is not
            something to change without being told: the model is converted as
            read, and every finding :meth:`validate` makes is logged as a
            warning instead.

            What it does NOT touch: a river whose stage is below its own bed.
            MODFLOW 6 refuses that at any ``icelltype``, and which of the two
            numbers is the wrong one is not something a converter can know.
        local_origin : bool, optional
            Default ``True`` writes the mesh relative to its own lower-left
            corner (the smallest vertex x and y) and declares that corner as the
            DISV ``xorigin``/``yorigin``. Leave this on: MODFLOW 6 derives
            conductances from the raw vertex coordinates, and on a State Plane
            grid (~1.34 million feet here) that arithmetic loses enough precision
            to produce a NaN budget while still reporting "Normal termination".
            The model stays georeferenced either way. ``False`` writes the raw
            coordinates with an origin of ``(0, 0)``.

        Returns
        -------
        SimulationSpec
            One GWF model on DISV with its TDIS (the DISU's period lengths,
            step counts and multipliers) and IMS. A list boundary's period that
            reused the previous one (``ITMP < 0``) is absent from its period
            dict, which MODFLOW 6 reads as reuse; a period that defined no
            groundwater record (``ITMP = 0``, or every record on a CLN node) is
            an EMPTY list, which switches the package off for that period as
            MODFLOW-USG did (ledger 367). A package with no groundwater
            record in any period is not written.

        Raises
        ------
        ValueError
            When the model has no grid (it was read without a ``.gsf``).
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

        **An entry whose ``blocks_run`` is false is not a harmless one.** For a
        boundary elevation below its own cell bottom, what MODFLOW 6 does splits
        on ``icelltype`` -- which a converted model inherits from the LPF
        ``LAYTYP``. A convertible layer is refused by name; a confined layer runs
        and over-discharges against a head below the cell, saying nothing in the
        listing file. Both are reported, with the marker saying which (ledger
        291). MODFLOW 6 also stops at the FIRST package it objects to, which is
        the argument for reading this list rather than running and seeing.
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
