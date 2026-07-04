"""Core simulation object shared by live builds and file-backed loaded runs."""

from __future__ import annotations

from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any

import flopy
import numpy as np
import pandas as pd

from myflopy.modflow.mf6.simulation.accessors import (
    build_choro,
    build_xsection,
    get_all_heads,
    get_budget,
    get_budget_cumulative,
    get_budget_incremental,
    get_hds,
    get_inputs,
    get_kstpkper,
    get_lak_output,
    get_outputs,
    get_packages,
    get_sfr_output,
    get_surface,
    get_uzf_output,
)
from myflopy.modflow.mf6.simulation.indexing import (
    build_idomain,
    build_model_times,
    build_ncpl_arr,
    build_node_to_lni,
    build_offsets,
    coerce_per_dates,
)
from myflopy.modflow.mf6.surface_water_validation import (
    SurfaceWaterValidationReport,
    validate_surface_water_configuration,
)
from myflopy.modflow.mf6.simulation.runtime import run_simulation as run_model_simulation
from myflopy.modflow.mf6.simulation.regions import (
    RegionGroup,
    RegionRegistry,
    filter_region_heads,
    get_region_cells,
    list_model_groups,
    list_model_regions,
)

if TYPE_CHECKING:
    import numpy as np

    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from myflopy.modflow.mf6.sfr import SFRBuilder
    from myflopy.modflow.mf6.interactive_plotting import ModelVisualization
    from myflopy.modflow.mf6.observations import TargetRegistry
    from myflopy.modflow.mf6.package_explorer import ModelPackages
    from myflopy.modflow.mf6.parallel import ParallelModelWorkflow
    from myflopy.modflow.mf6.prt import ParticleTracking
    from myflopy.modflow.mf6.simulation.accessors import ModelOutputs


class SimulationBase:
    """The legacy stateful model object: build, run, and explore one MF6 model.

    The original imperative entry point that predates the declarative spec API
    (:class:`~myflopy.specs.SimulationSpec`). An instance owns a model's FloPy
    simulation, its Voronoi grid, and -- after a run -- convenient access to
    outputs (heads, budgets, lists, observation CSVs via ``_OUTPUT_SUFFIXES``)
    plus plotting/visualization helpers. The OO builder classes
    (``SFRBuilder``/``LAKBuilder``/...) and ``canonical_example.py`` are built on
    this object.

    For new model assembly prefer the package-first API (``mf.gwf(...)`` on a
    :class:`~myflopy.workspace.Project`); reach for ``SimulationBase`` when working
    with existing legacy code paths or loading a previously built run to explore
    its results.
    """

    _OUTPUT_SUFFIXES = {".cbc", ".hds", ".lst", ".bud", ".obs.csv", ".grb"}

    def _initialize_model_state(
        self,
        *,
        name: str,
        vor: Vor | None = None,
        nper: int | None = 1,
        nlay: int | None = None,
        model_output_folder_path: Path,
        per_dates: list | pd.DatetimeIndex | None = None,
        idomain_path: Path | None = None,
        idomain: Any = None,
        grid_type_override: str | None = None,
    ) -> None:
        """Initialize model metadata and lazy caches shared by all live views."""

        self.name = name
        self.vor = vor
        self._vor = vor
        self.nlay = nlay
        self.nper = nper
        self.num_steps = None
        self.per_len = None
        self._obs = None
        self._times = None
        self.model_output_folder_path = Path(model_output_folder_path)
        self._master_celld = {}
        self._hds = None
        self._idomain_path = None
        self._bud = None
        self._per_dates = None
        self._idomain_gdf = None
        self._inactive_cells = None
        self._node_to_lni = None
        self._lni_to_node = None
        self._kstpkper = None
        self._idomain = idomain
        self._grid_type_override = grid_type_override
        self.regions = RegionRegistry(self)
        self._targets = None
        self._visualize = None
        self._particle_tracking = None
        self._parallel = None

        self.per_dates = per_dates
        self.idomain_path = idomain_path

    @classmethod
    def from_built_run(cls, run: Any, model_name: str) -> "SimulationBase":
        """Create a myflopy model view over a model in a live built run."""

        view = cls.__new__(cls)
        view._initialize_from_built_run(run, model_name)
        return view

    def _initialize_from_built_run(self, run: Any, model_name: str) -> None:
        """Hydrate this object from a model in a live built run."""

        if run.built is None:
            raise ValueError("A live model view requires a built run.")

        built_model = run.built.built_model(model_name)
        gwf = built_model.as_gwf()
        context = built_model.context

        self.sim = run.built.simulation
        self.gwf = gwf
        self.ims = None
        self._initialize_model_state(
            name=model_name,
            vor=context.grid,
            nper=self._simulation_nper(self.sim),
            nlay=self._model_nlay(gwf),
            model_output_folder_path=run.workspace,
            idomain=context.domain,
        )
        self.source = "built_run"

    @staticmethod
    def _simulation_nper(simulation: Any) -> int | None:
        """Return the simulation stress-period count when TDIS is available."""

        try:
            return int(simulation.tdis.nper.data)
        except Exception:
            return None

    @staticmethod
    def _model_nlay(model: Any) -> int | None:
        """Return the model layer count when grid metadata is available."""

        try:
            return int(model.modelgrid.nlay)
        except Exception:
            try:
                return int(np.asarray(model.modelgrid.botm).reshape(-1).size)
            except Exception:
                return None

    def __init__(
        self,
        name: str = 'mf6_model',
        mf_folder_path: Path = Path().home().joinpath('mf6'),
        nper: int = 1,
        vor: Vor = None,
        per_dates: list | pd.DatetimeIndex = None,
        idomain_path: Path = None,
        newton: bool = True,
        complexity: str = 'MODERATE',
    ):
        """Create a legacy single-GWF model helper."""
        self._initialize_model_state(
            name=name,
            vor=vor,
            nper=nper,
            nlay=None,
            model_output_folder_path=mf_folder_path.joinpath(f'{name}'),
            per_dates=per_dates,
            idomain_path=idomain_path,
        )

        newtonoptions = 'under_relaxation' if newton else None

        self.sim = flopy.mf6.MFSimulation(
            sim_name=self.name,
            exe_name="mf6",
            version="mf6",
            sim_ws=self.model_output_folder_path,
        )
        self.gwf = flopy.mf6.ModflowGwf(
            self.sim,
            modelname=self.name,
            model_nam_file=f"{self.name}.nam",
            print_flows=False,
            newtonoptions=newtonoptions,
            save_flows=True,
        )

        self.ims = flopy.mf6.modflow.mfims.ModflowIms(
            self.sim,
            print_option="SUMMARY",
            csv_inner_output_filerecord="ims_inner_convergence.csv",
            csv_outer_output_filerecord="ims_outer_convergence.csv",
            pname="ims_gwf",
            complexity=complexity,
            under_relaxation="DBD",
            under_relaxation_theta=0.7,
            under_relaxation_kappa=0.20,
            under_relaxation_momentum=0.001,
            backtracking_number=50,
            backtracking_tolerance=1.1,
            backtracking_reduction_factor=0.2,
            backtracking_residual_limit=10,
            outer_maximum=800,
            inner_maximum=300,
            outer_dvclose=1e-5,
            inner_dvclose=1e-6,
            linear_acceleration="BICGSTAB",
            scaling_method="L2NORM",
            reordering_method="RCM",
            preconditioner_levels=5,
            preconditioner_drop_tolerance=0,
            relaxation_factor=0.97,
        )

        self.sim.register_ims_package(self.ims, [self.name])

    @property
    def workspace(self) -> Path:
        """Workspace folder where MF6 input/output files for this model live."""

        return self.model_output_folder_path

    @property
    def pest_runs(self) -> list:
        """The PEST runs done on this model (``<workspace>/pest``), ready to review.

        Returns a list of :class:`~myflopy.modflow.mf6.pest.runs.PestRunHandle`;
        call ``.review()`` on one to open it as :class:`IesResults` (phi, ensemble
        vs. observations, forecasts, parameter-field maps). Empty until a
        :class:`~myflopy.modflow.mf6.pest.PestProject` has been built/run with its
        default workspace beside this model.
        """

        from myflopy.modflow.mf6.pest.runs import find_pest_runs

        return find_pest_runs(self.workspace / "pest", model_name=self.name, model=self)

    def pest(self, name: str, *, start_datetime: str | None = None, **kwargs):
        """Start a PEST/PEST++ calibration of this model.

        This is the write-side companion to :attr:`pest_runs` (the read side):
        :attr:`pest_runs` discovers calibrations already done on the model, while
        ``pest()`` creates a new one. It returns a fresh
        :class:`~myflopy.modflow.mf6.pest.PestProject` already bound to this model,
        with its template workspace defaulting to ``<workspace>/pest/<name>`` so the
        run is auto-discoverable afterwards. The whole loop then lives on the model::

            cal = model.pest("calib", start_datetime="2020-01-01")
            cal.parameterize("k", style="pilotpoints", pp_space=8)
            cal.observe(head_targets)
            cal.build(); cal.run_ies()
            model.pest_runs[-1].review().plot_phi()   # no imports needed

        Parameters
        ----------
        name
            Calibration name; used as the ``.pst`` stem and the default
            ``<workspace>/pest/<name>`` template directory.
        start_datetime
            Simulation start date (e.g. ``"2020-01-01"``) that pyEMU uses to place
            time-varying parameters/observations on the time axis. When omitted it
            is read from the model's TDIS ``start_date_time``; pass it explicitly if
            the model has none.
        **kwargs
            Forwarded to :class:`~myflopy.modflow.mf6.pest.PestProject`
            (``workspace``, ``spatial_reference``, ``zero_based``, ``longnames``).

        Returns
        -------
        PestProject
            A calibration project bound to this model.
        """

        from myflopy.modflow.mf6.pest.project import PestProject

        if start_datetime is None:
            start_datetime = self._tdis_start_datetime()
            if start_datetime is None:
                raise ValueError(
                    "start_datetime could not be read from the model's TDIS; "
                    "pass start_datetime=... explicitly (e.g. '2020-01-01')."
                )
        return PestProject(self, name, start_datetime=start_datetime, **kwargs)

    def _tdis_start_datetime(self) -> str | None:
        """The TDIS ``start_date_time`` as a string, or ``None`` if unset."""

        tdis = getattr(self.sim, "tdis", None)
        sdt = getattr(tdis, "start_date_time", None)
        if sdt is None:
            return None
        getter = getattr(sdt, "get_data", None)
        value = getter() if callable(getter) else sdt
        return str(value) if value else None

    @property
    def package_names(self) -> list[str]:
        """Sorted list of package names currently attached to the groundwater model."""

        return sorted(self.gwf.get_package_list())

    @property
    def targets(self) -> "TargetRegistry":
        """Model-bound registry for reusable calibration target sets."""

        if self._targets is None:
            from myflopy.modflow.mf6.observations import TargetRegistry

            self._targets = TargetRegistry(self)
        return self._targets

    @property
    def visualize(self) -> "ModelVisualization":
        """Model-bound standalone visualization and export helpers."""

        if self._visualize is None:
            from myflopy.modflow.mf6.interactive_plotting import ModelVisualization

            self._visualize = ModelVisualization(self)
        return self._visualize

    @property
    def particle_tracking(self) -> "ParticleTracking":
        """Model-bound MF6 PRT and MP3DU workflow entry point."""

        if self._particle_tracking is None:
            from myflopy.modflow.mf6.prt import ParticleTracking

            self._particle_tracking = ParticleTracking(self)
        return self._particle_tracking

    @property
    def parallel(self) -> "ParallelModelWorkflow":
        """Model-bound unified splitting and parallel execution workflow."""

        if self._parallel is None:
            from myflopy.modflow.mf6.parallel import ParallelModelWorkflow

            self._parallel = ParallelModelWorkflow(self)
        return self._parallel

    @property
    def grid_type(self) -> str:
        """Model grid flavor, typically ``disv`` or ``disu``."""

        override = getattr(self, "_grid_type_override", None)
        if override is not None:
            return override
        if getattr(self.gwf, "disu", None) is not None:
            return "disu"
        if getattr(self.gwf, "disv", None) is not None:
            return "disv"
        return "unknown"

    @grid_type.setter
    def grid_type(self, value: str | None):
        """Override the inferred grid type when rebuilding from external data."""

        self._grid_type_override = value

    @property
    def inactive_cells(self):
        """Flattened cell indices whose idomain value is inactive."""

        if self._inactive_cells is None:
            idomain = np.asarray(self.idomain).reshape(-1)
            self._inactive_cells = [idx for idx, value in enumerate(idomain) if int(value) == 0]
        return self._inactive_cells

    @property
    def idomain(self):
        """Flattened idomain array, building it from the Voronoi view when needed."""

        if self._idomain is None:
            self._idomain, self._idomain_gdf = build_idomain(self.vor, self.idomain_path)
        return self._idomain

    @property
    def idomain_gdf(self):
        """GeoDataFrame version of the idomain selection, built lazily."""

        _ = self.idomain
        return self._idomain_gdf

    @property
    def idomain_path(self):
        """Optional GIS/raster path used to build idomain selections."""

        return self._idomain_path

    @idomain_path.setter
    def idomain_path(self, idomain_path):
        """Set the optional path used when building idomain from external data."""

        if idomain_path is not None:
            assert isinstance(idomain_path, Path), 'idomain_path must be a Path object'
            self._idomain_path = idomain_path

    @property
    def per_dates(self):
        """Optional stress-period date labels."""

        return self._per_dates

    @per_dates.setter
    def per_dates(self, per_dates):
        """Coerce and store stress-period date labels."""

        self._per_dates = coerce_per_dates(per_dates)

    @property
    def times(self):
        """Model times extracted from the MF6 outputs."""

        if self._times is None:
            self._times = build_model_times(self.gwf)
        return self._times

    @property
    def modelgrid(self) -> flopy.discretization.vertexgrid.VertexGrid:
        """Underlying FloPy modelgrid object."""

        modelgrid: flopy.discretization.vertexgrid.VertexGrid = self.gwf.modelgrid
        return modelgrid

    @cached_property
    def node_to_lni(self) -> dict[int, tuple[int, int]]:
        """Map flattened node ids to ``(layer, index_in_layer)`` tuples."""

        return build_node_to_lni(self._ncpl_arr, self._offsets)

    @cached_property
    def _ncpl_arr(self) -> np.ndarray:
        """Per-layer cell counts cached from the model grid."""

        return build_ncpl_arr(self.modelgrid)

    @cached_property
    def _offsets(self) -> np.ndarray:
        """Prefix offsets used to map between flattened nodes and layers."""

        return build_offsets(self._ncpl_arr)

    def node_from_lni(self, layer: int, idx_in_layer: int) -> int:
        """Convert ``(layer, index_in_layer)`` back to a flattened node id."""

        return int(self._offsets[layer] + idx_in_layer)

    @property
    def cellids(self):
        """Convenience list of all flattened cell ids in the model."""

        return list(self.node_to_lni.keys())

    @property
    def hds(self):
        """FloPy heads reader for the current model outputs."""

        return get_hds(self)

    @property
    def all_heads(self):
        """Tabular heads view across available results."""

        return get_all_heads(self)

    def to_xugrid(self, *, layers=None, times=None, name: str = "head", masked: bool = True):
        """Export simulated heads as an xugrid ``(time, layer, cell)`` object.

        Convenience that stacks this model's saved heads across all layers and
        output times into one :class:`xugrid.UgridDataArray` on the Voronoi mesh
        -- so you never assemble a ``(nlay, ncpl)`` array by hand. The result
        plugs into xarray slicing, native unstructured plotting, and UGRID-NetCDF
        sharing (QGIS/ParaView). Delegates to
        :meth:`~myflopy.modflow.mf6.headsplus.HeadsPlus.to_xugrid`; see it for the
        full parameter list.

        Returns
        -------
        xugrid.UgridDataArray
            Dims ``("time", "layer", <face_dim>)`` with ``kstp``/``kper`` coords.

        Notes
        -----
        Keep the full ``(time, layer, cell)`` shape for analysis/NetCDF, but
        reduce to a single field (``.isel(time=-1, layer=0)``) before
        ``.ugrid.plot()`` -- it draws one value per cell and does not facet over
        ``time``/``layer``.

        Examples
        --------
        >>> uda = model.to_xugrid()
        >>> uda.isel(time=-1, layer=0).ugrid.plot()       # water table map
        >>> uda.isel(time=-1).ugrid.to_netcdf("heads.nc")  # share all layers
        """

        return self.hds.to_xugrid(layers=layers, times=times, name=name, masked=masked)

    @property
    def srf(self):
        """Surface/top-of-model values used by plotting and summaries."""

        return get_surface(self)

    def cor(
        self,
        kstpkper: tuple = None,
        per: int = None,
        per_timestep: int | str = "last",
        layer: int = 0,
        type: str = 'hds',
        custom_hover: dict = None,
        custom_zs: list = None,
        zmin: float | int = None,
        zmax: float | int = None,
        zoom: int = 13,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        show_layer_elevs: bool = True,
        show_mounding: bool = False,
        hover_heads: bool = True,
        hover_ks: bool = False,
        locs=None,
        rch_scale=None,
        bgs=False,
        hillshade_path: Path = None,
        colorscale: str = None,
        logscale: bool = False,
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list[float] = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        **kwargs,
    ):
        """Build a choropleth-style spatial plot from heads or other model values."""

        return build_choro(
            self,
            kstpkper=kstpkper,
            per=per,
            per_timestep=per_timestep,
            layer=layer,
            type=type,
            custom_hover=custom_hover,
            custom_zs=custom_zs,
            zmin=zmin,
            zmax=zmax,
            zoom=zoom,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            show_layer_elevs=show_layer_elevs,
            show_mounding=show_mounding,
            hover_heads=hover_heads,
            hover_ks=hover_ks,
            locs=locs,
            rch_scale=rch_scale,
            bgs=bgs,
            hillshade_path=hillshade_path,
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
            **kwargs,
        )

    def xs(
        self,
        per: int = None,
        kstpkper: tuple = None,
        layer: int = 0,
        cells: int | list[int] = None,
        x_or_y: str = None,
        spacing: int = 10,
        num_points: int = 100,
        extrapolate_beyond_section_ends: bool = False,
        interpolate: bool = False,
        use_rbf: bool = False,
        show_model_top=True,
        show_model_btm=False,
        animation_kstpkpers=None,
    ):
        """Build a cross-section style plot/view through the current model."""

        return build_xsection(
            self,
            per=per,
            kstpkper=kstpkper,
            layer=layer,
            cells=cells,
            x_or_y=x_or_y,
            spacing=spacing,
            num_points=num_points,
            extrapolate_beyond_section_ends=extrapolate_beyond_section_ends,
            interpolate=interpolate,
            use_rbf=use_rbf,
            show_model_top=show_model_top,
            show_model_btm=show_model_btm,
            animation_kstpkpers=animation_kstpkpers,
        )

    @property
    def inputs(self):
        """Convenience view of input packages and selected input data."""

        return get_inputs(self)

    @property
    def outputs(self) -> "ModelOutputs":
        """Namespace of package-specific output helpers."""

        return get_outputs(self)

    @property
    def packages(self) -> "ModelPackages":
        """Preferred package exploration namespace for inputs and maps.

        This is the higher-level package inspection surface intended to grow
        over time. It complements the older ``inputs`` helper with a more
        consistent API such as ``model.packages.rch.inputs.get()`` and
        ``model.packages.rch.inputs.map()``.
        """

        return get_packages(self)

    @property
    def lak_output(self):
        """Lake output accessor for the current run."""

        return get_lak_output(self)

    @property
    def uzf_output(self):
        """UZF output accessor for the current run."""

        return get_uzf_output(self)

    @property
    def sfr_output(self):
        """SFR output accessor for the current run."""

        return get_sfr_output(self)

    def package(self, package_names: str | list[str]):
        """Return one or more MF6 package objects attached to this model.

        Parameters
        ----------
        package_names
            One package name such as ``"rch"`` or a list of package names.
        """

        if isinstance(package_names, (list, tuple, set)):
            return {str(name).lower(): self.package(str(name)) for name in package_names}

        package_name = str(package_names).lower()
        package = self.gwf.get_package(package_name)
        if package is None and getattr(self, "sim", None) is not None:
            package = self.sim.get_package(package_name)
        if package is None:
            raise AttributeError(f"Package {package_name!r} is not attached to model '{self.name}'.")
        return package

    def load_all(self):
        """Eagerly load all model data.

        For live ``SimulationBase`` instances this is a no-op and simply
        returns ``self``. File-backed loaded runs override this to trigger a
        full FloPy load.
        """

        return self

    def _get_budget_reader(self):
        """Return the low-level cell-budget reader for this model.

        Live ``SimulationBase`` instances use FloPy's model-bound budget
        reader. File-backed runs override this to use the binary ``.cbc`` file
        directly so budget-only workflows can stay lighter.
        """

        if self._bud is None:
            self._bud = self.gwf.output.budget()
        return self._bud

    def _get_budget_kstpkper(self) -> list[tuple[int, int]]:
        """Return ``(kstp, kper)`` pairs available in the budget reader."""

        return [tuple(int(value) for value in pair) for pair in self._get_budget_reader().get_kstpkper()]

    @property
    def lak(self):
        """Return the MF6 LAK package attached to this model."""

        return self.package("lak")

    @property
    def uzf(self):
        """Return the MF6 UZF package attached to this model."""

        return self.package("uzf")

    @property
    def sfr(self):
        """Return the MF6 SFR package attached to this model."""

        return self.package("sfr")

    @property
    def rch(self):
        """Return the MF6 RCH package attached to this model."""

        return self.package("rch")

    def bud(self, package: str = None):
        """Return budget information, optionally filtered to one package."""

        return get_budget(self, package)

    @property
    def budget_cumulative(self):
        """Cumulative budget table for the current run."""

        return get_budget_cumulative(self)

    @property
    def budget_incremental(self):
        """Incremental budget table for the current run."""

        return get_budget_incremental(self)

    @property
    def master_celld(self):
        """Legacy storage used by some workflows for per-cell bookkeeping."""

        return self._master_celld

    @property
    def kstpkper(self):
        """Available ``(kstp, kper)`` output indices for this run."""

        return get_kstpkper(self)

    def add_region(
        self,
        name: str,
        *,
        cellids,
        layer: int | list[int] | None = None,
        category: str = "selection",
        package: str | None = None,
        tags: list[str] | None = None,
        geometry=None,
        metadata: dict | None = None,
        overwrite: bool = False,
    ):
        """Register a named region from explicit cell ids.

        Parameters
        ----------
        name
            Region name.
        cellids
            Cell ids to include in the region.
        layer
            Optional layer restriction or layer list.
        category, package, tags, geometry, metadata
            Optional descriptive metadata stored with the region.
        overwrite
            Whether an existing region of the same name may be replaced.
        """

        return self.regions.add_from_cells(
            name,
            cellids,
            layer=layer,
            category=category,
            package=package,
            tags=tags,
            geometry=geometry,
            metadata=metadata,
            overwrite=overwrite,
        )

    def add_region_from_cells(self, name: str, cellids, **kwargs):
        """Alias for :meth:`add_region` with explicit cell ids."""

        return self.add_region(name, cellids=cellids, **kwargs)

    def add_region_from_geometry(
        self,
        name: str,
        geometry,
        *,
        predicate: str = "intersects",
        layer: int | list[int] | None = None,
        category: str = "selection",
        package: str | None = None,
        tags: list[str] | None = None,
        metadata: dict | None = None,
        overwrite: bool = False,
    ):
        """Register a named region by selecting cells from a geometry."""

        return self.regions.add_from_geometry(
            name,
            geometry,
            predicate=predicate,
            layer=layer,
            category=category,
            package=package,
            tags=tags,
            metadata=metadata,
            overwrite=overwrite,
        )

    def get_region(self, name: str):
        """Return one registered region object by name."""

        return self.regions.get_region(name)

    def add_group(
        self,
        name: str,
        *,
        members=None,
        tags: list[str] | None = None,
        metadata: dict | None = None,
        overwrite: bool = False,
    ) -> RegionGroup:
        """Create a named region group from region/group members."""

        return self.regions.add_group(
            name,
            members=members,
            tags=tags,
            metadata=metadata,
            overwrite=overwrite,
        )

    def add_to_group(self, group_name: str, member_name: str):
        """Add one region or group to an existing region group."""

        return self.regions.add_to_group(group_name, member_name)

    def remove_from_group(self, group_name: str, member_name: str):
        """Remove one member from an existing region group."""

        return self.regions.remove_from_group(group_name, member_name)

    def get_group(self, name: str):
        """Return one registered region group by name."""

        return self.regions.get_group(name)

    def get_region_cells(self, name: str) -> list[int]:
        """Return the resolved cell ids for a named region or group."""

        return get_region_cells(self, name)

    def list_regions(self):
        """List registered regions on the model."""

        return list_model_regions(self)

    def list_groups(self):
        """List registered region groups on the model."""

        return list_model_groups(self)

    def resolve_region_cells(self, name: str) -> list[int]:
        """Resolve a region/group to its unique flattened cell ids."""

        return self.regions.resolve_cells(name)

    def resolve_region_cells_with_trace(self, name: str):
        """Resolve a region/group and return overlap/contribution trace data."""

        return self.regions.resolve_cells_with_trace(name)

    def remove_region(self, name: str):
        """Remove a registered region or region group by name."""

        return self.regions.remove(name)

    def region_heads(
        self,
        name: str,
        *,
        per: int | None = None,
        kstpkper: tuple | None = None,
        layer: int | list[int] | None = None,
    ):
        """Return head results filtered to a named region or region group."""

        return filter_region_heads(
            self,
            name,
            per=per,
            kstpkper=kstpkper,
            layer=layer,
        )

    def validate_surface_water(
        self,
        *,
        nlakes: int | None = None,
        lak_packagedata: list | None = None,
        lak_connectiondata: list | None = None,
        lak_perioddata: dict | None = None,
        sfr: "SFRBuilder" | None = None,
        maxmvr: int | None = None,
        maxpackages: int | None = None,
        mvr_packages: list | None = None,
        mvr_perioddata: dict | None = None,
        raise_on_error: bool = False,
    ) -> SurfaceWaterValidationReport:
        """Validate coupled LAK/SFR/MVR inputs against the current model.

        Parameters
        ----------
        nlakes, lak_packagedata, lak_connectiondata, lak_perioddata
            Optional LAK package inputs to validate.
        sfr
            Optional SFR builder instance to validate.
        maxmvr, maxpackages, mvr_packages, mvr_perioddata
            Optional MVR package sizing and period-data inputs to validate.
        raise_on_error
            Whether to raise ``ValueError`` immediately when any error-level
            findings are present in the validation report.
        """

        report = validate_surface_water_configuration(
            self,
            nlakes=nlakes,
            lak_packagedata=lak_packagedata,
            lak_connectiondata=lak_connectiondata,
            lak_perioddata=lak_perioddata,
            sfr=sfr,
            maxmvr=maxmvr,
            maxpackages=maxpackages,
            mvr_packages=mvr_packages,
            mvr_perioddata=mvr_perioddata,
        )
        if raise_on_error:
            report.raise_for_errors()
        return report

    def list_input_files(self) -> list[Path]:
        """List input-like files currently present in the model workspace."""

        return sorted(
            path
            for path in self.workspace.iterdir()
            if path.is_file() and path.suffix.lower() not in self._OUTPUT_SUFFIXES
        )

    def list_output_files(self) -> list[Path]:
        """List output-like files currently present in the model workspace."""

        return sorted(
            path
            for path in self.workspace.iterdir()
            if path.is_file() and path.suffix.lower() in self._OUTPUT_SUFFIXES
        )

    def summary(self) -> pd.DataFrame:
        """Return a one-row high-level summary of the model/run."""

        return pd.DataFrame(
            [
                {
                    "name": self.name,
                    "workspace": str(self.workspace),
                    "grid_type": self.grid_type,
                    "ncpl": self._summary_ncpl(),
                    "nlay": self._summary_nlay(),
                    "nper": self._summary_nper(),
                    "packages": self._summary_packages(),
                    "source": self._summary_source(),
                }
            ]
        )

    def diff(
        self,
        other,
        *others,
        crs: str = "EPSG:2927",
        shared_grid: bool = False,
        verbosity_level: int = 0,
    ):
        """Diff this model against one or more others (this model is the reference).

        ``other`` / ``others`` may be model objects or workspace paths (each is
        loaded lazily), and a single list/tuple is accepted too -- so
        ``model.diff(other)``, ``model.diff(a, b)``, and ``model.diff([a, b])``
        all work. Returns a
        :class:`~myflopy.project.model_diff.ModelDiff` whose reference is this
        model; every other model is compared against it (reference-star). For
        explicit reference control or a reusable, named comparison, build a
        :class:`~myflopy.project.model_group.ModelGroup` and call its ``.diff()``.
        """

        from myflopy.project.model_group import ModelGroup

        extra: list = []
        for item in (other, *others):
            if isinstance(item, (list, tuple)):
                extra.extend(item)
            else:
                extra.append(item)
        group = ModelGroup(
            [self, *extra],
            reference=self.name,
            crs=crs,
            shared_grid=shared_grid,
            verbosity_level=verbosity_level,
        )
        return group.diff()

    def _summary_ncpl(self) -> int | None:
        """Return the cell count shown by :meth:`summary`.

        File-backed subclasses may override this to avoid triggering geometry or
        FloPy package loads when a lightweight summary is enough.
        """

        return int(np.asarray(self.modelgrid.xcellcenters, dtype=float).reshape(-1).size)

    def _summary_nlay(self) -> int | None:
        """Return the layer count shown by :meth:`summary`."""

        return self.nlay

    def _summary_nper(self) -> int | None:
        """Return the stress-period count shown by :meth:`summary`."""

        return self.nper

    def _summary_packages(self) -> list[str]:
        """Return the package list shown by :meth:`summary`."""

        return self.package_names

    def _summary_source(self) -> str:
        """Return the data source label shown by :meth:`summary`."""

        return getattr(self, "source", "simulation")

    def file_summary(self) -> pd.DataFrame:
        """Return a file-level summary of the current workspace contents."""

        rows = []
        for path in sorted(self.workspace.iterdir()):
            if not path.is_file():
                continue
            rows.append(
                {
                    "name": path.name,
                    "suffix": path.suffix.lower(),
                    "category": "output" if path.suffix.lower() in self._OUTPUT_SUFFIXES else "input",
                    "size_bytes": int(path.stat().st_size),
                    "modified_at": pd.Timestamp(path.stat().st_mtime, unit="s"),
                }
            )
        return pd.DataFrame(rows)

    def package_summary(self) -> pd.DataFrame:
        """Return a table describing currently attached FloPy package objects."""

        rows = []
        for package in self.gwf.packagelist:
            package_name = getattr(package, "package_name", None)
            if isinstance(package_name, (list, tuple)):
                package_name = ",".join(str(value) for value in package_name)
            package_label = package_name or getattr(package, "pname", None) or package.__class__.__name__
            rows.append(
                {
                    "package": str(package_label).upper(),
                    "package_type": getattr(package, "package_type", package.__class__.__name__),
                    "filename": getattr(package, "filename", None),
                }
            )
        return pd.DataFrame(rows)

    def output_summary(self) -> pd.DataFrame:
        """Return a lightweight summary of key MF6 output files for this run."""

        heads_path = self.workspace / f"{self.name}.hds"
        budget_path = self.workspace / f"{self.name}.cbc"
        listing_path = self.workspace / "mfsim.lst"

        kstpkper_count = None
        final_time = None
        if heads_path.exists():
            try:
                hds = flopy.utils.HeadFile(heads_path)
                kstpkpers = hds.get_kstpkper()
                kstpkper_count = len(kstpkpers)
                times = hds.get_times()
                if times:
                    final_time = float(times[-1])
            except Exception:
                pass

        return pd.DataFrame(
            [
                {
                    "name": self.name,
                    "heads_file": heads_path.name if heads_path.exists() else None,
                    "budget_file": budget_path.name if budget_path.exists() else None,
                    "listing_file": listing_path.name if listing_path.exists() else None,
                    "has_heads": heads_path.exists(),
                    "has_budget": budget_path.exists(),
                    "has_listing": listing_path.exists(),
                    "kstpkper_count": kstpkper_count,
                    "final_time": final_time,
                }
            ]
        )

    def grid_summary(self) -> pd.DataFrame:
        """Return a lightweight summary of the current model grid."""

        modelgrid = self.modelgrid
        xcenters = np.asarray(modelgrid.xcellcenters, dtype=float).reshape(-1)
        ycenters = np.asarray(modelgrid.ycellcenters, dtype=float).reshape(-1)
        ncpl_arr = np.asarray(self._ncpl_arr, dtype=int)

        minx = maxx = miny = maxy = None
        if getattr(modelgrid, "verts", None) is not None:
            verts = np.asarray(modelgrid.verts, dtype=float)
            if verts.size:
                minx = float(np.min(verts[:, 0]))
                maxx = float(np.max(verts[:, 0]))
                miny = float(np.min(verts[:, 1]))
                maxy = float(np.max(verts[:, 1]))

        return pd.DataFrame(
            [
                {
                    "name": self.name,
                    "grid_type": self.grid_type,
                    "nlay": int(modelgrid.nlay),
                    "ncpl": int(xcenters.size),
                    "node_count": int(np.sum(ncpl_arr)),
                    "minx": minx,
                    "maxx": maxx,
                    "miny": miny,
                    "maxy": maxy,
                    "crs": getattr(self.vor, "crs", None) if hasattr(self, "_vor") and self._vor is not None else getattr(self, "_crs", None),
                }
            ]
        )

    def result_summary(self) -> pd.DataFrame:
        """Return a simple statistical summary of the final heads output."""

        heads_path = self.workspace / f"{self.name}.hds"
        if not heads_path.exists():
            return pd.DataFrame(
                [
                    {
                        "name": self.name,
                        "has_heads": False,
                        "cell_count": None,
                        "head_min": None,
                        "head_max": None,
                        "head_mean": None,
                    }
                ]
            )

        hds = flopy.utils.HeadFile(heads_path)
        data = np.asarray(hds.get_data()).astype(float).reshape(-1)
        finite = data[np.isfinite(data)]
        return pd.DataFrame(
            [
                {
                    "name": self.name,
                    "has_heads": True,
                    "cell_count": int(data.size),
                    "head_min": float(np.min(finite)) if finite.size else None,
                    "head_max": float(np.max(finite)) if finite.size else None,
                    "head_mean": float(np.mean(finite)) if finite.size else None,
                }
            ]
        )

    def plot(self, *args, **kwargs):
        """Delegate to FloPy's simulation-level plot helper."""

        return self.sim.plot(*args, **kwargs)

    def run_simulation(self):
        """Write and run the current simulation."""

        return run_model_simulation(self)
