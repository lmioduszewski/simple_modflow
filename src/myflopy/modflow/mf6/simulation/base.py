"""Core simulation object shared by live builds and file-backed loaded runs."""

from __future__ import annotations

import os
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any

import flopy
import numpy as np
import pandas as pd
from flopy.mf6.mfbase import MFDataException

from myflopy._logging import get_logger
from myflopy.areas import Area, live_areas, saved_areas
from myflopy.library import AoiLibrary
from myflopy.modflow.mf6.lak_implicit import warn_if_implicit_unsupported
from myflopy.modflow.mf6.package_registry import (
    _package_type_name,
    package_name_of,
    packages_of_type,
)
from myflopy.modflow.mf6.period_calendar import (
    PeriodCalendar,
    period_calendar,
    tdis_start_date_time,
)
from myflopy.modflow.mf6.simulation.accessors import (
    field_reader,
    get_all_conc,
    get_all_heads,
    get_all_temp,
    get_budget,
    get_budget_cumulative,
    get_budget_incremental,
    get_conc,
    get_hds,
    get_inputs,
    get_kstpkper,
    get_lak_output,
    get_model_budget,
    get_outputs,
    get_packages,
    get_sfr_output,
    get_temp,
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
from myflopy.modflow.mf6.simulation.regions import (
    RegionGroup,
    RegionRegistry,
    filter_region_heads,
    get_region_cells,
    list_model_groups,
    list_model_regions,
)
from myflopy.modflow.mf6.simulation.runtime import run_simulation as run_model_simulation
from myflopy.modflow.mf6.surface_water_validation import (
    SurfaceWaterValidationReport,
    validate_surface_water_configuration,
)
from myflopy.modflow.utils.datatypes.readers import as_path
from myflopy.plot import ModelPlots

if TYPE_CHECKING:
    import numpy as np

    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from myflopy.modflow.mf6.observations import TargetRegistry
    from myflopy.modflow.mf6.package_explorer import ModelPackages
    from myflopy.modflow.mf6.package_model import (
        HfbPackageExplorer,
        NpfPackageExplorer,
        PackageExplorer,
        StaticArrayPackageExplorer,
        UzfPackageExplorer,
    )
    from myflopy.modflow.mf6.package_surface_water import (
        LakPackageExplorer,
        SfrPackageExplorer,
        SurfaceWaterPackageExplorer,
    )
    from myflopy.modflow.mf6.parallel import ParallelModelWorkflow
    from myflopy.modflow.mf6.prt import ParticleTracking
    from myflopy.modflow.mf6.sfr import SFRBuilder
    from myflopy.modflow.mf6.simulation.accessors import ModelOutputs

logger = get_logger(__name__)


def _model_config():
    """:class:`~myflopy.project.model_config.ModelConfig`, imported at CALL time.

    The ONE deferred import of it in this module, shared by ``config`` and
    ``model_options``. Hoisting it to module level looks safe -- ``model_config``
    imports nothing from myflopy but ``_logging`` -- and is not: a bare
    ``import myflopy`` succeeds while the test suite's import order hits a circular
    import through ``myflopy.project.__init__``. Measured, not assumed. Two call sites
    would be two entries in the deferred-import ratchet, and a deferral that cannot be
    hoisted is worth stating once rather than counting twice.
    """

    from myflopy.project.model_config import ModelConfig

    return ModelConfig


def _type_then_name(typed: list, owner, key: str):
    """``model.package``'s answer on one FloPy model or simulation, or ``None``.

    ``typed`` is the owner's packages of TYPE ``key``: the one, or a list of several
    (FloPy's convention). With none, the owner's package NAMED ``key`` exactly (any
    case) -- never FloPy's filename or name-prefix matches (ledger 352).
    """

    if typed:
        return typed[0] if len(typed) == 1 else list(typed)
    named = [
        package
        for package in (owner.get_package() or ())
        if package_name_of(package).strip().lower() == key
    ]
    return named[0] if named else None



class FlopyPackages:
    """The model's FloPy package objects: ``model.fl.drn``, ``model.fl.npf``, ...

    ``model.drn`` is myflopy's view of a package (``model.drn.inputs.map()``,
    ``model.drn.results.q.get()``); ``model.fl.drn`` is the FloPy object the
    model was built with, for its data (``model.fl.drn.stress_period_data``).

    The ten registry packages are typed properties so an editor completes them,
    each found by TYPE alone: the model's one package of the type, or a list of
    every one, in model order, when it carries several -- whatever they are named
    (ledger 352: FloPy's own lookup answered a package NAMED ``lak`` before the
    type, so beside ``lak_2`` it returned ``lak`` alone). A model with none raises
    ``AttributeError`` naming it. Any other name -- ``npf``, ``ic``, ``sto``,
    ``disv``, ``oc``, ``tdis``, ``ims``, or a package's own name such as ``lak_2``
    -- resolves through :meth:`SimulationBase.package`: a type first, then a name.

    Parameters
    ----------
    model : SimulationBase
        The model whose packages are looked up, on every access (nothing is
        cached, so a package added or replaced later is seen).
    """

    def __init__(self, model: SimulationBase) -> None:
        self._model = model

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        if isinstance(getattr(type(self), name, None), property):
            # A typed property raised AttributeError -- the model has none of its type --
            # and Python retries the name HERE. The name lookup below would answer a
            # package merely NAMED `lak` for `model.fl.lak`; re-raise the type's refusal.
            return self._model._flopy_of_type(name)
        return self._model.package(name)

    def __dir__(self) -> list[str]:
        model = self._model
        return sorted(set(super().__dir__()) | set(model.package_types) | set(model.package_names))

    def __repr__(self) -> str:
        return f"<FloPy packages of {self._model.name!r}: {', '.join(self._model.package_types)}>"

    @property
    def chd(self) -> flopy.mf6.ModflowGwfchd | list[flopy.mf6.ModflowGwfchd]:
        """FloPy's CHD package -- a list when the model carries several."""

        return self._model._flopy_of_type("chd")

    @property
    def drn(self) -> flopy.mf6.ModflowGwfdrn | list[flopy.mf6.ModflowGwfdrn]:
        """FloPy's DRN package -- a list when the model carries several."""

        return self._model._flopy_of_type("drn")

    @property
    def ghb(self) -> flopy.mf6.ModflowGwfghb | list[flopy.mf6.ModflowGwfghb]:
        """FloPy's GHB package -- a list when the model carries several."""

        return self._model._flopy_of_type("ghb")

    @property
    def riv(self) -> flopy.mf6.ModflowGwfriv | list[flopy.mf6.ModflowGwfriv]:
        """FloPy's RIV package -- a list when the model carries several."""

        return self._model._flopy_of_type("riv")

    @property
    def wel(self) -> flopy.mf6.ModflowGwfwel | list[flopy.mf6.ModflowGwfwel]:
        """FloPy's WEL package -- a list when the model carries several."""

        return self._model._flopy_of_type("wel")

    @property
    def rch(self) -> flopy.mf6.ModflowGwfrch | list[flopy.mf6.ModflowGwfrch]:
        """FloPy's RCH package -- a list when the model carries several."""

        return self._model._flopy_of_type("rch")

    @property
    def evt(self) -> flopy.mf6.ModflowGwfevt | list[flopy.mf6.ModflowGwfevt]:
        """FloPy's EVT package -- a list when the model carries several."""

        return self._model._flopy_of_type("evt")

    @property
    def uzf(self) -> flopy.mf6.ModflowGwfuzf | list[flopy.mf6.ModflowGwfuzf]:
        """FloPy's UZF package -- a list when the model carries several."""

        return self._model._flopy_of_type("uzf")

    @property
    def sfr(self) -> flopy.mf6.ModflowGwfsfr | list[flopy.mf6.ModflowGwfsfr]:
        """FloPy's SFR package -- a list when the model carries several."""

        return self._model._flopy_of_type("sfr")

    @property
    def lak(self) -> flopy.mf6.ModflowGwflak | list[flopy.mf6.ModflowGwflak]:
        """FloPy's LAK package -- a list when the model carries several."""

        return self._model._flopy_of_type("lak")


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

    _OUTPUT_SUFFIXES = {".cbc", ".hds", ".ucn", ".lst", ".bud", ".obs.csv", ".grb"}

    #: MF6 model kind ("gwf6"/"gwt6"/"gwe6"); the legacy builder path is GWF, and
    #: ``_initialize_from_built_run`` overrides it from the built flopy model. It
    #: gates the dependent-variable readers (``hds``/``conc``/``temp``).
    model_type = "gwf6"

    #: The root of the project this model's run belongs to, recorded by
    #: ``Run.model``; ``None`` for a model outside any project.
    _project_root = None
    #: This model's areas when no LIVE project holds them (see :attr:`aoi`).
    _session_aoi = None

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
        self._particle_tracking = None
        self._parallel = None

        self.per_dates = per_dates
        self.idomain_path = idomain_path

    @classmethod
    def from_built_run(cls, run: Any, model_name: str) -> SimulationBase:
        """Create a myflopy model view over a model in a live built run.

        Skips ``__init__`` (no new FloPy simulation is created): the view shares
        the run's in-memory ``MFSimulation`` and the built FloPy model, takes its
        grid and idomain from that model's ``ModelContext``, its model kind from the
        FloPy model's ``model_type``, and its workspace from ``run.workspace``.

        Parameters
        ----------
        run : Run
            A :class:`~myflopy.workspace.Run` (anything with ``built`` and
            ``workspace`` attributes) that has been built; ``run.built`` being
            ``None`` raises ``ValueError``.
        model_name : str
            Name of a model in the built simulation; a name it does not hold
            raises ``KeyError``.

        Returns
        -------
        SimulationBase
            A view of ``cls`` bound to that model.

        Raises
        ------
        ValueError
            If ``run.built`` is ``None``.
        KeyError
            If the built simulation has no model named ``model_name``.
        """

        view = cls.__new__(cls)
        view._initialize_from_built_run(run, model_name)
        return view

    def _initialize_from_built_run(self, run: Any, model_name: str) -> None:
        """Hydrate this object from a model in a live built run."""

        if run.built is None:
            raise ValueError("A live model view requires a built run.")

        built_model = run.built.built_model(model_name)
        # Kind-neutral flopy model (GWF/GWT/GWE): ``.gwf`` stays the pervasive
        # handle everyone reads (modelgrid/output/get_package are kind-agnostic),
        # while ``.model_type`` gates the field readers (hds/conc/temp).
        flopy_model = built_model.model
        context = built_model.context

        self.sim = run.built.simulation
        self.gwf = flopy_model
        self.model_type = getattr(flopy_model, "model_type", "gwf6")
        self.ims = None
        self._initialize_model_state(
            name=model_name,
            vor=context.grid,
            nper=self._simulation_nper(self.sim),
            nlay=self._model_nlay(flopy_model),
            model_output_folder_path=run.workspace,
            idomain=context.domain,
        )
        self.source = "built_run"

    @staticmethod
    def _simulation_nper(simulation: Any) -> int | None:
        """Return the simulation stress-period count when TDIS is available."""

        try:
            return int(simulation.tdis.nper.data)
        except (AttributeError, TypeError, MFDataException):
            # AttributeError: no TDIS attached yet. TypeError: `nper.data` is
            # None on a TDIS declared but not populated. MFDataException: flopy
            # cannot read the block.
            logger.debug("no TDIS nper available", exc_info=True)
            return None

    @staticmethod
    def _model_nlay(model: Any) -> int | None:
        """Return the model layer count when grid metadata is available."""

        try:
            return int(model.modelgrid.nlay)
        except (AttributeError, TypeError, NotImplementedError, MFDataException):
            # NotImplementedError is the interesting one: flopy's BASE Grid
            # raises it from `nlay` when a model has no dis/disv/disu attached,
            # which is exactly the case this fallback is here for.
            logger.debug("modelgrid has no nlay; deriving it from botm")

        # `botm` is (nlay, ncpl) on DISV and (nlay, nrow, ncol) on DIS, so the
        # LAYER COUNT is its first axis. This used to read `.reshape(-1).size`,
        # which is the total number of botm ENTRIES -- on a 3-layer, 100-cell
        # grid it reported nlay as 300.
        try:
            botm = np.asarray(model.modelgrid.botm)
        except (AttributeError, ValueError, MFDataException):
            logger.debug("no botm to derive nlay from", exc_info=True)
            return None
        if botm.ndim == 0 or botm.size == 0:
            return None
        # A 1-D botm is DISU, which is a single-layer grid by construction.
        return int(botm.shape[0]) if botm.ndim > 1 else 1

    def __init__(
        self,
        name: str = 'mf6_model',
        mf_folder_path: Path = Path().home().joinpath('mf6'),
        nper: int = 1,
        vor: Vor = None,
        per_dates: list | pd.DatetimeIndex = None,
        idomain_path: str | os.PathLike | None = None,
        newton: bool = True,
        complexity: str = 'MODERATE',
    ):
        """Create a legacy single-GWF model helper.

        Builds an in-memory FloPy ``MFSimulation`` (``exe_name="mf6"``), one
        ``ModflowGwf`` named ``name`` (``save_flows=True``) and one IMS solution
        registered to it. Nothing is written to disk and no TDIS or discretization
        is added here: the caller adds them, as
        :func:`~myflopy.modflow.mf6.canonical_example.build_transient_model` does.

        Parameters
        ----------
        name : str, optional
            Simulation and model name (``sim_name``, ``modelname``, and the
            ``<name>.nam`` name file); also the leaf of the workspace. Default
            ``"mf6_model"``.
        mf_folder_path : pathlib.Path, optional
            Parent folder; the workspace is ``mf_folder_path / name``. Must be a
            ``Path`` (``joinpath`` is called on it; a ``str`` raises
            ``AttributeError``). Default ``~/mf6``, evaluated once when the module
            is imported.
        nper : int, optional
            Number of stress periods, stored as ``self.nper``. Default 1. This
            constructor writes no TDIS from it.
        vor : VoronoiGridPlus or None, optional
            The model's Voronoi grid, stored as ``self.vor``; the grid-based
            helpers (idomain, plotting, regions) read it. Default ``None``.
        per_dates : list or pandas.DatetimeIndex or None, optional
            Stress-period date labels, stored as :attr:`per_dates`. A list is
            parsed with ``pandas.to_datetime``; if that raises ``ValueError`` a
            warning is logged and the list is kept as given. A
            ``DatetimeIndex`` is kept as is; any other type fails an ``assert``.
            Default ``None``.
        idomain_path : str, os.PathLike or None, optional
            Layer-tagged polygon file used to build :attr:`idomain` lazily, the
            first time it is read, when no idomain array was supplied: each
            feature's ``layer`` tag says which layers it switches OFF. Stored as a
            ``pathlib.Path``; another type raises ``TypeError`` (ledger 457).
            Default ``None``.
        newton : bool, optional
            Default ``True`` writes the GWF name-file option ``NEWTON
            UNDER_RELAXATION`` (``newtonoptions="under_relaxation"``); ``False``
            writes no Newton option.
        complexity : str, optional
            IMS ``complexity`` passed to FloPy's ``ModflowIms`` (MF6 accepts
            ``"SIMPLE"``, ``"MODERATE"``, ``"COMPLEX"``). Default ``"MODERATE"``.
            The constructor also passes explicit IMS settings beside it (DBD
            under-relaxation, backtracking, ``outer_maximum=800``,
            ``inner_maximum=300``, ``outer_dvclose=1e-5``, ``inner_dvclose=1e-6``,
            BICGSTAB with RCM reordering and L2NORM scaling).
        """
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
        """The PEST runs done on this model (``<workspace>.pest``), ready to review.

        Returns a list of :class:`~myflopy.modflow.mf6.pest.runs.PestRunHandle`;
        call ``.review()`` on one to open it as :class:`IesResults` (phi, ensemble
        vs. observations, forecasts, parameter-field maps). Empty until a
        :class:`~myflopy.modflow.mf6.pest.PestProject` has been built/run with its
        default workspace beside this model. Calibrations built before
        2026-07-29, which lived in ``<workspace>/pest``, are found too.
        """

        from myflopy.modflow.mf6.pest.runs import find_model_pest_runs

        return find_model_pest_runs(self.workspace, model_name=self.name, model=self)

    def pest(
        self,
        name: str,
        *,
        start_datetime: str | None = None,
        workspace: str | Path | None = None,
        **kwargs,
    ):
        """Start a PEST/PEST++ calibration of this model.

        This is the write-side companion to :attr:`pest_runs` (the read side):
        :attr:`pest_runs` discovers calibrations already done on the model, while
        ``pest()`` creates a new one. It returns a fresh
        :class:`~myflopy.modflow.mf6.pest.PestProject` already bound to this model,
        with its template workspace defaulting to ``<workspace>.pest/<name>`` so the
        run is auto-discoverable afterwards. The whole loop then lives on the model::

            cal = model.pest("calib", start_datetime="2020-01-01")
            cal.parameterize("k", style="pilotpoints", pp_spacing=1500.0)
            cal.observe(head_targets)
            cal.build(); cal.run_ies()
            model.pest_runs[-1].review().plot_phi()   # no imports needed

        Parameters
        ----------
        name : str
            Calibration name (stored as ``str(name)``); the default ``.pst`` file is
            ``<name>.pst``. The default template directory is
            ``<workspace>.pest/<slug>``, where the slug is ``name`` with every run of
            characters outside ``[A-Za-z0-9_.-]`` replaced by ``_`` and leading and
            trailing ``_`` stripped (``"pest_run"`` if nothing is left).
        start_datetime : str or None, optional
            Simulation start date (e.g. ``"2020-01-01"``) that pyEMU uses to place
            time-varying parameters/observations on the time axis. Default ``None``
            reads it from the model's TDIS ``START_DATE_TIME`` (through
            ``period_calendar.tdis_start_date_time``, the one reader of it); when
            that is unset too, ``ValueError`` is raised, so pass it explicitly for a
            model without one.
        workspace : str or pathlib.Path or None, optional
            The template directory ``build()`` writes. Default ``None`` uses
            ``<workspace>.pest/<slug>``, beside the model directory, where
            :attr:`pest_runs` finds it. **A build deletes it, and the staging
            directory beside it** (``<template>.staging``), **only when each is
            absent, empty, or one a previous build wrote there** -- marked by the
            ``myflopy_pest_build.json`` the build writes into it -- so
            re-running ``build()`` rebuilds it, while an existing folder no build
            wrote (a Project root, a data folder) is refused with
            ``FileExistsError`` and left as it was (ledger 360). A ``workspace``
            that is, contains or sits inside this model's run directory is refused
            here with ``ValueError``, whatever it holds. See
            :class:`~myflopy.modflow.mf6.pest.PestProject` for the full rule.
        **kwargs
            Forwarded to :class:`~myflopy.modflow.mf6.pest.PestProject`
            (``spatial_reference``, ``zero_based``, ``longnames``).

        Returns
        -------
        PestProject
            A calibration project bound to this model.
        """

        from myflopy.modflow.mf6.pest.project import PestProject

        if start_datetime is None:
            start_datetime = tdis_start_date_time(self)
            if start_datetime is None:
                raise ValueError(
                    "start_datetime could not be read from the model's TDIS; "
                    "pass start_datetime=... explicitly (e.g. '2020-01-01')."
                )
        return PestProject(
            self, name, start_datetime=start_datetime, workspace=workspace, **kwargs
        )

    @property
    def calendar(self) -> PeriodCalendar:
        """The dates of every stress period: the TDIS read as a :class:`~myflopy.modflow.mf6.period_calendar.PeriodCalendar`.

        ``mf.period_calendar(model)``, read anew on each access: start, end and
        midpoint of each period, its calendar month (the MIDPOINT's) and its
        steady-state flag from the model's STO (ledger 529). ``.table()``
        shows them. Raises ``ValueError`` when the TDIS declares no
        ``START_DATE_TIME`` or no ``TIME_UNITS``, naming the ``mf.tdis(...)``
        argument that declares it.
        """

        return period_calendar(self)

    @property
    def package_names(self) -> list[str]:
        """Sorted list of package names currently attached to the groundwater model."""

        return sorted(self.gwf.get_package_list())

    @property
    def package_types(self) -> list[str]:
        """Each package TYPE once -- ``CHD``, ``NPF``, ``RCH`` (array recharge is RCH).

        :attr:`package_names` are the packages' own names (``CHD_0``, ``CHD_1``);
        this is what they are (ledger 228).
        """

        return sorted({_package_type_name(package.package_type) for package in self.gwf.packagelist})

    def _package_pairs(self) -> list[tuple[str, str]]:
        """``(name, type)`` for each of :attr:`package_names`, in its order.

        Driven by the names, so ``packages.summary()`` keeps its contract -- a row per
        listed package, in model order -- and gains each one's type. A listed name
        with no package behind it is its own type.
        """

        types = {
            str(package.package_name).upper(): _package_type_name(package.package_type)
            for package in self.gwf.packagelist
        }
        return [(str(name).upper(), types.get(str(name).upper(), str(name).upper())) for name in self.package_names]

    @property
    def targets(self) -> TargetRegistry:
        """Model-bound registry for reusable calibration target sets."""

        if self._targets is None:
            from myflopy.modflow.mf6.observations import TargetRegistry

            self._targets = TargetRegistry(self)
        return self._targets

    @property
    def particle_tracking(self) -> ParticleTracking:
        """Model-bound MF6 PRT and MP3DU workflow entry point."""

        if self._particle_tracking is None:
            from myflopy.modflow.mf6.prt import ParticleTracking

            self._particle_tracking = ParticleTracking(self)
        return self._particle_tracking

    @property
    def parallel(self) -> ParallelModelWorkflow:
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
        """Set the optional path used when building idomain from external data.

        A ``str`` or ``os.PathLike`` is stored as a ``pathlib.Path``, ``None`` clears
        it, and anything else raises ``TypeError`` naming it (ledger 457): a
        ``str`` used to fail an ``assert``, and ``None`` left the old path in place.
        """

        self._idomain_path = as_path(idomain_path, name="idomain_path")

    @property
    def per_dates(self):
        """Optional stress-period date labels, set by hand.

        A map's hover footer labels its period from these when they are set, and
        otherwise from the TDIS clock, :attr:`calendar` (the period's midpoint;
        ledger 584). Nothing else reads them.
        """

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
        """Convert ``(layer, index_in_layer)`` back to a flattened node id.

        The node is ``offsets[layer] + idx_in_layer``, where ``offsets`` is the
        cumulative per-layer cell count of the model grid with a leading 0. The
        inverse is :attr:`node_to_lni`.

        Parameters
        ----------
        layer : int
            Zero-based layer index; it indexes the offsets array directly, so a
            negative value counts from the end of that array (which has
            ``nlay + 1`` entries) rather than being refused.
        idx_in_layer : int
            Zero-based cell index within ``layer``. Added as-is: it is not checked
            against the layer's cell count, so an out-of-range value returns a
            node in another layer (or past the grid) without an error.

        Returns
        -------
        int
            The zero-based flattened node number.
        """

        return int(self._offsets[layer] + idx_in_layer)

    @property
    def cellids(self):
        """Convenience list of all flattened cell ids in the model."""

        return list(self.node_to_lni.keys())

    def _require_model_kind(self, model_type: str, accessor: str, field: str) -> None:
        """Raise a clear error if this view is not the model kind ``accessor`` needs."""

        actual = getattr(self, "model_type", "gwf6")
        if actual != model_type:
            raise AttributeError(
                f"model {self.name!r} is a {actual[:3].upper()} model; '.{accessor}' "
                f"({field}) is only available on {model_type[:3].upper()} models."
            )

    @property
    def hds(self):
        """FloPy heads reader for the current model outputs (GWF models)."""

        self._require_model_kind("gwf6", "hds", "heads")
        return get_hds(self)

    @property
    def conc(self):
        """FloPy concentration reader for the current model outputs (GWT models)."""

        self._require_model_kind("gwt6", "conc", "concentration")
        return get_conc(self)

    @property
    def temp(self):
        """FloPy temperature reader for the current model outputs (GWE models)."""

        self._require_model_kind("gwe6", "temp", "temperature")
        return get_temp(self)

    @property
    def all_heads(self):
        """Tabular heads view across available results."""

        return get_all_heads(self)

    @property
    def all_conc(self):
        """Tabular concentration view across available results (GWT models).

        Kind-GATED, unlike :attr:`all_heads`. The ungated spelling is a legacy
        the heads table carries because it predates the kind guard; a new
        transport table should not inherit it, or ``group.conc`` over a group
        that accidentally holds a flow model reads a missing ``.ucn`` and fails
        somewhere far from the cause.
        """

        self._require_model_kind("gwt6", "all_conc", "concentration")
        return get_all_conc(self)

    @property
    def all_temp(self):
        """Tabular temperature view across available results (GWE models).

        Kind-gated for the same reason as :attr:`all_conc`.
        """

        self._require_model_kind("gwe6", "all_temp", "temperature")
        return get_all_temp(self)

    @property
    def _field_reader(self):
        """This model's own dependent-variable reader, whatever kind it is.

        Private and deliberately UNGATED, unlike ``.hds``/``.conc``/``.temp``:
        it exists for kind-agnostic machinery that needs *a* time axis or field
        without knowing (or caring) which physics it is looking at.
        """

        return field_reader(self)

    def to_xugrid(
        self,
        *,
        layers=None,
        times=None,
        name: str = "head",
        masked: bool = True,
        surfaces: bool = False,
    ):
        """Export simulated heads as an xugrid ``(time, layer, cell)`` object.

        Convenience that stacks this model's saved heads across all layers and
        output times into one :class:`xugrid.UgridDataArray` on the Voronoi mesh
        -- so you never assemble a ``(nlay, ncpl)`` array by hand. The result
        plugs into xarray slicing, native unstructured plotting, and UGRID-NetCDF
        sharing (QGIS/ParaView). Delegates to
        :meth:`~myflopy.modflow.mf6.headsplus.HeadsPlus.to_xugrid` on :attr:`hds`,
        so it is available on a GWF model only (``hds`` raises ``AttributeError``
        on a GWT/GWE model).

        Parameters
        ----------
        layers : iterable of int or None, optional
            Zero-based layers to include, in the order given (each is passed
            through ``int``). Default ``None`` includes all ``nlay`` layers.
        times : iterable of (int, int) or None, optional
            ``(kstp, kper)`` keys to include, in the order given, spelled as the
            heads reader's own ``kstpkper`` lists them. Default ``None`` includes
            every saved time. A key that is not a saved time raises ``ValueError``.
        name : str, optional
            Variable name for the field. Default ``"head"``; an empty string falls
            back to the reader's own ``value_name``.
        masked : bool, optional
            Default ``True`` replaces MODFLOW's no-value sentinels (``|value| >=
            1e29``: +1e30 inactive, -1e30 dry) with ``NaN``; ``False`` keeps the
            raw file values.
        surfaces : bool, optional
            Default ``False``. ``True`` also attaches the grid's layer frame --
            ``top`` and ``botm_1`` ... ``botm_n``, one flat variable per layer --
            and returns a ``UgridDataset`` instead of a ``UgridDataArray``.

        Returns
        -------
        xugrid.UgridDataArray or xugrid.UgridDataset
            Dims ``("time", "layer", <face_dim>)`` with ``kstp``/``kper`` coords; a
            ``UgridDataset`` carrying the field plus the flat elevation variables
            when ``surfaces=True``.

        Raises
        ------
        ImportError
            If the optional ``xugrid`` / ``xarray`` packages are not installed.
        ValueError
            If a ``times`` key is not a saved output time.

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

        ``surfaces=True`` also attaches the grid's layer frame (``top``,
        ``botm_1`` ...), one flat variable per layer. For QGIS, flatten the field
        in the same call -- a mesh dataset there is one value per face, and a
        variable that still carries ``time`` or ``layer`` is dropped silently:

        >>> uds = model.to_xugrid(surfaces=True)
        >>> uds.isel(time=-1, layer=0).ugrid.to_netcdf("model.nc")
        """

        return self.hds.to_xugrid(
            layers=layers, times=times, name=name, masked=masked, surfaces=surfaces
        )

    @property
    def plot(self) -> ModelPlots:
        """All six plotting verbs for this model.

        ``map`` / ``section`` / ``surface`` / ``grid`` draw a picture of this
        model; ``mosaic`` / ``animate`` are **combinators** that compose pictures
        you already have, so their first argument is a collection rather than a
        subject.

        ``model.plot.map(layer=0)`` is exactly ``myflopy.plot.map(model, layer=0)``
        -- one implementation, two spellings. Everything returned is a
        :class:`~myflopy.viz.Picture`, so it renders itself in Jupyter and answers
        ``.fig`` / ``.show()`` / ``.save(path)`` / ``.html(path)``. There is never
        a trailing ``.plot()``.

        Each verb spells its parameters out rather than taking ``**kwargs``, so
        an editor can complete them; see :class:`~myflopy.plot.ModelPlots`.

        Replaces ``cor()`` (now ``plot.map()``), ``section()`` (``plot.section()``)
        and ``srf.hds()``/``srf.lyr()`` (``plot.surface()``), each of which was a
        hand-restated copy of the picture class's own signature. It also shadows
        the old one-line delegate to FloPy's ``MFSimulation.plot``; that renderer
        is still there as ``model.sim.plot(...)``.

        See Also
        --------
        myflopy.plot : the same six verbs as free functions.

        Examples
        --------
        >>> model.plot.map(layer=0, contours=True).show()
        >>> model.plot.section(cells=[1653, 651]).savefig("section.png")  # mpl default
        >>> model.plot.surface(layer=0)
        >>> model.plot.grid()                                    # the bare mesh
        >>> model.plot.grid(pathlines=tracks, backend="vtk")     # 3-D scene
        >>> model.plot.mosaic([a, b, c], ncols=2)
        """

        return ModelPlots(self)

    @property
    def inputs(self):
        """Convenience view of input packages and selected input data."""

        return get_inputs(self)

    @property
    def outputs(self) -> ModelOutputs:
        """Namespace of package-specific output helpers."""

        return get_outputs(self)

    @property
    def packages(self) -> ModelPackages:
        """Every package at once: ``summary()`` and the cross-package ``mosaic()``.

        One package is ``model.<pkg>`` (``model.rch.inputs.map()``), which
        returns the same view this namespace builds; FloPy's own package object
        is ``model.fl.<pkg>``.
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
        package_names : str or list of str or tuple of str or set of str
            One package TYPE or NAME such as ``"rch"`` or ``"lak_2"``, or a
            list/tuple/set of them. Each is lower-cased (through ``str``) and looked
            up by TYPE first -- every package of that MF6 type, array-form RCHA/EVTA
            counting as RCH/EVT -- then as a package's exact NAME (any case), on the
            model and then on its simulation (so ``"tdis"`` and ``"ims"`` resolve
            too). A collection returns a ``dict`` keyed by the lower-cased strings,
            each value looked up the same way.

        Returns
        -------
        flopy.mf6.mfpackage.MFPackage or list or dict
            The package; a list of every package of the type, in model order, when
            a TYPE matches several; a ``dict`` of those when ``package_names`` is a
            collection.

        Raises
        ------
        AttributeError
            If a string names no package type or package on the model or the
            simulation.

        Notes
        -----
        Asked by TYPE before NAME (ledger 352). FloPy's ``get_package`` asks the NAME
        first, so on a model holding LAK packages named ``lak`` and ``lak_2`` it
        answered ``"lak"`` with the one named ``lak`` -- and a type's other packages
        vanished from every reader built on this. FloPy's filename and name-PREFIX
        fallbacks are gone too: ``"dis"`` on a DISV model prefix-matched the DISV.
        The package NAMED like its type beside others is one element of the type's
        list (``model.package("lak")[0]``, or FloPy's
        ``model.gwf.get_package("lak", name_only=True)``); every other name is
        itself: ``model.package("lak_2")``.

        Examples
        --------
        >>> model.package("npf")                       # MF6 allows one NPF
        >>> [p.package_name for p in model.package("lak")]
        ['lak', 'lak_2']
        >>> model.package("lak_2").package_name
        'lak_2'
        >>> model.package(["npf", "tdis"])             # dict keyed by the strings
        """

        if isinstance(package_names, (list, tuple, set)):
            return {str(name).lower(): self.package(str(name)) for name in package_names}

        package_name = str(package_names).lower()
        package = _type_then_name(self._packages_of_type(package_name), self.gwf, package_name)
        simulation = getattr(self, "sim", None)
        if package is None and simulation is not None:
            package = _type_then_name(packages_of_type(simulation, package_name), simulation, package_name)
        if package is None:
            raise AttributeError(f"Package {package_name!r} is not attached to model '{self.name}'.")
        return package

    def _packages_of_type(self, package_type: str) -> list:
        """Every FloPy package of TYPE ``package_type`` in this model, in model order.

        The one lookup behind ``model.fl.<pkg>``, :meth:`package` and
        ``package_registry.model_packages_of_type`` (ledger 352). A file-backed run
        loads its packages first (:class:`~myflopy.project.run_model.LoadedMf6Run`).
        """

        return packages_of_type(self.gwf, package_type)

    def _flopy_of_type(self, package_type: str):
        """``model.fl.<pkg>``: the one package of the type, a list of several, by TYPE only.

        Never by name -- a package NAMED ``lak`` of another type is not a LAK -- and
        ``AttributeError`` naming the model when it has none.
        """

        found = self._packages_of_type(package_type)
        if not found:
            raise AttributeError(
                f"Model '{self.name}' has no {str(package_type).upper()} package."
            )
        return found[0] if len(found) == 1 else found

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

    # -- one package: model.<pkg> ---------------------------------------------
    # myflopy's view of each package; FloPy's own object is ``model.fl.<pkg>``.
    # Hand-written, not generated, for the reason 4.7.7 gives: a ``.pyi`` stub
    # replaces the WHOLE module for type checkers. ``test_package_descriptor.py``
    # asserts the list covers the registry.

    @property
    def fl(self) -> FlopyPackages:
        """The model's FloPy package objects: ``model.fl.drn`` is FloPy's DRN."""

        return FlopyPackages(self)

    @property
    def chd(self) -> PackageExplorer:
        """Constant-head package: ``inputs`` / ``results`` views and their maps. FloPy's object is ``model.fl.chd``."""

        return self.packages.chd

    @property
    def drn(self) -> PackageExplorer:
        """Drain package: ``inputs`` / ``results`` views and their maps. FloPy's object is ``model.fl.drn``."""

        return self.packages.drn

    @property
    def ghb(self) -> PackageExplorer:
        """General-head boundary package: ``inputs`` / ``results`` views and their maps. FloPy's object is ``model.fl.ghb``."""

        return self.packages.ghb

    @property
    def riv(self) -> PackageExplorer:
        """River package: ``inputs`` / ``results`` views and their maps. FloPy's object is ``model.fl.riv``."""

        return self.packages.riv

    @property
    def wel(self) -> PackageExplorer:
        """Well package: ``inputs`` / ``results`` views and their maps. FloPy's object is ``model.fl.wel``."""

        return self.packages.wel

    @property
    def rch(self) -> PackageExplorer:
        """Recharge package: ``inputs`` / ``results`` views and their maps. FloPy's object is ``model.fl.rch``."""

        return self.packages.rch

    @property
    def evt(self) -> PackageExplorer:
        """Evapotranspiration package: ``inputs`` / ``results`` views and their maps. FloPy's object is ``model.fl.evt``."""

        return self.packages.evt

    @property
    def uzf(self) -> UzfPackageExplorer:
        """Unsaturated-zone flow package views. FloPy's object is ``model.fl.uzf``."""

        return self.packages.uzf

    @property
    def sfr(self) -> SfrPackageExplorer:
        """Streamflow-routing package views. FloPy's object is ``model.fl.sfr``."""

        return self.packages.sfr

    @property
    def lak(self) -> LakPackageExplorer:
        """Lake package views. FloPy's object is ``model.fl.lak``."""

        return self.packages.lak

    @property
    def hfb(self) -> HfbPackageExplorer:
        """Horizontal-flow-barrier views (a barrier is a cell PAIR; no budget record)."""

        return self.packages.hfb

    @property
    def ic(self) -> StaticArrayPackageExplorer:
        """Initial-conditions array views (``strt``)."""

        return self.packages.ic

    @property
    def npf(self) -> NpfPackageExplorer:
        """Node-property-flow array views (``k``/``k22``/``k33``, ``transmissivity``)."""

        return self.packages.npf

    @property
    def sto(self) -> StaticArrayPackageExplorer:
        """Storage array views (``ss``/``sy``)."""

        return self.packages.sto

    @property
    def surface_water(self) -> SurfaceWaterPackageExplorer:
        """Combined SFR/LAK views."""

        return self.packages.surface_water

    def bud(self, package: str = None):
        """Return budget information, optionally filtered to one package.

        Parameters
        ----------
        package : str or None, optional
            A budget record to filter to, such as ``"drn"``, ``"ghb"`` or
            ``"rch"``. Matched case-insensitively as a SUBSTRING of the record
            names in this model's cell-budget file; ``"ssm"`` also reaches the
            ``SOURCE-SINK MIX`` record MF6 writes for it. A name matching no record
            raises ``ValueError`` listing the records present. Default ``None``
            returns the unfiltered helper, whose ``budget`` is the raw budget-file
            reader.

        Returns
        -------
        myflopy.modflow.mf6.budget.Budget
            The budget helper; with ``package`` set, its ``df`` is that record's
            per-cell table.

        Raises
        ------
        ValueError
            If ``package`` matches no record in the budget file.
        """

        return get_budget(self, package)

    @property
    def budget(self):
        """Every term in this model's own budget file, as a noun.

        ``model.budget.<term>.get()/summary()/plot()/map()`` -- the terms are
        discovered from the file, so they follow the model's kind and packages
        (``model.budget.source_sink_mix`` on a transport model,
        ``model.budget.drn`` or ``model.budget.sto_ss`` on a flow model). See
        :class:`~myflopy.modflow.mf6.package_results.ModelBudgetNamespace` for
        how this differs from ``bud()``, ``budget_cumulative``, and the
        package-output ``model.<pkg>.budget``.
        """

        return get_model_budget(self)

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

    def _project_roots(self) -> list[Path]:
        """Where this model's project may be: the project it physically LIVES IN,
        then the root its run recorded.

        The first is the nearest ancestor of the workspace that is a live project
        root or holds a ``project.json``. It comes FIRST because a copied project's
        runs still record the ORIGINAL root, so trusting the record gave a copy's
        models the original's areas -- and ``add_aoi`` wrote into the original
        (ledger 312). Nearest, not ``<workspace>/../..``: a run can live deeper
        (``runs/scenarios/x``, ``imports/usg/converted``), and a discovered run or
        a ``load_mf6_run`` model records no root at all.
        """

        roots: list[Path] = []
        workspace = getattr(self, "workspace", None)
        if workspace is not None:
            for ancestor in Path(workspace).resolve().parents:
                if live_areas(ancestor) is not None or (ancestor / "project.json").is_file():
                    roots.append(ancestor)
                    break
        if self._project_root is not None:
            recorded = Path(self._project_root).resolve()
            if recorded not in roots:
                roots.append(recorded)
        return roots

    @property
    def aoi(self) -> AoiLibrary:
        """This model's named areas of interest -- a dict of :class:`~myflopy.areas.Area`.

        What ``plot.map(extent="site")`` and ``select="site"`` look names up in.
        Defined here or on the project, and SHARED with the project (ledger 310):

        * a model whose project is LIVE in this session -- it came from
          ``prepare_run``, ``reopen_run``, ``project.runs[...]`` or ``load_run`` of
          a run inside it -- answers with THAT project's library, so
          ``model.add_aoi(...)`` is ``project.add_aoi(...)`` and ``project.save()``
          keeps it;
        * with no live project, the project's SAVED areas are read once, and
          anything added lasts for this session only (only ``project.save()``
          writes a project's files);
        * a model in no project keeps its areas in memory, like ``regions``.
        """

        roots = self._project_roots()
        for root in roots:
            live = live_areas(root)
            if live is not None:
                return live
        if self._session_aoi is None:
            library = AoiLibrary()
            for root in roots:
                saved = saved_areas(root)
                if saved is not None:
                    library.update(saved)
                    break
            self._session_aoi = library
        return self._session_aoi

    def add_aoi(
        self,
        name: str,
        geometry,
        *,
        crs=None,
        layer: str | None = None,
        note: str | None = None,
        overwrite: bool = False,
    ) -> Area:
        """Name an area of interest on this model -- and on its project, if it has one.

        ``geometry`` is four numbers ``(xmin, ymin, xmax, ymax)``, a run of
        ``(x, y)`` points, any shapely geometry, a GeoDataFrame, or a path to a
        vector file (``layer=`` picks a GeoPackage layer); ``crs=None`` means the
        grid's frame. See :attr:`aoi` for where it is kept. Refuses an existing
        name unless ``overwrite=True``; ``model.aoi[name] = ...`` replaces.

        Parameters
        ----------
        name : str
            The area's name (converted with ``str()``): the key in
            :attr:`aoi` and what ``extent=``/``select=`` accept.
        geometry : tuple, sequence of points, shapely geometry, GeoDataFrame, str or Path
            Four numbers ``(xmin, ymin, xmax, ymax)`` in shapely's order (a box),
            an ``(N, 2)`` run of ``(x, y)`` points (the box bounding them), a
            shapely geometry, a GeoDataFrame/GeoSeries (its rows unioned), an
            :class:`~myflopy.areas.Area`, or a path to a vector file (read with
            ``geopandas.read_file`` and unioned). Coerced by
            :meth:`~myflopy.areas.Area.coerce`, which raises ``TypeError`` for
            anything else and ``ValueError`` for an empty frame or an inverted
            box.
        crs : str, int, pyproj.CRS or None, optional
            The CRS ``geometry`` is in. Default ``None`` means the grid's frame for
            numbers, points and a bare geometry, or the CRS a file or frame
            carries. A ``crs`` that disagrees with the CRS a file, frame or
            :class:`~myflopy.areas.Area` already carries raises ``ValueError``.
        layer : str or None, optional
            Layer to read when ``geometry`` is a path to a multi-layer file (a
            GeoPackage); ignored for every other form. Default ``None`` reads the
            file's first layer.
        note : str or None, optional
            Free text stored on the area and shown in ``model.aoi.summary()``.
            Default ``None``.
        overwrite : bool, optional
            Default ``False`` raises ``KeyError`` when ``name`` is already an area
            in :attr:`aoi` (the project's library when the project is live);
            ``True`` replaces it.

        Returns
        -------
        Area
            The stored area.

        Examples
        --------
        >>> model.add_aoi("site", (500, 600, 1000, 1100))
        >>> model.plot.map(layer=0, extent="site", backend="mpl")
        >>> model.hds.mosaic(by="layer", extent="site")
        """

        library = self.aoi
        return library._define(
            name, geometry, crs=crs, layer=layer, note=note, overwrite=overwrite,
            where=library._where() if library._project_name else f"Model '{self.name}'",
        )

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

        The region is kept in this model's in-memory :attr:`regions` registry
        (not saved with a project). Regions and region groups share one
        namespace.

        Parameters
        ----------
        name : str
            Region name. A name already used by a region GROUP raises
            ``ValueError`` whatever ``overwrite`` says.
        cellids : iterable of int or tuple
            The cells. An ``int`` is an in-layer cell number as given, and a
            one-item tuple ``(cell,)`` the same -- on DISU that is MF6's own
            cellid ``(node,)``, which this model's heads table and cell polygons
            number as cell ``node`` of layer 0; an empty tuple is dropped. A
            longer tuple is the model grid's own cellid, read by the rule every
            package table uses (ledger 419): ``(layer, cell)`` on DISV,
            ``(layer, row, col)`` on DIS, whose cell number is FloPy's row-major
            ``row * ncol + col``; on DISU a longer tuple is refused naming
            ``(node,)``. A tuple of another length raises ``ValueError`` naming
            the expected form (ledger 376: a ``(layer, row, col)`` used to
            register the ROW as the cell), and so does any layer, row, column or
            cell off the grid -- a bare cell number included (ledger 419).
            The region stores the sorted unique cell numbers, plus the raw ids as
            given. :meth:`region_heads` matches them against the heads table's
            in-layer ``cell`` column, so give in-layer (per-layer) cell indices.
        layer : int or list of int or None, optional
            Zero-based layer (or layers) the region applies to; :meth:`region_heads`
            filters to it when called without its own ``layer``. A layer off the
            grid raises ``ValueError`` (on DISU the only layer is 0). Default
            ``None`` takes the layers named by tuple ``cellids`` (one layer as an
            ``int``, several as a sorted list), else leaves the region layer-free.
        category : str, optional
            Free label stored on the region and shown in :meth:`list_regions`.
            Default ``"selection"``.
        package : str or None, optional
            Package name stored on the region as a label and shown in
            :meth:`list_regions`; nothing is looked up. Default ``None``.
        tags : list of str or None, optional
            Labels stored (copied into a list) on the region. Default ``None``
            stores none.
        geometry : object, optional
            Any geometry to keep on the region for reference; not used to select
            cells. Default ``None``.
        metadata : dict or None, optional
            Free key/values stored (copied into a dict) on the region. Default
            ``None`` stores an empty dict.
        overwrite : bool, optional
            Default ``False`` raises ``ValueError`` when a region called ``name``
            exists; ``True`` replaces it (group memberships naming it are kept).

        Returns
        -------
        ModelRegion
            The registered region.

        Raises
        ------
        ValueError
            If ``name`` is a group's name, an existing region's without
            ``overwrite=True``, a tuple cellid is not of the grid's form, or a
            layer or cell is off the grid.
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
        """Alias for :meth:`add_region` with explicit cell ids.

        Parameters
        ----------
        name : str
            Region name, as in :meth:`add_region`.
        cellids : iterable of int or tuple
            The cells, as in :meth:`add_region` (here positional).
        **kwargs
            Forwarded to :meth:`add_region`: ``layer``, ``category``,
            ``package``, ``tags``, ``geometry``, ``metadata``, ``overwrite``.
            Anything else raises ``TypeError``.

        Returns
        -------
        ModelRegion
            The registered region.
        """

        return self.add_region(name, cellids=cellids, **kwargs)

    def add_region_from_geometry(
        self,
        name: str,
        geometry,
        *,
        predicate: str | None = None,
        layer: int | list[int] | None = None,
        category: str = "selection",
        package: str | None = None,
        tags: list[str] | None = None,
        metadata: dict | None = None,
        overwrite: bool = False,
    ):
        """Register a named region by selecting cells from a geometry.

        The cells are selected on this model's grid (``self.vor.gdf_vorPolys``),
        flattened into one list, and registered through :meth:`add_region`, so
        every rule there applies. The cells are the grid's in-layer cell indices.

        Parameters
        ----------
        name : str
            Region name, as in :meth:`add_region`.
        geometry : shapely geometry, GeoDataFrame, GeoSeries, str or pathlib.Path
            What selects the cells: any shapely geometry, a GeoDataFrame or
            GeoSeries, or a path to a vector file (its first layer is read). A
            geometry carrying a CRS is reprojected to the grid's; one carrying
            none is read in the grid's frame. Any other type raises
            ``ValueError``. It is also stored on the region.
        predicate : str or None, optional
            Default ``None`` selects by MEASURE: a polygon reaches the cells it
            covers by area, a line those it runs through, a point the one cell
            holding it (a feature reaching no cell is left out and named in a
            logged warning). A
            geopandas predicate (``"intersects"``, ``"contains"``, ...) uses the
            raw spatial join instead, which counts any touch, a shared face or a
            single vertex included.
        layer : int or list of int or None, optional
            Zero-based layer (or layers) stored on the region, as in
            :meth:`add_region`. Default ``None``.
        category : str, optional
            Free label stored on the region. Default ``"selection"``.
        package : str or None, optional
            Package name stored on the region as a label. Default ``None``.
        tags : list of str or None, optional
            Labels stored on the region. Default ``None`` stores none.
        metadata : dict or None, optional
            Free key/values stored on the region. Default ``None`` stores an
            empty dict.
        overwrite : bool, optional
            Default ``False`` raises ``ValueError`` when a region called ``name``
            exists; ``True`` replaces it.

        Returns
        -------
        ModelRegion
            The registered region.

        Raises
        ------
        ValueError
            If ``geometry`` is of an unsupported type or reaches no cell, or for
            a name collision as in :meth:`add_region`.
        """

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
        """Return one registered region object by name.

        Parameters
        ----------
        name : str
            A region's name; a group's name, or an unknown one, raises
            ``KeyError`` (use :meth:`get_group` for a group).

        Returns
        -------
        ModelRegion
            The registered region.
        """

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
        """Create a named region group from region/group members.

        A group holds region and group NAMES, resolved to cells only when read
        (:meth:`resolve_region_cells`), so a later change to a member region shows
        through.

        Parameters
        ----------
        name : str
            Group name. A name already used by a REGION raises ``ValueError``
            whatever ``overwrite`` says.
        members : iterable of str or None, optional
            Region or group names to add, in order, each through
            :meth:`add_to_group`: an unknown name raises ``KeyError`` and one
            that would make a cycle raises ``ValueError``. A refused member
            leaves the registry exactly as it was before the call -- no group
            ``name`` if there was none, the previous group under
            ``overwrite=True`` (ledger 375; it used to keep the group
            half-filled, the replaced one already gone). Default ``None``
            creates an empty group.
        tags : list of str or None, optional
            Labels stored (copied into a list) on the group. Default ``None``
            stores none.
        metadata : dict or None, optional
            Free key/values stored (copied into a dict) on the group. Default
            ``None`` stores an empty dict.
        overwrite : bool, optional
            Default ``False`` raises ``ValueError`` when a group called ``name``
            exists; ``True`` replaces it with the new group, once every member
            has been accepted.

        Returns
        -------
        RegionGroup
            The registered group.

        Raises
        ------
        ValueError
            If ``name`` is a region's name, an existing group's without
            ``overwrite=True``, or a member would create a cycle.
        KeyError
            If a ``members`` name is neither a region nor a group.
        """

        return self.regions.add_group(
            name,
            members=members,
            tags=tags,
            metadata=metadata,
            overwrite=overwrite,
        )

    def add_to_group(self, group_name: str, member_name: str):
        """Add one region or group to an existing region group.

        Parameters
        ----------
        group_name : str
            The group to add to; an unknown group raises ``KeyError``.
        member_name : str
            A registered region or group name; an unknown name raises
            ``KeyError``, ``group_name`` itself or a group that already contains
            ``group_name`` (a cycle) raises ``ValueError``. Adding a member already
            present changes nothing.

        Returns
        -------
        RegionGroup
            The updated group.
        """

        return self.regions.add_to_group(group_name, member_name)

    def remove_from_group(self, group_name: str, member_name: str):
        """Remove one member from an existing region group.

        Only the membership is removed; the member region or group itself stays
        registered (see :meth:`remove_region` to delete it).

        Parameters
        ----------
        group_name : str
            The group to remove from; an unknown group raises ``KeyError``.
        member_name : str
            The member to drop. A name that is not a member is ignored without
            an error.

        Returns
        -------
        RegionGroup
            The updated group.
        """

        return self.regions.remove_from_group(group_name, member_name)

    def get_group(self, name: str):
        """Return one registered region group by name.

        Parameters
        ----------
        name : str
            A group's name; a region's name, or an unknown one, raises
            ``KeyError`` (use :meth:`get_region` for a region).

        Returns
        -------
        RegionGroup
            The registered group.
        """

        return self.regions.get_group(name)

    def get_region_cells(self, name: str) -> list[int]:
        """Return the resolved cell ids for a named region or group.

        The same answer as :meth:`resolve_region_cells`.

        Parameters
        ----------
        name : str
            A region or group name; a group is resolved through its members,
            recursively. An unknown name raises ``KeyError``.

        Returns
        -------
        list of int
            The sorted unique cell numbers the regions store.
        """

        return get_region_cells(self, name)

    def list_regions(self):
        """List registered regions on the model."""

        return list_model_regions(self)

    def list_groups(self):
        """List registered region groups on the model."""

        return list_model_groups(self)

    def resolve_region_cells(self, name: str) -> list[int]:
        """Resolve a region/group to its unique cell ids.

        Parameters
        ----------
        name : str
            A region or group name; a group is resolved through its members,
            recursively (a cycle raises ``ValueError``). An unknown name raises
            ``KeyError``.

        Returns
        -------
        list of int
            The sorted union of the cell numbers the reached regions store --
            in-layer indices for a region built from geometry or from
            ``(layer, cell)`` tuples, the numbers as given otherwise.
        """

        return self.regions.resolve_cells(name)

    def resolve_region_cells_with_trace(self, name: str):
        """Resolve a region/group and return overlap/contribution trace data.

        Parameters
        ----------
        name : str
            A region or group name, resolved as in :meth:`resolve_region_cells`;
            an unknown name raises ``KeyError``, a cycle ``ValueError``.

        Returns
        -------
        tuple of (list of int, dict of int to list of str)
            The sorted unique cells, and for each cell the sorted names of the
            REGIONS that contributed it (a cell in several regions lists each).
        """

        return self.regions.resolve_cells_with_trace(name)

    def remove_region(self, name: str):
        """Remove a registered region or region group by name.

        The name is also dropped from every group's members.

        Parameters
        ----------
        name : str
            A region or group name; an unknown name raises ``KeyError``.

        Returns
        -------
        ModelRegion or RegionGroup
            The removed object.
        """

        return self.regions.remove(name)

    def region_heads(
        self,
        name: str,
        *,
        per: int | None = None,
        kstpkper: tuple | None = None,
        layer: int | list[int] | None = None,
    ):
        """Return head results filtered to a named region or region group.

        Reads :attr:`all_heads` (every saved time, layer and cell), keeps the rows
        whose in-layer ``cell`` is one of the region's resolved cells, and attaches
        each cell's polygon. It opens ``<workspace>/<name>.hds`` with no
        model-kind check (unlike :attr:`hds`), so it is meaningful on a GWF model
        only.

        Parameters
        ----------
        name : str
            A region or group name; an unknown name raises ``KeyError``.
        per : int or None, optional
            Zero-based stress period, read at its LAST saved step on the head
            file's clock -- the one time ``kstpkper`` would name, by ledger 337's
            rule (``period_end.last_saved_step``). It kept EVERY saved step of the
            period until ledger 348, so at ``nstp=3`` a region's table was three
            times as long and mixed three heads per cell. A period the head file
            saved nothing in raises ``ValueError``. Ignored when ``kstpkper`` is
            given. Default ``None`` keeps all periods.
        kstpkper : tuple of (int, int) or None, optional
            One zero-based ``(kstp, kper)`` saved time to keep; takes precedence
            over ``per``. A time not in the file gives an empty frame, not an
            error. Default ``None``.
        layer : int or list of int or None, optional
            Zero-based layer (or list of layers) to keep. Default ``None`` uses the
            layer stored on the region (a GROUP carries none, so every layer is
            kept), and keeps every layer when that is ``None`` too. A layer the
            model's views do not number -- ONE on a DISU, whatever FloPy infers
            from ``IHC`` -- raises ``ValueError`` naming the count (ledger
            496); it filtered to an empty frame.

        Returns
        -------
        geopandas.GeoDataFrame
            One row per kept ``(kstpkper, layer, cell)``, with the head column, the
            cell ``geometry`` in the grid's CRS, and a ``region`` column holding
            ``name``. A cell MF6 wrote as dry (-1e30) or inactive (+1e30) reads
            ``NaN``, its row kept, by the head file's rule (ledger 315); it read
            the number until ledger 495, so a region mean averaged it in.
        """

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
        sfr: SFRBuilder | None = None,
        maxmvr: int | None = None,
        maxpackages: int | None = None,
        mvr_packages: list | None = None,
        mvr_perioddata: dict | None = None,
        raise_on_error: bool = False,
    ) -> SurfaceWaterValidationReport:
        """Validate coupled LAK/SFR/MVR inputs against the current model.

        A thin wrapper over
        :func:`~myflopy.modflow.mf6.surface_water_validation.validate_surface_water_configuration`
        with this model as ``model``; that function states every rule checked.
        Each of the three checks runs only when its inputs are given: LAK when any
        of ``nlakes``/``lak_packagedata``/``lak_connectiondata``/``lak_perioddata``
        is not ``None``, SFR when ``sfr`` is, MVR when ``mvr_packages`` or
        ``mvr_perioddata`` is. Cell ids, lake ids, reach numbers and mover ids are
        zero-based.

        Parameters
        ----------
        nlakes : int or None, optional
            Number of lakes the LAK package declares. Default ``None`` is read as
            ONE lake whenever another LAK input is given, so pass it for a
            multi-lake package.
        lak_packagedata : list or None, optional
            LAK ``packagedata`` rows ``[lakeno, strt, nlakeconn, ...]``. Default
            ``None``.
        lak_connectiondata : list or None, optional
            LAK ``connectiondata`` rows ``[lakeno, iconn, cellid, claktype,
            bedleak, belev, telev, connlen, connwidth]``, ``cellid`` a
            ``(layer, cell)`` tuple. Default ``None``.
        lak_perioddata : dict or None, optional
            LAK ``perioddata`` as ``{period: [[lakeno, keyword, value, ...], ...]}``.
            Default ``None``.
        sfr : SFRBuilder or None, optional
            An :class:`~myflopy.modflow.mf6.sfr.SFRBuilder` checked through its
            derived package, connection and period data; reading those resolves
            the builder, which can raise. Default ``None`` skips the SFR check.
        maxmvr : int or None, optional
            MVR ``MAXMVR``. ``None`` (the default) is checked as 0, so give it
            whenever ``mvr_perioddata`` holds records.
        maxpackages : int or None, optional
            MVR ``MAXPACKAGES``. ``None`` (the default) is checked as 0, so give
            it whenever ``mvr_packages`` declares packages.
        mvr_packages : list or None, optional
            The MVR ``packages`` block as a list of rows whose first item is a
            package name, e.g. ``[["sfr"], ["lak"]]``. Default ``None``.
        mvr_perioddata : dict or None, optional
            MVR ``perioddata`` as ``{period: [(pname1, id1, pname2, id2, mvrtype,
            value), ...]}``. Default ``None``.
        raise_on_error : bool, optional
            Default ``False`` returns the report with its findings collected.
            ``True`` calls ``report.raise_for_errors()``, which raises
            ``ValueError`` listing every error-level finding; warnings alone never
            raise.

        Returns
        -------
        SurfaceWaterValidationReport
            Every finding; empty when nothing was given to check.

        Raises
        ------
        ValueError
            With ``raise_on_error=True``, when the report holds an error-level
            finding.
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
        crs: str | None = None,
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

        Parameters
        ----------
        other : SimulationBase, str, pathlib.Path, list or tuple
            The first model to compare: a model object (used as given), or a
            workspace path opened lazily with :func:`load_mf6_run` as a
            file-backed :class:`LoadedMf6Run`. A ``list``/``tuple`` is flattened
            one level into its items.
        *others : SimulationBase, str, pathlib.Path, list or tuple
            More models, in the same forms as ``other`` and flattened the same way;
            together with ``other`` they become the non-reference members of the
            :class:`~myflopy.project.model_group.ModelGroup` built for the diff.
            Each model's key in that group -- what ``diff.model_names``, the
            ``model`` column of every diff table and the map titles show -- is its
            ``name`` when every name differs. When two share one (two runs of one
            model, the usual case) every member is keyed by the folder that
            tells the workspaces apart: the workspace's own folder name when
            those differ (a project run's is its run name, a loaded run's the
            folder it was opened from), else its parent's (a built model's
            workspace is ``<folder>/<model name>``), and so on up; one workspace
            given twice is keyed ``"reference"`` and ``"other"`` (``"other_1"``,
            ``"other_2"``, ... for several). Ledger 373: a model named like
            this one used to REPLACE it, and the diff compared nothing.
        crs : str or None, optional
            CRS (any string pyproj accepts) given to :func:`load_mf6_run` for the
            workspace-path entries; model objects are not touched. Default
            ``None`` takes each run's CRS recorded in its ``run.json``, else
            :data:`~myflopy._defaults.DEFAULT_CRS` (``"EPSG:2927"``), as
            ``load_mf6_run(ws)`` does; a string given here wins (ledger 374).
        shared_grid : bool, optional
            Default ``False``. ``True`` requires every model to carry the same
            discretization (compared by grid signature, ``ValueError`` otherwise)
            and reuses this model's grid view for all of them.
        verbosity_level : int, optional
            Given to :func:`load_mf6_run` for the workspace-path entries. Default
            ``0`` is silent; above 0 prints the loader's progress.

        Returns
        -------
        ModelDiff
            A :class:`~myflopy.project.model_diff.ModelDiff` whose reference is
            this model, under the key ``other``/``others`` describes (its
            ``name`` unless a name repeats).
        """

        from myflopy.project.group.core import _group_for_diff

        extra: list = []
        for item in (other, *others):
            if isinstance(item, (list, tuple)):
                extra.extend(item)
            else:
                extra.append(item)
        group = _group_for_diff(
            self, extra, crs=crs, shared_grid=shared_grid, verbosity_level=verbosity_level
        )
        return group.diff()

    @property
    def config(self):
        """Return a normalized configuration view of this model.

        A :class:`~myflopy.project.model_config.ModelConfig` capturing the
        scalar settings that control how the run behaves -- TDIS timing
        (including per-period ``perlen``/``nstp``/``tsmult``), the IMS solver
        block, OC save/print records, the MODEL's own options block
        (``config.model`` -- where ``NEWTON`` lives, ledger 286) and every package's
        OPTIONS -- as one comparable ``(section, setting, value)`` table
        (``config.settings()``). This is what the
        :class:`~myflopy.project.model_diff.ModelDiff` config tier diffs across models.
        """

        return _model_config().from_model(self)

    @property
    def model_options(self) -> dict:
        """The MODEL-level options this model declares, as ``mf.gwf(...)`` keyword arguments.

        Rebuilding specs from a finished workspace reads each PACKAGE off the loaded
        model -- and the model's own options block belongs to none of them, so it was
        dropped in silence. On a real 18,107-cell model that lost ``NEWTON``: every cell
        is convertible, **90,408 cells converted to DRY** straight from the initial
        conditions, recharge on a dry column is skipped so the period-1 rate fell from
        2,854,736 to 1,645 ft3/d, every head-dependent boundary went to zero, and the
        solver diverged with a 200% budget discrepancy. Restoring the one option
        reproduced the original run's iteration trace digit for digit (ledger 286).

        Spread it into the rebuild::

            run = mf.load_mf6_run(workspace)
            model = mf.gwf(name, context=ctx, packages=[...], **run.model_options)

        Only options the model actually SETS are returned, so the result is also a
        readable answer to "what is switched on here".

        **Filename-valued options carry the SOURCE model's names** (``list``, the
        ``nc_*`` file records). Rebuilding under a new model name, pass the ones you
        want yourself -- they are returned rather than dropped, because a silent drop is
        the defect this exists to fix.
        """

        return _model_config().model_options(getattr(self, "gwf", None))

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

    def run_simulation(self):
        """Write and run the current simulation.

        A lake written with LAK IMPLICIT is checked against the executable first: one
        that does not read the option is WARNED about, and the run goes ahead (ledger
        326).
        """

        warn_if_implicit_unsupported(self.sim)
        return run_model_simulation(self)
