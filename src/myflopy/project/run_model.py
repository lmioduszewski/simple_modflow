"""File-backed model loader built on top of FloPy MF6 loading utilities."""

from __future__ import annotations

import hashlib
from pathlib import Path
from types import MethodType

import flopy
import numpy as np

from myflopy._defaults import DEFAULT_CRS
from myflopy.modflow.mf6.grid.domain import idomain_from_flopy
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus, grid_view_from_modelgrid
from myflopy.modflow.mf6.lak_implicit import keep_implicit, warn_if_implicit_unsupported
from myflopy.modflow.mf6.output_files import (
    _declared_models,
    _declared_output_path,
    _model_dir,
    recorded_crs,
)
from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS, _package_type_name
from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.simulation.regions import RegionRegistry


def _progress(verbosity_level: int, message: str):
    """Emit a simple progress message when verbose loading is enabled."""

    if verbosity_level > 0:
        print(f"[LoadedMf6Run] {message}")


def _infer_model_name(workspace: Path) -> str:
    """The FLOW model's name for an MF6 workspace: the first ``gwf6`` in ``mfsim.nam``.

    The alphabetically first model name file used to win, and in a coupled run that
    is as likely the GWT, GWE or PRT model -- ``load_mf6_run`` then opened transport
    as if it were flow. The name files decide only without a ``mfsim.nam``.
    """

    for model_type, name in _declared_models(workspace):
        if model_type == "gwf6":
            return name

    name_files = sorted(path for path in workspace.glob("*.nam") if path.name.lower() != "mfsim.nam")
    if not name_files:
        # Multi-model run: the name files live one level down, one per model.
        name_files = sorted(
            path
            for path in workspace.glob("*/*.nam")
            if path.name.lower() != "mfsim.nam"
        )
    if not name_files:
        raise FileNotFoundError(f"No model name file found in workspace: {workspace}")
    return name_files[0].stem


# Structural packages every model has, which carry no explorer metadata and so
# are absent from the registry. Kept explicit: the registry describes packages
# with per-package *behaviour* (fields, results, capabilities), and discretization
# or output control has none of that to describe.
_NON_REGISTRY_SUFFIXES = {
    # HFB describes a property of the CONNECTION between two cells rather than
    # either cell's own stress, so it has no registry descriptor -- but a reopened
    # run still has to recognise the file.
    "hfb": "HFB",
    "dis": "DIS",
    "disu": "DISU",
    "disv": "DISV",
    "ic": "IC",
    "mvr": "MVR",
    "npf": "NPF",
    "oc": "OC",
    "sto": "STO",
}

# Boundary/advanced packages come from the registry (plan 4.7.3); the union is
# what file discovery recognizes. A package added to the registry is discovered
# automatically -- it was possible before to add one and silently have its file
# ignored on reopen.
_PACKAGE_SUFFIX_TO_TYPE = dict(
    sorted(
        {
            **{
                spec.file_suffix: spec.name.upper()
                for spec in _PACKAGE_EXPLORER_SPECS.values()
                if spec.file_suffix
            },
            **_NON_REGISTRY_SUFFIXES,
        }.items()
    )
)
_PACKAGE_TYPE_TO_SUFFIX = {value: key for key, value in _PACKAGE_SUFFIX_TO_TYPE.items()}
_CORE_PACKAGE_TYPES = {"DIS", "DISU", "DISV", "OC"}


def _model_namefile(workspace: Path, model_name: str) -> Path | None:
    """The model's name file: from ``mfsim.nam``'s models block, else ``<model>.nam``."""

    try:
        lines = (Path(workspace) / "mfsim.nam").read_text(encoding="utf-8").splitlines()
    except OSError:
        lines = []
    inside = False
    for line in lines:
        stripped = line.split("#", 1)[0].strip()
        lowered = stripped.lower()
        if lowered.startswith("begin") and "models" in lowered:
            inside = True
            continue
        if lowered.startswith("end") and "models" in lowered:
            break
        fields = stripped.split()
        if inside and len(fields) >= 2:
            name = fields[2] if len(fields) >= 3 else Path(fields[1]).stem
            if name.lower() == str(model_name).lower():
                return Path(workspace) / fields[1]
    candidate = _model_dir(workspace, model_name) / f"{model_name}.nam"
    return candidate if candidate.exists() else None


def _repeats(prefix: str, package_type: str) -> bool:
    """Whether MF6 allows several packages of this type in a model -- FloPy's dfn says."""

    cls = getattr(flopy.mf6, f"Modflow{prefix}{package_type.lower()}", None)
    header = next((row for row in getattr(cls, "dfn", ()) if row and row[0] == "header"), ())
    return "multi-package" in header


def _discover_packages(workspace: Path, model_name: str, model_type: str) -> list[tuple[str, str]]:
    """``(type, name)`` for each package in the model's NAME FILE, in its order.

    The name file is the authority. Discovery used to match FILE SUFFIXES against the
    registry's ten packages, and missed every GWT/GWE/PRT package, RCHA/EVTA, MAW,
    CSUB and any custom filename (ledger 228). The type is the ftype without its
    ``6`` -- an array recharge package is ``RCH6`` there, MF6 having one recharge
    package. The name is the entry's own pname; an entry without one takes the name
    MF6 and FloPy give it: ``CHD-1``, ``CHD-2`` for a type that may repeat, the bare
    type otherwise.
    """

    namefile = _model_namefile(workspace, model_name)
    if namefile is None:
        return []
    try:
        lines = namefile.read_text(encoding="utf-8").splitlines()
    except OSError:
        return []
    prefix = str(model_type)[:3].capitalize()
    pairs: list[tuple[str, str]] = []
    counts: dict[str, int] = {}
    inside = False
    for line in lines:
        stripped = line.split("#", 1)[0].strip()
        lowered = stripped.lower()
        if lowered.startswith("begin") and "packages" in lowered:
            inside = True
            continue
        if lowered.startswith("end") and "packages" in lowered:
            break
        if not inside or not stripped:
            continue
        fields = stripped.split()
        package_type = fields[0].upper().removesuffix("6")
        if len(fields) >= 3:
            name = fields[2].upper()
        elif _repeats(prefix, package_type):
            counts[package_type] = counts.get(package_type, 0) + 1
            name = f"{package_type}-{counts[package_type]}"
        else:
            name = package_type
        pairs.append((package_type, name))
    return pairs


def _discover_package_types(workspace: Path, model_name: str) -> list[str]:
    """Return recognized MF6 package types from files in ``workspace``."""

    package_types: list[str] = []
    for path in sorted(_model_dir(workspace, model_name).glob(f"{model_name}.*")):
        suffix = path.suffix.lower().lstrip(".")
        package_type = _PACKAGE_SUFFIX_TO_TYPE.get(suffix)
        if package_type is not None and package_type not in package_types:
            package_types.append(package_type)
    return package_types


def _discover_model_type(workspace: Path, model_name: str) -> str:
    """The MF6 model kind (``gwf6``/``gwt6``/``gwe6``/``prt6``) from ``mfsim.nam``.

    Without this a reopened GWT model reports the class default ``'gwf6'`` and so
    offers ``.hds`` instead of ``.conc`` -- silently the wrong physics, and the
    kind guard cannot help because it believes the lie. ``mfsim.nam``'s ``models``
    block is the authority: ``<type>  <namefile>  <name>``.
    """

    namefile = Path(workspace) / "mfsim.nam"
    if not namefile.exists():
        return "gwf6"
    try:
        lines = namefile.read_text(encoding="utf-8").splitlines()
    except OSError:  # pragma: no cover - unreadable workspace
        return "gwf6"
    in_models = False
    for line in lines:
        stripped = line.strip()
        lowered = stripped.lower()
        if lowered.startswith("begin") and "models" in lowered:
            in_models = True
            continue
        if lowered.startswith("end") and "models" in lowered:
            break
        if not in_models or not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        # <type> <namefile> <name>; the name is optional in MF6, defaulting to
        # the namefile stem.
        if len(fields) >= 3 and fields[2].lower() == str(model_name).lower():
            return fields[0].lower()
        if len(fields) == 2 and Path(fields[1]).stem.lower() == str(model_name).lower():
            return fields[0].lower()
    return "gwf6"


def _discover_grid_type(package_types: list[str]) -> str:
    """Infer the grid type from discovered package types."""

    if "DISU" in package_types:
        return "disu"
    if "DISV" in package_types:
        return "disv"
    if "DIS" in package_types:
        return "dis"
    return "unknown"


def _grid_package_path(workspace: Path, model_name: str, grid_type: str) -> Path | None:
    """Return the discretization-package path for a workspace when present."""

    if grid_type not in {"dis", "disu", "disv"}:
        return None
    path = _model_dir(workspace, model_name) / f"{model_name}.{grid_type}"
    return path if path.exists() else None


def _file_sha256(path: Path) -> str:
    """Return a normalized SHA-256 hash of one text MF6 package file.

    Leading comment lines generated by FloPy contain timestamps and should not
    cause otherwise identical discretization files to compare as different.
    """

    digest = hashlib.sha256()
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            digest.update(stripped.encode("utf-8"))
            digest.update(b"\n")
    return digest.hexdigest()


def _grid_ncpl(modelgrid) -> int:
    """Return the number of cells per layer from an MF6 model grid."""

    return int(np.asarray(modelgrid.xcellcenters, dtype=float).reshape(-1).size)


def _recorded_crs(workspace: Path, model_name: str | None) -> str | None:
    """The CRS a myflopy run recorded for ``model_name``, if it has a manifest.

    ``run.json`` is ``workspace.RUN_MANIFEST_NAME``; the name is spelled here
    because ``workspace`` imports this module.
    """

    try:
        name = model_name or _infer_model_name(workspace)
    except (OSError, ValueError):
        return None
    return recorded_crs(workspace, name)


def load_mf6_run(
    workspace: Path,
    *,
    model_name: str | None = None,
    crs: str | None = None,
    verbosity_level: int = 0,
):
    """Open any MF6 workspace as a lazy, file-backed model for inspection.

    Unlike :func:`load_run` (which reopens a myflopy-managed :class:`Run` from its
    manifest), this loads a *plain* MF6 workspace -- any folder with ``mfsim.nam``,
    however it was created -- into a :class:`LoadedMf6Run` that lazily exposes the
    grid, heads, budgets and packages. Handy for post-processing external models or
    inside a PEST forward run.

    Parameters
    ----------
    workspace : str or pathlib.Path
        Directory containing the MF6 simulation files (``mfsim.nam``). A
        directory that does not exist raises ``FileNotFoundError``.
    model_name : str or None, optional
        Which model's files to open, by its MF6 model name. ``None`` (the
        default) takes the first ``gwf6`` model ``mfsim.nam`` declares; when it
        declares none (or there is no ``mfsim.nam``), the stem of the
        alphabetically first model name file in the workspace, else one level
        down; with no model name file at all it raises ``FileNotFoundError``.
    crs : str or None, optional
        CRS to attach to the grid view (any string pyproj accepts). ``None`` (the
        default) takes the one a myflopy run recorded for this model in its
        ``run.json``, else ``"EPSG:2927"``. A reopened run used to take 2927
        whatever CRS it had run in.
    verbosity_level : int, optional
        ``0`` (the default) is silent; above 0 prints ``[LoadedMf6Run]`` progress
        messages and is passed to FloPy's ``MFSimulation.load`` when the files
        are first read.

    Returns
    -------
    LoadedMf6Run

    Examples
    --------
    >>> model = mf.load_mf6_run("path/to/mf6_run")
    >>> model.load_all()
    >>> heads = model.hds.array()
    """

    if crs is None:
        crs = _recorded_crs(Path(workspace), model_name) or DEFAULT_CRS
    _progress(verbosity_level, f"Opening workspace: {workspace}")
    return LoadedMf6Run(
        workspace=workspace,
        model_name=model_name,
        crs=crs,
        verbosity_level=verbosity_level,
    )


def _simulation_plot_compat(sim, model_list=None, SelPackList=None, **kwargs):
    """Compatibility shim for FloPy simulation plotting."""

    from flopy.plot.plotutil import PlotUtilities

    return PlotUtilities._plot_simulation_helper(
        sim,
        model_list=model_list,
        SelPackList=SelPackList,
        **kwargs,
    )


def patch_simulation_plot(sim):
    """Replace a loaded simulation's ``.plot`` with a compatibility shim, in place.

    FloPy's ``MFSimulation.plot`` can raise on simulations reloaded from disk
    (e.g. unstructured grids or partially populated output). This binds a
    forgiving ``_simulation_plot_compat`` method onto ``sim`` so ``sim.plot(...)``
    works after a reload. Applied automatically by the run loaders; call it
    yourself only on a hand-loaded simulation.

    Parameters
    ----------
    sim : flopy.mf6.MFSimulation
        The FloPy simulation to patch (mutated in place): its ``plot`` attribute
        is replaced by a method that calls FloPy's
        ``PlotUtilities._plot_simulation_helper(sim, model_list=,
        SelPackList=, **kwargs)``.

    Returns
    -------
    flopy.mf6.MFSimulation
        The same ``sim``, for chaining.
    """
    sim.plot = MethodType(_simulation_plot_compat, sim)
    return sim


class LoadedMf6Run(SimulationBase):
    """File-backed exploration object built from MF6 input and output files.

    Construct it through :func:`load_mf6_run`, which also reads the CRS a
    myflopy run recorded; constructing the class directly does not.

    Parameters
    ----------
    workspace : str or pathlib.Path
        Path to a standard MF6 run directory. A directory that does not exist
        raises ``FileNotFoundError``.
    model_name : str or None, optional
        Which model's files to open, by its MF6 model name. ``None`` (the
        default) takes the first ``gwf6`` model ``mfsim.nam`` declares; when it
        declares none (or there is no ``mfsim.nam``), the stem of the
        alphabetically first model name file in the workspace, else one level
        down; with no model name file at all it raises ``FileNotFoundError``. A name
        FloPy cannot find in the simulation raises ``FileNotFoundError`` when the
        files are first loaded.
    crs : str, optional
        CRS assigned to the reconstructed ``VoronoiGridPlus`` view if one is
        built lazily later (any string pyproj accepts). Default
        :data:`~myflopy._defaults.DEFAULT_CRS` (``"EPSG:2927"``); the CRS a run
        recorded in ``run.json`` is read by :func:`load_mf6_run`, not here.
    verbosity_level : int, optional
        ``0`` (the default) is silent; above 0 prints ``[LoadedMf6Run]`` progress
        messages and is passed to FloPy's ``MFSimulation.load`` as its
        ``verbosity_level``.
    """

    def __init__(
        self,
        *,
        workspace: Path,
        model_name: str | None = None,
        crs: str = DEFAULT_CRS,
        verbosity_level: int = 0,
    ):
        """Prepare a lazy, file-backed model view over a completed MF6 run ``workspace``.

        Infers the model name/grid type and discovers package types from filenames
        without loading FloPy; the simulation, grid, and outputs load on first use.
        """

        workspace = Path(workspace)
        if not workspace.exists():
            raise FileNotFoundError(f"Run workspace does not exist: {workspace}")

        model_name = _infer_model_name(workspace) if model_name is None else model_name
        pairs = _discover_packages(workspace, model_name, _discover_model_type(workspace, model_name))
        if pairs:
            available_package_types = list(dict.fromkeys(package_type for package_type, _ in pairs))
        else:  # no name file to read: the file suffixes are all there is
            available_package_types = _discover_package_types(workspace, model_name)
            pairs = [(package_type, package_type) for package_type in available_package_types]
        self._discovered_pairs = tuple(pairs)
        grid_type = _discover_grid_type(available_package_types)

        self.name = model_name
        self.model_type = _discover_model_type(workspace, model_name)
        self._sim = None
        self._gwf = None
        self.ims = None
        self._crs = crs
        self._verbosity_level = verbosity_level
        self._vor = None
        self.nper = None
        self.nlay = None
        self.num_steps = None
        self.per_len = None
        self._obs = None
        self._times = None
        self.model_output_folder_path = workspace
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
        self._budget_kstpkper = None
        self._idomain = None
        self._grid_type_override = grid_type if grid_type != "unknown" else None
        self._targets = None
        self._particle_tracking = None
        self._parallel = None
        self.regions = RegionRegistry(self)
        self.source = "mf6_files"
        self._model_name = model_name
        self._available_package_types = tuple(available_package_types)
        self._available_package_types_lower = {value.lower() for value in available_package_types}
        self._loaded_package_types: set[str] = set()
        self._fully_loaded = False
        self._core_loaded = False
        self._shared_vor_source = None
        _progress(
            verbosity_level,
            (
                # The DISCOVERED grid type, not the `grid_type` property: the
                # property falls back to inspecting `self.gwf`, so building this
                # message loaded the entire FloPy simulation -- in a constructor
                # whose whole promise is that it does not, and whose f-string is
                # evaluated even at verbosity 0 where nothing is printed.
                f"Prepared lazy file-backed model '{self.name}' "
                f"(grid_type={self._grid_type_override or 'unknown'}, "
                f"packages={list(self.package_names)}). "
                "Grid, outputs, and package definitions will load only when needed."
            ),
        )

    @property
    def sim(self):
        """Lazily loaded FloPy simulation object."""

        self._ensure_core_loaded()
        return self._sim

    @sim.setter
    def sim(self, value):
        """Set the cached FloPy simulation object."""

        self._sim = value

    @property
    def gwf(self):
        """Lazily loaded FloPy groundwater-flow model object."""

        self._ensure_core_loaded()
        return self._gwf

    @gwf.setter
    def gwf(self, value):
        """Set the cached FloPy groundwater-flow model object."""

        self._gwf = value

    @property
    def package_names(self) -> list[str]:
        """The packages' own NAMES, as MF6 has them -- ``CHD_0``, ``RCHA_0``, ``NPF``.

        Read from the model's name file, so it stays cheap before any FloPy load and
        agrees with the same model live. It returned TYPES until 2026-09-11, so a
        reopened run and the live model disagreed (ledger 228); those are
        :attr:`package_types`.
        """

        return sorted(name for _, name in self._discovered_pairs)

    @property
    def package_types(self) -> list[str]:
        """Each package TYPE once -- ``CHD``, ``NPF``, ``RCH`` (array recharge is RCH)."""

        return sorted(set(self._available_package_types))

    def _package_pairs(self) -> list[tuple[str, str]]:
        """``(name, type)`` for each of :attr:`package_names`, in its order."""

        types = {name: package_type for package_type, name in self._discovered_pairs}
        return [(str(name).upper(), types.get(str(name).upper(), str(name).upper())) for name in self.package_names]

    def _summary_ncpl(self) -> int | None:
        """Return cached ``ncpl`` for a lightweight summary when possible."""

        if self._vor is not None:
            return int(self._vor.ncpl)
        if self._gwf is not None:
            return _grid_ncpl(self._gwf.modelgrid)
        return None

    def _summary_nlay(self) -> int | None:
        """Return cached layer count for :meth:`summary`."""

        return self.nlay

    def _summary_nper(self) -> int | None:
        """Return cached stress-period count for :meth:`summary`."""

        return self.nper

    @property
    def modelgrid(self):
        """Underlying FloPy modelgrid object, loaded only when needed."""

        self._ensure_core_loaded()
        return self._gwf.modelgrid

    def _normalize_package_names(self, package_names: str | list[str] | tuple[str, ...] | set[str]) -> list[str]:
        """Normalize package names into lower-case MF6 package identifiers."""

        if isinstance(package_names, (list, tuple, set)):
            return [str(name).lower() for name in package_names]
        return [str(package_names).lower()]

    def _load_simulation(
        self,
        *,
        load_only: list[str] | None,
        mark_full: bool,
        reason: str,
    ):
        """Load the FloPy simulation with the requested package subset."""

        # The caller's level, never floored: `load_mf6_run(verbosity_level=0)` is the
        # default and asks for silence, and a floor of 1 printed FloPy's load chatter
        # the first time a reopened run touched a package.
        display_verbosity = self._verbosity_level

        if load_only is None:
            _progress(display_verbosity, f"Loading full MF6 simulation from: {self.workspace}")
        else:
            _progress(
                display_verbosity,
                f"Loading MF6 packages for {reason}: {', '.join(package.upper() for package in load_only)}",
            )
        sim = flopy.mf6.MFSimulation.load(
            sim_ws=str(self.workspace),
            load_only=load_only,
            lazy_io=True,
            verbosity_level=display_verbosity,
        )
        # FloPy reads a LAK file's IMPLICIT without a word and drops it, so a rewrite
        # (run_simulation below writes first) would lose it (ledger 326).
        keep_implicit(sim)
        patch_simulation_plot(sim)
        _progress(display_verbosity, f"Locating GWF model '{self._model_name}'")
        gwf = sim.get_model(self._model_name)
        if gwf is None:
            raise FileNotFoundError(
                f"Could not load GWF model '{self._model_name}' from workspace: {self.workspace}"
            )

        self._sim = sim
        self._gwf = gwf
        self.nper = int(sim.tdis.nper.data)
        self.nlay = int(np.asarray(gwf.modelgrid.botm).reshape(-1).size / _grid_ncpl(gwf.modelgrid))
        if self._idomain is None:
            self._idomain = idomain_from_flopy(gwf)
        # TYPES, not names: `get_package_list()` answers `CHD_0` where a package is
        # named for FloPy's default, and every check below asks for `CHD`. Recording
        # names made each read of such a package reload the whole simulation.
        self._loaded_package_types = {_package_type_name(package.package_type) for package in gwf.packagelist}
        self._core_loaded = True
        self._fully_loaded = mark_full or self._loaded_package_types >= set(self.package_types)

    def _ensure_core_loaded(self):
        """Ensure core grid/output infrastructure is available."""

        if self._core_loaded and self._gwf is not None:
            return
        available = [name.lower() for name in self.package_types if name in _CORE_PACKAGE_TYPES]
        # The grid type DISCOVERED from filenames, never the `grid_type`
        # property: that property falls back to inspecting `self.gwf`, and
        # `self.gwf` calls this method -- an infinite recursion the moment
        # discovery comes up empty. Keeping the loader on the cheap, file-derived
        # answer makes the cycle structurally impossible rather than merely
        # unreached.
        discovered = getattr(self, "_grid_type_override", None)
        if discovered in {"dis", "disu", "disv"} and discovered not in available:
            available.insert(0, discovered)
        if "oc" not in available and "OC" in self.package_types:
            available.append("oc")
        self._load_simulation(load_only=available or None, mark_full=False, reason="core access")

    def _ensure_packages_loaded(self, package_names: list[str]):
        """Ensure one or more named packages are available on ``self.gwf``."""

        normalized = self._normalize_package_names(package_names)
        if all(name.upper() in self._loaded_package_types for name in normalized):
            return
        # A full load holds every package in the workspace, so one still missing is
        # ABSENT, not unloaded. Reloading for it replaced the in-memory simulation,
        # and `run_simulation` then wrote the files back without the caller's edits.
        if self._fully_loaded:
            return

        # One FULL load, whatever was asked for. A package handed out is one the
        # caller may edit, and `run_simulation` writes the full simulation -- so a
        # package served from a partial load was dropped, edits and all, by the
        # full load writing needs (an edited OC, measured). Growing `load_only` a
        # package at a time also cost more than one full load.
        self._load_simulation(load_only=None, mark_full=True, reason="package access")

    @property
    def vor(self) -> VoronoiGridPlus:
        """Build and cache a ``VoronoiGridPlus`` view on first access."""

        if self._shared_vor_source is not None:
            return self._shared_vor_source.vor
        if self._vor is None:
            self._ensure_core_loaded()
            _progress(self._verbosity_level, "Building Voronoi grid view from the MF6 modelgrid")
            vor = grid_view_from_modelgrid(self._gwf.modelgrid, crs=self._crs, name=self.name)
            self._vor = vor
            _progress(self._verbosity_level, f"Voronoi grid ready (ncpl={self._vor.ncpl})")
        return self._vor

    def grid_signature(self) -> str | None:
        """Return a stable signature for the run's discretization package.

        For file-backed runs this uses the exact ``.dis/.disv/.disu`` file
        bytes, which allows higher-level helpers such as ``ModelGroup`` to
        confirm grids match without forcing full grid reconstruction first.
        """

        grid_path = _grid_package_path(self.workspace, self.name, self.grid_type)
        if grid_path is None:
            return None
        return _file_sha256(grid_path)

    def package(self, package_names: str | list[str]):
        """Load and return one or more MF6 package objects on demand.

        The first request for a package not yet loaded loads the WHOLE
        simulation (once), so a package handed out is one ``run_simulation``
        will write back with any edits.

        Parameters
        ----------
        package_names : str or list, tuple or set of str
            A package TYPE or NAME (e.g. ``"chd"``, ``"npf"``, ``"lak_2"``), or
            several; lower-cased and looked up as
            :meth:`SimulationBase.package <myflopy.modflow.mf6.simulation.base.SimulationBase.package>`
            does -- by TYPE first, then a package's exact NAME, on the model and
            then on the simulation (so ``"tdis"`` and ``"ims"`` resolve too).

        Returns
        -------
        flopy.mf6.mfpackage.MFPackage or list or dict
            For one string, the package, or a list of every package of the type
            when a TYPE matches several (ledger 352). For two or more, a dict of
            lower-cased string to that result.

        Raises
        ------
        AttributeError
            When a string matches no package type or name in the model or the
            simulation.
        """

        normalized = self._normalize_package_names(package_names)
        if len(normalized) > 1:
            self._ensure_packages_loaded(normalized)
            return {name: self.package(name) for name in normalized}

        self._ensure_packages_loaded(normalized)
        return super().package(normalized[0])

    def _packages_of_type(self, package_type: str) -> list:
        """Every FloPy package of TYPE ``package_type``, loading the run's packages first."""

        self._ensure_packages_loaded([package_type])
        return super()._packages_of_type(package_type)

    def _get_budget_reader(self):
        """Return the binary cell-budget reader without loading FloPy core packages."""

        if self._bud is None:
            # MF6 always writes double-precision binaries. FloPy's default
            # precision="auto" tries single first and can misread a double file:
            # on Windows the bogus record length makes file.seek raise
            # OSError (Errno 22) -- which escapes FloPy's BudgetIndexError/EOFError
            # guard, so the double-precision fallback never runs. Read as double.
            self._bud = flopy.utils.CellBudgetFile(
                _declared_output_path(self.workspace, self.name, "budget"),
                precision="double",
            )
        return self._bud

    def _get_budget_kstpkper(self) -> list[tuple[int, int]]:
        """Return budget ``(kstp, kper)`` pairs directly from the model budget file.

        Cached in its OWN slot. It used to share ``_kstpkper`` with
        :func:`~myflopy.modflow.mf6.simulation.accessors.get_kstpkper`, which
        caches the HEAD file's clock -- two different clocks, one slot, so
        whichever accessor was touched first answered for both for the life of
        the object. Measured on a model with 6 head records and 2 budget records:
        head-first made this return the 6 head pairs, budget-first made
        ``model.kstpkper`` return only the 2 budget times (ledger 301).
        """

        if self._budget_kstpkper is None:
            self._budget_kstpkper = [
                tuple(int(value) for value in pair)
                for pair in self._get_budget_reader().get_kstpkper()
            ]
        return self._budget_kstpkper

    def load_all(self):
        """Fully load the MF6 simulation and all discovered packages."""

        if not self._fully_loaded:
            self._load_simulation(load_only=None, mark_full=True, reason="full load")
        return self

    def run_simulation(self):
        """Rewrite the loaded simulation files and run MF6 from the workspace."""

        self.load_all()
        self._sim.write_simulation()
        warn_if_implicit_unsupported(self._sim)
        success, buff = self._sim.run_simulation(silent=False, report=True)
        return success, buff
