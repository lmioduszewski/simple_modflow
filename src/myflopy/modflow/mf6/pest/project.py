"""Project orchestration for the first ``myflopy`` PEST slice."""

from __future__ import annotations

import contextlib
import io
import json
import re
import sys
import warnings
from pathlib import Path

from myflopy.modflow.mf6.observations import HeadTargets
from myflopy.modflow.mf6.pest.observations import (
    prepare_drn_flow_observations,
    finalize_observations,
    prepare_head_target_observations,
    prepare_lake_stage_observations,
    prepare_sfr_flow_observations,
    prepare_sfr_stage_observations,
    _observation_name,
)
from myflopy.modflow.mf6.pest.parameters import prepare_parameter_spec, register_parameter_spec
from myflopy.modflow.mf6.pest.native_parameters import NativeParameterSpec, add_native_parameter
from myflopy.modflow.mf6.pest.pilot_points import (
    add_pilot_point_parameter,
    register_pilot_point_parameters,
)
from myflopy.modflow.mf6.pest.summary import PestSettings
from myflopy.modflow.mf6.pest.specs import (
    HeadTargetObservationSpec,
    KPilotPointParameter,
    DrainElevationParameter,
    DrainConductanceParameter,
    DrnFlowObservationSpec,
    LakeStageObservationSpec,
    SfrFlowObservationSpec,
    SfrStageObservationSpec,
)

METADATA_FILENAME = "myflopy_pest_metadata.json"


def _import_pyemu():
    """Import pyEMU lazily with a workflow-oriented error message."""

    try:
        import pyemu
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "pyemu is required for PestProject workflows. Install pyemu in the "
            "active environment before calling PestProject.build_pst()."
        ) from exc
    return pyemu


def _pyemu_warning_class(pyemu_module):
    """Return the installed ``PyemuWarning`` class if available."""

    try:
        return pyemu_module.pyemu_warnings.PyemuWarning
    except AttributeError:
        return Warning


class PestProject:
    """Build a pyEMU/PEST++ calibration workspace from a ``myflopy`` model.

    ``PestProject`` is the single front door for calibration. The recommended
    (modern) workflow compiles a few readable declarations straight to native
    ``pyemu.utils.PstFrom`` and lets pyEMU drive the forward run::

        cal = PestProject(model, "calib", workspace="calib_template",
                          start_datetime="2020-01-01")
        cal.parameterize("k",        style="constant", bounds=(0.2, 5), physical=(1e-3, 100))
        cal.parameterize("recharge", bounds=(0.5, 1.5))
        cal.parameterize("ghb.cond", bounds=(0.1, 10))
        cal.observe(head_targets)          # history-matching targets
        cal.forecast(prediction_targets)   # predictions of interest (zero weight)
        print(cal.settings())              # review the resolved configuration
        pst = cal.build("calib.pst")       # -> native .pst + forward_run.py

    Use :meth:`parameterize`, :meth:`observe`, :meth:`forecast`, :meth:`build`
    and :meth:`settings` for new work. The legacy spec-object path
    (:meth:`add_parameter`/:meth:`add_observation`/:meth:`build_pst` with
    :class:`KPilotPointParameter` etc.) is still supported for Voronoi pilot
    points and named-series (lake/SFR/DRN) observations until those land on the
    native path.

    Parameters
    ----------
    model
        A built ``myflopy`` model (``SimulationBase`` or compatible) exposing
        ``.sim`` (the FloPy simulation), ``.name``, and ``.workspace``.
    name
        Project name; used as the default ``.pst`` filename stem.
    workspace
        Directory for the generated PEST template workspace. It is created
        (and cleared) when the control file is built.
    start_datetime
        Simulation start date (e.g. ``"2020-01-01"``). Required by pyEMU to
        place time-varying parameters and observations on the time axis.
    spatial_reference
        Optional cell spatial reference; only needed for spatial array
        parameters (pilot points / grid geostatistics).
    zero_based
        Whether index columns in list files are zero-based (default ``False``,
        matching pyEMU's convention for the indices pyEMU itself writes).
    longnames
        Use long PEST parameter/observation names (default ``True``).
    """

    def __init__(
        self,
        model,
        name: str,
        workspace: str | Path,
        start_datetime: str,
        spatial_reference=None,
        zero_based: bool = False,
        longnames: bool = True,
    ):
        """Create a calibration project around an existing model workspace.

        See the class docstring for the full parameter reference and workflow.
        """

        self.model = model
        self.name = str(name)
        self.original_workspace = Path(model.workspace)
        self.template_workspace = Path(workspace)
        self.start_datetime = start_datetime
        self.spatial_reference = spatial_reference
        self.zero_based = bool(zero_based)
        self.longnames = bool(longnames)
        self._parameter_specs: list = []
        self._observation_specs: list = []
        self._forecast_specs: list = []
        self._native_parameter_specs: list[NativeParameterSpec] = []
        self._native_parameter_frames: dict[str, object] = {}
        self._pilot_point_frames: dict[str, list] = {}
        self._capture_field_specs: dict[str, dict] = {}
        self._prepared_observations: list[dict] = []
        self._prepared_parameters: dict[str, dict] = {}
        self._parameter_frames: dict[str, object] = {}
        self._registered_geostructs: dict[str, object] = {}
        self._forward_run_config: dict = {}
        self.pyemu = None
        self.pf = None
        self.pst = None

    def add_parameter(self, spec):
        """Register a parameter specification for the project."""

        self._parameter_specs.append(spec)
        return spec

    def add_observation(self, spec):
        """Register an observation specification for the project."""

        self._observation_specs.append(spec)
        return spec

    # -- modern declarative facade ----------------------------------------
    #
    # These methods compile directly to native ``pyemu.utils.PstFrom`` calls
    # (see native_parameters.py) rather than the hand-rolled template path used
    # by ``add_parameter``/``build_pst``. Prefer ``parameterize``/``observe``/
    # ``forecast``/``build`` for new calibrations.

    def parameterize(
        self,
        target: str,
        *,
        style: str | None = None,
        bounds: tuple[float, float] = (0.5, 2.0),
        physical: tuple[float, float] | None = None,
        transform: str = "log",
        additive: bool | None = None,
        zones=None,
        layers=None,
        pp_space: int | None = None,
        pp_points=None,
        correlation: float | None = None,
        temporal: float | None = None,
        name: str | None = None,
        capture: bool = False,
        **extra,
    ) -> NativeParameterSpec:
        """Declare a calibration parameter that compiles to ``pf.add_parameters``.

        Each call expands, at :meth:`build` time, into the correct
        ``pyemu.utils.PstFrom.add_parameters`` invocation against the external
        MODFLOW input file(s) for ``target``. Call it once per property you want
        to adjust; call it multiple times for the same target to stack styles
        (e.g. a broad ``"constant"`` plus a fine ``"grid"``).

        Parameters
        ----------
        target
            Friendly name of the model property to calibrate. Supported targets
            (with aliases): ``"k"`` (``"npf.k"``/``"hk"``/``"kh"``),
            ``"k33"`` (``"kv"``), ``"recharge"`` (``"rch"``),
            ``"chd"``, ``"ghb.cond"``/``"ghb.bhead"`` (``"ghb"``),
            ``"drn.cond"``/``"drn.elev"`` (``"drn"``), and ``"wel"`` (``"pumping"``).
        style
            Spatial parameterization style:

            - ``"constant"`` -- a single multiplier for the whole target;
            - ``"zone"`` -- one multiplier per ``zones`` value (pass ``zones=``);
            - ``"grid"`` -- one multiplier per list entry / cell;
            - ``"pilotpoints"`` -- geostatistical pilot points.

            ``None`` uses the recipe default (``"constant"``). ``"grid"`` on
            *array* targets (K, K33) builds one geostatistically-correlated
            multiplier per Voronoi/DISV cell, drawn against the model's spatial
            reference (the modelgrid is used automatically); pass ``correlation=``
            for the variogram range and ``layers=`` to restrict which layers.
            ``"grid"`` also works for list targets. ``"pilotpoints"`` on array
            targets is not wired for Voronoi grids yet.
        bounds
            ``(lower, upper)`` bounds on the adjustable value. For a multiplier
            these are factors (e.g. ``(0.2, 5)`` = 5x down to 5x up); for an
            additive parameter they are offsets in model units.
        physical
            ``(ult_lbound, ult_ubound)`` -- hard limits clamped onto the *final*
            model value after all multipliers are applied. Strongly recommended
            for multipliers so calibration cannot push K (etc.) to nonphysical
            values. ``None`` leaves the input unclamped.
        transform
            ``"log"`` (default, recommended for strictly positive properties
            like K and recharge) or ``"none"``. Forced to ``"none"`` for
            additive parameters.
        additive
            Apply the parameter as an additive offset instead of a multiplier.
            ``None`` uses the recipe default (drain elevation is additive;
            everything else is multiplicative).
        zones
            Zone array for ``style="zone"`` (and to mask inactive cells). Values
            map one parameter per distinct zone.
        layers
            For multi-layer *array* targets, restrict the parameter to these
            zero-based model layers (e.g. ``layers=[0, 1]`` for the unconfined
            aquifers); other layers keep their input value. Defaults to all
            layers.
        correlation
            Variogram range (model length units) for ``grid``/``pilotpoints``
            spatial correlation. Ignored for ``constant``/``zone``.
        temporal
            Temporal correlation range (days) for time-varying list packages.
            Recorded now; wired in a later phase.
        name
            Parameter group / name base. Defaults to a slug of ``target``
            (e.g. ``"ghbcond"``). Useful when stacking multiple calls.
        **extra
            Advanced keyword arguments passed straight through to
            ``pf.add_parameters`` for cases the facade does not cover.

        Returns
        -------
        NativeParameterSpec
            The recorded specification (also stored on the project). Inspect the
            full set at any time with :meth:`settings`.

        Examples
        --------
        >>> cal.parameterize("k", style="constant", bounds=(0.2, 5), physical=(1e-3, 100))
        >>> cal.parameterize("recharge", bounds=(0.5, 1.5), physical=(0, 1e-2))
        >>> cal.parameterize("ghb.cond", bounds=(0.1, 10))
        >>> cal.parameterize("drn.elev", bounds=(-2, 2))   # additive offsets
        """

        spec = NativeParameterSpec(
            target=target,
            style=style,
            bounds=bounds,
            physical=physical,
            transform=transform,
            additive=additive,
            zones=zones,
            layers=tuple(layers) if layers is not None else None,
            pp_space=pp_space,
            pp_points=pp_points,
            correlation=correlation,
            temporal=temporal,
            name=name,
            capture=capture,
            extra=extra,
        )
        self._native_parameter_specs.append(spec)
        return spec

    def observe(self, targets, *, prefix: str | None = None):
        """Register history-matching observations from a ``myflopy`` target set.

        Observations are the measured data calibration tries to reproduce. At
        :meth:`build` time the targets are registered with pyEMU
        (``pf.add_observations``), their measured values and weights are written
        into the control file, and a post-processor is added to the forward run
        so simulated equivalents are regenerated on every model run.

        Parameters
        ----------
        targets
            A :class:`HeadTargets` object (locations + measured values), or a
            pre-built observation spec (e.g. :class:`HeadTargetObservationSpec`).
            Other target types (lake/SFR/DRN) are currently supported only via
            :meth:`add_parameter`-style specs on the legacy :meth:`build_pst`.
        prefix
            Observation-name prefix (default ``"hds"``). Use distinct prefixes
            when registering more than one head-target set.

        Returns
        -------
        The registered observation spec.
        """

        spec = self._coerce_observation_spec(targets, prefix=prefix, default_prefix="hds")
        self._observation_specs.append(spec)
        return spec

    def forecast(self, targets, *, prefix: str | None = None):
        """Register a prediction of interest as a zero-weight forecast.

        Forecasts use the *same* target objects as :meth:`observe`, but they are
        never history-matched: at :meth:`build` time their weight is set to zero
        and their observation names are recorded in
        ``pst.pestpp_options["forecasts"]``. The downstream uncertainty tools
        (IES ensembles, FOSM) report posterior uncertainty for exactly these
        quantities, so declare here whatever model output you ultimately care
        about predicting (a future head, a stream flux, a seepage rate).

        Parameters
        ----------
        targets
            A :class:`HeadTargets` (or observation spec) describing *where* and
            *when* the prediction is taken. Any values supplied are placeholders
            (e.g. an expected value); the weight is forced to zero regardless.
        prefix
            Observation-name prefix. Defaults to ``"fore1"``, ``"fore2"``, ...
            so forecast names never collide with calibration observations.

        Returns
        -------
        The registered forecast spec.
        """

        default = prefix or f"fore{len(self._forecast_specs) + 1}"
        spec = self._coerce_observation_spec(targets, prefix=default, default_prefix=default)
        self._forecast_specs.append(spec)
        return spec

    def _coerce_observation_spec(self, targets, *, prefix, default_prefix):
        """Turn a target object (or spec) into a registered observation spec."""

        observation_spec_types = (
            HeadTargetObservationSpec,
            LakeStageObservationSpec,
            SfrStageObservationSpec,
            SfrFlowObservationSpec,
            DrnFlowObservationSpec,
        )
        if isinstance(targets, observation_spec_types):
            return targets
        if isinstance(targets, HeadTargets):
            return HeadTargetObservationSpec(targets=targets, prefix=prefix or default_prefix)
        raise TypeError(
            "observe()/forecast() accept a HeadTargets object or a pre-built "
            f"observation spec; got {type(targets).__name__}. Other target types "
            "are supported through add_observation() with the legacy build."
        )

    def _ensure_original_workspace(self):
        """Ensure MF6 input files exist before creating the PEST template."""

        self.original_workspace.mkdir(parents=True, exist_ok=True)
        self.model.sim.write_simulation(silent=True)

    def _build_pstfrom(self):
        """Instantiate the underlying ``pyemu.utils.PstFrom`` object."""

        self.pyemu = _import_pyemu()
        from pyemu.utils.pst_from import PstFrom

        # Spatial reference enables geostatistical (grid / pilot-point)
        # parameters. Default to the model's (DISV/Voronoi) modelgrid -- pyEMU
        # reads cell centroids from it to build spatially-correlated parameters.
        spatial_reference = self.spatial_reference
        if spatial_reference is None:
            spatial_reference = getattr(self.model.gwf, "modelgrid", None)

        self.pf = PstFrom(
            original_d=self.original_workspace,
            new_d=self.template_workspace,
            remove_existing=True,
            longnames=self.longnames,
            spatial_reference=spatial_reference,
            zero_based=self.zero_based,
            start_datetime=self.start_datetime,
            echo=False,
        )

    @contextlib.contextmanager
    def _quiet_pyemu_context(self):
        """Suppress known harmless pyEMU chatter during workspace assembly."""

        pyemu_warning = _pyemu_warning_class(self.pyemu)
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=".*not being called directly.*",
                category=pyemu_warning,
            )
            warnings.filterwarnings(
                "ignore",
                message="no adjustable pars",
                category=pyemu_warning,
            )
            warnings.filterwarnings(
                "ignore",
                message="Setting an item of incompatible dtype is deprecated.*",
                category=FutureWarning,
            )
            with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(stderr_buffer):
                yield

    def _prepare_one_observation_spec(self, spec):
        """Create simulated files and register one observation spec with pyEMU."""

        if isinstance(spec, HeadTargetObservationSpec):
            return prepare_head_target_observations(self, spec)
        if isinstance(spec, LakeStageObservationSpec):
            return prepare_lake_stage_observations(self, spec)
        if isinstance(spec, SfrStageObservationSpec):
            return prepare_sfr_stage_observations(self, spec)
        if isinstance(spec, SfrFlowObservationSpec):
            return prepare_sfr_flow_observations(self, spec)
        if isinstance(spec, DrnFlowObservationSpec):
            return prepare_drn_flow_observations(self, spec)
        raise TypeError(f"Unsupported observation spec type: {type(spec).__name__}")

    def _prepare_observation_specs(self):
        """Create simulated observation files and register them with pyEMU."""

        prepared: list[dict] = []
        for spec in self._observation_specs:
            prepared.append(self._prepare_one_observation_spec(spec))
        for spec in self._forecast_specs:
            item = self._prepare_one_observation_spec(spec)
            item["is_forecast"] = True
            prepared.append(item)
        self._prepared_observations = prepared

    def _prepare_parameter_specs(self):
        """Create template/support files for parameter specs before ``build_pst()``."""

        self._prepared_parameters = {}
        for spec in self._parameter_specs:
            if not isinstance(
                spec,
                (KPilotPointParameter, DrainElevationParameter, DrainConductanceParameter),
            ):
                raise TypeError(f"Unsupported parameter spec type: {type(spec).__name__}")
            prepare_parameter_spec(self, spec)

    def _build_forward_run_config(self):
        """Build the JSON configuration consumed by the injected forward run."""

        head_target_outputs = [
            item["forward_run_config"]
            for item in self._prepared_observations
            if "forward_run_config" in item
        ]
        named_series_outputs = [
            item["named_series_forward_run_config"]
            for item in self._prepared_observations
            if "named_series_forward_run_config" in item
        ]
        k_specs = [
            prepared["config"]
            for prepared in self._prepared_parameters.values()
            if prepared["config"]["kind"] == "k_pilotpoints"
        ]
        drain_specs = [
            prepared["config"]
            for prepared in self._prepared_parameters.values()
            if prepared["config"]["kind"] in {"drn_elev", "drn_cond"}
        ]
        self._forward_run_config = {
            "model_name": self.model.name,
            "exe_name": getattr(self.model.sim, "exe_name", "mf6"),
            "k_specs": k_specs,
            "drain_specs": drain_specs,
            "head_target_outputs": head_target_outputs,
            "named_series_outputs": named_series_outputs,
        }
        config_path = self.template_workspace / "pest_forward_config.json"
        config_path.write_text(json.dumps(self._forward_run_config, indent=2), encoding="utf-8")
        return config_path

    def _write_project_metadata(self, *, filename: str | Path | None = None):
        """Persist reopen-friendly run metadata inside the generated workspace."""

        observation_sets = [
            item["metadata"]
            for item in self._prepared_observations
            if "metadata" in item and item["metadata"] is not None
        ]
        parameter_sets = []
        for spec in self._parameter_specs:
            prepared = self._prepared_parameters.get(spec.name, {})
            config = prepared.get("config", {})
            parameter_sets.append(
                {
                    "name": spec.name,
                    "type": type(spec).__name__,
                    "kind": config.get("kind"),
                    "transform": getattr(spec, "transform", None),
                    "bounds_mode": getattr(spec, "bounds_mode", None),
                    "bounds": list(spec.bounds) if getattr(spec, "bounds", None) is not None else None,
                    "parameter_space": getattr(spec, "parameter_space", None),
                    "parameter_style": getattr(spec, "parameter_style", None),
                    "files": {key: value for key, value in config.items() if key.endswith("_csv")},
                }
            )

        metadata = {
            "version": 1,
            "project_name": self.name,
            "model_name": self.model.name,
            "original_workspace": str(self.original_workspace.resolve()),
            "forward_run_config_file": "pest_forward_config.json",
            "pst_file": Path(filename).name if filename is not None else None,
            "observation_sets": observation_sets,
            "parameter_sets": parameter_sets,
            "capture_fields": list(self._capture_field_specs.values()),
        }
        metadata_path = self.template_workspace / METADATA_FILENAME
        metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        return metadata_path

    def _attach_forward_run(self):
        """Inject the reusable forward-run helper into pyEMU's ``forward_run.py``."""

        helper_path = Path(__file__).with_name("forward_run.py")
        helper_calls = [
            "_expand_cellid_columns(frame=None)",
            "_build_drn_frame(package=None)",
            "_rebuild_drn_stress_period_data(frame=None)",
            "_load_parameter_values(path=None)",
            "_apply_drain_specs(gwf=None, specs=None)",
            "_interpolate_idw(x=None, y=None, px=None, py=None, values=None)",
            "_apply_k_specs(gwf=None, specs=None)",
            "_write_head_target_csv(model_name=None, mapping_csv=None, output_csv=None)",
            "_read_saved_locations(path=None)",
            "_load_model_for_named_series(sim_ws='.')",
            "_write_named_series_target_csv(model=None, kind=None, locations_file=None, output_csv=None)",
            "_write_simulation_with_retry(sim=None)",
        ]
        for call in helper_calls:
            self.pf.add_py_function(str(helper_path), call, is_pre_cmd=None)
        self.pf.add_py_function(
            str(helper_path),
            "apply_pest_forward_run(config_path='pest_forward_config.json')",
            is_pre_cmd=True,
        )

    def _finalize_forward_run_script(self):
        """Trim noisy default pyEMU helpers from the generated script."""

        forward_run_path = self.template_workspace / "forward_run.py"
        text = forward_run_path.read_text(encoding="utf-8")
        filtered = []
        for line in text.splitlines():
            if "apply_list_and_array_pars(arr_par_file='mult2model_info.csv'" in line:
                continue
            if "print(r'error removing tmp file:" in line:
                filtered.append("       pass")
                continue
            filtered.append(line)
        forward_run_path.write_text("\n".join(filtered) + "\n", encoding="utf-8")

    def _register_parameter_specs(self):
        """Add template-file parameters to the built control file."""

        for spec in self._parameter_specs:
            if not isinstance(
                spec,
                (KPilotPointParameter, DrainElevationParameter, DrainConductanceParameter),
            ):
                raise TypeError(f"Unsupported parameter spec type: {type(spec).__name__}")
            register_parameter_spec(self, spec)

    def build_pst(self, filename: str | Path | None = None):
        """Build and return a ``pyemu.Pst`` control object."""

        self._ensure_original_workspace()
        self._build_pstfrom()
        self._prepare_observation_specs()
        self._prepare_parameter_specs()
        self._build_forward_run_config()
        self._write_project_metadata(filename=filename)
        with self._quiet_pyemu_context():
            self._attach_forward_run()
            self.pst = self.pf.build_pst(filename=filename)
        self.pst.model_command = [f'"{sys.executable}" forward_run.py']
        self._finalize_forward_run_script()
        with self._quiet_pyemu_context():
            self._register_parameter_specs()
        finalize_observations(self, self._prepared_observations)
        if filename is not None:
            with self._quiet_pyemu_context():
                self.pst.write(self.template_workspace / Path(filename).name)
        return self.pst

    def write(self, filename: str | Path | None = None):
        """Write the built control file to disk."""

        if self.pst is None:
            raise ValueError("Call build_pst() before write().")
        target = self.template_workspace / (Path(filename).name if filename else f"{self.name}.pst")
        self.pst.write(target)
        return target

    def draw_prior(self, num_reals: int = 1000, use_specsim: bool = True):
        """Draw a prior parameter ensemble using the underlying ``PstFrom``."""

        if self.pf is None:
            raise ValueError("Call build_pst() or build() before draw_prior().")
        return self.pf.draw(num_reals=num_reals, use_specsim=use_specsim)

    # -- native build path ------------------------------------------------

    def _ensure_external_model(self):
        """Write MF6 inputs as external array/list files for native PstFrom."""

        self.original_workspace.mkdir(parents=True, exist_ok=True)
        self.model.sim.set_all_data_external()
        self.model.sim.write_simulation(silent=True)

    def _resolve_exe(self) -> str:
        """Return a forward-run-usable MF6 command (absolute path when found)."""

        exe = getattr(self.model.sim, "exe_name", "mf6") or "mf6"
        exe_path = Path(exe)
        if not exe_path.is_absolute():
            candidate = (Path(self.model.workspace) / exe).resolve()
            if candidate.exists():
                return str(candidate)
        if exe_path.exists():
            return str(exe_path.resolve())
        return exe

    def _geostruct_for(self, spec: NativeParameterSpec):
        """Build a pyEMU geostruct from a parameter spec's correlation range."""

        pyemu = self.pyemu or _import_pyemu()
        vario = pyemu.geostats.ExpVario(contribution=1.0, a=float(spec.correlation))
        return pyemu.geostats.GeoStruct(variograms=[vario], transform=spec.resolved_transform)

    def _attach_native_observation_postprocessors(self):
        """Add post-model functions that regenerate simulated observation CSVs."""

        helper_path = Path(__file__).with_name("forward_run.py")
        for item in self._prepared_observations:
            config = item.get("forward_run_config")
            if config is None:
                raise NotImplementedError(
                    f"The native build() currently wires head-target observations "
                    f"only; observation set {item.get('prefix')!r} is not yet "
                    "supported. Use the legacy build_pst() for it, or wait for the "
                    "named-series observation phase."
                )
            self.pf.add_py_function(
                str(helper_path),
                "_write_head_target_csv("
                f"model_name='{self.model.name}', "
                f"mapping_csv='{config['mapping_csv']}', "
                f"output_csv='{config['output_csv']}')",
                is_pre_cmd=False,
            )

    def _forecast_observation_names(self) -> list[str]:
        """Return the pyEMU observation names registered as forecasts."""

        names: list[str] = []
        for item in self._prepared_observations:
            if not item.get("is_forecast"):
                continue
            target_frame = item.get("target_frame")
            if target_frame is None:
                continue
            for row in target_frame.itertuples(index=False):
                names.append(_observation_name(item["prefix"], row.col_label, row.row_label))
        return names

    def _apply_forecasts(self):
        """Zero forecast weights and register them as PEST++ forecasts."""

        names = [name for name in self._forecast_observation_names() if name in self.pst.observation_data.index]
        if not names:
            return
        self.pst.observation_data.loc[names, "weight"] = 0.0
        self.pst.pestpp_options["forecasts"] = ",".join(names)

    def _add_capture_field_observations(self, spec: NativeParameterSpec):
        """Record a parameter's resolved per-cell field as zero-weight observations."""

        recipe = spec.recipe
        base = f"{spec.name}field"
        # Each external array file is captured under its own prefix so that
        # multi-layer DISV K (one file per layer, each indexed 0..ncpl-1) does
        # not collide -- a single shared prefix would make layer 1 and layer 2
        # produce identical observation names.
        layer_prefixes: dict[int, str] = {}
        for filename in spec.resolved_files:
            if recipe.family == "array":
                match = re.search(r"_layer(\d+)\.txt$", str(filename))
                layer = int(match.group(1)) - 1 if match else 0
                prefix = f"{base}l{layer}" if match else base
                layer_prefixes[layer] = prefix
                self.pf.add_observations(filename, prefix=prefix)
            else:
                self.pf.add_observations(
                    filename, prefix=base, index_cols=[0, 1], use_cols=[recipe.use_col]
                )
        self._capture_field_specs[spec.name] = {
            "prefix": base,
            "layer_prefixes": layer_prefixes,
            "target": recipe.canonical,
            "family": recipe.family,
            "name": spec.name,
        }

    def _finalize_capture_fields(self):
        """Zero the weights of captured-field observations after the pst is built."""

        if not self._capture_field_specs:
            return
        obs = self.pst.observation_data
        index = obs.index.to_series()
        for info in self._capture_field_specs.values():
            prefixes = list(info.get("layer_prefixes", {}).values()) or [info["prefix"]]
            for prefix in prefixes:
                mask = index.str.contains(f"oname:{prefix.lower()}_", regex=False)
                obs.loc[mask.to_numpy(), "weight"] = 0.0

    def build(self, filename: str | Path | None = None, *, noptmax: int = 0):
        """Compile the declarations into a runnable PEST(++) control file.

        This is the modern build path. In order it: writes the model inputs as
        external array/list files (``set_all_data_external``); creates a
        ``pyemu.utils.PstFrom`` over the model workspace; registers every
        :meth:`parameterize` call natively (``pf.add_parameters``) and every
        :meth:`observe`/:meth:`forecast` set; adds the MF6 run command and the
        observation post-processors to ``forward_run.py``; and builds the
        ``.pst``. Crucially, it keeps pyEMU's own ``apply_list_and_array_pars``
        in the forward run (the legacy path deleted it), so the multiplier
        machinery is handled by pyEMU rather than hand-rolled code.

        After this call the workspace contains a self-contained PEST setup:
        the ``.pst``, ``forward_run.py``, template/instruction files, and the
        ``mult/`` multiplier files. Run it with PEST++ (e.g. ``pestpp-ies`` /
        ``pestpp-glm``) or, in a later phase, via ``cal.run_ies()``.

        Parameters
        ----------
        filename
            Output ``.pst`` filename (default ``"<name>.pst"``). Written into
            the template workspace.
        noptmax
            Initial ``NOPTMAX`` written to the control file. ``0`` (default)
            means "run the model once and compute residuals" -- the standard
            way to sanity-check a fresh setup before launching a real run.

        Returns
        -------
        pyemu.Pst
            The built control object (also available as ``cal.pst``). The
            underlying ``PstFrom`` is ``cal.pf`` for advanced edits.
        """

        target_name = Path(filename).name if filename else f"{self.name}.pst"
        self._ensure_external_model()
        self._build_pstfrom()
        self._prepare_observation_specs()
        with self._quiet_pyemu_context():
            for spec in self._native_parameter_specs:
                if spec.recipe.family == "array" and spec.style == "pilotpoints":
                    add_pilot_point_parameter(self, spec)
                else:
                    add_native_parameter(self, spec)
                if spec.capture:
                    self._add_capture_field_observations(spec)
        self.pf.mod_sys_cmds.append(self._resolve_exe())
        self._attach_native_observation_postprocessors()
        self._write_project_metadata(filename=target_name)
        with self._quiet_pyemu_context():
            self.pst = self.pf.build_pst(filename=target_name)
        if self._pilot_point_frames:
            register_pilot_point_parameters(self)
        # Pin the forward-run interpreter to this environment's Python so PEST++
        # workers use the venv that has numpy/flopy/pyemu/myflopy, not a bare
        # "python" that may resolve elsewhere.
        self.pst.model_command = [f'"{sys.executable}" forward_run.py']
        finalize_observations(self, self._prepared_observations)
        self._apply_forecasts()
        self._finalize_capture_fields()
        self.pst.control_data.noptmax = int(noptmax)
        with self._quiet_pyemu_context():
            self.pst.write(self.template_workspace / target_name, version=2)
        return self.pst

    def settings(self) -> PestSettings:
        """Return a printable snapshot of the resolved calibration configuration.

        The facade hides pyEMU boilerplate but never hides *state*: this lets
        you review exactly what was declared (every parameter with its style,
        bounds, physical limits, and transform; every observation and forecast)
        and -- once :meth:`build` has run -- the resulting control-file counts
        (number of parameters and groups, observations, nonzero-weight
        observations, forecasts, and ``NOPTMAX``).

        Returns
        -------
        PestSettings
            A dataclass with a readable ``str``/``repr`` (``print(cal.settings())``)
            plus ``.parameter_frame()`` / ``.observation_frame()`` accessors for
            the declared parameters and observations as ``pandas`` tables.

        Examples
        --------
        >>> print(cal.settings())          # before build: declared config
        >>> pst = cal.build("calib.pst")
        >>> print(cal.settings())          # after build: + control-file counts
        """

        parameters = [
            {
                "target": spec.recipe.canonical,
                "style": spec.style,
                "bounds": tuple(spec.bounds),
                "physical": tuple(spec.physical) if spec.physical is not None else None,
                "transform": spec.resolved_transform,
                "additive": bool(spec.additive),
                "name": spec.name,
            }
            for spec in self._native_parameter_specs
        ]

        def _obs_rows(specs, prepared_filter):
            rows = []
            prepared = {
                item.get("prefix"): item
                for item in self._prepared_observations
                if bool(item.get("is_forecast")) == prepared_filter
            }
            for spec in specs:
                prefix = getattr(spec, "prefix", "?")
                item = prepared.get(prefix)
                target_frame = item.get("target_frame") if item else None
                rows.append(
                    {
                        "prefix": prefix,
                        "kind": type(spec).__name__.replace("ObservationSpec", "").lower(),
                        "n": int(len(target_frame)) if target_frame is not None else None,
                    }
                )
            return rows

        observations = _obs_rows(self._observation_specs, prepared_filter=False)
        forecasts = _obs_rows(self._forecast_specs, prepared_filter=True)

        built = self.pst is not None
        npar = npar_groups = nobs = nnz_obs = n_forecasts = noptmax = None
        if built:
            par_data = self.pst.parameter_data
            npar = int(self.pst.npar)
            npar_groups = int(par_data["pargp"].nunique())
            nobs = int(self.pst.nobs)
            nnz_obs = int(self.pst.nnz_obs)
            forecast_option = self.pst.pestpp_options.get("forecasts", "")
            n_forecasts = len([name for name in forecast_option.split(",") if name]) if forecast_option else 0
            noptmax = int(self.pst.control_data.noptmax)

        return PestSettings(
            name=self.name,
            model_name=self.model.name,
            original_workspace=str(self.original_workspace),
            template_workspace=str(self.template_workspace),
            start_datetime=str(self.start_datetime),
            parameters=parameters,
            observations=observations,
            forecasts=forecasts,
            built=built,
            npar=npar,
            npar_groups=npar_groups,
            nobs=nobs,
            nnz_obs=nnz_obs,
            n_forecasts=n_forecasts,
            noptmax=noptmax,
        )

    # -- running PEST++ ---------------------------------------------------

    def _resolve_pestpp(self, exe: str):
        """Return ``(exe_name, exe_dir)`` for a PEST++ tool next to the MF6 binary."""

        name = exe[:-4] if str(exe).lower().endswith(".exe") else str(exe)
        mf6 = Path(self._resolve_exe())
        if mf6.is_absolute() and mf6.parent.exists() and (mf6.parent / f"{name}.exe").exists():
            return name, str(mf6.parent)
        return name, None

    @contextlib.contextmanager
    def _augmented_path(self, exe_dir):
        """Temporarily prepend ``exe_dir`` to ``PATH`` so PEST++ tools resolve."""

        import os

        original = os.environ.get("PATH", "")
        if exe_dir and exe_dir not in original:
            os.environ["PATH"] = exe_dir + os.pathsep + original
        try:
            yield
        finally:
            os.environ["PATH"] = original

    def run_ies(
        self,
        *,
        reals: int = 50,
        iterations: int = 3,
        workers: int | None = None,
        noise: bool = True,
        bad_phi_sigma: float | None = None,
        master_dir: str | Path | None = None,
        exe: str = "pestpp-ies",
        **pestpp_options,
    ):
        """Run iterative ensemble smoother history matching with PESTPP-IES.

        This is the one-line entry point to ensemble calibration and "free"
        uncertainty analysis. It sets sensible PESTPP-IES options on the built
        control file, launches the run (serial, or parallel across ``workers``
        agents), and returns an :class:`~myflopy.modflow.mf6.pest.ies.IesResults`
        for assessing the outcome (phi convergence, ensembles vs observations,
        posterior forecast distributions).

        Call :meth:`build` first. Reasonable defaults make ``cal.run_ies()`` a
        valid first call; tune from there.

        Parameters
        ----------
        reals
            Number of realizations in the ensemble (``ies_num_reals``). More is
            better for posterior coverage but costs more model runs; 50 is a
            reasonable starting point, production runs often use 100-300.
        iterations
            Number of smoother iterations (``NOPTMAX``). 3 is a common default;
            PESTPP-IES may stop earlier if it converges.
        workers
            Number of parallel agents. ``None``/``1`` runs serially in the
            template workspace; ``>1`` deploys a master + worker pool.
        noise
            Whether to carry measurement noise into the posterior (recommended).
            Set ``False`` to disable (``ies_no_noise``).
        bad_phi_sigma
            Optional adaptive rejection threshold (``ies_bad_phi_sigma``);
            values of 1.5 (aggressive) to 2.5 (tolerant) help highly nonlinear
            problems by dropping lagging realizations.
        master_dir
            Master directory for a parallel run (default ``<name>_ies_master``
            beside the template workspace).
        exe
            PEST++ IES executable name (resolved next to the MF6 binary, then on
            ``PATH``).
        **pestpp_options
            Any additional ``pst.pestpp_options`` to set (e.g.
            ``ies_localizer="loc.mat"``), passed straight through.

        Returns
        -------
        IesResults
            Reader/visualizer for the completed run.
        """

        from myflopy.modflow.mf6.pest.ies import IesResults

        results_dir = self._launch_pestpp_ies(
            reals=reals,
            noptmax=iterations,
            workers=workers,
            exe=exe,
            master_suffix="ies_master",
            master_dir=master_dir,
            noise=noise,
            bad_phi_sigma=bad_phi_sigma,
            pestpp_options=pestpp_options,
        )
        return IesResults(results_dir, case_name=self.name, model=self.model)

    def prior(
        self,
        *,
        reals: int = 50,
        workers: int | None = None,
        noise: bool = True,
        master_dir: str | Path | None = None,
        exe: str = "pestpp-ies",
        **pestpp_options,
    ):
        """Run the prior parameter ensemble once (prior Monte Carlo).

        Before history matching, this draws ``reals`` parameter sets from the
        prior and runs the model once for each -- "everything the model could
        plausibly do given only expert knowledge, without looking at the data."
        It is the cheap go/no-go check that the *prior brackets the observations*;
        if it does not, you have prior-data conflict (a model problem, not
        something history matching can fix), and history matching is premature.

        Mechanically this is PESTPP-IES with ``NOPTMAX=-1`` (evaluate the prior
        ensemble and stop). Call :meth:`build` first.

        Parameters
        ----------
        reals
            Number of prior realizations (``ies_num_reals``). 50 is plenty for a
            visual bracketing check.
        workers
            Parallel agents (``None``/``1`` runs serially).
        noise
            Generate the measurement-noise ensemble too (default ``True``).
        master_dir
            Master directory for a parallel run (default ``<name>_prior_master``).
        exe
            PEST++ IES executable (resolved next to the MF6 binary, then ``PATH``).
        **pestpp_options
            Additional ``pst.pestpp_options`` passed straight through.

        Returns
        -------
        IesResults
            Use :meth:`~...ies.IesResults.plot_prior_vs_obs` (grey prior spaghetti
            vs measured) and :meth:`~...ies.IesResults.conflict` to inspect it.
        """

        from myflopy.modflow.mf6.pest.ies import IesResults

        results_dir = self._launch_pestpp_ies(
            reals=reals,
            noptmax=-1,
            workers=workers,
            exe=exe,
            master_suffix="prior_master",
            master_dir=master_dir,
            noise=noise,
            bad_phi_sigma=None,
            pestpp_options=pestpp_options,
        )
        return IesResults(results_dir, case_name=self.name, model=self.model)

    def _inject_geostatistical_prior(self, reals: int) -> None:
        """Draw a geostatistically-correlated prior parameter ensemble.

        Grid / pilot-point parameters carry a variogram (``correlation=``). A
        plain bounds-based prior would treat each cell independently, so prior
        realizations come out as per-cell spatial noise and IES cannot recover a
        coherent field. Drawing from the prior covariance instead yields *smooth*
        prior realizations (handed to PESTPP-IES via ``ies_par_en``) -- this is
        what makes posterior property-pattern maps meaningful rather than
        laughable. No-op when there are no geostatistical parameters.
        """

        # Grid parameters live in pyEMU's PstFrom (pf.draw covers them). Pilot
        # points are PEST template parameters interpolated by IDW, so they are
        # not in pf -- their prior is drawn from bounds by PESTPP-IES, and the
        # IDW interpolation itself provides spatial smoothness.
        spatial = [
            spec
            for spec in self._native_parameter_specs
            if spec.style == "grid" and spec.correlation is not None
        ]
        if not spatial or self.pf is None:
            return
        with self._quiet_pyemu_context():
            ensemble = self.draw_prior(int(reals), use_specsim=False)
        ensemble_path = self.template_workspace / f"{self.name}.prior_par.csv"
        ensemble.to_csv(ensemble_path)
        self.pst.pestpp_options["ies_par_en"] = ensemble_path.name

    def _launch_pestpp_ies(
        self,
        *,
        reals: int,
        noptmax: int,
        workers: int | None,
        exe: str,
        master_suffix: str,
        master_dir: str | Path | None = None,
        noise: bool = True,
        bad_phi_sigma: float | None = None,
        pestpp_options: dict | None = None,
    ) -> Path:
        """Configure options, write the pst, run PESTPP-IES, and return the results dir."""

        if self.pst is None:
            raise ValueError("Call build() before running PESTPP-IES.")
        case = f"{self.name}.pst"
        self.pst.pestpp_options["ies_num_reals"] = int(reals)
        if not noise:
            self.pst.pestpp_options["ies_no_noise"] = True
        if bad_phi_sigma is not None:
            self.pst.pestpp_options["ies_bad_phi_sigma"] = float(bad_phi_sigma)
        for key, value in (pestpp_options or {}).items():
            self.pst.pestpp_options[key] = value
        self._inject_geostatistical_prior(int(reals))
        self.pst.control_data.noptmax = int(noptmax)
        with self._quiet_pyemu_context():
            self.pst.write(str(self.template_workspace / case), version=2)

        exe_name, exe_dir = self._resolve_pestpp(exe)
        pyemu = self.pyemu or _import_pyemu()
        if workers and int(workers) > 1:
            master = Path(master_dir) if master_dir else self.template_workspace.parent / f"{self.name}_{master_suffix}"
            with self._augmented_path(exe_dir):
                pyemu.os_utils.start_workers(
                    str(self.template_workspace),
                    exe_name,
                    case,
                    num_workers=int(workers),
                    worker_root=str(master.parent),
                    master_dir=str(master),
                )
            return master

        with self._augmented_path(exe_dir):
            pyemu.os_utils.run(f"{exe_name} {case}", cwd=str(self.template_workspace))
        return self.template_workspace
