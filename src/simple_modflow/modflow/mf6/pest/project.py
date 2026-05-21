"""Project orchestration for the first ``simple_modflow`` PEST slice."""

from __future__ import annotations

import contextlib
import io
import json
import sys
import warnings
from pathlib import Path

from simple_modflow.modflow.mf6.pest.observations import (
    finalize_observations,
    prepare_head_target_observations,
    prepare_lake_stage_observations,
)
from simple_modflow.modflow.mf6.pest.parameters import prepare_parameter_spec, register_parameter_spec
from simple_modflow.modflow.mf6.pest.specs import (
    HeadTargetObservationSpec,
    KPilotPointParameter,
    DrainElevationParameter,
    DrainConductanceParameter,
    LakeStageObservationSpec,
)

METADATA_FILENAME = "simple_modflow_pest_metadata.json"


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
    """Build a reusable pyEMU/Pest workspace from a ``SimulationBase`` model.

    The first implementation slice supports:

    - head-target observations
    - lake-stage observations
    - pilot-point K parameter metadata and forward-run application
    - drain elevation/conductance parameter metadata and forward-run application
    - head-target output regeneration during the forward run
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
        """Create a calibration project around an existing model workspace."""

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

    def _ensure_original_workspace(self):
        """Ensure MF6 input files exist before creating the PEST template."""

        self.original_workspace.mkdir(parents=True, exist_ok=True)
        self.model.sim.write_simulation(silent=True)

    def _build_pstfrom(self):
        """Instantiate the underlying ``pyemu.utils.PstFrom`` object."""

        self.pyemu = _import_pyemu()
        from pyemu.utils.pst_from import PstFrom

        self.pf = PstFrom(
            original_d=self.original_workspace,
            new_d=self.template_workspace,
            remove_existing=True,
            longnames=self.longnames,
            spatial_reference=self.spatial_reference,
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

    def _prepare_observation_specs(self):
        """Create simulated observation files and register them with pyEMU."""

        prepared: list[dict] = []
        for spec in self._observation_specs:
            if isinstance(spec, HeadTargetObservationSpec):
                prepared.append(prepare_head_target_observations(self, spec))
            elif isinstance(spec, LakeStageObservationSpec):
                prepared.append(prepare_lake_stage_observations(self, spec))
            else:
                raise TypeError(f"Unsupported observation spec type: {type(spec).__name__}")
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
            raise ValueError("Call build_pst() before draw_prior().")
        return self.pf.draw(num_reals=num_reals, use_specsim=use_specsim)
