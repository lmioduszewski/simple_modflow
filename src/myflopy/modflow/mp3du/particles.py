from __future__ import annotations

import json
import os
import shutil
import subprocess
from collections import Counter
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import geopandas as gpd
import numpy as np
import pandas as pd
from flopy.mf6.mfbase import FlopyException, MFDataException

from myflopy._deprecation import deprecated_module_getattr, warn_deprecated
from myflopy._logging import get_logger
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase

logger = get_logger(__name__)


__all__ = [
    "ParticleTrackingInput",
    "ParticleTrackingResult",
    "prepare_particle_tracking",
    "run_particle_tracking",
]


_MODULE_DIR = Path(__file__).resolve().parent
_REPO_TOOLS_DIR = _MODULE_DIR.parents[3] / "tools" / "mp3du"
_LEGACY_PRT_NAMES = {"PRT", "PrtMip", "PrtOc", "PrtPrp", "PrtDisv", "PrtFmi"}


def resolve_mp3du_executable(name: str, explicit: Path | str | None = None) -> Path:
    """Resolve an MP3DU executable path (implementation plan 2.4).

    Resolution order: explicit argument -> ``MYFLOPY_MP3DU_DIR`` environment
    variable -> the repo's untracked ``tools/mp3du/`` -> the deprecated
    package-directory fallback (executables no longer ship inside ``src/``;
    see ``tools/mp3du/README.md`` for where to put them).
    """

    if explicit is not None:
        return Path(explicit)
    env_dir = os.environ.get("MYFLOPY_MP3DU_DIR")
    if env_dir:
        return Path(env_dir) / name
    repo_candidate = _REPO_TOOLS_DIR / name
    if repo_candidate.exists():
        return repo_candidate
    legacy = _MODULE_DIR / name
    if legacy.exists():
        # Behavior deprecation (not an alias, so not in __compatibility__).
        warn_deprecated(
            f"myflopy.modflow.mp3du: loading {name} from inside the package",
            "tools/mp3du/ or MYFLOPY_MP3DU_DIR",
            since="0.1.0",
        )
        return legacy
    return repo_candidate


@dataclass(frozen=True)
class ParticleTrackingResult(Mapping[str, Any]):
    """Structured result from preparing and optionally running one MP3DU job."""

    json_file: Path
    path_file: Path
    output_json: Path | None = None
    start_cell_diagnostics: dict[str, Any] | None = None
    endpoint_summary: dict[str, int] | None = None
    diagnostics_file: Path | None = None

    def to_dict(self) -> dict[str, Any]:
        """The result fields as a plain dict (for mapping-style access)."""

        return {
            "json_file": self.json_file,
            "path_file": self.path_file,
            "output_json": self.output_json,
            "start_cell_diagnostics": self.start_cell_diagnostics,
            "endpoint_summary": self.endpoint_summary,
            "diagnostics_file": self.diagnostics_file,
        }

    def __getitem__(self, key: str) -> Any:
        """Mapping-style access to a result field by name."""

        return self.to_dict()[key]

    def __iter__(self) -> Iterator[str]:
        """Iterate the result field names."""

        return iter(self.to_dict())

    def __len__(self) -> int:
        """The number of result fields."""

        return len(self.to_dict())


class ParticleTrackingInput:
    """Assemble (and optionally run) mod-PATH3DU particle inputs for one MF6 model.

    The legacy, external-tool path to particle tracking: it translates a flow
    model and a set of release locations into the input files mod-PATH3DU expects
    -- mapping GIS attribute columns to release fields via ``_DEFAULT_FIELD_ALIASES``
    and assigning boundary-package IFACE codes via ``_DEFAULT_IFACE_OVERRIDES`` --
    then can invoke the mod-PATH3DU executable. Prefer the native MF6
    :class:`~myflopy.modflow.mf6.prt.PRTProject` (``model.particle_tracking.prt(...)``)
    for new work; use this when you specifically need mod-PATH3DU.

    Construction resolves the executable paths, stores the options and creates
    ``output_path``; no input file is read, checked or written until a method
    needs it (:meth:`run`, :meth:`validate_inputs`, ...). The three tools are run
    as subprocesses with ``output_path`` as the working directory.

    Parameters
    ----------
    model : SimulationBase, optional
        The flow model. Its ``name`` names the default MF6 output files and the
        MP3DU outputs, its ``model_output_folder_path`` is where the MF6 outputs
        are copied from (and the default ``output_path``), its ``vor`` locates
        particles on cells, and its ``gwf.modelgrid``/``gwf.package_dict`` give the
        layer count, the idomain and the boundary packages. ``None`` (the default)
        only constructs when ``output_path`` is given (otherwise
        ``AttributeError``), and :meth:`validate_inputs` and :meth:`run` then raise
        ``ValueError``.
    writep3dgsf_path : str or pathlib.Path, optional
        The writeP3DGSF executable, which turns the MF6 binary grid file into the
        grid specification (GSF) file mod-PATH3DU reads; every :meth:`run` runs
        it. ``None`` (the default) resolves ``writep3dgsf.exe`` from the
        ``MYFLOPY_MP3DU_DIR`` directory, else the repository's ``tools/mp3du/``,
        else (deprecated, with a warning) the package directory. A missing file
        raises ``FileNotFoundError`` in :meth:`validate_inputs`.
    mp3du_path : str or pathlib.Path, optional
        The mod-PATH3DU executable (``mp3du.exe``), resolved the same way.
        Required only when :meth:`run` is asked to ``execute``.
    writep3doutput_path : str or pathlib.Path, optional
        The writeP3DOutput executable (``writep3doutput.exe``), resolved the same
        way, which converts mod-PATH3DU's binary pathline output into shapefiles
        and a DBF table. Required only when :meth:`run` executes AND
        ``convert_output``.
    model_output_files : dict of str to str, optional
        File names of the MF6 files mod-PATH3DU reads, under the keys ``"grb"``,
        ``"tdis"``, ``"hds"``, ``"cbc"`` and ``"gsf"``. ``None`` (the default)
        uses ``<model name>.disv.grb``, ``.tdis``, ``.hds``, ``.cbc`` and ``.gsf``
        -- names derived from the model name, not read from its OC or name file.
        A dict replaces the defaults wholesale, so give all five keys. Each file
        but the GSF (which :meth:`run` generates) must be in ``output_path`` or in
        the model's output folder, whence it is copied; ``FileNotFoundError``
        otherwise.
    output_path : str or pathlib.Path, optional
        Workspace for every MP3DU input and output, created (with its parents) at
        construction. ``None`` (the default) uses the model's output folder.
    porosities_by_layer : Sequence[float], optional
        One porosity per model layer, written as a ``CONSTANT`` ``POROSITY`` line
        per layer in the PATH file (``mp3du.p3d``). Required: :meth:`validate_inputs`
        raises ``ValueError`` when it is ``None`` and ``AssertionError`` unless it
        holds ``int``/``float`` instances, exactly one per layer.
    particle_shp : str or pathlib.Path, optional
        A ``.shp`` or ``.gpkg`` (first layer) of release locations. It is passed to
        mod-PATH3DU as-is when every feature is a point or multipoint and it has
        cell-id, release-time and z-location columns under one of the default names
        (``cells``/``Node``/``node``/``P3D_CellID``/``CELLID``/``cellid``;
        ``TimeRel``/``TREL``/``time``/``Time``/``trel``; ``ZLoc``/``ZLOC``/``zloc``).
        Otherwise its features are mapped onto model cells
        (``model.vor.get_vor_cells_as_series``) and a point file with one particle
        per distinct cell is generated (``mp3du_particles_from_vector.shp``), as
        for ``particle_cells``. Exactly one of ``particle_shp`` and
        ``particle_cells`` is required (``ValueError`` otherwise;
        ``FileNotFoundError`` when the file is missing).
    particle_cells : list of int, optional
        Zero-based cell (CELL2D) numbers; one particle is generated at each cell's
        polygon centroid, into ``mp3du_particle_cells.shp`` in ``output_path``
        (columns ``P3D_CellID`` = cell + 1, ``TimeRel``, ``ZLoc``, ``LocName`` and
        ``Cell0``). A number outside ``0 .. ncpl - 1`` raises ``ValueError``.
    particle_field_map : dict, optional
        Which columns of the particle file supply mod-PATH3DU's release fields:
        ``"CELLID_ATTR"``, ``"TIME_ATTR"`` and ``"ZLOC_ATTR"`` (one column name
        each) and ``"ADDTL_ATTR"`` (a name or list of extra columns, written to the
        JSON as ``ADDTL_ATTR``). A key left out is found among the default names above
        (``ADDTL_ATTR`` among ``LocName``/``locname``/``name``/``Name``), and a
        named or required column that is missing raises ``ValueError``. For a
        generated particle file an empty map (the default) uses the generated
        columns. It is NOT consulted when deciding whether a ``particle_shp`` can
        be used as-is (see above).
    cellid_index_base : int, default 0
        The base of the cell ids in the particle file's cell-id column: subtracted
        to compare them with the cell each geometry falls in, for the start-cell
        diagnostics only (the file reaches mod-PATH3DU unchanged). Forced to 1
        whenever the particle file is generated (``particle_cells``, or a
        ``particle_shp`` that had to be regenerated), whose ``P3D_CellID`` is
        one-based.
    iface_overrides : dict of str to int, optional
        IFACE codes by boundary package type (keys upper-cased), merged over the
        defaults CHD 2, DRN 7, GHB 2, RCH 6, RIV 6, SFR 6, EVT 6 and WEL 0. Only
        the types of the model's packages that name cells (in
        ``stress_period_data``/``packagedata``/``connectiondata``) are written to
        the JSON's ``IFACE`` list -- all of them when no such package is found or
        the packages cannot be read.
    direction : str, default "FORWARD"
        Tracking direction, upper-cased and written unchecked as the pathline
        block's ``DIRECTION``.
    simulation_end_time : float, optional
        Written as the pathline block's ``SIMULATION_END_TIME`` (model time
        units). ``None`` (the default) leaves it out.
    flow_thread_count : int, default 10
        Written as the flow-model block's ``THREAD_COUNT``.
    pathline_thread_count : int, default 4
        Written as the pathline block's ``THREAD_COUNT``.
    initial_stepsize : float, default 0.1
        Written as the pathline block's ``INITIAL_STEPSIZE``.
    euler_dt : float, default 1.0e-4
        Written as the pathline block's ``EULER_DT``.
    adaptive_step_error : float, default 1.0e-6
        Written as the pathline block's ``ADAPTIVE_STEP_ERROR``.
    capture_radius : float, default 10.0
        Written as the pathline block's ``CAPTURE_RADIUS``.
    tracking_options : list of str, optional
        Written as the pathline block's ``OPTIONS``. ``None`` (the default) writes
        ``["TRACK_TO_TERMINATION"]``.
    generated_particle_zloc : float, default 0.95
        The ``ZLoc`` value (what mod-PATH3DU reads as the z location) stamped on
        every GENERATED particle; ignored for a particle file used as-is.
    generated_particle_release_time : float, default 0.0
        The ``TimeRel`` (release time) stamped on every generated particle.
    generated_particle_label_prefix : str, default "cell_"
        Generated particles are labelled ``<prefix><cell>`` in ``LocName``; a
        regenerated ``particle_shp`` replaces the prefix with
        ``<file stem>_cell_``.
    """

    _DEFAULT_FIELD_ALIASES = {
        "CELLID_ATTR": ("cells", "Node", "node", "P3D_CellID", "CELLID", "cellid"),
        "TIME_ATTR": ("TimeRel", "TREL", "time", "Time", "trel"),
        "ZLOC_ATTR": ("ZLoc", "ZLOC", "zloc"),
        "ADDTL_ATTR": ("LocName", "locname", "name", "Name"),
    }
    _DEFAULT_IFACE_OVERRIDES = {
        "CHD": 2,
        "DRN": 7,
        "GHB": 2,
        "RCH": 6,
        "RIV": 6,
        "SFR": 6,
        "EVT": 6,
        "WEL": 0,
    }
    _BOUNDARY_PACKAGE_NAMES = {
        "CHD",
        "DRN",
        "EVT",
        "GHB",
        "LAK",
        "RCH",
        "RIV",
        "SFR",
        "UZF",
        "WEL",
    }

    def __init__(
        self,
        model: SimulationBase = None,
        writep3dgsf_path: Path = None,
        mp3du_path: Path = None,
        writep3doutput_path: Path = None,
        model_output_files: dict | None = None,
        output_path: Path | None = None,
        porosities_by_layer: list[float] | None = None,
        particle_shp: Path | None = None,
        particle_cells: list[int] | None = None,
        particle_field_map: dict[str, Any] | None = None,
        cellid_index_base: int = 0,
        iface_overrides: dict[str, int] | None = None,
        direction: str = "FORWARD",
        simulation_end_time: float | None = None,
        flow_thread_count: int = 10,
        pathline_thread_count: int = 4,
        initial_stepsize: float = 0.1,
        euler_dt: float = 1.0e-4,
        adaptive_step_error: float = 1.0e-6,
        capture_radius: float = 10.0,
        tracking_options: list[str] | None = None,
        generated_particle_zloc: float = 0.95,
        generated_particle_release_time: float = 0.0,
        generated_particle_label_prefix: str = "cell_",
    ):
        """Configure a mod-PATH3DU particle-tracking input build for ``model``.

        Resolves executable paths and stores the release source (a particle
        shapefile or explicit cells), field mapping, IFACE overrides, and the
        numeric tracking options; see the class docstring for the workflow.
        """

        self.model = model
        self.writep3dgsf_path = resolve_mp3du_executable("writep3dgsf.exe", writep3dgsf_path)
        self.mp3du_path = resolve_mp3du_executable("mp3du.exe", mp3du_path)
        self.writep3doutput_path = resolve_mp3du_executable("writep3doutput.exe", writep3doutput_path)
        self._model_output_files = model_output_files
        self._output_path = Path(output_path) if output_path is not None else None
        self._porosities_by_layer = porosities_by_layer
        self.particle_shp = None if particle_shp is None else Path(particle_shp)
        self.particle_cells = None if particle_cells is None else [int(cell) for cell in particle_cells]
        self._particle_field_map_input = dict(particle_field_map or {})
        self._resolved_particle_field_map = None
        self._resolved_particle_input_path = None
        self._generated_from_vector = False
        self.cellid_index_base = 1 if (particle_cells is not None and particle_shp is None) else int(cellid_index_base)
        self._iface_overrides = self._coerce_iface_overrides(iface_overrides)
        self.direction = direction.upper()
        self.simulation_end_time = simulation_end_time
        self.flow_thread_count = int(flow_thread_count)
        self.pathline_thread_count = int(pathline_thread_count)
        self.initial_stepsize = float(initial_stepsize)
        self.euler_dt = float(euler_dt)
        self.adaptive_step_error = float(adaptive_step_error)
        self.capture_radius = float(capture_radius)
        self.tracking_options = ["TRACK_TO_TERMINATION"] if tracking_options is None else list(tracking_options)
        self.generated_particle_zloc = float(generated_particle_zloc)
        self.generated_particle_release_time = float(generated_particle_release_time)
        self.generated_particle_label_prefix = str(generated_particle_label_prefix)
        self._variables = None
        self.path_file_path = self.output_path / "mp3du.p3d"

    @staticmethod
    def _coerce_iface_overrides(iface_overrides: dict[str, int] | None) -> dict[str, int]:
        """Merge user IFACE overrides (upper-cased) onto the package defaults."""

        overrides = dict(ParticleTrackingInput._DEFAULT_IFACE_OVERRIDES)
        if iface_overrides:
            for key, value in iface_overrides.items():
                overrides[str(key).upper()] = int(value)
        return overrides

    @property
    def model_output_files(self):
        """The MF6 output filenames MP3DU reads (grb/tdis/hds/cbc/gsf), inferred if unset (cached)."""

        if self._model_output_files is None:
            model_name = self.model.name
            self._model_output_files = {
                "grb": f"{model_name}.disv.grb",
                "tdis": f"{model_name}.tdis",
                "hds": f"{model_name}.hds",
                "cbc": f"{model_name}.cbc",
                "gsf": f"{model_name}.gsf",
            }
        return self._model_output_files

    @property
    def output_path(self) -> Path:
        """The directory MP3DU inputs/outputs are written to (the model's output folder by default)."""

        if self._output_path is None:
            self._output_path = self.model.model_output_folder_path
        self._output_path.mkdir(parents=True, exist_ok=True)
        return self._output_path

    @property
    def porosities_by_layer(self):
        """Per-layer aquifer porosities used to convert flow to velocity."""

        return self._porosities_by_layer

    @porosities_by_layer.setter
    def porosities_by_layer(self, val):
        """Set per-layer porosities, asserting one numeric value per model layer."""

        assert isinstance(val, list), "porosities_by_layer must be a list of porosities, one per layer"
        assert len(val) == self.model.gwf.modelgrid.nlay, "number of porosities must match number of layers"
        for por in val:
            assert isinstance(por, (int, float)), f"porosity value: {por} must be an integer or a float (decimal)"
        self._porosities_by_layer = val

    @property
    def variables(self):
        """The MP3DU transport variables block (velocity method, porosity, retardation, dispersion)."""

        if self._variables is None:
            porosities = self.porosities_by_layer
            self._variables = {
                "VELOCITY METHOD LAYER": 3,
                "POROSITY": porosities,
                "RETARDATION": 1,
                "DispH": 0,
                "DISPT": 0,
                "DISPV": 0,
            }
        return self._variables

    @property
    def particle_field_map(self) -> dict[str, Any]:
        """The resolved mapping of MP3DU release fields to source-attribute columns (cached)."""

        if self._resolved_particle_field_map is None:
            self._resolved_particle_field_map = self._resolve_particle_field_map()
        return self._resolved_particle_field_map

    @property
    def particle_input_path(self) -> Path:
        """The path to the resolved particle-release input file (cached)."""

        if self._resolved_particle_input_path is None:
            self._resolved_particle_input_path = self._resolve_particle_input_path()
        return self._resolved_particle_input_path

    def _require_model(self):
        """Raise if no flow model is bound (required for input generation)."""

        if self.model is None:
            raise ValueError("model is required for MP3DU input generation.")

    def _require_particle_source(self):
        """Validate exactly one release source is given (a particle file or explicit cells) and exists."""

        if self.particle_shp is None and self.particle_cells is None:
            raise ValueError("Either particle_shp or particle_cells is required for MP3DU input generation.")
        if self.particle_shp is not None and self.particle_cells is not None:
            raise ValueError("Specify only one particle source: particle_shp or particle_cells.")
        if self.particle_shp is not None and not self.particle_shp.exists():
            raise FileNotFoundError(f"Particle file not found: {self.particle_shp}")

    def _resolve_particle_input_path(self) -> Path:
        """The particle-release input file path: the given vector, or one generated from cells."""

        self._require_particle_source()
        if self.particle_shp is not None:
            return self._prepare_particle_vector_input()
        return self._write_generated_particle_file()

    def _prepare_particle_vector_input(self) -> Path:
        """Use the release vector directly if it is point features with the required fields, else regenerate.

        A non-point or field-incomplete vector is mapped to grid cells and a
        proper point-release file is generated from them.
        """

        particle_data = read_shp_gpkg(self.particle_shp)
        columns = set(particle_data.columns)
        has_required = all(
            any(alias in columns for alias in self._DEFAULT_FIELD_ALIASES[key])
            for key in ("CELLID_ATTR", "TIME_ATTR", "ZLOC_ATTR")
        )
        geom_types = {
            str(geom_type).upper()
            for geom_type in getattr(particle_data, "geom_type", pd.Series(dtype=str)).dropna().tolist()
        }
        is_point_input = all(geom_type in {"POINT", "MULTIPOINT"} for geom_type in geom_types) if geom_types else False
        if has_required and is_point_input:
            return self.particle_shp
        self._generated_from_vector = True
        self.cellid_index_base = 1
        return self._write_generated_particle_file_from_vector(particle_data)

    def _write_generated_particle_file_from_vector(self, particle_data: pd.DataFrame) -> Path:
        """Map vector features to unique grid cells and write a generated point-release file."""

        self._require_model()
        cells_series = self.model.vor.get_vor_cells_as_series(particle_data)
        ordered_cells: list[int] = []
        seen: set[int] = set()
        for value in cells_series.tolist():
            for cell in self._flatten_cellid(value):
                if cell in seen:
                    continue
                seen.add(cell)
                ordered_cells.append(int(cell))
        if not ordered_cells:
            raise ValueError(f"No model cells were found from particle vector input: {self.particle_shp}")
        self.particle_cells = ordered_cells
        self.generated_particle_label_prefix = f"{self.particle_shp.stem}_cell_"
        return self._write_generated_particle_file(filename="mp3du_particles_from_vector.shp")

    def _grid_cell_count(self) -> int:
        """The model's cells-per-layer count, from the Voronoi grid or model grid (raises if unknown)."""

        vor = getattr(self.model, "vor", None)
        if vor is not None:
            ncpl = getattr(vor, "ncpl", None)
            if ncpl is not None:
                return int(ncpl)
            gdf_vor_polys = getattr(vor, "gdf_vorPolys", None)
            if gdf_vor_polys is not None:
                return int(len(gdf_vor_polys))
        modelgrid = getattr(getattr(self.model, "gwf", None), "modelgrid", None)
        if modelgrid is not None:
            ncpl = getattr(modelgrid, "ncpl", None)
            if ncpl is not None:
                return int(ncpl)
        raise AttributeError("Could not determine grid cell count for particle generation.")

    def _write_generated_particle_file(self, filename: str = "mp3du_particle_cells.shp") -> Path:
        """Write a point-release shapefile with one particle at each ``particle_cells`` centroid.

        Validates the cell ids are in range and stamps release time, z-location, and
        a label on each generated point.
        """

        self._require_model()
        if self.particle_cells is None:
            raise ValueError("particle_cells are required to generate a particle shapefile.")

        grid_cell_count = self._grid_cell_count()
        invalid_cells = [cell for cell in self.particle_cells if cell < 0 or cell >= grid_cell_count]
        if invalid_cells:
            raise ValueError(f"particle_cells includes out-of-range cell IDs: {invalid_cells[:10]}")

        centroids = self.model.vor.gdf_vorPolys.geometry.centroid
        particle_stem = Path(filename).stem
        particle_path = self.output_path / filename
        for sidecar in self.output_path.glob(f"{particle_stem}.*"):
            sidecar.unlink()

        particle_gdf = gpd.GeoDataFrame(
            {
                "P3D_CellID": [cell + 1 for cell in self.particle_cells],
                "TimeRel": [self.generated_particle_release_time for _ in self.particle_cells],
                "ZLoc": [self.generated_particle_zloc for _ in self.particle_cells],
                "LocName": [f"{self.generated_particle_label_prefix}{cell}" for cell in self.particle_cells],
                "Cell0": self.particle_cells,
            },
            geometry=[centroids.iloc[cell] for cell in self.particle_cells],
            crs=self.model.vor.gdf_vorPolys.crs,
        )
        particle_gdf.to_file(particle_path, driver="ESRI Shapefile")
        return particle_path

    def _read_particle_locations(self) -> pd.DataFrame:
        """Read the resolved particle-release input file into a GeoDataFrame."""

        return read_shp_gpkg(self.particle_input_path)

    def _resolve_particle_field_map(self) -> dict[str, Any]:
        """Resolve which source columns supply each MP3DU release field.

        Uses generated-file defaults for cell-based/regenerated input, else matches
        the configured aliases against the file's columns.
        """

        particle_data = self._read_particle_locations()
        columns = set(particle_data.columns)
        explicit = dict(self._particle_field_map_input)
        particle_name = self.particle_input_path.name

        if (self.particle_cells is not None or self._generated_from_vector) and not explicit:
            return {
                "CELLID_ATTR": "P3D_CellID",
                "TIME_ATTR": "TimeRel",
                "ZLOC_ATTR": "ZLoc",
                "ADDTL_ATTR": ["LocName", "Cell0"],
            }

        resolved: dict[str, Any] = {}
        for key in ("CELLID_ATTR", "TIME_ATTR", "ZLOC_ATTR"):
            explicit_value = explicit.get(key)
            if explicit_value is not None:
                if explicit_value not in columns:
                    raise ValueError(
                        f"{key} '{explicit_value}' was not found in {particle_name}. "
                        f"Available fields: {sorted(columns)}"
                    )
                resolved[key] = explicit_value
                continue

            alias = next((name for name in self._DEFAULT_FIELD_ALIASES[key] if name in columns), None)
            if alias is None:
                raise ValueError(
                    f"Could not infer {key} from {particle_name}. "
                    f"Available fields: {sorted(columns)}"
                )
            resolved[key] = alias

        explicit_additional = explicit.get("ADDTL_ATTR")
        if explicit_additional is None:
            resolved["ADDTL_ATTR"] = [name for name in self._DEFAULT_FIELD_ALIASES["ADDTL_ATTR"] if name in columns]
        else:
            additional = [explicit_additional] if isinstance(explicit_additional, str) else list(explicit_additional)
            missing = [name for name in additional if name not in columns]
            if missing:
                raise ValueError(
                    f"ADDTL_ATTR fields {missing} were not found in {particle_name}. "
                    f"Available fields: {sorted(columns)}"
                )
            resolved["ADDTL_ATTR"] = additional
        return resolved

    @staticmethod
    def _normalize_int(value: Any) -> int | None:
        """Coerce a value to a Python int, or ``None`` for null/NaN."""

        if value is None or pd.isna(value):
            return None
        if isinstance(value, np.generic):
            value = value.item()
        return int(value)

    def _normalize_declared_cellid(self, value: Any) -> int | None:
        """A declared cell id normalized to a zero-based index (subtracting ``cellid_index_base``)."""

        normalized = self._normalize_int(value)
        if normalized is None:
            return None
        return normalized - self.cellid_index_base

    @classmethod
    def _flatten_cellid(cls, value: Any) -> list[int]:
        """Flatten any cellid shape (scalar, ``(layer, cell)``, nested list/array) to zero-based cell ints."""

        if value is None or (not isinstance(value, (list, tuple, np.ndarray)) and pd.isna(value)):
            return []
        if isinstance(value, np.ndarray):
            return cls._flatten_cellid(value.tolist())
        if isinstance(value, tuple):
            if len(value) == 2 and isinstance(value[0], (int, np.integer)) and isinstance(value[1], (int, np.integer)):
                return [int(value[1])]
            cells: list[int] = []
            for item in value:
                cells.extend(cls._flatten_cellid(item))
            return cells
        if isinstance(value, list):
            cells: list[int] = []
            for item in value:
                cells.extend(cls._flatten_cellid(item))
            return cells
        return [int(value)]

    @classmethod
    def _aligned_cells_from_series(cls, frame: pd.DataFrame, cells_series: pd.Series) -> list[int | None]:
        """One cell (or ``None``) per row of ``frame``, taken from the geometry-derived ``cells_series``."""

        aligned: list[int | None] = []
        for idx in range(len(frame)):
            if idx not in cells_series.index:
                aligned.append(None)
                continue
            flattened = cls._flatten_cellid(cells_series.loc[idx])
            aligned.append(flattened[0] if flattened else None)
        return aligned

    def get_start_cell_diagnostics(self) -> dict[str, Any]:
        """Diagnose each particle's start cell: mapped counts, declared-vs-geometry mismatches, and boundary/inactive/IFACE coverage.

        Returns a dict summarizing how many particles mapped to cells, how many fall
        on boundary packages or inactive cells, and which packages lack an IFACE override.
        """

        self._require_model()
        particle_data = self._read_particle_locations()
        declared_cells = [
            self._normalize_declared_cellid(value)
            for value in particle_data[self.particle_field_map["CELLID_ATTR"]].tolist()
        ]
        geometry_cells = self._aligned_cells_from_series(particle_data, self.model.vor.get_vor_cells_as_series(particle_data))
        start_cells = [declared if declared is not None else geometry for declared, geometry in zip(declared_cells, geometry_cells, strict=False)]

        mismatch_count = sum(
            declared is not None and geometry is not None and declared != geometry
            for declared, geometry in zip(declared_cells, geometry_cells, strict=False)
        )
        boundary_cells = self.collect_boundary_cell_sets()
        inactive_cells = self.collect_inactive_cells()

        mapped_cells = [cell for cell in start_cells if cell is not None]
        total_particles = len(start_cells)
        mapped_particles = len(mapped_cells)
        boundary_summary = {}
        for package_name, cells in boundary_cells.items():
            count = sum(cell in cells for cell in mapped_cells)
            boundary_summary[package_name] = {
                "count": count,
                "pct_mapped": round((count / mapped_particles) * 100.0, 2) if mapped_particles else 0.0,
            }

        boundary_packages_seen = sorted(boundary_cells)
        missing_iface = sorted(
            package_name
            for package_name in boundary_packages_seen
            if package_name in self._BOUNDARY_PACKAGE_NAMES and package_name not in self._iface_overrides
        )
        inactive_count = sum(cell in inactive_cells for cell in mapped_cells)

        return {
            "particle_file": str(self.particle_input_path),
            "source_particle_file": None if self.particle_shp is None else str(self.particle_shp),
            "field_map": self.particle_field_map,
            "cellid_index_base": self.cellid_index_base,
            "total_particles": total_particles,
            "mapped_particles": mapped_particles,
            "unmapped_particles": total_particles - mapped_particles,
            "declared_vs_geometry_mismatches": mismatch_count,
            "inactive_particles": inactive_count,
            "inactive_pct_mapped": round((inactive_count / mapped_particles) * 100.0, 2) if mapped_particles else 0.0,
            "boundary_packages": boundary_summary,
            "iface_overrides": dict(sorted(self._iface_overrides.items())),
            "packages_without_iface_override": missing_iface,
        }

    def collect_inactive_cells(self) -> set[int]:
        """The set of inactive (idomain==0) cell ids for the model's top layer."""

        try:
            inactive = getattr(self.model, "inactive_cells", None)
        except (TypeError, KeyError, ValueError, AssertionError, OSError,
                RuntimeError):
            # `getattr` with a default absorbs AttributeError, so anything that
            # arrives here came from `inactive_cells` being a PROPERTY that does
            # work: reading a GeoPackage (pyogrio's DataSourceError is a
            # RuntimeError), asserting its inputs are Paths (AssertionError), or
            # coercing them. Falling through to the idomain array below is the
            # documented second source.
            logger.debug("no inactive_cells on the model; reading idomain instead",
                         exc_info=True)
            inactive = None
        if inactive is not None:
            return {int(cell) for cell in inactive}

        idomain = getattr(self.model.gwf.modelgrid, "idomain", None)
        if idomain is None:
            return set()
        idomain_arr = np.asarray(idomain)
        if idomain_arr.ndim == 1:
            layer_zero = idomain_arr
        else:
            layer_zero = idomain_arr.reshape(idomain_arr.shape[0], -1)[0]
        return {int(idx) for idx, value in enumerate(layer_zero) if int(value) == 0}

    @classmethod
    def _package_type_name(cls, name: str, package) -> str:
        """The uppercase three-letter package type (e.g. ``"CHD"``) from a package or its name."""

        package_type = getattr(package, "package_type", None) or getattr(package, "_package_type", None) or name
        return str(package_type).upper().split("_")[0]

    @classmethod
    def _collect_cells_from_source(cls, source: Any) -> set[int]:
        """Recursively extract all cell ids from a package data source (arrays, dicts, records)."""

        cells: set[int] = set()
        if source is None:
            return cells

        if isinstance(source, dict):
            for value in source.values():
                cells.update(cls._collect_cells_from_source(value))
            return cells

        if not isinstance(source, (np.ndarray, list, tuple)):
            array = getattr(source, "array", None)
            if array is not None:
                cells.update(cls._collect_cells_from_source(array))
                return cells

            data = getattr(source, "data", None)
            if data is not None and data is not source:
                cells.update(cls._collect_cells_from_source(data))
                return cells
        else:
            array = source

        if isinstance(array, np.ndarray) and array.dtype.names:
            cellid_field = next((name for name in array.dtype.names if name.lower() == "cellid"), None)
            if cellid_field is not None:
                for row in array:
                    cells.update(cls._flatten_cellid(row[cellid_field]))
            return cells

        if isinstance(array, (list, tuple)):
            for item in array:
                if isinstance(item, np.void) and item.dtype.names:
                    cellid_field = next((name for name in item.dtype.names if name.lower() == "cellid"), None)
                    if cellid_field is not None:
                        cells.update(cls._flatten_cellid(item[cellid_field]))
                        continue
                if isinstance(item, (list, tuple, np.ndarray)) and len(item) > 1:
                    first_value = item[0]
                    if isinstance(first_value, (list, tuple, np.ndarray)):
                        cells.update(cls._flatten_cellid(first_value))
                        continue
                cells.update(cls._flatten_cellid(item))
            return cells

        return cells

    def collect_boundary_cell_sets(self) -> dict[str, set[int]]:
        """Map each boundary package type to the set of cells it applies to across the model."""

        self._require_model()
        boundary_cells: dict[str, set[int]] = {}
        package_dict = getattr(self.model.gwf, "package_dict", {})
        for package_name, package in package_dict.items():
            package_type = self._package_type_name(package_name, package)
            cells = set()
            for attr_name in ("stress_period_data", "packagedata", "connectiondata"):
                cells.update(self._collect_cells_from_source(getattr(package, attr_name, None)))
            if cells:
                boundary_cells[package_type] = boundary_cells.get(package_type, set()).union(cells)
        return dict(sorted(boundary_cells.items()))

    def summarize_endpoint_output(self, endpoint_path: Path | None = None) -> dict[str, int]:
        """Count particle termination reasons (``PTERM``) in the endpoint output shapefile.

        Parameters
        ----------
        endpoint_path : str or pathlib.Path, optional
            The endpoint shapefile (or GeoPackage) writeP3DOutput produced. ``None``
            (the default) reads ``<model name>_endpoint.shp`` in ``output_path``.

        Returns
        -------
        dict of str to int
            Particle count per ``PTERM`` value (as text, a missing value counted
            as ``"None"``), sorted by that text.

        Raises
        ------
        FileNotFoundError
            When the endpoint file does not exist.
        ValueError
            When it has no ``PTERM`` column.
        """

        path = self.output_path / self.output_names["ENDPOINT"] if endpoint_path is None else Path(endpoint_path)
        if not path.exists():
            raise FileNotFoundError(f"Endpoint output not found: {path}")
        endpoint = read_shp_gpkg(path)
        if "PTERM" not in endpoint.columns:
            raise ValueError(f"Endpoint output {path.name} does not include a PTERM field.")
        counts = Counter(str(value) for value in endpoint["PTERM"].fillna("None"))
        return dict(sorted(counts.items()))

    @property
    def output_names(self) -> dict[str, str]:
        """The MP3DU output filenames (pathline table/shapes, points-in-time, endpoint), prefixed by model name."""

        prefix = self.model.name
        return {
            "DBF_TABLE": f"{prefix}_pathline_table.dbf",
            "PATHLINE_WHOLE": f"{prefix}_pathline_whole.shp",
            "PATHLINE_PARTS": f"{prefix}_pathline_parts.shp",
            "POINTS_IN_TIME": f"{prefix}_points_in_time.shp",
            "ENDPOINT": f"{prefix}_endpoint.shp",
        }

    def validate_inputs(self, *, execute: bool, convert_output: bool):
        """Validate everything needed to build (and optionally run/convert) an MP3DU run.

        Checks the model, particle source, porosities, and the required executables
        exist, and eagerly resolves the particle input path + field map (which may
        write a generated particle file into ``output_path``).

        Parameters
        ----------
        execute : bool
            Whether mod-PATH3DU will be run: ``True`` also requires the
            ``mp3du_path`` executable to exist. The writeP3DGSF executable is
            required either way.
        convert_output : bool
            Whether its output will be converted: with ``execute``, also requires
            the ``writep3doutput_path`` executable. Ignored without ``execute``.

        Raises
        ------
        ValueError
            No model, no (or two) particle sources, ``porosities_by_layer`` is
            ``None``, a ``particle_cells`` number out of range, a ``particle_shp``
            that has to be regenerated but maps onto no cell, or a particle file
            that lacks a required field.
        FileNotFoundError
            The particle file or a required executable is missing.
        AssertionError
            ``porosities_by_layer`` is not one number per model layer.
        """

        self._require_model()
        self._require_particle_source()
        if self.porosities_by_layer is None:
            raise ValueError("porosities_by_layer is required for MP3DU input generation.")
        required_executables = [self.writep3dgsf_path]
        if execute:
            required_executables.append(self.mp3du_path)
        if execute and convert_output:
            required_executables.append(self.writep3doutput_path)
        for executable in required_executables:
            if not Path(executable).exists():
                raise FileNotFoundError(f"Required MP3DU executable not found: {executable}")
        self.porosities_by_layer = list(self.porosities_by_layer)
        _ = self.particle_input_path
        _ = self.particle_field_map

    def get_active_iface_overrides(self) -> dict[str, int]:
        """The IFACE overrides restricted to boundary packages actually present in the model."""

        try:
            boundary_packages = set(self.collect_boundary_cell_sets())
        except (MFDataException, FlopyException, AttributeError, TypeError,
                ValueError, KeyError, IndexError):
            # Walking every boundary package on the model reaches flopy's own
            # two exception classes (neither subclasses a builtin) as well as
            # the usual shape errors. With no packages resolved, ALL configured
            # iface overrides are returned unfiltered -- which is what the
            # caller below does when the set is empty.
            logger.debug("could not resolve boundary packages; returning every "
                         "iface override", exc_info=True)
            boundary_packages = set()
        if not boundary_packages:
            return dict(sorted(self._iface_overrides.items()))
        return {
            name: value
            for name, value in sorted(self._iface_overrides.items())
            if name in boundary_packages
        }

    def create_gsf_file(self):
        """Create the mod-PATH3DU grid specification file from the MF6 GRB."""

        gsf_json = {
            "FLOW_MODEL_TYPE": {
                "USGS_HFWK": {
                    "GRB_FILE": self.model_output_files["grb"],
                    "GSF_FILE": {
                        "TYPE": "HFWK_GRB_V.1.0.0"
                    },
                }
            },
            "OUTPUT_FILENAME": self.model_output_files["gsf"],
        }
        gsf_json_file_path = self.output_path / "grb_to_gsf.json"
        with open(gsf_json_file_path, "w", encoding="utf-8") as f:
            json.dump(gsf_json, f, indent=4)

        cmd = [self.writep3dgsf_path.as_posix(), gsf_json_file_path.name, "colorcode"]
        run = subprocess.run(
            cmd,
            cwd=self.output_path,
            capture_output=True,
            text=True,
        )

        if run.returncode != 0:
            raise RuntimeError(f"Error running writeP3DGSF.exe: {run.stderr}")
        return gsf_json_file_path

    def create_modflow_input_files(self):
        """Ensure the required MF6 output files are present in the MP3DU workspace."""

        self.output_path.mkdir(parents=True, exist_ok=True)
        for typ, file_name in self.model_output_files.items():
            destination = self.output_path / file_name
            if destination.exists():
                continue
            original_file = self.model.model_output_folder_path / file_name
            if original_file.exists():
                if original_file.resolve() == destination.resolve():
                    continue
                shutil.copyfile(original_file, destination)
            elif typ == "gsf":
                continue
            else:
                raise FileNotFoundError(f"Required file {original_file} not found.")

    def create_path_file(self):
        """Create the per-cell property file used by mod-PATH3DU."""

        with open(self.path_file_path, "w", encoding="utf-8") as f:
            f.write("# PATH3D input file\n\n")
            for variable in self.variables:
                for layer in range(self.model.gwf.modelgrid.nlay):
                    if variable == "POROSITY":
                        f.write(f"  CONSTANT    {self.variables[variable][layer]}   POROSITY {layer + 1}\n")
                    else:
                        f.write(f"  CONSTANT    {self.variables[variable]}   {variable} {layer + 1}\n")
        return self.path_file_path

    def run_mp3du(self, json_file_path: Path):
        """Run mod-PATH3DU in the configured workspace.

        Parameters
        ----------
        json_file_path : str or pathlib.Path
            The mod-PATH3DU JSON configuration (:meth:`create_json_file` writes
            ``mp3du_input.json``). Only its file NAME is passed to the tool, which
            runs in ``output_path``, so the file must be there.

        Returns
        -------
        subprocess.CompletedProcess
            The finished process, with its captured text output.

        Raises
        ------
        RuntimeError
            When the tool exits non-zero; the message carries its stderr.
        """

        cmd = [self.mp3du_path.as_posix(), Path(json_file_path).name, "colorcode"]
        result = subprocess.run(
            cmd,
            cwd=self.output_path,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"Error running mp3du.exe: {result.stderr}")
        return result

    def create_json_file(self):
        """Create the primary mod-PATH3DU JSON configuration file."""

        particle_fields = self.particle_field_map
        flow_model_block = {
            "GRB_FILE": self.model_output_files["grb"],
            "TDIS_FILE": self.model_output_files["tdis"],
            "PATH_FILE": self.path_file_path.name,
            "HDS_FILE": self.model_output_files["hds"],
            "CBB_FILE": self.model_output_files["cbc"],
            "GSF_FILE": {
                "TYPE": "GSF_V.1.1.0",
                "FILE_NAME": self.model_output_files["gsf"],
            },
            "OUTPUT_PRECISION": "DOUBLE",
            "IFACE": [{name: value} for name, value in self.get_active_iface_overrides().items()],
            "THREAD_COUNT": self.flow_thread_count,
        }
        pathline_block = {
            "NAME": self.model.name,
            "DIRECTION": self.direction,
            "THREAD_COUNT": self.pathline_thread_count,
            "INITIAL_STEPSIZE": self.initial_stepsize,
            "EULER_DT": self.euler_dt,
            "ADAPTIVE_STEP_ERROR": self.adaptive_step_error,
            "CAPTURE_RADIUS": self.capture_radius,
            "OPTIONS": self.tracking_options,
            "PARTICLE_START_LOCATIONS": {
                "SHAPEFILE": {
                    "FILE_NAME": self.particle_input_path.as_posix(),
                    "CELLID_ATTR": particle_fields["CELLID_ATTR"],
                    "TIME_ATTR": particle_fields["TIME_ATTR"],
                    "ZLOC_ATTR": particle_fields["ZLOC_ATTR"],
                }
            },
        }
        if self.simulation_end_time is not None:
            pathline_block["SIMULATION_END_TIME"] = float(self.simulation_end_time)
        if particle_fields["ADDTL_ATTR"]:
            pathline_block["PARTICLE_START_LOCATIONS"]["SHAPEFILE"]["ADDTL_ATTR"] = particle_fields["ADDTL_ATTR"]

        json_data = {
            "FLOW_MODEL_TYPE": {"USGS_HFWK": flow_model_block},
            "SIMULATIONS": [{"PATHLINE": pathline_block}],
        }

        json_file_path = self.output_path / "mp3du_input.json"
        with open(json_file_path, "w", encoding="utf-8") as f:
            json.dump(json_data, f, indent=4)
        return json_file_path

    def create_output_json_file(self) -> Path:
        """Write the MP3DU output-conversion JSON (which pathline/endpoint files to emit) and return its path."""

        output_json = {
            "MP3DU_BIN": f"{self.model.name}_PATHLINE.bin",
            "OUTPUTS": [
                {"SUMMARY": {}},
                {"DBF_TABLE": {"FILE_NAME": self.output_names["DBF_TABLE"]}},
                {"PATHLINE_WHOLE": {"FILE_NAME": self.output_names["PATHLINE_WHOLE"]}},
                {"PATHLINE_PARTS": {"FILE_NAME": self.output_names["PATHLINE_PARTS"]}},
                {"POINTS_IN_TIME": {"FILE_NAME": self.output_names["POINTS_IN_TIME"]}},
                {"ENDPOINT": {"FILE_NAME": self.output_names["ENDPOINT"]}},
            ],
        }
        output_json_path = self.output_path / "P3DOutput_json.json"
        with open(output_json_path, "w", encoding="utf-8") as f:
            json.dump(output_json, f, indent=4)
        return output_json_path

    def get_output_json(self):
        """Backwards-compatible alias for writing the output-conversion JSON."""

        return self.create_output_json_file()

    def run_output_conversion(self, output_json_path: Path):
        """Run the ``writep3doutput`` tool: convert the binary pathline output to shapefiles/tables.

        Parameters
        ----------
        output_json_path : str or pathlib.Path
            The conversion configuration (:meth:`create_output_json_file` writes
            ``P3DOutput_json.json``). Only its file NAME is passed to the tool,
            which runs in ``output_path``, so the file must be there.

        Returns
        -------
        subprocess.CompletedProcess
            The finished process, with its captured text output.

        Raises
        ------
        RuntimeError
            When the tool exits non-zero; the message carries its stderr.
        """

        cmd = [self.writep3doutput_path.as_posix(), Path(output_json_path).name, "colorcode"]
        result = subprocess.run(
            cmd,
            cwd=self.output_path,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"Error running writeP3DOutput.exe: {result.stderr}")
        return result

    def write_diagnostics_file(
        self,
        *,
        start_cells: dict[str, Any] | None = None,
        endpoint_summary: dict[str, int] | None = None,
    ) -> Path:
        """Write the start-cell and endpoint diagnostics to a JSON file and return its path.

        The file is ``mp3du_diagnostics.json`` in ``output_path``, overwritten.

        Parameters
        ----------
        start_cells : dict, optional
            Written under the ``"start_cells"`` key -- normally
            :meth:`get_start_cell_diagnostics`' result. Must be JSON-serializable;
            ``None`` (the default) writes ``null``.
        endpoint_summary : dict of str to int, optional
            Written under the ``"endpoint_summary"`` key -- normally
            :meth:`summarize_endpoint_output`' result. ``None`` (the default)
            writes ``null``.

        Returns
        -------
        pathlib.Path
        """

        diagnostics = {
            "start_cells": start_cells,
            "endpoint_summary": endpoint_summary,
        }
        diagnostics_path = self.output_path / "mp3du_diagnostics.json"
        with open(diagnostics_path, "w", encoding="utf-8") as f:
            json.dump(diagnostics, f, indent=4)
        return diagnostics_path

    def run(
        self,
        *,
        execute: bool = True,
        convert_output: bool = True,
        write_diagnostics: bool = True,
    ) -> ParticleTrackingResult:
        """Create MP3DU inputs and optionally execute the full run.

        In order: :meth:`validate_inputs`, copy the MF6 output files into
        ``output_path``, compute the start-cell diagnostics, run writeP3DGSF to make
        the GSF file (always -- even without ``execute``), write the PATH file and
        ``mp3du_input.json``, then optionally run and convert.

        Parameters
        ----------
        execute : bool, default True
            Run mod-PATH3DU on the inputs just written. ``False`` stops after
            writing them.
        convert_output : bool, default True
            With ``execute``, write ``P3DOutput_json.json``, run writeP3DOutput to
            produce the pathline, points-in-time and endpoint shapefiles and the
            DBF table, and summarize the endpoint file when it exists. Ignored
            without ``execute``.
        write_diagnostics : bool, default True
            Write ``mp3du_diagnostics.json`` (:meth:`write_diagnostics_file`) with
            the start-cell diagnostics and the endpoint summary.

        Returns
        -------
        ParticleTrackingResult
            The written file paths and diagnostics; ``output_json`` and
            ``endpoint_summary`` are ``None`` for a step that did not run, and
            ``diagnostics_file`` is ``None`` without ``write_diagnostics``.

        Raises
        ------
        ValueError, FileNotFoundError, AssertionError
            From :meth:`validate_inputs`, or a missing MF6 output file.
        RuntimeError
            When one of the three tools exits non-zero.
        """

        self.validate_inputs(execute=bool(execute), convert_output=bool(convert_output))
        self.create_modflow_input_files()
        start_cell_diagnostics = self.get_start_cell_diagnostics()
        self.create_gsf_file()
        self.create_path_file()
        json_file_path = self.create_json_file()

        output_json_path = None
        endpoint_summary = None
        if execute:
            self.run_mp3du(json_file_path)
            if convert_output:
                output_json_path = self.create_output_json_file()
                self.run_output_conversion(output_json_path)
                endpoint_path = self.output_path / self.output_names["ENDPOINT"]
                if endpoint_path.exists():
                    endpoint_summary = self.summarize_endpoint_output(endpoint_path)

        diagnostics_path = None
        if write_diagnostics:
            diagnostics_path = self.write_diagnostics_file(
                start_cells=start_cell_diagnostics,
                endpoint_summary=endpoint_summary,
            )

        return ParticleTrackingResult(
            json_file=json_file_path,
            path_file=self.path_file_path,
            output_json=output_json_path,
            start_cell_diagnostics=start_cell_diagnostics,
            endpoint_summary=endpoint_summary,
            diagnostics_file=diagnostics_path,
        )


def prepare_particle_tracking(
    *,
    model,
    particles: Path | str | list[int] | tuple[int, ...],
    porosity: float | list[float] = 0.2,
    output_path: Path | str | None = None,
    execute: bool = False,
    convert_output: bool = False,
    write_diagnostics: bool = True,
    **kwargs,
) -> tuple[ParticleTrackingInput, ParticleTrackingResult]:
    """Create and optionally run a MP3DU particle-tracking job from a vector file or cell IDs.

    ``prepare_particle_tracking`` defaults to ``execute=False`` and
    ``convert_output=False`` so it can be used as a setup-first workflow before the
    final ``tracker.run(...)`` call. It still runs writeP3DGSF, because
    :meth:`ParticleTrackingInput.run` always makes the GSF file.

    Parameters
    ----------
    model : SimulationBase
        The flow model (``SimulationBase``-like: ``name``,
        ``model_output_folder_path``, ``vor`` and ``gwf``). See
        :class:`ParticleTrackingInput`.
    particles : str, pathlib.Path, or Sequence[int]
        The release source. A ``str`` or ``Path`` is a ``.shp``/``.gpkg`` file,
        passed as ``particle_shp``; anything else is iterated as zero-based cell
        (CELL2D) numbers, passed as ``particle_cells``.
    porosity : float or Sequence[float], default 0.2
        An ``int``/``float`` is repeated for every model layer
        (``model.gwf.modelgrid.nlay``); anything else is iterated, one value per
        layer. Becomes ``porosities_by_layer`` (a ``porosities_by_layer`` passed in
        ``kwargs`` replaces it).
    output_path : str or pathlib.Path, optional
        Workspace for MP3DU inputs and outputs. Defaults to
        ``model.model_output_folder_path``.
    execute : bool, default False
        Passed to :meth:`ParticleTrackingInput.run`: run mod-PATH3DU after writing
        its inputs.
    convert_output : bool, default False
        Passed to :meth:`ParticleTrackingInput.run`: with ``execute``, convert the
        binary output to shapefiles and summarize the endpoints.
    write_diagnostics : bool, default True
        Passed to :meth:`ParticleTrackingInput.run`: write
        ``mp3du_diagnostics.json``.
    **kwargs
        Forwarded to the :class:`ParticleTrackingInput` constructor -- any of its
        keywords, such as ``direction``, ``simulation_end_time``,
        ``iface_overrides``, ``generated_particle_zloc``, the executable paths and
        the threading controls. A ``particle_shp``/``particle_cells`` here of the
        same kind as ``particles`` is overwritten by it; one of the other kind is
        kept beside it, and the run then raises ``ValueError`` for two particle
        sources. An unknown keyword raises ``TypeError``.

    Returns
    -------
    tuple[ParticleTrackingInput, ParticleTrackingResult]
        The configured tracker and the structured result from ``run()``.
    """

    if isinstance(porosity, (int, float)):
        porosities_by_layer = [float(porosity) for _ in range(model.gwf.modelgrid.nlay)]
    else:
        porosities_by_layer = [float(value) for value in porosity]

    tracker_kwargs = dict(
        model=model,
        output_path=None if output_path is None else Path(output_path),
        porosities_by_layer=porosities_by_layer,
    )
    tracker_kwargs.update(kwargs)

    if isinstance(particles, (str, Path)):
        tracker_kwargs["particle_shp"] = Path(particles)
    else:
        tracker_kwargs["particle_cells"] = [int(cell) for cell in particles]

    tracker = ParticleTrackingInput(**tracker_kwargs)
    result = tracker.run(
        execute=bool(execute),
        convert_output=bool(convert_output),
        write_diagnostics=bool(write_diagnostics),
    )
    return tracker, result


def run_particle_tracking(
    *,
    model,
    particles: Path | str | list[int] | tuple[int, ...],
    porosity: float | list[float] = 0.2,
    output_path: Path | str | None = None,
    execute: bool = True,
    convert_output: bool = True,
    write_diagnostics: bool = True,
    **kwargs,
) -> ParticleTrackingResult:
    """Run MP3DU end-to-end from either a vector selection or explicit cell IDs.

    This is the shortest public entry point for common use:

    ``run_particle_tracking(model=model, particles=Path('my.gpkg'))``
    ``run_particle_tracking(model=model, particles=[9151, 9152, 9747])``

    It is :func:`prepare_particle_tracking` with ``execute`` and ``convert_output``
    defaulting to ``True``, returning only the result.

    Parameters
    ----------
    model : SimulationBase
        The flow model. See :class:`ParticleTrackingInput`.
    particles : str, pathlib.Path, or Sequence[int]
        A ``.shp``/``.gpkg`` path of release locations, or zero-based cell (CELL2D)
        numbers. See :func:`prepare_particle_tracking`.
    porosity : float or Sequence[float], default 0.2
        One porosity for every layer, or one per layer.
    output_path : str or pathlib.Path, optional
        Workspace for MP3DU inputs and outputs. Defaults to
        ``model.model_output_folder_path``.
    execute : bool, default True
        Run mod-PATH3DU after writing its inputs (requires the ``mp3du``
        executable).
    convert_output : bool, default True
        With ``execute``, convert the binary output to shapefiles with
        writeP3DOutput and summarize the endpoints.
    write_diagnostics : bool, default True
        Write ``mp3du_diagnostics.json``.
    **kwargs
        Forwarded through :func:`prepare_particle_tracking` to the
        :class:`ParticleTrackingInput` constructor; an unknown keyword raises
        ``TypeError``.

    Returns
    -------
    ParticleTrackingResult
        The structured result of :meth:`ParticleTrackingInput.run`.
    """

    _, result = prepare_particle_tracking(
        model=model,
        particles=particles,
        porosity=porosity,
        output_path=output_path,
        execute=execute,
        convert_output=convert_output,
        write_diagnostics=write_diagnostics,
        **kwargs,
    )
    return result


# Legacy FloPy/MODFLOW PRT helpers moved to mp3du.legacy_prt; the supported
# MP3DU API is ParticleTrackingInput / prepare_particle_tracking /
# run_particle_tracking. Warned + hidden from dir()/completion per D12.
_DEPRECATED = {
    name: (
        f"myflopy.modflow.mp3du.legacy_prt:{name}",
        f"myflopy.modflow.mp3du.legacy_prt.{name}",
        "0.1.0",
    )
    for name in _LEGACY_PRT_NAMES
}
__getattr__, __dir__ = deprecated_module_getattr(_DEPRECATED, __name__)
