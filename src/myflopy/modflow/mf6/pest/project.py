"""Project orchestration for the first ``myflopy`` PEST slice."""

from __future__ import annotations

import contextlib
import io
import json
import os
import re
import sys
import warnings
from pathlib import Path

from myflopy._optional import require
from myflopy.modflow.mf6.observations import (
    ConcTargets,
    DrnFlowTargets,
    HeadTargets,
    LakeStageTargets,
    SfrFlowTargets,
    SfrStageTargets,
)
from myflopy.modflow.mf6.period_end import period_lengths
from myflopy.modflow.mf6.pest.geostats import build_geostruct
from myflopy.modflow.mf6.pest.native_parameters import (
    NativeParameterSpec,
    add_native_parameter,
    relayer_array_target,
)
from myflopy.modflow.mf6.pest.observations import (
    finalize_observations,
    prepare_conc_observations,
    prepare_drn_flow_observations,
    prepare_head_target_observations,
    prepare_lake_stage_observations,
    prepare_sfr_flow_observations,
    prepare_sfr_stage_observations,
)
from myflopy.modflow.mf6.pest.pilot_points import (
    add_pilot_point_parameter,
    register_pilot_point_parameters,
)
from myflopy.modflow.mf6.pest.runs import default_pest_root
from myflopy.modflow.mf6.pest.specs import (
    ConcObservationSpec,
    DrnFlowObservationSpec,
    ExpGeoStruct,
    HeadTargetObservationSpec,
    LakeStageObservationSpec,
    SfrFlowObservationSpec,
    SfrStageObservationSpec,
)
from myflopy.modflow.mf6.pest.summary import PestSettings

METADATA_FILENAME = "myflopy_pest_metadata.json"

# Pre-built observation specs accepted by observe()/forecast() as-is.
_OBSERVATION_SPEC_TYPES = (
    HeadTargetObservationSpec,
    ConcObservationSpec,
    LakeStageObservationSpec,
    SfrStageObservationSpec,
    SfrFlowObservationSpec,
    DrnFlowObservationSpec,
)

# High-level target sets and the observation spec each is wrapped in, so
# observe()/forecast() accept head, lake-stage, SFR stage/flow and DRN seepage
# targets uniformly (each spec carries its own default prefix).
_TARGET_SPEC_TYPES = (
    # ConcTargets subclasses HeadTargets, so it MUST come first -- an isinstance
    # ladder would otherwise wrap every concentration target as a head spec.
    (ConcTargets, ConcObservationSpec),
    (HeadTargets, HeadTargetObservationSpec),
    (LakeStageTargets, LakeStageObservationSpec),
    (SfrStageTargets, SfrStageObservationSpec),
    (SfrFlowTargets, SfrFlowObservationSpec),
    (DrnFlowTargets, DrnFlowObservationSpec),
)


def _pest_run_slug(name: str) -> str:
    """Filesystem-safe directory name for a PEST run's default workspace."""

    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name)).strip("_")
    return slug or "pest_run"


def build_forward_run_command(template_workspace: Path, interpreter: str) -> str:
    """Return the PEST++ model command that runs ``forward_run.py``.

    The forward-run interpreter is pinned to ``interpreter`` (this
    environment's Python) so PEST++ workers use the venv that has
    numpy/flopy/pyemu/myflopy, not a bare ``python`` that may resolve
    elsewhere. On Windows the quoted absolute path works directly. On POSIX,
    PEST++'s run manager neither shell-parses quotes nor accepts absolute
    command paths (it mangles the leading ``/`` and ``execv`` fails), so the
    command is a ``./run_forward.sh`` wrapper written into the template
    workspace that ``exec``s the absolute interpreter; PEST++ worker dirs are
    copies of the template, so the wrapper travels with them.
    """

    if os.name == "nt":
        return f'"{interpreter}" forward_run.py'
    wrapper = Path(template_workspace) / "run_forward.sh"
    wrapper.write_text(f'#!/bin/sh\nexec "{interpreter}" forward_run.py "$@"\n')
    wrapper.chmod(0o755)
    return "./run_forward.sh"


def _import_pyemu():
    """Import pyEMU lazily with a workflow-oriented error message."""

    return require("pyemu", feature="PestProject workflows")


def _pyemu_warning_class(pyemu_module):
    """Return the installed ``PyemuWarning`` class if available."""

    try:
        return pyemu_module.pyemu_warnings.PyemuWarning
    except AttributeError:
        return Warning


class PestProject:
    """Build a pyEMU/PEST++ calibration workspace from a ``myflopy`` model.

    It compiles a few readable declarations straight to native
    ``pyemu.utils.PstFrom`` and lets pyEMU drive the forward run.

    **The front door is** ``model.pest("calib")``, which returns one of these
    with ``start_datetime`` filled in from the model's TDIS. Constructing
    ``PestProject(...)`` directly is the advanced path: everything else is the
    same -- ``workspace`` still defaults to ``<model workspace>.pest/<name>``, so
    the run is still discoverable through ``model.pest_runs`` -- but
    ``start_datetime`` has no TDIS fallback and you must supply it::

        cal = model.pest("calib", start_datetime="2020-01-01")
        cal.parameterize("k",        style="pilotpoints", pp_space=8, physical=(1e-3, 100))
        cal.parameterize("recharge", bounds=(0.5, 1.5))
        cal.parameterize("ghb.cond", bounds=(0.1, 10))
        cal.observe(head_targets)          # history-matching targets
        cal.forecast(prediction_targets)   # predictions of interest (zero weight)
        print(cal.settings())              # review the resolved configuration
        pst = cal.build("calib.pst")       # -> native .pst + forward_run.py
        cal.run_ies()                      # PESTPP-IES ensemble smoother

    Use :meth:`parameterize` (styles ``constant`` / ``zone`` / ``grid`` /
    ``pilotpoints`` across every target), :meth:`observe`, :meth:`forecast`,
    :meth:`build`, :meth:`run_ies`/:meth:`prior` and :meth:`settings`. The runs
    land beside the model and are reopened for review via ``model.pest_runs``
    (or ``run.pest_runs``) -> :meth:`~...runs.PestRunHandle.review`
    (:class:`~...ies.IesResults`).

    Parameters
    ----------
    model : SimulationBase
        A built ``myflopy`` model (``SimulationBase`` or compatible) exposing
        ``.sim`` (the FloPy simulation), ``.gwf``, ``.name``,
        ``.workspace`` and, for ``pp_space=`` pilot points, any ``zones=`` and
        head or concentration targets given by cell id, ``.vor``. Its workspace
        is what ``PstFrom``
        copies: :meth:`build` first REWRITES the simulation there with every
        array and list stored in external files
        (``set_all_data_external``), and observation targets read the model's
        current output at build time.
    name : str
        Project name (``str()`` is applied). It names the default template
        directory (each run of characters outside ``A-Za-z0-9_.-`` becomes one
        ``_``, leading and trailing ``_`` are stripped), the
        default ``.pst`` stem, the ``<name>.pst`` that :meth:`run_ies` and
        :meth:`prior` always write and run, their master directories
        (``<name>_ies_master`` / ``<name>_prior_master``), and the run as
        ``model.pest_runs`` lists it.
    workspace : str or pathlib.Path or None, optional
        Directory for the generated PEST template workspace. Default ``None``
        uses ``<model workspace>.pest/<name>``, a sibling of the model
        directory where ``model.pest_runs`` finds it. It is DELETED and
        recreated when :meth:`build` runs. A directory inside the model
        workspace is refused by :meth:`build` (``RuntimeError``) unless its
        parent directory is empty once it is removed, because ``PstFrom``
        would copy the workspace into itself.
    start_datetime : str or pandas.Timestamp
        Simulation start date (e.g. ``"2020-01-01"``), passed to
        ``PstFrom(start_datetime=...)``, which parses it to a datetime to place
        time-varying parameters and observations on the time axis. Required
        here, keyword-only; ``model.pest(name)`` fills it from the model's
        TDIS ``start_date_time``.
    spatial_reference : object, optional
        Cell spatial reference passed to ``PstFrom(spatial_reference=...)``,
        which ``grid``-style parameters and ``correlation=`` geostatistics
        locate cells with. Default ``None`` uses ``model.gwf.modelgrid`` at
        :meth:`build`.
    zero_based : bool, optional
        Passed to ``PstFrom(zero_based=...)``: whether the index columns of the
        model's list files count from zero. Default ``False``, because MF6
        writes its ``layer cell`` columns one-based (pyEMU's own default is
        ``True``).
    longnames : bool, optional
        Passed to ``PstFrom(longnames=...)``: allow parameter and observation
        names longer than classic PEST's limits. Default ``True``.
    """

    def __init__(
        self,
        model,
        name: str,
        *,
        workspace: str | Path | None = None,
        start_datetime: str,
        spatial_reference=None,
        zero_based: bool = False,
        longnames: bool = True,
    ):
        """Create a calibration project around an existing model workspace.

        ``workspace`` is the PEST template directory. When omitted it defaults to
        ``<model workspace>.pest/<name>`` -- a SIBLING of the model directory, so
        the calibration lives beside the model it calibrates, is
        auto-discoverable via ``model.pest_runs`` / ``run.pest_runs``, and is
        not inside the tree ``PstFrom`` copies (ledger 107). See the class
        docstring for the full reference.
        """

        self.model = model
        self.name = str(name)
        self.original_workspace = Path(model.workspace)
        if workspace is None:
            # A SIBLING of the model directory, not a child: `PstFrom` copies the
            # model workspace wholesale via a bare `shutil.copytree` in a private
            # function, with no way to exclude a subtree, so a template inside it
            # gets copied into itself (ledger 107).
            workspace = default_pest_root(self.original_workspace) / _pest_run_slug(name)
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

    # -- declarative facade -----------------------------------------------
    #
    # These methods compile directly to native ``pyemu.utils.PstFrom`` calls
    # (see native_parameters.py): ``parameterize``/``observe``/``forecast``
    # declare the calibration, ``build`` writes the ``.pst`` and forward run.

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
        anisotropy: float = 1.0,
        bearing: float = 0.0,
        nugget: float = 0.0,
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
        target : str
            Friendly name of the model property to calibrate (case and
            surrounding spaces ignored). Supported targets
            (with aliases): ``"k"`` (``"npf.k"``/``"hk"``/``"kh"``),
            ``"k33"`` (``"kv"``), ``"recharge"`` (``"rch"``),
            ``"chd"``, ``"ghb.cond"``/``"ghb.bhead"`` (``"ghb"``),
            ``"drn.cond"``/``"drn.elev"`` (``"drn"``), ``"wel"`` (``"pumping"``),
            ``"uzf.vks"`` (``"vks"``) — UZF saturated vertical K, the one UZF
            quantity heads respond to — and ``"porosity"``
            (``"mst.porosity"``/``"n"``), which lives on the GWT sibling of a
            coupled transport simulation and is resolved there automatically.
            ``k``, ``k33`` and ``porosity`` are ARRAY targets; the rest scale a
            column of a LIST package's files. An unknown name raises
            ``KeyError`` here, listing the known ones. Each target scales the
            external files FloPy writes for it (``recharge``, for instance,
            reads ``<model>.rch_stress_period_data_*.txt``), and :meth:`build`
            raises ``FileNotFoundError`` when the model has none.

            Porosity is absent from the flow equation, so it is identifiable
            only from CONCENTRATION data — and since transport velocity is
            ``v = Ki/n``, K and porosity are near-collinear from concentration
            alone (measured cosine 0.98 on the canonical model). Estimate both
            only with heads in the mix: heads pin K, and concentration then pins
            porosity.
        style : str or None, optional
            Spatial parameterization style:

            - ``"constant"`` -- a single multiplier for the whole target;
            - ``"zone"`` -- one multiplier per ``zones`` value (pass ``zones=``);
            - ``"grid"`` -- one multiplier per list entry / cell;
            - ``"pilotpoints"`` -- geostatistical pilot points.

            ``None`` uses the recipe default (``"constant"`` for every target).
            The value is lower-cased and not checked against these four: it
            reaches pyEMU as ``par_type``, except ``"pilotpoints"`` on an array
            target, which myflopy builds itself. ``"grid"`` on
            *array* targets (K, K33) builds one multiplier per Voronoi/DISV
            cell, and an independent set per layer file, located through the
            model's spatial reference (the modelgrid is used automatically);
            pass ``correlation=`` to correlate them and ``layers=`` to restrict
            which layers. ``"grid"`` also works for list targets.
            ``"pilotpoints"`` on an array target places a pilot-point net
            (``pp_space=`` spacing or explicit ``pp_points=``) and interpolates
            the multipliers to cells by inverse-distance weighting in the
            forward run -- the supported sparse-field style on Voronoi/DISV
            grids.
        bounds : tuple of (float, float), optional
            ``(lower, upper)`` bounds on the adjustable value, pyEMU's
            ``lower_bound``/``upper_bound`` (each pilot point's ``parlbnd`` /
            ``parubnd``). Default ``(0.5, 2.0)``. For a multiplier
            these are factors (e.g. ``(0.2, 5)`` = 5x down to 5x up); for an
            additive parameter they are offsets in model units.
        physical : tuple of (float, float) or None, optional
            ``(ult_lbound, ult_ubound)`` -- hard limits clamped onto the *final*
            model value after all multipliers are applied (pyEMU's
            ``ult_lbound``/``ult_ubound``; the pilot-point interpolation clamps
            to them itself). Strongly recommended
            for multipliers so calibration cannot push K (etc.) to nonphysical
            values. ``None`` (default) leaves the input unclamped.
        transform : str, optional
            ``"log"`` (default, recommended for strictly positive properties
            like K and recharge) or ``"none"``; pyEMU also accepts
            ``"fixed"``, but not together with ``correlation=`` (the variogram
            conversion asserts ``"none"`` or ``"log"``). Forced to ``"none"``
            for additive parameters.
        additive : bool or None, optional
            Apply the parameter as an additive offset instead of a multiplier
            (pyEMU ``par_style="a"``). ``None`` (default) uses the recipe
            default (drain elevation is additive; everything else is
            multiplicative). Array-target pilot points still MULTIPLY in their
            forward run; for them it only forces the transform to ``"none"``.
        zones : ZoneSpec or array_like of int or None, optional
            Per-cell integer zone labels, normalized to the shape pyEMU wants
            for the target's family at :meth:`build` (see :class:`ZoneSpec`,
            which also reads a raster or polygon layer) and passed as pyEMU's
            ``zone_array``. With ``style="zone"`` each distinct label gets one
            parameter. With any style, an array target leaves cells labelled
            ``<= 0`` unparameterized; a list target makes zone ``0`` a real
            parameter and warns. A length that is not one label per cell raises
            ``ValueError`` at build. Default ``None`` passes no zones.
            Array-target pilot points ignore it.
        layers : iterable of int or None, optional
            For multi-layer *array* targets, restrict the parameter to these
            zero-based model layers (e.g. ``layers=[0, 1]`` for the unconfined
            aquifers); other layers keep their input value. ``None`` (default)
            covers all layers. Read only for array targets -- a list target
            ignores it. :meth:`build` raises ``ValueError`` when no requested
            layer has its own file, or the target is written as one whole-grid
            file.
        pp_space : int or None, optional
            Pilot-point spacing for ``style="pilotpoints"`` on an array target,
            in multiples of the median cell width (the square root of the median
            ``model.vor`` cell area): the grid is binned at that spacing and the
            cell centre nearest each occupied bin's middle becomes a point.
            Truncated with ``int``; ``None`` (default) or ``0`` means ``5``.
            Ignored when ``pp_points`` is given.
        pp_points : geopandas.GeoDataFrame or array_like or None, optional
            Explicit pilot-point locations for ``style="pilotpoints"`` on an
            array target: a GeoDataFrame of Point geometries, or ``(x, y)``
            pairs (anything reshapable to ``(n, 2)``), in the model grid's own
            coordinates -- a GeoDataFrame's CRS is not read or reprojected.
            Takes precedence over ``pp_space``. Default ``None``.
        correlation : float or None, optional
            Variogram range (model length units) of an exponential variogram
            (see :class:`~myflopy.modflow.mf6.pest.specs.ExpGeoStruct`). With
            ``style="grid"`` it becomes the parameters' pyEMU ``geostruct``, and
            :meth:`run_ies`/:meth:`prior` draw a spatially correlated prior
            ensemble from it (``ies_par_en``). With ``style="pilotpoints"`` on an
            array target it shapes only that prior draw -- the pilot values are
            drawn correlated, while the interpolation to cells stays
            inverse-distance. ``None`` (default) builds no variogram; ignored
            for ``"constant"``/``"zone"``.
        anisotropy : float, optional
            Anisotropy ratio of the variogram -- correlated this many
            times farther along the major axis than across it (``1.0`` =
            isotropic, the default). For an alluvial valley, e.g. ``5`` makes K
            correlate farther down-valley than across it. Acts only when
            ``correlation`` is given; pyEMU asserts it is positive.
        bearing : float, optional
            Azimuth of the anisotropy major axis, in degrees east of north.
            Default ``0.0``. Ignored when ``anisotropy`` is ``1.0``, and
            without ``correlation``.
        nugget : float, optional
            Nugget (unresolved short-scale variance) of the variogram, added to
            the covariance diagonal on top of the variogram's contribution,
            which is fixed at ``1.0`` here (so the total sill is
            ``1.0 + nugget``). Default ``0.0``. Acts only when ``correlation``
            is given.
        temporal : float or None, optional
            Temporal correlation range (days) for time-varying list packages.
            Recorded on the returned spec and read by nothing yet, so it has no
            effect on the build. Default ``None``.
        name : str or None, optional
            Parameter group / name base (pyEMU's ``par_name_base`` and
            ``pargp``). Default ``None`` uses the canonical target with dots
            removed (``"k"``, ``"ghbcond"``, ``"uzfvks"``). A multi-layer
            array ``"grid"`` parameter appends ``l1``, ``l2``, ... per selected
            layer file, pilot points append ``l<layer>`` (zero-based), and a
            captured field is named ``<name>field``. Give each call a distinct
            name when stacking calls on one target: the project's
            per-parameter records are keyed by it, so a second call under the
            same name overwrites the first's.
        capture : bool, optional
            Also record the *resolved* model-input field (after multipliers are
            applied) as zero-weight observations named from ``<name>field``
            (one prefix per layer file for a layered array), so each
            realization's per-cell property field is carried in the ensembles.
            Enables spatial parameter maps via ``IesResults.plot_field`` /
            ``IesResults.field``. Default ``False``.
        **extra
            Advanced keyword arguments passed straight through to pyEMU's
            ``PstFrom.add_parameters`` for cases the facade does not cover.
            They are applied LAST, so a key here overrides what the facade
            derived (``par_type``, the bounds, ``transform``, ``zone_array``,
            ``geostruct``, ...). Not forwarded for array-target pilot points,
            which do not go through ``add_parameters``. Stored on the returned
            spec as ``extra``.

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
            anisotropy=anisotropy,
            bearing=bearing,
            nugget=nugget,
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

        What the control file holds is exactly the targets (ledger 340). Each is
        matched to the simulated row whose index EQUALS its ``time`` read as a
        zero-based stress period (``1.0`` is period 1, ``1.5`` is refused), and a
        target no row holds -- a period the output did not save, a name with no
        simulated column, a head or concentration target on no active cell
        (ledger 341) -- makes :meth:`build` raise, naming it. A blank value
        is unmeasured. Every simulated row no target measures stays in the
        control file at weight 0, its ``obsval`` the base run's own value, so it
        is reported but never history-matched; it used to carry weight 1, which
        pulled the calibration toward the uncalibrated model at every period
        nobody measured. A set given no values at all registers every row at
        weight 0: its simulated values are captured, nothing is fitted.

        A FLOW is written to the control file in MF6's own sign, on both routes
        (ledger 350): an SFR reach's streamflow and a DRN zone's discharge are
        negative, water leaving the reach or the aquifer, so a gauged ``Q`` is
        entered as ``-Q`` -- the targets refuse a value above zero when they are
        built, which scored ``-2Q`` at a perfect fit.

        Parameters
        ----------
        targets : target set or observation spec
            A high-level target set -- :class:`HeadTargets`,
            :class:`ConcTargets`, :class:`LakeStageTargets`,
            :class:`SfrStageTargets`, :class:`SfrFlowTargets` or
            :class:`DrnFlowTargets` -- which is wrapped in its matching
            ``*ObservationSpec``, or a pre-built observation spec (e.g.
            :class:`HeadTargetObservationSpec`), registered as it is. The
            matching simulated values are regenerated by a post-processor on
            every run. Anything else raises ``TypeError``.
        prefix : str or None, optional
            Observation-name prefix. ``None`` (default) uses the target kind's
            own prefix (``"hds"``, ``"conc"``, ``"stage"``, ``"sfr_stage"``,
            ``"sfr_flow"``, ``"drn_flow"``). Pass distinct prefixes to register
            more than one set of the same kind. IGNORED when ``targets`` is
            already a spec: the spec's own ``prefix`` is used.

        Returns
        -------
        observation spec
            The registered spec -- the one passed in, or the
            ``*ObservationSpec`` (e.g. :class:`HeadTargetObservationSpec`)
            the target set was wrapped in.
        """

        spec = self._coerce_observation_spec(targets, prefix=prefix)
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

        A forecast is the rows its targets NAME, matched as :meth:`observe`
        matches them (a target no row holds makes :meth:`build` raise); its value
        may be blank. The set's other rows -- the periods nobody forecast -- are
        zero-weighted too and are not forecasts; they used to stay at weight 1,
        history-matched to the base run's own values (ledger 340). A set given
        NO values forecasts every row of its table: its locations at every
        period.

        Parameters
        ----------
        targets : target set or observation spec
            Any target set :meth:`observe` accepts (:class:`HeadTargets`,
            :class:`ConcTargets`, :class:`LakeStageTargets`,
            :class:`SfrStageTargets`, :class:`SfrFlowTargets`,
            :class:`DrnFlowTargets`) or a pre-built observation spec,
            describing *where* and *when* the prediction is taken. Any values
            supplied are placeholders (e.g. an expected value), written as the
            observation's value; a blank keeps the base run's. The weight is
            forced to zero regardless. Anything else raises ``TypeError``.
        prefix : str or None, optional
            Observation-name prefix. ``None`` (default, and an empty string
            too) gives ``"fore1"``, ``"fore2"``, ... in call order, so forecast
            names never collide with calibration observations. IGNORED when
            ``targets`` is already a spec: the spec's own ``prefix`` is used,
            which for a default-built spec is the same as a calibration set's
            (``"hds"``, ...).

        Returns
        -------
        observation spec
            The registered forecast spec -- the one passed in, or the
            ``*ObservationSpec`` the target set was wrapped in.
        """

        default = prefix or f"fore{len(self._forecast_specs) + 1}"
        spec = self._coerce_observation_spec(targets, prefix=default)
        self._forecast_specs.append(spec)
        return spec

    def _coerce_observation_spec(self, targets, *, prefix):
        """Turn a target set (or spec) into a registered observation spec.

        Accepts a pre-built observation spec, or any high-level target set
        (head / lake stage / SFR stage / SFR flow / DRN seepage), wrapping it in
        the matching spec. When ``prefix`` is ``None`` the spec's own default
        prefix is used (``"hds"``, ``"stage"``, ``"sfr_flow"``, ...); pass distinct
        prefixes to register more than one set of the same kind.
        """

        if isinstance(targets, _OBSERVATION_SPEC_TYPES):
            return targets
        for target_type, spec_type in _TARGET_SPEC_TYPES:
            if isinstance(targets, target_type):
                if prefix is None:
                    return spec_type(targets=targets)
                return spec_type(targets=targets, prefix=prefix)
        raise TypeError(
            "observe()/forecast() accept a head/lake/SFR/DRN target set or a "
            f"pre-built observation spec; got {type(targets).__name__}."
        )

    def _ensure_original_workspace(self):
        """Ensure MF6 input files exist before creating the PEST template."""

        self.original_workspace.mkdir(parents=True, exist_ok=True)
        self.model.sim.write_simulation(silent=True)

    def _clear_stale_template(self):
        """Make sure the PEST template directory is not copied into itself.

        ``PstFrom`` copies ``original_workspace`` wholesale into
        ``template_workspace``. If the template lives INSIDE the workspace being
        copied, ``shutil.copytree`` descends into the destination it is
        currently filling -- ``pest/<name>/pest/<name>/pest/...`` -- until it
        dies of recursion. Measured 2026-07-28: re-running a calibration cell
        left a 17 GB tree 40 levels deep. ``PstFrom``'s own
        ``remove_existing=True`` does not help -- it clears the destination,
        which is not what recurses. Nor is deleting just this run's template
        enough: an EMPTY subdirectory still gets walked.

        **The default no longer lands there** (ledger 107): it is now
        ``<model workspace>.pest/<name>``, a sibling of the copied tree, so this
        returns early and several named calibrations coexist freely. What
        remains is the guard for an explicit ``workspace=`` that the caller put
        inside the model directory. There it removes this run's template and, if
        that leaves the parent empty, removes the parent too. Other calibrations
        are never deleted -- their IES master directories hold finished results
        -- so if any exist the copy cannot be made safe here and this raises
        rather than exploding later.
        """

        import shutil

        template = self.template_workspace
        if template.exists():
            shutil.rmtree(template, ignore_errors=True)

        pest_root = template.parent
        try:
            inside = pest_root.resolve().is_relative_to(
                self.original_workspace.resolve()
            )
        except (OSError, ValueError):  # pragma: no cover - defensive
            inside = False
        if not inside or not pest_root.is_dir():
            return

        siblings = sorted(entry.name for entry in pest_root.iterdir())
        if not siblings:
            # Empty now -- remove it so the copy sees no `pest/` at all.
            shutil.rmtree(pest_root, ignore_errors=True)
            return

        raise RuntimeError(
            f"Cannot build this calibration at {template}: it sits inside the "
            f"model workspace that PstFrom copies, and {pest_root} still holds "
            f"{len(siblings)} other PEST run(s) ({', '.join(siblings[:5])}"
            f"{' ...' if len(siblings) > 5 else ''}). Copying would nest the "
            "workspace inside itself until it exhausts the disk.\n\n"
            "Drop the explicit workspace= to use the default, which sits beside "
            f"the model rather than inside it:\n"
            f"    model.pest({self.name!r}, ...)   # -> "
            f"{default_pest_root(self.original_workspace) / _pest_run_slug(self.name)}\n"
            "or point workspace= anywhere outside the model directory."
        )

    def _build_pstfrom(self):
        """Instantiate the underlying ``pyemu.utils.PstFrom`` object."""

        self._clear_stale_template()
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

        if isinstance(spec, ConcObservationSpec):
            return prepare_conc_observations(self, spec)
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

    def _refuse_shared_prefixes(self):
        """Refuse two observation or forecast sets under one prefix, before anything is written.

        The prefix names every file a set writes, so the second overwrote the
        first's table, instruction file and map, and the control file the build
        handed back could not run -- nor could its targets be told apart by the
        positional matching of ledger 340.
        """

        prefixes = [
            str(spec.prefix).strip().lower()
            for spec in [*self._observation_specs, *self._forecast_specs]
        ]
        repeated = sorted({prefix for prefix in prefixes if prefixes.count(prefix) > 1})
        if repeated:
            raise ValueError(
                f"Two observation or forecast sets share the prefix {repeated[0]!r}"
                + (f" (and {', '.join(map(repr, repeated[1:]))})" if len(repeated) > 1 else "")
                + ". The prefix names every file a set writes, so the second would "
                "overwrite the first. Pass observe(targets, prefix=...) a distinct "
                "prefix per set."
            )

    def _prepare_observation_specs(self):
        """Create simulated observation files and register them with pyEMU."""

        self._refuse_shared_prefixes()
        prepared: list[dict] = []
        for spec in self._observation_specs:
            prepared.append(self._prepare_one_observation_spec(spec))
        for spec in self._forecast_specs:
            item = self._prepare_one_observation_spec(spec)
            item["is_forecast"] = True
            prepared.append(item)
        self._prepared_observations = prepared

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
            "pst_file": Path(filename).name if filename is not None else None,
            "observation_sets": observation_sets,
            "parameter_sets": parameter_sets,
            "capture_fields": list(self._capture_field_specs.values()),
        }
        metadata_path = self.template_workspace / METADATA_FILENAME
        metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        return metadata_path

    def write(self, filename: str | Path | None = None):
        """Write the built control file to disk.

        Writes ``cal.pst`` -- including any edits made to it since
        :meth:`build` -- in the format pyEMU's ``Pst.write`` picks when given
        no version: version 2 when there are more than 10,000 parameters or
        observations or a version-2-only column is present, else version 1.
        :meth:`build`, :meth:`run_ies` and :meth:`prior` always write
        version 2.

        Parameters
        ----------
        filename : str or pathlib.Path or None, optional
            Control-file name. Only its final component is used: the file
            always lands in the template workspace. ``None`` (default) writes
            ``<name>.pst``, the file :meth:`run_ies`/:meth:`prior` run.

        Returns
        -------
        pathlib.Path
            The file written.

        Raises
        ------
        ValueError
            When :meth:`build` has not run yet.
        """

        if self.pst is None:
            raise ValueError("Call build() before write().")
        target = self.template_workspace / (Path(filename).name if filename else f"{self.name}.pst")
        self.pst.write(target)
        return target

    def draw_prior(self, num_reals: int = 1000, use_specsim: bool = True):
        """Draw a prior parameter ensemble using the underlying ``PstFrom``.

        Calls ``cal.pf.draw(num_reals=..., use_specsim=...)``, which draws every
        parameter in ``cal.pst`` from its bounds, correlated through any
        ``correlation=`` geostructs registered with ``PstFrom``. Array-target
        pilot points are added to ``cal.pst`` after ``PstFrom`` builds and have
        no geostruct there, so this draw includes them but draws each one
        independently; their ``correlation=`` is honoured only by the separate
        draw :meth:`run_ies`/:meth:`prior` make.

        Parameters
        ----------
        num_reals : int, optional
            Number of realizations to draw. Default ``1000``.
        use_specsim : bool, optional
            Use pyEMU's spectral simulation for grid-scale parameters. Default
            ``True`` here, where pyEMU's own default is ``False``;
            :meth:`run_ies`/:meth:`prior` call this with ``False``. Measured on
            a DISV (``VertexGrid``) model with a ``grid`` parameter, ``True``
            raises ``AttributeError`` from inside pyEMU; pass ``False`` there.

        Returns
        -------
        pyemu.ParameterEnsemble or None
            The drawn ensemble, one row per realization; ``None`` (pyEMU's
            answer) when the control file has no adjustable parameter.

        Raises
        ------
        ValueError
            When :meth:`build` has not run yet.
        """

        if self.pf is None:
            raise ValueError("Call build() before draw_prior().")
        return self.pf.draw(num_reals=num_reals, use_specsim=use_specsim)

    # -- native build path ------------------------------------------------

    def _drop_unused_multiplier_apply(self):
        """Remove pyEMU's multiplier-apply command when there are no multipliers.

        ``PstFrom.__init__`` inserts ``apply_list_and_array_pars`` into
        ``pre_py_cmds`` unconditionally (pst_from.py:296-302), but only WRITES
        the ``mult2model_info.csv`` manifest it reads when ``par_dfs`` is
        non-empty (pst_from.py:757-762). A calibration whose only parameter is
        pilot points never calls ``pf.add_parameters`` -- pilot points register
        through template files after ``build_pst`` -- so the manifest is absent
        and every forward run dies with ``FileNotFoundError:
        mult2model_info.csv`` from inside pyEMU, naming nothing the user wrote.

        An empty manifest is not an option: ``apply_list_and_array_pars`` asserts
        ``ddf.shape[0] > 0`` (helpers.py:1953) and has no zero-parameter path. So
        the command is dropped instead, keyed on the same ``par_dfs`` signal
        pyEMU itself branches on.
        """

        if getattr(self.pf, "par_dfs", None):
            return
        self.pf.pre_py_cmds = [
            command for command in self.pf.pre_py_cmds
            if "apply_list_and_array_pars" not in str(command)
        ]

    def _ensure_external_model(self):
        """Write MF6 inputs as external array/list files for native PstFrom."""

        self.original_workspace.mkdir(parents=True, exist_ok=True)
        # Before externalizing, not after: a griddata array supplied as a scalar
        # writes one whole-grid file that pyEMU cannot index against a DISV
        # spatial reference. See `relayer_array_target`.
        for spec in self._native_parameter_specs:
            relayer_array_target(self, spec.recipe)
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
        """Build the pyEMU geostruct for a spatial parameter spec.

        Routes through the single :func:`~...geostats.build_geostruct` builder so
        there is one place that turns a variogram description into a pyEMU
        ``GeoStruct``. The spec's ``correlation`` is the range; ``anisotropy`` and
        ``bearing`` shape the ellipse (e.g. K correlated farther down an alluvial
        valley than across it); ``nugget`` is the unresolved short-scale variance.
        """

        return build_geostruct(
            ExpGeoStruct(
                range=float(spec.correlation),
                anisotropy=float(spec.anisotropy),
                bearing=float(spec.bearing),
                nugget=float(spec.nugget),
                transform=spec.resolved_transform,
            )
        )

    def _attach_native_observation_postprocessors(self):
        """Add post-model functions that regenerate simulated observation CSVs."""

        helper_path = Path(__file__).with_name("forward_run.py")
        # Each period's length, so the forward run can check that the record it
        # reads ENDS its period, as the base run's build did (ledger 336). OC is
        # not a parameter, but ATS can move a STEPS save off a period's end in one
        # realization and not another.
        perlen = period_lengths(self.model)
        for item in self._prepared_observations:
            head_config = item.get("forward_run_config")
            if head_config is not None:
                self.pf.add_py_function(
                    str(helper_path),
                    "_write_head_target_csv("
                    f"model_name='{self.model.name}', "
                    f"mapping_csv='{head_config['mapping_csv']}', "
                    f"output_csv='{head_config['output_csv']}', "
                    f"perlen={perlen!r})",
                    is_pre_cmd=False,
                )
                continue
            conc_config = item.get("conc_forward_run_config")
            if conc_config is not None:
                self.pf.add_py_function(
                    str(helper_path),
                    "_write_conc_target_csv("
                    f"model_name='{conc_config['model_name']}', "
                    f"mapping_csv='{conc_config['mapping_csv']}', "
                    f"output_csv='{conc_config['output_csv']}', "
                    f"perlen={perlen!r})",
                    is_pre_cmd=False,
                )
                continue
            series_config = item.get("named_series_forward_run_config")
            if series_config is not None:
                locations_file = series_config.get("locations_file")
                if locations_file is None:
                    raise ValueError(
                        f"Named-series observation set {item.get('prefix')!r} has no "
                        "saved location definition to regenerate simulated values from."
                    )
                # The calibrated model by NAME: without it the forward run reloads
                # the simulation's first GWF model, another model's lake or zone
                # in a simulation of several (ledger 340).
                self.pf.add_py_function(
                    str(helper_path),
                    "write_named_series_targets("
                    f"kind='{series_config['kind']}', "
                    f"locations_file='{locations_file}', "
                    f"output_csv='{series_config['output_csv']}', "
                    f"model_name='{self.model.name}', "
                    f"package={series_config.get('package')!r})",
                    is_pre_cmd=False,
                )
                continue
            raise NotImplementedError(
                f"Observation set {item.get('prefix')!r} has no forward-run "
                "post-processor; supported kinds are head-target and named-series "
                "(lake stage, SFR stage/flow, DRN seepage) observations."
            )

    def _forecast_observation_names(self) -> list[str]:
        """The observations the forecast sets name, as ``finalize_observations`` read them.

        Each forecast target's observation is the one the instruction file reads at
        its row and column -- not a name rebuilt from its key's spelling, which
        dropped a forecast keyed ``1.0`` where pyEMU wrote ``per:1`` (ledger 340).
        """

        names: list[str] = []
        for item in self._prepared_observations:
            if item.get("is_forecast"):
                names.extend(item.get("forecast_names", []))
        return list(dict.fromkeys(names))

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
        #
        # `layer_prefixes` stays EMPTY for a whole-grid array (MST porosity is
        # one file of nlay*ncpl values, with no LAYERED keyword). That is what
        # selects the flat `layer * ncpl + cell` reader in `IesResults._cell_map`
        # -- recording it as `{0: base}` instead would label every layer's
        # values "layer 0" and hand `plot_field` nlay*ncpl values for ncpl cells.
        layer_prefixes: dict[int, str] = {}
        for filename in spec.resolved_files:
            if recipe.family == "array":
                match = re.search(r"_layer(\d+)\.txt$", str(filename))
                prefix = f"{base}l{int(match.group(1)) - 1}" if match else base
                if match:
                    layer_prefixes[int(match.group(1)) - 1] = prefix
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

        In order it: writes the model inputs as external array/list files
        (``set_all_data_external``); creates a ``pyemu.utils.PstFrom`` over the
        model workspace; registers every :meth:`parameterize` call natively
        (``pf.add_parameters``) and every :meth:`observe`/:meth:`forecast` set;
        adds the MF6 run command and the observation post-processors to
        ``forward_run.py``; and builds the ``.pst``. The multiplier machinery is
        pyEMU's own ``apply_list_and_array_pars`` in the forward run, so parameter
        application is handled by pyEMU rather than hand-rolled code.

        After this call the workspace contains a self-contained PEST setup:
        the ``.pst``, ``forward_run.py``, template/instruction files, and the
        ``mult/`` multiplier files. Run it with PEST++ (e.g. ``pestpp-ies`` /
        ``pestpp-glm``) or via :meth:`run_ies` / :meth:`prior`.

        Parameters
        ----------
        filename : str or pathlib.Path or None, optional
            Output ``.pst`` filename (default ``"<name>.pst"``). Only its final
            component is used: it is written into the template workspace and
            recorded in the run's metadata, whose stem is the case
            ``model.pest_runs[...].review()`` opens. :meth:`run_ies` and
            :meth:`prior` always write and run ``<name>.pst`` instead, so a
            stem other than ``name`` leaves ``review()`` reading the wrong
            case.
        noptmax : int, optional
            Initial ``NOPTMAX`` written to the control file. ``0`` (default)
            means "run the model once and compute residuals" -- the standard
            way to sanity-check a fresh setup before launching a real run.
            :meth:`run_ies`/:meth:`prior` overwrite it.

        Returns
        -------
        pyemu.Pst
            The built control object (also available as ``cal.pst``). The
            underlying ``PstFrom`` is ``cal.pf`` for advanced edits.

        Raises
        ------
        ValueError
            When an observation or forecast target matches no row of its
            simulated table (naming it and the rows the table has), or a head or
            concentration target's location lands on no active cell (ledger
            341); when a target's ``time`` is not a whole number, or two values
            name one target in one period; when two sets share a prefix; or when
            a head or concentration ``simulated_values`` table's columns are not
            the forward run's, in its order (ledger 340).
        """

        target_name = Path(filename).name if filename else f"{self.name}.pst"
        self._refuse_shared_prefixes()
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
        self._drop_unused_multiplier_apply()
        self._attach_native_observation_postprocessors()
        self._write_project_metadata(filename=target_name)
        with self._quiet_pyemu_context():
            self.pst = self.pf.build_pst(filename=target_name)
        if self._pilot_point_frames:
            register_pilot_point_parameters(self)
        self.pst.model_command = [
            build_forward_run_command(self.template_workspace, sys.executable)
        ]
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
            """Summary rows for observation/forecast ``specs``, merged with prepared build info."""

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

        ``noise``, ``bad_phi_sigma`` and ``**pestpp_options`` are written onto
        ``cal.pst.pestpp_options`` and STAY there: ``noise=False``, a
        ``bad_phi_sigma`` or an extra option from an earlier
        :meth:`run_ies`/:meth:`prior` call is still set on the next one, whose
        defaults do not clear it.

        Parameters
        ----------
        reals : int, optional
            Number of realizations in the ensemble (``ies_num_reals``). More is
            better for posterior coverage but costs more model runs; 50 (the
            default) is a reasonable starting point, production runs often use
            100-300. Also the size of the correlated prior drawn for
            ``correlation=`` parameters.
        iterations : int, optional
            Number of smoother iterations (``NOPTMAX``). Default ``3``, a
            common choice; PESTPP-IES may stop earlier if it converges.
        workers : int or None, optional
            Number of parallel agents. ``None`` (default) or ``1`` runs
            serially in the template workspace, where the results then are;
            ``>1`` starts that many workers with pyEMU's ``start_workers``
            and collects the results in ``master_dir``.
        noise : bool, optional
            Whether to carry measurement noise into the posterior (recommended,
            the default ``True``, which sets nothing). ``False`` sets
            ``ies_no_noise``.
        bad_phi_sigma : float or None, optional
            Adaptive rejection threshold (``ies_bad_phi_sigma``);
            values of 1.5 (aggressive) to 2.5 (tolerant) help highly nonlinear
            problems by dropping lagging realizations. ``None`` (default) sets
            nothing.
        master_dir : str or pathlib.Path or None, optional
            Master directory for a parallel run; ``None`` (default) uses
            ``<name>_ies_master`` beside the template workspace, which is where
            ``model.pest_runs`` looks. The workers are created in its parent
            directory. Unused by a serial run.
        exe : str, optional
            PEST++ IES executable name, default ``"pestpp-ies"`` (a trailing
            ``.exe`` is dropped). When the model's MF6 executable resolves to an
            absolute path and a file ``<exe>.exe`` sits beside it, that
            directory is prepended to ``PATH`` for the run (unless ``PATH``
            already contains it); otherwise the name is resolved on ``PATH`` as
            it is.
        **pestpp_options
            Any additional ``pst.pestpp_options`` to set (e.g.
            ``ies_localizer="loc.mat"``), passed straight through and applied
            after the options above, so they override them. ``ies_par_en`` is
            then replaced by the drawn prior when a ``grid`` parameter, or an
            array-target ``pilotpoints`` parameter, carries ``correlation=``.

        Returns
        -------
        IesResults
            Reader/visualizer for the completed run, opened on the directory
            holding its results with case ``name``.

        Raises
        ------
        ValueError
            When :meth:`build` has not run, or ``cal.pst`` has been put in
            ``regularization`` mode, which PESTPP-IES cannot parse.
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

        Options persist on ``cal.pst`` between calls exactly as in
        :meth:`run_ies` (a ``bad_phi_sigma`` set there is still set here).

        Parameters
        ----------
        reals : int, optional
            Number of prior realizations (``ies_num_reals``). Default ``50``,
            plenty for a visual bracketing check. Also the size of the
            correlated prior drawn for ``correlation=`` parameters.
        workers : int or None, optional
            Parallel agents. ``None`` (default) or ``1`` runs serially in the
            template workspace; ``>1`` starts that many workers and collects the
            results in ``master_dir``.
        noise : bool, optional
            Generate the measurement-noise ensemble too (default ``True``,
            which sets nothing). ``False`` sets ``ies_no_noise``.
        master_dir : str or pathlib.Path or None, optional
            Master directory for a parallel run. ``None`` (default) uses
            ``<name>_prior_master`` beside the template workspace, which is where
            ``model.pest_runs`` looks. Unused by a serial run.
        exe : str, optional
            PEST++ IES executable, default ``"pestpp-ies"``, resolved as in
            :meth:`run_ies`.
        **pestpp_options
            Additional ``pst.pestpp_options``, passed straight through and
            applied after ``noise``; ``ies_par_en`` is replaced as in
            :meth:`run_ies`.

        Returns
        -------
        IesResults
            Use :meth:`~...ies.IesResults.plot_prior_vs_obs` (grey prior spaghetti
            vs measured) and :meth:`~...ies.IesResults.conflict` to inspect it.

        Raises
        ------
        ValueError
            When :meth:`build` has not run, or ``cal.pst`` has been put in
            ``regularization`` mode.
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

    def _refuse_regularization_mode(self) -> None:
        """Refuse a control file PESTPP-IES cannot parse, before launching it.

        Tikhonov regularization is a PESTPP-GLM concept. Applying pyEMU's
        ``zero_order_tikhonov``/``first_order_pearson_tikhonov`` to ``cal.pst``
        flips ``pestmode`` to ``regularization`` and adds prior-information
        equations, and PESTPP-IES does two things with that: it prints
        ``prior information equations not supported in ensemble methods,
        ignoring``, and -- because myflopy writes version-2 control files -- the
        ``* regularization`` keywords land in the control-data keyword block and
        the run dies with ``control file parsing error`` before a single model
        run. Measured 2026-07-30 on PEST++ 5.2.16.

        So the failure is not "regularization is ignored", it is "the run does
        not start", from an error naming a keyword the caller never typed. IES
        has its own knob for this -- ``run_ies(ies_reg_factor=...)`` -- which
        needs none of the above.
        """

        mode = str(getattr(self.pst.control_data, "pestmode", "")).strip().lower()
        if not mode.startswith("regul"):
            return
        raise ValueError(
            "This control file is in 'regularization' mode, which PESTPP-IES "
            "cannot parse in the version-2 format myflopy writes -- the run "
            "would fail with a control-file parsing error before any model run. "
            "Tikhonov regularization is a PESTPP-GLM concept; ensemble methods "
            "ignore prior-information equations entirely.\n\n"
            "For IES, the equivalent knob is a regularization FACTOR on the "
            "ensemble update:\n"
            "    cal.run_ies(reals=..., ies_reg_factor=0.1)\n"
            "and spatial structure belongs in the prior itself -- pass "
            "`correlation=` to a grid or pilot-point parameter."
        )

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

        # Grid parameters live in pyEMU's PstFrom, so `pf.draw` covers them.
        # Pilot points do NOT: they are PEST template parameters registered after
        # `build_pst`, so they need their own draw (see
        # `_draw_pilot_point_prior`).
        grid_specs = [
            spec
            for spec in self._native_parameter_specs
            if spec.style == "grid" and spec.correlation is not None
        ]
        pilot_specs = [
            spec
            for spec in self._native_parameter_specs
            if spec.style == "pilotpoints" and spec.correlation is not None
        ]
        if self.pf is None or not (grid_specs or pilot_specs):
            return

        ensemble = None
        if grid_specs:
            with self._quiet_pyemu_context():
                ensemble = self.draw_prior(int(reals), use_specsim=False)
            ensemble = ensemble._df if hasattr(ensemble, "_df") else ensemble
        for spec in pilot_specs:
            drawn = self._draw_pilot_point_prior(spec, int(reals))
            if drawn is None:
                continue
            ensemble = drawn if ensemble is None else ensemble.join(drawn, how="inner")
        if ensemble is None or ensemble.empty:
            return

        ensemble_path = self.template_workspace / f"{self.name}.prior_par.csv"
        ensemble.to_csv(ensemble_path)
        self.pst.pestpp_options["ies_par_en"] = ensemble_path.name

    def _draw_pilot_point_prior(self, spec, reals: int):
        """A geostatistically correlated prior for one pilot-point spec.

        ``correlation=`` was documented for ``style="pilotpoints"`` since the
        module was written (``pilot_points.py:14``) and read nowhere in it: the
        pilot VALUES were drawn independently from their bounds by PESTPP-IES,
        and the code's stated reasoning was that IDW interpolation supplies the
        spatial structure. It supplies smoothness, but not the CORRELATION LENGTH
        the caller asked for -- the field's structure came from pilot spacing
        instead, so a 200 m and a 2000 m variogram produced the same prior.

        Drawing the pilot values from the variogram's covariance fixes that. It
        matters most where the data are thin: with a correlated prior, IES moves
        neighbouring points together and unobserved areas keep a plausible field,
        rather than each point wandering on its own.
        """

        frames = self._pilot_point_frames.get(spec.name)
        if not frames:
            return None
        import pandas as pd

        pyemu = self.pyemu or _import_pyemu()
        geostruct = self._geostruct_for(spec)
        columns = []
        for frame, _ in frames:
            # One covariance per LAYER's point set: the layers are separate
            # parameter groups here, and correlating across them would assert a
            # vertical structure the variogram never described.
            names = [str(name) for name in frame["parnme"]]
            with self._quiet_pyemu_context():
                cov = geostruct.covariance_matrix(
                    frame["x"].to_numpy(dtype=float),
                    frame["y"].to_numpy(dtype=float),
                    names=names,
                )
                drawn = pyemu.ParameterEnsemble.from_gaussian_draw(
                    pst=self.pst, cov=cov, num_reals=int(reals), by_groups=False,
                    fill=False,
                )
            values = drawn._df if hasattr(drawn, "_df") else drawn
            columns.append(values.loc[:, [name for name in names if name in values.columns]])
        if not columns:
            return None
        return pd.concat(columns, axis=1)

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
        self._refuse_regularization_mode()
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
