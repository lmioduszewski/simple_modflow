"""Dataclass specifications for the first ``myflopy`` PEST slice."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from myflopy._defaults import DEFAULT_EPSG
from myflopy.modflow.mf6.observations import (
    ConcDifferences,
    ConcTargets,
    DrnFlowTargets,
    HeadDifferences,
    HeadTargets,
    LakeStageTargets,
    RivFlowTargets,
    SfrFlowTargets,
    SfrStageTargets,
)


@dataclass
class VectorParameterSource:
    """Describe a GIS vector source used to define calibration parameters.

    A leftover of the retired ``add_parameter``/``build_pst`` path (deleted
    2026-06-22). No calibration route takes one today: ``cal.parameterize``
    has no ``source=`` argument, and the only reader of this record,
    ``myflopy.modflow.mf6.pest.gis.load_vector_parameter_source``, is called
    by nothing in the library. It is still exported, so it is documented here
    as it stands; to calibrate from polygons use
    ``cal.parameterize(..., style="zone", zones=ZoneSpec.from_polygons(...))``.

    Parameters
    ----------
    path : str or pathlib.Path
        Vector file (shapefile, GeoPackage, ...) opened with
        ``geopandas.read_file``, which reads its first layer.
        ``load_vector_parameter_source`` refuses a file holding an empty
        geometry.
    value_column : str
        Base-value column such as ``"k"``, ``"elev"``, or ``"cond"``.
        ``load_vector_parameter_source`` raises ``ValueError`` when the file has
        no such column.
    feature_id_column : str or None, optional
        Unique feature-identifier column. Default ``"name"``; when not
        ``None`` it must exist in the file (same check as ``value_column``).
        ``None`` skips that check.
    zone_column : str or None, optional
        Spatial zone column. Default ``None``. Stored only; nothing reads it.
    group_column : str or None, optional
        Logical group column. Default ``None``. Stored only; nothing reads it.
    layer_column : str or None, optional
        Layer column. Default ``None``. Stored only; nothing reads it.
    crs : int or None, optional
        EPSG code of the features. Default ``myflopy._defaults.DEFAULT_EPSG``
        (2927). Stored only: ``load_vector_parameter_source`` neither checks it
        nor reprojects.
    lower_bound_column : str or None, optional
        Column of per-feature lower bounds. Default ``None``. Stored only;
        nothing reads it.
    upper_bound_column : str or None, optional
        Column of per-feature upper bounds. Default ``None``. Stored only;
        nothing reads it.
    """

    path: str | Path
    value_column: str
    feature_id_column: str | None = "name"
    zone_column: str | None = None
    group_column: str | None = None
    layer_column: str | None = None
    crs: int | None = DEFAULT_EPSG
    lower_bound_column: str | None = None
    upper_bound_column: str | None = None


@dataclass
class ExpGeoStruct:
    """Simple exponential geostatistical structure definition.

    This wrapper stores values in a ``myflopy``-friendly form and can be
    converted to a real ``pyemu.geostats.GeoStruct`` lazily when pyEMU is
    available, by :func:`myflopy.modflow.mf6.pest.geostats.build_geostruct`.
    Nothing is validated here; pyEMU asserts ``range``, ``contribution`` and
    ``anisotropy`` are positive and ``transform`` is ``"none"`` or ``"log"``
    when that conversion runs (an ``AssertionError``).

    For each ``grid``/``pilotpoints`` parameter given a ``correlation=``, the
    project builds one of these from what ``cal.parameterize(...,
    correlation=, anisotropy=, bearing=, nugget=)`` recorded (at ``build()``
    and when drawing the prior), with ``contribution`` left at ``1.0`` and
    the parameter's transform; you only construct one yourself to call
    ``build_geostruct`` directly.

    Parameters
    ----------
    range : float
        Variogram range in model length units, passed to pyEMU as
        ``ExpVario(a=range)``. pyEMU's exponential covariance is
        ``contribution * exp(-h / range)`` at separation ``h``.
    contribution : float, optional
        Sill (variance contribution) of the exponential variogram, pyEMU's
        ``ExpVario(contribution=...)``. Default ``1.0``.
    anisotropy : float, optional
        Ratio of the correlation length along ``bearing`` to that across it,
        pyEMU's ``ExpVario(anisotropy=...)``. Default ``1.0`` (isotropic).
    bearing : float, optional
        Direction of the long axis, in degrees east of north, pyEMU's
        ``ExpVario(bearing=...)``. Default ``0.0``. No effect while
        ``anisotropy`` is ``1.0``.
    nugget : float, optional
        Nugget variance added to the covariance diagonal, pyEMU's
        ``GeoStruct(nugget=...)``. Default ``0.0``. The project's prior draws
        scale ``contribution`` and ``nugget`` to sum to ``1.0``
        (:func:`~myflopy.modflow.mf6.pest.geostats.unit_sill`) before handing the
        structure to pyEMU, whose draw replaces the diagonal with the
        bound-derived variance and so would discard a nugget on top of a sill-1
        variogram (ledger 489); :func:`build_geostruct` itself passes both
        through unchanged.
    transform : str, optional
        ``"none"`` (default) or ``"log"``: the transform of the property the
        structure describes, pyEMU's ``GeoStruct(transform=...)``.
    """

    range: float
    contribution: float = 1.0
    anisotropy: float = 1.0
    bearing: float = 0.0
    nugget: float = 0.0
    transform: str = "none"


@dataclass
class HeadTargetObservationSpec:
    """Pre-built head observations for the PEST build, wrapping :class:`HeadTargets`.

    Carries a :class:`~myflopy.modflow.mf6.observations.HeadTargets` plus the
    pyEMU observation-name ``prefix`` and, optionally, an explicit
    ``simulated_values`` table. ``cal.observe(...)`` accepts a bare ``HeadTargets``
    and wraps it for you, so you only construct this spec when you need to pin the
    prefix or supply simulated values directly; the project's native build turns
    it into pyEMU observations + the forward-run head post-processor.

    Parameters
    ----------
    targets : HeadTargets
        The head observation locations and target heads. Its rows are matched
        to model cells at build time; a row whose target matched no cell gets
        no PEST observation.
    simulated_values : str or pathlib.Path or pandas.DataFrame or None, optional
        The base-run simulated-head table the build writes as
        ``<prefix>_simulated_heads.csv`` and builds pyEMU's instruction file
        from. Default ``None`` has the forward run's own writer
        (``forward_run._write_head_target_csv``) write it from the model's
        current head file, so the base run and every forward run are one table
        in one order (ledger 340). A path is read with ``pandas.read_csv``; a
        DataFrame is copied -- and either must hold the forward run's columns in
        its order, ``per`` then each observation label in the order of the
        target map, or the build raises ``ValueError``: PEST reads the file by
        position. The forward run rewrites the file on every model run, so this
        only seeds the build.
    prefix : str, optional
        pyEMU observation-name prefix, default ``"hds"``. It also names the
        files the build writes into the template workspace
        (``<prefix>_simulated_heads.csv``, ``<prefix>_head_target_map.csv``,
        ``<prefix>_target_locations.gpkg``, ``<prefix>_target_values.csv``), so
        two sets need two prefixes.

    See Also
    --------
    myflopy.modflow.mf6.observations.HeadTargets
    """

    targets: HeadTargets
    simulated_values: str | Path | pd.DataFrame | None = None
    prefix: str = "hds"


@dataclass
class ConcObservationSpec:
    """Pre-built CONCENTRATION observations, wrapping :class:`ConcTargets`.

    The transport twin of :class:`HeadTargetObservationSpec`. ``cal.observe(...)``
    accepts a bare ``ConcTargets`` and wraps it for you.

    Parameters
    ----------
    targets : ConcTargets
        Concentration observation locations and measured concentrations,
        matched to model cells at build time.
    simulated_values : str or pathlib.Path or pandas.DataFrame or None, optional
        The base-run simulated-concentration table the build writes as
        ``<prefix>_simulated_conc.csv`` and builds pyEMU's instruction file
        from. Default ``None`` has the forward run's own writer
        (``forward_run._write_conc_target_csv``) write it from the transport
        model's ``.ucn`` in the model workspace (ledger 340). A path is read with
        ``pandas.read_csv``; a DataFrame is copied -- and either must hold the
        forward run's columns in its order (``per``, then each label in target-map
        order), or the build raises ``ValueError``. The forward run rewrites the
        file on every model run, so this only seeds the build.
    prefix : str, optional
        pyEMU observation-name prefix, default ``"conc"``. It also names the
        files the build writes (``<prefix>_simulated_conc.csv``,
        ``<prefix>_conc_target_map.csv``, ``<prefix>_target_locations.gpkg``,
        ``<prefix>_target_values.csv``), so two sets need two prefixes.
    transport_model_name : str or None, optional
        Name of the GWT model whose ``.ucn`` holds the simulated field. Left
        ``None`` it is discovered from the simulation -- the calibration hangs
        off the FLOW model (that is where ``k``/``recharge`` live), so the
        transport sibling has to be found rather than assumed. The discovery
        raises ``ValueError`` at build when the simulation holds no GWT model or
        more than one; name it here in the second case.

    See Also
    --------
    myflopy.modflow.mf6.observations.ConcTargets
    """

    targets: ConcTargets
    simulated_values: str | Path | pd.DataFrame | None = None
    prefix: str = "conc"
    transport_model_name: str | None = None


@dataclass
class LakeStageObservationSpec:
    """Pre-built lake-stage observations, wrapping :class:`LakeStageTargets`.

    ``cal.observe(...)`` / ``cal.forecast(...)`` accept a bare
    ``LakeStageTargets`` and wrap it; build this spec explicitly to pin a
    non-default ``prefix`` or a uniform ``weight``. Although ``targets`` has a
    default, the build REQUIRES it: a spec without one raises ``ValueError``
    when ``cal.build()`` prepares the observations.

    The first four fields predate ``targets=`` and are read by nothing: the
    build takes names, values and times from ``targets`` alone. They are kept
    only so old call sites still construct.

    Parameters
    ----------
    lake_names : list of str or None, optional
        Ignored (see above). Default ``None``.
    values : str or pathlib.Path or pandas.DataFrame or None, optional
        Ignored (see above). Default ``None``.
    time_column : str, optional
        Ignored (see above). Default ``"time"``.
    value_column : str, optional
        Ignored (see above). Default ``"stage"``.
    weight : float or None, optional
        The set's observation weight. When not ``None`` every target row holding
        a value and no weight of its own gets it (cast with ``float``): a
        reading's weight in the targets' values table wins (ledger
        504; this used to overwrite it). Default ``None``: those rows get
        weight ``1.0``. Every row no target measures gets weight 0 either way
        (ledger 340), and a forecast's weights are zero regardless.
    prefix : str, optional
        pyEMU observation-name prefix, default ``"lak_stage"``; it is the set's
        PEST observation group (ledger 358). It also names the files the build
        writes (``<prefix>_simulated_lake_stage.csv``,
        ``<prefix>_target_values.csv``, ``<prefix>_target_locations.csv``), so
        two sets need two prefixes. The default was ``"stage"`` until ledger
        372: PESTPP-IES matches an ``ies_phi_factor_file`` tag anywhere in
        a group's name, and ``stage`` is inside an SFR set's ``sfr_stage``, so
        no tag could give a lake set's group a phi factor of its own.
    targets : LakeStageTargets or None, optional
        The named lake-stage targets. Required in practice: the default
        ``None`` constructs, and ``cal.build()`` then raises ``ValueError``.
    """

    lake_names: list[str] | None = None
    values: str | Path | pd.DataFrame | None = None
    time_column: str = "time"
    value_column: str = "stage"
    weight: float | None = None
    prefix: str = "lak_stage"
    targets: LakeStageTargets | None = None


@dataclass
class SfrStageObservationSpec:
    """Pre-built SFR-stage observations, wrapping :class:`SfrStageTargets`.

    Pairs a :class:`~myflopy.modflow.mf6.observations.SfrStageTargets` with the
    pyEMU observation-name ``prefix``. ``cal.observe(...)`` accepts a bare
    ``SfrStageTargets`` and wraps it for you; build this spec explicitly to pin a
    non-default prefix or a set ``weight``. The native build registers the observations and the
    forward run regenerates the simulated SFR-stage CSV.

    Parameters
    ----------
    targets : SfrStageTargets
        The named-reach stage targets. Their simulated series is read from the
        model's output at build time and regenerated by the forward run.
    prefix : str, optional
        pyEMU observation-name prefix, default ``"sfr_stage"``. It also names
        the files the build writes (``<prefix>_simulated_sfr_stage.csv``,
        ``<prefix>_target_values.csv``, ``<prefix>_target_locations.csv``), so
        two sets need two prefixes.
    weight : float or None, optional
        The set's observation weight (ledger 504, as a lake set's): every
        target row holding a value and no weight of its own in the targets'
        values table gets it (cast with ``float``). Default ``None``: those rows
        get weight ``1.0``. Rows no target measures get weight 0 either way, and
        a forecast's weights are zero regardless. ``cal.observe(targets)`` builds
        the spec without one; build this spec to set it. PESTPP-IES draws a
        reading's noise with standard deviation 1/weight, so a flow in
        ft3/d at weight 1 carries 1 ft3/d of noise.
    """

    targets: SfrStageTargets
    prefix: str = "sfr_stage"
    weight: float | None = None


@dataclass
class SfrFlowObservationSpec:
    """Pre-built SFR-flow observations, wrapping :class:`SfrFlowTargets`.

    Pairs a :class:`~myflopy.modflow.mf6.observations.SfrFlowTargets` with the
    pyEMU observation-name ``prefix``. ``cal.observe(...)`` / ``cal.forecast(...)``
    accept a bare ``SfrFlowTargets`` and wrap it; build this spec explicitly to
    pin a non-default prefix or a set ``weight``. Typical for gaged-discharge calibration targets.

    Parameters
    ----------
    targets : SfrFlowTargets
        The named-reach flow targets. Their simulated series is read from the
        model's output at build time and regenerated by the forward run.
    prefix : str, optional
        pyEMU observation-name prefix, default ``"sfr_flow"``. It also names
        the files the build writes (``<prefix>_simulated_sfr_flow.csv``,
        ``<prefix>_target_values.csv``, ``<prefix>_target_locations.csv``), so
        two sets need two prefixes.
    weight : float or None, optional
        The set's observation weight (ledger 504, as a lake set's): every
        target row holding a value and no weight of its own in the targets'
        values table gets it (cast with ``float``). Default ``None``: those rows
        get weight ``1.0``. Rows no target measures get weight 0 either way, and
        a forecast's weights are zero regardless. ``cal.observe(targets)`` builds
        the spec without one; build this spec to set it. PESTPP-IES draws a
        reading's noise with standard deviation 1/weight, so a flow in
        ft3/d at weight 1 carries 1 ft3/d of noise.
    """

    targets: SfrFlowTargets
    prefix: str = "sfr_flow"
    weight: float | None = None


@dataclass
class DrnFlowObservationSpec:
    """Pre-built DRN-seepage observations, wrapping :class:`DrnFlowTargets`.

    Pairs a :class:`~myflopy.modflow.mf6.observations.DrnFlowTargets` with the
    pyEMU observation-name ``prefix``. ``cal.observe(...)`` / ``cal.forecast(...)``
    accept a bare ``DrnFlowTargets`` and wrap it; build this spec explicitly only
    to pin a non-default prefix. Used for spring/seep discharge aggregated over
    a group of drain cells.

    Parameters
    ----------
    targets : DrnFlowTargets
        The seepage-zone flow targets, each aggregating a group of drain
        cells. Their simulated series is read from the model's output at build
        time and regenerated by the forward run.
    prefix : str, optional
        pyEMU observation-name prefix, default ``"drn_flow"``. It also names
        the files the build writes (``<prefix>_simulated_drn_flow.csv``,
        ``<prefix>_target_values.csv``, ``<prefix>_target_locations.*``), so
        two sets need two prefixes.
    """

    targets: DrnFlowTargets
    prefix: str = "drn_flow"


@dataclass
class RivFlowObservationSpec:
    """Pre-built RIV exchange-flow (baseflow) observations, wrapping :class:`RivFlowTargets`.

    Pairs a :class:`~myflopy.modflow.mf6.observations.RivFlowTargets` with the
    pyEMU observation-name ``prefix``. ``cal.observe(...)`` / ``cal.forecast(...)``
    accept a bare ``RivFlowTargets`` and wrap it; build this spec explicitly only
    to pin a non-default prefix. Used for a creek's gauged baseflow, summed over
    the river cells of a reach.

    Parameters
    ----------
    targets : RivFlowTargets
        The river-zone exchange-flow targets, each aggregating a group of river
        cells, in MF6's sign (a gaining reach's baseflow negative). Their
        simulated series is read from the model's output at build time and
        regenerated by the forward run.
    prefix : str, optional
        pyEMU observation-name prefix, default ``"riv_flow"``. It also names
        the files the build writes (``<prefix>_simulated_riv_flow.csv``,
        ``<prefix>_target_values.csv``, ``<prefix>_target_locations.*``), so
        two sets need two prefixes; and it is the PEST observation group of the
        zones given no ``group`` label (ledger 358).
    """

    targets: RivFlowTargets
    prefix: str = "riv_flow"


@dataclass
class HeadDifferenceObservationSpec:
    """Pre-built head-DIFFERENCE observations, wrapping :class:`HeadDifferences` (ledger 565).

    ``cal.observe(head_targets.differences())`` wraps the set for you; build this spec
    to pin a prefix or supply simulated values. Every measured difference is in ONE
    PEST group, the prefix; every other row in ``<prefix>_unmeasured`` at weight 0.

    Parameters
    ----------
    targets : HeadDifferences
        The differenced set (``HeadTargets.differences(...)``).
    simulated_values : str or pathlib.Path or pandas.DataFrame or None, optional
        The base-run table of simulated DIFFERENCES the build writes as
        ``<prefix>_simulated_differences.csv``. Default ``None`` has the forward run's
        own writer (``forward_run._write_target_difference_csv``) write it from the
        model's head file, so the base run and every forward run are one table in
        one order (ledger 340). A path is read with ``pandas.read_csv``; a DataFrame
        is copied -- either must hold ``per`` then each kept label in the target-map
        order, or the build raises ``ValueError``.
    prefix : str, optional
        pyEMU observation-name prefix and the set's PEST group, default ``"hdiff"``.
        It names the files the build writes (``<prefix>_simulated_differences.csv``,
        ``<prefix>_difference_map.csv``, ``<prefix>_difference_periods.csv``,
        ``<prefix>_target_locations.gpkg``, ``<prefix>_target_values.csv``), so two
        sets need two prefixes.

    See Also
    --------
    myflopy.modflow.mf6.observations.HeadDifferences
    """

    targets: HeadDifferences
    simulated_values: str | Path | pd.DataFrame | None = None
    prefix: str = "hdiff"


@dataclass
class ConcDifferenceObservationSpec:
    """Pre-built concentration-DIFFERENCE observations, wrapping :class:`ConcDifferences` (ledger 565).

    The transport twin of :class:`HeadDifferenceObservationSpec`;
    ``cal.observe(conc_targets.differences())`` wraps the set for you.

    Parameters
    ----------
    targets : ConcDifferences
        The differenced concentration set (``ConcTargets.differences(...)``).
    simulated_values : str or pathlib.Path or pandas.DataFrame or None, optional
        The base-run table of simulated differences, as
        :class:`HeadDifferenceObservationSpec` reads it, from the transport model's
        concentration file. Default ``None``: the forward run's writer writes it.
    prefix : str, optional
        pyEMU observation-name prefix and the set's PEST group, default ``"cdiff"``;
        it names the set's files, as :class:`HeadDifferenceObservationSpec`'s does.
    transport_model_name : str or None, optional
        The GWT model whose concentration file is read; ``None`` (default) finds the
        simulation's one GWT model, raising ``ValueError`` at build when it holds none
        or several.

    See Also
    --------
    myflopy.modflow.mf6.observations.ConcDifferences
    """

    targets: ConcDifferences
    simulated_values: str | Path | pd.DataFrame | None = None
    prefix: str = "cdiff"
    transport_model_name: str | None = None
