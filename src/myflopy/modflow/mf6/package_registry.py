"""Registry metadata for package input and result explorers."""

from __future__ import annotations

from dataclasses import dataclass, field

# Suffix (kind.ext) for each output filerecord an advanced package (uzf/lak/sfr)
# emits. The full record is always ``{model_name}_{name}_{suffix}`` -- see
# advanced_output_filerecords.
_ADVANCED_OUTPUT_FILERECORD_SUFFIXES: dict[str, dict[str, str]] = {
    "uzf": {
        "budget_filerecord": "budget.uzf",
        "budgetcsv_filerecord": "budget.csv",
        "package_convergence_filerecord": "package_convergence.csv",
    },
    "lak": {
        "stage_filerecord": "stage.lak",
        "budget_filerecord": "budget.lak",
        "budgetcsv_filerecord": "budget.csv",
        "package_convergence_filerecord": "package_convergence.csv",
    },
    "sfr": {
        "stage_filerecord": "stage.sfr",
        "budget_filerecord": "budget.sfr",
    },
}


def advanced_output_filerecords(
    package_type: str, name: str, model_name: str = "{model_name}"
) -> dict[str, str]:
    """The MF6 output filerecords for one advanced package (uzf/lak/sfr).

    Every record is uniquely named ``{model_name}_{name}_{kind}.{ext}``. This is
    the ONE source of truth shared by the ``*_spec`` factories (build), the
    legacy imperative ``.uzf()`` path, and
    :func:`myflopy.project.components.apply_package_artifact` (restore), so that
    (a) two packages of the same type never collide on an output file -- e.g.
    ``uzf`` budget used to ignore ``name`` and two UZF packages both wrote
    ``{model}_budget.uzf`` -- and (b) a restored package's output files match the
    ones the original build wrote.

    Lives here (layer 0, the single source of per-package truth) so every build
    and restore path can import it downward at module level. ``model_name``
    defaults to the literal ``"{model_name}"`` template the spec factories
    resolve later via ``str.format``; the restore path passes a concrete name.
    """

    suffixes = _ADVANCED_OUTPUT_FILERECORD_SUFFIXES[package_type]
    return {key: f"{model_name}_{name}_{suffix}" for key, suffix in suffixes.items()}


#: The kinds a :class:`FieldSpec` may declare. A closed set, checked at construction,
#: because a typo in a kind would silently fall through to the intensive default and
#: quietly change how a field aggregates.
_FIELD_KINDS = frozenset({"extensive", "elevation", "intensive"})


def agg_for_kind(kind: str) -> str:
    """How several values of one quantity combine, from what the quantity IS.

    THE rule of ledger 297, in one place: only an ``"extensive"`` value -- a flux or
    a conductance, which MF6 superposes -- ADDS; an elevation or an intensive value
    never does. :attr:`FieldSpec.agg` returns this for its own ``kind``, and the
    observation targets' reduction (``observations._shared._aggregate_named_series``)
    asks it before it decides whether rows add, so a stage is never summed there
    either (ledger 331: a pond's stage read 763 times itself).

    Parameters
    ----------
    kind : str
        What the number is: ``"extensive"``, ``"elevation"`` or ``"intensive"``,
        the closed set a :class:`FieldSpec` declares.

    Returns
    -------
    str
        ``"sum"`` for ``"extensive"``, ``"mean"`` for the other two.

    Raises
    ------
    ValueError
        If ``kind`` is not one of the three: a typo would otherwise fall through
        to the non-additive branch and quietly change how a value combines.

    Examples
    --------
    >>> agg_for_kind("extensive")
    'sum'
    >>> agg_for_kind("elevation")
    'mean'
    """

    if kind not in _FIELD_KINDS:
        raise ValueError(f"kind {kind!r} is not one of {sorted(_FIELD_KINDS)}.")
    return "sum" if kind == "extensive" else "mean"


@dataclass(frozen=True)
class FieldSpec:
    """Registry metadata for one package input field."""

    name: str
    label: str | None = None
    colorscale: str | None = None
    fill_value: float = 0.0
    #: What the number IS, which decides how several records on one cell combine and
    #: how the field would survive a regrid. ONE declaration, three consumers: `agg`
    #: below, and `geopackage`'s elevation and extensive sets (ledger 297).
    #:
    #: ``"extensive"``  a flux or a conductance -- MF6 superposes them, so they ADD;
    #: ``"elevation"``  a water level or an invert, which `mf.Spread` is refused on;
    #: ``"intensive"``  a length, a fraction or a coefficient -- a level, not an amount.
    kind: str = "intensive"

    def __post_init__(self) -> None:
        """Refuse a kind nothing downstream knows how to treat."""

        if self.kind not in _FIELD_KINDS:
            raise ValueError(
                f"FieldSpec {self.name!r}: kind {self.kind!r} is not one of "
                f"{sorted(_FIELD_KINDS)}."
            )

    @property
    def agg(self) -> str:
        """How several records on ONE cell combine for a map or table.

        DERIVED from :attr:`kind` rather than declared beside it, because two fields
        that must agree are two fields that drift. Every entry read ``"sum"`` before,
        elevations included -- so a cell holding two CHD records of 100.0 and 102.0
        rendered **202.0** on `packages.chd.inputs.head.map()` (ledger 297). Summing a
        flux is right; summing a water level is not.

        A caller who wants something else passes ``agg=`` at the call, which every map
        and table already accepts -- this is the DEFAULT, not a restriction. The rule
        itself is :func:`agg_for_kind`, which the observation targets ask too.
        """

        return agg_for_kind(self.kind)


@dataclass(frozen=True)
class ResultSpec:
    """Registry metadata for one package result term.

    ``reference_frame`` declares whose point of view MF6 wrote this record from.
    myflopy keeps MF6's RAW sign everywhere and never negates ``q``; the
    ambiguity that a gaining stream and a gaining lake report opposite signs is
    resolved by naming the frame -- both in this field (so plots orient their
    colours) and in the emitted column (``q_gwf`` vs ``q_lake``).

    * ``"gwf"`` -- an aquifer-referenced cell-by-cell (.cbc) record, written as
      "flow FROM the feature TO the GWF cell". A GAINING feature is therefore
      NEGATIVE (water leaving the aquifer). This is SFR's ``SFR`` record and
      every list BC (CHD/DRN/GHB/RIV/WEL/RCH/EVT). Column: ``q_gwf``.
    * ``"feature"`` -- a feature-referenced package budget file, written from the
      FEATURE's perspective, so a positive value means water ENTERING the
      feature (gaining). This is LAK's ``GWF`` record -- named ``GWF`` in the
      file but written lake's-perspective, which is the trap. Column: ``q_lake``.

    Verified against the canonical model: its perched lake leaks downward and
    reports ``GWF`` in the range -5888..0 (feature frame: negative = losing),
    while its losing stream reaches report positive ``q`` (gwf frame: positive =
    losing). Same physical direction, opposite signs -- exactly why the frame
    must be declared rather than assumed.

    Unlike the ``raw_gaining_sign`` field this replaced, ``reference_frame``
    cannot go stale: it describes what MF6 *wrote*, and the data is never
    flipped out from under it. The old field named both the raw sign and the
    effective sign, which are equal only before a normalization step -- so it
    silently inverted the SFR map on 2026-07-19 when consulted after that step.
    There is no normalization step any more.
    """

    name: str
    budget_text: str
    value_name: str = "q"
    label: str | None = None
    colorscale: str | None = None
    diverging: bool = True
    reference_frame: str = "gwf"


@dataclass(frozen=True)
class PackageCapabilities:
    """What MF6/FloPy will accept for one package.

    Mirrors the live FloPy constructor -- ``tests/test_package_descriptor.py``
    asserts each flag against ``inspect.signature`` rather than against a
    literal, so a FloPy upgrade that changes a package's options fails loudly
    instead of leaving the registry quietly wrong.

    ``edges_only`` is the odd one out: it is a *myflopy* capability, not an MF6
    one -- whether ``GeoPackageSource.<pkg>`` can restrict mapped features to
    perimeter cells. Only the head-dependent boundaries expose it.
    """

    mover: bool = False
    auxiliary: bool = True
    boundnames: bool = True
    observations: bool = True
    edges_only: bool = False


@dataclass(frozen=True)
class PackageTiers:
    """Which cross-cutting subsystems cover one package.

    These are deliberately separate booleans rather than one "supported" flag:
    the subsystems genuinely disagree about which packages they handle, and
    several of those disagreements are unintended gaps (plan 4.7). Recording
    them per tier is what makes the gaps visible instead of implied.
    """

    #: input-diffable row-by-row (``model_diff._DIFF_PACKAGES``)
    diffable: bool = False
    #: input-diffed by connection/reach geometry (``_CONNECTION_PACKAGES``)
    connection_diffable: bool = False
    #: results-diffable per cell (``model_results_diff._CELL_BUDGET_PACKAGES``)
    results_diffable: bool = False
    #: capturable/restorable as a run artifact (``project.components``)
    artifact_serializable: bool = False
    #: relative order artifacts are re-applied in; None when not serializable
    artifact_apply_order: int | None = None
    #: has a direct ``SimulationBase.<pkg>`` accessor (several do NOT -- 4.7)
    model_accessor: bool = False


@dataclass(frozen=True)
class PackageExplorerSpec:
    """Registry metadata for one MF6 package explorer.

    Plan 4.7.2 grew this from explorer-only metadata into the single source of
    per-package truth. The fields below the original five are the knowledge
    that was previously restated at ~92 sites across ``geopackage.py``,
    ``package_api.py``, ``advanced.py``, ``run_model.py``, ``budget_tables.py``,
    ``model_diff.py`` and ``components.py``. Nothing consumes them yet -- 4.7.3
    deletes the hardcoded lists one at a time -- but every value is pinned by a
    test that reads the ORIGINAL source, so the descriptor cannot drift away
    from the code it is about to replace.
    """

    name: str
    kind: str = "cell_stress"
    default_input: str | None = None
    colorscale: str | None = None
    inputs: dict[str, FieldSpec] = field(default_factory=dict)
    results: dict[str, ResultSpec] = field(default_factory=dict)

    # -- 4.7.2: the knowledge the other sites currently restate ---------------
    #: FloPy class NAME, not the class. Kept as a string so this module stays
    #: import-free and the deferred-import ratchet is undisturbed; callers do
    #: ``getattr(flopy.mf6, spec.flopy_class)``.
    flopy_class: str | None = None
    #: stress-period record fields AFTER the cellid, in MF6 order. Empty for
    #: the advanced packages, whose input is packagedata, not a flat record.
    record_fields: tuple[str, ...] = ()
    #: ``GeoPackageSource.<pkg>`` parameter -> the column name it defaults to.
    #: Empty when the package has no GeoPackage resolver.
    gpkg_defaults: dict[str, str] = field(default_factory=dict)
    capabilities: PackageCapabilities = field(default_factory=PackageCapabilities)
    tiers: PackageTiers = field(default_factory=PackageTiers)
    #: the MF6 input-file suffix (``model.<suffix>``)
    file_suffix: str | None = None
    #: Whether this package's budget ``node2`` is zero-based by myflopy. MISNAMED
    #: as of 2026-07-27 -- it no longer governs ``node`` at all (rename deferred,
    #: ledger 93).
    #:
    #: It used to, on the belief that "the advanced packages report feature
    #: numbers instead" of cells. **That belief was measured FALSE.** In the
    #: *model* budget file EVERY record's ``node`` is a 1-based model cell,
    #: ``SFR``/``LAK``/``UZF-GWRCH`` included -- there it is ``node2`` that holds
    #: the feature id. (The feature-first layout the belief described is real but
    #: belongs to the separate *package-output* budget file; the two were
    #: conflated.) So ``node`` is now zero-based once, universally, in
    #: ``package_explorer_utils._model_budget_record_frame``, and this flag is
    #: left governing only the package-specific ``node2``.
    #:
    #: Got the ``node`` half wrong for years -- see the budget off-by-one fixed
    #: in 4.7.1, and its SFR/LAK/UZF tail in ledger 92.
    zero_base_budget_nodes: bool = False
    #: one hand-written sentence naming what the package IS, hydrologically.
    #: Not generated: prose stays human-written (see the deprecation policy's
    #: companion rule for docstrings).
    blurb: str = ""


_PACKAGE_EXPLORER_SPECS: dict[str, PackageExplorerSpec] = {
    # Colorscale policy: diverging red/white/blue is reserved for signed,
    # gaining/losing "q"-like fields and diff maps; everything else uses the
    # house brown-to-blue "earth" scale (the mounding-figure default).
    "rch": PackageExplorerSpec(
        name="rch",
        flopy_class="ModflowGwfrch",
        record_fields=("recharge",),
        gpkg_defaults={"recharge": "recharge"},
        capabilities=PackageCapabilities(mover=False, edges_only=False),
        tiers=PackageTiers(
            diffable=True,
            results_diffable=True,
            artifact_serializable=True,
            artifact_apply_order=30,
            model_accessor=True,
        ),
        file_suffix="rch",
        zero_base_budget_nodes=True,
        blurb=("Areally distributed recharge applied to the top active cell."),
        default_input="recharge",
        colorscale="earth",
        inputs={
            "recharge": FieldSpec("recharge", label="Recharge", colorscale="earth", kind="extensive"),
        },
        results={
            "q": ResultSpec("q", budget_text="RCH", value_name="q_gwf", colorscale="RdBu"),
        },
    ),
    "chd": PackageExplorerSpec(
        name="chd",
        flopy_class="ModflowGwfchd",
        record_fields=("head",),
        gpkg_defaults={"head": "head"},
        capabilities=PackageCapabilities(mover=False, edges_only=True),
        tiers=PackageTiers(
            diffable=True,
            results_diffable=True,
            artifact_serializable=True,
            artifact_apply_order=40,
            model_accessor=True,
        ),
        file_suffix="chd",
        zero_base_budget_nodes=True,
        blurb=(
            "Constant-head cells that hold a prescribed head, sourcing or sinking whatever flow that requires."
        ),
        default_input="head",
        colorscale="earth",
        inputs={
            "head": FieldSpec("head", label="Constant head", colorscale="earth", kind="elevation"),
        },
        results={
            "q": ResultSpec("q", budget_text="CHD", value_name="q_gwf", colorscale="RdBu"),
        },
    ),
    "drn": PackageExplorerSpec(
        name="drn",
        flopy_class="ModflowGwfdrn",
        record_fields=("elev", "cond"),
        gpkg_defaults={"elevation": "elevation", "conductance": "conductance"},
        capabilities=PackageCapabilities(mover=True, edges_only=True),
        tiers=PackageTiers(
            diffable=True,
            results_diffable=True,
            artifact_serializable=True,
            artifact_apply_order=50,
            model_accessor=True,
        ),
        file_suffix="drn",
        zero_base_budget_nodes=True,
        blurb=(
            "Drains that remove water only while head stands above the drain elevation, and never add any."
        ),
        default_input="elev",
        colorscale="earth",
        inputs={
            "elev": FieldSpec("elev", label="Drain elevation", colorscale="earth", kind="elevation"),
            "cond": FieldSpec("cond", label="Drain conductance", colorscale="earth", kind="extensive"),
        },
        results={
            "q": ResultSpec("q", budget_text="DRN", value_name="q_gwf", colorscale="RdBu"),
        },
    ),
    "ghb": PackageExplorerSpec(
        name="ghb",
        flopy_class="ModflowGwfghb",
        record_fields=("bhead", "cond"),
        gpkg_defaults={"head": "head", "conductance": "conductance"},
        capabilities=PackageCapabilities(mover=True, edges_only=True),
        tiers=PackageTiers(
            diffable=True,
            results_diffable=True,
            artifact_serializable=True,
            artifact_apply_order=60,
            model_accessor=True,
        ),
        file_suffix="ghb",
        zero_base_budget_nodes=True,
        blurb=(
            "General-head boundary: flow proportional to the difference between a distant boundary head and the cell head, through a fixed conductance."
        ),
        default_input="bhead",
        colorscale="earth",
        inputs={
            "bhead": FieldSpec("bhead", label="Boundary head", colorscale="earth", kind="elevation"),
            "cond": FieldSpec("cond", label="Boundary conductance", colorscale="earth", kind="extensive"),
        },
        results={
            "q": ResultSpec("q", budget_text="GHB", value_name="q_gwf", colorscale="RdBu"),
        },
    ),
    "riv": PackageExplorerSpec(
        name="riv",
        flopy_class="ModflowGwfriv",
        record_fields=("stage", "cond", "rbot"),
        gpkg_defaults={"stage": "stage", "conductance": "conductance", "rbot": "rbot"},
        capabilities=PackageCapabilities(mover=True, edges_only=True),
        tiers=PackageTiers(
            diffable=True,
            results_diffable=True,
            artifact_serializable=True,
            artifact_apply_order=65,
            model_accessor=True,
        ),
        file_suffix="riv",
        zero_base_budget_nodes=True,
        blurb=(
            "A river reach exchanging flow with the aquifer through a streambed conductance, limited by the bed bottom once the aquifer falls below it."
        ),
        default_input="stage",
        colorscale="earth",
        inputs={
            "stage": FieldSpec("stage", label="River stage", colorscale="earth", kind="elevation"),
            "cond": FieldSpec("cond", label="Riverbed conductance", colorscale="earth", kind="extensive"),
            "rbot": FieldSpec("rbot", label="River bottom", colorscale="earth", kind="elevation"),
        },
        results={
            "q": ResultSpec("q", budget_text="RIV", value_name="q_gwf", colorscale="RdBu"),
        },
    ),
    "wel": PackageExplorerSpec(
        name="wel",
        flopy_class="ModflowGwfwel",
        record_fields=("q",),
        gpkg_defaults={"rate": "rate"},
        capabilities=PackageCapabilities(mover=True, edges_only=False),
        tiers=PackageTiers(
            diffable=True,
            results_diffable=True,
            artifact_serializable=True,
            artifact_apply_order=45,
            model_accessor=True,
        ),
        file_suffix="wel",
        zero_base_budget_nodes=True,
        blurb=("Wells injecting or withdrawing a specified volumetric rate, independent of head."),
        default_input="q",
        colorscale="RdBu",
        inputs={
            "q": FieldSpec("q", label="Well flow", colorscale="RdBu", kind="extensive"),
        },
        results={
            "q": ResultSpec("q", budget_text="WEL", value_name="q_gwf", colorscale="RdBu"),
        },
    ),
    "evt": PackageExplorerSpec(
        name="evt",
        flopy_class="ModflowGwfevt",
        record_fields=("surface", "rate", "depth"),
        gpkg_defaults={"surface": "surface", "rate": "rate", "depth": "depth"},
        capabilities=PackageCapabilities(mover=False, edges_only=False),
        tiers=PackageTiers(
            diffable=True,
            results_diffable=True,
            artifact_serializable=True,
            artifact_apply_order=35,
            model_accessor=True,
        ),
        file_suffix="evt",
        zero_base_budget_nodes=True,
        blurb=(
            "Evapotranspiration drawn from the water table at a rate that falls from a maximum at the ET surface to zero at the extinction depth."
        ),
        default_input="rate",
        colorscale="earth",
        inputs={
            "surface": FieldSpec("surface", label="ET surface", colorscale="earth", kind="elevation"),
            "rate": FieldSpec("rate", label="Maximum ET rate", colorscale="earth", kind="extensive"),
            "depth": FieldSpec("depth", label="Extinction depth", colorscale="earth"),
        },
        results={
            "q": ResultSpec("q", budget_text="EVT", value_name="q_gwf", colorscale="RdBu"),
        },
    ),
    "uzf": PackageExplorerSpec(
        name="uzf",
        flopy_class="ModflowGwfuzf",
        record_fields=(),
        gpkg_defaults={},
        capabilities=PackageCapabilities(mover=True, edges_only=False),
        tiers=PackageTiers(
            diffable=False,
            results_diffable=False,
            artifact_serializable=True,
            artifact_apply_order=70,
            model_accessor=True,
        ),
        file_suffix="uzf",
        zero_base_budget_nodes=False,
        blurb=(
            "Unsaturated-zone flow routing infiltration through a vadose column before it reaches the water table, with optional vadose ET."
        ),
        kind="uzf",
        inputs={
            "finf": FieldSpec("finf", label="UZF infiltration", colorscale="earth", kind="extensive"),
            "pet": FieldSpec("pet", label="Potential evapotranspiration", colorscale="earth", kind="extensive"),
            "extdp": FieldSpec("extdp", label="ET extinction depth", colorscale="earth"),
            "extwc": FieldSpec("extwc", label="ET extinction water content", colorscale="earth"),
            "ha": FieldSpec("ha", label="Surface depression storage depth", colorscale="earth"),
            "hroot": FieldSpec("hroot", label="Root zone thickness", colorscale="earth"),
            "rootact": FieldSpec("rootact", label="Root activity", colorscale="earth"),
        },
        results={
            "gwrch": ResultSpec(
                "gwrch",
                budget_text="UZF-GWRCH",
                value_name="gwrch",
                label="UZF groundwater recharge",
                colorscale="earth",
                diverging=False,
            ),
            "sat": ResultSpec(
                "sat",
                budget_text="DATA-SAT",
                value_name="sat",
                label="UZF saturation",
                colorscale="earth",
                diverging=False,
            ),
        },
    ),
    "sfr": PackageExplorerSpec(
        name="sfr",
        flopy_class="ModflowGwfsfr",
        record_fields=(),
        gpkg_defaults={},
        capabilities=PackageCapabilities(mover=True, edges_only=False),
        tiers=PackageTiers(
            diffable=False,
            connection_diffable=True,
            results_diffable=True,
            artifact_serializable=True,
            artifact_apply_order=90,
            model_accessor=True,
        ),
        file_suffix="sfr",
        zero_base_budget_nodes=False,
        blurb=(
            "Streamflow routing through connected reaches, exchanging with the aquifer along the way."
        ),
        kind="surface_water",
        results={
            "q": ResultSpec("q", budget_text="SFR", value_name="q_gwf", colorscale="RdBu"),
        },
    ),
    "lak": PackageExplorerSpec(
        name="lak",
        flopy_class="ModflowGwflak",
        record_fields=(),
        gpkg_defaults={},
        capabilities=PackageCapabilities(mover=True, edges_only=False),
        tiers=PackageTiers(
            diffable=False,
            connection_diffable=True,
            results_diffable=True,
            artifact_serializable=True,
            artifact_apply_order=80,
            model_accessor=True,
        ),
        file_suffix="lak",
        zero_base_budget_nodes=False,
        blurb=(
            "A lake whose stage responds to its own water balance while it exchanges with the aquifer through lakebed connections."
        ),
        kind="surface_water",
        results={
            # LAK's GWF record comes from the LAK package budget file, written
            # from the LAKE's perspective: positive = water entering the lake.
            # Opposite frame to SFR's cell record -- see ResultSpec.reference_frame.
            "q": ResultSpec(
                "q",
                budget_text="GWF",
                value_name="q_lake",
                colorscale="RdBu",
                reference_frame="feature",
            ),
        },
    ),
}

# The policy for grouped difference maps: diverging, negative red / positive
# blue. Explicit STOPS, not the name "RdBu" this used to be: the stops carry
# the ORIENTATION (which end is red), which a bare name cannot. (A retired hand
# table once also sent the name to matplotlib's REVERSED colormap -- ledger
# 69/70, closed by 253 -- so this was doubly necessary then.)
#
# Stated here rather than imported from `package_plotting`: this registry is the
# low-level source of per-package truth and must not depend on a drawing module
# (that import inverts the layering and pushes half the package up a layer).
# `test_the_two_diverging_orientations_do_not_drift` pins these against
# `_red_white_blue_diverging_colorscale()` so the duplication cannot rot.
_GROUP_COMPARE_DEFAULT_COLORSCALE = [
    [0.0, "#d62728"],
    [0.5, "#ffffff"],
    [1.0, "#1f77b4"],
]


def get_default_package_value_column(package_name: str) -> str | None:
    """Return the preferred primary numeric input column for one package."""

    spec = get_package_explorer_spec(package_name)
    return spec.default_input if spec is not None else None


def get_default_package_colorscale(package_name: str) -> str | None:
    """Return the preferred choropleth colorscale for one package or field."""

    normalized = str(package_name).lower()
    if "_" in normalized:
        package, field_name = normalized.split("_", 1)
        field_spec = get_package_input_field_spec(package, field_name)
        if field_spec is not None and field_spec.colorscale is not None:
            return field_spec.colorscale

    spec = get_package_explorer_spec(normalized)
    return spec.colorscale if spec is not None else None


def get_default_budget_term(package_name: str) -> tuple[str, str] | None:
    """Return the preferred budget text and public value name for a package."""

    normalized = str(package_name).lower()
    if "_" in normalized:
        package, result_name = normalized.split("_", 1)
        result_spec = get_package_result_spec(package, result_name)
    else:
        result_spec = get_package_result_spec(normalized, "q")
    if result_spec is None:
        return None
    return result_spec.budget_text, result_spec.value_name


def get_package_explorer_spec(package_name: str) -> PackageExplorerSpec | None:
    """Return registry metadata for one package, if it is known."""

    return _PACKAGE_EXPLORER_SPECS.get(str(package_name).lower())


def get_package_input_field_spec(package_name: str, field_name: str) -> FieldSpec | None:
    """Return registry metadata for one input field, if it is known."""

    spec = get_package_explorer_spec(package_name)
    if spec is None:
        return None
    return spec.inputs.get(str(field_name).lower())


#: The LAK/SFR input fields that are TOTALS -- an area, a length, a volumetric rate --
#: and so ADD over the records they combine; every other one (an elevation, a width, a
#: leakance, a roughness, a starting stage) averages. The registry declares no LAK/SFR
#: input field, so this set is their ``FieldSpec.kind`` (ledger 475).
_SURFACE_WATER_EXTENSIVE_INPUTS = frozenset({"connection_area", "rlen", "inflow", "runoff"})

#: The packages whose input fields :data:`_SURFACE_WATER_EXTENSIVE_INPUTS` classifies.
_SURFACE_WATER_INPUT_PACKAGES = frozenset({"lak", "sfr"})


def input_field_agg(package_name: str, field_name: str, agg: str | None = None) -> str:
    """How several records of one package input field combine, from what it IS.

    The default ``agg`` of every view that reduces an input field's records -- a
    map's cell, a series' period point, a group or diff comparison -- through
    :attr:`FieldSpec.agg` and so :func:`agg_for_kind` (ledger 297): ``"sum"`` for
    an extensive field (a conductance, a rate), ``"mean"`` for an elevation or an
    intensive one. The series panels and the group tier used to take ``"sum"``
    for every field, so ``model.drn.inputs.elev.plot()`` drew the SUM of the
    drains' elevations (ledger 472). A LAK or SFR input field, which the
    registry does not declare, takes its kind from one set instead: ``"sum"``
    for ``connection_area``, ``rlen``, ``inflow`` and ``runoff``, ``"mean"`` for
    every other column (ledger 475).

    Parameters
    ----------
    package_name : str
        The package type, any case (``"drn"``, ``"uzf"``, ``"lak"``).
    field_name : str
        The input field, any case, as the registry names it (``"elev"``,
        ``"cond"``, ``"finf"``), or a column of the LAK/SFR inputs table
        (``"belev"``, ``"rlen"``).
    agg : str or None, optional
        The caller's own choice, any pandas aggregation name; returned unchanged
        when given. Default ``None``: the field's kind decides.

    Returns
    -------
    str
        ``agg`` when given, else ``"sum"`` or ``"mean"`` by the field's kind.

    Raises
    ------
    ValueError
        When ``agg`` is ``None``, the package is not LAK or SFR, and the registry
        declares no such field (an auxiliary column, say): what its records ARE
        is unknown, so the caller must say how they combine rather than have one
        guessed.
    """

    if agg is not None:
        return str(agg)
    if str(package_name).lower() in _SURFACE_WATER_INPUT_PACKAGES:
        # The registry declares no LAK/SFR input field: their table's columns are
        # read off the package, so the kind is this set's (ledger 475).
        extensive = str(field_name).lower() in _SURFACE_WATER_EXTENSIVE_INPUTS
        return agg_for_kind("extensive" if extensive else "intensive")
    spec = get_package_input_field_spec(package_name, field_name)
    if spec is None:
        raise ValueError(
            f"{package_name}.{field_name} is not a declared input field, so how its "
            f"records on one cell combine is not known; pass agg='sum' (a flux or a "
            f"conductance) or agg='mean' (an elevation or a property)."
        )
    return spec.agg


def get_package_input_field_names(package_name: str) -> list[str]:
    """Return the known input field names for one package (empty if unknown)."""

    spec = get_package_explorer_spec(package_name)
    return list(spec.inputs) if spec is not None else []


def get_package_result_spec(package_name: str, result_name: str) -> ResultSpec | None:
    """Return registry metadata for one result field, if it is known."""

    spec = get_package_explorer_spec(package_name)
    if spec is None:
        return None
    return spec.results.get(str(result_name).lower())


def get_default_group_compare_colorscale() -> list[list[object]]:
    """Return the default diverging colorscale for grouped difference maps.

    Explicit STOPS, not the name ``'RdBu'`` this used to return: the stops state
    which end is red, which is the policy, and a bare name cannot carry it.
    (Until ledger 253 a hand table also sent the name to matplotlib's REVERSED
    colormap, so a named diff map read "lower" on one backend and "higher" on
    the other -- ledger 69/70. Names now render from Plotly's own stops on both
    backends; the orientation reason stands.)
    """

    return [list(stop) for stop in _GROUP_COMPARE_DEFAULT_COLORSCALE]


__all__ = [
    "FieldSpec",
    "PackageExplorerSpec",
    "ResultSpec",
    "SeveralPackagesError",
    "get_default_budget_term",
    "get_default_group_compare_colorscale",
    "get_default_package_colorscale",
    "get_default_package_value_column",
    "get_package_explorer_spec",
    "get_package_input_field_names",
    "get_package_input_field_spec",
    "get_package_result_spec",
]


def _package_type_name(package_type: str) -> str:
    """An MF6 package type as a name file spells it: ``RCHA``/``EVTA`` are RCH/EVT.

    FloPy gives array-form recharge and ET classes of their own; MF6 has ONE recharge
    and ONE ET package, whose array form is an option inside the file.
    """

    upper = str(package_type).upper()
    return {"RCHA": "RCH", "EVTA": "EVT"}.get(upper, upper)


# -- a model's packages, found by TYPE (ledger 352) ----------------------------------
# FloPy's `get_package(name)` answers a package NAMED like its argument before it looks
# at types, so on a model holding LAK packages named `lak` and `lak_2` a lookup of "lak"
# returned `lak` alone -- `model.fl.lak`, `model.package("lak")`, and through them every
# LAK view, which read `lak`'s lake as the model's in silence. Every "the packages of
# this type" read goes through `packages_of_type` / `model_packages_of_type`, and every
# view that reads ONE package through `one_package`.


class SeveralPackagesError(LookupError):
    """A view that reads ONE package was asked of a model carrying several of its type.

    Raised by the LAK, SFR and UZF views (``model.lak.results.stage.get()``, ...) and
    their readers (``model.outputs.lak``) on a model with two or more packages of the
    type and none named. The message names every package and the spelling that picks
    one -- ``model.lak(package="lak_2")``.

    A ``LookupError``, deliberately NOT a ``ValueError``: the tiers that tolerate a
    MISSING package catch ``ValueError`` among others
    (``package_tables.PACKAGE_TABLE_UNAVAILABLE``, the group and diff tiers), and one
    catching this would compare a two-LAK model as having no LAK at all, in silence.

    Examples
    --------
    >>> try:
    ...     model.lak.results.stage.get()
    ... except mf.SeveralPackagesError as error:
    ...     print(error)          # names 'lak', 'lak_2' and model.lak(package=...)
    """


def packages_of_type(owner, package_type: str) -> list:
    """Every package of TYPE ``package_type`` a FloPy model or simulation holds.

    Asked by type ALONE -- never by name first, as FloPy's ``get_package`` does -- and
    by MF6's type: FloPy's array-form ``rcha``/``evta`` classes are RCH/EVT
    (:func:`_package_type_name`), on either side, so ``"rch"`` finds a list RCH and an
    array one both.

    Parameters
    ----------
    owner : flopy.mf6.MFModel or flopy.mf6.MFSimulation
        What holds the packages; its ``get_package()`` (no name) lists them all.
    package_type : str
        An MF6 package type, any case -- ``"lak"``, ``"RCH"``, ``"rcha"``, ``"ims"``.

    Returns
    -------
    list
        The FloPy packages of that type, in the order the owner holds them; ``[]``
        when it holds none.

    Examples
    --------
    >>> [p.package_name for p in packages_of_type(model.gwf, "lak")]
    ['lak', 'lak_2']
    """

    wanted = _package_type_name(package_type)
    held = owner.get_package() or ()
    return [
        package
        for package in held
        if _package_type_name(getattr(package, "package_type", "") or "") == wanted
    ]


def model_packages_of_type(model, package_type: str) -> list:
    """Every FloPy package of TYPE ``package_type`` in a myflopy model; ``[]`` when none.

    The one rule behind ``model.fl.<pkg>``, ``model.package(<type>)``, the
    one-package views (:func:`one_package`), the tables that read every package of a
    type, and the observation targets' package lookup.

    Parameters
    ----------
    model : SimulationBase or duck-typed model
        A myflopy model: its ``_packages_of_type`` loads a file-backed run's packages
        first. A test double without it answers through ``model.package(type)``
        (one package, a list, or ``AttributeError`` for none).
    package_type : str
        An MF6 package type, any case (see :func:`packages_of_type`).

    Returns
    -------
    list
        The model's FloPy packages of that type, in model order.
    """

    own = getattr(model, "_packages_of_type", None)
    if callable(own):
        return list(own(package_type))
    try:
        found = model.package(package_type)
    except AttributeError:
        return []
    if found is None:
        return []
    return list(found) if isinstance(found, (list, tuple)) else [found]


def package_name_of(package) -> str:
    """A FloPy package's own name (``package_name``), as FloPy holds it; ``""`` for none."""

    return str(getattr(package, "package_name", "") or "")


def package_named(model, package_type: str, package: str, *, everyone: list | None = None):
    """The model's package of TYPE ``package_type`` NAMED ``package`` (any case).

    Parameters
    ----------
    model : SimulationBase or duck-typed model
        The model (see :func:`model_packages_of_type`).
    package_type : str
        The MF6 package type the name must belong to.
    package : str
        A package NAME, FloPy's ``package_name``, compared stripped and lower-cased.
        Never a type: on a model holding ``lak`` and ``lak_2``, ``"lak"`` is the
        package named ``lak``.
    everyone : list or None, optional
        The model's packages of the type when the caller already has them. Default
        ``None`` looks them up.

    Returns
    -------
    flopy.mf6.mfpackage.MFPackage
        The one package of the type bearing that name.

    Raises
    ------
    ValueError
        When no package of the type bears the name, naming the ones that do exist.
    """

    everyone = model_packages_of_type(model, package_type) if everyone is None else everyone
    wanted = str(package).strip().lower()
    for candidate in everyone:
        if package_name_of(candidate).strip().lower() == wanted:
            return candidate
    held = ", ".join(repr(package_name_of(p)) for p in everyone) or "none"
    raise ValueError(
        f"package={package!r} names no {str(package_type).upper()} package of model "
        f"{getattr(model, 'name', '?')!r}; its {str(package_type).upper()} packages: {held}."
    )


def checked_package_name(model, package_type: str, package, *, spelling: str) -> str:
    """A ``package=`` argument checked where it is given: the named package's own name.

    Behind ``model.lak(package=...)`` and ``model.outputs.lak(package=...)`` (ledger
    352), so a typo raises at the call rather than at the first read.

    Parameters
    ----------
    model : SimulationBase or duck-typed model
        The model (see :func:`model_packages_of_type`).
    package_type : str
        The MF6 package type the name must belong to.
    package : str
        The argument: a package NAME, any case.
    spelling : str
        The call, for the ``TypeError`` -- ``"model.lak(package=...)"``.

    Returns
    -------
    str
        The package's ``package_name`` as FloPy holds it.

    Raises
    ------
    TypeError
        When ``package`` is not a non-empty string.
    ValueError
        When no package of the type bears the name (:func:`package_named`).
    """

    if not isinstance(package, str) or not package.strip():
        raise TypeError(
            f"{spelling} takes a {str(package_type).upper()} package's name, not {package!r}."
        )
    return package_name_of(package_named(model, package_type, package))


def one_package(model, package_type: str, package: str | None = None, *, spelling: str | None = None):
    """The ONE FloPy package of ``package_type`` a one-package view reads.

    A LAK, SFR or UZF view reads one package's lakes, reaches or UZF cells -- they are
    numbered within their package -- so on a model with several it must be told which.
    It never picks one itself (ledger 352).

    Parameters
    ----------
    model : SimulationBase or duck-typed model
        The model (see :func:`model_packages_of_type`).
    package_type : str
        ``"lak"``, ``"sfr"`` or ``"uzf"`` -- any MF6 package type, any case.
    package : str or None, optional
        The package's NAME (FloPy's ``package_name``), any case: that package, or
        ``ValueError`` naming the ones the model has. Default ``None``: the model's
        only package of the type.
    spelling : str or None, optional
        How a caller names the package, quoted by the refusal with ``{name}`` standing
        for each package's name. Default ``"model.<type>(package={name!r})"`` -- a
        view's spelling; ``model.outputs`` passes its own.

    Returns
    -------
    flopy.mf6.mfpackage.MFPackage
        The package.

    Raises
    ------
    AttributeError
        When the model has no package of the type -- what FloPy's lookup raised, so a
        tier that tolerates an absent package still does.
    SeveralPackagesError
        When ``package`` is ``None`` and the model carries several packages of the
        type, naming each and the spelling that picks it.
    ValueError
        When ``package`` names no package of the type.

    Examples
    --------
    >>> one_package(model, "lak").package_name                  # one LAK
    'lak'
    >>> one_package(model, "lak", "LAK_2").package_name         # any case
    'lak_2'
    """

    everyone = model_packages_of_type(model, package_type)
    if package is not None:
        return package_named(model, package_type, package, everyone=everyone)
    kind = str(package_type).upper()
    if not everyone:
        raise AttributeError(
            f"Model {getattr(model, 'name', '?')!r} has no {kind} package."
        )
    if len(everyone) == 1:
        return everyone[0]
    spelling = spelling or f"model.{str(package_type).lower()}(package={{name!r}})"
    names = [package_name_of(p) for p in everyone]
    features = _NUMBERED_FEATURES.get(kind, "features")
    raise SeveralPackagesError(
        f"Model {getattr(model, 'name', '?')!r} carries {len(everyone)} {kind} packages "
        f"({', '.join(repr(name) for name in names)}), and this reads ONE package's "
        f"{features}, which each package numbers from 1. Name the package: "
        + " or ".join(spelling.format(name=name) for name in names)
        + "."
    )


def package_outputs(model, package_type: str, package: str | None = None):
    """``model.outputs.<type>`` for ONE package: bound to ``package`` when one is named.

    The unbound namespace reads the model's only package of the type, so a reader
    that has already refused several (:func:`one_package`, whose refusal names the
    VIEW's spelling) reads the files exactly as it did before ledger 352.

    Parameters
    ----------
    model : SimulationBase or duck-typed model
        The model; its ``outputs`` namespace is read.
    package_type : str
        ``"lak"``, ``"sfr"`` or ``"uzf"``.
    package : str or None, optional
        The package's NAME, or ``None`` for the model's only one.

    Returns
    -------
    LakOutputData or SFROutputData or UzfOutputData
        ``model.outputs.<type>``, or ``model.outputs.<type>(package=package)``.
    """

    namespace = getattr(model.outputs, package_type)
    return namespace if package is None else namespace(package=package)


#: What a one-package view's rows are, per type, for :func:`one_package`'s refusal.
_NUMBERED_FEATURES = {"LAK": "lakes", "SFR": "reaches", "UZF": "UZF cells"}


def feature_boundnames(packagedata) -> dict[str, list[int]]:
    """Which features of a LAK, SFR or UZF package each boundname names, from its ``packagedata``.

    The ONE reading of a boundname (ledger 558): stripped and lower-cased,
    as MF6 compares one, a blank or missing one naming nothing. An observation's
    boundname (``observations._shared``) and a ``style="feature"`` parameter's
    (``pest.native_parameters``) are both resolved through it.

    Parameters
    ----------
    packagedata : numpy.recarray or pandas.DataFrame or None
        The package's ``packagedata`` as FloPy's ``get_data()`` returns it (a record
        array) or as a table; its FIRST field is the feature's ZERO-based number (a
        lake's ``ifno``/``lakeno``, a reach's ``ifno``), and a ``boundname`` field, when
        present, its name. ``None`` names nothing.

    Returns
    -------
    dict of str to list of int
        Each lower-cased, stripped boundname to the zero-based numbers of the
        features carrying it, in ``packagedata`` order; empty when there is no
        ``boundname`` field.
    """

    if packagedata is None:
        return {}
    import pandas as pd

    table = packagedata if isinstance(packagedata, pd.DataFrame) else pd.DataFrame(packagedata)
    if table.empty or "boundname" not in table.columns:
        return {}
    named: dict[str, list[int]] = {}
    for number, bound in zip(table.iloc[:, 0], table["boundname"], strict=True):
        if bound is None or (isinstance(bound, float) and bound != bound):
            continue
        text = str(bound).strip().lower()
        if text:
            named.setdefault(text, []).append(int(number))
    return named
