"""Registry metadata for package input and result explorers."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class FieldSpec:
    """Registry metadata for one package input field."""

    name: str
    label: str | None = None
    colorscale: str | None = None
    fill_value: float = 0.0
    agg: str = "sum"


@dataclass(frozen=True)
class ResultSpec:
    """Registry metadata for one package result term."""

    name: str
    budget_text: str
    value_name: str = "q"
    label: str | None = None
    colorscale: str | None = None
    diverging: bool = True


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
    #: MF6 reports this package's budget node numbers 1-based, so readers must
    #: subtract one. True for the cell-stress list BCs; the advanced packages
    #: report feature numbers instead. Got this wrong for years -- see the
    #: budget off-by-one fixed in 4.7.1.
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
            "recharge": FieldSpec("recharge", label="Recharge", colorscale="earth"),
        },
        results={
            "q": ResultSpec("q", budget_text="RCH", value_name="q", colorscale="RdBu"),
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
            "head": FieldSpec("head", label="Constant head", colorscale="earth"),
        },
        results={
            "q": ResultSpec("q", budget_text="CHD", value_name="q", colorscale="RdBu"),
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
            "elev": FieldSpec("elev", label="Drain elevation", colorscale="earth"),
            "cond": FieldSpec("cond", label="Drain conductance", colorscale="earth"),
        },
        results={
            "q": ResultSpec("q", budget_text="DRN", value_name="q", colorscale="RdBu"),
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
            "bhead": FieldSpec("bhead", label="Boundary head", colorscale="earth"),
            "cond": FieldSpec("cond", label="Boundary conductance", colorscale="earth"),
        },
        results={
            "q": ResultSpec("q", budget_text="GHB", value_name="q", colorscale="RdBu"),
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
            "stage": FieldSpec("stage", label="River stage", colorscale="earth"),
            "cond": FieldSpec("cond", label="Riverbed conductance", colorscale="earth"),
            "rbot": FieldSpec("rbot", label="River bottom", colorscale="earth"),
        },
        results={
            "q": ResultSpec("q", budget_text="RIV", value_name="q", colorscale="RdBu"),
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
            "q": FieldSpec("q", label="Well flow", colorscale="RdBu"),
        },
        results={
            "q": ResultSpec("q", budget_text="WEL", value_name="q", colorscale="RdBu"),
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
            "surface": FieldSpec("surface", label="ET surface", colorscale="earth"),
            "rate": FieldSpec("rate", label="Maximum ET rate", colorscale="earth"),
            "depth": FieldSpec("depth", label="Extinction depth", colorscale="earth"),
        },
        results={
            "q": ResultSpec("q", budget_text="EVT", value_name="q", colorscale="RdBu"),
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
            "finf": FieldSpec("finf", label="UZF infiltration", colorscale="earth"),
            "pet": FieldSpec("pet", label="Potential evapotranspiration", colorscale="earth"),
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
            "q": ResultSpec("q", budget_text="SFR", value_name="q", colorscale="RdBu"),
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
            "q": ResultSpec("q", budget_text="GWF", value_name="q", colorscale="RdBu"),
        },
    ),
}

_GROUP_COMPARE_DEFAULT_COLORSCALE = "RdBu"


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


def get_default_group_compare_colorscale() -> str:
    """Return the default diverging colorscale for grouped difference maps."""

    return _GROUP_COMPARE_DEFAULT_COLORSCALE


__all__ = [
    "FieldSpec",
    "PackageExplorerSpec",
    "ResultSpec",
    "get_default_budget_term",
    "get_default_group_compare_colorscale",
    "get_default_package_colorscale",
    "get_default_package_value_column",
    "get_package_explorer_spec",
    "get_package_input_field_names",
    "get_package_input_field_spec",
    "get_package_result_spec",
]
