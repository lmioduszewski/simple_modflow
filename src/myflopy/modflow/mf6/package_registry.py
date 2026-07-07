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
class PackageExplorerSpec:
    """Registry metadata for one MF6 package explorer."""

    name: str
    kind: str = "cell_stress"
    default_input: str | None = None
    colorscale: str | None = None
    inputs: dict[str, FieldSpec] = field(default_factory=dict)
    results: dict[str, ResultSpec] = field(default_factory=dict)


_PACKAGE_EXPLORER_SPECS: dict[str, PackageExplorerSpec] = {
    # Colorscale policy: diverging red/white/blue is reserved for signed,
    # gaining/losing "q"-like fields and diff maps; everything else uses the
    # house brown-to-blue "earth" scale (the mounding-figure default).
    "rch": PackageExplorerSpec(
        name="rch",
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
    "wel": PackageExplorerSpec(
        name="wel",
        default_input="q",
        colorscale="RdBu",
        inputs={
            "q": FieldSpec("q", label="Well flow", colorscale="RdBu"),
        },
        results={
            "q": ResultSpec("q", budget_text="WEL", value_name="q", colorscale="RdBu"),
        },
    ),
    "uzf": PackageExplorerSpec(
        name="uzf",
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
        kind="surface_water",
        results={
            "q": ResultSpec("q", budget_text="SFR", value_name="q", colorscale="RdBu"),
        },
    ),
    "lak": PackageExplorerSpec(
        name="lak",
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
