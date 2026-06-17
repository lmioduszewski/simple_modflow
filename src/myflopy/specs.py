"""Composable specifications for building MODFLOW 6 simulations.

The classes in this module are intentionally small. A spec stores a named
builder and its keyword arguments, while the builder remains ordinary Python.
This keeps model assembly explicit and makes package alternatives easy to
replace without introducing a framework around FloPy.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum
from html import escape
from pathlib import Path
from typing import Any, Callable, Iterable, TypeAlias, TypeVar, cast

import flopy


Builder = Callable[..., Any]
Hook = Callable[[Any, dict[str, Any], "ModelContext"], Any]
Mf6Simulation: TypeAlias = flopy.mf6.MFSimulation
GwfModel: TypeAlias = flopy.mf6.ModflowGwf
GwtModel: TypeAlias = flopy.mf6.ModflowGwt
GweModel: TypeAlias = flopy.mf6.ModflowGwe
PrtModel: TypeAlias = flopy.mf6.ModflowPrt
Mf6Model: TypeAlias = GwfModel | GwtModel | GweModel | PrtModel
_ModelT = TypeVar("_ModelT", bound=Mf6Model)


class ModelType(str, Enum):
    """MODFLOW 6 model types supported by the default model builders."""

    GWF = "gwf"
    GWT = "gwt"
    GWE = "gwe"
    PRT = "prt"


_DEFAULT_MODEL_BUILDERS: dict[ModelType, Builder] = {
    ModelType.GWF: flopy.mf6.ModflowGwf,
    ModelType.GWT: flopy.mf6.ModflowGwt,
    ModelType.GWE: flopy.mf6.ModflowGwe,
    ModelType.PRT: flopy.mf6.ModflowPrt,
}


def _require_unique_names(items: Iterable[Any], *, item_type: str) -> None:
    """Raise when a collection contains duplicate ``name`` values."""

    names = [item.name for item in items]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"Duplicate {item_type} names: {', '.join(duplicates)}")


def _builder_label(builder: Builder | None) -> str:
    """Return a short readable label for a builder callable."""

    if builder is None:
        return "default"
    module = getattr(builder, "__module__", "")
    name = getattr(builder, "__name__", builder.__class__.__name__)
    if module.startswith("flopy."):
        return f"flopy.{name}"
    return name


def _preview(values: Iterable[Any], *, limit: int = 4) -> str:
    """Return a short comma-separated preview of an iterable."""

    items = list(values)
    shown = ", ".join(repr(value) for value in items[:limit])
    if len(items) > limit:
        shown = f"{shown}, ..."
    return shown


def _shape_label(value: Any) -> str | None:
    """Return a compact shape label for array/table-like values."""

    shape = getattr(value, "shape", None)
    if shape is None:
        return None
    class_name = value.__class__.__name__
    return f"{class_name}(shape={tuple(shape)!r})"


def _summarize_mapping(value: dict[Any, Any]) -> str:
    """Return a compact summary for mappings, including stress-period data."""

    keys = list(value.keys())
    if not keys:
        return "dict(len=0)"
    if all(isinstance(key, int) for key in keys):
        record_counts = [
            len(records)
            for records in value.values()
            if isinstance(records, (list, tuple))
        ]
        if record_counts:
            return (
                "period_data("
                f"periods={len(keys)}, "
                f"records={sum(record_counts)}, "
                f"max_per_period={max(record_counts)}"
                ")"
            )
    return f"dict(len={len(keys)}, keys=[{_preview(keys)}])"


def _summarize_value(value: Any) -> str:
    """Return a compact, notebook-friendly value summary."""

    shape = _shape_label(value)
    if shape is not None:
        return shape
    if isinstance(value, dict):
        return _summarize_mapping(value)
    if isinstance(value, (list, tuple)):
        if not value:
            return f"{type(value).__name__}(len=0)"
        return f"{type(value).__name__}(len={len(value)}, preview=[{_preview(value)}])"
    if isinstance(value, str):
        return repr(value if len(value) <= 80 else f"{value[:77]}...")
    text = repr(value)
    return text if len(text) <= 80 else f"{text[:77]}..."


def _html_table(rows: Iterable[tuple[str, str]], *, empty: str = "None") -> str:
    """Return a small HTML table for notebook rich display."""

    rows = list(rows)
    if not rows:
        return f"<p>{escape(empty)}</p>"
    body = "\n".join(
        "<tr>"
        f"<td><code>{escape(key)}</code></td>"
        f"<td>{escape(value)}</td>"
        "</tr>"
        for key, value in rows
    )
    return (
        "<table>"
        "<thead><tr><th>Name</th><th>Summary</th></tr></thead>"
        f"<tbody>{body}</tbody>"
        "</table>"
    )


def _require_model_type(model: Any, model_type: type[_ModelT], label: str) -> _ModelT:
    """Return ``model`` as ``model_type`` or explain the mismatch clearly."""

    if not isinstance(model, model_type):
        raise TypeError(
            f"Built model is not a {label}. "
            f"Found {model.__class__.__module__}.{model.__class__.__name__}."
        )
    return cast(_ModelT, model)


@dataclass(frozen=True, slots=True)
class PackageSpec:
    """Reusable instructions for building one package on a model or simulation.

    Builders receive the target as their first positional argument followed by
    the values in ``options``. A builder can be a FloPy package class or a
    small project-specific function.
    """

    name: str
    builder: Builder
    options: dict[str, Any] = field(default_factory=dict)
    enabled: bool = True
    requires: tuple[str, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "requires", tuple(self.requires))

    def build(self, target: Any) -> Any | None:
        """Build the package on ``target``, or return ``None`` when disabled."""

        if not self.enabled:
            return None
        return self.builder(target, **self.options)

    def with_options(self, **overrides: Any) -> PackageSpec:
        """Return a copy with selected options added or replaced."""

        return replace(self, options={**self.options, **overrides})

    def disabled(self) -> PackageSpec:
        """Return a disabled copy of this package specification."""

        return replace(self, enabled=False)

    def with_metadata(self, **updates: Any) -> PackageSpec:
        """Return a copy with selected provenance metadata added or replaced."""

        return replace(self, metadata={**self.metadata, **updates})

    @property
    def option_summary(self) -> dict[str, str]:
        """Return a compact summary of package options for display."""

        return {
            key: _summarize_value(value)
            for key, value in self.options.items()
        }

    def __repr__(self) -> str:
        status = "enabled" if self.enabled else "disabled"
        pieces = [
            f"name={self.name!r}",
            f"builder={_builder_label(self.builder)!r}",
            f"status={status!r}",
            f"options={self.option_summary!r}",
        ]
        if self.requires:
            pieces.append(f"requires={self.requires!r}")
        if self.metadata:
            pieces.append(f"metadata_keys={tuple(self.metadata)!r}")
        return f"PackageSpec({', '.join(pieces)})"

    def _repr_html_(self) -> str:
        """Return a compact rich representation for notebooks."""

        metadata = ""
        if self.metadata:
            metadata = (
                "<p><strong>Metadata:</strong> "
                f"{escape(', '.join(self.metadata.keys()))}</p>"
            )
        requires = ""
        if self.requires:
            requires = (
                "<p><strong>Requires:</strong> "
                f"{escape(', '.join(self.requires))}</p>"
            )
        return (
            "<div>"
            f"<h4>PackageSpec: <code>{escape(self.name)}</code></h4>"
            f"<p><strong>Builder:</strong> {escape(_builder_label(self.builder))} "
            f"<strong>Status:</strong> {'enabled' if self.enabled else 'disabled'}</p>"
            f"{requires}"
            f"{metadata}"
            f"{_html_table(self.option_summary.items(), empty='No options')}"
            "</div>"
        )


@dataclass(frozen=True, slots=True)
class ModelContext:
    """Domain information carried beside a model specification.

    Context is intentionally separate from FloPy package options. It gives
    data builders and post-build hooks access to geometry, surfaces, dates,
    domain information, and project-specific metadata.
    """

    grid: Any = None
    surfaces: Any = None
    dates: Any = None
    domain: Any = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def with_metadata(self, **updates: Any) -> ModelContext:
        """Return a copy with selected metadata added or replaced."""

        return replace(self, metadata={**self.metadata, **updates})


@dataclass(frozen=True, slots=True)
class PostBuildHook:
    """Named callback run after every package on a model has been built."""

    name: str
    callback: Hook

    def run(self, model: Any, packages: dict[str, Any], context: ModelContext) -> Any:
        """Run the callback with the completed FloPy model and its context."""

        return self.callback(model, packages, context)


@dataclass(frozen=True, slots=True)
class BuiltModel:
    """A built FloPy model and the package objects created from its spec."""

    model: Mf6Model
    packages: dict[str, Any]
    context: ModelContext = field(default_factory=ModelContext)
    hook_results: dict[str, Any] = field(default_factory=dict)
    simulation: Mf6Simulation | None = None

    @property
    def sim(self) -> Mf6Simulation:
        """Return the parent FloPy simulation for this built model."""

        if self.simulation is not None:
            return self.simulation
        model_simulation = getattr(self.model, "sim", None)
        if isinstance(model_simulation, Mf6Simulation):
            return model_simulation
        raise AttributeError("This built model does not have a parent MFSimulation.")

    @property
    def gwf(self) -> GwfModel:
        """Return the model as a typed FloPy GWF model."""

        return self.as_gwf()

    @property
    def gwt(self) -> GwtModel:
        """Return the model as a typed FloPy GWT model."""

        return self.as_gwt()

    @property
    def gwe(self) -> GweModel:
        """Return the model as a typed FloPy GWE model."""

        return self.as_gwe()

    @property
    def prt(self) -> PrtModel:
        """Return the model as a typed FloPy PRT model."""

        return self.as_prt()

    def as_gwf(self) -> GwfModel:
        """Return the model as ``flopy.mf6.ModflowGwf``."""

        return _require_model_type(self.model, GwfModel, "GWF model")

    def as_gwt(self) -> GwtModel:
        """Return the model as ``flopy.mf6.ModflowGwt``."""

        return _require_model_type(self.model, GwtModel, "GWT model")

    def as_gwe(self) -> GweModel:
        """Return the model as ``flopy.mf6.ModflowGwe``."""

        return _require_model_type(self.model, GweModel, "GWE model")

    def as_prt(self) -> PrtModel:
        """Return the model as ``flopy.mf6.ModflowPrt``."""

        return _require_model_type(self.model, PrtModel, "PRT model")

    def package(self, name: str) -> Any:
        """Return one built package by name."""

        try:
            return self.packages[name]
        except KeyError as error:
            raise KeyError(f"Built model has no package named '{name}'.") from error


@dataclass(frozen=True, slots=True)
class ModelSpec:
    """Definition of one GWF, GWT, GWE, or PRT model.

    Package names are unique within the model. Calling :meth:`with_package`
    replaces an existing package with the same name, which makes scenario
    variants explicit and inexpensive to create.
    """

    name: str
    model_type: ModelType | str
    packages: tuple[PackageSpec, ...] = ()
    options: dict[str, Any] = field(default_factory=dict)
    builder: Builder | None = None
    context: ModelContext = field(default_factory=ModelContext)
    hooks: tuple[PostBuildHook, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "model_type", ModelType(self.model_type))
        object.__setattr__(self, "packages", tuple(self.packages))
        object.__setattr__(self, "hooks", tuple(self.hooks))
        _require_unique_names(self.packages, item_type="package")
        _require_unique_names(self.hooks, item_type="post-build hook")
        available = {package.name for package in self.packages if package.enabled}
        package_positions = {
            package.name: index
            for index, package in enumerate(self.packages)
            if package.enabled
        }
        for package in self.packages:
            missing = sorted(set(package.requires) - available)
            if package.enabled and missing:
                raise ValueError(
                    f"Package '{package.name}' requires missing packages: {', '.join(missing)}"
                )
            if package.enabled:
                later = sorted(
                    name
                    for name in package.requires
                    if package_positions[name] > package_positions[package.name]
                )
                if later:
                    raise ValueError(
                        f"Package '{package.name}' must be declared after required packages: "
                        f"{', '.join(later)}"
                    )

    def with_package(self, package: PackageSpec) -> ModelSpec:
        """Return a copy with ``package`` added or replaced by name."""

        packages = list(self.packages)
        for index, existing in enumerate(packages):
            if existing.name == package.name:
                packages[index] = package
                break
        else:
            packages.append(package)
        return replace(self, packages=tuple(packages))

    def package(self, name: str) -> PackageSpec:
        """Return one package specification by name."""

        for package in self.packages:
            if package.name == name:
                return package
        raise KeyError(f"Model '{self.name}' has no package named '{name}'.")

    def without_package(self, name: str) -> ModelSpec:
        """Return a copy without the named package."""

        return replace(
            self,
            packages=tuple(package for package in self.packages if package.name != name),
        )

    def with_options(self, **overrides: Any) -> ModelSpec:
        """Return a copy with selected model options added or replaced."""

        return replace(self, options={**self.options, **overrides})

    def with_context(self, context: ModelContext) -> ModelSpec:
        """Return a copy carrying the supplied model context."""

        return replace(self, context=context)

    def with_hook(self, hook: PostBuildHook) -> ModelSpec:
        """Return a copy with a post-build hook added or replaced by name."""

        hooks = list(self.hooks)
        for index, existing in enumerate(hooks):
            if existing.name == hook.name:
                hooks[index] = hook
                break
        else:
            hooks.append(hook)
        return replace(self, hooks=tuple(hooks))

    def build(self, simulation: Any) -> BuiltModel:
        """Build the model and all enabled packages in declaration order."""

        builder = self.builder or _DEFAULT_MODEL_BUILDERS[self.model_type]
        model = builder(simulation, modelname=self.name, **self.options)
        model.myflopy_context = self.context
        built_packages = {
            package.name: built
            for package in self.packages
            if (built := package.build(model)) is not None
        }
        hook_results = {
            hook.name: hook.run(model, built_packages, self.context)
            for hook in self.hooks
        }
        return BuiltModel(
            model=model,
            packages=built_packages,
            context=self.context,
            hook_results=hook_results,
            simulation=simulation,
        )

    @property
    def package_names(self) -> tuple[str, ...]:
        """Return enabled package names in declaration order."""

        return tuple(package.name for package in self.packages if package.enabled)

    def __repr__(self) -> str:
        return (
            "ModelSpec("
            f"name={self.name!r}, "
            f"model_type={self.model_type.value!r}, "
            f"packages={self.package_names!r}, "
            f"options={_summarize_mapping(self.options)!r}, "
            f"hooks={tuple(hook.name for hook in self.hooks)!r}"
            ")"
        )

    def _repr_html_(self) -> str:
        """Return a compact rich representation for notebooks."""

        package_rows = [
            (
                package.name,
                f"{_builder_label(package.builder)}; "
                f"{len(package.options)} options"
                + ("; disabled" if not package.enabled else "")
                + (f"; requires {', '.join(package.requires)}" if package.requires else ""),
            )
            for package in self.packages
        ]
        hooks = ", ".join(hook.name for hook in self.hooks) or "None"
        return (
            "<div>"
            f"<h4>ModelSpec: <code>{escape(self.name)}</code></h4>"
            f"<p><strong>Type:</strong> {escape(self.model_type.value)} "
            f"<strong>Packages:</strong> {len(self.packages)} "
            f"<strong>Hooks:</strong> {escape(hooks)}</p>"
            f"{_html_table(package_rows, empty='No packages')}"
            "</div>"
        )


@dataclass(frozen=True, slots=True)
class ExchangeSpec:
    """Definition of an exchange between two or more models.

    The exchange builder receives the simulation, a tuple containing the built
    model objects named by ``models``, and then the configured keyword options.
    """

    name: str
    builder: Builder
    models: tuple[str, ...]
    options: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "models", tuple(self.models))
        if len(self.models) < 2:
            raise ValueError("An exchange must reference at least two models.")

    def build(self, simulation: Any, built_models: dict[str, BuiltModel]) -> Any:
        """Build the exchange after resolving its model names."""

        missing = [name for name in self.models if name not in built_models]
        if missing:
            raise KeyError(
                f"Exchange '{self.name}' references unknown models: {', '.join(missing)}"
            )
        models = tuple(built_models[name].model for name in self.models)
        return self.builder(simulation, models, **self.options)

    def with_options(self, **overrides: Any) -> ExchangeSpec:
        """Return a copy with selected exchange options added or replaced."""

        return replace(self, options={**self.options, **overrides})

    def __repr__(self) -> str:
        return (
            "ExchangeSpec("
            f"name={self.name!r}, "
            f"builder={_builder_label(self.builder)!r}, "
            f"models={self.models!r}, "
            f"options={_summarize_mapping(self.options)!r}"
            ")"
        )


@dataclass(frozen=True, slots=True)
class BuiltSimulation:
    """A built FloPy simulation and the objects created from its specification."""

    simulation: Mf6Simulation
    models: dict[str, BuiltModel]
    packages: dict[str, Any]
    exchanges: dict[str, Any]

    @property
    def sim(self) -> Mf6Simulation:
        """Return the built FloPy ``MFSimulation``."""

        return self.simulation

    def model(self, name: str) -> Mf6Model:
        """Return one built FloPy model by name."""

        return self.built_model(name).model

    def built_model(self, name: str) -> BuiltModel:
        """Return one ``BuiltModel`` wrapper by name."""

        try:
            return self.models[name]
        except KeyError as error:
            raise KeyError(f"Built simulation has no model named '{name}'.") from error

    def gwf(self, name: str) -> GwfModel:
        """Return a named model as ``flopy.mf6.ModflowGwf``."""

        return self.built_model(name).as_gwf()

    def gwt(self, name: str) -> GwtModel:
        """Return a named model as ``flopy.mf6.ModflowGwt``."""

        return self.built_model(name).as_gwt()

    def gwe(self, name: str) -> GweModel:
        """Return a named model as ``flopy.mf6.ModflowGwe``."""

        return self.built_model(name).as_gwe()

    def prt(self, name: str) -> PrtModel:
        """Return a named model as ``flopy.mf6.ModflowPrt``."""

        return self.built_model(name).as_prt()

    def package(self, name: str) -> Any:
        """Return one built simulation-level package by name."""

        try:
            return self.packages[name]
        except KeyError as error:
            raise KeyError(
                f"Built simulation has no simulation-level package named '{name}'."
            ) from error


@dataclass(frozen=True, slots=True)
class SimulationSpec:
    """Definition of a complete simulation containing coupled model specs."""

    name: str
    models: tuple[ModelSpec, ...]
    packages: tuple[PackageSpec, ...] = ()
    exchanges: tuple[ExchangeSpec, ...] = ()
    options: dict[str, Any] = field(default_factory=dict)
    builder: Builder = flopy.mf6.MFSimulation
    workspace: Path | str | None = None
    executable: str = "mf6"
    run_name: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "models", tuple(self.models))
        object.__setattr__(self, "packages", tuple(self.packages))
        object.__setattr__(self, "exchanges", tuple(self.exchanges))
        _require_unique_names(self.models, item_type="model")
        _require_unique_names(self.packages, item_type="simulation package")
        _require_unique_names(self.exchanges, item_type="exchange")

    def with_model(self, model: ModelSpec) -> SimulationSpec:
        """Return a copy with ``model`` added or replaced by name."""

        models = list(self.models)
        for index, existing in enumerate(models):
            if existing.name == model.name:
                models[index] = model
                break
        else:
            models.append(model)
        return replace(self, models=tuple(models))

    def model(self, name: str) -> ModelSpec:
        """Return one model specification by name."""

        for model in self.models:
            if model.name == name:
                return model
        raise KeyError(f"Simulation '{self.name}' has no model named '{name}'.")

    def with_package(self, package: PackageSpec) -> SimulationSpec:
        """Return a copy with a simulation-level package added or replaced."""

        packages = list(self.packages)
        for index, existing in enumerate(packages):
            if existing.name == package.name:
                packages[index] = package
                break
        else:
            packages.append(package)
        return replace(self, packages=tuple(packages))

    def with_exchange(self, exchange: ExchangeSpec) -> SimulationSpec:
        """Return a copy with ``exchange`` added or replaced by name."""

        exchanges = list(self.exchanges)
        for index, existing in enumerate(exchanges):
            if existing.name == exchange.name:
                exchanges[index] = exchange
                break
        else:
            exchanges.append(exchange)
        return replace(self, exchanges=tuple(exchanges))

    def with_workspace(self, workspace: Path | str) -> SimulationSpec:
        """Return a copy with a default run workspace."""

        return replace(self, workspace=workspace)

    def build(self, workspace: Path | str | None = None, *, run_name: str | None = None):
        """Build this specification as a workflow ``Run``.

        Use :meth:`build_flopy` when you need only the raw in-memory FloPy
        objects without a run lifecycle wrapper.
        """

        from myflopy.workspace import Run

        run_workspace = workspace if workspace is not None else self.workspace
        if run_workspace is None:
            run_workspace = Path(self.run_name or self.name)
        run = Run(
            name=run_name or self.run_name or self.name,
            workspace=run_workspace,
            spec=self,
            executable=self.executable,
        )
        run.build()
        return run

    def build_flopy(self, workspace: Path | str | None = None) -> BuiltSimulation:
        """Build and return the raw FloPy simulation objects."""

        options = dict(self.options)
        run_workspace = workspace if workspace is not None else self.workspace
        if run_workspace is not None:
            options["sim_ws"] = str(run_workspace)
        simulation = self.builder(sim_name=self.name, **options)

        simulation_packages = {
            package.name: built
            for package in self.packages
            if (built := package.build(simulation)) is not None
        }
        built_models = {model.name: model.build(simulation) for model in self.models}
        built_exchanges = {
            exchange.name: exchange.build(simulation, built_models)
            for exchange in self.exchanges
        }
        return BuiltSimulation(
            simulation=simulation,
            models=built_models,
            packages=simulation_packages,
            exchanges=built_exchanges,
        )

    @property
    def model_names(self) -> tuple[str, ...]:
        """Return model names in declaration order."""

        return tuple(model.name for model in self.models)

    @property
    def package_names(self) -> tuple[str, ...]:
        """Return simulation-level package names in declaration order."""

        return tuple(package.name for package in self.packages if package.enabled)

    def __repr__(self) -> str:
        return (
            "SimulationSpec("
            f"name={self.name!r}, "
            f"models={self.model_names!r}, "
            f"packages={self.package_names!r}, "
            f"exchanges={tuple(exchange.name for exchange in self.exchanges)!r}"
            ")"
        )

    def _repr_html_(self) -> str:
        """Return a compact rich representation for notebooks."""

        model_rows = [
            (
                model.name,
                f"{model.model_type.value}; packages={len(model.packages)}",
            )
            for model in self.models
        ]
        package_rows = [
            (
                package.name,
                f"{_builder_label(package.builder)}; {len(package.options)} options",
            )
            for package in self.packages
            if package.enabled
        ]
        exchange_rows = [
            (
                exchange.name,
                f"{_builder_label(exchange.builder)}; models={', '.join(exchange.models)}",
            )
            for exchange in self.exchanges
        ]
        return (
            "<div>"
            f"<h4>SimulationSpec: <code>{escape(self.name)}</code></h4>"
            f"<p><strong>Models:</strong> {len(self.models)} "
            f"<strong>Simulation Packages:</strong> {len(self.packages)} "
            f"<strong>Exchanges:</strong> {len(self.exchanges)}</p>"
            "<h5>Models</h5>"
            f"{_html_table(model_rows, empty='No models')}"
            "<h5>Simulation Packages</h5>"
            f"{_html_table(package_rows, empty='No simulation packages')}"
            "<h5>Exchanges</h5>"
            f"{_html_table(exchange_rows, empty='No exchanges')}"
            "</div>"
        )
