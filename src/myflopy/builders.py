"""Small builders for common FloPy 3.10 simulation-level objects."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any, Protocol

import flopy

from myflopy.specs import PackageSpec


class PackageBuilder(Protocol):
    """Contract for configured objects that produce reusable package specs.

    A package builder receives all configuration when it is created. Its
    ``build`` method only validates that configuration and returns the spec.
    """

    def build(self) -> PackageSpec:
        """Return a package spec from the builder's stored configuration."""


def build_ims(simulation, *, models: Iterable[str], **options):
    """Build an IMS solver package and register it with the named models.

    Each solution needs a distinct file. Without one, multiple IMS packages
    collapse into a single solution, which MF6 rejects for coupled simulations
    (for example GWF-GWT). The filename defaults to the package name. This is the
    engine under ``mf.ims(...)``.

    Parameters
    ----------
    simulation
        The owning FloPy ``MFSimulation`` to attach the solver to.
    models : Iterable[str]
        Names of the models this IMS solves (registered to each).
    **options
        ``flopy.mf6.ModflowIms`` options (``complexity``, ``outer_maximum``,
        ``linear_acceleration``, ...); ``filename`` defaults to ``"<pname>.ims"``.

    Returns
    -------
    flopy.mf6.ModflowIms
        The constructed and registered IMS package.
    """

    models = list(models)
    base = options.get("pname") or (models[0] if models else "ims")
    options.setdefault("filename", f"{base}.ims")
    ims = flopy.mf6.ModflowIms(simulation, **options)
    simulation.register_ims_package(ims, models)
    return ims


# Model-type dispatch for shared core packages (plan 5.3A). ``ic``/``oc``/``disv``
# are the same package on GWF, GWT and GWE models but map to different FloPy
# classes; these module-level builders resolve the class from the built model's
# ``model_type`` (``"gwf6"``/``"gwt6"``/``"gwe6"``). They MUST stay module-level
# functions: ``PackageSpec`` serializes its builder by importable reference
# (``specs._callable_ref`` rejects lambdas/closures), and only the function is
# persisted -- the class tables below are re-read at build time, never serialized.
# PRT is intentionally absent (MF6 has no ``ModflowPrtic``; PRT builds its own dis
# via ``mf.prt``); ``npf``/``sto`` are GWF-only in MF6 and are not dispatched.
_IC_CLASSES = {
    "gwf6": flopy.mf6.ModflowGwfic,
    "gwt6": flopy.mf6.ModflowGwtic,
    "gwe6": flopy.mf6.ModflowGweic,
}
_OC_CLASSES = {
    "gwf6": flopy.mf6.ModflowGwfoc,
    "gwt6": flopy.mf6.ModflowGwtoc,
    "gwe6": flopy.mf6.ModflowGweoc,
}
_DISV_CLASSES = {
    "gwf6": flopy.mf6.ModflowGwfdisv,
    "gwt6": flopy.mf6.ModflowGwtdisv,
    "gwe6": flopy.mf6.ModflowGwedisv,
}


def _dispatch_core_class(table: dict, model, package_label: str):
    """Return the model-kind-specific FloPy class for a shared core package.

    Resolves the class from the built model's ``model_type`` (the FloPy instance
    attribute, ``"gwf6"``/``"gwt6"``/``"gwe6"``), raising a clear error for any
    kind the package does not support (e.g. PRT, which has no IC package).
    """

    model_type = getattr(model, "model_type", None)
    try:
        return table[model_type]
    except KeyError:
        supported = ", ".join(sorted(key.removesuffix("6") for key in table))
        raise ValueError(
            f"mf.{package_label} is not available for model type {model_type!r}; "
            f"supported model types: {supported}."
        ) from None


def build_ic(model, **options):
    """Build the IC package for the model's kind (GWF/GWT/GWE). Engine under ``mf.ic``."""

    return _dispatch_core_class(_IC_CLASSES, model, "ic")(model, **options)


def build_oc(model, **options):
    """Build the OC package for the model's kind (GWF/GWT/GWE). Engine under ``mf.oc``."""

    return _dispatch_core_class(_OC_CLASSES, model, "oc")(model, **options)


def build_disv(model, **options):
    """Build the DISV package for the model's kind (GWF/GWT/GWE). Engine under ``mf.disv``."""

    return _dispatch_core_class(_DISV_CLASSES, model, "disv")(model, **options)


def _build_two_model_exchange(constructor, simulation, models: tuple[Any, ...], **options):
    """Build a standard two-model MF6 exchange from resolved model objects."""

    if len(models) != 2:
        raise ValueError("A standard MF6 exchange requires exactly two models.")
    model_a, model_b = models
    return constructor(
        simulation,
        exgmnamea=model_a.name,
        exgmnameb=model_b.name,
        **options,
    )


def build_gwf_gwf_exchange(simulation, models: tuple[Any, ...], **options):
    """Build a GWF-GWF exchange coupling two flow models (e.g. for LGR/nested grids).

    Wires :class:`flopy.mf6.ModflowGwfgwf` between exactly two GWF models so they
    exchange flow across a shared interface. ``models`` is the resolved
    ``(model_a, model_b)`` pair; ``**options`` (e.g. ``exchangedata``, ``nexg``)
    pass through to FloPy. Usually emitted for you from an
    :class:`~myflopy.specs.ExchangeSpec` rather than called directly.

    Parameters
    ----------
    simulation
        The owning FloPy simulation.
    models
        The two flow models to couple, as ``(model_a, model_b)``.
    **options
        Exchange options forwarded to the FloPy constructor.
    """

    return _build_two_model_exchange(flopy.mf6.ModflowGwfgwf, simulation, models, **options)


def build_gwf_gwt_exchange(simulation, models: tuple[Any, ...], **options):
    """Build a GWF-GWT exchange coupling a flow model to a transport model.

    Wires :class:`flopy.mf6.ModflowGwfgwt` so the GWT (solute transport) model
    consumes the flow field from its paired GWF model. ``models`` is the resolved
    ``(gwf, gwt)`` pair. Usually emitted from an
    :class:`~myflopy.specs.ExchangeSpec`; see ``mf.gwt(...)`` for the model side.

    Parameters
    ----------
    simulation
        The owning FloPy simulation.
    models
        The flow/transport pair, as ``(gwf_model, gwt_model)``.
    **options
        Exchange options forwarded to the FloPy constructor.
    """

    return _build_two_model_exchange(flopy.mf6.ModflowGwfgwt, simulation, models, **options)


def build_gwf_gwe_exchange(simulation, models: tuple[Any, ...], **options):
    """Build a GWF-GWE exchange coupling a flow model to an energy model.

    Wires :class:`flopy.mf6.ModflowGwfgwe` so the GWE (heat/energy transport)
    model consumes the flow field from its paired GWF model. ``models`` is the
    resolved ``(gwf, gwe)`` pair. Usually emitted from an
    :class:`~myflopy.specs.ExchangeSpec`; see ``mf.gwe(...)`` for the model side.

    Parameters
    ----------
    simulation
        The owning FloPy simulation.
    models
        The flow/energy pair, as ``(gwf_model, gwe_model)``.
    **options
        Exchange options forwarded to the FloPy constructor.
    """

    return _build_two_model_exchange(flopy.mf6.ModflowGwfgwe, simulation, models, **options)


def build_gwf_prt_exchange(simulation, models: tuple[Any, ...], **options):
    """Build a GWF-PRT exchange coupling a flow model to a particle-tracking model.

    Wires :class:`flopy.mf6.ModflowGwfprt` so the PRT (particle tracking) model
    advects particles through its paired GWF model's flow field. ``models`` is
    the resolved ``(gwf, prt)`` pair. Usually emitted from an
    :class:`~myflopy.specs.ExchangeSpec`; see ``mf.prt(...)`` for the model side.

    Parameters
    ----------
    simulation
        The owning FloPy simulation.
    models
        The flow/particle pair, as ``(gwf_model, prt_model)``.
    **options
        Exchange options forwarded to the FloPy constructor.
    """

    return _build_two_model_exchange(flopy.mf6.ModflowGwfprt, simulation, models, **options)
