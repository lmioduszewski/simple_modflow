"""Small builders for common FloPy 3.10 simulation-level objects."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any, Protocol

import flopy

from myflopy._length_units import mf6_length_option
from myflopy.specs import PackageSpec


class PackageBuilder(Protocol):
    """Contract for configured objects that produce reusable package specs.

    A package builder receives all configuration when it is created. Its
    ``build`` method only validates that configuration and returns the spec.

    Parameters
    ----------
    *args, **kwargs
        Accepted by nothing: this is a :class:`typing.Protocol`, a structural type
        for annotations, and calling ``PackageBuilder(...)`` raises ``TypeError:
        Protocols cannot be instantiated``. A class satisfies it by defining
        ``build(self) -> PackageSpec``; its own constructor takes the
        configuration.
    """

    def build(self) -> PackageSpec:
        """Return a package spec from the builder's stored configuration."""


def build_ims(simulation, *, models: Iterable[str] | None = None, **options):
    """Build an IMS solver package and register it with the named models.

    Each solution needs a distinct file. Without one, multiple IMS packages
    collapse into a single solution, which MF6 rejects for coupled simulations
    (for example GWF-GWT). The filename defaults to the package name. This is the
    engine under ``mf.ims(...)``.

    Parameters
    ----------
    simulation : flopy.mf6.MFSimulation
        The owning FloPy simulation to attach the solver to. The package is
        registered on it with ``simulation.register_ims_package(ims, models)``.
    models : iterable of str, optional
        Names of the models this IMS solves; each is registered to it. Required in
        practice: ``None`` (the default) raises, because it means the ``mf.ims(...)``
        declaration never reached an ``mf.simulation`` that fills in the models the
        other solvers left unclaimed. An empty iterable is accepted and registers
        the solver to no model.
    **options
        Forwarded to :class:`flopy.mf6.ModflowIms` (``complexity``,
        ``outer_maximum``, ``linear_acceleration``, ``pname``, ...); a name it does
        not take raises there. ``filename`` defaults to ``"<base>.ims"``, where
        ``<base>`` is ``pname`` when given, else the first model's name, else
        ``"ims"`` -- one file per solution, which is what keeps several IMS
        packages from collapsing into one.

    Returns
    -------
    flopy.mf6.ModflowIms
        The constructed and registered IMS package.

    Raises
    ------
    ValueError
        ``models`` is ``None``.
    """

    if models is None:
        raise ValueError(
            "mf.ims(...) was declared with no models= and never reached a\n"
            "simulation that could fill it in. Pass models=(...), or declare it\n"
            "inside mf.simulation(model, packages=[...]) which claims the rest."
        )

    models = list(models)
    base = options.get("pname") or (models[0] if models else "ims")
    options.setdefault("filename", f"{base}.ims")
    ims = flopy.mf6.ModflowIms(simulation, **options)
    simulation.register_ims_package(ims, models)
    return ims


def build_ems(simulation, *, models: Iterable[str] | None = None, **options):
    """Build an EMS (explicit) solver package and register it with the named models.

    PRT models are *explicit* in MF6 -- they must sit in an EMS6 solution group,
    never IMS6 ("Model X is an explicit model and cannot be added to an IMS6
    solution"). Registration is mandatory: an EMS merely *constructed* is dropped
    from ``mfsim.nam`` at write time and flopy silently attaches the PRT model to
    an existing IMS group, which MF6 then rejects. This is the engine under
    ``mf.ems(...)``.

    Parameters
    ----------
    simulation : flopy.mf6.MFSimulation
        The owning FloPy simulation to attach the solver to. The package is
        registered on it with ``simulation.register_solution_package(ems, models)``.
    models : iterable of str, optional
        Names of the (PRT) models this EMS solves; each is registered to it.
        ``None`` (the default) raises, for the same reason as :func:`build_ims`.
    **options
        Forwarded to :class:`flopy.mf6.ModflowEms`. ``filename`` defaults to
        ``"<base>.ems"``, where ``<base>`` is ``pname`` when given, else the first
        model's name, else ``"ems"``.

    Returns
    -------
    flopy.mf6.ModflowEms
        The constructed and registered EMS package.

    Raises
    ------
    ValueError
        ``models`` is ``None``.
    """

    if models is None:
        raise ValueError(
            "mf.ems(...) was declared with no models= and never reached a\n"
            "simulation that could fill it in. Pass models=(...), or declare it\n"
            "inside mf.simulation(model, packages=[...]) which claims the rest."
        )

    models = list(models)
    base = options.get("pname") or (models[0] if models else "ems")
    options.setdefault("filename", f"{base}.ems")
    ems = flopy.mf6.ModflowEms(simulation, **options)
    simulation.register_solution_package(ems, models)
    return ems


# Model-type dispatch for shared core packages (plan 5.3A, PRT rows 6.3A).
# ``ic``/``oc``/``dis(v/u)`` are the same package across model kinds but map to
# different FloPy classes; these module-level builders resolve the class from the
# built model's ``model_type`` (``"gwf6"``/``"gwt6"``/``"gwe6"``/``"prt6"``). They
# MUST stay module-level functions: ``PackageSpec`` serializes its builder by
# importable reference (``specs._callable_ref`` rejects lambdas/closures), and only
# the function is persisted -- the class tables below are re-read at build time,
# never serialized. PRT is dispatched for ``dis``/``disv``/``oc`` only: MF6 has no
# ``ModflowPrtic`` (PRT needs no initial condition) and no ``ModflowPrtdisu``.
# ``npf``/``sto`` are GWF-only in MF6 and are not dispatched.
_IC_CLASSES = {
    "gwf6": flopy.mf6.ModflowGwfic,
    "gwt6": flopy.mf6.ModflowGwtic,
    "gwe6": flopy.mf6.ModflowGweic,
}
_OC_CLASSES = {
    "gwf6": flopy.mf6.ModflowGwfoc,
    "gwt6": flopy.mf6.ModflowGwtoc,
    "gwe6": flopy.mf6.ModflowGweoc,
    "prt6": flopy.mf6.ModflowPrtoc,
}
_DISV_CLASSES = {
    "gwf6": flopy.mf6.ModflowGwfdisv,
    "gwt6": flopy.mf6.ModflowGwtdisv,
    "gwe6": flopy.mf6.ModflowGwedisv,
    "prt6": flopy.mf6.ModflowPrtdisv,
}
_DIS_CLASSES = {
    "gwf6": flopy.mf6.ModflowGwfdis,
    "gwt6": flopy.mf6.ModflowGwtdis,
    "gwe6": flopy.mf6.ModflowGwedis,
    "prt6": flopy.mf6.ModflowPrtdis,
}
_DISU_CLASSES = {
    "gwf6": flopy.mf6.ModflowGwfdisu,
    "gwt6": flopy.mf6.ModflowGwtdisu,
    "gwe6": flopy.mf6.ModflowGwedisu,
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
    """Build the IC package for the model's kind (GWF/GWT/GWE). Engine under ``mf.ic``.

    Parameters
    ----------
    model : flopy.mf6.MFModel
        The built FloPy model the package is added to. Its ``model_type``
        (``"gwf6"``, ...) picks the FloPy class.
    **options
        Forwarded to that class (``flopy.mf6.ModflowGwfic``, ``ModflowGwtic`` or ``ModflowGweic``).

    Raises
    ------
    ValueError
        ``model.model_type`` is not one of GWF/GWT/GWE.
    """

    return _dispatch_core_class(_IC_CLASSES, model, "ic")(model, **options)


def build_oc(model, **options):
    """Build the OC package for the model's kind (GWF/GWT/GWE/PRT). Engine under ``mf.oc``.

    Parameters
    ----------
    model : flopy.mf6.MFModel
        The built FloPy model the package is added to. Its ``model_type``
        (``"gwf6"``, ...) picks the FloPy class.
    **options
        Forwarded to that class (``flopy.mf6.ModflowGwfoc``, ``ModflowGwtoc``, ``ModflowGweoc`` or
        ``ModflowPrtoc``).

    Raises
    ------
    ValueError
        ``model.model_type`` is not one of GWF/GWT/GWE/PRT.
    """

    return _dispatch_core_class(_OC_CLASSES, model, "oc")(model, **options)


def _mf6_length_word(options: dict[str, Any]) -> dict[str, Any]:
    """``options`` with ``length_units`` as the word MF6 reads (ledger 432).

    Every route reaches the builder -- a frozen snapshot and a saved project
    included -- so this is where a spec holding ``"ft"`` (written ``FT``, read by
    MF6 as unit 0 before) is written ``FEET``; an unknown word raises. The spec
    itself is not changed, so its digest is not either.
    """

    if options.get("length_units") is None:
        return options
    return {**options, "length_units": mf6_length_option(options["length_units"])}


def build_disv(model, **options):
    """Build the DISV package for the model's kind (GWF/GWT/GWE/PRT). Engine under ``mf.disv``.

    Parameters
    ----------
    model : flopy.mf6.MFModel
        The built FloPy model the package is added to. Its ``model_type``
        (``"gwf6"``, ...) picks the FloPy class.
    **options
        Forwarded to that class (``flopy.mf6.ModflowGwfdisv``, ``ModflowGwtdisv``,
        ``ModflowGwedisv`` or ``ModflowPrtdisv``). ``length_units`` is written as MF6's word: one
        of MF6's four (``feet``/``meters``/``centimeters``/``unknown``, any
        case) as given, an abbreviation such as ``"ft"`` as ``FEET`` -- which
        MF6 read as undefined when written ``FT`` (ledger 432).

    Raises
    ------
    ValueError
        ``model.model_type`` is not one of GWF/GWT/GWE/PRT; or ``length_units``
        is a word that is neither MF6's nor a spelling it maps.
    """

    return _dispatch_core_class(_DISV_CLASSES, model, "disv")(model, **_mf6_length_word(options))


def build_dis(model, **options):
    """Build the DIS package for the model's kind (GWF/GWT/GWE/PRT). Engine under ``mf.dis``.

    Parameters
    ----------
    model : flopy.mf6.MFModel
        The built FloPy model the package is added to. Its ``model_type``
        (``"gwf6"``, ...) picks the FloPy class.
    **options
        Forwarded to that class (``flopy.mf6.ModflowGwfdis``, ``ModflowGwtdis``, ``ModflowGwedis``
        or ``ModflowPrtdis``). ``length_units`` is written as MF6's word: one
        of MF6's four (``feet``/``meters``/``centimeters``/``unknown``, any
        case) as given, an abbreviation such as ``"ft"`` as ``FEET`` -- which
        MF6 read as undefined when written ``FT`` (ledger 432).

    Raises
    ------
    ValueError
        ``model.model_type`` is not one of GWF/GWT/GWE/PRT; or ``length_units``
        is a word that is neither MF6's nor a spelling it maps.
    """

    return _dispatch_core_class(_DIS_CLASSES, model, "dis")(model, **_mf6_length_word(options))


def build_disu(model, **options):
    """Build the DISU package for the model's kind (GWF/GWT/GWE). Engine under ``mf.disu``.

    Parameters
    ----------
    model : flopy.mf6.MFModel
        The built FloPy model the package is added to. Its ``model_type``
        (``"gwf6"``, ...) picks the FloPy class.
    **options
        Forwarded to that class (``flopy.mf6.ModflowGwfdisu``, ``ModflowGwtdisu`` or
        ``ModflowGwedisu``). ``length_units`` is written as MF6's word: one
        of MF6's four (``feet``/``meters``/``centimeters``/``unknown``, any
        case) as given, an abbreviation such as ``"ft"`` as ``FEET`` -- which
        MF6 read as undefined when written ``FT`` (ledger 432).

    Raises
    ------
    ValueError
        ``model.model_type`` is not one of GWF/GWT/GWE; or ``length_units``
        is a word that is neither MF6's nor a spelling it maps.
    """

    return _dispatch_core_class(_DISU_CLASSES, model, "disu")(model, **_mf6_length_word(options))


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
    simulation : flopy.mf6.MFSimulation
        The owning FloPy simulation the exchange is added to.
    models : sequence of two FloPy models
        The two flow models to couple, as ``(model_a, model_b)``. Any length-2 sequence
        works; only each model's ``.name`` is read, as ``exgmnamea`` and
        ``exgmnameb`` in that order. Any other length raises ``ValueError``.
    **options
        Forwarded to :class:`flopy.mf6.ModflowGwfgwf` (``filename``, ``pname``, ...).
        ``exgmnamea``/``exgmnameb`` are set from ``models`` and so cannot be passed
        here: the call raises ``TypeError`` (got multiple values for keyword argument).

    Raises
    ------
    ValueError
        ``models`` does not hold exactly two models.
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
    simulation : flopy.mf6.MFSimulation
        The owning FloPy simulation the exchange is added to.
    models : sequence of two FloPy models
        The flow/transport pair, as ``(gwf_model, gwt_model)``. Any length-2 sequence
        works; only each model's ``.name`` is read, as ``exgmnamea`` and
        ``exgmnameb`` in that order. Any other length raises ``ValueError``.
    **options
        Forwarded to :class:`flopy.mf6.ModflowGwfgwt` (``filename``, ``pname``, ...).
        ``exgmnamea``/``exgmnameb`` are set from ``models`` and so cannot be passed
        here: the call raises ``TypeError`` (got multiple values for keyword argument).

    Raises
    ------
    ValueError
        ``models`` does not hold exactly two models.
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
    simulation : flopy.mf6.MFSimulation
        The owning FloPy simulation the exchange is added to.
    models : sequence of two FloPy models
        The flow/energy pair, as ``(gwf_model, gwe_model)``. Any length-2 sequence
        works; only each model's ``.name`` is read, as ``exgmnamea`` and
        ``exgmnameb`` in that order. Any other length raises ``ValueError``.
    **options
        Forwarded to :class:`flopy.mf6.ModflowGwfgwe` (``filename``, ``pname``, ...).
        ``exgmnamea``/``exgmnameb`` are set from ``models`` and so cannot be passed
        here: the call raises ``TypeError`` (got multiple values for keyword argument).

    Raises
    ------
    ValueError
        ``models`` does not hold exactly two models.
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
    simulation : flopy.mf6.MFSimulation
        The owning FloPy simulation the exchange is added to.
    models : sequence of two FloPy models
        The flow/particle pair, as ``(gwf_model, prt_model)``. Any length-2 sequence
        works; only each model's ``.name`` is read, as ``exgmnamea`` and
        ``exgmnameb`` in that order. Any other length raises ``ValueError``.
    **options
        Forwarded to :class:`flopy.mf6.ModflowGwfprt` (``filename``, ``pname``, ...).
        ``exgmnamea``/``exgmnameb`` are set from ``models`` and so cannot be passed
        here: the call raises ``TypeError`` (got multiple values for keyword argument).

    Raises
    ------
    ValueError
        ``models`` does not hold exactly two models.
    """

    return _build_two_model_exchange(flopy.mf6.ModflowGwfprt, simulation, models, **options)


# --- what each wrapper builder ultimately constructs ---------------------------
#
# A `PackageSpec` whose builder is a raw FloPy class can be inspected directly;
# a wrapper cannot, so it says here what it targets. `specs._accepted_options`
# reads this to reject a misspelled or misplaced option AT THE CALL rather than
# five frames into FloPy at build time -- which is how `model_nam_file=` reached
# `mf.ims` and surfaced only as `FlopyException: Extraneous kwargs`.
#
# A dispatching builder lists every kind it may construct; the accepted set is
# their union, because the model's kind is not known until build time.

build_ims.mf6_targets = (flopy.mf6.ModflowIms,)
build_ems.mf6_targets = (flopy.mf6.ModflowEms,)
build_ic.mf6_targets = tuple(_IC_CLASSES.values())
build_oc.mf6_targets = tuple(_OC_CLASSES.values())
build_disv.mf6_targets = tuple(_DISV_CLASSES.values())
build_dis.mf6_targets = tuple(_DIS_CLASSES.values())
build_disu.mf6_targets = tuple(_DISU_CLASSES.values())
