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
    """Build an IMS package and register it with the named models."""

    ims = flopy.mf6.ModflowIms(simulation, **options)
    simulation.register_ims_package(ims, list(models))
    return ims


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
    """Build a GWF-GWF exchange."""

    return _build_two_model_exchange(flopy.mf6.ModflowGwfgwf, simulation, models, **options)


def build_gwf_gwt_exchange(simulation, models: tuple[Any, ...], **options):
    """Build a GWF-GWT exchange."""

    return _build_two_model_exchange(flopy.mf6.ModflowGwfgwt, simulation, models, **options)


def build_gwf_gwe_exchange(simulation, models: tuple[Any, ...], **options):
    """Build a GWF-GWE exchange."""

    return _build_two_model_exchange(flopy.mf6.ModflowGwfgwe, simulation, models, **options)


def build_gwf_prt_exchange(simulation, models: tuple[Any, ...], **options):
    """Build a GWF-PRT exchange."""

    return _build_two_model_exchange(flopy.mf6.ModflowGwfprt, simulation, models, **options)
