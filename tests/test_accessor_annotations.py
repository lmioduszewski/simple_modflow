"""`run.model().fl.chd.` used to autocomplete LIST methods.

The ten FloPy package accessors -- on `SimulationBase` until ledger 339 moved them to
`model.fl` (`FlopyPackages`), `model.chd` being myflopy's view now -- carried no
return annotation, and
PyCharm and Pylance are STATIC -- they read the `def` line and never run the
module. With nothing there, `run.model().chd.` offered `append`, `count`, `sort`.

The annotation is a UNION on purpose. FloPy's `get_package` returns a list when a
model carries more than one package of a type -- two CHDs give a two-element
list, measured -- so promising a single package would be a lie an editor would
then propagate.
"""

from __future__ import annotations

import inspect
import typing

import pytest

from myflopy.modflow.mf6.simulation.base import FlopyPackages

ACCESSORS = {
    "chd": "ModflowGwfchd", "drn": "ModflowGwfdrn", "evt": "ModflowGwfevt",
    "ghb": "ModflowGwfghb", "lak": "ModflowGwflak", "rch": "ModflowGwfrch",
    "riv": "ModflowGwfriv", "sfr": "ModflowGwfsfr", "uzf": "ModflowGwfuzf",
    "wel": "ModflowGwfwel",
}


@pytest.mark.parametrize("verb, cls", sorted(ACCESSORS.items()))
def test_every_package_accessor_declares_what_it_returns(verb, cls):
    accessor = getattr(FlopyPackages, verb)
    annotation = inspect.signature(accessor.fget).return_annotation

    assert annotation is not inspect.Signature.empty, (
        f"{verb} has no return annotation, so an editor falls back to list methods"
    )
    assert cls in annotation


@pytest.mark.parametrize("verb", sorted(ACCESSORS))
def test_the_annotation_admits_the_list_flopy_can_return(verb):
    """A model with two CHDs really does hand back a list."""

    annotation = inspect.signature(getattr(FlopyPackages, verb).fget).return_annotation
    assert "list[" in annotation


def test_no_package_accessor_was_missed():
    """An accessor added later without an annotation fails here by name."""

    unannotated = [
        name
        for name in ACCESSORS
        if inspect.signature(getattr(FlopyPackages, name).fget).return_annotation
        is inspect.Signature.empty
    ]
    assert unannotated == []
