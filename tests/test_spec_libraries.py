"""The project's three spec libraries (design step 3).

``project.packages``/``grids``/``simulations`` were plain dicts, so an editor
could offer ``.keys()`` and nothing else and there was no way to see what a
project held without already knowing the keys. They are now ``dict``
**subclasses**, which is the whole bet: every existing dict operation must keep
working unchanged, or this is a breaking change wearing a feature's clothes.
"""

from __future__ import annotations

import inspect
import typing
from pathlib import Path

import numpy as np
import pytest

import myflopy as mf
from myflopy.library import (
    GridLibrary,
    PackageLibrary,
    SimulationLibrary,
    SpecLibrary,
)
from myflopy.workspace import Project


def _project(tmp_path: Path) -> Project:
    """A project carrying one of everything, keyed the way the user keys things."""

    project = mf.Project(tmp_path / "tentrails", name="tentrails")
    project.add_package("npf/base", mf.npf(k=1.0))
    project.add_package("rch_pond", mf.rch.array(recharge={0: np.zeros(4)}, irch=None))
    project.add_grid("grid/base", mf.GridSpec.python("g.py", function="build"))
    project.add_simulation(
        mf.SimulationSpec(
            "baseline", models=[mf.gwf("flow", packages=[mf.ref("npf/base")])]
        )
    )
    return project


# --- the bet: nothing a dict did stops working ------------------------------


def test_a_library_is_still_a_dict_in_every_way_that_matters(tmp_path):
    """If any of this fails, the subclass is a breaking change, not a feature."""

    project = _project(tmp_path)
    packages = project.packages

    assert isinstance(packages, dict)
    assert dict(packages) == {k: packages[k] for k in packages}
    assert "npf/base" in packages
    assert len(packages) == 2
    assert sorted(packages) == ["npf/base", "rch_pond"]
    assert packages.get("nope") is None            # never routes through __missing__
    assert packages.get("nope", "fallback") == "fallback"
    assert [key for key, _ in packages.items()] == list(packages)

    packages["extra"] = mf.ic(strt=1.0)            # plain assignment still registers
    assert project.package("extra").name == "ic"
    packages.pop("extra")
    assert "extra" not in packages


def test_copy_and_merge_stay_libraries(tmp_path):
    """``dict.copy`` and ``dict.__or__`` return a PLAIN dict -- measured -- which
    would silently strip the verbs. Both are overridden."""

    packages = _project(tmp_path).packages

    assert type(packages.copy()) is PackageLibrary
    assert type(packages | {}) is PackageLibrary
    assert type({} | packages) is PackageLibrary
    assert packages.copy()._project_name == "tentrails"
    assert (packages | {"x": mf.ic(strt=1.0)})._project_name == "tentrails"


def test_project_name_is_not_a_constructor_argument():
    """It is a class-level default on purpose: a required ``__init__`` argument on
    a dict subclass is dropped by ``copy()``/``|`` and by ``dict``-based revival."""

    # No custom __init__ at all -- it inherits dict's, so there is no argument
    # to lose and `PackageLibrary(some_dict)` revives cleanly.
    assert PackageLibrary.__init__ is dict.__init__
    assert SpecLibrary._project_name is None
    assert PackageLibrary({"a": mf.ic(strt=1.0)})._project_name is None


# --- discover / inspect ------------------------------------------------------


def test_a_missing_key_names_the_near_miss_and_everything_defined(tmp_path):
    project = _project(tmp_path)

    with pytest.raises(KeyError) as caught:
        project.packages["npf/bass"]

    message = caught.value.args[0]
    assert "Project 'tentrails' has no package 'npf/bass'." in message
    assert "Did you mean: 'npf/base'?" in message
    assert "npf/base, rch_pond" in message

    # The accessor raises the same thing -- one message, two spellings.
    with pytest.raises(KeyError, match="npf/base"):
        project.package("npf/bass")


def test_an_empty_library_says_so_rather_than_suggesting_nothing(tmp_path):
    project = mf.Project(tmp_path / "empty", name="empty")

    with pytest.raises(KeyError, match="No packages are defined."):
        project.packages["anything"]

    # ...and its summary still has columns, so `.empty` is the test.
    frame = project.packages.summary()
    assert frame.empty
    assert list(frame.columns) == list(PackageLibrary._columns)


def test_summaries_describe_each_tier(tmp_path):
    project = _project(tmp_path)

    packages = project.packages.summary().set_index("key")
    assert packages.loc["npf/base", "package"] == "npf"
    # A list-BC/array helper builds through functools.partial; the label unwraps
    # it to the FloPy class rather than saying "partial".
    assert packages.loc["rch_pond", "builder"] == "flopy.ModflowGwfrcha"

    grids = project.grids.summary().set_index("key")
    assert grids.loc["grid/base", "method"] == "python"
    assert bool(grids.loc["grid/base", "built"]) is False
    assert grids.loc["grid/base", "cells"] == "—"      # never resolves a recipe

    simulations = project.simulations.summary().set_index("name")
    assert simulations.loc["baseline", "models"] == "flow"


def test_persistence_is_a_method_not_a_summary_column(tmp_path):
    """Classifying persistence means attempting ``to_dict()`` -- 119 ms on a
    100k-record list BC, and Jupyter calls ``_repr_html_`` on every bare project.
    So it is opt-in, and ``summary()`` must not carry it."""

    project = _project(tmp_path)

    assert "persists" not in project.packages.summary().columns

    report = project.packages.persistence().set_index("key")
    assert report.loc["npf/base", "persists"] == "json"
    assert report.loc["rch_pond", "persists"] == "pickle"
    assert "ndarray" in report.loc["rch_pond", "reason"]


def test_project_repr_touches_no_disk(tmp_path, monkeypatch):
    """PyCharm calls ``__repr__`` on every debugger step, and ``Run.reopen`` does
    an unguarded ``json.loads(...)["run"]`` -- so one truncated manifest must not
    be able to make ``repr(project)`` raise."""

    project = _project(tmp_path)

    def explode(*args, **kwargs):
        raise AssertionError("repr must not scan the filesystem")

    monkeypatch.setattr(Path, "glob", explode)
    monkeypatch.setattr(Path, "rglob", explode)
    monkeypatch.setattr(Path, "iterdir", explode)

    assert repr(project) == (
        f"Project('tentrails', root={str(project.root)!r}, "
        "packages=2, grids=1, simulations=1)"
    )
    assert "tentrails" in project._repr_html_()
    assert "runs" not in repr(project)          # runs live on disk; not counted here


def test_project_summary_covers_all_three_libraries(tmp_path):
    frame = _project(tmp_path).summary()

    assert list(frame.columns) == ["library", "key", "detail"]
    assert set(frame["library"]) == {"package", "grid", "simulation"}
    assert set(frame["key"]) == {"npf/base", "rch_pond", "grid/base", "baseline"}


# --- the libraries survive persistence --------------------------------------


def test_libraries_survive_save_and_load(tmp_path):
    project = _project(tmp_path)
    project.save()

    loaded = mf.Project.load(project.root)

    assert type(loaded.packages) is PackageLibrary
    assert type(loaded.grids) is GridLibrary
    assert type(loaded.simulations) is SimulationLibrary
    assert loaded.packages.names() == ["npf/base", "rch_pond"]
    assert loaded.package("npf/base").options["k"] == 1.0
    # the reconstructed libraries know their project, so misses read the same
    with pytest.raises(KeyError, match="Project 'tentrails'"):
        loaded.packages["nope"]


# --- static discoverability ratchet -----------------------------------------


def test_the_annotations_an_editor_reads_stay_concrete():
    """PyCharm and Pylance never run the module -- only the annotation reaches
    them. ``Any`` here would undo the whole point (ledger 171)."""

    source = inspect.getsource(Project.__init__)
    assert "self.packages: PackageLibrary" in source
    assert "self.grids: GridLibrary" in source
    assert "self.simulations: SimulationLibrary" in source

    # The element type is fixed at the class statement, so `library[key]` infers.
    assert PackageLibrary.__orig_bases__[0] == SpecLibrary[mf.PackageSpec]
    assert SimulationLibrary.__orig_bases__[0] == SpecLibrary[mf.SimulationSpec]

    accessors = {
        "package": mf.PackageSpec,
        "simulation": mf.SimulationSpec,
    }
    for name, expected in accessors.items():
        assert typing.get_type_hints(getattr(Project, name))["return"] is expected


def test_every_library_answers_the_same_verbs():
    """One grammar over the tiers: adding a verb to one and forgetting the others
    is the asymmetry that produced this work in the first place."""

    verbs = {"names", "summary", "copy", "_row", "_repr_html_"}
    for library in (PackageLibrary, GridLibrary, SimulationLibrary):
        assert verbs <= set(dir(library)), library.__name__
        assert library._kind != SpecLibrary._kind
        assert library._columns != SpecLibrary._columns
