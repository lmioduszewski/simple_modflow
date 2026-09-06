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


# --- keys that name the same file ------------------------------------------


def test_two_keys_that_map_to_one_file_are_refused_before_writing(tmp_path):
    """``_slug`` rewrites every unsafe character to ``_``, so distinct keys can
    name one file. Writing both silently loses one -- measured before the fix:

        in memory : {'npf/high k': 111.0, 'npf/high_k': 999.0}
        after load: {'npf/high k': 999.0, 'npf/high_k': 999.0}

    Unlike a dangling reference, this cannot be written harmlessly and reported
    afterwards, so it is checked BEFORE anything is written.
    """

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("npf/high k", mf.npf(k=111.0))
    project.add_package("npf/high_k", mf.npf(k=999.0))

    with pytest.raises(ValueError) as caught:
        project.save()

    message = str(caught.value)
    assert "'npf/high k' and 'npf/high_k'" in message
    assert "specs/packages/npf/high_k.json" in message
    assert not (project.root / "specs").exists()      # nothing written


def test_grid_keys_collide_more_readily_than_package_keys(tmp_path):
    """A grid key is slugged FLAT -- every ``/`` collapses -- so ``base/coarse``
    and ``base_coarse`` name one file, where the package layout nests instead."""

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_grid("base/coarse", mf.GridSpec.python("g.py", function="b"))
    project.add_grid("base_coarse", mf.GridSpec.python("g.py", function="b"))

    with pytest.raises(ValueError, match="base_coarse"):
        project.save()

    # ...while the same shapes as PACKAGE keys nest and do not collide.
    ok = mf.Project(tmp_path / "ok", name="ok")
    ok.add_package("base/coarse", mf.ic(strt=1.0))
    ok.add_package("base_coarse", mf.ic(strt=2.0))
    ok.save()
    reloaded = mf.Project.load(ok.root)
    assert reloaded.package("base/coarse").options["strt"] == 1.0
    assert reloaded.package("base_coarse").options["strt"] == 2.0


# --- forking and editing (design step 4) ------------------------------------


def test_a_fork_is_built_not_patched(tmp_path):
    """The obvious implementation is silently wrong. ``PackageSpec`` is
    frozen+slots and validates in ``__post_init__``, but ``copy.deepcopy``
    reconstructs it through ``__reduce_ex__`` -- so deep-copying the spec and
    writing into ``copy.options`` accepts a typo that only surfaces at MF6 build.
    Forking goes through ``dataclasses.replace``, the ordinary path.
    """

    project = _project(tmp_path)

    with pytest.raises(TypeError, match="does not accept 'kk_typo'"):
        project.packages.fork("npf/base", "npf/typo", kk_typo=1.0)

    assert "npf/typo" not in project.packages       # nothing registered on failure


def test_a_fork_does_not_alias_the_options_it_did_not_touch(tmp_path):
    """Where the real hazard lives: ``with_options`` rebuilds the options dict but
    not its VALUES, so a "copy" can write through into the parent's arrays."""

    project = _project(tmp_path)
    project.add_package(
        "npf/arrays", mf.npf(k=np.array([1.0, 2.0]), k33=np.array([0.1, 0.2]))
    )

    forked = project.packages.fork("npf/arrays", "npf/high", k=np.array([9.0, 9.0]))
    forked.options["k33"][0] = -99.0                # an option the fork never set

    assert project.package("npf/arrays").options["k33"][0] == 0.1
    assert forked.options["k33"] is not project.package("npf/arrays").options["k33"]

    # ...and the sharing is available on purpose when asked for.
    shared = project.packages.fork("npf/arrays", "npf/shared", deep=False)
    assert shared.options["k33"] is project.package("npf/arrays").options["k33"]


def test_fork_registers_stamps_and_refuses_to_clobber(tmp_path):
    project = _project(tmp_path)

    forked = project.packages.fork(
        "rch_pond", "rch_pond_vol", name="rch_pond_vol", note="max-volume pond"
    )

    assert project.package("rch_pond_vol") is forked
    assert forked.name == "rch_pond_vol"
    stamp = forked.metadata["forked_from"]
    assert stamp["key"] == "rch_pond" and stamp["note"] == "max-volume pond"

    with pytest.raises(KeyError, match="already has a package 'rch_pond_vol'"):
        project.packages.fork("rch_pond", "rch_pond_vol")
    assert project.packages.fork("rch_pond", "rch_pond_vol", overwrite=True)


def test_edit_changes_what_a_key_means(tmp_path):
    """``fork`` leaves the original alone; ``edit`` changes what an existing key
    MEANS, which is the operation that can invalidate a run's history."""

    project = _project(tmp_path)

    edited = project.packages.edit("npf/base", k=50.0, note="post-calibration")

    assert project.package("npf/base") is edited
    assert edited.options["k"] == 50.0
    assert edited.metadata["edited"]["note"] == "post-calibration"
    assert "forked_from" not in edited.metadata     # a different stamp on purpose

    with pytest.raises(TypeError, match="does not accept"):
        project.packages.edit("npf/base", kk_typo=1.0)
    assert project.package("npf/base").options["k"] == 50.0   # unchanged by the failure


def _coupled(tmp_path) -> Project:
    """A GWF+GWT project: `ic` sits on BOTH models, which is the ambiguous case."""

    project = mf.Project(tmp_path / "coupled", name="coupled")
    project.add_package("ic/flow", mf.ic(strt=1.0))
    project.add_package("ic/tport", mf.ic(strt=0.0))
    dis = dict(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0)
    project.add_simulation(
        mf.SimulationSpec(
            "coupled",
            models=[
                mf.gwf("flow", packages=[mf.dis(**dis), mf.ic(strt=45.0)]),
                mf.gwt("tport", packages=[mf.dis(**dis), mf.ic(strt=0.0)]),
            ],
            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
        )
    )
    return project


def _labels(spec, model):
    return [getattr(e, "key", None) or e.name for e in spec.model(model).packages]


def test_a_flat_swap_reaches_every_model_carrying_the_label(tmp_path):
    project = _coupled(tmp_path)

    forked = project.simulations.fork("coupled", "both", swap={"ic": "ic/flow"})

    assert _labels(forked, "flow") == ["dis", "ic/flow"]
    assert _labels(forked, "tport") == ["dis", "ic/flow"]


def test_a_model_qualified_swap_reaches_exactly_one(tmp_path):
    """A coupled simulation carries ic/oc/dis on both models, so the flat form
    cannot express "just the transport one"."""

    project = _coupled(tmp_path)

    forked = project.simulations.fork(
        "coupled", "one", swap={"tport": {"ic": "ic/tport"}}
    )

    assert _labels(forked, "flow") == ["dis", "ic"]
    assert _labels(forked, "tport") == ["dis", "ic/tport"]


def test_a_swap_matching_nothing_raises_naming_the_source(tmp_path):
    """A typo that silently swapped nothing would hand back a "variant" identical
    to its parent. The message names the simulation the caller asked about, not
    the derived name they have not seen yet."""

    project = _coupled(tmp_path)

    with pytest.raises(KeyError) as caught:
        project.simulations.fork("coupled", "typo", swap={"rch/nope": "x"})

    assert "Simulation 'coupled' has no package 'rch/nope'" in caught.value.args[0]
    assert "typo" not in caught.value.args[0]
    assert "typo" not in project.simulations

    with pytest.raises(KeyError, match="has no model 'nope'"):
        project.simulations.fork("coupled", "bad", swap={"nope": {"ic": "ic/flow"}})


def test_a_simulation_fork_carries_its_lineage(tmp_path):
    project = _project(tmp_path)
    project.add_package("npf/high", mf.npf(k=50.0))

    forked = project.simulations.fork(
        "baseline", "high_k", swap={"npf/base": "npf/high"}
    )

    assert project.simulations["high_k"] is forked
    assert forked.derived_from == "baseline"
    assert [step["operation"] for step in forked.lineage] == [
        "derive",
        "replace_package",
    ]
    # the parent stays registered, or Project.save would refuse the fork
    assert "baseline" in project.simulations
    project.save()
