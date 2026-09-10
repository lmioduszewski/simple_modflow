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
from myflopy.specs import fork_package
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
        "with_packages",
    ]
    # the parent stays registered, or Project.save would refuse the fork
    assert "baseline" in project.simulations
    project.save()


def test_a_renamed_fork_does_not_collapse_onto_its_parent(tmp_path):
    """A package name is not just ``PackageSpec.name``. The helpers bake it into
    ``options`` as ``pname``/``filename``, so a fork that renamed only the field
    left those stale -- and MF6 kept ONE package of two under a FloPy "Replacing
    existing package" warning. Measured before the fix: two CHD specs on one
    model wrote a single ``flow.chd`` holding only the second.
    """

    base = mf.chd(stress_period_data={0: [[(0, 0, 0), 10.0]]})
    forked = fork_package(
        base, name="chd_east", stress_period_data={0: [[(0, 0, 2), 8.0]]}
    )

    assert forked.options["pname"] == "chd_east"
    assert forked.options["filename"] == "{model_name}.chd_east"

    simulation = mf.SimulationSpec(
        "s",
        models=[
            mf.gwf(
                "flow",
                packages=[
                    mf.dis(nlay=1, nrow=1, ncol=3, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                    mf.ic(strt=5.0),
                    mf.npf(k=1.0),
                    base,
                    forked,
                ],
            )
        ],
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["flow"])],
    )
    built = simulation.build_flopy(tmp_path / "ws")
    built.simulation.write_simulation(silent=True)

    names = [p.package_name for p in built.simulation.get_model("flow").packagelist]
    assert "chd" in names and "chd_east" in names
    written = sorted(f.name for f in (tmp_path / "ws").iterdir() if ".chd" in f.name)
    assert written == ["flow.chd", "flow.chd_east"]


def test_a_rename_that_cannot_be_carried_safely_is_refused(tmp_path):
    """uzf/lak/sfr write ``{model_name}_{name}_<suffix>`` while mvr writes
    ``{model_name}.{name}.bud``, and ``PackageSpec`` has no package-type field to
    dispatch on -- so a filerecord still naming the old package raises rather
    than being guessed at. A ``str.replace`` would corrupt
    ``{model_name}_sfr_stage.sfr`` when renaming a package called ``sfr``."""

    spec = mf.sfr.flopy(
        nreaches=1,
        packagedata=[],
        connectiondata=[],
        perioddata={},
        pname="sfr",
        budget_filerecord="{model_name}_sfr_budget.sfr",
    )

    with pytest.raises(ValueError) as caught:
        fork_package(spec, name="sfr_north")

    assert "budget_filerecord" in str(caught.value)
    assert "would be a guess" in str(caught.value)

    # ...and renaming is fine once the caller supplies the corrected values.
    # The check is keyed on whether the CALLER handled the key, not on whether
    # the value contains the old name: "{model_name}_sfr_north_budget.sfr" still
    # contains "sfr" twice, once because the new name embeds the old and once
    # because ".sfr" is the file suffix.
    renamed = fork_package(
        spec,
        name="sfr_north",
        budget_filerecord="{model_name}_sfr_north_budget.sfr",
        stage_filerecord="{model_name}_sfr_north_stage.sfr",
    )
    assert renamed.options["pname"] == "sfr_north"

    # Handling only some of them is still refused, naming the ones left.
    with pytest.raises(ValueError, match="stage_filerecord"):
        fork_package(spec, name="sfr_north", budget_filerecord="ok.bud")


def test_an_unrelated_pname_is_left_alone():
    """Only values still holding what the helper wrote are rewritten -- a pname
    the user set deliberately is theirs."""

    spec = mf.npf(k=1.0).with_options(pname="my_own_name")

    forked = fork_package(spec, name="npf_high")

    assert forked.name == "npf_high"
    assert forked.options["pname"] == "my_own_name"


# --- defects found by adversarial verification of step 4 --------------------


def _wells(tmp_path) -> Project:
    """Two wells keyed <type>/<variant> -- the project's own documented style,
    and the shape that made swap= replace the wrong one."""

    project = mf.Project(tmp_path / "wells", name="wells")
    # Distinct names, or _require_unique_names refuses the pair at build -- which
    # is the same guard commit 043f55c exists to keep working.
    project.add_package(
        "wel/north", mf.wel(stress_period_data={0: [((0, 0, 2), -100.0)]}, name="wel_north")
    )
    project.add_package(
        "wel/south", mf.wel(stress_period_data={0: [((0, 0, 7), -100.0)]}, name="wel_south")
    )
    project.add_package(
        "wel/south_off", mf.wel(stress_period_data={0: [((0, 0, 7), 0.0)]}, name="wel_south_off")
    )
    dis = dict(nlay=1, nrow=1, ncol=10, delr=10.0, delc=10.0, top=20.0, botm=0.0)
    project.add_simulation(
        mf.SimulationSpec(
            "baseline",
            models=[
                mf.gwf(
                    "flow",
                    packages=[
                        mf.dis(**dis),
                        mf.ic(strt=10.0),
                        mf.npf(k=10.0),
                        mf.ref("wel/north"),
                        mf.ref("wel/south"),
                    ],
                )
            ],
            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["flow"])],
        )
    )
    return project


def test_swap_replaces_the_named_package_not_its_slot_mate(tmp_path):
    """A slot is only the key's FIRST SEGMENT for a ref, so wel/north and
    wel/south both answered to "wel" and the swap landed on whichever came
    first. Measured before the fix: swapping wel/south deleted wel/north, left
    wel/south pumping, and MF6 reported normal termination on the wrong model."""

    project = _wells(tmp_path)

    forked = project.simulations.fork(
        "baseline", "south_off", swap={"wel/south": "wel/south_off"}
    )

    labels = [getattr(e, "key", None) or e.name for e in forked.model("flow").packages]
    assert labels == ["dis", "ic", "npf", "wel/north", "wel/south_off"]

    # ...and the built model carries the north well and the OFF south one.
    run = project.prepare_run("south_off", "south_off")
    built = sorted(
        n for n in run.flopy_model("flow").get_package_list() if n.startswith("WEL")
    )
    assert built == ["WEL_NORTH", "WEL_SOUTH_OFF"]


def test_the_model_qualified_swap_validates_its_labels_too(tmp_path):
    """The guarantee the flat form makes -- a label matching nothing raises --
    was absent from exactly the form the docs recommend for coupled models, so a
    typo registered a "variant" identical to its parent."""

    project = _coupled(tmp_path)

    with pytest.raises(KeyError) as caught:
        project.simulations.fork("coupled", "typo", swap={"tport": {"rch/nope": "ic/flow"}})

    assert "Model 'tport' of simulation 'coupled' has no package 'rch/nope'" in caught.value.args[0]
    assert "typo" not in project.simulations


def test_a_swap_naming_one_target_twice_is_refused(tmp_path):
    """Mixing the flat and qualified forms over the same entry made the result
    depend on the order the mapping happened to be written."""

    project = _coupled(tmp_path)

    with pytest.raises(ValueError, match="twice"):
        project.simulations.fork(
            "coupled", "v", swap={"ic": "ic/flow", "flow": {"ic": "ic/tport"}}
        )


def test_a_simulation_level_package_says_so_rather_than_denying_it_exists(tmp_path):
    """tdis/ims cannot be swapped, but the old message claimed the simulation had
    no such package while listing an inventory that omitted them."""

    project = _coupled(tmp_path)

    with pytest.raises(KeyError, match="simulation-level package"):
        project.simulations.fork("coupled", "v", swap={"tdis": "x"})


def test_a_non_string_swap_value_is_named(tmp_path):
    project = _coupled(tmp_path)

    with pytest.raises(TypeError, match="must be a library key"):
        project.simulations.fork("coupled", "v", swap={"ic": mf.ic(strt=1.0)})


def test_package_digests_are_cached_until_the_library_changes(tmp_path):
    """Digesting costs pickling PYTHON OBJECTS, not bytes: a 72x9,090-record EVT
    measures 1.3 s. ``RunLibrary.summary`` fingerprints the whole library per
    call and Jupyter calls it on every bare ``project.runs``, so without a cache
    the browse surface stalls on exactly the projects worth browsing."""

    project = _project(tmp_path)

    first = project.packages.digests()
    assert project.packages.digests() is first          # cached, not recomputed

    project.packages.edit("npf/base", k=50.0)
    assert project.packages._digest_cache is None       # dropped on mutation
    assert project.packages.digests()["npf/base"] != first["npf/base"]


@pytest.mark.parametrize(
    "mutate",
    [
        lambda lib: lib.__setitem__("new", mf.ic(strt=1.0)),
        lambda lib: lib.pop("npf/base"),
        lambda lib: lib.update({"another": mf.ic(strt=2.0)}),
        lambda lib: lib.setdefault("defaulted", mf.ic(strt=3.0)),
        lambda lib: lib.__delitem__("npf/base"),
        lambda lib: lib.clear(),
        # |= does NOT route through __or__ or update on a dict subclass, so it
        # slipped past every other hook and left a stale fingerprint behind.
        lambda lib: lib.__ior__({"merged": mf.ic(strt=4.0)}),
    ],
    ids=["setitem", "pop", "update", "setdefault", "delitem", "clear", "ior"],
)
def test_every_mutating_route_drops_the_digest_cache(tmp_path, mutate):
    """The cache is keyed on the library, so replacing an entry through ANY
    route must invalidate it -- a stale fingerprint reports a changed package
    as unchanged, which is the one thing drift must never do."""

    packages = _project(tmp_path).packages
    packages.digests()
    assert packages._digest_cache is not None

    mutate(packages)

    assert packages._digest_cache is None


def test_save_report_counts_pickled_grids_as_its_docstring_says(tmp_path):
    """``pickled`` documented "packages and grids" and counted only packages, so
    a project whose grid is a built mesh under-reported what it had cached."""

    from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("npf/base", mf.npf(k=1.0))          # JSON
    project.add_grid("grid/base", irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0))

    report = project.save()

    assert report.grids == 1
    assert report.pickled == 1        # the built mesh, which cannot be JSON


def test_prepare_run_reuses_the_library_digest_cache(tmp_path):
    """Digesting is expensive for record-list packages (1.3 s for the real EVT),
    so prepare_run must not re-pay it per run."""

    project = _project(tmp_path)
    project.add_simulation(
        mf.SimulationSpec(
            "s",
            models=[
                mf.gwf(
                    "flow",
                    packages=[
                        mf.dis(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                        mf.ref("npf/base"),
                    ],
                )
            ],
            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
        )
    )
    project.packages.digests()
    cached = project.packages._digest_cache

    project.prepare_run("r", "s", build=False)

    assert project.packages._digest_cache is cached   # reused, not recomputed
