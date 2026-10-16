"""`with_grid(remap=True)` carries a library-referenced npf onto a new mesh (ledger 610).

Ledger 597 closed the silent half of the defect: a move without `library=` refuses a
model holding an `mf.ref`, and with it the ref is classified as the package it names.
Measured at 55947911 on the route a calibration's write-back takes -- an `npf` built
by `substitute_values`, registered and `mf.ref`'d beside an EAGER `mf.disv(<result>)`,
which is the shape of every Ten Trails run -- four smaller defects remained:

* `remap=True` asked the overlay engine to average the DISV's own `top`/`botm`/
  `idomain` and refused the whole move (`RemapRefused`), inline or through a ref;
* a package named in `grid_free=` was remapped anyway, and the report asked with the
  same `grid_free=` then called it `unchanged`;
* a stamped array moved to a mesh of ANOTHER size died inside FloPy
  (`MFDataException: Unable to set data "k" layer 0`), never naming the move;
* a model a SAVE stripped of its mesh (ledger 257) was taken for a model with none:
  the simulation tier skipped it and returned it unmoved, and the model tier moved it
  unstamped, so it BUILT the old mesh's K on a new mesh of the same size;

and the new mesh could not be read through the library: a `grid_ref` target raised a
`TypeError` naming `project.grid(key)`, which then raised `FileNotFoundError` anywhere
but the project root, because a reloaded library grid is a pickle with a relative path.
"""

from __future__ import annotations

import os
import warnings

import numpy as np
import pytest
from flopy.mf6.mfbase import MFDataException

import myflopy as mf
from myflopy.grid_identity import (
    INLINED_FROM_REF,
    MOVED_FROM,
    REMAPPABLE_ARRAYS,
    REMAPPED_FROM,
    grid_pin,
    moved_package_complaints,
)
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.modflow.mf6.grid.remap import CONDUCTIVITY_ARRAYS, INTENSIVE_ARRAYS
from myflopy.specs import GridRef, ModelContext, PackageRef, substitute_values
from myflopy.workspace import Project

NLAY = 2


def _rect(nrow, ncol):
    return rectangular_voronoi(CanonicalModelConfig(nrow=nrow, ncol=ncol, nlay=NLAY, nper=1))


@pytest.fixture(scope="module")
def meshes():
    """Two meshes of 28 cells (the shape no range check sees) and one of 81."""

    a, b, c = _rect(4, 7), _rect(7, 4), _rect(9, 9)
    assert a.ncpl == b.ncpl == 28 and c.ncpl == 81
    return a, b, c


def _layering(grid):
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("upper", thickness=20.0)
    stack.add("lower", thickness=30.0)
    return stack.build(grid)


def _k(grid):
    return np.linspace(1.0, 50.0, NLAY * grid.ncpl).reshape(NLAY, grid.ncpl)


def _calibrated(k):
    """An npf as `CalibratedInputs.packages()` registers one: `substitute_values` over a base."""

    return substitute_values(
        mf.npf(k=np.ones_like(k), icelltype=1),
        {"k": k},
        provenance={"pest": "cal", "realization": "base", "source_key": "npf/base", "run": "ss"},
        metadata={"source_path": "/x/cal.par", "source_sha": "deadbeef"},
    )


def _sim(model):
    return mf.simulation(
        model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims()]
    )


def _model(grid, *packages, **gwf):
    context = ModelContext(grid=grid, domain=np.ones((NLAY, grid.ncpl), dtype=int))
    return mf.gwf(
        "m",
        context=context,
        packages=[*packages, mf.ic(strt=90.0), mf.oc(head_filerecord="m.hds")],
        **gwf,
    )


def _rows(spec, library=None, **kwargs):
    return {row["package"]: row for row in spec.move_report(library=library, **kwargs)}


def _built_k(run):
    return run.built.built_model("m").model.get_package("npf").k.array


@pytest.fixture
def quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


# --- (a) a discretization is never remapped -------------------------------------- #


@pytest.mark.parametrize("form", ["inline", "ref"])
def test_remap_leaves_an_eager_discretization_for_with_packages(tmp_path, meshes, quiet, form):
    """`RemapRefused: disv: 'botm', 'idomain', 'top' ...` refused the whole move before.

    The DISV IS the old mesh; its arrays are not data on it. It reads `stale`, and the
    caller swaps in the new mesh's discretization.
    """

    a, b, _c = meshes
    project = Project(tmp_path / "p", name="p")
    dis = mf.disv(_layering(a))
    project.add_package("dis/old", dis)
    entry = mf.ref("dis/old") if form == "ref" else dis
    model = _model(a, entry, mf.npf(k=_k(a), icelltype=1))

    moved = model.with_grid(b, remap=True, library=project)

    rows = _rows(moved, project)
    dis_label = "dis/old" if form == "ref" else "disv"
    assert rows[dis_label]["how"] == "stale"
    assert rows["npf"]["how"] == "remapped"
    kept = next(p for p in moved.packages if (getattr(p, "key", None) or p.name) == dis_label)
    if form == "ref":
        assert isinstance(kept, PackageRef), "a discretization a remap leaves alone stays a ref"
    else:
        assert REMAPPED_FROM not in (kept.metadata or {})
        assert np.array_equal(kept.options["top"], dis.options["top"])


# --- the route end to end: a calibrated library npf onto another mesh -------------- #


@pytest.mark.parametrize("target", [1, 2], ids=["same-size", "81-cells"])
def test_a_calibrated_ref_k_moves_to_a_mesh_of_the_same_and_another_size(
    tmp_path, meshes, quiet, target
):
    a, new = meshes[0], meshes[target]
    k_old = _k(a)
    project = Project(tmp_path / "p", name="p")
    project.add_package("npf/base/cal", _calibrated(k_old))
    project.add_package("dis/old", mf.disv(_layering(a)))
    project.add_package("dis/new", mf.disv(_layering(new)))
    by_ref = _model(a, mf.ref("dis/old"), mf.ref("npf/base/cal"))
    inline = _model(a, mf.disv(_layering(a)), _calibrated(k_old))

    moved = by_ref.with_grid(new, remap=True, library=project)
    twin = inline.with_grid(new, remap=True)

    row = _rows(moved, project)["npf"]
    assert row["how"] == "remapped"
    npf = next(p for p in moved.packages if getattr(p, "name", None) == "npf")
    assert npf.metadata[INLINED_FROM_REF] == "npf/base/cal"
    swapped = moved.with_packages(mf.ref("dis/new"), slot="dis/old")
    run = project.prepare_run("moved", _sim(swapped))
    k = _built_k(run)
    assert k.shape == (NLAY, new.ncpl)
    assert not (k.shape == k_old.shape and np.array_equal(k, k_old)), "the OLD mesh's K built"
    twin_npf = next(p for p in twin.packages if getattr(p, "name", None) == "npf")
    np.testing.assert_allclose(k, twin_npf.options["k"])
    assert np.array_equal(project.package("npf/base/cal").options["k"], k_old)


# --- (a') grid_free= is never remapped --------------------------------------------- #


def test_grid_free_is_never_remapped(meshes, quiet):
    """Measured before: carried, and the report asked with the same `grid_free=` said `unchanged`."""

    a, b, _c = meshes
    k = _k(a)
    model = _model(a, mf.disv(_layering(a)), mf.npf(k=k, icelltype=1))

    moved = model.with_grid(b, remap=True, grid_free=("npf",))

    npf = next(p for p in moved.packages if getattr(p, "name", None) == "npf")
    assert np.array_equal(npf.options["k"], k)
    assert REMAPPED_FROM not in (npf.metadata or {})
    assert MOVED_FROM not in (npf.metadata or {})


# --- (b) FloPy's refusal of a stamped array names the move ------------------------- #


@pytest.mark.parametrize("form", ["inline", "ref"])
def test_a_stamped_array_another_mesh_cannot_hold_is_refused_naming_the_move(
    tmp_path, meshes, quiet, form
):
    a, _b, c = meshes
    project = Project(tmp_path / "p", name="p")
    project.add_package("npf/k", mf.npf(k=_k(a), icelltype=1))
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("upper", thickness=20.0)
    stack.add("lower", thickness=30.0)
    npf = mf.ref("npf/k") if form == "ref" else mf.npf(k=_k(a), icelltype=1)
    model = _model(a, mf.disv(stack), npf)

    moved = model.with_grid(c, library=project)

    spelling = (
        "remap=True, library=project" if form == "ref" else r"with_grid\(\.\.\., remap=True\)"
    )
    with pytest.raises(ValueError, match="was moved to this grid by `with_grid`") as raised:
        project.prepare_run("r", _sim(moved))
    assert __import__("re").search(spelling, str(raised.value)), str(raised.value)
    assert isinstance(raised.value.__cause__, MFDataException)


def test_the_move_refusal_names_remap_only_for_arrays_it_carries(meshes):
    a, _b, _c = meshes
    k_array = mf.npf(k=_k(a), icelltype=1).with_metadata(**{MOVED_FROM: "x"})
    well = mf.wel(stress_period_data={0: [[(0, 3), -50.0]]}).with_metadata(**{MOVED_FROM: "x"})
    recharge = mf.rch.array(recharge=np.full(a.ncpl, 1e-3), irch=0).with_metadata(
        **{MOVED_FROM: "x"}
    )

    (array_says,) = moved_package_complaints([k_array])
    (well_says,) = moved_package_complaints([well])
    (recharge_says,) = moved_package_complaints([recharge])

    assert "remap=True" in array_says and "(k)" in array_says
    assert "remap=True" not in well_says, "a typed cellid has no value to weight"
    assert "remap=True" not in recharge_says, "the engine refuses an extensive array"


def test_the_refusal_states_what_the_engine_carries():
    """`REMAPPABLE_ARRAYS` is restated in a lower import layer; the two must not drift."""

    assert REMAPPABLE_ARRAYS == INTENSIVE_ARRAYS | CONDUCTIVITY_ARRAYS


# --- (c) a model a save stripped of its mesh -------------------------------------- #


def _saved_and_reloaded(tmp_path, a, *, gridded=True):
    """A project holding a library npf (an array) on mesh ``a``, saved and loaded."""

    root = tmp_path / "saved"
    project = Project(root, name="p")
    project.add_package("npf/cal", mf.npf(k=_k(a), icelltype=1))
    project.add_package("dis/old", mf.disv(_layering(a)))
    if gridded:
        model = _model(a, mf.ref("dis/old"), mf.ref("npf/cal"))
    else:
        model = mf.gwf("m", packages=[mf.ref("dis/old"), mf.ref("npf/cal"), mf.ic(strt=90.0)])
    project.add_simulation(_sim(model))
    report = project.save()
    assert report.simulations == 1, report
    return Project.load(root)


_STRIPPED_CALLS = {
    "simulation with_grid": lambda sim, b, q: sim.with_grid(b, library=q),
    "simulation with_grid remap": lambda sim, b, q: sim.with_grid(b, remap=True, library=q),
    "model with_grid": lambda sim, b, q: sim.model("m").with_grid(b, library=q),
    "model with_grid remap": lambda sim, b, q: sim.model("m").with_grid(b, remap=True, library=q),
    "model move_report": lambda sim, b, q: sim.model("m").move_report(library=q),
    "simulation move_report": lambda sim, b, q: sim.move_report(library=q),
}


@pytest.mark.parametrize("call", sorted(_STRIPPED_CALLS))
def test_a_model_a_save_stripped_of_its_mesh_is_refused_not_skipped(tmp_path, meshes, quiet, call):
    """Before: the simulation tier returned it unmoved, the model tier built the old K."""

    a, b, _c = meshes
    q = _saved_and_reloaded(tmp_path, a)
    sim = q.simulations["s"]
    assert sim.model("m").context.grid is None

    with pytest.raises(ValueError, match=rf"declared on mesh {grid_pin(a)}.*mf\.grid_ref"):
        _STRIPPED_CALLS[call](sim, b, q)


def test_with_grid_gives_a_reloaded_model_back_its_own_mesh(tmp_path, meshes, quiet):
    a, b, _c = meshes
    q = _saved_and_reloaded(tmp_path, a)
    q.add_grid("old", mf.GridSpec.from_object(a))
    q.add_grid("new", mf.GridSpec.from_object(b))
    stripped = q.simulations["s"].model("m")

    with pytest.raises(ValueError, match=rf"declared on mesh {grid_pin(a)}"):
        stripped.with_grid(mf.grid_ref("new"), library=q)
    with pytest.raises(ValueError, match="nothing to carry"):
        stripped.with_grid(mf.grid_ref("old"), remap=True, library=q)

    attached = stripped.with_grid(mf.grid_ref("old"), library=q)

    assert all(e.moved_from is None for e in attached.packages if isinstance(e, PackageRef))
    assert "grid_pin" not in (attached.context.metadata or {}), "the old record went with it"
    assert _rows(attached, q)["npf/cal"]["how"] == "unverifiable", "measured on its mesh again"
    moved = attached.with_grid(mf.grid_ref("new"), remap=True, library=q)
    assert _rows(moved, q)["npf"]["how"] == "remapped"


def test_a_never_gridded_model_takes_its_first_grid_after_a_reload(tmp_path, meshes, quiet):
    """The guard reads the RECORD, so a model that never had a mesh is not over-refused."""

    a, b, _c = meshes
    q = _saved_and_reloaded(tmp_path, a, gridded=False)
    model = q.simulations["s"].model("m")
    assert "grid_pin" not in (model.context.metadata or {})

    moved = q.simulations["s"].with_grid(b)

    assert moved.lineage[-1]["models"] == [], "a gridless model is still left alone"
    assert model.with_grid(b).grid is not None


# --- (d) the new grid is read through library= ------------------------------------ #


@pytest.fixture(params=["rectangular", "triangulated"])
def two_meshes(request, meshes):
    """Two meshes of one CRS: rectangular_voronoi's, and real `GridSpec.voronoi` grids,
    whose `cell2d` mixes numpy.int64 vertex indices into what is pickled (ledgers 263/265)."""

    if request.param == "rectangular":
        return meshes[0], meshes[1]
    return (
        request.getfixturevalue("triangulated_grid"),
        request.getfixturevalue("triangulated_grid_finer"),
    )


def test_remap_reads_the_new_grid_through_the_library_after_a_reload(
    tmp_path, two_meshes, quiet, monkeypatch
):
    """`project.grid(key)` raised FileNotFoundError off the project root; a grid_ref, TypeError."""

    a, b = two_meshes
    root = tmp_path / "p"
    project = Project(root, name="p")
    project.add_grid("old", mf.GridSpec.from_object(a))
    project.add_grid("new", mf.GridSpec.from_object(b))
    project.add_package("npf/cal", mf.npf(k=_k(a), icelltype=1))
    project.add_package("dis/old", mf.disv(_layering(a)))
    project.add_package("dis/new", mf.disv(_layering(b)))
    project.add_simulation(
        _sim(
            mf.gwf(
                "m",
                grid=mf.grid_ref("old"),
                packages=[mf.ref("dis/old"), mf.ref("npf/cal"), mf.ic(strt=90.0)],
            )
        )
    )
    project.save()
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    q = Project.load(root)
    model = q.simulations["s"].model("m")
    assert q.grid("new").method == "pickle"

    by_spec = model.with_grid(q.grid("new"), remap=True, library=q)
    by_ref = model.with_grid(mf.grid_ref("new"), remap=True, library=q)

    assert isinstance(by_ref.grid, GridRef) and by_ref.grid.key == "new"
    npf_spec = next(p for p in by_spec.packages if getattr(p, "name", None) == "npf")
    npf_ref = next(p for p in by_ref.packages if getattr(p, "name", None) == "npf")
    np.testing.assert_array_equal(npf_ref.options["k"], npf_spec.options["k"])
    assert npf_ref.options["k"].shape == (NLAY, b.ncpl)

    q.add_package("npf/cal/new", npf_ref)
    on_new = by_ref.with_packages(mf.ref("npf/cal/new"), slot="npf").with_packages(
        mf.ref("dis/new"), slot="dis/old"
    )
    q.add_simulation(_sim(on_new).derive("s_new"))
    assert q.save().skipped == ()
    again = Project.load(root)
    run = again.prepare_run("r", again.simulation("s_new"))
    np.testing.assert_allclose(_built_k(run), npf_ref.options["k"])
    assert (run.workspace / ".myflopy" / "frozen.npz").is_file()
    assert os.listdir(elsewhere) == [], "a move writes no mesh into the working directory"


def test_a_grid_ref_target_without_the_library_is_refused_naming_it(meshes):
    a, _b, _c = meshes
    model = _model(a, mf.disv(_layering(a)), mf.npf(k=_k(a), icelltype=1))
    with pytest.raises(TypeError, match=r"library=project.*not GridRef"):
        model.with_grid(mf.grid_ref("new"), remap=True)
