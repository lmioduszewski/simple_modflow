"""Numpy types at the serialization, hashing and snapshot boundaries.

Ledger 263 fixed `np.int64` at the `to_dict` boundary and drew the lesson "use a
real grid fixture". Measured afterwards (ledger 265), that was the smaller half of
the lesson: a real grid would not have found any of the six defects below, because
five of them need no grid at all and the sixth needs a `LayerBuildResult` rather
than a mesh. The durable guard is TYPE-LEVEL -- push the numpy types through each
boundary directly -- because it does not depend on which grid a future test picks.

Each test here names the symptom a user would actually have seen.
"""

from __future__ import annotations

import dataclasses
import json

import numpy as np
import pytest

import myflopy as mf
from myflopy.project.freeze import FrozenSnapshotError, freeze_simulation, thaw_simulation
from myflopy.specs import ModelContext, _canonical, _json_value, package_digest

SCALARS = [
    (np.int64(7), 7),
    (np.int32(7), 7),
    (np.float64(1.5), 1.5),
    (np.float32(1.5), 1.5),
    (np.array(3), 3),
    (np.array(1.5), 1.5),
]


@pytest.mark.parametrize("raw, expected", SCALARS)
def test_a_numpy_scalar_lands_as_an_exact_python_type(raw, expected):
    """`isinstance` is the wrong check here, and using it is what hid the bug.

    `np.float64` SUBCLASSES `float`, so it satisfied the primitive passthrough and
    was returned unconverted while `np.int64`, which subclasses nothing, raised.
    Converting only the half that raised left numpy floats in the output: json
    accepts them, PyYAML does not.
    """

    converted = _json_value(raw)
    assert converted == pytest.approx(expected)
    assert type(converted) is type(expected)


def test_a_numpy_string_lands_as_an_exact_str():
    """`np.str_` subclasses `str`, so it took the same free ride as `np.float64`."""

    assert type(_json_value(np.str_("HEAD"))) is str


def test_a_numpy_bool_lands_as_a_bool():
    """`np.bool_` does NOT subclass `bool`, which is the asymmetry in miniature."""

    assert _json_value(np.bool_(True)) is True


def test_a_scalar_that_cannot_unwrap_is_refused_rather_than_recursed():
    """`np.longdouble.item()` returns a `np.longdouble`, so recursing never ends.

    The user-visible difference is a `RecursionError` naming nothing against this
    module's own message naming the type.
    """

    with pytest.raises(ValueError, match="numpy.longdouble"):
        _json_value(np.longdouble(1.5))


def test_a_numpy_bool_in_enabled_still_serializes():
    """`enabled` was written raw, so it bypassed every conversion above.

    numpy 2 names its own type `bool`, so the failure read "Object of type bool is
    not JSON serializable" -- a message that sends you looking at the wrong type.
    """

    spec = dataclasses.replace(mf.npf(k=1.0), enabled=np.bool_(True))
    assert json.loads(json.dumps(spec.to_dict()))["enabled"] is True


def test_a_numpy_string_is_interned_like_any_other():
    """A pickle records identity, so an un-interned repeat reads as an edit.

    `np.str_` IS a `numpy.generic`, so it returned from the numpy branch before
    reaching the interning that ledger 222 added -- and `np.unique` on a boundname
    column is exactly how a GeoPackage produces one.
    """

    assert _canonical(np.str_("north_reach")) is _canonical("north_reach")


def test_a_zero_dimensional_array_canonicalizes_like_its_scalar():
    """`_json_value` normalized it and `_canonical` did not, so they disagreed.

    The disagreement IS the bug: the value the digest sees changes across a save,
    and `project.runs.drift()` reports a package nobody touched as edited.
    """

    assert type(_canonical(np.array(2.0))) is float
    assert _canonical(np.array(2.0)) == _canonical(np.float64(2.0))
    assert isinstance(_canonical(np.array([1.0, 2.0])), np.ndarray), "only 0-d unwraps"


def test_a_package_with_numpy_boundnames_digests_as_one_with_python_ones():
    """The end-to-end shape of the two tests above: a spurious `edited` verdict."""

    # The SAME name twice, deliberately. A pickle records identity, so the defect is a
    # missing back-reference and appears only when one object REPEATS -- with each name
    # used once this test passed with the fix fully reverted, certifying rather than
    # checking (ledger 266).
    names = np.unique(np.array(["north_reach", "north_reach", "south_reach"]))
    numpy_side = mf.chd(
        stress_period_data={0: [[(0, 0), 10.0, names[0]], [(0, 1), 10.0, names[0]]]},
        boundnames=True,
    )
    python_side = mf.chd(
        stress_period_data={0: [[(0, 0), 10.0, "north_reach"], [(0, 1), 10.0, "north_reach"]]},
        boundnames=True,
    )
    assert package_digest(numpy_side) == package_digest(python_side)


def test_a_numpy_valued_package_does_not_drift_across_a_save(tmp_path):
    """The symptom itself, not its mechanism: `runs.drift()` said EDITED.

    The two tests above pin the causes -- an un-interned `np.str_` and a 0-d array
    the two normalizers disagreed about. This is what the user saw: a package nobody
    touched, reported as changed after a plain `Project.save()` + `load()`, because
    "the key still exists and means something else".
    """

    from myflopy.workspace import Project

    names = np.unique(np.array(["north_reach", "north_reach", "south_reach"]))
    project = Project(tmp_path / "p", name="p")
    project.add_package(
        "chd_np",
        mf.chd(
            stress_period_data={
                0: [
                    [(0, 0), 10.0, names[0]],
                    [(0, 1), 10.0, names[0]],
                    [(0, 2), 10.0, names[1]],
                ]
            },
            boundnames=True,
        ),
    )
    project.add_package("npf_0d", mf.npf(k=np.squeeze(np.array([12.5]))))
    before = {key: package_digest(project.package(key)) for key in ("chd_np", "npf_0d")}
    project.save()

    reloaded = Project.load(tmp_path / "p")

    after = {key: package_digest(reloaded.package(key)) for key in ("chd_np", "npf_0d")}
    assert after == before


def _sim(grid, *, domain, surfaces=None, name="n"):
    """A minimal buildable simulation on ``grid``."""

    model = mf.gwf(
        "m",
        context=ModelContext(grid=grid, domain=domain, surfaces=surfaces),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
            mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    return mf.simulation(model, name=name, packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])


def test_an_object_dtype_domain_is_refused_at_freeze_time(tmp_path, triangulated_grid):
    """It froze "successfully" and the WHOLE snapshot was then unreadable.

    `ModelContext.domain` is typed `Any` and validated nowhere, so a
    `LayerBuildResult` passed in place of its `.idomain` -- or a list holding a
    `None` -- reached `np.asarray` and became an object array. `savez` pickles one
    and `load` refuses a pickle, so every package in the snapshot went with it,
    while `frozen_summary` said `readable: False` and named no cause.
    """

    spec = _sim(triangulated_grid, domain=[1, None, 1])
    with pytest.raises(FrozenSnapshotError, match=r"domain for 'm'"):
        freeze_simulation(spec, tmp_path, run_name="n", contexts={"m": spec.models[0].context})


def test_the_elevation_table_keeps_its_integer_column_labels(tmp_path, triangulated_grid):
    """The table's contract is integer columns, and the snapshot stringified them.

    `_elevation_table_from_layers` builds `{0: top, 1..nlay: bottoms}` and
    `GeoPackageSource.surface_value` looks them up by integer, so a thawed run
    rebuilt its packages and then died on the first `CellSurfaceOffset`, SFR or LAK
    elevation read -- the one failure `_surfaces_block` exists to prevent.
    """

    surfaces = triangulated_grid.gdf_vorPolys[["geometry"]].copy()
    for layer, elevation in enumerate([10.0, 5.0, 0.0]):
        surfaces[layer] = elevation
    spec = _sim(
        triangulated_grid,
        domain=np.ones((1, triangulated_grid.ncpl), dtype=int),
        surfaces=surfaces,
    )

    freeze_simulation(spec, tmp_path, run_name="n", contexts={"m": spec.models[0].context})
    thawed = thaw_simulation(tmp_path)

    assert list(thawed.models[0].context.surfaces.columns) == ["geometry", 0, 1, 2]


def test_a_column_label_the_converter_cannot_carry_degrades_to_a_string():
    """Keeping integer labels must not make an exotic one WORSE than before.

    A `MultiIndex` label converts to a tagged tuple, which is a dict, and a dict is
    not hashable as a column name; a label the converter refuses outright would raise
    inside `_surfaces_block`'s broad `except` and drop the surfaces entirely. Both
    are worse than the stringified label they replaced, so both fall back to it.
    """

    from myflopy.project.freeze import _column_label

    assert _column_label(0) == 0
    assert _column_label(np.int64(2)) == 2 and type(_column_label(np.int64(2))) is int
    assert _column_label("top") == "top"
    assert _column_label(("a", "b")) == "('a', 'b')"
    assert _column_label(1.5) == "1.5"


# --- the paths that did not go through the converter at all ----------------- #
#
# Fixing `_json_value` only helps a value that REACHES it. These four routes into
# JSON did not, and each was found after the first five fixes were already in --
# which is the shape of this whole entry: the boundary is wider than the function.


def test_a_numpy_value_inside_a_source_spec_serializes():
    """`_json_value`'s FIRST branch handed a source spec straight to its own `to_dict`.

    A source writes its fields raw, and `PackageSpec.inputs` is exactly where source
    specs live, so this was the one branch above the normalization that could still
    return numpy. Measured: `Project.validate()` returned `[]` and `Project.save()`
    then died inside `json.dumps` having already written the project and simulation
    files, so the next `Project.load` raised `FileNotFoundError` on the package it
    never reached.
    """

    from myflopy.sources import LiteralSource
    from myflopy.specs import PackageSpec

    package = PackageSpec(
        name="npf",
        builder=mf.npf(k=1.0).builder,
        options={},
        inputs={"k": LiteralSource(value=[np.float64(12.5), np.int64(3)])},
    )

    written = package.to_dict()

    assert {type(v).__name__ for v in written["inputs"]["k"]["value"]} == {"float", "int"}
    assert json.loads(json.dumps(written))


def test_a_project_with_a_numpy_valued_source_saves_and_reloads(tmp_path):
    """The half-written project directory is the part that made this expensive."""

    from myflopy.sources import LiteralSource
    from myflopy.specs import PackageSpec
    from myflopy.workspace import Project

    project = Project(tmp_path / "p", name="p")
    project.add_package(
        "np/src",
        PackageSpec(
            name="npf", builder=mf.npf(k=1.0).builder, options={},
            inputs={"k": LiteralSource(value=[np.int64(3)])},
        ),
    )
    report = project.save()

    # `pickled == 0` is the load-bearing half. Without it this passed with the fix
    # reverted, because the `TypeError` fallback quietly pickled the package instead --
    # a reload still returns something, so `is not None` cannot tell the two apart. It
    # is also the only assertion anywhere that reaches that fallback (ledger 266).
    assert report.pickled == 0, "the package fell back to a pickle instead of writing JSON"
    assert Project.load(tmp_path / "p").package("np/src") is not None


def test_a_package_json_cannot_write_falls_back_to_a_pickle(tmp_path, monkeypatch):
    """The fallback is reached for a `json.dumps` failure, not only a `to_dict` one.

    It guarded `to_dict()` alone and caught `ValueError`, but `json.dumps` raises
    `TypeError` on the line AFTER the guard -- so for the commonest way a package
    cannot be written as JSON the exception escaped `Project.save()` and broke its
    documented "persist what it can rather than refusing the whole project" contract,
    leaving a half-written directory (ledger 265).

    The failure is INJECTED because the numpy routes into it are now closed; this is a
    safety net, and a net nothing can reach is still worth knowing works.
    """

    from myflopy import workspace
    from myflopy.workspace import Project

    real = workspace._write_json
    calls = {"n": 0}

    def refuse_once(path, doc):
        if path.name == "npf.json" and calls["n"] == 0:
            calls["n"] += 1
            raise TypeError("Object of type Whatever is not JSON serializable")
        return real(path, doc)

    monkeypatch.setattr(workspace, "_write_json", refuse_once)
    project = Project(tmp_path / "p", name="p")
    project.add_package("npf", mf.npf(k=1.0))

    report = project.save()

    assert calls["n"] == 1, "the injected refusal never fired; the test proves nothing"
    assert report.pickled == 1
    assert not report.skipped
    assert Project.load(tmp_path / "p").package("npf") is not None


def test_numpy_metadata_does_not_stop_a_run_being_prepared(tmp_path):
    """`run.json` is written with `json.dumps` and no `default=`, by a path of its own.

    `PackageSpec.to_dict` normalizes `metadata`; `_spec_summary` does not call it, so
    both metadata dicts reached the manifest raw. `np.True_` in a package's metadata
    or a `np.int64` in the context's made `prepare_run` raise outright.
    """

    from myflopy.workspace import Project

    model = mf.gwf(
        "m",
        context=ModelContext(metadata={"zones": np.int64(4)}),
        packages=[
            mf.dis(nlay=1, nrow=2, ncol=2, delr=10.0, delc=10.0, top=10.0, botm=[0.0]),
            mf.npf(k=1.0).with_metadata(calibrated=np.True_),
            mf.ic(strt=1.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    spec = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    project = Project(tmp_path / "p", name="p")
    project.add_simulation(spec)

    run = project.prepare_run("r", "s")

    manifest = json.loads((run.workspace / "run.json").read_text())
    assert manifest["spec"]["models"][0]["context_metadata"] == {"zones": 4}


def test_a_numpy_crs_on_a_grid_spec_serializes(tmp_path):
    """`crs` documents itself as taking `2927` as well as `"EPSG:2927"`.

    pandas hands back `numpy.int64` for an integer column, so a CRS read from a
    config table is the ordinary producer of the value that broke the save.
    """

    from myflopy.specs import GridSpec

    written = GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=tmp_path / "b.gpkg"), crs=np.int64(2927)
    ).to_dict()

    assert type(written["crs"]) is int
    assert json.loads(json.dumps(written))


@pytest.mark.parametrize(
    "field, value",
    [("enabled", np.bool_(True)), ("name", np.str_("npf")), ("requires", (np.str_("dis"),))],
)
def test_a_numpy_scalar_field_does_not_drift_across_a_save(tmp_path, field, value):
    """`package_digest` read three fields raw while `from_dict` normalized them.

    `np.bool_` is what any numpy predicate returns -- `arr.any()`, `mask[i]`, a
    comparison -- so this is reachable without anyone meaning to. Before `enabled`
    went through the converter the package could not be saved at all, which means
    that fix is what turned a loud failure into a silent permanent false drift; both
    halves had to move together.
    """

    from myflopy.workspace import Project

    project = Project(tmp_path / "p", name="p")
    project.add_package("np", dataclasses.replace(mf.npf(k=1.0), **{field: value}))
    before = package_digest(project.package("np"))
    project.save()

    assert package_digest(Project.load(tmp_path / "p").package("np")) == before
