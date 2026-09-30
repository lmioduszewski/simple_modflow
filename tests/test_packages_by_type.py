"""Packages are found by TYPE, and a view that reads one package is told which (ledger 352).

FloPy's ``get_package(name)`` answers a package NAMED like its argument before it
looks at types. So on a model holding LAK packages named ``lak`` and ``lak_2`` --
what a second ``mf.lak`` beside a default-named first one builds -- ``model.fl.lak``
and ``model.package("lak")`` returned ``lak`` alone, and ``model.lak.results.stage``
read ``lak``'s lake as the model's while ``lak_2``'s vanished. With no package named
like its type (``lak_a``, ``lak_b``) FloPy returned a list and every LAK view died on
``'list' object has no attribute 'nlakes'``.

Beside it, the budget-file reader under every cell-budget view kept the FIRST record
of each saved time: a model's second RIV package -- Ten Trails' ``crisp`` beside
``riv``, 53% of its river exchange -- was absent from ``model.riv.results.q``,
``model.budget.riv`` and every map of them.

Real MF6 runs of the model in ``tests/test_obs_route_declared_identity.py`` (6 x 6
irregular Voronoi, one layer, CHD on the west and east edges, three 10-day periods)
carrying two packages each of LAK (cells 8-9, 25-26), SFR (13-16, 31-32), RIV (20;
21-22) and UZF (1-2, 3-4). ``named_like_the_type`` names the first of each like its
type (``lak``, ``sfr``, ``riv``, ``uzf``); ``named_apart`` names none so (``lak_a``,
``sfr_a`` ...). What each file holds is read with FloPy directly, never through the
code under test.
"""

from __future__ import annotations

import flopy
import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from tests.test_obs_route_declared_identity import DRYING, HELD, _lake, _run

REACHES = {"first": [13, 14, 15, 16], "second": [31, 32]}
RIVERS = {"first": [20], "second": [21, 22]}
UZF_CELLS = {"first": [1, 2], "second": [3, 4]}


def _stream(name, cells, inflow):
    """An SFR package of its own name and files: a straight reach chain on ``cells``."""

    last = len(cells) - 1
    packagedata = [
        [i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
         1 if i in (0, last) else 2, 1.0, 0]
        for i, c in enumerate(cells)
    ]
    connectiondata = [
        [i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
        for i in range(len(cells))
    ]
    return mf.sfr.flopy(
        name=name, filename=f"m.{name}", packagedata=packagedata, connectiondata=connectiondata,
        perioddata={0: [[0, "INFLOW", inflow]]},
        stage_filerecord=f"m.{name}.stg", budget_filerecord=f"m.{name}.bud")


def _river(name, cells, stage):
    return mf.riv.flopy(name=name, filename=f"m.{name}", save_flows=True,
                        stress_period_data={0: [[(0, c), stage, 50.0, stage - 2.0] for c in cells]})


def _uzf(name, cells, finf):
    return mf.uzf.flopy(
        name=name, filename=f"m.{name}", save_flows=True,
        packagedata=[[i, (0, c), 1, -1, 0.1, 0.1, 0.05, 0.30, 0.15, 4.0] for i, c in enumerate(cells)],
        perioddata={0: [[i, finf, 0.0, 1.0, 0.1, 0.0, 0.0, 0.0] for i in range(len(cells))]},
    )


def _two_of_each(root, names):
    """LAK, SFR, RIV and UZF twice each, the ``(first, second)`` names given per type."""

    model, _ = _run(root, [
        _lake(name=names["lak"][0], lakes=((DRYING, 114.0, "pa", 100.0),)),
        _lake(name=names["lak"][1], lakes=((HELD, 116.0, "pb", 100.0),)),
        _stream(names["sfr"][0], REACHES["first"], 50.0),
        _stream(names["sfr"][1], REACHES["second"], 20.0),
        _river(names["riv"][0], RIVERS["first"], 101.0),
        _river(names["riv"][1], RIVERS["second"], 100.5),
        _uzf(names["uzf"][0], UZF_CELLS["first"], 0.05),
        _uzf(names["uzf"][1], UZF_CELLS["second"], 0.02),
    ])
    return model


LIKE_THE_TYPE = {"lak": ("lak", "lak_2"), "sfr": ("sfr", "sfr_2"), "riv": ("riv", "crisp"),
                 "uzf": ("uzf", "uzf_2")}
APART = {"lak": ("lak_a", "lak_b"), "sfr": ("sfr_a", "sfr_b"), "riv": ("riv_a", "riv_b"),
         "uzf": ("uzf_a", "uzf_b")}


@pytest.fixture(scope="module")
def named_like_the_type(tmp_path_factory):
    return _two_of_each(tmp_path_factory.mktemp("like_the_type"), LIKE_THE_TYPE)


@pytest.fixture(scope="module")
def named_apart(tmp_path_factory):
    return _two_of_each(tmp_path_factory.mktemp("apart"), APART)


def _records(model, text):
    """Every record of budget term ``text`` in the model's budget file, read with FloPy."""

    cbc = flopy.utils.CellBudgetFile(f"{model.workspace}/m.cbc", precision="double")
    headers = cbc.recordarray
    frames = []
    for index in np.flatnonzero(np.char.strip(headers["text"].astype(str)) == text):
        record = cbc.get_record(int(index))
        frames.append(pd.DataFrame({
            "per": int(headers["kper"][index]) - 1,
            "pname": headers["paknam2"][index].decode().strip(),
            "cell": np.asarray(record["node"], dtype=np.int64) - 1,  # one layer: the node is the cell
            "q": np.asarray(record["q"]),
        }))
    return pd.concat(frames, ignore_index=True).sort_values(["per", "pname", "cell"], ignore_index=True)


def _stage_file(model, name):
    """One LAK or SFR package's stage file as MF6 wrote it: a row per period, a column per feature."""

    stages = flopy.utils.HeadFile(f"{model.workspace}/m.{name}.stg", text="STAGE", precision="double")
    return stages.get_alldata().reshape(len(stages.get_kstpkper()), -1)


def _read(frame, value):
    """A view's rows as ``(per, pname, cell, q)``, sorted like :func:`_records`."""

    return (
        frame.loc[:, ["per", "pname", "cell", value]]
        .rename(columns={value: "q"})
        .sort_values(["per", "pname", "cell"], ignore_index=True)
    )


# -- the accessors answer by TYPE ---------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("fixture", ["named_like_the_type", "named_apart"])
def test_model_fl_and_package_answer_by_type(fixture, request):
    """``model.fl.<pkg>`` and ``model.package(<type>)`` are every package of the type,
    in model order, on a live model and a reloaded one alike. FloPy's name-first lookup
    answered ``lak`` alone beside ``lak_2``."""

    live = request.getfixturevalue(fixture)
    names = LIKE_THE_TYPE if fixture == "named_like_the_type" else APART
    for model in (live, mf.load_mf6_run(live.workspace, verbosity_level=0)):
        for kind, pair in names.items():
            assert [p.package_name for p in getattr(model.fl, kind)] == list(pair), kind
            assert [p.package_name for p in model.package(kind)] == list(pair), kind
            assert model.package(pair[1]).package_name == pair[1]
        assert model.fl.chd.package_name == "chd"  # one of its type: the package itself


@pytest.mark.slow
def test_a_string_is_a_type_or_a_name_never_a_prefix(named_like_the_type):
    """FloPy's last resort matched a name PREFIX: ``"dis"`` answered the DISV."""

    model = named_like_the_type
    with pytest.raises(AttributeError, match="'dis'"):
        model.package("dis")
    assert model.package("disv").package_type == "disv"
    assert model.package("tdis").package_type == "tdis"  # the simulation's, as before
    with pytest.raises(AttributeError, match="no GHB package"):
        model.fl.ghb  # noqa: B018 -- the property is the read


@pytest.mark.slow
def test_a_package_merely_named_like_a_type_is_not_that_type(tmp_path):
    """A GHB NAMED ``lak`` on a model with no LAK. ``model.fl.lak`` is type-only, but a
    property raising ``AttributeError`` makes Python retry the name through
    ``FlopyPackages.__getattr__`` -- whose name lookup answered the GHB."""

    model, _ = _run(tmp_path, [mf.ghb.flopy(name="lak", filename="m.ghb_named_lak",
                                            stress_period_data={0: [[(0, 7), 101.0, 10.0]]})])
    with pytest.raises(AttributeError, match="no LAK package"):
        model.fl.lak  # noqa: B018 -- the property is the read
    assert model.fl.ghb.package_name == "lak"
    assert model.package("lak").package_type == "ghb"  # a NAME, once no type matches
    with pytest.raises(AttributeError, match="no LAK package"):
        model.lak.results.stage.get()


# -- a view that reads ONE package refuses several, naming the spelling -------------

#: Every LAK, SFR and UZF view that reads one package's numbered features.
ONE_PACKAGE_VIEWS = {
    "lak.results.stage.get": lambda m: m.lak.results.stage.get(),
    "lak.results.stage_change.get": lambda m: m.lak.results.stage_change.get(),
    "lak.results.q.get": lambda m: m.lak.results.q.get(),
    "lak.results.q.budget.get": lambda m: m.lak.results.q.budget.get(),
    "lak.budget.get": lambda m: m.lak.budget.get(),
    "lak.budget.types": lambda m: m.lak.budget.types,
    "lak.connections.get": lambda m: m.lak.connections.get(),
    "lak.inputs.get": lambda m: m.lak.inputs.get(),
    "sfr.results.stage.get": lambda m: m.sfr.results.stage.get(),
    "sfr.results.q.get": lambda m: m.sfr.results.q.get(),
    "sfr.results.profile.get": lambda m: m.sfr.results.profile.get(),
    "sfr.budget.get": lambda m: m.sfr.budget.get(),
    "sfr.inputs.get": lambda m: m.sfr.inputs.get(),
    "uzf.inputs.finf.get": lambda m: m.uzf.inputs.finf.get(),
}


@pytest.mark.slow
@pytest.mark.parametrize("fixture", ["named_like_the_type", "named_apart"])
def test_a_one_package_view_refuses_several_naming_them_and_the_spelling(fixture, request):
    """Beside ``lak_2`` every LAK view read ``lak``'s lake as the model's, in silence;
    beside nothing named ``lak`` they died on FloPy's list. Each refuses now, naming
    both packages and ``model.lak(package=...)``."""

    live = request.getfixturevalue(fixture)
    names = LIKE_THE_TYPE if fixture == "named_like_the_type" else APART
    # the live run and the reloaded one -- `mf.load_mf6_run(ws).lak.results.stage.get()`
    # is how the defect was reported
    for model in (live, mf.load_mf6_run(live.workspace, verbosity_level=0)):
        for label, read in ONE_PACKAGE_VIEWS.items():
            kind = label.split(".")[0]
            first, second = names[kind]
            with pytest.raises(mf.SeveralPackagesError) as caught:
                read(model)
            message = str(caught.value)
            assert f"({first!r}, {second!r})" in message, label
            assert f"model.{kind}(package={second!r})" in message, label
        for kind in ("lak", "sfr"):
            first, second = names[kind]
            with pytest.raises(mf.SeveralPackagesError, match=rf"model\.outputs\.{kind}\(package={second!r}\)"):
                getattr(model.outputs, kind).stage.get()


@pytest.mark.slow
def test_each_named_package_reads_its_own_files(named_like_the_type):
    """``model.lak(package=...)`` / ``model.sfr(package=...)`` / ``model.uzf(package=...)``
    read that package's own stage, budget and input -- each against the file MF6 wrote.
    The package NAMED like its type is picked by its name, never taken for the type."""

    model = named_like_the_type
    cells = {"lak": (DRYING, HELD), "sfr": (REACHES["first"], REACHES["second"])}
    for name, lake_cells in zip(LIKE_THE_TYPE["lak"], cells["lak"], strict=True):
        view = model.lak(package=name)
        stage = view.results.stage.get().drop_duplicates(["per", "lake"]).sort_values("per")
        np.testing.assert_array_equal(stage["stage"].to_numpy(), _stage_file(model, name)[:, 0])
        assert sorted(view.connections.get()["cell"]) == sorted(lake_cells)
        lake_budget = flopy.utils.CellBudgetFile(f"{model.workspace}/m.{name}.bud", precision="double")
        exchange = np.concatenate([record["q"] for record in lake_budget.get_data(text="GWF")])
        np.testing.assert_array_equal(np.sort(view.results.q.get()["q_lake"].to_numpy()), np.sort(exchange))
        assert set(view.budget.get(term="GWF")["cell"]) == set(lake_cells)
    for name, reach_cells in zip(LIKE_THE_TYPE["sfr"], cells["sfr"], strict=True):
        view = model.sfr(package=name)
        stage = view.results.stage.get().pivot(index="per", columns="reach", values="stage")
        np.testing.assert_array_equal(stage.to_numpy(), _stage_file(model, name))
        truth = _records(model, "SFR").query("pname == @name.upper()").reset_index(drop=True)
        np.testing.assert_array_equal(_read(view.results.q.get().assign(pname=name.upper()), "q_gwf")["q"],
                                      truth["q"])
        assert sorted(view.results.profile.get(per=2)["cell"]) == sorted(reach_cells)
    for name, uzf_cells in zip(LIKE_THE_TYPE["uzf"], UZF_CELLS.values(), strict=True):
        finf = model.uzf(package=name).inputs.finf.get(per=0)
        assert sorted(finf["cell"]) == sorted(uzf_cells)
        truth = _records(model, "UZF-GWRCH").query("pname == @name.upper()").reset_index(drop=True)
        pd.testing.assert_frame_equal(_read(model.uzf(package=name).results.gwrch.get(), "gwrch"), truth)
    # the package NAMED like its type, by its name: `lak`, never both LAKs
    assert set(model.lak(package="LAK").results.stage.get()["cell"]) == set(DRYING)
    loaded = mf.load_mf6_run(model.workspace, verbosity_level=0)
    stage = loaded.lak(package="lak_2").results.stage.get().drop_duplicates(["per", "lake"])
    np.testing.assert_array_equal(stage.sort_values("per")["stage"].to_numpy(), _stage_file(model, "lak_2")[:, 0])


@pytest.mark.slow
def test_a_package_name_is_checked_at_the_call(named_like_the_type):
    model = named_like_the_type
    with pytest.raises(ValueError, match=r"package='nope' names no LAK package.*'lak', 'lak_2'"):
        model.lak(package="nope")
    with pytest.raises(ValueError, match=r"package='sfr' names no LAK package"):
        model.lak(package="sfr")  # a package, of another type
    with pytest.raises(TypeError, match=r"model\.lak\(package=\.\.\.\) takes a LAK package's name"):
        model.lak(package=2)
    with pytest.raises(ValueError, match=r"package='nope' names no SFR package"):
        model.outputs.sfr(package="nope")
    with pytest.raises(ValueError, match=r"package='nope' names no UZF package"):
        model.uzf(package="nope")


# -- a view whose rows are CELLS reads every package, each row tagged ---------------


@pytest.mark.slow
def test_every_package_of_a_cell_term_is_read_from_the_budget_file(named_like_the_type):
    """The reader kept the FIRST record of each saved time: ``crisp``'s flows were
    absent from ``model.riv.results.q`` and ``model.budget.riv`` beside ``riv`` (Ten
    Trails: 53% of its river exchange). Every record is read now, tagged ``pname``."""

    model = named_like_the_type
    truth = _records(model, "RIV")
    assert set(truth["pname"]) == {"RIV", "CRISP"}
    pd.testing.assert_frame_equal(_read(model.riv.results.q.get(), "q_gwf"), truth)
    pd.testing.assert_frame_equal(_read(model.budget.riv.get(), "q"), truth)
    wide = model.riv.results.q.wide()
    assert sorted(wide["cell"]) == RIVERS["first"] + RIVERS["second"]
    # UZF's cell terms, unbound, are every UZF package's too.
    pd.testing.assert_frame_equal(_read(model.uzf.results.gwrch.get(), "gwrch"),
                                  _records(model, "UZF-GWRCH"))


@pytest.mark.slow
def test_the_combined_surface_water_view_reads_every_lak_and_sfr_package(named_like_the_type):
    frame = named_like_the_type.surface_water.results.q.get(per=2)
    assert set(frame.query("source == 'lak'")["cell"]) == set(DRYING) | set(HELD)
    assert set(frame.query("source == 'sfr'")["cell"]) == set(REACHES["first"]) | set(REACHES["second"])


@pytest.mark.slow
def test_the_all_packages_table_and_mosaic_take_each_package(named_like_the_type):
    """``summary(detail='data')`` looked explorers up by package NAME: ``CRISP`` read
    blank and ``RIV`` counted ``CRISP``'s records too; a second LAK has its own row."""

    model = named_like_the_type
    counts = model.packages.summary(detail="data").set_index("package")["records"]
    assert counts["RIV"] == 1 and counts["CRISP"] == 2       # one period of 1 and 2 cells
    assert counts["LAK"] == counts["LAK_2"] == 6              # 3 periods x 2 connections
    assert counts["SFR"] == 12 and counts["SFR_2"] == 6       # 3 periods x 4 and 2 reaches
    figure = model.packages.mosaic(packages=["lak", "lak_2", "sfr_2", "riv", "crisp"])
    titles = [annotation.text for annotation in figure.layout.annotations]
    assert titles == ["lak", "lak_2", "sfr_2", "riv"], titles  # `riv` is every RIV's panel


# -- the other readers that asked FloPy for "the" package of a type -----------------


@pytest.mark.slow
def test_attach_flopy_obs_takes_the_named_package_or_refuses_several(named_like_the_type):
    """With no ``package=`` the OBS went on FloPy's ``model.gwf.lak`` -- ``lak`` beside
    ``lak_2``, in silence. Attached in memory to a reloaded copy; nothing is written."""

    fresh = mf.load_mf6_run(named_like_the_type.workspace, verbosity_level=0)
    with pytest.raises(ValueError, match=r"2 LAK packages \('lak', 'lak_2'\).*attaches to ONE"):
        mf.LakeStageTargets(locations={"x": 0}).attach_flopy_obs(fresh, filename="x.obs")
    with pytest.raises(ValueError, match=r"2 RIV|2 SFR packages \('sfr', 'sfr_2'\)"):
        mf.SfrStageTargets(locations={"x": 0}).attach_flopy_obs(fresh, filename="y.obs")
    attached = mf.LakeStageTargets(locations={"x": 0}, package="lak_2").attach_flopy_obs(
        fresh, filename="z.obs")
    assert attached.parent_file.package_name == "lak_2"


@pytest.mark.slow
def test_an_artifact_is_one_package_by_name_and_a_type_of_several_is_refused(named_like_the_type, named_apart):
    """An artifact is captured by the package's exact NAME first (the builders capture
    the package they just built by its name); a type of several is refused. FloPy's list
    reached the type check as package type ``''``."""

    from myflopy.project.components import build_package_artifact

    def rows(artifact):
        return sum(len(records) for records in artifact.package_data["stress_period_data"].values())

    like = mf.load_mf6_run(named_like_the_type.workspace, verbosity_level=0)
    assert rows(build_package_artifact(like, artifact_id="a", package_name="crisp")) == 2
    assert rows(build_package_artifact(like, artifact_id="b", package_name="riv")) == 1  # by NAME
    apart = mf.load_mf6_run(named_apart.workspace, verbosity_level=0)
    with pytest.raises(ValueError, match=r"2 RIV packages \('riv_a', 'riv_b'\).*ONE package"):
        build_package_artifact(apart, artifact_id="c", package_name="riv")


@pytest.mark.slow
def test_split_lake_groups_cover_every_lak_package(named_like_the_type):
    """The splitter keeps each lake's cells in one part; it grouped only ``lak``'s."""

    from myflopy.modflow.mf6.parallel import _lake_column_groups

    groups = sorted(sorted(int(cell) for cell in group) for group in _lake_column_groups(named_like_the_type.gwf))
    assert groups == [sorted(DRYING), sorted(HELD)]


@pytest.mark.slow
def test_the_diff_tier_refuses_several_rather_than_reading_none(named_like_the_type, named_apart):
    """The diff and group tiers tolerate a MISSING package by catching ``ValueError``
    among others (``PACKAGE_TABLE_UNAVAILABLE``). The refusal is a ``LookupError`` so
    that no such handler compares a two-LAK model as having no LAK at all."""

    from myflopy.project.model_group import ModelGroup

    diff = ModelGroup({"like": named_like_the_type, "apart": named_apart}, reference="like").diff()
    with pytest.raises(mf.SeveralPackagesError, match=r"2 LAK packages \('lak', 'lak_2'\)"):
        diff.packages.lak.inputs.connections()
    assert not issubclass(mf.SeveralPackagesError, ValueError)
