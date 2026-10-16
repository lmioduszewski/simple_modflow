"""A target's group is its PEST observation group (ledger 358).

Every target class carried a ``group`` label that reached no PEST file: pyEMU gave
each simulated column a group of its own (``oname:hds_otype:lst_usecol:w1``), so
every well was a group of one and PEST++'s phi by group -- and the review's group
summaries -- could not see the user's grouping. Now a MEASURED observation is in its
target's label (made a PEST name), else in its set's prefix; every row no target
measures is in ``<prefix>_unmeasured`` and a forecast in ``<prefix>``, so neither
ever joins a user's group. A label that cannot be a PEST group is refused, naming it.

The review drew one time-series panel per GROUP, which was one location only while
pyEMU gave every column a group; it draws one per target LOCATION now, and a run
built before this commit -- one group per column -- reviews the same.

Every test builds a real control file on the real MF6 run of
``test_pest_handed_targets.py``'s toy (6 x 6 Voronoi, three 10-day periods, heads,
a stream, a three-cell DRN zone and a lake); the review tests run PESTPP-IES on it.
"""

from __future__ import annotations

import re
import shutil
import warnings

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.pest.ies import IesResults
from tests.test_pest_handed_targets import DRN_CELLS, _build, run_toy

pyemu = pytest.importorskip("pyemu")


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """``test_pest_handed_targets.py``'s model, run."""

    return run_toy(tmp_path_factory.mktemp("pest_groups"))

HEADS = pd.DataFrame(
    {"time": [0, 1, 2], "W2": [91.0, 92.0, None], "W1": [93.0, None, 94.0], "W10": [95.0, 96.0, 97.0]}
)


def _heads(groups, *, values=HEADS):
    locations = pd.DataFrame(
        {"name": ["W2", "W1", "W10"], "cell": [1, 2, 34], "layer": 0, "group": groups}
    )
    return mf.HeadTargets(locations=locations, values=values)


def _zone(group):
    return mf.DrnFlowTargets(
        locations=pd.DataFrame({"name": ["zone"], "cells": [DRN_CELLS], "group": [group]}),
        values=pd.DataFrame({"time": [2.0], "zone": [-3.0]}),
    )


def _members(obs: pd.DataFrame, group: str) -> list[str]:
    return sorted(obs.index[obs["obgnme"] == group])


def _measured(obs: pd.DataFrame, column: str, prefix: str = "hds") -> list[str]:
    stem = f"oname:{prefix}_otype:lst_usecol:{column}_"
    return sorted(name for name in obs.index if name.startswith(stem) and obs.loc[name, "weight"] > 0)


# -- the build: which group each observation is in ----------------------------


def test_a_measured_observation_is_in_its_targets_group_and_nothing_else_joins_it(toy, tmp_path):
    """The user's group holds exactly what was measured under its label: the
    unmeasured periods of a labeled well, and a forecast at a labeled well, do not."""

    heads = _heads([" North Wells", "North Wells  ", "South"])
    pond = mf.LakeStageTargets(locations={"pond": 0}, values=pd.DataFrame({"time": [1.0], "pond": [99.0]}))
    flow = mf.SfrFlowTargets(locations={"a": 0}, values=pd.DataFrame({"time": [1.0], "a": [-7.0]}))
    # A forecast at a well labeled like the measured ones, with a placeholder value.
    fore = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["P"], "cell": [23], "layer": 0, "group": ["North Wells"]}),
        values=pd.DataFrame({"time": [2], "P": [100.0]}),
    )
    cal, pst = _build(toy, tmp_path / "pest", observe=(heads, _zone("Seeps"), pond, flow), forecast=(fore,))
    obs = pst.observation_data

    # The label, made a PEST name -- and " North Wells" / "North Wells  " are ONE
    # label (surrounding whitespace is not part of it): no refusal.
    assert _members(obs, "north_wells") == sorted(_measured(obs, "w2") + _measured(obs, "w1"))
    assert len(_members(obs, "north_wells")) == 4
    assert _members(obs, "south") == _measured(obs, "w10")
    assert _members(obs, "seeps") == _measured(obs, "zone", "drn_flow")
    # No label: the set's prefix.
    assert _members(obs, "lak_stage") == _measured(obs, "pond", "lak_stage")
    assert _members(obs, "sfr_flow") == _measured(obs, "a", "sfr_flow")
    # Every other row is its set's `_unmeasured`, at weight 0.
    for prefix in ("hds", "drn_flow", "lak_stage", "sfr_flow", "fore1"):
        rows = obs.loc[obs["obgnme"] == f"{prefix}_unmeasured"]
        assert len(rows) and (rows["weight"] == 0).all(), prefix
    assert sorted(obs.loc[obs["obgnme"] == "hds_unmeasured"].index) == [
        "oname:hds_otype:lst_usecol:w1_per:1", "oname:hds_otype:lst_usecol:w2_per:2",
    ]
    # The forecast is its prefix's, not its label's.
    assert _members(obs, "fore1") == ["oname:fore1_otype:lst_usecol:p_per:2"]
    assert pst.pestpp_options["forecasts"] == "oname:fore1_otype:lst_usecol:p_per:2"
    # The values snapshot carries each target's group and column for the review.
    snapshot = pd.read_csv(cal.template_workspace / "hds_target_values.csv")
    assert snapshot.set_index("obsnme").loc["oname:hds_otype:lst_usecol:w10_per:0", "obgnme"] == "south"
    assert set(snapshot["column"]) == {"W1", "W2", "W10"}


def test_a_target_with_no_label_or_a_blank_one_is_in_its_sets_prefix(toy, tmp_path):
    """Two head sets, one unlabeled and one with a blank and a missing label: phi
    splits by set out of the box, each set under the prefix it was registered with."""

    deep = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["D1", "D2"], "cell": [7, 28], "layer": 0, "group": ["  ", None]}),
        values=pd.DataFrame({"time": [0], "D1": [90.0], "D2": [91.0]}),
    )
    cal = toy.pest("prefixes", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(_heads([None, None, None]))
    cal.observe(deep, prefix="deep")
    pst = cal.build("prefixes.pst")
    obs = pst.observation_data

    assert _members(obs, "hds") == sorted(
        _measured(obs, "w1") + _measured(obs, "w2") + _measured(obs, "w10")
    )
    assert _members(obs, "deep") == ["oname:deep_otype:lst_usecol:d1_per:0", "oname:deep_otype:lst_usecol:d2_per:0"]
    assert set(obs.loc[obs["weight"] > 0, "obgnme"]) == {"hds", "deep"}


# -- the refusals ---------------------------------------------------------------


def test_two_labels_that_become_one_group_name_are_refused_naming_both(toy, tmp_path):
    with pytest.raises(ValueError, match=r"'North Wells' and 'north_wells' both become 'north_wells'"):
        _build(toy, tmp_path / "pest", observe=(_heads(["North Wells", "north_wells", "south"]),))


@pytest.mark.parametrize("label", ["Greater Seattle", "lessard", "l_upper", "G_deep"])
def test_a_label_pestpp_ies_reads_as_an_inequality_is_refused(toy, tmp_path, label):
    """PESTPP-IES reads a group starting ``l_``/``less``/``g_``/``greater`` as an
    inequality: measured, a ``greater...`` group of heads above their targets scored
    phi 0 where ``north`` scored 424.864 (``docs/reproducers/pest_observation_groups.py``)."""

    with pytest.raises(ValueError, match="INEQUALITY"):
        _build(toy, tmp_path / "pest", observe=(_heads([label, label, "south"]),))


def test_a_set_prefix_pestpp_ies_reads_as_an_inequality_is_refused(toy, tmp_path):
    cal = toy.pest("ineq", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(_heads([None, None, None]), prefix="less_wells")
    with pytest.raises(ValueError, match=r"observation set 'less_wells'.*INEQUALITY"):
        cal.build("ineq.pst")


def _default_prefixes():
    """Every observation spec's DEFAULT prefix, found from ``pest.specs`` rather than listed."""

    import dataclasses

    from myflopy.modflow.mf6.pest import specs

    found = {}
    for name in dir(specs):
        cls = getattr(specs, name)
        if isinstance(cls, type) and dataclasses.is_dataclass(cls) and name.endswith("ObservationSpec"):
            fields = {field.name: field for field in dataclasses.fields(cls)}
            if "prefix" in fields:
                found[name] = fields["prefix"].default
    return found


def _tags_matching(group: str, tags) -> list[str]:
    """The ``ies_phi_factor_file`` tags PESTPP-IES maps ``group`` to: each tag found
    ANYWHERE in the group's name, ignoring case (measured with 5.2.16, ledger 358)."""

    return [tag for tag in tags if tag.lower() in group.lower()]


def test_no_default_prefix_is_inside_another_sets_group():
    """A phi-factor tag matches anywhere in a group's name, and a lake set's default
    ``stage`` was inside an SFR set's ``sfr_stage``, so no tag could give the lake its
    own factor: PESTPP-IES stopped, "group 'SFR_STAGE' mapped to multiple tags"
    (ledger 359). The lake default is ``lak_stage`` now (ledger 372), and no
    default prefix is inside another family's groups, ``<prefix>_unmeasured`` included."""

    defaults = _default_prefixes()
    assert defaults["LakeStageObservationSpec"] == "lak_stage"
    assert defaults["RivFlowObservationSpec"] == "riv_flow"
    # The difference sets (ledger 565) are one group each, checked too.
    assert defaults["HeadDifferenceObservationSpec"] == "hdiff"
    assert defaults["ConcDifferenceObservationSpec"] == "cdiff"
    assert len(defaults) == 9 and len(set(defaults.values())) == 9, defaults
    for owner, prefix in defaults.items():
        for other, theirs in defaults.items():
            if other == owner:
                continue
            for group in (theirs, f"{theirs}_unmeasured"):
                assert prefix not in group, (owner, prefix, other, group)


def test_the_default_groups_take_one_phi_factor_tag_each(toy, tmp_path):
    """Heads, a lake and a reach-stage set under their default prefixes: the control
    file's groups, and each group matched by exactly the one tag of its own set."""

    pond = mf.LakeStageTargets(locations={"pond": 0}, values=pd.DataFrame(
        {"time": [1.0, 2.0], "pond": [99.0, 98.0]}))
    gage = mf.SfrStageTargets(locations={"r1": 1}, values=pd.DataFrame(
        {"time": [0.0, 1.0], "r1": [95.0, 96.0]}))
    _cal, pst = _build(toy, tmp_path / "pest", observe=(_heads([None, None, None]), pond, gage))
    obs = pst.observation_data
    measured = sorted(set(obs.loc[obs["weight"] > 0, "obgnme"]))
    assert measured == ["hds", "lak_stage", "sfr_stage"]
    assert sorted(set(obs["obgnme"])) == [
        "hds", "hds_unmeasured", "lak_stage", "lak_stage_unmeasured", "sfr_stage",
        "sfr_stage_unmeasured"]
    assert _members(obs, "lak_stage") == _measured(obs, "pond", "lak_stage")
    for group in measured:
        assert _tags_matching(group, measured) == [group], group
    # What the old default did: the tag `stage` matched the reach set's group as well.
    assert _tags_matching("sfr_stage", ["hds", "stage", "sfr_stage"]) == ["stage", "sfr_stage"]


@pytest.mark.parametrize("label", [1, 2.0, True])
def test_a_label_that_is_not_text_is_refused_naming_it(toy, tmp_path, label):
    with pytest.raises(ValueError, match=r"not text: 'W2' in set 'hds': " + repr(label).replace(".", r"\.")):
        _build(toy, tmp_path / "pest", observe=(_heads([label, "a", "b"]),))


@pytest.mark.parametrize("label", ["lak_stage", "hds_unmeasured", "Fore1"])
def test_a_label_taking_a_group_the_build_gives_a_set_is_refused(toy, tmp_path, label):
    pond = mf.LakeStageTargets(locations={"pond": 0}, values=pd.DataFrame({"time": [1.0], "pond": [99.0]}))
    fore = mf.HeadTargets(locations={"P": 23}, values=pd.DataFrame({"time": [2], "P": [None]}))
    with pytest.raises(ValueError, match="take a name the build gives another group"):
        _build(toy, tmp_path / "pest", observe=(_heads([label, "a", "b"]), pond), forecast=(fore,))


@pytest.mark.parametrize("label", ["nan", "None", "<NA>"])
def test_the_label_astype_str_writes_for_a_blank_is_refused(toy, tmp_path, label):
    """``astype(str)`` -- the spelling the not-text refusal used to advise -- writes
    a missing label as ``'nan'``/``'None'``/``'<NA>'``, which made a group of that
    name of every well left unlabelled (ledger 359)."""

    with pytest.raises(ValueError, match=(
            r"1 target\(s\) carry the group label 'nan', 'None' or '<NA>': 'W10' in set "
            r"'hds': '" + re.escape(label) + r"'.*what `astype\(str\)` writes for a MISSING")):
        _build(toy, tmp_path / "pest", observe=(_heads(["north", "north", label]),))


def test_the_not_text_refusals_advice_keeps_a_blank_label_in_the_sets_group(toy, tmp_path):
    """A zone id column read with a blank is floats. The refusal names a spelling,
    and that spelling must keep the blank well in its set's group: the advice was
    ``astype(str)``, which put it in a group named ``nan`` (ledger 359)."""

    locations = pd.DataFrame(
        {"name": ["W2", "W1", "W10"], "cell": [1, 2, 34], "layer": 0, "group": [1.0, 1.0, np.nan]}
    )
    with pytest.raises(ValueError, match=(
            r"^2 target\(s\) carry a group label that is not text: 'W2' in set 'hds': 1\.0 "
            r"\(float\); 'W1' in set 'hds': 1\.0 \(float\)\..*"
            r"`\.astype\('Int64'\)\.astype\('string'\)`")):
        _build(toy, tmp_path / "refused", observe=(mf.HeadTargets(locations=locations, values=HEADS),))

    advised = locations.assign(group=locations["group"].astype("Int64").astype("string"))
    _cal, pst = _build(toy, tmp_path / "advised", observe=(mf.HeadTargets(locations=advised, values=HEADS),))
    obs = pst.observation_data
    assert _members(obs, "1") == sorted(_measured(obs, "w2") + _measured(obs, "w1"))
    assert _members(obs, "hds") == _measured(obs, "w10")


def test_a_zone_given_two_labels_is_refused(toy, tmp_path):
    zone = mf.DrnFlowTargets(
        locations=pd.DataFrame({
            "name": ["zone", "zone"], "cells": [DRN_CELLS[:2], DRN_CELLS[2:]], "group": ["seeps", "springs"],
        }),
        values=pd.DataFrame({"time": [2.0], "zone": [-3.0]}),
    )
    with pytest.raises(ValueError, match=r"'zone' \('seeps', 'springs'\).*DIFFERENT group labels"):
        _build(toy, tmp_path / "pest", observe=(zone,))


def test_a_zone_whose_first_row_is_unlabelled_reads_one_label_in_compare_and_in_pest(toy, tmp_path):
    """Rows with no label take the others' in the build; ``compare()`` read the FIRST
    row's, so it showed ``<NA>`` -- and no weight -- for a zone PEST put in ``seeps``
    at weight 2 (ledger 359)."""

    zone = mf.DrnFlowTargets(
        locations=pd.DataFrame({
            "name": ["zone", "zone"], "cells": [DRN_CELLS[:2], DRN_CELLS[2:]],
            "group": [None, "Seeps"], "weight": [np.nan, 2.0],
        }),
        values=pd.DataFrame({"time": [2.0], "zone": [-3.0]}),
    )
    compared = zone.compare(toy)
    assert list(compared["group"]) == ["Seeps"]
    assert list(compared["weight"]) == [2.0]
    _cal, pst = _build(toy, tmp_path / "pest", observe=(zone,))
    measured = pst.observation_data.loc[_measured(pst.observation_data, "zone", "drn_flow")]
    assert list(measured["obgnme"]) == ["seeps"]
    assert list(measured["weight"]) == [2.0]


def test_a_concentration_targets_label_is_its_group(tmp_path):
    """``ConcTargets`` shares ``_point_target_frame`` with heads, through the same
    assignment: a labelled well's measured observations are in its label, an
    unlabelled one's in ``conc``, and the forecast set's in its prefix."""

    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
    from myflopy.modflow.mf6.canonical_transport import build_canonical_transport_calibration_demo

    with warnings.catch_warnings():
        # The canonical lake is known unexcavated (tests/conftest.py, ledger 325).
        warnings.filterwarnings("ignore", message=r"mf\.lak: lake\(s\) NOT excavated")
        demo = build_canonical_transport_calibration_demo(
            tmp_path / "model", config=CanonicalModelConfig.testing()
        )
    source = demo.conc_targets
    locations = source.locations_gdf.loc[:, ["name", "layer", "cell"]].reset_index(drop=True)
    plume = set(locations["name"].iloc[::2])
    locations["group"] = [("Plume Core" if name in plume else None) for name in locations["name"]]
    cal = demo.model.pest("conc_groups", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.05, 2.0), physical=(0.01, 300.0))
    cal.observe(mf.ConcTargets(locations=locations, values=source.to_wide()))
    cal.forecast(demo.forecast_targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obs = cal.build("conc_groups.pst").observation_data

    measured = obs.loc[obs["weight"] > 0]
    assert set(measured["obgnme"]) == {"plume_core", "conc"}
    for name, group in zip(measured.index, measured["obgnme"], strict=True):
        well = re.search(r"usecol:(.+)_per:", name).group(1)
        assert group == ("plume_core" if well in {w.lower() for w in plume} else "conc"), name
    assert set(obs.loc[obs["weight"] == 0, "obgnme"]) <= {"conc_unmeasured", "fore1", "fore1_unmeasured"}
    assert (obs.loc[cal.pst.pestpp_options["forecasts"].split(","), "obgnme"] == "fore1").all()


# -- PESTPP-IES and the review ---------------------------------------------------


@pytest.fixture(scope="module")
def prior_run(toy, tmp_path_factory):
    """W2 and W1 measured under ``north``, W10 under ``south``, run through PESTPP-IES's prior."""

    workspace = tmp_path_factory.mktemp("groups_ies")
    cal = toy.pest("grouped", workspace=workspace / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(_heads(["north", "north", "south"]))
    cal.build("grouped.pst")
    return cal, cal.prior(reals=3)


def _base_residual_phi(cal, columns) -> float:
    """phi of the base run's own heads against the targets, by hand."""

    simulated = pd.read_csv(cal.template_workspace / "hds_simulated_heads.csv").set_index("per")
    targets = HEADS.set_index("time")
    return float(sum(
        (simulated.loc[per, column] - value) ** 2
        for column in columns for per, value in targets[column].dropna().items()
    ))


def test_pestpp_ies_reports_phi_under_the_users_groups(prior_run):
    cal, ies = prior_run
    by_group = pd.read_csv(ies.workspace / "grouped.phi.group.csv")
    base = by_group.loc[by_group["par_realization"].astype(str) == "base"].iloc[0]

    groups = [column for column in by_group.columns if column not in
              ("iteration", "total_runs", "obs_realization", "par_realization")]
    assert sorted(groups) == ["hds_unmeasured", "north", "south"]
    assert base["north"] == pytest.approx(_base_residual_phi(cal, ["W2", "W1"]), rel=1e-4)
    assert base["south"] == pytest.approx(_base_residual_phi(cal, ["W10"]), rel=1e-4)
    assert base["hds_unmeasured"] == 0.0
    # The review's group summaries read the same groups.
    assert set(ies.phi_contributions().index) == {"north", "south"}
    assert set(ies.conflict()["obgnme"]) == {"north", "south"}


def _panels(figure):
    """Per subplot title: the y of its 'measured' markers and how many lines it draws."""

    titles = [annotation.text for annotation in figure.layout.annotations]
    axes = {}
    for trace in figure.data:
        axes.setdefault(trace.yaxis, []).append(trace)
    out = {}
    for title, (_axis, traces) in zip(titles, sorted(axes.items(), key=lambda kv: int(kv[0][1:] or 1)), strict=True):
        measured = [list(trace.y) for trace in traces if trace.name == "measured"]
        lines = [trace for trace in traces if trace.mode == "lines"]
        out[title] = (measured[0], len(lines), len(lines[0].y))
    return out


def test_the_review_draws_one_series_per_location_not_per_group(prior_run):
    """``north`` holds two wells; drawn per group they were one line zig-zagging
    between them, their measured markers merged."""

    _cal, ies = prior_run
    for figure in (ies.plot_vs_obs(), ies.plot_prior_vs_obs(show_conflict=False)):
        panels = _panels(figure)
        assert set(panels) == {"hds: W2 [north]", "hds: W1 [north]", "hds: W10 [south]"}
        # Each well's own measurements, and a line through all three periods.
        assert panels["hds: W2 [north]"][0] == [91.0, 92.0]
        assert panels["hds: W1 [north]"][0] == [93.0, 94.0]
        assert panels["hds: W10 [south]"][0] == [95.0, 96.0, 97.0]
        assert {length for _measured, _count, length in panels.values()} == {3}
    assert set(_panels(ies.plot_vs_obs(groups=["NORTH"]))) == {"hds: W2 [north]", "hds: W1 [north]"}
    assert len(_panels(ies.plot_vs_obs(max_locations=1))) == 1
    static = ies.plot_vs_obs(backend="matplotlib")
    assert {ax.get_title(loc="left") for ax in static.axes} == set(panels)


@pytest.mark.parametrize("era", ["before_358", "before_340"])
def test_a_run_built_before_this_commit_reviews_per_location_too(prior_run, tmp_path, era):
    """The layout the commits before wrote: every observation in pyEMU's per-column
    group, and a values snapshot with no ``column``/``obgnme`` (before ledger 340, no
    ``obsnme`` either). One panel per column, each holding its own measurements."""

    _cal, ies = prior_run
    old = tmp_path / "old"
    shutil.copytree(ies.workspace, old)
    pst = pyemu.Pst(str(old / "grouped.pst"))
    obs = pst.observation_data
    obs["obgnme"] = [name.rsplit("_", 1)[0] for name in obs.index]
    pst.write(str(old / "grouped.pst"), version=2)
    snapshot = pd.read_csv(old / "hds_target_values.csv").drop(columns=["column", "obgnme"])
    if era == "before_340":
        snapshot = snapshot.drop(columns=["obsnme"])
    snapshot.to_csv(old / "hds_target_values.csv", index=False)

    review = IesResults(old, case_name="grouped")
    panels = _panels(review.plot_vs_obs())
    base = "oname:hds_otype:lst_usecol:{}"
    titles = (
        {"hds: W2": "w2", "hds: W1": "w1", "hds: W10": "w10"} if era == "before_358"
        else {base.format(c): c for c in ("w2", "w1", "w10")}
    )
    assert set(panels) == set(titles)
    measured = {titles[title]: values for title, (values, _n, _len) in panels.items()}
    assert measured == {"w2": [91.0, 92.0], "w1": [93.0, 94.0], "w10": [95.0, 96.0, 97.0]}
    assert set(review.phi_contributions().index) == {base.format(c) for c in ("w2", "w1", "w10")}
    assert np.isfinite(review.phi_contributions()).all()
