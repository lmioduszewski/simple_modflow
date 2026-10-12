"""Pilot points stack, interpolate in log space from their nearest points, and honour zones.

Ledgers 525 (one command per array file, multiplying what pyEMU wrote),
526 (log10 IDW from the ``pp_nearest`` nearest points of the set, weights
precomputed at build) and 527 (``pp_spacing`` in model length units, ``zones=``
honoured, one set per layer and zone).

Before them a pilot-point set REBUILT its layer file from the base captured at build,
after pyEMU's ``apply_list_and_array_pars``: a constant x2 under pilot points at 1.0
wrote 1.0x the base, and a second pilot-point set on a layer erased the first --
both measured with a real forward run, exit 0. ``zones=`` was accepted and ignored.
"""

from __future__ import annotations

import json
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
from myflopy.modflow.mf6.pest import forward_run
from myflopy.modflow.mf6.pest.pilot_points import (
    PILOT_POINT_MANIFEST,
    idw_factors,
)

pytest.importorskip("pyemu")

#: The canonical testing profile's cells are ~100 m wide on a 2,100 m domain.
SPACING = 600.0


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    """The canonical testing model, run once, with four head targets."""

    from myflopy.modflow.mf6.canonical_calibration import build_canonical_calibration_demo

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build_canonical_calibration_demo(
            tmp_path_factory.mktemp("ppz") / "model",
            config=CanonicalModelConfig.testing(), n_head_wells=4,
        )


def _base_k(demo, layer=0):
    return np.asarray(demo.model.gwf.npf.k.get_data(), dtype=float)[layer]


def _xc(demo):
    return np.asarray(demo.model.gwf.modelgrid.xcellcenters, dtype=float).reshape(-1)


def _west_east(demo, west=1, east=0):
    xc = _xc(demo)
    return np.where(xc < np.median(xc), west, east)


def _build(demo, tmp_path, declarations, name="ppz"):
    cal = demo.model.pest(name, workspace=tmp_path / name, start_datetime="2024-01-01")
    for args, kwargs in declarations:
        cal.parameterize(*args, **kwargs)
    cal.observe(demo.head_targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build(f"{name}.pst", noptmax=0)
    return cal


def _write_inputs(cal):
    """``Pst.write_input_files``, serially: pyEMU's opens a process pool, whose temp
    directory outlives the test (ledger 335)."""

    from pyemu.pst.pst_utils import write_to_template

    par = cal.pst.parameter_data
    values = par["parval1"] * par["scale"] + par["offset"]
    values.index = values.index.str.lower()
    for tpl, model_file in zip(cal.pst.template_files, cal.pst.input_files, strict=True):
        write_to_template(values, str(cal.template_workspace / tpl),
                          str(cal.template_workspace / model_file))


def _set_values(cal, values: dict[str, float]):
    """Write the template's model input files with each named group at a value."""

    par = cal.pst.parameter_data
    for group, value in values.items():
        rows = par["pargp"] == group
        assert rows.any(), f"no parameter group {group!r}: {sorted(par['pargp'].unique())}"
        par.loc[rows, "parval1"] = value
    _write_inputs(cal)


def _pre_commands(cal, monkeypatch):
    """Run the generated forward run's pre-model commands, in order, without MF6."""

    import pyemu

    text = (cal.template_workspace / "forward_run.py").read_text()
    main = text[text.index("def main():"):]
    calls = [line.strip() for line in main.splitlines()
             if "apply_list_and_array_pars(" in line or "apply_pilotpoints_to_array(" in line]
    assert calls, "the forward run has no pre-model commands"
    monkeypatch.chdir(cal.template_workspace)
    namespace = {"pyemu": pyemu,
                 "apply_pilotpoints_to_array": forward_run.apply_pilotpoints_to_array}
    for call in calls:
        exec(call, namespace)  # noqa: S102 -- the template's own forward-run lines
    return calls


def _written(cal, demo, layer=1):
    return np.loadtxt(cal.template_workspace / f"{demo.model.name}.npf_k_layer{layer}.txt")


# --- (a) stacking ---------------------------------------------------------------------


def test_a_constant_under_pilot_points_reaches_mf6_in_a_real_forward_run(demo, tmp_path):
    """Constant x3 with every pilot point at 1.0: the file MF6 reads is 3x the base, and
    the generated forward run runs pyEMU's apply BEFORE the pilot-point command. Before
    ledger 525 this read 1.0000000000041."""

    cal = _build(demo, tmp_path, [
        (("k",), dict(style="constant", layers=[0], name="kc", bounds=(0.1, 10.0))),
        (("k",), dict(style="pilotpoints", layers=[0], name="kp", pp_spacing=SPACING,
                      bounds=(0.1, 10.0))),
    ])
    _set_values(cal, {"kc": 3.0})
    text = (cal.template_workspace / "forward_run.py").read_text()
    main = text[text.index("def main():"):]
    assert main.index("apply_list_and_array_pars(") < main.index("apply_pilotpoints_to_array(")
    result = subprocess.run([sys.executable, "forward_run.py"], cwd=cal.template_workspace,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout[-1500:] + result.stderr[-1500:]
    assert np.allclose(_written(cal, demo) / _base_k(demo), 3.0, rtol=1e-6)


def test_two_pilot_point_sets_on_one_layer_multiply(demo, tmp_path, monkeypatch):
    """Two calls on layer 0, every point of one at 2.0 and of the other at 1.5: 3.0x. The
    second set used to erase the first (ratio 1.0)."""

    cal = _build(demo, tmp_path, [
        (("k",), dict(style="pilotpoints", layers=[0], name="ka", pp_spacing=SPACING,
                      bounds=(0.1, 10.0))),
        (("k",), dict(style="pilotpoints", layers=[0], name="kb", pp_spacing=SPACING / 2,
                      bounds=(0.1, 10.0))),
    ])
    _set_values(cal, {"kal0": 2.0, "kbl0": 1.5})
    calls = _pre_commands(cal, monkeypatch)
    assert sum("apply_pilotpoints_to_array(" in call for call in calls) == 1, (
        "one command per array file")
    assert np.allclose(_written(cal, demo) / _base_k(demo), 3.0, rtol=1e-8)


@pytest.mark.parametrize("with_pyemu", [True, False], ids=["pyemu-also-writes", "captured-base"])
def test_a_reused_worker_directory_writes_the_same_file_every_run(demo, tmp_path, monkeypatch,
                                                                  with_pyemu):
    """A worker runs realization after realization in one directory. With pyEMU also
    parameterizing the file the pilot step multiplies what pyEMU rewrote from org/;
    without, it starts from the captured base -- either way, twice is the same as once."""

    declarations = [(("k",), dict(style="pilotpoints", layers=[0], name="kp",
                                  pp_spacing=SPACING, bounds=(0.1, 10.0)))]
    if with_pyemu:
        declarations.insert(0, (("k",), dict(style="constant", layers=[0], name="kc",
                                             bounds=(0.1, 10.0))))
    cal = _build(demo, tmp_path, declarations)
    manifest = json.loads((cal.template_workspace / PILOT_POINT_MANIFEST).read_text())
    assert [entry["base_from_file"] for entry in manifest["files"]] == [with_pyemu]
    values = {"kpl0": 1.7, **({"kc": 2.0} if with_pyemu else {})}
    _set_values(cal, values)
    _pre_commands(cal, monkeypatch)
    once = _written(cal, demo)
    _pre_commands(cal, monkeypatch)
    assert np.array_equal(_written(cal, demo), once)
    assert np.allclose(once / _base_k(demo), 1.7 * (2.0 if with_pyemu else 1.0), rtol=1e-8)


# --- (b) zones -------------------------------------------------------------------------


def test_zones_are_honoured_and_cells_outside_keep_their_value(demo, tmp_path, monkeypatch):
    """West half zone 1, east half 0, every pilot value 2: the west is exactly 2x, the
    east untouched. ``zones=`` was ignored: 221 east cells moved."""

    zones = _west_east(demo)
    cal = _build(demo, tmp_path, [(("k",), dict(
        style="pilotpoints", layers=[0], name="kp", pp_spacing=SPACING, zones=zones,
        bounds=(0.1, 10.0)))])
    points = pd.read_csv(cal.template_workspace / "kpl0z1_pp.points.csv")
    assert (points["x"] < np.median(_xc(demo))).all(), "a point was laid outside its zone"
    _set_values(cal, {"kpl0z1": 2.0})
    _pre_commands(cal, monkeypatch)
    ratio = _written(cal, demo) / _base_k(demo)
    assert np.allclose(ratio[zones == 1], 2.0, rtol=1e-8)
    assert np.allclose(ratio[zones == 0], 1.0, rtol=1e-9, atol=0.0)


def test_nothing_blends_across_a_zone_edge(demo, tmp_path, monkeypatch):
    """Zones 1 and 2 at 2.0 and 0.5: each zone exactly its own value, at its edge too."""

    zones = _west_east(demo, west=1, east=2)
    cal = _build(demo, tmp_path, [(("k",), dict(
        style="pilotpoints", layers=[0], name="kp", pp_spacing=SPACING, zones=zones,
        bounds=(0.1, 10.0)))])
    assert set(cal.pst.parameter_data["pargp"]) >= {"kpl0z1", "kpl0z2"}
    _set_values(cal, {"kpl0z1": 2.0, "kpl0z2": 0.5})
    _pre_commands(cal, monkeypatch)
    ratio = _written(cal, demo) / _base_k(demo)
    assert np.allclose(ratio[zones == 1], 2.0, rtol=1e-8)
    assert np.allclose(ratio[zones == 2], 0.5, rtol=1e-8)


def test_per_zone_settings_are_separate_calls_with_masked_zones(demo, tmp_path, monkeypatch):
    """The one spelling for different settings per zone: two calls, each a masked zones=
    and its own name; disjoint sets on one layer each reach their own cells."""

    zones = _west_east(demo, west=1, east=2)
    cal = _build(demo, tmp_path, [
        (("k",), dict(style="pilotpoints", layers=[0], name="kw", pp_spacing=SPACING,
                      zones=np.where(zones == 1, 1, 0), bounds=(0.1, 10.0))),
        (("k",), dict(style="pilotpoints", layers=[0], name="ke", pp_spacing=SPACING / 2,
                      zones=np.where(zones == 2, 2, 0), bounds=(0.2, 5.0))),
    ])
    _set_values(cal, {"kwl0z1": 3.0, "kel0z2": 0.25})
    _pre_commands(cal, monkeypatch)
    ratio = _written(cal, demo) / _base_k(demo)
    assert np.allclose(ratio[zones == 1], 3.0, rtol=1e-8)
    assert np.allclose(ratio[zones == 2], 0.25, rtol=1e-8)


def test_set_names_carry_layer_and_zone_and_unzoned_names_are_unchanged(demo, tmp_path):
    """Unzoned: ``<name>l<layer>pp0000``; zoned: ``<name>l<layer>z<zone>pp0000``. An
    upper-case name= is lower-cased as PEST reads it -- it died with a KeyError."""

    zones = _west_east(demo, west=3, east=0)
    cal = _build(demo, tmp_path, [
        (("k",), dict(style="pilotpoints", layers=[0, 1], name="kU", pp_spacing=SPACING,
                      bounds=(0.1, 10.0))),
        (("k",), dict(style="pilotpoints", layers=[1], name="KZ", pp_spacing=SPACING,
                      zones=zones, bounds=(0.1, 10.0))),
    ])
    par = cal.pst.parameter_data
    assert set(par["pargp"]) == {"kul0", "kul1", "kzl1z3"}
    assert "kul0pp0000" in par.index and "kzl1z3pp0000" in par.index
    sets = cal._native_parameter_specs[1].pilot_point_sets
    assert [(item["name"], item["layer"], item["zone"]) for item in sets] == [("kzl1z3", 1, 3)]
    assert sets[0]["points"] == int((par["pargp"] == "kzl1z3").sum())
    assert "kzl1z3" in str(cal.settings())


def test_two_calls_on_separate_layers_keep_both_sets_under_one_name(demo, tmp_path):
    """The natural per-layer spelling, default name both times: both sets are in the
    control file. The second call used to drop the first set's parameters while its
    forward-run command kept applying them at 1.0."""

    cal = _build(demo, tmp_path, [
        (("k",), dict(style="pilotpoints", layers=[0], pp_spacing=SPACING, bounds=(0.1, 10.0))),
        (("k",), dict(style="pilotpoints", layers=[1], pp_spacing=SPACING / 2,
                      bounds=(0.1, 10.0))),
    ])
    groups = cal.pst.parameter_data["pargp"].value_counts()
    assert {"kl0", "kl1"} <= set(groups.index)
    assert groups["kl1"] > groups["kl0"]


def test_two_sets_under_one_name_on_one_layer_are_refused(demo, tmp_path):
    with pytest.raises(ValueError, match="Two pilot-point sets are named 'kl0'"):
        _build(demo, tmp_path, [
            (("k",), dict(style="pilotpoints", layers=[0], pp_spacing=SPACING)),
            (("k",), dict(style="pilotpoints", layers=[0], pp_spacing=SPACING / 2)),
        ])


# --- (c) interpolation and its factors --------------------------------------------------


def test_interpolation_is_in_log_space_from_the_nearest_points():
    """A cell midway between a point at 10 and one at 0.1 reads 1.0 (their geometric
    mean), where linear IDW read 5.05; a far point at 100 is not among the 2 nearest
    and moves nothing."""

    index, weight = idw_factors([0.5, 0.0], [0.0, 0.0], [0.0, 1.0, 50.0], [0.0, 0.0, 0.0],
                                nearest=2)
    values = np.array([10.0, 0.1, 100.0])
    field = 10.0 ** np.sum(weight * np.log10(values)[index], axis=1)
    assert field[0] == pytest.approx(1.0)
    assert field[1] == pytest.approx(10.0), "a cell at a point takes that point's value"
    assert not (index == 2).any()


def _manifest(tmp_path, values, *, interpolation, base_from_file=False, limits=(None, None)):
    """A one-set manifest over three cells at x = 0, 0.5 and 1 between points at 0 and 1."""

    pd.DataFrame({"parnme": ["sp0000", "sp0001"], "value": values}).to_csv(
        tmp_path / "s_pp.csv", index=False)
    pd.DataFrame({"parnme": ["sp0000", "sp0001"], "x": [0.0, 1.0], "y": [0.0, 0.0],
                  "layer": 0, "zone": 0}).to_csv(tmp_path / "s_pp.points.csv", index=False)
    cell_x = np.array([0.0, 0.5, 1.0])
    index, weight = idw_factors(cell_x, np.zeros(3), [0.0, 1.0], [0.0, 0.0])
    np.savez(tmp_path / "s_pp.factors.npz", cells=np.arange(3), x=cell_x, y=np.zeros(3),
             index=index, weight=weight)
    np.savetxt(tmp_path / "a.txt", np.full((3, 1), 2.0))
    np.savetxt(tmp_path / "a.txt.pp_base", np.full((3, 1), 2.0))
    (tmp_path / PILOT_POINT_MANIFEST).write_text(json.dumps({
        "version": 1,
        "files": [{"array_file": "a.txt", "base_from_file": base_from_file,
                   "base_file": None if base_from_file else "a.txt.pp_base", "sets": ["s"]}],
        "sets": [{"name": "s", "array_file": "a.txt", "parameter_csv": "s_pp.csv",
                  "points_csv": "s_pp.points.csv", "factors": "s_pp.factors.npz",
                  "interpolation": interpolation, "lower_limit": limits[0],
                  "upper_limit": limits[1]}],
    }))


@pytest.mark.parametrize("interpolation,midway", [("log", 2.0), ("linear", 10.1)])
def test_the_forward_run_interpolates_a_log_parameter_in_log_space(tmp_path, monkeypatch,
                                                                   interpolation, midway):
    """Points at 10 and 0.1, base 2: midway reads 2 x 1.0 (the geometric mean) for a
    log-transformed parameter -- linear IDW read 2 x 5.05 -- and the cells AT the points
    take their values."""

    _manifest(tmp_path, [10.0, 0.1], interpolation=interpolation)
    monkeypatch.chdir(tmp_path)
    forward_run.apply_pilotpoints_to_array("a.txt")
    assert np.loadtxt("a.txt") == pytest.approx([20.0, midway, 0.2])


def test_the_forward_run_clamps_and_starts_from_the_file_pyemu_wrote(tmp_path, monkeypatch):
    _manifest(tmp_path, [10.0, 10.0], interpolation="log", base_from_file=True,
              limits=(0.0, 30.0))
    np.savetxt(tmp_path / "a.txt", np.array([[1.0], [2.0], [4.0]]))
    monkeypatch.chdir(tmp_path)
    forward_run.apply_pilotpoints_to_array("a.txt")
    assert np.loadtxt("a.txt") == pytest.approx([10.0, 20.0, 30.0])


def test_the_forward_run_applies_the_weights_stored_at_build(demo, tmp_path, monkeypatch):
    """The interpolation is a gather through the factors the build wrote, not distances
    recomputed per cell every run: rewriting them changes the field."""

    cal = _build(demo, tmp_path, [(("k",), dict(
        style="pilotpoints", layers=[0], name="kp", pp_spacing=SPACING, bounds=(0.1, 10.0)))])
    stored = dict(np.load(cal.template_workspace / "kpl0_pp.factors.npz"))
    assert stored["index"].shape == (stored["cells"].size, 8)
    assert np.allclose(stored["weight"].sum(axis=1), 1.0)
    par = cal.pst.parameter_data
    par.loc[par["pargp"] == "kpl0", "parval1"] = 1.0
    par.loc["kpl0pp0000", "parval1"] = 4.0
    _write_inputs(cal)
    stored["index"][:] = 0
    stored["weight"][:] = 0.0
    stored["weight"][:, 0] = 1.0
    np.savez(cal.template_workspace / "kpl0_pp.factors.npz", **stored)
    _pre_commands(cal, monkeypatch)
    assert np.allclose(_written(cal, demo) / _base_k(demo), 4.0, rtol=1e-8)


def test_a_block_of_points_reaches_its_value_between_them(demo, tmp_path, monkeypatch):
    """Points within 500 m of the domain's centre at 5, the rest at 1: a cell inside the
    block reads 5 between its points. Global IDW read 2.3 inside such a block."""

    cal = _build(demo, tmp_path, [(("k",), dict(
        style="pilotpoints", layers=[0], name="kp", pp_spacing=200.0, bounds=(0.1, 10.0)))])
    points = pd.read_csv(cal.template_workspace / "kpl0_pp.points.csv").set_index("parnme")
    xc, yc = _xc(demo), np.asarray(demo.model.gwf.modelgrid.ycellcenters, float).reshape(-1)
    centre = np.array([np.median(xc), np.median(yc)])
    inside = np.hypot(points["x"] - centre[0], points["y"] - centre[1]) < 500.0
    par = cal.pst.parameter_data
    par.loc[points.index, "parval1"] = np.where(inside, 5.0, 1.0)
    _write_inputs(cal)
    _pre_commands(cal, monkeypatch)
    ratio = _written(cal, demo) / _base_k(demo)
    core = np.hypot(xc - centre[0], yc - centre[1]) < 150.0
    assert core.any()
    assert np.allclose(ratio[core], 5.0, rtol=1e-6)


def test_the_net_spacing_is_in_model_length_units(demo, tmp_path):
    """Neighbouring points of a net are about its spacing apart, and halving the spacing
    roughly quadruples the count."""

    counts = {}
    for spacing in (SPACING, SPACING / 2):
        cal = _build(demo, tmp_path / str(int(spacing)), [(("k",), dict(
            style="pilotpoints", layers=[0], name="kp", pp_spacing=spacing))])
        points = pd.read_csv(cal.template_workspace / "kpl0_pp.points.csv")
        xy = points[["x", "y"]].to_numpy()
        gaps = np.hypot(*(xy[:, None, :] - xy[None, :, :]).transpose(2, 0, 1))
        np.fill_diagonal(gaps, np.inf)
        nearest = gaps.min(axis=1)
        # Each point is the cell centre nearest its square's middle, ~100 m cells.
        assert 0.6 * spacing < np.median(nearest) < 1.2 * spacing
        counts[spacing] = len(points)
    assert 3.0 < counts[SPACING / 2] / counts[SPACING] < 5.0


def test_a_covered_cell_is_clamped_to_the_tightest_covering_limits(demo, tmp_path, monkeypatch):
    zones = _west_east(demo, west=1, east=0)
    base = _base_k(demo)
    cap = float(np.median(base))
    cal = _build(demo, tmp_path, [
        (("k",), dict(style="pilotpoints", layers=[0], name="ka", pp_spacing=SPACING,
                      bounds=(0.1, 10.0), physical=(1e-6, 1e6))),
        (("k",), dict(style="pilotpoints", layers=[0], name="kb", pp_spacing=SPACING,
                      zones=zones, bounds=(0.1, 10.0), physical=(1e-6, cap))),
    ])
    _set_values(cal, {"kal0": 5.0, "kbl0z1": 1.0})
    _pre_commands(cal, monkeypatch)
    written = _written(cal, demo)
    assert written[zones == 1].max() <= cap * (1 + 1e-9)
    assert np.allclose(written[zones == 0], 5.0 * base[zones == 0], rtol=1e-8)


# --- (d) refusals ----------------------------------------------------------------------


def _points(demo, cells):
    grid = demo.model.gwf.modelgrid
    return np.c_[np.asarray(grid.xcellcenters).reshape(-1)[cells],
                 np.asarray(grid.ycellcenters).reshape(-1)[cells]]


def test_explicit_points_take_their_cells_zone_and_are_refused_outside_one(demo, tmp_path):
    zones = _west_east(demo)
    west, east = int(np.flatnonzero(zones == 1)[0]), int(np.flatnonzero(zones == 0)[0])
    with pytest.raises(ValueError, match="whose zone is 0"):
        _build(demo, tmp_path / "a", [(("k",), dict(
            style="pilotpoints", layers=[0], pp_points=_points(demo, [west, east]),
            zones=zones))])
    cal = _build(demo, tmp_path / "b", [(("k",), dict(
        style="pilotpoints", layers=[0], pp_points=_points(demo, [west]), zones=zones))])
    assert list(cal.pst.parameter_data.index) == ["kl0z1pp0000"]


def test_a_zone_left_without_an_explicit_point_is_refused(demo, tmp_path):
    zones = _west_east(demo, west=1, east=2)
    west = int(np.flatnonzero(zones == 1)[0])
    with pytest.raises(ValueError, match="Zone 2 .* holds none of the explicit pilot points"):
        _build(demo, tmp_path, [(("k",), dict(
            style="pilotpoints", layers=[0], pp_points=_points(demo, [west]), zones=zones))])


def test_zones_labelling_nothing_are_refused(demo, tmp_path):
    with pytest.raises(ValueError, match="parameterizes nothing"):
        _build(demo, tmp_path, [(("k",), dict(
            style="pilotpoints", layers=[0], pp_spacing=SPACING,
            zones=np.zeros_like(_west_east(demo))))])


def test_an_explicit_point_off_the_grid_or_given_twice_is_refused(demo, tmp_path):
    xc = _xc(demo)
    with pytest.raises(ValueError, match="lie in no cell of the grid"):
        _build(demo, tmp_path / "a", [(("k",), dict(
            style="pilotpoints", pp_points=[(xc.max() + 1e5, 0.0)]))])
    twice = np.vstack([_points(demo, [3]), _points(demo, [3])])
    with pytest.raises(ValueError, match="given more than once"):
        _build(demo, tmp_path / "b", [(("k",), dict(style="pilotpoints", pp_points=twice))])


def test_physical_limits_that_share_no_value_on_one_cell_are_refused(demo, tmp_path):
    with pytest.raises(ValueError, match="share no value"):
        _build(demo, tmp_path, [
            (("k",), dict(style="pilotpoints", layers=[0], name="ka", pp_spacing=SPACING,
                          physical=(10.0, 20.0))),
            (("k",), dict(style="pilotpoints", layers=[0], name="kb", pp_spacing=SPACING,
                          physical=(1.0, 5.0))),
        ])


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"style": "pilotpoints"}, "give exactly one"),
        ({"style": "pilotpoints", "pp_spacing": 100.0, "pp_points": [(0.0, 0.0)]},
         "give exactly one"),
        ({"style": "pilotpoints", "pp_spacing": 0.0}, "finite length above zero"),
        ({"style": "pilotpoints", "pp_spacing": float("nan")}, "finite length above zero"),
        ({"style": "pilotpoints", "pp_spacing": 100.0, "pp_nearest": 0}, "1 or more"),
        ({"style": "pilotpoints", "pp_spacing": 100.0, "pp_nearest": 2.5}, "1 or more"),
        ({"style": "grid", "pp_spacing": 100.0}, "pp_spacing place and interpolate"),
        ({"style": "constant", "pp_nearest": 4}, "pp_nearest place and interpolate"),
    ],
)
def test_pilot_point_arguments_that_cannot_act_are_refused_at_the_call(kwargs, match):
    from myflopy.modflow.mf6.pest.native_parameters import NativeParameterSpec

    with pytest.raises(ValueError, match=match):
        NativeParameterSpec(target="k", **kwargs)


def test_pp_space_is_gone_with_no_alias():
    from myflopy.modflow.mf6.pest.native_parameters import NativeParameterSpec
    from myflopy.modflow.mf6.pest.project import PestProject

    with pytest.raises(TypeError):
        NativeParameterSpec(target="k", style="pilotpoints", pp_space=5)
    import inspect

    assert "pp_space" not in inspect.signature(PestProject.parameterize).parameters


def test_the_manifest_names_each_set_for_a_write_back(demo, tmp_path):
    """Everything a write-back needs to rebuild a set's field is in the template."""

    zones = _west_east(demo, west=2, east=0)
    cal = _build(demo, tmp_path, [(("k",), dict(
        style="pilotpoints", layers=[0], name="kp", pp_spacing=SPACING, zones=zones,
        physical=(1e-3, 300.0)))])
    manifest = json.loads((cal.template_workspace / PILOT_POINT_MANIFEST).read_text())
    (record,) = manifest["sets"]
    assert record["name"] == "kpl0z2" and record["layer"] == 0 and record["zone"] == 2
    assert record["interpolation"] == "log" and record["nearest"] == 8
    assert record["lower_limit"] == 1e-3 and record["upper_limit"] == 300.0
    for key in ("parameter_csv", "points_csv", "factors", "template"):
        assert Path(cal.template_workspace / record[key]).is_file(), key
    assert record["array_file"] == f"{demo.model.name}.npf_k_layer1.txt"
