"""``uzf.vks`` as a calibration target (plan §5.8 item 1).

UZF is the first list package whose external file does NOT open with the cell
identity: it numbers its own cells, so a row reads ``iuzno layer icell2d ...``.
That one difference is the whole item, and it fails silently — pyEMU reads the
LAYER as the cell number, places every parameter at one cell's coordinates, and
the control file, the forward run and MODFLOW are all perfectly happy.
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest
from pyemu.pst.pst_utils import write_to_template

from myflopy.modflow.mf6.canonical_calibration import build_canonical_calibration_demo
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
from myflopy.modflow.mf6.pest.native_parameters import resolve_target


def _uzf_rows(path):
    """The external UZF packagedata table as a float array."""

    return np.array(
        [line.split() for line in path.read_text().splitlines() if line.strip()],
        dtype=float,
    )


def test_uzf_vks_declares_the_columns_that_make_it_different():
    """Two numbers carry this target, and both were measured against a written
    model rather than read off the MF6 manual: ``vks`` is column 6, and the cell
    identity is at columns (1, 2) because ``iuzno`` occupies column 0."""

    recipe = resolve_target("uzf.vks")
    assert (recipe.family, recipe.use_col) == ("list", 6)
    assert recipe.index_cols == (1, 2)
    assert resolve_target("vks") is recipe
    assert resolve_target("uzf") is recipe

    # Every other list target keeps the ordinary (layer, cell) opening.
    for target in ("recharge", "chd", "ghb.cond", "drn.elev", "wel"):
        assert resolve_target(target).index_cols == (0, 1)


@pytest.mark.slow
def test_uzf_parameters_land_on_their_own_cells(tmp_path):
    """The silent-wrongness case. With the default (0, 1), pyEMU takes the LAYER
    column as the cell number, so all 189 UZF parameters get cell 0's
    coordinates. Nothing fails at build: the .pst is valid, the multipliers are
    applied to the right rows, and MODFLOW runs. Only the geostatistics are
    garbage, and that surfaces much later — as ``error inverting cov`` from the
    prior draw, an error naming parameters but not the cause.
    """

    demo = build_canonical_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing(), n_head_wells=4
    )
    model = demo.model
    cal = model.pest("uzfcal", start_datetime="2024-01-01")
    cal.parameterize("uzf.vks", style="grid", correlation=400.0,
                     bounds=(0.2, 5.0), physical=(1e-4, 10.0))
    cal.observe(demo.head_targets)
    cal.build("uzfcal.pst", noptmax=0)

    frame = cal._native_parameter_frames["uzfvks"]
    assert frame["x"].nunique() == len(frame), (
        "UZF parameters share coordinates -- they were geolocated by layer"
    )

    # ...and the coordinates are the RIGHT cells, not merely distinct ones.
    uzf_cells = [
        int(row[1][1]) for row in model.gwf.get_package("uzf").packagedata.get_data()
    ]
    centers = np.asarray(model.gwf.modelgrid.xcellcenters, dtype=float).reshape(-1)
    assert np.allclose(
        sorted(frame["x"].astype(float)), sorted(centers[uzf_cells])
    )

    # The prior draw is what actually blew up on the wrong coordinates.
    cal.pf.build_prior(fmt="none")


@pytest.mark.slow
def test_the_uzf_multiplier_scales_vks_and_nothing_else(tmp_path):
    """`use_col` off by one would scale `surfdep` or `thtr` instead — a model
    that still runs and still calibrates, against the wrong property."""

    demo = build_canonical_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing(), n_head_wells=4
    )
    model = demo.model
    cal = model.pest("uzfmult", start_datetime="2024-01-01")
    cal.parameterize("uzf.vks", style="constant", bounds=(0.2, 5.0),
                     physical=(1e-4, 10.0))
    cal.observe(demo.head_targets)
    pst = cal.build("uzfmult.pst", noptmax=0)

    template = cal.template_workspace
    packagedata = template / f"{model.name}.uzf_packagedata.txt"

    subprocess.run([sys.executable, "forward_run.py"], cwd=template,
                   check=True, capture_output=True, text=True)
    baseline = _uzf_rows(packagedata)
    assert np.allclose(np.unique(baseline[:, 6]), [0.25]), (
        "unit multipliers must reproduce the model's own vks"
    )

    parameters = pst.parameter_data
    selected = parameters.index[parameters["pargp"].str.startswith("uzfvks")]
    parameters.loc[selected, "parval1"] = 2.0
    for tpl, inp in zip(pst.template_files, pst.input_files, strict=False):
        write_to_template(parameters["parval1"], str(template / tpl), str(template / inp))
    subprocess.run([sys.executable, "forward_run.py"], cwd=template,
                   check=True, capture_output=True, text=True)

    doubled = _uzf_rows(packagedata)
    assert np.allclose(np.unique(doubled[:, 6]), [0.5]), "vks did not scale"
    untouched = [0, 1, 2, 3, 4, 5, 7, 8, 9, 10]
    assert np.allclose(doubled[:, untouched], baseline[:, untouched]), (
        "a column other than vks moved -- use_col is wrong"
    )


@pytest.mark.slow
def test_boundnames_do_not_shift_the_vks_column(tmp_path):
    """`use_col=6` is a positional promise about a file MODFLOW writes, so it is
    only safe if optional fields cannot move it. `boundname` is declared LAST in
    the MF6 dfn; this asserts that empirically rather than trusting the ordering,
    because a shifted column would silently calibrate `thtr` instead."""

    from myflopy.modflow.mf6.canonical_example import build_canonical_model

    model = build_canonical_model(
        tmp_path / "model", config=CanonicalModelConfig.testing()
    )
    uzf = model.gwf.get_package("uzf")
    rows = [list(row) for row in uzf.packagedata.get_data()]
    uzf.boundnames = True
    uzf.packagedata = [
        tuple(row[:-1]) + (f"uz{index:03d}",) for index, row in enumerate(rows)
    ]
    model.sim.set_all_data_external()
    model.sim.write_simulation(silent=True)

    path = next((tmp_path / "model").rglob(f"{model.name}.uzf_packagedata.txt"))
    first = path.read_text().splitlines()[0].split()

    assert len(first) == 12, "boundname should append one column"
    assert float(first[6]) == pytest.approx(0.25), (
        "vks moved when boundnames were enabled; use_col=6 is unsafe"
    )
