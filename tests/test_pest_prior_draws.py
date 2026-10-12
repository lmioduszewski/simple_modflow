"""The prior PESTPP-IES is handed, and ``cal.draw_prior()``, on the grids myflopy builds.

Ledger 385: ``run_ies()``/``prior()`` draw the prior themselves when a ``grid`` or
pilot-point parameter carries ``correlation=``. They died before PESTPP-IES started
in two layouts: a correlated ``grid`` parameter beside correlated pilot points (the
grid draw covers every parameter of ``cal.pst``, pilot points included, so joining the
pilot draw onto it overlapped on every pilot column and pandas refused the join), and
correlated pilot points beside any other parameter with no correlated ``grid`` one
(the ensemble held the pilot points alone, and PESTPP-IES refuses an ensemble missing
an adjustable parameter). The ensemble is now the whole-control-file draw with each
pilot set's columns replaced by its correlated draw.

Ledger 384: ``draw_prior()`` defaulted to ``use_specsim=True``, and pyEMU's
spectral simulation needs a regular structured grid; on a DISV model -- the toolkit's
primary grid -- every default call raised ``AttributeError`` from inside pyEMU. The
default is ``False`` (pyEMU's own, and the draw ``run_ies``/``prior`` make), and
``True`` off a regular structured grid is refused naming it.

The model is the canonical calibration demo on ``CanonicalModelConfig.testing()``
(21 x 21 DISV, 4 layers, four head wells).
"""

from __future__ import annotations

import types
import warnings
from pathlib import Path

import flopy
import numpy as np
import pandas as pd
import pytest

from myflopy.modflow.mf6.canonical_calibration import build_canonical_calibration_demo
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
from myflopy.modflow.mf6.pest.project import PestProject

pyemu = pytest.importorskip("pyemu")


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build_canonical_calibration_demo(
            tmp_path_factory.mktemp("prior_draws") / "model",
            config=CanonicalModelConfig.testing(), n_head_wells=4,
        )


def _calibration(demo, name: str, *, k33: str):
    """``k`` as correlated pilot points, and ``k33`` as ``grid`` (correlated) or ``constant``."""

    cal = demo.model.pest(name, start_datetime="2024-01-01")
    cal.parameterize("k", style="pilotpoints", pp_spacing=500.0, correlation=3000.0,
                     bounds=(0.1, 10.0), physical=(1e-3, 300.0))
    if k33 == "grid":
        cal.parameterize("k33", style="grid", correlation=3000.0, layers=[0])
    else:
        cal.parameterize("k33", style="constant")
    cal.observe(demo.head_targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build(f"{name}.pst", noptmax=0)
    return cal


def _prior_csv(cal) -> pd.DataFrame:
    return pd.read_csv(cal.template_workspace / cal.pst.pestpp_options["ies_par_en"], index_col=0)


@pytest.mark.parametrize("k33", ["grid", "constant"])
def test_the_prior_holds_every_parameter_with_the_correlated_pilot_draw(demo, monkeypatch, k33):
    cal = _calibration(demo, f"whole_{k33}", k33=k33)
    pilot_draws = []
    draw = PestProject._draw_pilot_point_prior

    def _record(self, spec, reals):
        drawn = draw(self, spec, reals)
        pilot_draws.append(drawn)
        return drawn

    monkeypatch.setattr(PestProject, "_draw_pilot_point_prior", _record)
    cal._inject_geostatistical_prior(6)

    ensemble = _prior_csv(cal)
    assert ensemble.shape[0] == 6
    assert sorted(ensemble.columns) == sorted(cal.pst.par_names)
    # The pilot columns are the CORRELATED draw, not the whole-file draw's own.
    (pilot,) = pilot_draws
    assert len(pilot.columns) == 100
    np.testing.assert_allclose(ensemble[pilot.columns].to_numpy(), pilot.to_numpy(), rtol=1e-12)
    others = [name for name in cal.pst.par_names if name not in set(pilot.columns)]
    assert others and not ensemble[others].isna().any().any()


@pytest.mark.slow
@pytest.mark.parametrize("k33", ["grid", "constant"])
def test_pestpp_ies_runs_the_prior_it_is_handed(demo, k33):
    """A real PESTPP-IES prior run (NOPTMAX=-1): it starts, and its realizations are the
    rows of the ensemble the calibration drew, held to each parameter's bounds as
    PESTPP-IES holds them (PESTPP-IES writes six digits)."""

    cal = _calibration(demo, f"ran_{k33}", k33=k33)
    results = cal.prior(reals=3)
    ran = pd.read_csv(Path(results.workspace) / f"ran_{k33}.0.par.csv", index_col=0)
    ran.columns = [column.lower() for column in ran.columns]
    bounds = cal.pst.parameter_data
    drawn = _prior_csv(cal)
    drawn = drawn.clip(lower=bounds.loc[drawn.columns, "parlbnd"],
                       upper=bounds.loc[drawn.columns, "parubnd"], axis=1)
    for realization in ("0", "1"):
        row = ran.loc[[label for label in ran.index if str(label) == realization]].iloc[0]
        np.testing.assert_allclose(row[list(drawn.columns)].to_numpy(dtype=float),
                                   drawn.loc[int(realization)].to_numpy(dtype=float), rtol=1e-5)
    obs = pd.read_csv(Path(results.workspace) / f"ran_{k33}.0.obs.csv", index_col=0)
    assert len(obs) == 3 and np.isfinite(obs.to_numpy(dtype=float)).all()


def test_draw_prior_draws_by_default_on_a_disv_model_and_refuses_specsim_there(demo):
    cal = demo.model.pest("specsim", start_datetime="2024-01-01")
    cal.parameterize("k", style="grid", correlation=3000.0, layers=[0])
    cal.observe(demo.head_targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("specsim.pst", noptmax=0)
    assert isinstance(cal.pf.spatial_reference, flopy.discretization.VertexGrid)

    drawn = cal.draw_prior(5)
    frame = drawn._df if hasattr(drawn, "_df") else drawn
    assert frame.shape == (5, cal.pst.npar)
    with pytest.raises(ValueError, match=r"a VertexGrid, which has no rows or columns.*use_specsim=False"):
        cal.draw_prior(5, use_specsim=True)


def _stub(reference):
    """Just enough of a built calibration for ``draw_prior``: ``pf`` records its draw."""

    calls = []
    pf = types.SimpleNamespace(
        spatial_reference=reference,
        draw=lambda num_reals, use_specsim, sigma_range, scale_offset: calls.append(
            (num_reals, use_specsim, sigma_range, scale_offset)) or "drawn",
    )
    stub = types.SimpleNamespace(pf=pf, pyemu=None)
    stub._refuse_specsim_off_a_regular_grid = (
        lambda: PestProject._refuse_specsim_off_a_regular_grid(stub))
    return stub, calls


def _structured(delr):
    return flopy.discretization.StructuredGrid(
        delc=np.full(4, 10.0), delr=np.asarray(delr, dtype=float),
        top=np.full((4, len(delr)), 10.0), botm=np.zeros((1, 4, len(delr))),
    )


def test_specsim_is_passed_to_pyemu_only_on_a_regular_structured_grid():
    stub, calls = _stub(_structured([10.0] * 5))
    assert PestProject.draw_prior(stub, 7, use_specsim=True) == "drawn"
    assert PestProject.draw_prior(stub, 7) == "drawn"
    # Both at the prior's one sigma range (ledger 424).
    assert calls == [(7, True, 6, True), (7, False, 6, True)]

    stub, calls = _stub(_structured([10.0, 10.0, 20.0, 10.0, 10.0]))
    with pytest.raises(ValueError, match="not all one width.*use_specsim=False"):
        PestProject.draw_prior(stub, 7, use_specsim=True)
    stub, calls = _stub(None)
    with pytest.raises(ValueError, match="a NoneType, which has no rows or columns"):
        PestProject.draw_prior(stub, 7, use_specsim=True)
    assert calls == []
