from __future__ import annotations

import flopy
import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy import EVTBuilder, ModelContext, ModelSpec, PackageSpec, SimulationSpec
from myflopy.geopackage import CellSurfaceOffset
from myflopy.modflow.mf6.areal import _ArealBuilder


def _surfaces(top=(100.0, 90.0, 80.0)) -> pd.DataFrame:
    return pd.DataFrame({"top": list(top), "bottom": [0.0] * len(top)})


def _context() -> ModelContext:
    # Same idomain as the RCH builder test: top-active cells are (0,0),(1,1),(0,2).
    return ModelContext(domain=np.array([[1, 0, 1], [1, 1, 0]]), surfaces=_surfaces())


def test_evt_builder_shares_the_areal_base():
    assert issubclass(EVTBuilder, _ArealBuilder)


def test_evt_builder_selects_top_active_cells_and_resolves_model_top_surface():
    builder = EVTBuilder(context=_context(), nper=2, rate=2.0e-3, depth=2.5)

    assert builder.selected_cells == ((0, 0), (1, 1), (0, 2))
    # surface defaults to the model top (land surface) per column: 100/90/80.
    assert builder.stress_period_data == {
        0: [[(0, 0), 100.0, 2.0e-3, 2.5], [(1, 1), 90.0, 2.0e-3, 2.5], [(0, 2), 80.0, 2.0e-3, 2.5]],
        1: [[(0, 0), 100.0, 2.0e-3, 2.5], [(1, 1), 90.0, 2.0e-3, 2.5], [(0, 2), 80.0, 2.0e-3, 2.5]],
    }


def test_evt_builder_broadcasts_rate_and_depth_across_shapes():
    period_scalars = EVTBuilder(
        context=_context(), nper=2, rate={0: 1.0e-3, 1: 2.0e-3}, depth=2.0
    )
    by_cell = EVTBuilder(
        context=_context(),
        nper=1,
        rate={(0, 0): 1.0e-3, (1, 1): 2.0e-3, (0, 2): 3.0e-3},
        depth={(0, 0): 1.0, (1, 1): 2.0, (0, 2): 3.0},
    )
    per_cell_seq = EVTBuilder(
        context=_context(), nper=1, rate=[1.0e-3, 2.0e-3, 3.0e-3], depth=1.5
    )

    # rate broadcasts exactly like RCH recharge (index 2 of the EVT record).
    assert [record[2] for record in period_scalars.stress_period_data[1]] == [2.0e-3] * 3
    assert [record[2] for record in by_cell.stress_period_data[0]] == [1.0e-3, 2.0e-3, 3.0e-3]
    # depth (index 3) uses the same machinery.
    assert [record[3] for record in by_cell.stress_period_data[0]] == [1.0, 2.0, 3.0]
    assert [record[2] for record in per_cell_seq.stress_period_data[0]] == [1.0e-3, 2.0e-3, 3.0e-3]


def test_evt_builder_resolves_surface_offset_constant_and_explicit_forms():
    offset = EVTBuilder(
        context=_context(), nper=1, rate=1.0e-3, depth=1.0,
        surface=CellSurfaceOffset("model_top", offset=-2.0),
    )
    constant = EVTBuilder(context=_context(), nper=1, rate=1.0e-3, depth=1.0, surface=55.0)
    explicit = EVTBuilder(
        context=_context(), nper=1, rate=1.0e-3, depth=1.0, surface=[10.0, 20.0, 30.0]
    )

    assert [record[1] for record in offset.stress_period_data[0]] == [98.0, 88.0, 78.0]
    assert [record[1] for record in constant.stress_period_data[0]] == [55.0, 55.0, 55.0]
    assert [record[1] for record in explicit.stress_period_data[0]] == [10.0, 20.0, 30.0]


def test_evt_builder_accepts_explicit_cells_and_boundnames():
    builder = EVTBuilder(
        context=_context(),
        nper=1,
        cells=[0, (1, 1)],
        rate=[1.0e-3, 2.0e-3],
        depth=1.0,
        surface=[12.0, 34.0],
        boundnames=True,
        name_by_cell={0: "upland", (1, 1): "lowland"},
    )

    assert builder.stress_period_data[0] == [
        [(0, 0), 12.0, 1.0e-3, 1.0, "upland"],
        [(1, 1), 34.0, 2.0e-3, 1.0, "lowland"],
    ]


def test_evt_builder_build_is_immutable_and_stamps_metadata():
    base = EVTBuilder(context=_context(), nper=1, rate=1.0e-3, depth=2.0)
    wetter = base.with_updates(rate=2.0e-3)

    assert base.stress_period_data[0][0][2] == 1.0e-3
    assert wetter.stress_period_data[0][0][2] == 2.0e-3
    built = wetter.build()
    assert built.metadata["builder"] == "EVTBuilder"
    assert built.options["nseg"] == 1
    with pytest.raises(TypeError):
        wetter.build(name="other")


def test_evt_builder_validates_inputs():
    with pytest.raises(ValueError, match="missing cell"):
        EVTBuilder(context=_context(), nper=1, rate={(0, 0): 1.0e-3}, depth=1.0).build()
    with pytest.raises(ValueError, match="rate must be nonnegative"):
        EVTBuilder(context=_context(), nper=1, rate=-1.0e-3, depth=1.0).build()
    with pytest.raises(ValueError, match="depth must be nonnegative"):
        EVTBuilder(context=_context(), nper=1, rate=1.0e-3, depth=-1.0).build()
    with pytest.raises(ValueError, match="inactive"):
        EVTBuilder(context=_context(), nper=1, cells=[(0, 1)], rate=1.0e-3, depth=1.0).build()
    # File-less builder cannot resolve a field-name (string) surface offset.
    with pytest.raises(ValueError, match="numeric CellSurfaceOffset"):
        EVTBuilder(
            context=_context(), nper=1, rate=1.0e-3, depth=1.0,
            surface=CellSurfaceOffset("model_top", offset="fld"),
        )


def test_package_api_evt_uses_builder_by_default_and_flopy_for_direct_data():
    built = mf.evt(context=_context(), nper=1, rate=2.0e-3, depth=2.5)
    assert built.metadata["builder"] == "EVTBuilder"

    direct = mf.evt(stress_period_data={0: [[(0, 3), 100.0, 2.0e-3, 2.5]]})
    assert direct.metadata == {}

    with pytest.raises(TypeError, match="context=, nper=, rate=, and depth="):
        mf.evt(context=_context(), nper=1, rate=2.0e-3)  # depth missing


def test_evt_builder_writes_real_flopy_package(tmp_path):
    builder = EVTBuilder(
        context=ModelContext(
            domain=np.array([[1, 1]]),
            surfaces=pd.DataFrame({"top": [10.0, 10.0], "bottom": [0.0, 0.0]}),
        ),
        nper=1,
        rate={(0, 0): 1.0e-3, (0, 1): 2.0e-3},
        depth=2.0,
    )
    flow = ModelSpec(
        "flow",
        "gwf",
        packages=(
            PackageSpec(
                "disv",
                flopy.mf6.ModflowGwfdisv,
                {
                    "nlay": 1,
                    "ncpl": 2,
                    "nvert": 6,
                    "vertices": [
                        [0, 0.0, 0.0],
                        [1, 1.0, 0.0],
                        [2, 1.0, 1.0],
                        [3, 0.0, 1.0],
                        [4, 2.0, 0.0],
                        [5, 2.0, 1.0],
                    ],
                    "cell2d": [
                        [0, 0.5, 0.5, 4, 0, 1, 2, 3],
                        [1, 1.5, 0.5, 4, 1, 4, 5, 2],
                    ],
                    "top": 10.0,
                    "botm": 0.0,
                },
            ),
            builder.build(),
        ),
    )

    built = SimulationSpec(
        "evt",
        models=(flow,),
        packages=(
            PackageSpec(
                "tdis",
                flopy.mf6.ModflowTdis,
                {"nper": 1, "perioddata": [(1.0, 1, 1.0)]},
            ),
        ),
    ).build_flopy(tmp_path)
    built.simulation.write_simulation(silent=True)

    assert (tmp_path / "flow.evt").exists()
