from __future__ import annotations

import flopy
import numpy as np
import pytest

from myflopy import ModelContext, ModelSpec, PackageSpec, RCHBuilder, SimulationSpec


def _context() -> ModelContext:
    return ModelContext(domain=np.array([[1, 0, 1], [1, 1, 0]]))


def test_rch_builder_selects_top_active_cells_and_applies_flat_rate():
    builder = RCHBuilder(context=_context(), nper=2, recharge=1.0e-4)

    assert builder.selected_cells == ((0, 0), (1, 1), (0, 2))
    assert builder.stress_period_data == {
        0: [[(0, 0), 1.0e-4], [(1, 1), 1.0e-4], [(0, 2), 1.0e-4]],
        1: [[(0, 0), 1.0e-4], [(1, 1), 1.0e-4], [(0, 2), 1.0e-4]],
    }


def test_rch_builder_accepts_period_values_and_cell_keyed_values():
    period_scalars = RCHBuilder(
        context=_context(),
        nper=2,
        recharge={0: 1.0e-4, 1: 2.0e-4},
    )
    by_cell = RCHBuilder(
        context=_context(),
        nper=2,
        recharge={(0, 0): 1.0e-4, (1, 1): 2.0e-4, (0, 2): 3.0e-4},
    )
    transient_by_cell = RCHBuilder(
        context=_context(),
        nper=2,
        recharge={
            0: {(0, 0): 1.0e-4, (1, 1): 2.0e-4, (0, 2): 3.0e-4},
            1: {(0, 0): 4.0e-4, (1, 1): 5.0e-4, (0, 2): 6.0e-4},
        },
    )

    assert [record[1] for record in period_scalars.stress_period_data[1]] == [2.0e-4] * 3
    assert [record[1] for record in by_cell.stress_period_data[0]] == [1.0e-4, 2.0e-4, 3.0e-4]
    assert [record[1] for record in transient_by_cell.stress_period_data[1]] == [
        4.0e-4,
        5.0e-4,
        6.0e-4,
    ]


def test_rch_builder_accepts_explicit_cells_and_boundnames():
    builder = RCHBuilder(
        context=_context(),
        nper=1,
        cells=[0, (1, 1)],
        recharge=[1.0e-4, 2.0e-4],
        boundnames=True,
        name_by_cell={0: "upland", (1, 1): "lowland"},
    )

    assert builder.stress_period_data[0] == [
        [(0, 0), 1.0e-4, "upland"],
        [(1, 1), 2.0e-4, "lowland"],
    ]


def test_rch_builder_build_is_immutable_and_accepts_no_configuration():
    base = RCHBuilder(context=_context(), nper=1, recharge=1.0e-4)
    wet = base.with_updates(recharge=2.0e-4)

    assert base.stress_period_data[0][0][1] == 1.0e-4
    assert wet.stress_period_data[0][0][1] == 2.0e-4
    assert wet.build().metadata["builder"] == "RCHBuilder"
    with pytest.raises(TypeError):
        wet.build(name="other")


def test_rch_builder_validates_recharge_inputs():
    with pytest.raises(ValueError, match="missing cell"):
        RCHBuilder(context=_context(), nper=1, recharge={(0, 0): 1.0e-4}).build()
    with pytest.raises(ValueError, match="nonnegative"):
        RCHBuilder(context=_context(), nper=1, recharge=-1.0e-4).build()
    with pytest.raises(ValueError, match="inactive"):
        RCHBuilder(context=_context(), nper=1, cells=[(0, 1)], recharge=1.0e-4).build()


def test_rch_builder_writes_real_flopy_310_package(tmp_path):
    builder = RCHBuilder(
        context=ModelContext(domain=np.array([[1, 1]])),
        nper=1,
        recharge={(0, 0): 1.0e-4, (0, 1): 2.0e-4},
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
        "rch",
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

    assert (tmp_path / "flow.rch").exists()
