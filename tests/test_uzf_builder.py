from __future__ import annotations

import numpy as np
import pytest

from myflopy import ModelContext, ModelSpec, PackageSpec, SimulationSpec, UZFBuilder


def _context() -> ModelContext:
    return ModelContext(
        domain=np.array(
            [
                [1, 0, 1],
                [1, 1, 1],
                [1, 1, 0],
            ]
        )
    )


def _builder(**updates) -> UZFBuilder:
    values = {
        "context": _context(),
        "nper": 2,
        "vks": 0.1,
        "thtr": 0.05,
        "thts": 0.30,
        "thti": 0.15,
        "finf": {0: [0.001, 0.002, 0.003], 1: [0.004, 0.005, 0.006]},
    }
    values.update(updates)
    return UZFBuilder(**values)


def test_uzf_builder_creates_vertical_chains_and_surface_period_data():
    builder = _builder()

    assert builder.uzf_cells == (
        (0, 0),
        (1, 0),
        (2, 0),
        (1, 1),
        (2, 1),
        (0, 2),
        (1, 2),
    )
    assert builder.surface_cells == ((0, 0), (1, 1), (0, 2))

    packagedata = builder.packagedata
    assert [record[2] for record in packagedata] == [1, 0, 0, 1, 0, 1, 0]
    assert [record[3] for record in packagedata] == [1, 2, -1, 4, -1, 6, -1]

    assert len(builder.perioddata[0]) == 3
    assert [record[0] for record in builder.perioddata[0]] == [0, 3, 5]
    assert [record[1] for record in builder.perioddata[1]] == [0.004, 0.005, 0.006]


def test_uzf_builder_stops_vertical_chain_at_inactive_gap():
    context = ModelContext(domain=np.array([[1], [0], [1]]))

    builder = _builder(context=context, finf=0.001)

    assert builder.uzf_cells == ((0, 0),)
    assert builder.packagedata[0][3] == -1


def test_uzf_builder_supports_surface_only_and_explicit_cells():
    surface = _builder(cells="surface_only")
    explicit = _builder(cells=[(0, 0), (1, 0), (2, 0)], finf=0.001)

    assert surface.uzf_cells == ((0, 0), (1, 1), (0, 2))
    assert [record[3] for record in surface.packagedata] == [-1, -1, -1]
    assert [record[3] for record in explicit.packagedata] == [1, 2, -1]


def test_uzf_builder_is_immutable_and_build_accepts_no_configuration():
    base = _builder()
    dry = base.with_updates(finf=0.0001)

    assert base.perioddata[0][0][1] == 0.001
    assert dry.perioddata[0][0][1] == 0.0001
    assert dry.build().metadata["builder"] == "UZFBuilder"
    with pytest.raises(TypeError):
        dry.build(name="other")


def test_uzf_builder_validates_physical_properties():
    with pytest.raises(ValueError, match="thtr must be less than thts"):
        _builder(thtr=0.3, thts=0.3).build()
    with pytest.raises(ValueError, match="finf must be nonnegative"):
        _builder(finf=-0.1).build()
    with pytest.raises(ValueError, match="inactive"):
        _builder(cells=[(0, 1)]).build()


def test_uzf_builder_allows_transient_time_series_names():
    builder = _builder(finf="rainfall_series", pet="pet_series")

    builder.validate()

    assert builder.perioddata[0][0][1:3] == ["rainfall_series", "pet_series"]


def test_uzf_builder_writes_real_flopy_310_package(tmp_path):
    builder = UZFBuilder(
        context=ModelContext(domain=np.array([[1, 1]])),
        nper=1,
        vks=0.1,
        thtr=0.05,
        thts=0.30,
        thti=0.15,
        finf=[0.001, 0.002],
    )
    flow = ModelSpec(
        "flow",
        "gwf",
        packages=(
            PackageSpec(
                "dis",
                __import__("flopy").mf6.ModflowGwfdisv,
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
        "uzf",
        models=(flow,),
        packages=(
            PackageSpec(
                "tdis",
                __import__("flopy").mf6.ModflowTdis,
                {"nper": 1, "perioddata": [(1.0, 1, 1.0)]},
            ),
        ),
    ).build_flopy(tmp_path)
    built.simulation.write_simulation(silent=True)

    assert (tmp_path / "flow.uzf").exists()
