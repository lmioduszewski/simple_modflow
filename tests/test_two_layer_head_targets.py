"""A well screened in two layers is two observations, not one that overwrites itself.

The PEST side keyed everything on the bare target NAME. A name in two layers
therefore had one simulated column, written twice by the forward run (the last
mapping row won) and read once by the build (the pivot kept the first): 1005.0
against 5.0 on a two-layer model. The MF6 records had always disambiguated such
a name as ``<name>_lay<layer>``; that label is now the one definition.
"""

from __future__ import annotations

from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point

from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.observations.heads import HeadTargets
from myflopy.modflow.mf6.pest.forward_run import _write_head_target_csv


def _grid(n: int = 4, length: float = 100.0) -> VoronoiGridPlus:
    """n x n square cells, row-major from the bottom left."""

    d = length / n
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def _model(grid, nlay: int = 2):
    """A model-like object with a DISV idomain and a head table (head = 1000*layer + cell)."""

    ncpl = len(grid.gdf_vorPolys)
    heads = pd.DataFrame(
        [
            {"kstpkper": (0, 0), "layer": layer, "cell": cell, "elev": 1000.0 * layer + cell}
            for layer in range(nlay)
            for cell in range(ncpl)
        ]
    )
    return SimpleNamespace(
        vor=grid,
        gwf=SimpleNamespace(
            disv=SimpleNamespace(idomain=SimpleNamespace(array=np.ones((nlay, ncpl), dtype=int)))
        ),
        all_heads=heads,
    )


def _targets(named):
    """Point targets as (name, layer, point) triples."""

    frame = gpd.GeoDataFrame(
        {"name": [name for name, _, _ in named], "layer": [layer for _, layer, _ in named]},
        geometry=[point for _, _, point in named],
        crs="EPSG:2927",
    )
    names = list(dict.fromkeys(name for name, _, _ in named))
    values = pd.DataFrame({"per": [0], **{name: [1.0] for name in names}})
    return HeadTargets(locations=frame, values=values, time_column="per")


@pytest.fixture
def grid():
    return _grid()


def test_a_unique_name_is_its_own_label(grid):
    """The common case must not be renamed: existing pst files name observations
    after these columns."""

    labels = _targets([("W1", 0, Point(37.5, 37.5)), ("W2", 0, Point(62.5, 62.5))])
    labels = labels.observation_labels(_model(grid))

    assert labels["obsname"].tolist() == ["W1", "W2"]


def test_a_name_in_two_layers_becomes_two_labels(grid):
    labels = _targets(
        [("W1", 0, Point(37.5, 37.5)), ("W1", 1, Point(37.5, 37.5))]
    ).observation_labels(_model(grid))

    assert labels["obsname"].tolist() == ["W1_lay0", "W1_lay1"]
    assert labels["cell"].tolist() == [5, 5]


def test_the_flopy_records_and_the_labels_are_one_definition(grid):
    targets = _targets([("W1", 0, Point(37.5, 37.5)), ("W1", 1, Point(37.5, 37.5))])
    model = _model(grid)

    records = targets.to_flopy_obs(model)["head_targets.csv"]

    assert [name for name, _kind, _cellid in records] == \
        targets.observation_labels(model)["obsname"].tolist()


def test_each_layer_gets_its_own_simulated_column(grid):
    """The defect: one column, holding layer 0's head, against a forward run that
    wrote layer 1's into it."""

    targets = _targets([("W1", 0, Point(37.5, 37.5)), ("W1", 1, Point(37.5, 37.5))])

    simulated = targets.simulated_heads(_model(grid))

    assert list(simulated.columns) == ["per", "W1_lay0", "W1_lay1"]
    assert simulated.loc[0, "W1_lay0"] == 5.0
    assert simulated.loc[0, "W1_lay1"] == 1005.0


def test_the_build_and_the_forward_run_agree_on_the_columns(grid, tmp_path, monkeypatch):
    """The forward run writes what the build registered, for every layer, IN ORDER.

    Names deliberately out of sorted order: PEST reads the file by position and
    names each value after the BASE copy's header, and the base copy was a pivot,
    whose columns pandas sorts -- ``W1_lay0, W1_lay1, W2`` against the forward
    run's ``W2, W1_lay0, W1_lay1``, so every realization handed PEST W2's head
    as W1's. Names already in sorted order could not see it (ledger 340)."""

    from myflopy.modflow.mf6.pest.observations import prepare_head_target_observations
    from myflopy.modflow.mf6.pest.specs import HeadTargetObservationSpec

    targets = _targets([
        ("W2", 0, Point(62.5, 62.5)),
        ("W1", 0, Point(37.5, 37.5)),
        ("W1", 1, Point(37.5, 37.5)),
    ])
    model = _model(grid)
    model.name, model.workspace = "flow", tmp_path
    _stub_head_file(monkeypatch)
    prepare_head_target_observations(
        SimpleNamespace(model=model, template_workspace=tmp_path,
                        pf=SimpleNamespace(add_observations=lambda *args, **kwargs: None)),
        HeadTargetObservationSpec(targets=targets, prefix="hds"),
    )
    base = pd.read_csv(tmp_path / "hds_simulated_heads.csv")

    written = _run_forward(tmp_path, monkeypatch, tmp_path / "hds_head_target_map.csv")

    assert list(written.columns) == list(base.columns) == ["per", "W2", "W1_lay0", "W1_lay1"]
    assert written.loc[0, "W1_lay0"] == 1005.0      # layer 0 of the fake head file
    assert written.loc[0, "W1_lay1"] == 2005.0
    assert written.loc[0, "W2"] == 1010.0


def test_a_mapping_written_before_the_label_still_works(tmp_path, monkeypatch):
    """A run prepared earlier has no label column; its name still serves."""

    mapping = tmp_path / "map.csv"
    pd.DataFrame({"name": ["W1"], "layer": [0], "cell": [5]}).to_csv(mapping, index=False)

    written = _run_forward(tmp_path, monkeypatch, mapping)

    assert list(written.columns) == ["per", "W1"]


def _stub_head_file(monkeypatch):
    """Stand a two-layer, 16-cell fake in for FloPy's head file (layer 0 = 1000 + cell)."""

    import flopy

    class _FakeHeadFile:
        # The file's own clock, one-based as on disk: the forward run picks each
        # period's END from it (ledger 336).
        recordarray = np.array(
            [(1, 1, 1.0)], dtype=[("kstp", int), ("kper", int), ("pertim", float)]
        )

        def __init__(self, *_args, **_kwargs):
            pass

        def get_kstpkper(self):
            return [(0, 0)]

        def get_data(self, kstpkper=None):
            return np.array([[1000.0 + cell for cell in range(16)],
                             [2000.0 + cell for cell in range(16)]])

        def close(self):
            pass

    monkeypatch.setattr(flopy.utils, "HeadFile", _FakeHeadFile)


def _run_forward(tmp_path, monkeypatch, mapping) -> pd.DataFrame:
    """Run the injected forward-run writer against a stub head file."""

    _stub_head_file(monkeypatch)
    monkeypatch.chdir(tmp_path)
    output = tmp_path / "simulated.csv"
    _write_head_target_csv("flow", str(mapping), str(output))
    return pd.read_csv(output)


def test_the_mapping_and_the_pst_rows_carry_the_label(grid, tmp_path, monkeypatch):
    """pyEMU names an observation after the simulated CSV's COLUMN, so the mapping
    the forward run reads and the pst rows must carry the same label it does.

    The build writes that CSV with the forward run's own writer (ledger 340), so
    the model double needs a head file: a stub, carrying MF6's clock."""

    from myflopy.modflow.mf6.pest.observations import prepare_head_target_observations
    from myflopy.modflow.mf6.pest.specs import HeadTargetObservationSpec

    _stub_head_file(monkeypatch)
    targets = _targets([("W1", 0, Point(37.5, 37.5)), ("W1", 1, Point(37.5, 37.5))])
    model = _model(grid)
    model.name, model.workspace = "flow", tmp_path
    project = SimpleNamespace(
        model=model,
        template_workspace=tmp_path,
        pf=SimpleNamespace(add_observations=lambda *args, **kwargs: None),
    )

    prepared = prepare_head_target_observations(
        project, HeadTargetObservationSpec(targets=targets, prefix="hds")
    )

    mapping = pd.read_csv(tmp_path / "hds_head_target_map.csv")
    simulated = pd.read_csv(tmp_path / "hds_simulated_heads.csv")
    assert mapping["obsname"].tolist() == ["W1_lay0", "W1_lay1"]
    assert sorted(set(prepared["target_frame"]["column"])) == ["W1_lay0", "W1_lay1"]
    assert set(simulated.columns) - {"per"} == set(mapping["obsname"])
