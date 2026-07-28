"""Budget node ids are zero-based on EVERY path, for every cell-based BC.

Regression guard for the coupled off-by-one found by the 2026-07-18 review:

* ``budget_tables._zero_base_budget_frame`` normalized a hardcoded
  ``{"drn", "ghb", "rch"}`` only, so ``model.bud("chd"|"riv"|"wel"|"evt").df``
  handed back MF6's raw 1-based node ids.
* ``GroupBudget._normalize_nodes`` then re-normalized with an
  ``if node.min() >= 1: subtract 1`` heuristic. That heuristic is correct on RAW
  MF6 records (whose nodes are always >= 1) but was being applied to output that
  ``budget_df`` had ALREADY zero-based -- so ``group.bud("drn"|"ghb")`` came back
  one cell low, while ``rch`` escaped only because it covers node 0.

The two bugs cancelled for chd/riv/wel/evt, so fixing either alone traded one
off-by-one for another. These tests pin both paths against ground truth taken
from the packages' own stress-period data.
"""

from __future__ import annotations

import pandas as pd
import pytest

from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS
from myflopy.project.model_group import ModelGroup

CELL_BC_PACKAGES = sorted(
    name for name, spec in _PACKAGE_EXPLORER_SPECS.items() if spec.kind == "cell_stress"
)

#: The advanced packages, which the cell-BC parametrization above never reached.
ADVANCED_BUDGET_PACKAGES = ("sfr", "lak", "uzf")


def _truth_nodes(model, package: str) -> set[int]:
    """Zero-based flat node ids a package actually occupies, from its own inputs."""

    ncpl = int(model.vor.ncpl)
    records = model.gwf.get_package(package).stress_period_data.get_data(key=0)
    return {int(row[0][0]) * ncpl + int(row[0][1]) for row in records}


def _truth_cells_from_cellid(model, package: str) -> set[int]:
    """Zero-based flat node ids an ADVANCED package occupies, from its own inputs.

    SFR and UZF list their cells in ``packagedata``; LAK lists them per
    connection in ``connectiondata``.
    """

    ncpl = int(model.vor.ncpl)
    owner = getattr(model, package)
    table = owner.connectiondata if package == "lak" else owner.packagedata
    rows = pd.DataFrame(table.get_data())
    return {int(cellid[0]) * ncpl + int(cellid[1]) for cellid in rows["cellid"]}


def test_cell_bc_registry_covers_the_packages_under_test():
    # if a new cell BC is added, it is covered here automatically -- and if the
    # registry ever loses one, this fails rather than silently shrinking coverage
    assert set(CELL_BC_PACKAGES) == {"chd", "drn", "evt", "ghb", "rch", "riv", "wel"}


@pytest.mark.parametrize("package", CELL_BC_PACKAGES)
def test_model_budget_nodes_are_zero_based(canonical_run, package):
    """``model.bud(pkg).df`` node ids must be zero-based model cells."""

    frame = canonical_run.bud(package).df.reset_index()
    column = "node" if "node" in frame.columns else frame.columns[0]
    nodes = {int(value) for value in frame[column].dropna().unique()}

    assert nodes, f"{package} budget produced no nodes"
    assert nodes <= _truth_nodes(canonical_run, package), (
        f"{package}: model.bud() node ids are not the package's zero-based cells "
        f"(smallest returned {min(nodes)}, smallest true {min(_truth_nodes(canonical_run, package))})"
    )


@pytest.mark.parametrize("package", CELL_BC_PACKAGES)
def test_group_budget_nodes_match_the_single_model_path(canonical_run, package):
    """``group.bud(pkg)`` must agree with ``model.bud(pkg)`` -- no re-normalization."""

    group = ModelGroup({"a": canonical_run}, reference="a")
    grouped = {int(v) for v in group.bud(package).get()["node"].dropna().unique()}

    assert grouped, f"{package} group budget produced no nodes"
    assert grouped <= _truth_nodes(canonical_run, package), (
        f"{package}: group.bud() node ids drifted from the package's zero-based "
        f"cells (smallest returned {min(grouped)})"
    )

    frame = canonical_run.bud(package).df.reset_index()
    column = "node" if "node" in frame.columns else frame.columns[0]
    single = {int(value) for value in frame[column].dropna().unique()}
    assert grouped == single, (
        f"{package}: the group and single-model budget paths disagree "
        f"(group-only {sorted(grouped - single)[:5]}, model-only {sorted(single - grouped)[:5]})"
    )


@pytest.mark.parametrize("package", ADVANCED_BUDGET_PACKAGES)
def test_advanced_package_budget_nodes_are_zero_based(canonical_run, package):
    """In the MODEL budget, ``node`` is a CELL for SFR/LAK/UZF too (2026-07-27).

    The registry marked these three ``zero_base_budget_nodes=False`` on the
    stated belief that "their records carry feature ids, not model cells". That
    was measured false: in the *model* budget file the ``SFR``/``LAK``/
    ``UZF-GWRCH`` records put the 1-based model cell in ``node`` and the feature
    id in ``node2``. (The feature-first layout the belief described is real, but
    it belongs to the separate *package-output* budget file, which
    ``model.outputs.<pkg>.bud`` reads.)

    So ``model.bud("sfr"|"lak"|"uzf").df`` handed back node ids exactly one cell
    high -- the same off-by-one this module exists to prevent, escaping only
    because the parametrization above covers ``cell_stress`` packages alone.
    """

    frame = canonical_run.bud(package).df.reset_index()
    nodes = {int(value) for value in frame["node"].dropna().unique()}
    truth = _truth_cells_from_cellid(canonical_run, package)

    assert nodes, f"{package} budget produced no nodes"
    assert nodes <= truth, (
        f"{package}: model.bud() node ids are not the package's zero-based cells "
        f"(smallest returned {min(nodes)}, smallest true {min(truth)})"
    )


def test_a_full_array_budget_record_builds_a_real_per_cell_table(canonical_run):
    """imeth=1 records are positional arrays, and must NOT be shifted (2026-07-27).

    ``STO-SS`` is written as a plain ``(nlay, 1, ncpl)`` float array with no
    ``node`` column at all. It used to reach ``pd.DataFrame(record)`` and die on
    pandas' opaque ``"Must pass 2-d input"``. Now it builds one row per model
    node -- and because a full array is *positional*, node 0 is index 0, so the
    1-based shift the recarray branch applies must not happen here.
    """

    frame = canonical_run.bud("STO-SS").df.reset_index()
    nodes = sorted({int(value) for value in frame["node"].dropna().unique()})

    assert nodes == list(range(len(canonical_run.node_to_lni))), (
        "a full-array record covers every model node exactly once, zero-based"
    )
    assert "q" in frame.columns and frame["q"].notna().all()


def test_a_connection_indexed_record_says_so_instead_of_returning_nothing(canonical_run):
    """FLOW-JA-FACE cannot be a cell table, and must refuse rather than go quiet.

    It is a full array indexed by cell CONNECTION, so there is no honest mapping
    onto cells. The bug being fixed was silence, so the replacement must be an
    error that names the mismatch -- not an empty frame.
    """

    with pytest.raises(ValueError, match="cannot be mapped to cells"):
        canonical_run.bud("FLOW-JA-FACE").df

    # the alias route must refuse too: MF6's own record name is not what a
    # caller necessarily types, and "flow" reaches FLOW-JA-FACE by substring
    with pytest.raises(ValueError, match="cannot be mapped to cells"):
        canonical_run.bud("flow").df


def test_flow_ja_face_is_refused_by_name_even_when_its_length_looks_like_cells(tmp_path):
    """The connection guard must key on the record NAME, never on its length.

    Under an idomain reduction MF6 expands cell arrays back to ``nodesuser``
    while FLOW-JA-FACE stays at the reduced ``nja``, so the two counts are
    independent and CAN collide. This model is the smallest case where they do:
    a 1-layer DISV with ``ncpl=4`` and two active adjacent cells gives
    ``nja == nodesuser == 4``, and both records arrive with the identical shape
    ``(1, 1, 4)`` -- so neither size nor shape can tell them apart.

    A size-only guard (the first version of this fix, 2026-07-27) accepted that
    FLOW-JA-FACE as a four-cell table and reported flow through the two
    ``idomain=0`` cells, with no error at all -- strictly worse than the empty
    frame it replaced, and on the very failure mode the fix existed to remove.
    Every other test here uses an all-ones idomain, where ``nja > nodesuser`` is
    guaranteed, so none of them can see this.
    """

    import numpy as np

    import myflopy as mf
    from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
    from myflopy.specs import ModelContext

    vor = irregular_voronoi_grid(nrow=2, ncol=2, cell_size=100.0)
    gridprops = vor.get_disv_gridprops()
    ncpl = gridprops["ncpl"]
    idomain = np.zeros((1, ncpl), dtype=int)
    idomain[0, 0] = 1
    idomain[0, 1] = 1

    model_spec = mf.gwf(
        "m",
        context=ModelContext(grid=vor, domain=idomain),
        packages=[
            mf.disv(
                nlay=1, ncpl=ncpl, nvert=len(gridprops["vertices"]),
                vertices=gridprops["vertices"], cell2d=gridprops["cell2d"],
                top=10.0, botm=0.0, idomain=idomain,
            ),
            mf.ic(strt=9.0),
            mf.npf(k=1.0),
            mf.chd(stress_period_data={0: [[(0, 0), 10.0], [(0, 1), 8.0]]}),
            mf.oc(
                head_filerecord="m.hds", budget_filerecord="m.cbc",
                saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
            ),
        ],
    )
    sim = mf.SimulationSpec(
        "collide",
        models=[model_spec],
        packages=[
            mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
            mf.ims(models=["m"], complexity="SIMPLE"),
        ],
    )
    run = mf.Project(tmp_path / "collide", name="p").prepare_run(
        "r", sim, overwrite=True
    )
    success, report = run.execute()
    assert success, "\n".join(report[-15:])
    model = run.model("m")

    # the collision is the whole point of this fixture -- assert it really holds,
    # or the test silently stops covering anything
    reader = model._get_budget_reader()
    record = reader.get_data(
        text="FLOW-JA-FACE", kstpkper=model._get_budget_kstpkper()[0]
    )[0]
    assert np.asarray(record).size == len(model.node_to_lni), (
        "fixture no longer collides, so it cannot catch a size-only guard"
    )

    with pytest.raises(ValueError, match="indexed by cell CONNECTION"):
        model.bud("FLOW-JA-FACE").df


@pytest.mark.parametrize("package", CELL_BC_PACKAGES)
def test_modern_explorer_path_agrees_with_the_legacy_budget_path(canonical_run, package):
    """The registry-backed explorer and the legacy ``bud()`` path must not disagree."""

    ncpl = int(canonical_run.vor.ncpl)
    results = getattr(canonical_run.packages, package).results.q.get()
    explorer = {
        int(row.layer) * ncpl + int(row.cell) for row in results.itertuples()
    }

    frame = canonical_run.bud(package).df.reset_index()
    column = "node" if "node" in frame.columns else frame.columns[0]
    legacy = {int(value) for value in frame[column].dropna().unique()}

    assert explorer and legacy
    assert explorer == legacy, (
        f"{package}: explorer and legacy budget paths disagree by "
        f"{sorted(explorer ^ legacy)[:5]}"
    )
