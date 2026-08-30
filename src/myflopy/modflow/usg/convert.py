"""Map an imported MODFLOW-USG model onto a MODFLOW 6 :class:`SimulationSpec`.

Every judgment this module makes is recorded in :func:`conversion_notes`, which
is what :meth:`UsgModel.report` prints. The rule is that a conversion may be
lossy but never silently lossy: anything approximated or omitted is counted and
named, with the cells it touched.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

import myflopy.package_api as api
from myflopy._logging import get_logger
from myflopy.specs import ModelContext, PackageSpec, SimulationSpec

if TYPE_CHECKING:
    from myflopy.modflow.usg.model import UsgModel

logger = get_logger(__name__)

__all__ = ["Finding", "conversion_notes", "mf6_findings", "to_mf6"]

#: USG list packages and the MODFLOW 6 factory that receives them.
_LIST_TARGETS = {"CHD": "chd", "DRN": "drn", "GHB": "ghb", "RIV": "riv", "WEL": "wel"}


@dataclass(frozen=True, slots=True)
class Finding:
    """One way in which MODFLOW 6 is stricter than MODFLOW-USG about this model.

    MODFLOW-USG accepts inputs MODFLOW 6 refuses -- most often a head-type
    boundary placed below the bottom of the cell it sits in, which USG treats as
    a strong sink and MODFLOW 6 rejects outright. Finding these before the run
    is worth more than discovering them from a failed run, because the listing
    file names six of them and stops.
    """

    code: str
    severity: str
    count: int
    detail: str
    fix: str

    @property
    def blocks_run(self) -> bool:
        """True when MODFLOW 6 will refuse to run until this is resolved."""

        return self.severity == "reject"


def mf6_findings(model: UsgModel) -> list[Finding]:
    """Return everything MODFLOW 6 will object to in this model, before running it."""

    findings: list[Finding] = []
    for ftype, column, label in (("GHB", 1, "bhead"), ("RIV", 1, "stage"), ("DRN", 1, "elev")):
        records = model.boundaries.get(ftype)
        if records is None:
            continue
        offenders = 0
        worst = 0.0
        share = 0.0
        total = 0.0
        for block in records.periods.values():
            if not len(block):
                continue
            cellid = model.to_cellid(block[:, 0].astype(np.int64))
            valid = cellid[:, 0] >= 0
            if not valid.any():
                continue
            bottom = model.botm[cellid[valid, 0], cellid[valid, 1]]
            below = block[valid, column] < bottom
            offenders = max(offenders, int(below.sum()))
            if below.any():
                worst = max(worst, float((bottom - block[valid, column])[below].max()))
                if block.shape[1] > column + 1:
                    share = max(share, float(block[valid, column + 1][below].sum()))
                    total = max(total, float(block[valid, column + 1].sum()))
        if offenders and ftype != "DRN":  # a DRN below its cell bottom is simply inactive
            portion = f", carrying {share:.3g} of {total:.3g} total conductance" if total else ""
            findings.append(
                Finding(
                    code=f"{ftype.lower()}_below_cell_bottom",
                    severity="reject",
                    count=offenders,
                    detail=(
                        f"{offenders} {ftype} record(s) have {label} below the bottom of "
                        f"their cell (by up to {worst:.2f}){portion}. MODFLOW-USG allows "
                        "this; MODFLOW 6 rejects it and stops."
                    ),
                    fix=f"to_mf6(fix_for_mf6=True) raises each {label} to its cell bottom",
                )
            )

    chd = model.boundaries.get("CHD")
    if chd is not None:
        worst_count = 0
        for block in chd.periods.values():
            if not len(block):
                continue
            cellid = model.to_cellid(block[:, 0].astype(np.int64))
            valid = cellid[:, 0] >= 0
            if not valid.any():
                continue
            bottom = model.botm[cellid[valid, 0], cellid[valid, 1]]
            head = block[valid, 2] if block.shape[1] > 2 else block[valid, 1]
            worst_count = max(worst_count, int((head < bottom).sum()))
        if worst_count:
            findings.append(
                Finding(
                    code="chd_below_cell_bottom",
                    severity="reject",
                    count=worst_count,
                    detail=(
                        f"up to {worst_count} CHD record(s) per period hold a head below their "
                        "cell bottom. MODFLOW 6 rejects that outright. MODFLOW-USG accepts it, "
                        "but the cell then has zero saturated thickness and passes no water -- "
                        "so the boundary was already inert wherever this happens."
                    ),
                    fix=(
                        "to_mf6(fix_for_mf6=True) omits each such record in the periods where "
                        "it is below the bottom, and keeps it in the periods where it is not"
                    ),
                )
            )
    return findings


def to_mf6(
    model: UsgModel,
    *,
    name: str = "usg",
    crs: str | None = None,
    newton: bool | None = None,
    under_relaxation: bool = True,
    start_date_time: str | None = None,
    save_budget: bool = True,
    complexity: str = "COMPLEX",
    include: tuple[str, ...] | None = None,
    fix_for_mf6: bool = False,
    local_origin: bool = True,
) -> SimulationSpec:
    """Build a MODFLOW 6 simulation spec from a USG model. See :meth:`UsgModel.to_mf6`."""

    if model.grid is None:
        raise ValueError(
            "conversion needs the grid geometry: a USG DISU carries no coordinates, so "
            "read_usg() must be given the model's .gsf file"
        )

    blocking = [f for f in mf6_findings(model) if f.blocks_run]
    if blocking and not fix_for_mf6:
        for finding in blocking:
            logger.warning("MODFLOW 6 will reject this model: %s %s", finding.detail, finding.fix)

    wanted = {p.lower() for p in include} if include is not None else None
    keep = lambda pkg: wanted is None or pkg in wanted  # noqa: E731

    grid_props = model.grid.get_disv_gridprops()
    vertices, cell2d, xorigin, yorigin = _localize(grid_props, shift=local_origin)
    packages: list[PackageSpec] = [
        api.disv(
            nlay=model.nlay,
            ncpl=grid_props["ncpl"],
            nvert=len(vertices),
            vertices=vertices,
            cell2d=cell2d,
            top=model.top,
            botm=model.botm,
            idomain=model.idomain,
            xorigin=xorigin,
            yorigin=yorigin,
        ),
        api.npf(k=model.k, k33=model.k33, icelltype=model.icelltype, save_flows=True),
        api.ic(strt=model.strt),
    ]

    # Pass a period map only when it has entries. An EMPTY dict is not the same
    # as omitting it: FloPy treats `steady_state={}` as a supplied period map and
    # writes `BEGIN period / END period` with no keyword inside, which silently
    # loses the TRANSIENT flags as well -- MODFLOW 6 then defaults every period
    # to steady state, and the model solves to NaN rather than failing loudly.
    steady_flags = np.asarray(model.disu.steady, dtype=bool)
    transient = {i: True for i, flag in enumerate(steady_flags) if not flag}
    steady = {i: True for i, flag in enumerate(steady_flags) if flag}
    packages.append(
        api.sto(
            ss=model.ss,
            sy=model.sy,
            iconvert=np.broadcast_to(model.icelltype[:, None], (model.nlay, model.ncpl)).copy(),
            transient=transient or None,
            steady_state=steady or None,
        )
    )

    for ftype, target in _LIST_TARGETS.items():
        records = model.boundaries.get(ftype)
        if records is None or not keep(target):
            continue
        data = _list_period_data(model, records, fix_for_mf6=fix_for_mf6)
        if not data:
            logger.debug("%s has no groundwater records after CLN filtering; skipped", ftype)
            continue
        factory = getattr(api, target)
        packages.append(factory.flopy(stress_period_data=data, name=target))

    if model.rch is not None and keep("rch"):
        packages.append(_recharge_spec(model))
    if model.ets is not None and keep("evt"):
        packages.append(_evt_spec(model))
    if model.hfb.size and keep("hfb"):
        spec = _hfb_spec(model)
        if spec is not None:
            packages.append(spec)

    packages.append(
        api.oc(
            head_filerecord=f"{name}.hds",
            budget_filerecord=f"{name}.cbc" if save_budget else None,
            saverecord=[("HEAD", "ALL")] + ([("BUDGET", "ALL")] if save_budget else []),
        )
    )

    use_newton = newton if newton is not None else bool(model.sms and model.sms.is_newton)
    newton_options = None
    if use_newton:
        newton_options = "NEWTON UNDER_RELAXATION" if under_relaxation else "NEWTON"

    flow = api.gwf(
        name,
        packages=packages,
        context=ModelContext(grid=model.grid),
        newtonoptions=newton_options,
        save_flows=True,
    )

    perioddata = [
        (float(p), int(n), float(t))
        for p, n, t in zip(model.disu.perlen, model.disu.nstp, model.disu.tsmult, strict=True)
    ]
    simulation_packages = [
        api.tdis(
            nper=len(perioddata),
            perioddata=perioddata,
            time_units=model.time_units,
            start_date_time=start_date_time or model.start_date_time,
        ),
        _ims_spec(model, name=name, complexity=complexity),
    ]

    logger.info(
        "converted %s to MODFLOW 6: %d packages on a %d x %d DISV grid, %d periods",
        model.name_file.path.name,
        len(packages),
        model.nlay,
        model.ncpl,
        len(perioddata),
    )
    if crs:
        model.grid.crs = crs
    return SimulationSpec(name, models=(flow,), packages=tuple(simulation_packages))


# ---------------------------------------------------------------- pieces


def _localize(grid_props: dict, *, shift: bool):
    """Move the mesh next to the origin and declare where it came from.

    MODFLOW 6 builds DISV connection geometry -- shared edge lengths and the
    centre-to-edge distances that become conductances -- from the raw vertex
    coordinates. On a State Plane grid those are around 1.34 MILLION feet, and
    the differences that matter are a few feet, so the subtraction loses most of
    its significant digits. For the smallest cells here (under 3 square feet)
    the result degrades to NaN, and MODFLOW 6 then reports a NaN budget and
    "Normal termination" rather than an error: the run looks like it worked.

    Measured on the Ten Trails grid: identical input, run twice, differing only
    in this shift -- 0 of 9405 cells finite as-is, 9405 of 9405 finite shifted.

    Shifting the mesh to its own lower-left corner and passing that corner as
    ``xorigin``/``yorigin`` keeps the model georeferenced (FloPy's model grid and
    every export put the cells back where they belong) while the arithmetic
    happens near zero.
    """

    vertices = grid_props["vertices"]
    cell2d = grid_props["cell2d"]
    if not shift:
        return vertices, cell2d, 0.0, 0.0

    xorigin = float(min(v[1] for v in vertices))
    yorigin = float(min(v[2] for v in vertices))
    moved_vertices = [[v[0], v[1] - xorigin, v[2] - yorigin] for v in vertices]
    moved_cells = [
        [r[0], r[1] - xorigin, r[2] - yorigin, r[3], *r[4 : 4 + r[3]]] for r in cell2d
    ]
    logger.debug("shifted the mesh by (%.1f, %.1f) and declared it as the DISV origin", xorigin, yorigin)
    return moved_vertices, moved_cells, xorigin, yorigin


def _list_period_data(model: UsgModel, records, *, fix_for_mf6: bool = False) -> dict[int, list[Any]]:
    """Convert one list boundary condition to a MODFLOW 6 period dict.

    A USG period whose ``ITMP`` was negative reused the previous period. MODFLOW
    6 expresses exactly that by omitting the period from the dict, so a reused
    period is simply absent here -- no data is duplicated and no period is
    accidentally emptied. (An *empty list* would mean something else entirely in
    MODFLOW 6: delete every record.)
    """

    ftype = records.ftype
    data: dict[int, list[Any]] = {}
    clamped = [0]
    dropped = [0]
    for period, block in sorted(records.periods.items()):
        if len(block) == 0:
            continue
        nodes = block[:, 0].astype(np.int64)
        cellid = model.to_cellid(nodes)
        valid = cellid[:, 0] >= 0
        values_block = block[valid, 1:]
        if fix_for_mf6 and ftype == "CHD":
            # A constant head below its cell bottom is inert in MODFLOW-USG
            # (zero saturated thickness) and fatal in MODFLOW 6, so omitting it
            # reproduces the USG behaviour exactly. It has to be decided PER
            # PERIOD: these heads follow a seasonal cycle, so the same record is
            # below the bottom in a low-stage month and above it in a high one.
            bottom = model.botm[cellid[valid, 0], cellid[valid, 1]]
            head = values_block[:, 1] if values_block.shape[1] > 1 else values_block[:, 0]
            guard = np.maximum(np.abs(bottom) * 1e-6, 1e-6)
            keep_row = head > bottom + guard
            if not keep_row.all():
                dropped[0] += int((~keep_row).sum())
                cellid = cellid[valid][keep_row]
                values_block = values_block[keep_row]
                valid = np.ones(len(cellid), dtype=bool)

        if fix_for_mf6 and ftype in ("GHB", "RIV"):
            bottom = model.botm[cellid[valid, 0], cellid[valid, 1]]
            below = values_block[:, 0] < bottom
            if below.any():
                values_block = values_block.copy()
                # Clamp to just ABOVE the bottom, not to it. MODFLOW 6 compares
                # the two numbers as it read them back from text, and both are
                # written with finite precision -- so a head set exactly equal to
                # the bottom still trips "HEAD IS LESS THAN CELL BOTTOM" once the
                # last digit rounds the wrong way. The nudge is relative so it
                # stays negligible at any elevation (5e-5 ft at 50 ft).
                target = bottom[below]
                values_block[below, 0] = target + np.maximum(np.abs(target) * 1e-6, 1e-6)
                clamped[0] += int(below.sum())
        rows = []
        for (layer, cell), values in zip(cellid[valid], values_block, strict=True):
            rows.append(_bc_record(ftype, (int(layer), int(cell)), values))
        if rows:
            data[period] = rows
    if dropped[0]:
        logger.warning(
            "%s: omitted %d record-period(s) whose head sits below the cell bottom -- inert in "
            "MODFLOW-USG, rejected by MODFLOW 6 (fix_for_mf6=True)",
            ftype,
            dropped[0],
        )
    if clamped[0]:
        logger.warning(
            "%s: raised %d boundary head(s) to their cell bottom so MODFLOW 6 accepts them "
            "(fix_for_mf6=True)",
            ftype,
            clamped[0],
        )
    return data


def _bc_record(ftype: str, cellid: tuple[int, int], values: np.ndarray) -> list[Any]:
    """Build one MODFLOW 6 boundary record from a USG record's trailing values.

    ``CHD`` is the one that loses information: MODFLOW-USG ramps a constant head
    from ``shead`` to ``ehead`` across the stress period, while MODFLOW 6 holds
    one value for the whole period. The end-of-period head is taken, because
    that is the value MODFLOW-USG itself reaches by the end of the period and so
    the one a following period continues from.
    """

    if ftype == "CHD":
        return [cellid, float(values[1] if values.size > 1 else values[0])]
    return [cellid, *(float(v) for v in values)]


def _recharge_spec(model: UsgModel) -> PackageSpec:
    """Build MODFLOW 6 ``RCHA`` from the USG ``RCH`` arrays.

    USG recharge is already one array per period, so the array form is the
    faithful shape as well as the compact one -- on the Ten Trails model, 0.88 s
    end to end against 5.35 s for the equivalent list.

    ``NRCHOP = 1`` (top grid layer) and ``NRCHOP = 3`` (highest ACTIVE cell) both
    resolve through ``irch``: MODFLOW 6 without an IRCH array applies recharge to
    layer 1 unconditionally and silently skips any column whose layer 1 is
    inactive, which is not what either option means.
    """

    recharge = {period: np.asarray(values) for period, values in sorted(model.rch.rech.items())}
    return api.rch.array(
        recharge=recharge,
        context=ModelContext(grid=model.grid, domain=model.idomain),
        name="rch",
    )


def _evt_spec(model: UsgModel) -> PackageSpec:
    """Build MODFLOW 6 ``EVT`` from the USG ``ETS`` package.

    Two mismatches are resolved here, both forced:

    ``NETSOP = 3`` applies ET to the highest *active* cell in each column.
    MODFLOW 6 has no such option, so the cell is resolved once from IBOUND --
    which is static for the whole run, making this exact rather than an
    approximation.

    ``NETSEG > 1`` (segmented ET) cannot use the array-based ``EVTA`` at all,
    because ``READASARRAYS`` and segments are mutually exclusive in MODFLOW 6.
    The package therefore becomes list-based, one record per active column per
    period, which is a large but honest translation of the same physics.
    """

    ets = model.ets
    columns = np.flatnonzero(model.has_active_column)
    layers = model.uppermost_active[columns]
    nseg = max(int(ets.netseg), 1)

    def latest(store: dict[int, Any], period: int, default=None):
        """Value for ``period``, falling back to the most recent one defined."""

        if period in store:
            return store[period]
        earlier = [p for p in store if p < period]
        return store[max(earlier)] if earlier else default

    data: dict[int, list[Any]] = {}
    for period in range(model.nper):
        rate = latest(ets.rate, period)
        if rate is None:
            continue
        surf = latest(ets.surf, period, np.zeros(model.ncpl))
        depth = latest(ets.depth, period, np.ones(model.ncpl))
        pxdp = latest(ets.pxdp, period, [])
        petm = latest(ets.petm, period, [])
        rows = []
        for layer, cell in zip(layers, columns, strict=True):
            record: list[Any] = [
                (int(layer), int(cell)),
                float(surf[cell]),
                float(rate[cell]),
                float(depth[cell]),
            ]
            record += [float(a[cell]) for a in pxdp[: nseg - 1]]
            record += [float(a[cell]) for a in petm[: nseg - 1]]
            rows.append(record)
        data[period] = rows

    logger.info(
        "ETS -> list-based EVT: %d periods x %d columns = %d records (nseg=%d)",
        len(data),
        columns.size,
        sum(len(v) for v in data.values()),
        nseg,
    )
    return api.evt.flopy(stress_period_data=data, nseg=nseg, name="evt")


def _hfb_spec(model: UsgModel) -> PackageSpec | None:
    """Build MODFLOW 6 ``HFB`` from USG node pairs.

    Vertical barriers ARE legal on DISV as of MODFLOW 6 6.7.0, so a cross-layer
    pair is no longer dropped on sight -- it is handed to the same validator every
    other route uses, which accepts a pair in adjacent layers and rejects one that
    skips a layer. (FloPy 3.10's embedded definition still asserts the old
    same-layer rule and is out of date.)

    ``mf.hfb`` also dedupes: MODFLOW-USG tolerates a repeated barrier face, while
    MODFLOW 6 silently applies its series formula twice and then leaves that face
    permanently wrong.
    """

    if not len(model.hfb):
        return None

    first = model.to_cellid(model.hfb[:, 0].astype(np.int64))
    second = model.to_cellid(model.hfb[:, 1].astype(np.int64))
    on_grid = (first[:, 0] >= 0) & (second[:, 0] >= 0)
    if not on_grid.any():
        logger.warning("every HFB barrier references a CLN node; HFB omitted")
        return None
    if not on_grid.all():
        logger.warning(
            "%d of %d HFB barriers reference CLN nodes and were dropped with the CLN network",
            int((~on_grid).sum()),
            on_grid.size,
        )

    pairs = [
        ((int(a[0]), int(a[1])), (int(b[0]), int(b[1])))
        for a, b in zip(first[on_grid], second[on_grid], strict=True)
    ]
    return api.hfb(
        pairs=pairs,
        hydchr=[float(v) for v in model.hfb[on_grid, 2]],
        context=ModelContext(grid=model.grid),
        name="hfb",
    )


def _ims_spec(model: UsgModel, *, name: str, complexity: str) -> PackageSpec:
    """Translate the USG ``SMS`` solver onto MODFLOW 6 ``IMS``.

    Only the convergence criteria and iteration caps carry over literally. The
    rest of SMS -- its under-relaxation flavour, backtracking and linear
    acceleration -- has no term-for-term MODFLOW 6 equivalent, so a complexity
    preset supplies those and the result is a *comparable*, not identical,
    solver.

    The preset defaults to ``COMPLEX`` rather than ``MODERATE`` deliberately. A
    USG model that reached for SMS with a Newton formulation is a hard nonlinear
    problem -- convertible layers, wetting and drying, thin cells -- and on the
    Ten Trails model ``MODERATE`` does not merely converge poorly: MODFLOW 6
    dies with SIGFPE during the first solve. ``COMPLEX`` runs the same model to
    normal termination, with or without the SMS-derived tolerances.
    """

    sms = model.sms
    if sms is None:
        return api.ims(models=(name,), complexity=complexity)
    # Carry only what means the same thing in both solvers: the OUTER (nonlinear)
    # head-change closure and the iteration caps. SMS's HICLOSE is not MODFLOW 6's
    # INNER_DVCLOSE -- the inner solvers differ, and SMS's remaining inner controls
    # (IACL, NORDER, LEVEL, NORTH, RCLOSEPCGU) have no MODFLOW 6 counterpart at
    # all. Copying half of SMS's numbers on top of half a preset produces a solver
    # that is neither, so the preset owns the inner solve outright.
    values: dict[str, Any] = {
        "complexity": complexity,
        "outer_dvclose": sms.hclose,
        "outer_maximum": sms.mxiter,
        "inner_maximum": sms.iter1,
    }
    if sms.linear_acceleration:
        values["linear_acceleration"] = sms.linear_acceleration
    if sms.is_newton and sms.theta is not None:
        # SMS's delta-bar-delta under-relaxation and backtracking map field for
        # field onto MODFLOW 6's, and they are not decoration: they are how the
        # original model was made to converge. Dropping them and keeping only the
        # tolerances leaves MODFLOW 6 taking 100,000-foot Newton steps on a model
        # whose heads span 400 feet.
        values.update(
            under_relaxation="DBD",
            under_relaxation_theta=sms.theta,
            under_relaxation_kappa=sms.akappa,
            under_relaxation_gamma=sms.gamma,
            under_relaxation_momentum=sms.amomentum,
        )
        if sms.numtrack:
            values.update(
                backtracking_number=sms.numtrack,
                backtracking_tolerance=sms.btol,
                backtracking_reduction_factor=sms.breduc,
                backtracking_residual_limit=sms.reslim,
            )
    return api.ims(models=(name,), **values)


# ---------------------------------------------------------------- report


def conversion_notes(model: UsgModel) -> str:
    """Return the human-readable account of the conversion. See :meth:`UsgModel.report`."""

    lines: list[str] = []
    add = lines.append

    add(f"MODFLOW-USG model: {model.name_file.path.name}")
    add(f"  grid       {model.nlay} layers x {model.ncpl} cells = {model.nodes} nodes")
    add(f"  time       {model.nper} stress periods, {model.disu.perlen.sum():.0f} {model.time_units}")
    add(f"  units      length={model.length_units}, time={model.time_units}")
    if model.length_units == "unknown":
        add("             ^ DISU LENUNI is 0; MODFLOW 6 will be told 'unknown' unless you say otherwise")
    if model.start_date_time:
        add(f"  starts     {model.start_date_time} (inferred from the DISU period labels)")
    active = int((model.idomain != 0).sum())
    add(f"  active     {active} of {model.idomain.size} cells")
    add("")

    add("CONVERTED")
    add("  DISU  -> DISV      layered grid confirmed; geometry from the .gsf")
    add(f"  BAS6  -> IDOMAIN+IC  starting heads {model.strt.min():.1f}-{model.strt.max():.1f}")
    convertible = "convertible" if model.icelltype.any() else "confined"
    add(f"  LPF   -> NPF+STO   icelltype {convertible}; K {model.k.min():.3g}-{model.k.max():.3g}")
    if model.sms:
        add("  SMS   -> IMS       approximate: criteria carry over, solver internals do not")
    for ftype in ("CHD", "DRN", "GHB", "WEL", "RIV"):
        records = model.boundaries.get(ftype)
        if records is None:
            continue
        if records.n_records == 0 and records.n_cln_records:
            add(f"  {ftype:<5} -> (none)     every record is on a CLN node; see below")
            continue
        note = f"{records.n_defined} period(s) defined, {len(records.reused)} reused"
        add(f"  {ftype:<5} -> {ftype:<9} {records.n_records} records, {note}")
    if model.rch is not None:
        add(f"  RCH   -> RCHA      {len(model.rch.rech)} arrays, NRCHOP={model.rch.nrchop}")
    if model.ets is not None:
        add(f"  ETS   -> EVT       segmented (NETSEG={model.ets.netseg}), list-based")
    if model.hfb.size:
        add(f"  HFB6  -> HFB       {len(model.hfb)} barriers")
    add("")

    findings = mf6_findings(model)
    if findings:
        add("MODFLOW 6 IS STRICTER")
        for finding in findings:
            marker = "REJECTS" if finding.blocks_run else "allows "
            add(f"  {marker}  {finding.detail}")
            add(f"           fix: {finding.fix}")
        add("")

    approximations, omissions = _caveats(model)
    if approximations:
        add("APPROXIMATED")
        for line in approximations:
            add(f"  - {line}")
        add("")
    if omissions:
        add("NOT CONVERTED")
        for line in omissions:
            add(f"  - {line}")
        add("")
    return "\n".join(lines)


def _caveats(model: UsgModel) -> tuple[list[str], list[str]]:
    """Return the approximations and the omissions, each as report lines."""

    approximations: list[str] = []
    omissions: list[str] = []

    chd = model.boundaries.get("CHD")
    if chd is not None and chd.periods:
        drift = [
            float(np.abs(block[:, 2] - block[:, 1]).max())
            for block in chd.periods.values()
            if block.shape[1] > 2 and len(block)
        ]
        if drift and max(drift) > 0:
            approximations.append(
                f"CHD ramps from shead to ehead within each period (max {max(drift):.3f} "
                "of head); MODFLOW 6 holds one value, so the end-of-period head is used"
            )

    if model.ets is not None:
        if model.ets.netsop == 3:
            approximations.append(
                "ETS applied to the highest ACTIVE cell (NETSOP=3), which MODFLOW 6 cannot "
                "express; resolved once from IBOUND, which is static -- so this is exact"
            )
        if model.ets.netseg > 1:
            columns = int(model.has_active_column.sum())
            approximations.append(
                f"ETS is segmented (NETSEG={model.ets.netseg}), and MODFLOW 6 cannot combine "
                f"segments with array input, so EVT becomes list-based: about "
                f"{columns * model.nper:,} records"
            )

    if model.sms is not None:
        approximations.append(
            "SMS -> IMS keeps the convergence criteria; under-relaxation, backtracking and "
            "linear acceleration have no term-for-term equivalent and come from the preset"
        )

    if "DRAWDOWN SAVE" in model.oc_words:
        omissions.append("OC requested DRAWDOWN output, which MODFLOW 6 does not produce")

    if model.cln is not None:
        cells = model.cln.cells_touched(model.ncpl)
        bodies = model.cln.waterbodies
        streams = model.cln.streams
        omissions.append(
            f"CLN ({model.cln.nclnnds} nodes) has no MODFLOW 6 counterpart and is NOT "
            f"converted. It touches {cells.size} groundwater cells in "
            f"{len(model.cln.features)} features "
            f"({len(bodies)} waterbody, {len(streams)} stream). "
            "Rebuild them as LAK/SFR from model.cln_polygons()."
        )
        for feature in model.cln.features:
            omissions.append(f"    [{feature.index}] {feature.describe()}")

    for ftype, records in model.boundaries.items():
        if records.n_cln_records:
            omissions.append(
                f"{ftype}: {records.n_cln_records} records address CLN nodes, not groundwater "
                f"cells, and are dropped with the CLN network"
            )
    return approximations, omissions
