"""A transport (GWT) sibling for the canonical valley model.

Plan §6.1/6.2 item 6. The canonical model is the one model every notebook and
test shares; this attaches solute transport to it rather than introducing a
second model family.

**Why a sibling on the existing simulation, not a ``Project`` spec.** Ledger 56
recorded the transport fixture as blocked on "the multi-model spec path". That
was measured wrong twice over (2026-07-28):

* ``SimulationBase`` already owns a real ``MFSimulation``; a GWT model, a second
  IMS and a ``GWF6-GWT6`` exchange attach to it with no library change, and
  ``CANONICAL_MODEL_CONTRACT.validate()`` is entirely ``model.gwf``-scoped so it
  keeps passing unchanged.
* The spec path would actively BREAK calibration. ``Project.prepare_run`` gives
  each model its own subdirectory, so external arrays land at
  ``<ws>/flow/flow.npf_k.txt`` -- while PEST's parameter-file resolution globs
  the workspace ROOT. A flat sibling keeps every existing PEST plumbing working.
  (Binary *outputs* stay at the root either way, which is why observations alone
  would have survived the spec path.)

The contaminant source is a CNC package on the transport model itself, so the
flow model is untouched -- no auxiliary concentration columns threaded through
CHD, and the canonical contract sees exactly the model it always saw.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import flopy
import numpy as np

from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    build_canonical_model,
)
from myflopy.modflow.mf6.simulation.base import SimulationBase

#: Suffix appended to the flow model's name to name its transport sibling.
TRANSPORT_NAME_SUFFIX = "_t"

#: MF6 rejects a MODELNAME longer than this ("exceeds maximum allowed length").
#: The canonical flow model is already 14 characters, so the sibling's name has
#: to be built to fit rather than simply concatenated.
MAX_MF6_MODEL_NAME = 16


def transport_model_name(flow_name: str, suffix: str = TRANSPORT_NAME_SUFFIX) -> str:
    """Return a valid GWT sibling name for ``flow_name``, truncated if needed."""

    suffix = str(suffix)
    keep = MAX_MF6_MODEL_NAME - len(suffix)
    if keep < 1:
        raise ValueError(f"Transport name suffix {suffix!r} leaves no room for a name.")
    return f"{str(flow_name)[:keep]}{suffix}"


@dataclass
class CanonicalTransportModel:
    """A canonical flow model plus its attached GWT sibling.

    Attributes
    ----------
    model
        The canonical GWF model, exactly as :func:`build_canonical_model`
        returns it -- the contract still validates, and ``model.gwf`` is still
        the flow model.
    transport_name
        Name of the GWT model inside the same simulation. Read its results with
        ``transport_view()``.
    source_cells
        Zero-based cells carrying the constant-concentration source.
    source_concentration
        Concentration held at those cells.
    """

    model: SimulationBase
    transport_name: str
    source_cells: tuple[int, ...]
    source_concentration: float

    @property
    def workspace(self) -> Path:
        """The (flat) simulation workspace holding both models."""

        return Path(self.model.sim.sim_path)

    def transport_view(self):
        """Return a kind-aware view of the GWT sibling (``.conc``, ``.budget``, ...).

        ``SimulationBase`` exposes only its own ``model.gwf``, so the sibling is
        reached by reopening the written workspace, which
        :func:`~myflopy.project.run_model.load_mf6_run` resolves by model name.
        """

        from myflopy.project.run_model import load_mf6_run

        return load_mf6_run(
            self.workspace, model_name=self.transport_name, verbosity_level=0
        )


def _disv_gridprops_from(model: SimulationBase) -> dict:
    """Read the flow model's DISV mesh back out, to build the sibling on it.

    The two models must share one grid: MF6 requires it for a GWF-GWT exchange,
    and every myflopy view assumes one ``vor`` behind both.
    """

    disv = model.gwf.disv
    return {
        "nlay": int(disv.nlay.get_data()),
        "ncpl": int(disv.ncpl.get_data()),
        "nvert": int(disv.nvert.get_data()),
        "top": disv.top.get_data(),
        "botm": disv.botm.get_data(),
        "vertices": disv.vertices.get_data(),
        "cell2d": disv.cell2d.get_data(),
        "idomain": disv.idomain.get_data(),
    }


def _upgradient_source_cells(model: SimulationBase, count: int) -> tuple[int, ...]:
    """Pick source cells at the up-valley (high-head) end of the domain.

    Chosen by normalized position rather than by a named region so the source
    sits where the plume has the whole valley to travel down -- a source in the
    discharge zone would leave nothing to observe.
    """

    centers = np.asarray(model.vor.points, dtype=float)
    xs = centers[:, 0]
    span = xs.max() - xs.min()
    upgradient = np.where(xs <= xs.min() + 0.18 * span)[0]
    if upgradient.size == 0:  # pragma: no cover - defensive
        upgradient = np.argsort(xs)[:count]
    # spread them across the inflow face rather than clustering
    ys = centers[upgradient, 1]
    ordered = upgradient[np.argsort(ys)]
    step = max(1, len(ordered) // max(count, 1))
    return tuple(int(c) for c in ordered[::step][:count])


def attach_transport_model(
    model: SimulationBase,
    *,
    porosity: float = 0.25,
    longitudinal_dispersivity: float = 10.0,
    source_cells: tuple[int, ...] | list[int] | None = None,
    n_source_cells: int = 6,
    source_concentration: float = 1.0,
    name: str | None = None,
) -> CanonicalTransportModel:
    """Attach a GWT sibling to ``model``'s simulation and return both.

    The simulation stays FLAT (one directory, both models), which is what keeps
    PEST's root-globbing parameter resolution working -- see the module
    docstring.

    Parameters
    ----------
    porosity
        Mobile porosity for MST. The single strongest control on plume timing.
    longitudinal_dispersivity
        DSP ``alh``; transverse is a tenth of it, the usual convention.
    source_cells, n_source_cells
        Zero-based layer-0 cells held at ``source_concentration`` by a CNC
        package. Defaults to ``n_source_cells`` cells across the up-valley face.
    source_concentration
        Concentration held at the source cells (relative units).
    """

    simulation = model.sim
    flow_name = model.name
    transport_name = name or transport_model_name(flow_name)
    grid = _disv_gridprops_from(model)

    cells = (
        _upgradient_source_cells(model, n_source_cells)
        if source_cells is None
        else tuple(int(c) for c in source_cells)
    )

    gwt = flopy.mf6.ModflowGwt(
        simulation,
        modelname=transport_name,
        # save_flows on the MODEL is what makes MF6 write a non-empty .cbc for a
        # transport model; OC asking for a budget is not sufficient (ledger 90).
        save_flows=True,
    )
    flopy.mf6.ModflowGwtdisv(
        gwt,
        nlay=grid["nlay"],
        ncpl=grid["ncpl"],
        nvert=grid["nvert"],
        top=grid["top"],
        botm=grid["botm"],
        vertices=grid["vertices"],
        cell2d=grid["cell2d"],
        idomain=grid["idomain"],
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtdsp(
        gwt,
        alh=float(longitudinal_dispersivity),
        ath1=float(longitudinal_dispersivity) / 10.0,
    )
    flopy.mf6.ModflowGwtmst(gwt, porosity=float(porosity))

    # SSM is MANDATORY once the flow model has any boundary package ("Flow model
    # has boundary packages, but there is no SSM package"). `sources=None` means
    # no package carries an auxiliary concentration, so every boundary inflow
    # enters at concentration zero -- clean water, with the CNC below as the only
    # source. That is exactly the intended conceptual model.
    flopy.mf6.ModflowGwtssm(gwt, sources=None, pname="ssm")

    # Every ADVANCED flow package needs its transport counterpart, or MF6 refuses
    # to run ("GWF water mover is active but the GWT MVT package has not been
    # specified"). Each is keyed to its flow package by name and starts clean.
    flow_packages = {str(name).lower() for name in model.gwf.package_names}
    if "lak" in flow_packages:
        flopy.mf6.ModflowGwtlkt(
            gwt, flow_package_name="lak", boundnames=False, pname="lkt",
            packagedata=[(i, 0.0) for i in range(int(model.lak.nlakes.get_data()))],
        )
    if "sfr" in flow_packages:
        flopy.mf6.ModflowGwtsft(
            gwt, flow_package_name="sfr", boundnames=False, pname="sft",
            packagedata=[(i, 0.0) for i in range(int(model.sfr.nreaches.get_data()))],
        )
    if "uzf" in flow_packages:
        flopy.mf6.ModflowGwtuzt(
            gwt, flow_package_name="uzf", boundnames=False, pname="uzt",
            packagedata=[(i, 0.0) for i in range(int(model.uzf.nuzfcells.get_data()))],
        )
    if "mvr" in flow_packages:
        # MVT carries solute along the SFR->LAK mover connections; it needs LKT
        # and SFT to exist, which is why it is registered last.
        flopy.mf6.ModflowGwtmvt(gwt, pname="mvt")
    flopy.mf6.ModflowGwtcnc(
        gwt,
        stress_period_data={
            0: [[(0, int(cell)), float(source_concentration)] for cell in cells]
        },
        pname="cnc",
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{transport_name}.ucn",
        budget_filerecord=f"{transport_name}.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    # A transport model needs its own solver: the two are solved sequentially,
    # and BICGSTAB is the usual choice for the non-symmetric transport matrix.
    # The EXPLICIT filename is load-bearing, not cosmetic. Without it FloPy gives
    # the new solver the "<sim>_0.ims" auto-name, which lands it BEFORE the flow
    # model's solver in mfsim.nam -- and MF6 rejects that outright: "the IMS
    # specified for GWF must be listed in mfsim.nam before the IMS for GWT".
    # Registration order does not fix it; the filename does.
    transport_ims = flopy.mf6.ModflowIms(
        simulation,
        pname="ims_gwt",
        filename=f"{transport_name}.ims",
        complexity="SIMPLE",
        linear_acceleration="BICGSTAB",
        outer_maximum=100,
        inner_maximum=200,
    )
    simulation.register_solution_package(transport_ims, [transport_name])

    flopy.mf6.ModflowGwfgwt(
        simulation,
        exgtype="GWF6-GWT6",
        exgmnamea=flow_name,
        exgmnameb=transport_name,
        filename=f"{flow_name}_{transport_name}.gwfgwt",
    )

    return CanonicalTransportModel(
        model=model,
        transport_name=transport_name,
        source_cells=cells,
        source_concentration=float(source_concentration),
    )


def build_canonical_transport_model(
    workspace: str | Path,
    *,
    config: CanonicalModelConfig | None = None,
    run: bool = True,
    **transport_kwargs,
) -> CanonicalTransportModel:
    """Build the canonical valley model with a solute-transport sibling.

    Parameters
    ----------
    workspace
        Directory for the (flat) simulation.
    config
        Canonical profile; defaults to :meth:`CanonicalModelConfig.validation`.
    run
        Run MF6 after building. Set ``False`` to inspect or modify first.
    **transport_kwargs
        Passed to :func:`attach_transport_model` (``porosity``,
        ``longitudinal_dispersivity``, ``source_cells``, ...).
    """

    config = config or CanonicalModelConfig.validation()
    built = attach_transport_model(
        build_canonical_model(Path(workspace), config=config), **transport_kwargs
    )
    if run:
        success, report = built.model.run_simulation()
        if not success:
            raise RuntimeError(
                "Canonical transport model failed to run:\n" + "\n".join(report[-25:])
            )
    return built


__all__ = [
    "CanonicalTransportModel",
    "MAX_MF6_MODEL_NAME",
    "TRANSPORT_NAME_SUFFIX",
    "attach_transport_model",
    "build_canonical_transport_model",
    "transport_model_name",
]
