"""Finding the sibling model a calibration target or observation lives on.

A ``PestProject`` hangs off the FLOW model -- that is where ``k``, ``recharge``
and the boundary packages are -- but a coupled simulation puts transport
properties (``mst.porosity``) and transport state (concentration) on a GWT
sibling. Both the parameter side (``native_parameters``) and the observation
side (``observations``) need the same lookup, and this module is deliberately
dependency-free so the low-layer parameter modules can import it without
dragging in the observation stack.
"""

from __future__ import annotations

__all__ = ["resolve_transport_model_name"]


def resolve_transport_model_name(project) -> str:
    """Find the GWT sibling of the model this calibration hangs off.

    A transport calibration estimates FLOW parameters (``k``/``recharge`` live on
    the GWF model) from CONCENTRATION data (which lives on the GWT model), so the
    project is built on the flow model and the transport sibling has to be found.
    The same lookup names the transport model's external input files when a
    TRANSPORT property (``porosity``) is the calibration target.
    """

    simulation = project.model.sim
    transport = [
        name
        for name in simulation.model_names
        if str(getattr(simulation.get_model(name), "model_type", "")).lower().startswith("gwt")
    ]
    if not transport:
        raise ValueError(
            "No GWT model found in this simulation, so there is nothing for a "
            "concentration observation to read or a transport parameter "
            "(e.g. porosity) to adjust. Attach one with "
            "`myflopy.modflow.mf6.canonical_transport.attach_transport_model`, or "
            "name it explicitly via ConcObservationSpec(transport_model_name=...)."
        )
    if len(transport) > 1:
        raise ValueError(
            f"Simulation has several GWT models ({sorted(transport)}); name the one "
            "to observe via ConcObservationSpec(transport_model_name=...)."
        )
    return str(transport[0])
