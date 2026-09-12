"""The active domain of a model, read from whatever carries it.

MF6 keeps idomain on the discretization package, so a model loaded from disk
carries it there and nowhere else. Two callers need the same read -- a reopened
run's package views and observation targets -- and the first one's copy sat in
``myflopy.project.run_model``, the top import layer, where nothing else can
reach it. That is the whole reason this is its own module.
"""

from __future__ import annotations

import numpy as np


def idomain_from_flopy(gwf) -> list[int] | None:
    """Return a flattened idomain list from DISV or DISU packages when present."""

    for package_name in ("disu", "disv"):
        package = getattr(gwf, package_name, None)
        if package is None or getattr(package, "idomain", None) is None:
            continue
        values = np.asarray(package.idomain.array).reshape(-1)
        if values.size and all(value is not None for value in values):
            return values.astype(int).tolist()
    return None


def model_domain(model, ncpl: int | None = None):
    """The model's idomain as an array, or ``None`` when nothing declares one.

    Read in the order MF6 itself would: the discretization package first, since
    that is what MF6 reads and the only place a loaded model keeps it, then the
    myflopy context's ``domain``.

    ``SimulationBase.idomain`` is deliberately NOT consulted. That property
    BUILDS a domain from an ``idomain_path`` and raises ``TypeError`` when there
    is none, so asking it turns "no domain declared" into an error; a caller that
    swallowed the error would then treat an inactive cell as active. ``None`` here
    means "nothing says otherwise", which :func:`domain_active` reads as active.
    """

    gwf = getattr(model, "gwf", None)
    if gwf is None and (hasattr(model, "disv") or hasattr(model, "disu")):
        gwf = model
    flat = idomain_from_flopy(gwf) if gwf is not None else None
    if flat is not None:
        values = np.asarray(flat, dtype=int)
        return values.reshape(-1, ncpl) if ncpl else values
    context = getattr(model, "myflopy_context", None)
    domain = getattr(context, "domain", None) if context is not None else None
    return None if domain is None else np.asarray(domain)
