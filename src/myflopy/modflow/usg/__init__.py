"""Read a MODFLOW-USG model and convert it to MODFLOW 6.

MODFLOW-USG models are the ones myflopy could not previously reach at all: its
``DISU`` carries connectivity and geometric measures but no coordinates, its
grid lives in a separate ``.gsf``, and packages it supports (``CLN``, ``ETS``)
have no MODFLOW 6 counterpart.

Start with :func:`~myflopy.modflow.usg.reader.read_usg`::

    usg = mf.read_usg("flow-tt01_USE.nam", gsf="flow-tt01.gsf")
    print(usg.report())          # what mapped, what did not, and where
    sim = usg.to_mf6("tentrails")   # -> SimulationSpec, ready for Project.prepare_run
"""

from __future__ import annotations

from myflopy.modflow.usg.model import UsgModel
from myflopy.modflow.usg.reader import read_usg

__all__ = ["UsgModel", "read_usg"]
