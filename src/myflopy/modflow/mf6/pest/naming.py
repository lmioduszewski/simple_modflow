"""The one rule that turns a label into a PEST name, for parameters and observations alike.

An observation group label (ledger 358) and a lake's or stream's boundname, made a
parameter name (``style="feature"``, ledger 558), are slugged by the
same function, :func:`pest_group_name`. It lives here, with no import of the rest of
the PEST engine, so the parameter side (``native_parameters``) and the observation
side (``observations``) both import it without either importing the other.
"""

from __future__ import annotations

import re

#: Every run of characters a PEST group name does not keep -- anything but a
#: lower-case letter, a digit, ``_``, ``.``, ``:`` and ``-`` -- becomes one ``_``.
#: PEST++ 5.2.16 refuses the version-2 control file over a comma or a double
#: quote in a group name, and whitespace splits a classic (version-1) one
#: (``docs/reproducers/pest_observation_groups.py``).
_GROUP_NAME_DROPS = re.compile(r"[^a-z0-9_.:\-]+")


def pest_group_name(label: str) -> str:
    """The PEST name a group label, a set prefix or a feature's boundname makes.

    Lower-cased (pyEMU lower-cases every name it reads), stripped, and every run
    of characters other than letters, digits, ``_``, ``.``, ``:`` and ``-``
    turned into one ``_`` -- so ``"North Wells"`` is ``north_wells``. No length
    limit is applied: PEST++ 5.2.16 read a 5000-character group name from the
    version-2 control file myflopy writes, and pyEMU's own per-column groups
    (``oname:hds_otype:lst_usecol:w1``) are already past classic PEST's 12.

    Parameters
    ----------
    label : str
        The text: a target's ``group`` label, an observation set's prefix, or a
        lake's or stream's boundname naming a ``style="feature"`` parameter. Any
        other object is read through ``str()``.

    Returns
    -------
    str
        The name. Two labels may give one name (``"North"``, ``"north"``); the
        observation groups (``assign_observation_groups``) and the feature
        parameters (``native_parameters``) each refuse that.
    """

    return _GROUP_NAME_DROPS.sub("_", str(label).strip().lower())
