"""Where a warning is attributed, decided once.

A warning pinned on library code is one Python's stock filters suppress in a notebook
or a script, which makes it a warning nobody hears. This module holds the one helper
that finds the first frame outside ``myflopy``, so every caller attributes alike.

It lived in ``modflow/utils/datatypes/choros.py`` until ledger 288 needed it in
``geopackage.py`` as well; a second copy is what the shared-rules consolidation of
ledgers 210-212 exists to prevent.
"""

from __future__ import annotations


def stacklevel_outside_package() -> int:
    """The ``stacklevel`` that pins a warning on the first frame outside ``myflopy``.

    Counted from the CALLER of this helper. The number of frames differs per route --
    a deprecation reaches ``Choro`` through a different count from every scope, and a
    ``.gpkg`` miss is ~5 frames deep when resolved eagerly and 10+ and variable when a
    lazy declaration resolves it inside ``prepare_run`` -- so no constant is right.
    """

    import sys

    level = 1
    frame = sys._getframe(1)
    while frame is not None:
        module = frame.f_globals.get("__name__", "")
        if not (module == "myflopy" or module.startswith("myflopy.")):
            return level
        frame = frame.f_back
        level += 1
    return level
