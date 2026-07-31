"""The one logger factory, and the convention for reporting a swallowed error.

A library must never configure logging -- that is the application's call -- so
this module does exactly two things: attach a ``NullHandler`` to the ``myflopy``
root (which is what stops "No handlers could be found" noise in a program that
never opts in), and hand out child loggers.

Getting a logger::

    from myflopy._logging import get_logger

    logger = get_logger(__name__)

Turning it on, from a notebook or a script::

    import logging
    logging.basicConfig(level=logging.DEBUG)

**The swallow convention.** Plan 7.3 narrowed ~60 broad ``except Exception``
handlers down to the exceptions each guarded call can actually raise. Narrowing
alone is only half the fix: a fallback that fires silently is still invisible,
so a caller who gets a blank hover label or a default value has no way to learn
why. Every such handler logs at DEBUG, on one line, in the same shape::

    except (KeyError, AttributeError) as error:
        logger.debug("<what could not be done>; <what happens instead>: %s", error)

DEBUG rather than WARNING because these fallbacks are *expected* -- a partial
model legitimately has no TDIS, a static field legitimately has no time axis.
They are diagnostics for the person asking "why is this blank?", not events the
user needs to be told about. The message names the degradation, not just the
error, because ``KeyError: 'per'`` on its own explains nothing.
"""

from __future__ import annotations

import logging

__all__ = ["get_logger"]

#: The package root every myflopy logger hangs off.
ROOT_LOGGER_NAME = "myflopy"

logging.getLogger(ROOT_LOGGER_NAME).addHandler(logging.NullHandler())


def get_logger(name: str | None = None) -> logging.Logger:
    """Return the ``myflopy`` logger, or the child named for a module.

    Pass ``__name__`` from inside the package and you get the matching child
    (``myflopy.modflow.mf6.prt_maps``), so an application can raise or lower the
    volume of one subsystem. A name from outside the package -- or none at all --
    returns the package root rather than silently creating a sibling tree.
    """

    if not name or not name.startswith(ROOT_LOGGER_NAME):
        return logging.getLogger(ROOT_LOGGER_NAME)
    return logging.getLogger(name)
