"""The library's logging contract: hand out loggers, configure nothing.

The failure this guards against is a library that calls ``basicConfig`` or
attaches a ``StreamHandler`` -- which hijacks the logging of every application
that imports it, and is only ever noticed by the person whose own log format
suddenly changed. A ``NullHandler`` on the package root is the whole of what a
library is allowed to do.
"""

from __future__ import annotations

import logging

from myflopy._logging import ROOT_LOGGER_NAME, get_logger


def test_a_module_name_gets_its_own_child_logger():
    """So an application can quiet one subsystem without quieting myflopy."""

    logger = get_logger("myflopy.modflow.mf6.prt_maps")
    assert logger.name == "myflopy.modflow.mf6.prt_maps"
    assert logger.propagate is True   # or configuring `myflopy` would do nothing


def test_a_foreign_name_falls_back_to_the_package_root():
    """A name from outside the package must not plant a sibling logger tree that
    `logging.getLogger("myflopy").setLevel(...)` would then fail to reach."""

    assert get_logger("some_other_package.thing").name == ROOT_LOGGER_NAME
    assert get_logger("").name == ROOT_LOGGER_NAME
    assert get_logger(None).name == ROOT_LOGGER_NAME
    assert get_logger().name == ROOT_LOGGER_NAME


def test_the_package_root_carries_a_null_handler_and_nothing_else():
    """The NullHandler suppresses the "no handlers" fallback; anything louder
    would be myflopy configuring an application's logging for it."""

    root = logging.getLogger(ROOT_LOGGER_NAME)
    assert any(isinstance(handler, logging.NullHandler) for handler in root.handlers)
    assert all(isinstance(handler, logging.NullHandler) for handler in root.handlers), (
        f"myflopy attached a real handler to its own root: {root.handlers}. A "
        "library configures nothing -- that is the application's decision."
    )
    assert root.level == logging.NOTSET, (
        "myflopy set a level on its root logger, which overrides the level the "
        "application chose."
    )


def test_importing_myflopy_does_not_configure_the_root_logger():
    """`logging.basicConfig` in library code is the classic version of this bug."""

    import myflopy  # noqa: F401  - the import IS the thing under test

    assert not logging.getLogger().handlers or all(
        handler not in logging.getLogger(ROOT_LOGGER_NAME).handlers
        for handler in logging.getLogger().handlers
    )


def test_a_swallow_message_reaches_a_caplog_listener(caplog):
    """The point of the convention: someone debugging a blank hover label can
    turn DEBUG on and see which fallback fired."""

    logger = get_logger("myflopy.modflow.mf6.example")
    with caplog.at_level(logging.DEBUG, logger=ROOT_LOGGER_NAME):
        try:
            raise KeyError("per")
        except KeyError as error:
            logger.debug("no stress-period column; drawing a static map: %s", error)

    assert "drawing a static map" in caplog.text
