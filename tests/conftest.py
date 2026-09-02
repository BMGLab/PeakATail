"""Shared pytest fixtures.

The autouse fixture here exists because PeakATail configures logging
*globally*: ``ema.logging_config.setup_logging`` clears handlers on both the
``ema`` and root loggers and sets ``ema.propagate = False``.  That is correct
for the CLI and hostile to a test suite -- once one test calls it, every later
test in the same process inherits the mutated state, and anything listening on
the root logger (notably pytest's ``caplog``) goes silent.

Symptom this fixes: ``test_polya_two_tier_integration.py::
test_clip_rate_warning_fires_on_clip_free_bam`` passes in isolation and fails
in a full run, because ``test_logging_config.py`` sorts earlier and leaves
``propagate`` False behind.
"""
import logging

import pytest


@pytest.fixture(autouse=True)
def _restore_logging_state():
    """Snapshot and restore global logger state around every test."""
    loggers = (logging.getLogger("ema"), logging.getLogger())
    saved = [
        (lg, list(lg.handlers), lg.level, lg.propagate)
        for lg in loggers
    ]
    try:
        yield
    finally:
        for lg, handlers, level, propagate in saved:
            lg.handlers[:] = handlers
            lg.setLevel(level)
            lg.propagate = propagate
