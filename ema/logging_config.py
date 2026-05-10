"""Centralised logging setup for PeakATail.

Uses stdlib `logging` + Rich for pretty console output + a per-run
timestamped file log + a multiprocessing queue listener so spawn workers
can route their records through the parent process's handlers.
"""
from __future__ import annotations

import logging
import multiprocessing
import sys
from datetime import datetime
from logging.handlers import QueueHandler, QueueListener
from pathlib import Path
from typing import Optional

from rich.logging import RichHandler

# Module-level globals (process-wide singletons)
_LOG_QUEUE: Optional[multiprocessing.Queue] = None
_QUEUE_LISTENER: Optional[QueueListener] = None
_FILE_HANDLER: Optional[logging.FileHandler] = None
_CONSOLE_HANDLER: Optional[RichHandler] = None


def get_log_queue() -> Optional[multiprocessing.Queue]:
    """Return the process-wide log queue (used by spawn workers).

    Workers should attach a `QueueHandler(get_log_queue())` to their loggers
    so records flow back to the main process's listener.
    """
    return _LOG_QUEUE


def setup_logging(
    level: str = "INFO",
    output_dir: Optional[Path] = None,
    quiet: bool = False,
    verbose_count: int = 0,
    no_log_file: bool = False,
    per_logger_overrides: Optional[dict[str, str]] = None,
) -> Optional[Path]:
    """Configure logging for the current process.

    Args:
        level: Default level for `ema.*` loggers (resolved before -v/-q apply).
        output_dir: Directory to write the per-run log file into.
        quiet: If True, force WARNING on `ema.*` regardless of `level`.
        verbose_count: -v count. 1 → DEBUG for ema.*. 2 → DEBUG everywhere.
        no_log_file: If True, no file handler is attached.
        per_logger_overrides: {"ema.countmatrix": "WARNING", ...} per-logger levels.

    Returns:
        The resolved log-file path, or None when no_log_file is True or
        output_dir is None.
    """
    global _LOG_QUEUE, _QUEUE_LISTENER, _FILE_HANDLER, _CONSOLE_HANDLER

    # Reset any previous state (test-safe).
    teardown_logging()

    ema_logger = logging.getLogger("ema")
    root_logger = logging.getLogger()

    # Resolve effective ema-level
    if quiet:
        ema_level = logging.WARNING
    elif verbose_count >= 1:
        ema_level = logging.DEBUG
    else:
        ema_level = getattr(logging, level.upper(), logging.INFO)

    ema_logger.setLevel(ema_level)
    ema_logger.handlers.clear()
    ema_logger.propagate = False  # don't double-emit via root

    # Root level: WARNING by default; DEBUG only when -vv
    root_logger.setLevel(logging.DEBUG if verbose_count >= 2 else logging.WARNING)
    root_logger.handlers.clear()

    # Console handler (Rich)
    _CONSOLE_HANDLER = RichHandler(
        show_time=True,
        show_level=True,
        show_path=False,
        markup=False,
        rich_tracebacks=True,
        tracebacks_suppress=[logging],
    )
    _CONSOLE_HANDLER.setLevel(ema_level)
    ema_logger.addHandler(_CONSOLE_HANDLER)
    root_logger.addHandler(_CONSOLE_HANDLER)

    # File handler
    log_path: Optional[Path] = None
    if output_dir is not None and not no_log_file:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        ts = datetime.now().strftime("%Y-%m-%d_%H%M%S")
        log_path = output_dir / f"peakatail_{ts}.log"
        _FILE_HANDLER = logging.FileHandler(log_path)
        _FILE_HANDLER.setLevel(logging.DEBUG)  # file always captures everything
        _FILE_HANDLER.setFormatter(
            logging.Formatter(
                "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        )
        ema_logger.addHandler(_FILE_HANDLER)
        root_logger.addHandler(_FILE_HANDLER)

    # Per-logger overrides
    for name, lvl in (per_logger_overrides or {}).items():
        logging.getLogger(name).setLevel(getattr(logging, lvl.upper(), logging.INFO))

    # Queue listener for spawn workers
    _LOG_QUEUE = multiprocessing.Manager().Queue()
    handlers = [h for h in (_CONSOLE_HANDLER, _FILE_HANDLER) if h is not None]
    _QUEUE_LISTENER = QueueListener(_LOG_QUEUE, *handlers, respect_handler_level=True)
    _QUEUE_LISTENER.start()

    return log_path


def teardown_logging() -> None:
    """Stop the queue listener and remove handlers. Idempotent. Safe in tests."""
    global _LOG_QUEUE, _QUEUE_LISTENER, _FILE_HANDLER, _CONSOLE_HANDLER
    if _QUEUE_LISTENER is not None:
        try:
            _QUEUE_LISTENER.stop()
        except Exception:
            pass
        _QUEUE_LISTENER = None
    _LOG_QUEUE = None
    for log in (logging.getLogger("ema"), logging.getLogger()):
        for h in list(log.handlers):
            try:
                h.flush()
                h.close()
            except Exception:
                pass
            log.removeHandler(h)
    _FILE_HANDLER = None
    _CONSOLE_HANDLER = None


def setup_worker_logging(queue: multiprocessing.Queue, level: str = "DEBUG") -> None:
    """Configure a spawn worker to forward all log records to the parent's queue.

    Call this once at the top of every worker function that runs under
    multiprocessing.get_context('spawn').Pool. After this call, any
    `logging.getLogger(__name__).info(...)` in the worker reaches the parent's
    handlers.

    Args:
        queue: The multiprocessing queue obtained from ``get_log_queue()``
            in the parent and passed to the worker as a JobSpec field.
        level: Worker-side floor; DEBUG is the right default so the parent
            handlers' own levels decide what gets shown.
    """
    h = QueueHandler(queue)
    h.setLevel(getattr(logging, level.upper(), logging.DEBUG))
    root = logging.getLogger()
    # Replace, don't append — workers always start fresh.
    for old in list(root.handlers):
        root.removeHandler(old)
    root.addHandler(h)
    root.setLevel(getattr(logging, level.upper(), logging.DEBUG))
