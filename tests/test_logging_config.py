"""Tests for the rewritten logging_config — Rich console + queue listener."""
import logging
import multiprocessing
from pathlib import Path

import pytest

from ema.logging_config import setup_logging, teardown_logging, get_log_queue


@pytest.fixture(autouse=True)
def _clean_loggers():
    yield
    teardown_logging()


def test_setup_returns_log_path_when_output_dir_given(tmp_path: Path):
    log_path = setup_logging(level="INFO", output_dir=tmp_path)
    assert log_path is not None
    assert log_path.parent == tmp_path
    assert log_path.suffix == ".log"


def test_setup_returns_none_when_no_log_file(tmp_path: Path):
    log_path = setup_logging(level="INFO", output_dir=tmp_path, no_log_file=True)
    assert log_path is None


def test_quiet_overrides_to_warning(tmp_path: Path):
    setup_logging(level="INFO", output_dir=tmp_path, quiet=True)
    assert logging.getLogger("ema").level == logging.WARNING


def test_verbose_count_one_sets_debug_for_ema(tmp_path: Path):
    setup_logging(level="INFO", output_dir=tmp_path, verbose_count=1)
    assert logging.getLogger("ema").level == logging.DEBUG
    # third-party still at default WARNING (NOTSET=0 means inherits WARNING from root)
    assert logging.getLogger("pysam").getEffectiveLevel() >= logging.WARNING


def test_verbose_count_two_sets_debug_everywhere(tmp_path: Path):
    setup_logging(level="INFO", output_dir=tmp_path, verbose_count=2)
    assert logging.getLogger("ema").level == logging.DEBUG
    assert logging.getLogger().level == logging.DEBUG  # root


def test_per_logger_override(tmp_path: Path):
    setup_logging(
        level="INFO", output_dir=tmp_path,
        per_logger_overrides={"ema.countmatrix": "WARNING"},
    )
    assert logging.getLogger("ema.countmatrix").level == logging.WARNING
    assert logging.getLogger("ema").level == logging.INFO


def test_log_file_receives_records(tmp_path: Path):
    log_path = setup_logging(level="DEBUG", output_dir=tmp_path)
    log = logging.getLogger("ema.test_logging")
    log.info("hello-from-test")
    teardown_logging()  # flush handlers
    text = log_path.read_text()
    assert "hello-from-test" in text


def test_queue_handler_in_spawn_workers(tmp_path: Path):
    """A child spawned process can post a record through the QueueHandler
    and the parent's listener picks it up."""
    log_path = setup_logging(level="DEBUG", output_dir=tmp_path)
    from ema.logging_config import get_log_queue
    queue = get_log_queue()
    assert queue is not None

    ctx = multiprocessing.get_context("spawn")
    proc = ctx.Process(target=_child_emit_log, args=(queue,))
    proc.start()
    proc.join(timeout=10)
    assert proc.exitcode == 0
    teardown_logging()
    text = log_path.read_text()
    assert "from-spawn-child" in text


def _child_emit_log(queue):
    """Top-level (picklable) helper for the spawn worker test."""
    import logging
    from logging.handlers import QueueHandler
    h = QueueHandler(queue)
    logger = logging.getLogger("ema.test_child")
    logger.addHandler(h)
    logger.setLevel(logging.DEBUG)
    logger.info("from-spawn-child")


def test_setup_worker_logging_replaces_handlers(tmp_path: Path):
    log_path = setup_logging(level="DEBUG", output_dir=tmp_path)
    queue = get_log_queue()
    # Pretend we are in a worker now:
    from ema.logging_config import setup_worker_logging
    setup_worker_logging(queue)
    root = logging.getLogger()
    assert len(root.handlers) == 1
    assert root.handlers[0].__class__.__name__ == "QueueHandler"
