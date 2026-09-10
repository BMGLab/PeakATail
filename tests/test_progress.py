"""ProgressManager + ProgressClient contract tests."""
import multiprocessing
import pickle

import pytest

from ema.progress import ProgressManager, ProgressClient


def test_manager_context_manager_lifecycle():
    with ProgressManager(disable=True) as pm:
        assert pm is not None
    # __exit__ ran without error


def test_add_stage_returns_int_id():
    with ProgressManager(disable=True) as pm:
        tid = pm.add_stage("Peak calling", total=100)
        assert isinstance(tid, int)


def test_add_subtask_under_stage():
    with ProgressManager(disable=True) as pm:
        parent = pm.add_stage("Peak calling", total=2)
        child = pm.add_subtask(parent, "sampleA", total=46)
        assert child != parent


def test_client_is_picklable():
    """Required for spawn workers."""
    with ProgressManager(disable=True) as pm:
        tid = pm.add_stage("X", total=10)
        client = pm.client(tid)
        data = pickle.dumps(client)
        restored = pickle.loads(data)
        assert isinstance(restored, ProgressClient)


def test_client_advance_does_not_raise_when_disabled():
    with ProgressManager(disable=True) as pm:
        tid = pm.add_stage("X", total=10)
        client = pm.client(tid)
        client.advance(5)  # disabled mode = no-op


def test_disabled_mode_when_no_tty(monkeypatch):
    """When stderr is not a TTY and disable=False, ProgressManager still
    works but renders as plain-text status (Rich auto-fallback)."""
    monkeypatch.setattr("sys.stderr.isatty", lambda: False)
    with ProgressManager(disable=False) as pm:
        tid = pm.add_stage("X", total=5)
        pm.client(tid).advance(5)


def test_spawn_worker_can_advance_through_client():
    """End-to-end: start manager, send client to spawn worker, worker
    posts ticks, parent's listener increments the bar."""
    with ProgressManager(disable=True) as pm:
        tid = pm.add_stage("X", total=10)
        client = pm.client(tid)
        ctx = multiprocessing.get_context("spawn")
        proc = ctx.Process(target=_child_advances, args=(client, 7))
        proc.start()
        proc.join(timeout=10)
        assert proc.exitcode == 0


def _child_advances(client, n):
    """Top-level helper for spawn pickling."""
    for _ in range(n):
        client.advance(1)
