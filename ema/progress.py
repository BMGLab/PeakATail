"""Spawn-safe Rich progress wrapper for PeakATail.

Single Progress instance lives in the main process. Workers receive a
ProgressClient (picklable) that posts tick messages over a multiprocessing
queue. A listener thread in the main process drains the queue and calls
progress.advance() on the right task.
"""
from __future__ import annotations

import multiprocessing
import sys
import threading
from dataclasses import dataclass
from typing import Optional

from rich.console import Console
from rich.progress import (
    Progress,
    BarColumn,
    ProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    MofNCompleteColumn,
    SpinnerColumn,
)
from rich.text import Text


class ThickBarColumn(ProgressColumn):
    """Bar column with a medium-weight glyph.

    Filled: U+2586 LOWER THREE QUARTERS BLOCK (75 % cell height) — visibly
    chunkier than Rich's default U+2501 (heavy horizontal, ~50 %) but not
    as overwhelming as U+2588 FULL BLOCK.  Empty: U+2591 LIGHT SHADE.
    """

    FILLED_CHAR = "▆"
    EMPTY_CHAR = "░"

    def __init__(self, bar_width: int | None = None, reserve: int = 35) -> None:
        super().__init__()
        self.bar_width = bar_width
        self.reserve = reserve  # cells left for the other columns

    def render(self, task) -> Text:
        # Auto-width to terminal if bar_width is None.
        if self.bar_width is None:
            try:
                width = self._table.console.width - self.reserve  # type: ignore[attr-defined]
            except Exception:
                width = 60
            width = max(20, width)
        else:
            width = self.bar_width
        total = task.total or 0
        if total <= 0:
            return Text(self.EMPTY_CHAR * width, style="bar.back")
        completed = min(task.completed, total)
        filled = int(width * completed / total)
        return Text(
            self.FILLED_CHAR * filled + self.EMPTY_CHAR * (width - filled),
            style="bar.complete",
        )


_TICK_ADVANCE = "advance"
_TICK_TOTAL = "set_total"
_TICK_STOP = "stop"


@dataclass
class ProgressClient:
    """Picklable handle. Workers call .advance(n) / .set_total(n)."""
    queue: multiprocessing.Queue
    task_id: int

    def advance(self, n: int = 1) -> None:
        try:
            self.queue.put((_TICK_ADVANCE, self.task_id, n))
        except Exception:
            pass  # never crash a worker over a tick

    def set_total(self, n: int) -> None:
        try:
            self.queue.put((_TICK_TOTAL, self.task_id, n))
        except Exception:
            pass


class ProgressManager:
    """Owns the Rich Progress instance + the tick queue listener.

    Use as a context manager:
        with ProgressManager() as pm:
            tid = pm.add_stage("Peak calling", total=2)
            sub = pm.add_subtask(tid, "sampleA", total=46)
            client = pm.client(sub)
            # ... pass client to a spawn worker
    """

    def __init__(self, disable: bool = False) -> None:
        self.disable = disable
        self._console = Console(file=sys.stderr)
        self._progress: Optional[Progress] = None
        self._queue: Optional[multiprocessing.Queue] = None
        self._listener: Optional[threading.Thread] = None
        self._stop_event = threading.Event()
        self._task_meta: dict[int, dict] = {}

    def __enter__(self) -> "ProgressManager":
        if self.disable:
            return self
        self._progress = Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            ThickBarColumn(bar_width=None),
            MofNCompleteColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=self._console,
            transient=False,
        )
        self._progress.start()
        # Use a Manager queue so it's spawn-shareable
        self._queue = multiprocessing.Manager().Queue()
        self._listener = threading.Thread(target=self._drain_queue, daemon=True)
        self._listener.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        try:
            if self._queue is not None:
                self._queue.put((_TICK_STOP, 0, 0))
            if self._listener is not None:
                self._listener.join(timeout=2)
        finally:
            if self._progress is not None:
                self._progress.stop()

    def add_stage(self, name: str, total: Optional[int] = None) -> int:
        if self.disable or self._progress is None:
            tid = len(self._task_meta) + 1
            self._task_meta[tid] = {"name": name, "total": total}
            return tid
        tid = self._progress.add_task(name, total=total)
        self._task_meta[int(tid)] = {"name": name, "total": total}
        return int(tid)

    def add_subtask(self, parent: int, name: str, total: Optional[int] = None) -> int:
        # Rich doesn't have nested tasks natively — render as indented child task.
        return self.add_stage(f"  └ {name}", total=total)

    def client(self, task_id: int) -> ProgressClient:
        # Even when disabled we still hand out a Client (no-op via empty queue).
        if self._queue is None:
            # Disabled or pre-enter: build a dummy queue. .put() succeeds, never read.
            self._queue = multiprocessing.Manager().Queue()
        return ProgressClient(queue=self._queue, task_id=task_id)

    def finish(self, task_id: int) -> None:
        """Force a stage bar to 100 % regardless of advances received.

        Use when the apriori ``total`` is an over-estimate — e.g. peak
        calling counts BAM references with mapped reads but inner filters
        drop a chunk, so the advance count never reaches the estimate.
        """
        if self.disable or self._progress is None:
            return
        try:
            task = self._progress.tasks[task_id]
            if task.total is not None:
                self._progress.update(task_id, completed=task.total)
        except Exception:
            pass

    # ------------------------------------------------------------------ #
    # internals                                                           #
    # ------------------------------------------------------------------ #
    def _drain_queue(self) -> None:
        assert self._queue is not None
        while not self._stop_event.is_set():
            try:
                msg = self._queue.get(timeout=0.2)
            except Exception:
                continue
            kind, tid, val = msg
            if kind == _TICK_STOP:
                self._stop_event.set()
                break
            if self._progress is None:
                continue
            try:
                if kind == _TICK_ADVANCE:
                    self._progress.advance(tid, val)
                elif kind == _TICK_TOTAL:
                    self._progress.update(tid, total=val)
            except Exception:
                pass
