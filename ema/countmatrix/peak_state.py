"""Per-run peak-calling state, encapsulated to avoid module-level globals.

Replaces Peak.pasnumber class attribute. Each peak_calling() invocation
gets its own PeakCallingState instance -- no implicit sharing across
parallel workers.

These module-level functions wrap a process-global BarcodeIndex singleton.
They work correctly under spawn-based multiprocessing (each worker has its
own copy) but ARE NOT THREAD-SAFE. New code should create explicit
BarcodeIndex instances and pass them through the call chain.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class PeakCallingState:
    """Per-invocation state for one peak_calling() call.

    pasnumber starts at 0; bump_pasnumber() returns the next ID and
    increments the counter atomically (within one process; for multi-process
    use one instance per worker).

    Attributes:
        pasnumber: Current peak counter. Starts at 0; after bump_pasnumber()
            returns 1, pasnumber is 1.
    """

    pasnumber: int = 0

    def bump_pasnumber(self) -> int:
        """Increment the counter and return the new (post-increment) value.

        Returns:
            The new pasnumber after incrementing (1-based).
        """
        self.pasnumber += 1
        return self.pasnumber

    def reset(self) -> None:
        """Reset the counter back to 0."""
        self.pasnumber = 0
