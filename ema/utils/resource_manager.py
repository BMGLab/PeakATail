"""Device-aware resource management for parallel work.

Computes safe `n_jobs` and batch sizes given the system's *current* free RAM
and CPU count, not the totals. Avoids over-subscribing when the box is already
loaded (e.g. by other heavy processes, or by an earlier stage of the pipeline).

Conservative defaults:
- max workers = min(physical_cpu - 2, free_ram_mb / per_worker_mb, hard_cap)
- never returns less than 1
- always leaves headroom (default 30% of free RAM)
"""

from __future__ import annotations
import os

try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False


class ResourceManager:
    """Singleton-ish accessor for resource decisions across the pipeline."""

    def __init__(
        self,
        ram_safety_fraction: float = 0.7,
        cpu_reserve: int = 2,
        hard_n_jobs_cap: int = 16,
    ):
        self.ram_safety_fraction = ram_safety_fraction
        self.cpu_reserve = cpu_reserve
        self.hard_n_jobs_cap = hard_n_jobs_cap

    # ---- system probes ----

    def free_ram_mb(self) -> int:
        if _HAS_PSUTIL:
            return int(psutil.virtual_memory().available // (1024 * 1024))
        # /proc/meminfo fallback
        try:
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemAvailable:"):
                        kb = int(line.split()[1])
                        return kb // 1024
        except OSError:
            pass
        return 4096  # 4 GB pessimistic fallback

    def physical_cpu_count(self) -> int:
        if _HAS_PSUTIL:
            n = psutil.cpu_count(logical=False)
            if n:
                return n
        return os.cpu_count() or 1

    def logical_cpu_count(self) -> int:
        return os.cpu_count() or 1

    # ---- decisions ----

    def get_n_jobs(self, per_worker_mb: int = 200) -> int:
        """How many parallel workers can we safely spin up?

        Args:
            per_worker_mb: Estimated peak RAM per worker, in MB.

        Returns:
            Worker count >= 1. Bounded by:
              - (free_ram * safety_fraction) / per_worker_mb
              - physical_cpu_count - cpu_reserve
              - hard_n_jobs_cap
        """
        ram_budget_mb = int(self.free_ram_mb() * self.ram_safety_fraction)
        ram_workers = max(1, ram_budget_mb // max(1, per_worker_mb))
        cpu_workers = max(1, self.physical_cpu_count() - self.cpu_reserve)
        return max(1, min(ram_workers, cpu_workers, self.hard_n_jobs_cap))

    def get_batch_size(
        self,
        per_item_mb: float,
        min_batch: int = 100,
        max_batch: int = 100_000,
    ) -> int:
        """Choose a batch size for streaming so a single batch fits in RAM budget."""
        ram_budget_mb = self.free_ram_mb() * self.ram_safety_fraction
        if per_item_mb <= 0:
            return max_batch
        n = int(ram_budget_mb / per_item_mb)
        return max(min_batch, min(n, max_batch))

    def report(self) -> dict:
        return {
            "free_ram_mb": self.free_ram_mb(),
            "physical_cpus": self.physical_cpu_count(),
            "logical_cpus": self.logical_cpu_count(),
            "ram_safety_fraction": self.ram_safety_fraction,
            "cpu_reserve": self.cpu_reserve,
            "hard_n_jobs_cap": self.hard_n_jobs_cap,
        }
