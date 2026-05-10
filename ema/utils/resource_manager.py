"""Device-aware resource management for parallel work.

Computes safe `n_jobs` and batch sizes given the system's *current* free RAM
and CPU count, not the totals. Avoids over-subscribing when the box is already
loaded (e.g. by other heavy processes, or by an earlier stage of the pipeline).

Conservative defaults:
- max workers = min(physical_cpu - 2, free_ram_mb / per_worker_mb, hard_cap)
- never returns less than 1
- always leaves headroom (default 30% of free RAM)
- user_max_threads (from --threads CLI flag) is an absolute ceiling
"""

from __future__ import annotations
import logging
import os

logger = logging.getLogger(__name__)

try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False


class ResourceManager:
    """Process-wide authority for parallelism decisions across the pipeline.

    All parallel call sites must query this class via
    ``get_resource_manager().get_n_jobs(...)`` rather than calling
    ``os.cpu_count()`` or hardcoding ``n_jobs=-1``.

    Args:
        user_max_threads: Absolute ceiling on worker count imposed by the
            ``--threads`` CLI flag.  ``None`` means no user override.
        user_per_worker_mb: Optional per-worker RAM override (unused by
            default; reserved for future per-stage config).
        ram_safety_fraction: Fraction of free RAM that can be allocated to
            workers (default 0.7, i.e. keep 30% headroom).
        cpu_reserve: Number of logical CPUs to leave idle (default 2).
        hard_n_jobs_cap: Absolute maximum even if RAM/CPU budget would allow
            more (default 16).
    """

    def __init__(
        self,
        user_max_threads: int | None = None,
        user_per_worker_mb: int | None = None,
        ram_safety_fraction: float = 0.7,
        cpu_reserve: int = 2,
        hard_n_jobs_cap: int = 16,
    ) -> None:
        self.user_max_threads = user_max_threads
        self.user_per_worker_mb = user_per_worker_mb
        self.ram_safety_fraction = ram_safety_fraction
        self.cpu_reserve = cpu_reserve
        self.hard_n_jobs_cap = hard_n_jobs_cap

    # ---- system probes ----

    def free_ram_mb(self) -> int:
        """Return available RAM in MB (conservative, live reading)."""
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
        """Return the number of physical (non-hyperthreaded) CPU cores."""
        if _HAS_PSUTIL:
            n = psutil.cpu_count(logical=False)
            if n:
                return n
        return os.cpu_count() or 1

    def logical_cpu_count(self) -> int:
        """Return the number of logical CPU threads (including HT)."""
        return os.cpu_count() or 1

    # ---- decisions ----

    def _compute_budget(self, per_worker_mb: int) -> int:
        """Compute raw worker budget from RAM + CPU constraints.

        Does NOT apply the user_max_threads ceiling — that is done in
        ``get_n_jobs`` so callers that need the raw budget can access it.

        Args:
            per_worker_mb: Estimated peak RAM per worker, in MB.

        Returns:
            Worker count bounded by RAM budget, CPU budget, and hard cap.
            Always >= 1.
        """
        effective_mb = (
            self.user_per_worker_mb
            if self.user_per_worker_mb is not None
            else per_worker_mb
        )
        ram_budget_mb = int(self.free_ram_mb() * self.ram_safety_fraction)
        ram_workers = max(1, ram_budget_mb // max(1, effective_mb))
        cpu_workers = max(1, self.physical_cpu_count() - self.cpu_reserve)
        return max(1, min(ram_workers, cpu_workers, self.hard_n_jobs_cap))

    def get_n_jobs(self, per_worker_mb: int = 200, stage: str = "default") -> int:
        """How many parallel workers can we safely spin up?

        Respects the ``user_max_threads`` ceiling set from ``--threads`` flag.
        Default behavior changes from "use all CPUs" to "use what RM permits".

        Args:
            per_worker_mb: Estimated peak RAM per worker, in MB.
            stage: Human-readable label for logging (e.g. ``"nb_pairwise"``).
                Does not affect the computation.

        Returns:
            Worker count >= 1. Bounded by:
              - (free_ram * safety_fraction) / per_worker_mb
              - physical_cpu_count - cpu_reserve
              - hard_n_jobs_cap
              - user_max_threads (absolute ceiling when set)
        """
        budget = self._compute_budget(per_worker_mb)
        if self.user_max_threads is not None:
            budget = min(budget, self.user_max_threads)
        return max(1, budget)

    def split_jobs(
        self,
        n_outer: int,
        per_inner_mb: int = 200,
    ) -> tuple[int, int]:
        """Allocate workers for two nested parallel stages.

        For nested parallelism (e.g. ``ema_switch`` pair × PAS), given
        ``n_outer`` parallel outer-stage units, distributes the total worker
        budget between outer and inner (per-PAS within each pair) so that
        ``n_outer_workers * n_inner_per_outer <= total_budget``.

        Args:
            n_outer: Number of parallel units in the outer stage (e.g. number
                of cluster pairs to test simultaneously).
            per_inner_mb: Estimated RAM per inner worker, in MB.

        Returns:
            Tuple ``(n_outer_workers, n_inner_per_outer)`` where:
                - ``n_outer_workers`` is how many outer units run in parallel.
                - ``n_inner_per_outer`` is how many inner workers each outer
                  unit may spawn.
                - ``n_outer_workers * n_inner_per_outer <= total_budget``.
        """
        total = self.get_n_jobs(per_inner_mb)
        if n_outer >= total:
            return (total, 1)
        outer = min(n_outer, max(1, total // 2))
        inner = max(1, total // outer)
        return (outer, inner)

    def get_batch_size(
        self,
        per_item_mb: float,
        min_batch: int = 100,
        max_batch: int = 100_000,
    ) -> int:
        """Choose a batch size for streaming so a single batch fits in RAM budget.

        Args:
            per_item_mb: Estimated memory footprint of one batch item, in MB.
            min_batch: Minimum batch size regardless of RAM.
            max_batch: Maximum batch size regardless of RAM.

        Returns:
            Batch size between ``min_batch`` and ``max_batch``.
        """
        ram_budget_mb = self.free_ram_mb() * self.ram_safety_fraction
        if per_item_mb <= 0:
            return max_batch
        n = int(ram_budget_mb / per_item_mb)
        return max(min_batch, min(n, max_batch))

    def get_tile_size(
        self,
        bam_path: str,
        n_workers: int,
        target_per_worker_mb: int = 300,
        sample_region_mb: int = 1,
        min_tile_bp: int = 5_000_000,
        max_tile_bp: int = 100_000_000,
    ) -> int:
        """Compute a RAM-adaptive tile size in base pairs for tiled peak calling.

        Samples the first ``sample_region_mb`` Mb of the first chromosome in the
        BAM to estimate read density, then derives a tile size that keeps each
        worker's peak memory within ``target_per_worker_mb`` MB.

        Memory model per worker (conservative):
        - Base overhead: 60 MB (Python interpreter + pysam init)
        - pysam BGZF buffer: 32 MB
        - cb_dict read storage: ~50 bytes/read

        Formula::

            budget_mb = min(target_per_worker_mb, free_ram_mb * 0.7 / n_workers)
            bytes_available = (budget_mb - 92) * 1024 * 1024
            tile_size_bp = bytes_available / 50 / density_reads_per_bp

        Result is clamped to [min_tile_bp, max_tile_bp].

        Args:
            bam_path: Path to the indexed BAM file.
            n_workers: Number of parallel tile workers.
            target_per_worker_mb: Desired peak RAM per worker in MB (default 300).
            sample_region_mb: How many Mb of the first chromosome to sample for
                density estimation (default 1).
            min_tile_bp: Minimum returned tile size in bp (default 5 Mb).
            max_tile_bp: Maximum returned tile size in bp (default 100 Mb).

        Returns:
            Tile size in bp, clamped to [min_tile_bp, max_tile_bp].
        """
        try:
            import pysam

            with pysam.AlignmentFile(bam_path, "rb") as bam:
                if bam.nreferences == 0:
                    logger.warning("get_tile_size: BAM has no references, using min_tile_bp")
                    return min_tile_bp

                first_chrom = bam.get_reference_name(0)
                first_chrom_len = bam.get_reference_length(first_chrom)
                sample_end = min(first_chrom_len, sample_region_mb * 1_000_000)

                read_count = sum(
                    1 for _ in bam.fetch(first_chrom, 0, sample_end)
                )

        except Exception as exc:
            logger.warning(
                "get_tile_size: BAM density sampling failed (%s); using min_tile_bp",
                exc,
            )
            return min_tile_bp

        sample_bp = max(1, sample_end)
        density = read_count / sample_bp  # reads per bp

        free_mb = self.free_ram_mb()
        budget_per_worker_mb = min(
            target_per_worker_mb,
            int(free_mb * self.ram_safety_fraction / max(1, n_workers)),
        )
        # Clamp budget to at least 1 MB above overhead so we don't get 0
        overhead_mb = 92
        usable_mb = max(1, budget_per_worker_mb - overhead_mb)
        bytes_available = usable_mb * 1024 * 1024

        if density <= 0:
            tile_size_bp = max_tile_bp
        else:
            tile_size_bp = int(bytes_available / 50 / density)

        result = max(min_tile_bp, min(tile_size_bp, max_tile_bp))
        logger.info(
            "get_tile_size: bam=%s density=%.6f reads/bp budget_per_worker=%dMB "
            "computed_tile=%d clamped_tile=%d",
            bam_path,
            density,
            budget_per_worker_mb,
            tile_size_bp,
            result,
        )
        return result

    def report(self) -> dict:
        """Return a summary dict of current resource state for logging."""
        return {
            "free_ram_mb": self.free_ram_mb(),
            "physical_cpus": self.physical_cpu_count(),
            "logical_cpus": self.logical_cpu_count(),
            "ram_safety_fraction": self.ram_safety_fraction,
            "cpu_reserve": self.cpu_reserve,
            "hard_n_jobs_cap": self.hard_n_jobs_cap,
            "user_max_threads": self.user_max_threads,
        }
