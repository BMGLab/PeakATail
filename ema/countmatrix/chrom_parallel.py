"""Per-(contig, strand) parallel peak calling with a deterministic merge.

The default caller streams the whole BAM twice (one pass per strand) in a
single pure-Python thread: on the PBMC 10k v3 BAM (734M records) that is
3h01m of a 3h46m run, whatever ``--threads`` says (``--threads`` only fed
the ``--tiles`` pool and the multi-dataset downstream pool).  Peak calling
is embarrassingly parallel over contigs, though: the streaming state machine
in :func:`~ema.countmatrix.peackcalling.peak_calling` is reset at every
chromosome change and the clip-seeded emitter flushes per chromosome, so a
job that region-fetches one whole contig on one strand sees exactly the
read sequence the sequential pass saw for that contig and produces exactly
its records.

What makes the result **byte-identical** to the sequential two-pass run
(verified on the chr19+21 PBMC slice, ``results/perf/M4_region_parallel``,
and pinned by ``tests/test_chrom_parallel_identity.py``):

* one job per contig with mapped reads (BAI index statistics), both strands,
  each running ``peak_calling(direction, region=(contig, 0, length))`` in a
  spawned worker with a fresh :class:`~ema.countmatrix.indexing.BarcodeIndex`
  and :class:`~ema.countmatrix.peak_state.PeakCallingState`;
* the merge walks the jobs in the order the sequential run would have
  emitted them — ``+`` strand first, contigs in BAM header order, then
  ``-`` strand — renumbering ``pas_id`` across both strands, continuing
  from ``Peak.pasnumber`` (0 for the first BAM of a run, the previous
  BAM's last id afterwards, exactly like the legacy loop);
* the shared cell-barcode index is rebuilt in first-write order: each job's
  local ``cb.tsv`` lists its barcodes in local first-write order, so
  appending the unseen ones while walking the jobs in emission order
  reproduces the sequential singleton's assignment exactly;
* the per-job matrices and support sidecars are re-keyed through those two
  maps; nothing else is touched.

Per-worker memory is the caller's per-chromosome working set (~1.5 GB on
PBMC chr19 (+); the sequential pass holds the same), so 16 workers stay
well under 40 GB.  Wall time is bounded by the largest contig (chr1,
~10 min on PBMC) instead of the whole BAM.

Configuration reaches the workers explicitly through :class:`ChromJob`
(``barcode_tag`` / ``cb_len`` / ``seq_len`` / ``ignore_chro`` / thresholds)
— the spawned interpreter re-parses the parent's ``sys.argv`` through the
legacy ``ema.config`` shim, which does not know the Click flags, so relying
on module state in the child is what crashes ``--tiles``.
"""
from __future__ import annotations

import logging
import multiprocessing
import os
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pysam

from ema.countmatrix.paswrite import SUPPORT_COLUMNS, support_path_for

log = logging.getLogger(__name__)

#: Default per-worker RAM estimate handed to ``ResourceManager.get_n_jobs``.
#: Measured: 1.48 GB for PBMC chr19 (+) (24M accepted reads) as a region
#: job; 2.5 GB leaves room for chr1 of a deeper library.
PER_WORKER_MB = 2500

#: Rows per chunk when re-keying the per-job MatrixMarket triplets.
_MTX_CHUNK_ROWS = 4_000_000


@dataclass
class ChromJob:
    """Everything one spawned worker needs to call one (contig, strand)."""

    job_id: int
    dataset_id: str
    bam_path: str
    contig: str
    length: int
    direction: bool          # False = '+' (forward), True = '-' (reverse)
    mapped_reads: int
    workdir: str
    # --- variable_config in the child (spawn re-imports ema.config with
    #     the defaults, so every value the caller reads is passed along) ---
    barcode_tag: str = "CB"
    cb_len: int = 16
    seq_len: int = 150
    ignore_chro: tuple = ("MT", "mt")
    default_threshold: int = 5
    merge_len: int = 100
    # peakAtail-prime read acceptance geometry.  Workers are SPAWNED, so the
    # legacy globals come back at their module defaults ("fixed", 0) in the
    # child -- these must travel with the job or the parallel path silently
    # runs v2 geometry while the monolithic path runs the branch default.
    read_geometry: str = "fixed"
    read_exclude_flags: int = 0
    # --- strategy (re-instantiated in the child) ---
    strategy_name: str = "original"
    strategy_kwargs: dict = field(default_factory=dict)
    # --- peak_calling kwargs other than the ones above (min_pas_spacing
    #     already resolved by the dispatcher) ---
    peak_kwargs: dict = field(default_factory=dict)
    bam_threads: int = 1
    log_queue: Any = None

    @property
    def strand(self) -> str:
        return "-" if self.direction else "+"

    def out(self, suffix: str) -> str:
        return os.path.join(self.workdir, f"job{self.job_id:04d}.{suffix}")


# ---------------------------------------------------------------------------
# Worker (spawned process)
# ---------------------------------------------------------------------------


def chrom_worker(job: ChromJob) -> dict[str, Any]:
    """Call peaks for one (contig, strand); write bed/mtx/support/cb files."""
    import resource

    if job.log_queue is not None:
        from ema.logging_config import setup_worker_logging
        setup_worker_logging(job.log_queue)

    from ema import config as _cfg
    from ema.countmatrix.indexing import get_mapping, reset_index
    from ema.countmatrix.peackcalling import peak_calling
    from ema.countmatrix.peak import Peak
    from ema.countmatrix.read import set_default_sample_id
    from ema.strategies import get_strategy

    vc = _cfg.variable_config
    vc.barcode_tag = job.barcode_tag
    vc.cb_len = job.cb_len
    vc.seqlen = job.seq_len
    vc.ignore_chro = list(job.ignore_chro)
    vc.default_threshold = job.default_threshold
    vc.merge_len = job.merge_len
    vc.read_geometry = job.read_geometry
    vc.read_exclude_flags = job.read_exclude_flags

    reset_index()
    Peak.reset_pasnumber()
    set_default_sample_id(job.dataset_id)
    strategy = get_strategy(job.strategy_name, **job.strategy_kwargs)

    bed, mtx, cb = job.out("bed"), job.out("mtx"), job.out("cb.tsv")
    t0 = time.monotonic()
    peak_calling(
        job.direction,
        bedfilepath=bed,
        matrixpath=mtx,
        bamfile_dir=job.bam_path,
        default_threshold=job.default_threshold,
        merge_len=job.merge_len,
        strategy=strategy,
        bam_threads=job.bam_threads,
        region=(job.contig, 0, job.length),
        **job.peak_kwargs,
    )
    # Local barcode list in local column order (== local first-write order).
    mapping = get_mapping()
    with open(cb, "w") as fh:
        for cb_str, _ in sorted(mapping.items(), key=lambda kv: kv[1]):
            fh.write(cb_str + "\n")
    return {
        "job_id": job.job_id,
        "contig": job.contig,
        "direction": job.direction,
        "bed": bed,
        "mtx": mtx,
        "support": support_path_for(bed),
        "cb": cb,
        "n_pas": Peak.pasnumber,
        "n_cb": len(mapping),
        "wall_s": time.monotonic() - t0,
        "maxrss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
    }


# ---------------------------------------------------------------------------
# Dispatcher helpers
# ---------------------------------------------------------------------------


def bam_has_index(bam_path: str) -> bool:
    p = str(bam_path)
    return os.path.exists(p + ".bai") or os.path.exists(p[:-4] + ".bai")


def can_run_parallel(bam_path: str) -> bool:
    """True when *bam_path* has a BAM index with per-contig mapped-read
    statistics (what the job plan and the region fetches need)."""
    if not bam_has_index(bam_path):
        return False
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            bam.get_index_statistics()
        return True
    except (ValueError, OSError) as exc:
        log.debug("BAM index statistics unavailable for %s: %s", bam_path, exc)
        return False


def contigs_with_reads(bam_path: str, ignore_chro: Iterable[str]) -> list[tuple[str, int, int]]:
    """``[(contig, length, mapped_reads)]`` in BAM header order, for contigs
    that carry mapped reads and are not excluded by ``ignore_chro``."""
    skip = set(ignore_chro or ())
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        stats = {s.contig: s.mapped for s in bam.get_index_statistics()}
        out = []
        for name, length in zip(bam.references, bam.lengths):
            if name in skip:
                continue
            if stats.get(name, 0) > 0:
                out.append((name, int(length), int(stats[name])))
    return out


def worker_bam_threads(bam_threads: int, n_workers: int, thread_budget: int | None) -> int:
    """BGZF decompression threads per worker so ``workers x threads`` stays
    inside the user's ``--threads`` budget (1 when the budget is used up by
    the workers themselves)."""
    if thread_budget is None or thread_budget <= 0:
        return max(1, int(bam_threads))
    return max(1, min(int(bam_threads), thread_budget // max(1, n_workers)))


def _remap_mtx(src: str, dst, pas_map: np.ndarray, col_map: np.ndarray) -> int:
    """Append *src* triplets to the open file *dst* with both ids re-keyed.
    Returns the number of rows written.  Same text format as
    :func:`~ema.countmatrix.paswrite.matrix_write` (``"pas col count\\n"``)."""
    import pandas as pd

    if not os.path.exists(src) or os.path.getsize(src) == 0:
        return 0
    n = 0
    for chunk in pd.read_csv(
        src, sep=" ", header=None, names=("p", "c", "n"), dtype=np.int64,
        engine="c", chunksize=_MTX_CHUNK_ROWS,
    ):
        out = pd.DataFrame({
            "p": pas_map[chunk["p"].to_numpy()],
            "c": col_map[chunk["c"].to_numpy()],
            "n": chunk["n"].to_numpy(),
        })
        out.to_csv(dst, sep=" ", header=False, index=False, lineterminator="\n")
        n += len(out)
    return n


def merge_chrom_results(
    results: list[dict[str, Any]],
    contig_order: list[str],
    pos_bed: str, neg_bed: str, pos_mtx: str, neg_mtx: str, cb_tsv: str,
    write_support: bool = True,
    pas_start: int = 0,
) -> tuple[int, list[str]]:
    """Merge per-(contig, strand) outputs into the sequential layout.

    ``pas_start`` is the ``pas_id`` the previous dataset ended on: the legacy
    loop seeds every ``peak_calling()`` from the class-level
    ``Peak.pasnumber`` and never resets it between BAMs, so in a multi-BAM
    run the second dataset's ids continue from the first one's last id.

    Returns ``(last_pas_id, cb_list)`` — the last id written (``pas_start``
    plus the number of PAS merged) and the shared barcode list in column
    order (what the sequential run's ``BarcodeIndex.mapping`` would have
    held).
    """
    by_key = {(r["contig"], r["direction"]): r for r in results}
    cb_global: dict[str, int] = {}
    cb_list: list[str] = []
    pasnum = int(pas_start)

    for direction, bed_out_path, mtx_out_path in (
        (False, pos_bed, pos_mtx), (True, neg_bed, neg_mtx),
    ):
        sup_out_path = support_path_for(bed_out_path)
        with open(bed_out_path, "w") as bed_out, open(mtx_out_path, "w") as mtx_out:
            sup_out = open(sup_out_path, "w") if write_support else None
            if sup_out is not None:
                sup_out.write("\t".join(SUPPORT_COLUMNS) + "\n")
            try:
                for contig in contig_order:
                    r = by_key.get((contig, direction))
                    if r is None:
                        continue
                    # --- barcodes: local column -> global column --------------
                    with open(r["cb"]) as fh:
                        local_cbs = [ln.rstrip("\n") for ln in fh]
                    col_map = np.zeros(len(local_cbs) + 1, dtype=np.int64)
                    for i, cb in enumerate(local_cbs, 1):
                        g = cb_global.get(cb)
                        if g is None:
                            cb_list.append(cb)
                            g = len(cb_list)
                            cb_global[cb] = g
                        col_map[i] = g
                    # --- BED: local pas -> global pas -------------------------
                    local_to_global: dict[int, int] = {}
                    with open(r["bed"]) as fh:
                        for line in fh:
                            if not line.strip():
                                continue
                            f = line.split("\t")
                            pasnum += 1
                            local_to_global[int(f[3])] = pasnum
                            f[3] = str(pasnum)
                            bed_out.write("\t".join(f))
                    pas_map = np.zeros(max(local_to_global, default=0) + 1, dtype=np.int64)
                    for lp, gp in local_to_global.items():
                        pas_map[lp] = gp
                    # --- matrix ---------------------------------------------
                    _remap_mtx(r["mtx"], mtx_out, pas_map, col_map)
                    # --- support sidecar ------------------------------------
                    sp = r.get("support")
                    if sup_out is not None and sp and os.path.exists(sp):
                        with open(sp) as fh:
                            first = fh.readline()
                            if first and not first.startswith("pas_id"):
                                fh.seek(0)
                            for line in fh:
                                if not line.strip():
                                    continue
                                lp, rest = line.split("\t", 1)
                                sup_out.write(f"{local_to_global[int(lp)]}\t{rest}")
            finally:
                if sup_out is not None:
                    sup_out.close()

    with open(cb_tsv, "w") as fh:
        for cb in cb_list:
            fh.write(cb + "\n")
    return pasnum, cb_list


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


def run_chrom_parallel(
    *,
    dataset_id: str,
    bam_path: str,
    pos_bed: str, neg_bed: str, pos_mtx: str, neg_mtx: str, cb_tsv: str,
    n_workers: int,
    strategy_name: str,
    strategy_kwargs: dict | None,
    peak_kwargs: dict,
    workdir: str,
    thread_budget: int | None = None,
    progress_client=None,
    log_queue=None,
) -> dict[str, Any]:
    """Call peaks on both strands of *bam_path* with one job per (contig,
    strand) and merge deterministically into the five sequential-layout
    files.  ``peak_kwargs`` is the dict ``main.py`` passes to
    :func:`~ema.countmatrix.peackcalling.peak_calling` (the ``strategy``
    object is replaced by *strategy_name* / *strategy_kwargs* so it can be
    rebuilt in the child).

    Returns a summary dict (``n_pas``, ``n_cb``, ``cb_list``, per-job
    timings).
    """
    from ema.config import variable_config as vc
    from ema.countmatrix.peak import Peak

    peak_kwargs = dict(peak_kwargs)
    peak_kwargs.pop("strategy", None)
    bam_threads = int(peak_kwargs.pop("bam_threads", 4))
    # Parameters the worker sets from the job itself.
    for k in ("region", "progress_client", "index", "state", "sample_id",
              "use_tiles", "use_pipeline", "default_threshold", "merge_len"):
        peak_kwargs.pop(k, None)

    ignore_chro = tuple(vc.ignore_chro or ())
    contigs = contigs_with_reads(bam_path, ignore_chro)
    contig_order = [c for c, _, _ in contigs]

    # Once-per-BAM work the sequential peak_calling() does on its full-scan
    # path and region jobs skip: the poly(A) clip-rate QC and the median
    # read length for the PAS merger's distance tier.
    polya_on = bool(peak_kwargs.get("polya_enabled", True)) or strategy_name == "clip_seeded"
    if polya_on:
        from ema.countmatrix.polya import check_clip_rate
        check_clip_rate(
            str(bam_path),
            min_clip=peak_kwargs.get("polya_min_clip", 6),
            min_purity=peak_kwargs.get("polya_min_purity", 0.8),
            barcode_tag=vc.barcode_tag or "CB",
        )
    if int(peak_kwargs.get("min_pas_spacing", -1)) < 0:
        from ema.countmatrix.bam_utils import infer_median_read_length
        cache = vc.dataset_read_lengths
        key = str(bam_path)
        if key not in cache:
            cache[key] = infer_median_read_length(key)
            log.info(
                "Auto-detected median read length %d bp for %s "
                "(min_pas_spacing distance tier)", cache[key], bam_path,
            )
        peak_kwargs["min_pas_spacing"] = cache[key]

    os.makedirs(workdir, exist_ok=True)
    jobs: list[ChromJob] = []
    for direction in (False, True):
        for contig, length, mapped in contigs:
            jobs.append(ChromJob(
                job_id=len(jobs), dataset_id=dataset_id, bam_path=str(bam_path),
                contig=contig, length=length, direction=direction,
                mapped_reads=mapped, workdir=workdir,
                barcode_tag=vc.barcode_tag or "CB", cb_len=int(vc.cb_len or 16),
                seq_len=int(vc.seqlen), ignore_chro=ignore_chro,
                default_threshold=int(vc.default_threshold),
                merge_len=int(vc.merge_len),
                read_geometry=str(vc.read_geometry),
                read_exclude_flags=int(vc.read_exclude_flags),
                strategy_name=strategy_name,
                strategy_kwargs=dict(strategy_kwargs or {}),
                peak_kwargs=peak_kwargs, log_queue=log_queue,
            ))
    n_jobs = len(jobs)
    n_workers = max(1, min(int(n_workers), n_jobs)) if n_jobs else 1
    per_worker_threads = worker_bam_threads(bam_threads, n_workers, thread_budget)
    for j in jobs:
        j.bam_threads = per_worker_threads
    # Longest-processing-time first: the biggest contigs bound the wall time.
    order = sorted(jobs, key=lambda j: (-j.mapped_reads, j.job_id))

    log.info(
        "peak calling: %d (contig, strand) jobs on %s through %d worker(s) "
        "(%d BGZF thread(s) each); largest %s%s with %d mapped reads",
        n_jobs, bam_path, n_workers, per_worker_threads,
        order[0].contig if order else "-", order[0].strand if order else "",
        order[0].mapped_reads if order else 0,
    )
    if progress_client is not None:
        try:
            progress_client.set_total(max(1, n_jobs))
        except Exception:
            pass

    results: list[dict[str, Any]] = []
    t0 = time.monotonic()
    try:
        if n_jobs:
            ctx = multiprocessing.get_context("spawn")
            with ctx.Pool(processes=n_workers) as pool:
                try:
                    for res in pool.imap_unordered(chrom_worker, order, chunksize=1):
                        results.append(res)
                        log.info(
                            "  %s (%s): %d PAS, %d cells, %.0f s, %.2f GB RSS "
                            "[%d/%d done, %.0f s elapsed]",
                            res["contig"], "-" if res["direction"] else "+",
                            res["n_pas"], res["n_cb"], res["wall_s"],
                            res["maxrss_mb"] / 1024.0, len(results), n_jobs,
                            time.monotonic() - t0,
                        )
                        if progress_client is not None:
                            try:
                                progress_client.advance(1)
                            except Exception:
                                pass
                except Exception:
                    pool.terminate()
                    raise
        t_merge = time.monotonic()
        # Legacy numbering: peak_calling() seeds its state from the
        # class-level Peak.pasnumber and writes the final value back, and
        # main.py never resets it between BAMs -- so dataset 2 continues
        # where dataset 1 stopped.  Reproduce that here (worker processes
        # start at 0; only the merged ids matter).
        pas_start = int(Peak.pasnumber)
        last_pas, cb_list = merge_chrom_results(
            results, contig_order, pos_bed, neg_bed, pos_mtx, neg_mtx, cb_tsv,
            write_support=polya_on, pas_start=pas_start,
        )
        Peak.pasnumber = last_pas
        n_pas = last_pas - pas_start
        log.info(
            "peak calling merged: %d PAS (ids %d..%d), %d cells; jobs %.0f s, "
            "merge %.0f s",
            n_pas, pas_start + 1, last_pas, len(cb_list), t_merge - t0,
            time.monotonic() - t_merge,
        )
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    return {
        "n_pas": n_pas,
        "n_cb": len(cb_list),
        "cb_list": cb_list,
        "n_workers": n_workers,
        "bam_threads_per_worker": per_worker_threads,
        "jobs": [
            {k: r[k] for k in ("contig", "direction", "n_pas", "n_cb", "wall_s", "maxrss_mb")}
            for r in sorted(results, key=lambda r: r["job_id"])
        ],
    }


def load_index_from_cb_list(cb_list: list[str]) -> None:
    """Rebuild the module-level :class:`BarcodeIndex` singleton from the
    merged barcode list so ``get_mapping()`` (used by the CB filter when no
    explicit list is passed) sees the same column order the sequential run
    would have produced."""
    from ema.countmatrix.indexing import _index, reset_index

    reset_index()
    for cb in cb_list:
        _index.get_index(cb)
