"""Read-level poly(A) soft-clip evidence (Stage 1 of the caller fix plan).

This module is the single home for every poly(A)-clip primitive the caller
uses:

* :func:`clip_site` — strand-aware terminal soft-clip A/T detector.  Given a
  ``pysam.AlignedSegment`` it returns the inferred cleavage coordinate (the
  0-based genomic position of the last aligned base on the read's 3' side)
  or ``None``.  Semantics match the measured prototype exactly
  (1.152% of CB reads on the PBMC BAM, 92x wrong-end specificity, 73.9% of
  sites within 100 bp of PolyASite 2.0).
* :func:`read_umi` — best-effort UMI extraction (10x ``UB`` tag).
* :class:`ClipAccumulator` — per-chromosome accumulation of clip evidence
  (reads, distinct molecules, per-CB counts) keyed by cleavage site.
* :func:`cluster_clip_sites` — 25 bp single-linkage clustering with a
  read-weighted modal position, ported from the scoring prototype
  (F1@100 0.329 vs 0.177 shipped on chr19+21).
* :class:`ClipSeeder` — the Phase-2 two-tier emitter: clip-site clusters are
  PRIMARY PAS candidates; coverage peaks that do not overlap any cluster are
  emitted as a second tier tagged ``polya_support=0`` (BED score column 0).
* :func:`check_clip_rate` — R2 chemistry mitigation: loud warning when the
  observed clip rate is below 0.3% of CB reads (poly(A) trimmed upstream,
  wrong chemistry, etc.), instead of silently emitting an unsupported set.

Coordinate conventions (identical to the prototype's clip.awk / measure2.py):

* forward read: terminal soft clip at the END of the alignment; clipped
  bases are non-templated A's; cleavage site = ``reference_end - 1``
  (0-based last aligned base).
* reverse read: terminal soft clip at the START of the alignment; clipped
  bases are T's (revcomp of the tail); cleavage site = ``reference_start``.

A site ``s`` corresponds to the 1-bp BED interval ``[s, s+1)``.
"""

from __future__ import annotations

import logging
from bisect import bisect_left, bisect_right

log = logging.getLogger(__name__)

_BAM_CSOFT_CLIP = 4

# Observed-clip-rate warnings are emitted once per BAM per process.
_clip_rate_warned: dict[str, float] = {}

# Below this fraction of CB reads carrying a poly(A) clip, the evidence
# channel is considered destroyed (e.g. poly(A) trimmed before alignment).
LOW_CLIP_RATE = 0.003


def clip_site(read, min_clip: int = 6, min_purity: float = 0.8):
    """Return the inferred cleavage coordinate for a poly(A)-clipped read.

    Strand-aware terminal soft-clip A/T detector.  A read qualifies when its
    3'-side terminal CIGAR operation is a soft clip of at least *min_clip*
    bases, the clipped bases are at least *min_purity* A (forward) / T
    (reverse), and the homopolymer run adjacent to the alignment boundary is
    itself at least *min_clip* long (this last check is what buys the 92x
    wrong-end specificity — a templated A-rich clip rarely starts flush at
    the alignment edge).

    Args:
        read: ``pysam.AlignedSegment``.  Unmapped / CIGAR-less / SEQ-less
            reads return ``None`` (the caller normally invokes this only
            after ``read_check`` succeeds, but the guards keep the function
            safe standalone).
        min_clip: Minimum soft-clip length AND minimum adjacent A/T run.
        min_purity: Minimum A (or T) fraction over the whole clipped segment.

    Returns:
        int cleavage coordinate (0-based, last aligned base on the 3' side:
        ``reference_end - 1`` for forward reads, ``reference_start`` for
        reverse reads) or ``None`` when the read carries no qualifying clip.
    """
    cigar = read.cigartuples
    if not cigar:  # unmapped or CIGAR-less
        return None

    if read.is_reverse:
        op, n = cigar[0]
        if op != _BAM_CSOFT_CLIP or n < min_clip:
            return None
        seq = read.query_sequence
        if not seq:
            return None
        clip = seq[:n]
        if clip.count("T") / n < min_purity:
            return None
        # T-run adjacent to the alignment boundary (end of the clip)
        pure = len(clip) - len(clip.rstrip("T"))
        if pure < min_clip:
            return None
        return read.reference_start

    op, n = cigar[-1]
    if op != _BAM_CSOFT_CLIP or n < min_clip:
        return None
    seq = read.query_sequence
    if not seq:
        return None
    clip = seq[-n:]
    if clip.count("A") / n < min_purity:
        return None
    # A-run adjacent to the alignment boundary (start of the clip)
    pure = len(clip) - len(clip.lstrip("A"))
    if pure < min_clip:
        return None
    ref_end = read.reference_end
    if ref_end is None:
        return None
    return ref_end - 1


def read_umi(read):
    """Return the read's UMI (10x ``UB`` tag) or ``None``.

    Only called for the ~1% of reads that carry a qualifying clip, so the
    extra ``get_tag`` never touches the hot path.
    """
    try:
        return read.get_tag("UB")
    except KeyError:
        return None


class ClipAccumulator:
    """Accumulates clip evidence for one chromosome/strand pass.

    ``sites`` maps ``site -> [n_reads, umi_set, n_reads_without_umi,
    {cb: n_reads}]``.  Distinct-molecule support for a site is
    ``len(umi_set) + n_reads_without_umi`` (a read without a UB tag counts
    as its own molecule — the conservative choice for thresholds).
    """

    __slots__ = ("sites",)

    def __init__(self) -> None:
        self.sites: dict[int, list] = {}

    def add(self, site: int, cb, umi=None) -> None:
        rec = self.sites.get(site)
        if rec is None:
            rec = [0, set(), 0, {}]
            self.sites[site] = rec
        rec[0] += 1
        if umi is None:
            rec[2] += 1
        else:
            rec[1].add((cb, umi))
        cbs = rec[3]
        cbs[cb] = cbs.get(cb, 0) + 1

    def __bool__(self) -> bool:
        return bool(self.sites)


def cluster_clip_sites(sites: dict, seed_window: int = 25) -> list[dict]:
    """Single-linkage clustering of clip sites (prototype logic, productionised).

    Adjacent sites at most *seed_window* bp apart join the same cluster.
    Each cluster's position is the read-weighted mode: the member site with
    the highest read count, ties broken toward the lowest coordinate
    (identical to the prototype's ``pos[argmax(w)]`` on ascending
    positions — deterministic).

    Args:
        sites: ``{site: [n_reads, umi_set, n_no_umi, {cb: n}]}`` as built by
            :class:`ClipAccumulator`.
        seed_window: Single-linkage gap in bp (``--polya-seed-window``).

    Returns:
        List of dicts sorted by mode position, each with keys ``mode``,
        ``nreads`` (total clip reads), ``numis`` (distinct molecules),
        ``cb_counts`` (per-CB clip-read counts across the cluster) and
        ``sites`` (the member coordinates, ascending).
    """
    if not sites:
        return []
    positions = sorted(sites)
    groups: list[list[int]] = []
    current = [positions[0]]
    for p in positions[1:]:
        if p - current[-1] > seed_window:
            groups.append(current)
            current = [p]
        else:
            current.append(p)
    groups.append(current)

    out: list[dict] = []
    for grp in groups:
        # read-weighted mode; ties -> lowest coordinate (prototype argmax)
        mode = max(grp, key=lambda p: (sites[p][0], -p))
        nreads = 0
        n_no_umi = 0
        umis: set = set()
        cb_counts: dict = {}
        for p in grp:
            rec = sites[p]
            nreads += rec[0]
            n_no_umi += rec[2]
            umis |= rec[1]
            for cb, c in rec[3].items():
                cb_counts[cb] = cb_counts.get(cb, 0) + c
        out.append({
            "mode": mode,
            "nreads": nreads,
            "numis": len(umis) + n_no_umi,
            "cb_counts": cb_counts,
            "sites": grp,
        })
    out.sort(key=lambda d: d["mode"])
    return out


class _SupportIndex:
    """O(log n) clip-read counts within a window of a position."""

    __slots__ = ("pos", "cum")

    def __init__(self, sites: dict) -> None:
        self.pos = sorted(sites)
        cum = [0]
        for p in self.pos:
            cum.append(cum[-1] + sites[p][0])
        self.cum = cum

    def reads_within(self, center: int, window: int) -> int:
        lo = bisect_left(self.pos, center - window)
        hi = bisect_right(self.pos, center + window)
        return self.cum[hi] - self.cum[lo]


class ClipSeeder:
    """Two-tier PAS emission for the ``clip_seeded`` strategy (Phase 2).

    Collects, per chromosome:

    * every qualifying clip read (via :meth:`add_clip`), regardless of
      whether a coverage peak fired — this is the point of seeding: ~60% of
      the clip evidence lies outside every coverage-peak window;
    * every coverage-strategy PAS candidate (via :meth:`add_coverage_pas`).

    :meth:`flush` then produces the merged, coordinate-sorted record list:

    * TIER 1 — clip clusters clearing ``min_reads`` distinct molecules,
      emitted as the 1-bp interval ``[mode, mode+1)`` with the cluster's
      raw clip-read count in the BED score column and the cluster's own
      per-CB clip-read counts as the count-matrix column.
    * TIER 2 — coverage PAS that overlap no kept cluster, emitted with their
      clip-read support near the strand-aware 3' base as the score.  A peak
      counts as overlapping when a kept cluster's mode falls in
      ``[start - window, end + window]`` **or** any of that cluster's member
      sites lies within ``+/-window`` of the peak's 3' base — the second
      clause is what makes the score an exact tier tag: at the default
      ``min_reads=1`` no cluster is ever rejected, so every surviving tier-2
      peak scores exactly 0.  ``score == 0`` therefore means coverage_only,
      and pasbed.bed stays plain BED6.  (Raise ``min_reads`` and a tier-2
      peak may score >0 — read that as "has clip evidence, but not enough
      for a call of its own", which is exactly what ``--polya-mode filter``
      should keep.)

    The caller owns pasnumber assignment and the actual writes, so the same
    class serves the monolithic loop, the 3-stage pipeline writer, and the
    tile workers.
    """

    __slots__ = ("direction", "seed_window", "min_reads", "window",
                 "accum", "coverage")

    def __init__(self, direction: bool, seed_window: int = 25,
                 min_reads: int = 1, window: int = 100) -> None:
        self.direction = direction
        self.seed_window = seed_window
        self.min_reads = min_reads
        self.window = window
        self.accum = ClipAccumulator()
        self.coverage: list[tuple[int, int, dict]] = []

    def add_clip(self, site: int, cb, umi=None) -> None:
        self.accum.add(site, cb, umi)

    def add_coverage_pas(self, pas_1: int, pas_2: int, cb_dict: dict) -> None:
        self.coverage.append((pas_1, pas_2, cb_dict))

    def load_sites(self, sites: dict) -> None:
        """Adopt a pre-built site map (pipeline writer receives the finder's)."""
        self.accum.sites = sites

    def flush(self) -> list[tuple[int, int, int, dict]]:
        """Emit and reset. Returns ``[(bed_start, bed_end, score, cb_dict)]``
        sorted by coordinate."""
        sites = self.accum.sites
        clusters = cluster_clip_sites(sites, self.seed_window)
        kept = [c for c in clusters if c["numis"] >= self.min_reads]

        records: list[tuple[int, int, int, dict]] = [
            (c["mode"], c["mode"] + 1, c["nreads"], c["cb_counts"])
            for c in kept
        ]

        if self.coverage:
            modes = [c["mode"] for c in kept]  # already sorted
            support = _SupportIndex(sites)
            # Sites belonging to KEPT clusters: a coverage peak drawing its
            # support from one of these is tier-1 territory, so suppressing
            # on them keeps `score == 0` an exact coverage_only tag.
            kept_sites = _SupportIndex(
                {s: sites[s] for c in kept for s in c["sites"]}
            )
            w = self.window
            for pas_1, pas_2, cb_dict in self.coverage:
                bed_start = min(pas_1, pas_2)
                bed_end = max(pas_1, pas_2)
                three = bed_start if self.direction else bed_end - 1
                # suppressed when a kept cluster's mode falls within the
                # window-expanded interval...
                lo = bisect_left(modes, bed_start - w)
                if lo < len(modes) and modes[lo] <= bed_end + w:
                    continue
                # ...or when a kept cluster's own sites are what would be
                # scored as this peak's support
                if kept_sites.reads_within(three, w):
                    continue
                records.append((bed_start, bed_end,
                                support.reads_within(three, w), cb_dict))

        records.sort(key=lambda r: (r[0], r[1]))
        self.accum = ClipAccumulator()
        self.coverage = []
        return records


def check_clip_rate(bam_path: str, min_clip: int = 6, min_purity: float = 0.8,
                    barcode_tag: str = "CB", max_reads: int = 200_000):
    """Estimate the poly(A) clip rate on the first *max_reads* CB reads and
    warn loudly when it is below :data:`LOW_CLIP_RATE`.

    R2 mitigation (plan section 5): a pipeline that trims poly(A) before
    alignment destroys the evidence channel entirely; the caller must say so
    at startup rather than silently emitting an unsupported call set.

    Cached per BAM per process, so tile/pipeline workers and repeat strand
    passes never rescan.

    Returns the estimated rate (float) or ``None`` when the BAM could not be
    sampled.
    """
    key = str(bam_path)
    if key in _clip_rate_warned:
        return _clip_rate_warned[key]

    try:
        import pysam
        n_cb = 0
        n_clip = 0
        with pysam.AlignmentFile(key, "rb") as bam:
            for read in bam:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if not read.has_tag(barcode_tag):
                    continue
                n_cb += 1
                if clip_site(read, min_clip, min_purity) is not None:
                    n_clip += 1
                if n_cb >= max_reads:
                    break
        rate = (n_clip / n_cb) if n_cb else 0.0
    except Exception as exc:  # never fail a run over the QC estimate
        log.warning("poly(A) clip-rate estimate failed for %s: %s", key, exc)
        _clip_rate_warned[key] = None
        return None

    _clip_rate_warned[key] = rate
    if rate < LOW_CLIP_RATE:
        log.warning(
            "LOW POLY(A) CLIP RATE: %.4f%% of the first %d CB reads in %s "
            "carry a poly(A) soft clip (expected >~0.3%% for 10x cDNA kept "
            "untrimmed; PBMC v3 reference: 1.15%%). The clip evidence channel "
            "looks destroyed (poly(A) trimmed upstream, unusual chemistry, or "
            "aligner hard-clipping). clip_seeded calls and BED score-column "
            "support values from this BAM are NOT reliable.",
            rate * 100.0, n_cb, key,
        )
    else:
        log.info("poly(A) clip rate for %s: %.4f%% of %d sampled CB reads",
                 key, rate * 100.0, n_cb)
    return rate
