"""Read-level poly(A) soft-clip evidence (Stage 1 of the caller fix plan).

This module is the single home for every poly(A)-clip primitive the caller
uses:

* :func:`clip_site` — strand-aware terminal soft-clip A/T detector.  Given a
  ``pysam.AlignedSegment`` it returns the inferred cleavage coordinate (the
  0-based genomic position of the last aligned base on the read's 3' side)
  or ``None``.  Semantics match the measured prototype exactly
  (92x wrong-end specificity, 73.9% of sites within 100 bp of PolyASite 2.0).

  **The qualifying clip rate is 0.5730 % of accepted CB reads genome-wide on
  the PBMC 10k v3 BAM** (3,195,067 / 557,564,408), 0.5364 % on
  :func:`check_clip_rate`'s own denominator, and 0.3669 % (chr21) to 0.8179 %
  (chr19) per chromosome.  The **1.152 %** this docstring used to state was
  measured on the head of a coordinate-sorted BAM and is wrong by ~2x; it was
  independently re-derived from scratch in gawk over ``samtools view`` and
  agreed to six decimals (``results/algo_headroom/VERIFY``, §5).
* :func:`read_umi` — best-effort UMI extraction (10x ``UB`` tag).
* :func:`clip_read_ok` — ``samtools view -F 3844`` alignment filter
  (unmapped / secondary / supplementary / qcfail / duplicate) applied to
  every read that feeds the clip-evidence channel.
* :class:`ClipAccumulator` — per-chromosome accumulation of clip evidence
  keyed by cleavage site.  Support is counted in DISTINCT (CB, UMI)
  MOLECULES — the unit ``--polya-min-umis`` documents — both over every
  accepted clip read and over the ``-F 3844`` subset alone.
* :func:`cluster_clip_sites` — 25 bp single-linkage clustering with a
  read-weighted modal position, ported from the scoring prototype
  (F1@100 0.329 vs 0.177 shipped on chr19+21).
* :class:`ClipSeeder` — the Phase-2 two-tier emitter: clip-site clusters are
  PRIMARY PAS candidates; coverage peaks that do not overlap any cluster are
  emitted as a second tier tagged ``polya_support=0`` (BED score column 0).
* :func:`check_clip_rate` — R2 chemistry mitigation: loud warning when the
  observed clip rate is below 0.3% of CB reads (poly(A) trimmed upstream,
  wrong chemistry, etc.), instead of silently emitting an unsupported set.
  ``--clip-rate-sampling strided`` (the peakAtail-prime default) spreads the
  sample across the whole BAM; ``head`` is v2's first-200k-CB-reads scan,
  which on a coordinate-sorted BAM samples the head of chr1 and reported
  2.2565 % where the truth on the same library is 0.5364 %.

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
from array import array
from bisect import bisect_left, bisect_right
from collections import Counter

log = logging.getLogger(__name__)

_BAM_CSOFT_CLIP = 4

#: ``samtools view -F 3844``: unmapped (4) + secondary (256) + qcfail (512)
#: + duplicate (1024) + supplementary (2048).  ``read_check`` applies none of
#: these — it is the coverage path's contract and widening it would move
#: every peak — so the filter lives here, on the clip-evidence channel only.
CLIP_EXCLUDE_FLAGS = 3844

#: Accepted values of ``--polya-clip-filter``.
CLIP_FILTERS = ("none", "f3844")

# Observed-clip-rate warnings are emitted once per BAM per process.
_clip_rate_warned: dict[str, float] = {}

# Below this fraction of CB reads carrying a poly(A) clip, the evidence
# channel is considered destroyed (e.g. poly(A) trimmed before alignment).
# Reference point for the bar: the genome-wide truth on PBMC 10k v3 is
# 0.5364 % on this function's own denominator (0.3669 %-0.8179 % per
# chromosome), so 0.3 % is roughly "half of what a healthy 10x cDNA library
# gives".  See the module docstring for the provenance of that number.
LOW_CLIP_RATE = 0.003

#: Accepted values of ``--clip-rate-sampling``.
#:
#: ``"head"``
#:     v2: scan the BAM from record 0 until ``max_reads`` CB reads have been
#:     seen, which on a coordinate-sorted BAM means the head of the first
#:     contig.  On PBMC 10k v3 it returns 2.2565 % where the whole-file rate on
#:     the same denominator is 0.5364 %.
#: ``"strided"``
#:     Sample coordinate-uniform windows across every mapped contig, taking
#:     every read whose start falls in one.  Unbiased in construction and
#:     cheap (3 s), but MEASURED TO BE UNRELIABLE at an affordable budget:
#:     poly(A) clips are rare and concentrated at the 3' ends of expressed
#:     genes, so the estimate is dominated by which windows happen to hit one.
#:     On the full PBMC BAM against a 0.5364 % truth it returns 1.5258 %
#:     (200,000 reads) and 1.8068 % (1,000,000); on the chr19+21 slice against
#:     a 0.6169 % truth, 0.44 %-0.79 %.  Kept because it is the "sample across
#:     the BAM" design and its failure is worth being able to reproduce.
#: ``"pass"``
#:     **The peakAtail-prime default, and the only one that is not an
#:     estimate.**  Skip the startup scan entirely and COUNT, during the peak
#:     calling pass the run does anyway, every read the caller accepts and
#:     every qualifying poly(A) clip among them.  Exact, per (contig, strand),
#:     for zero extra I/O -- ``clip_site`` is already called on each of those
#:     reads.  Its denominator is the caller's own accepted-read denominator
#:     (0.5730 % genome-wide on PBMC 10k v3), which is the evidence the caller
#:     actually has, rather than ``check_clip_rate``'s slightly wider one.
CLIP_RATE_SAMPLINGS = ("head", "strided", "pass")

#: The v2-compatibility value of ``--clip-rate-sampling``.
V2_CLIP_RATE_SAMPLING = "head"

#: Windows per contig used by ``"strided"``.  Every read whose start falls in
#: a window is taken, so dense regions contribute proportionally more reads --
#: see :func:`_sample_strided` for why an equal-READS-per-stratum design is
#: biased and this one is not.
CLIP_RATE_STRATA = 10

#: Minimum window width (bp).  A narrow window is expensive rather than wrong:
#: ``fetch`` yields every read OVERLAPPING it, and a spliced alignment spans
#: its introns, so a 2 kb window in a gene-dense region yields ~685 reads of
#: which ~5 START inside it.  Measured on the full PBMC BAM: at 100 windows
#: per contig (2 kb wide) collecting 1 M in-window reads iterated ~100x that
#: many and took >10 min of CPU.  Widening the window amortises the span
#: overhead; the number of windows is reduced to keep the sampled FRACTION of
#: each contig unchanged, so the estimator is the same one.
CLIP_RATE_MIN_WINDOW = 25_000

#: CB-read budget for ``"strided"``.  Larger than v2's 200,000 because poly(A)
#: clips are RARE (0.57 % genome-wide) and CLUSTERED at 3' ends, so a
#: coordinate-window sample either hits a pileup or misses it: measured on the
#: PBMC chr19+21 slice against a whole-file truth of 0.6169 %, the estimate
#: lands at 0.44 % with a 200,000-read budget and 0.61 % with 1,000,000
#: (9 s).  The estimator is an ALARM, not a measurement -- across the budget
#: and window settings probed it lands within 0.85x-1.28x of truth, against
#: 4.2x for v2's head scan (2.2565 % vs 0.5364 % on the full PBMC BAM).
CLIP_RATE_STRIDED_READS = 1_000_000

#: CB-read budget for ``"head"`` -- v2's number, kept exactly.
CLIP_RATE_HEAD_READS = 200_000


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


def clip_read_ok(read) -> bool:
    """Return ``True`` when *read* passes ``samtools view -F 3844``.

    Mirrors the read definition the recall-ceiling measurement used
    (``samtools view -F 3844``): unmapped, secondary, supplementary,
    QC-fail and PCR/optical-duplicate alignments are not poly(A) evidence.

    Note what this is NOT for: PCR duplicates cannot inflate a *molecule*
    count in the first place, because duplicate reads share their
    ``(CB, UMI)`` key and collapse into one molecule.  Excluding them here
    additionally discards molecules whose unflagged representative carries
    no qualifying clip — which is why the filter is selectable
    (``--polya-clip-filter``) rather than unconditional.  See
    :meth:`ClipAccumulator.add`.
    """
    flag = getattr(read, "flag", None)
    if flag is not None:
        return not (flag & CLIP_EXCLUDE_FLAGS)
    return not (                                   # pragma: no cover - mocks
        read.is_unmapped or read.is_secondary or read.is_supplementary
        or read.is_qcfail or read.is_duplicate
    )


def molecule_cb(cb: str) -> str:
    """Barcode half of the composite ``"{RG}_{barcode}"`` cell id.

    ``read_check`` keys cells by ``RG + "_" + barcode`` (and
    :func:`~ema.countmatrix.indexing.split_cb` splits on the LAST
    underscore, the barcode being a fixed-length ACGTN string).  A MOLECULE,
    though, is a ``(cell barcode, UMI)`` pair within a library: a 10x BAM
    with one read group per flow-cell lane would otherwise report the same
    molecule once per lane — measured at +4.7% on the PBMC chr21 (+) slice,
    2 read groups (13,100 vs 12,483 molecules over identical clusters).

    Cells stay keyed by the composite id everywhere else (matrix columns are
    unchanged); only the molecule identity drops the read-group prefix.
    """
    i = cb.rfind("_")
    return cb[i + 1:] if i >= 0 else cb


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

    ``sites`` maps ``site -> record``, a plain list so the hot path stays
    free of attribute lookups::

        [0] n_reads          every accepted clip read
        [1] umis             set of (barcode, umi) molecules over those
                             reads (see :func:`molecule_cb`)
        [2] n_no_umi         reads without a UB tag (one molecule each)
        [3] cb_counts        {cb: n_reads}
        [4] ends             {end1: {cb: n}} — the peak caller's read-end
                             coordinate of each clip read, needed to keep
                             the tier-1 clip fallback from counting a read
                             a neighbouring cluster already counted
        [5] n_reads_f3844    clip reads passing :func:`clip_read_ok`
        [6] umis_f3844       set of (cb, umi) over those reads
        [7] n_no_umi_f3844   of those, reads without a UB tag

    Distinct-molecule support for a site is ``len(umis) + n_no_umi`` (a read
    without a UB tag counts as its own molecule — the conservative choice
    for thresholds), and the same arithmetic on slots 6/7 gives the
    ``-F 3844`` molecule count.  Both are always tracked; the emitted
    support unit is chosen by ``--polya-clip-filter``.
    """

    __slots__ = ("sites",)

    def __init__(self) -> None:
        self.sites: dict[int, list] = {}

    def add(self, site: int, cb, umi=None, end1=None, primary: bool = True) -> None:
        """Record one clip read.

        Args:
            site: Cleavage coordinate from :func:`clip_site`.
            cb: Composite cell barcode.
            umi: UB tag or ``None`` (``None`` == its own molecule).
            end1: The read's ``end1`` as the peak caller tracks it
                (``start1 + seq_len``).  Optional — without it the clip
                fallback cannot exclude reads counted by a neighbour.
            primary: Result of :func:`clip_read_ok` for this read
                (``-F 3844``).  Defaults to ``True`` so hand-built test
                fixtures and legacy callers keep their old meaning.
        """
        rec = self.sites.get(site)
        if rec is None:
            rec = [0, set(), 0, {}, {}, 0, set(), 0]
            self.sites[site] = rec
        rec[0] += 1
        mol = None if umi is None else (molecule_cb(cb), umi)
        if mol is None:
            rec[2] += 1
        else:
            rec[1].add(mol)
        cbs = rec[3]
        cbs[cb] = cbs.get(cb, 0) + 1
        if end1 is not None:
            at_end = rec[4].get(end1)
            if at_end is None:
                rec[4][end1] = {cb: 1}
            else:
                at_end[cb] = at_end.get(cb, 0) + 1
        if primary:
            rec[5] += 1
            if mol is None:
                rec[7] += 1
            else:
                rec[6].add(mol)

    def __bool__(self) -> bool:
        return bool(self.sites)


def site_molecules(rec: list, strict: bool = False) -> int:
    """Distinct molecules recorded at one site (``-F 3844`` subset when
    *strict*)."""
    if strict:
        return len(rec[6]) + rec[7]
    return len(rec[1]) + rec[2]


def cluster_clip_sites(sites: dict, seed_window: int = 25) -> list[dict]:
    """Single-linkage clustering of clip sites (prototype logic, productionised).

    Adjacent sites at most *seed_window* bp apart join the same cluster.
    Each cluster's position is the read-weighted mode: the member site with
    the highest read count, ties broken toward the lowest coordinate
    (identical to the prototype's ``pos[argmax(w)]`` on ascending
    positions — deterministic).

    Cluster geometry (membership and the modal position) is a pure function
    of the ACCEPTED clip reads, i.e. it does not depend on the alignment
    filter or on the support unit: changing ``--polya-clip-filter`` or
    ``--polya-min-umis`` can drop a cluster, but it never MOVES one.

    Args:
        sites: ``{site: record}`` as built by :class:`ClipAccumulator`.
        seed_window: Single-linkage gap in bp (``--polya-seed-window``).

    Returns:
        List of dicts sorted by mode position, each with keys ``mode``,
        ``nreads`` (total clip reads), ``numis`` (distinct molecules),
        ``nreads_f3844`` / ``numis_f3844`` (the same two counts over the
        ``-F 3844`` subset), ``cb_counts`` (per-CB clip-read counts across
        the cluster) and ``sites`` (the member coordinates, ascending).
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
        nreads_f = 0
        n_no_umi_f = 0
        umis_f: set = set()
        for p in grp:
            rec = sites[p]
            nreads += rec[0]
            n_no_umi += rec[2]
            umis |= rec[1]
            for cb, c in rec[3].items():
                cb_counts[cb] = cb_counts.get(cb, 0) + c
            nreads_f += rec[5]
            n_no_umi_f += rec[7]
            umis_f |= rec[6]
        out.append({
            "mode": mode,
            "nreads": nreads,
            "numis": len(umis) + n_no_umi,
            "nreads_f3844": nreads_f,
            "numis_f3844": len(umis_f) + n_no_umi_f,
            "cb_counts": cb_counts,
            "sites": grp,
        })
    out.sort(key=lambda d: d["mode"])
    return out


class _SupportIndex:
    """O(log n) clip support within a window of a position.

    Reads are a prefix-sum lookup; molecules need the ``(cb, umi)`` sets
    unioned over the window (summing per-site molecule counts would count a
    molecule twice when it clipped at two sites inside the same window), so
    :meth:`molecules_within` walks the — always short — window slice.
    """

    __slots__ = ("pos", "cum", "cum_f", "recs")

    def __init__(self, sites: dict) -> None:
        self.pos = sorted(sites)
        self.recs = [sites[p] for p in self.pos]
        cum = [0]
        cum_f = [0]
        for rec in self.recs:
            cum.append(cum[-1] + rec[0])
            cum_f.append(cum_f[-1] + rec[5])
        self.cum = cum
        self.cum_f = cum_f

    def _slice(self, center: int, window: int) -> tuple[int, int]:
        return (bisect_left(self.pos, center - window),
                bisect_right(self.pos, center + window))

    def reads_within(self, center: int, window: int, strict: bool = False) -> int:
        lo, hi = self._slice(center, window)
        cum = self.cum_f if strict else self.cum
        return cum[hi] - cum[lo]

    def sites_within(self, center: int, window: int) -> tuple[int, int]:
        """``(n distinct clip positions, bp span)`` in the window.

        peakAtail-prime ``--pas-features on``: the tier-2 counterpart of a
        tier-1 cluster's ``len(sites)`` / member span.  Two bisects and two
        list reads -- the same slice ``reads_within`` already computes.
        """
        lo, hi = self._slice(center, window)
        if hi <= lo:
            return 0, 0
        return hi - lo, self.pos[hi - 1] - self.pos[lo]

    def molecules_within(self, center: int, window: int,
                         strict: bool = False) -> int:
        lo, hi = self._slice(center, window)
        if hi <= lo:
            return 0
        if hi - lo == 1:
            return site_molecules(self.recs[lo], strict)
        umis: set = set()
        n_no_umi = 0
        i_set, i_no = (6, 7) if strict else (1, 2)
        for k in range(lo, hi):
            rec = self.recs[k]
            umis |= rec[i_set]
            n_no_umi += rec[i_no]
        return len(umis) + n_no_umi


def parse_count_window(spec) -> tuple[int, int]:
    """Parse ``--polya-count-window`` into ``(upstream, downstream)`` bp.

    The tier-1 count window is expressed in *transcript* orientation around
    the cluster's cleavage site: ``upstream`` bp toward the transcript 5' end
    (where R2 3' ends pile up) and ``downstream`` bp past the site.

    Accepted forms: ``"auto,25"`` (default; ``auto`` == the run's
    ``--seq-len``), ``"98,25"``, a lone ``"25"`` (downstream only, upstream
    auto), or an already-parsed ``(up, down)`` pair where ``up == -1`` means
    auto.
    """
    if isinstance(spec, (tuple, list)):
        if len(spec) != 2:
            raise ValueError(f"polya_count_window needs (up, down), got {spec!r}")
        up, down = int(spec[0]), int(spec[1])
    else:
        parts = [p.strip().lower() for p in str(spec).split(",")]
        if len(parts) == 1:
            parts = ["auto", parts[0]]
        if len(parts) != 2:
            raise ValueError(
                f"--polya-count-window expects 'UP,DOWN' (e.g. 'auto,25'), got {spec!r}"
            )
        up = -1 if parts[0] in ("auto", "", "-1") else int(parts[0])
        down = int(parts[1])
    if up < -1 or down < 0:
        raise ValueError(
            f"--polya-count-window: upstream must be >= 0 (or auto), downstream >= 0, got {spec!r}"
        )
    return up, down


def _merge_intervals(intervals: list[tuple[int, int]]) -> tuple[list[int], list[int]]:
    """Sort + merge inclusive integer intervals -> parallel (starts, ends)."""
    starts: list[int] = []
    ends: list[int] = []
    for s, e in sorted(intervals):
        if ends and s <= ends[-1] + 1:
            if e > ends[-1]:
                ends[-1] = e
        else:
            starts.append(s)
            ends.append(e)
    return starts, ends


def _subtract_intervals(lo: int, hi: int, starts: list[int], ends: list[int]):
    """Yield the sub-ranges of inclusive ``[lo, hi]`` not covered by the merged
    inclusive intervals ``(starts, ends)``."""
    if hi < lo:
        return
    cur = lo
    k = bisect_left(ends, lo)
    while k < len(starts) and starts[k] <= hi:
        if starts[k] > cur:
            yield cur, starts[k] - 1
        if ends[k] + 1 > cur:
            cur = ends[k] + 1
        k += 1
    if cur <= hi:
        yield cur, hi


def _point_covered(x: int, starts: list[int], ends: list[int]) -> bool:
    """True when ``x`` lies inside one of the merged inclusive intervals."""
    if not starts:
        return False
    k = bisect_right(starts, x) - 1
    return k >= 0 and x <= ends[k]


def _merge_counts(into: dict, src: dict) -> None:
    for cb, n in src.items():
        into[cb] = into.get(cb, 0) + n


class _CandidateSlice:
    """Compact positional view of one coverage candidate.

    The ``(end1, cb, n)`` triples of ``peak.cb_positions`` restricted to the
    candidate's inclusive ``[bed_start, bed_end]``, sorted by position, as
    three parallel ``int32`` arrays (cell barcodes interned by the owning
    :class:`ClipSeeder`).  12 bytes per (position, cell) pair instead of the
    nested dicts — keeping the ``Peak`` objects themselves alive until the
    chromosome flush cost ~4.5 GB on a single gene-dense human chromosome.

    :meth:`counts` over an inclusive range is exactly
    :func:`~ema.strategies.utils.reconstruct_cb_dict` on the same positions,
    which is what every positional coverage strategy's
    ``get_cb_dict_for_pas`` delegates to.
    """

    __slots__ = ("pos", "cbid", "n")

    def __init__(self, cb_positions: dict, bed_start: int, bed_end: int,
                 intern) -> None:
        pos = array("i")
        cbid = array("i")
        n = array("i")
        for p in sorted(q for q in cb_positions if bed_start <= q <= bed_end):
            for cb, c in cb_positions[p].items():
                pos.append(p)
                cbid.append(intern(cb))
                n.append(c)
        self.pos, self.cbid, self.n = pos, cbid, n

    def __bool__(self) -> bool:
        return bool(self.pos)

    def counts(self, lo: int, hi: int, names: list) -> dict:
        out: dict = {}
        if hi < lo:
            return out
        i = bisect_left(self.pos, lo)
        j = bisect_right(self.pos, hi)
        cbid, n = self.cbid, self.n
        for k in range(i, j):
            cb = names[cbid[k]]
            out[cb] = out.get(cb, 0) + n[k]
        return out


class ClipStream:
    """Streaming-side poly(A) evidence for one chromosome/strand pass.

    Two things are recorded as the reads stream by:

    * every qualifying clip read, keyed by cleavage site
      (:class:`ClipAccumulator`, via :meth:`add_clip`);
    * every accepted read end, as two parallel ``int32`` arrays —
      ``ends[i]`` is the read's ``end1`` as tracked by the peak caller and
      ``cbids[i]`` an interned cell-barcode id (via :meth:`add_read`).

    Under ``--read-geometry fixed`` (v2) ``read_check`` pads/drops reads so
    that ``end1 - start1 == seq_len`` for every accepted read; the BAM is
    coordinate-sorted, so ``ends`` is non-decreasing in arrival order and any
    ``[lo, hi]`` window is an ``O(log n)`` bisect plus a ``Counter`` over the
    slice.  Memory is 8 bytes per accepted read (plus one interned string per
    distinct barcode), i.e. tens of MB for the deepest chromosome of a
    200M-read BAM — and it is released at every chromosome flush.

    **Under the peakAtail-prime geometries that invariant is gone**, so two
    things change and both are explicit rather than inferred:

    * ``seq_len`` is no longer learned from the first read (``end1 - start1``
      is not constant any more) — it MUST be supplied by the caller from
      ``--seq-len``.  Constructing a non-``fixed`` stream without one raises,
      rather than silently shifting every minus-strand count window.
    * the key stored per read is the read's TRANSCRIPT 3'-MOST coordinate
      expressed in the caller's ``end1`` space (:meth:`three_key`), not the
      raw ``end1``.  On ``+`` that is ``end1`` itself (the true 3'-most
      aligned base, which under ``true`` geometry is no longer fabricated).
      On ``-`` the transcript 3' end is ``start1`` — exact in every geometry —
      and it is stored as ``start1 + seq_len`` so that it keeps sitting
      exactly on the cluster anchor ``mode + seq_len`` the rest of
      :class:`ClipSeeder` works in.  Under ``fixed`` both reduce to ``end1``,
      so v2 bytes are untouched.

    Sortedness is therefore preserved by construction everywhere EXCEPT
    ``true`` geometry on the ``+`` strand, where a soft-clipped read's true
    end can fall behind its predecessor's by up to ``seq_len``.  That case is
    detected while streaming and repaired once, at :meth:`finalize`, with a
    stable sort — ``count_ends`` bisects and must never see an unsorted array.

    The object is picklable (the 3-stage pipeline hands one per chromosome
    from the finder to the writer); the intern dict is dropped on pickling
    because the receiving side only reads.
    """

    __slots__ = ("accum", "ends", "cbids", "cb_names", "seq_len", "_cb_ids",
                 "geometry", "direction", "_sorted")

    def __init__(self, seq_len: int | None = None, geometry: str = "fixed",
                 direction: bool = False) -> None:
        if geometry != "fixed" and not seq_len:
            raise ValueError(
                "ClipStream(geometry=%r) needs an explicit seq_len: outside "
                "v2 geometry `end1 - start1` is not constant, so it cannot be "
                "learned from the first read." % (geometry,)
            )
        self.accum = ClipAccumulator()
        self.ends = array("i")
        self.cbids = array("i")
        self.cb_names: list = []
        self.seq_len: int | None = seq_len
        self.geometry = geometry
        self.direction = bool(direction)
        self._sorted = True
        self._cb_ids: dict = {}

    # -- coordinate space --------------------------------------------------
    def three_key(self, start1: int, end1: int) -> int:
        """The read's transcript 3'-most coordinate in ``end1`` space.

        ``fixed`` -> ``end1`` (v2, and ``end1 == start1 + seq_len`` anyway).
        ``+`` strand -> ``end1``.  ``-`` strand -> ``start1 + seq_len``.
        """
        if self.geometry == "fixed" or not self.direction:
            return end1
        return start1 + self.seq_len

    # -- streaming hooks ---------------------------------------------------
    def add_clip(self, site: int, cb, umi=None, end1=None,
                 primary: bool = True) -> None:
        self.accum.add(site, cb, umi, end1, primary)

    def add_read(self, start1: int, end1: int, cb) -> None:
        if self.seq_len is None:
            self.seq_len = end1 - start1
        key = end1 if self.geometry == "fixed" else self.three_key(start1, end1)
        cid = self._cb_ids.get(cb)
        if cid is None:
            cid = len(self.cb_names)
            self._cb_ids[cb] = cid
            self.cb_names.append(cb)
        ends = self.ends
        if self._sorted and ends and key < ends[-1]:
            self._sorted = False
        ends.append(key)
        self.cbids.append(cid)

    def finalize(self) -> None:
        """Restore the non-decreasing ``ends`` invariant ``count_ends`` needs.

        A no-op (and a single flag test) in every geometry that keeps the
        invariant while streaming, which includes all of v2 — so this cannot
        move a v2 byte.  Where it does fire, the sort is STABLE so that equal
        coordinates keep arrival order and the per-cell dict built by
        ``count_ends`` is insertion-ordered deterministically.
        """
        if self._sorted or not self.ends:
            return
        import numpy as np
        # Round-trip through BYTES, never through a Python list.  `.tolist()`
        # on a 19.8 M-element int32 array materialises 19.8 M int objects
        # (~28 B each) plus the list of pointers -- measured at +1.1 GB on the
        # chr19 (+) worker of the PBMC slice, which on its own pushed peak RSS
        # to 2.21x v2 and through the 1.5x compute guard rail of
        # manuscript/24 3.2.4.  frombuffer/tobytes keeps the transient at
        # ~24 B per read (int64 argsort + two int32 copies + two array copies).
        ends = np.frombuffer(self.ends, dtype=np.int32)
        cbids = np.frombuffer(self.cbids, dtype=np.int32)
        order = np.argsort(ends, kind="stable")
        new_ends = array("i")
        new_ends.frombytes(ends[order].tobytes())
        new_cbids = array("i")
        new_cbids.frombytes(cbids[order].tobytes())
        self.ends = new_ends
        self.cbids = new_cbids
        self._sorted = True

    @property
    def sites(self) -> dict:
        return self.accum.sites

    def __len__(self) -> int:
        return len(self.ends)

    # -- queries -------------------------------------------------------------
    def count_ends(self, lo: int, hi: int, into: dict) -> int:
        """Add per-CB counts of read ends with ``lo <= end1 <= hi`` into
        *into*; returns the number of reads added."""
        if hi < lo or not self.ends:
            return 0
        i = bisect_left(self.ends, lo)
        j = bisect_right(self.ends, hi)
        if j <= i:
            return 0
        names = self.cb_names
        for cid, n in Counter(self.cbids[i:j]).items():
            cb = names[cid]
            into[cb] = into.get(cb, 0) + n
        return j - i

    # -- pickling (pipeline hand-off) ---------------------------------------
    def __getstate__(self):
        return (self.accum.sites, self.ends, self.cbids, self.cb_names,
                self.seq_len, self.geometry, self.direction, self._sorted)

    def __setstate__(self, state):
        (sites, self.ends, self.cbids, self.cb_names, self.seq_len,
         self.geometry, self.direction, self._sorted) = state
        self.accum = ClipAccumulator()
        self.accum.sites = sites
        self._cb_ids = {}


class ClipSeeder:
    """Two-tier PAS emission for the ``clip_seeded`` strategy (Phase 2).

    Collects, per chromosome:

    * every qualifying clip read (via :meth:`add_clip`), regardless of
      whether a coverage peak fired — this is the point of seeding: ~60% of
      the clip evidence lies outside every coverage-peak window;
    * every accepted read end (via :meth:`add_read`) — the count source for
      tier-1 PAS that sit outside every coverage candidate;
    * every coverage-strategy PAS candidate together with the ``Peak`` it
      came from (via :meth:`add_coverage_pas`).

    :meth:`flush` then produces the merged, coordinate-sorted record list:

    * TIER 1 — clip clusters clearing ``min_umis`` distinct molecules,
      emitted as the 1-bp interval ``[mode, mode+1)`` with the cluster's
      DISTINCT MOLECULE count in the BED score column (a SUPPORT
      annotation, in the same unit the gate uses; the raw clip-read count
      and the ``-F 3844`` counts go to the ``pas_support.tsv`` sidecar via
      ``flush(support_out=...)``).
      The count-matrix column is NOT the clip reads (that was the Stage-2
      regression: ~1% of the reads, a ~10x mass loss) but every accepted
      read end that piles up at the cleavage site:

      - a cluster that suppresses one or more coverage candidates takes
        those candidates' accumulated ``cb_positions`` counts, restricted to
        its partition of each candidate.  A candidate claimed by several
        clusters is divided among them by midpoint splitting between the
        clusters' anchors — the same
        :func:`~ema.strategies.utils.partition_peak_region` rule
        ``find_pas`` applies to multi-PAS peaks — each share counted with
        :func:`~ema.strategies.utils.reconstruct_cb_dict` semantics on the
        candidate's own positions (what every positional coverage
        strategy's ``get_cb_dict_for_pas`` delegates to; a non-positional
        strategy such as ``original`` is detected by the shares not summing
        to the candidate and the candidate then goes whole to the nearest
        cluster), so the candidate's mass is conserved exactly (shipped
        counting, re-keyed to the cleavage site);
      - every cluster additionally counts the read ends within the count
        window around its cleavage site, ``[site - up, site + down]`` in
        transcript orientation (``--polya-count-window``, default
        ``auto,25`` == ``[site - seq_len, site + 25]``: R2 3' ends pile up
        just upstream of cleavage), EXCLUDING every coverage candidate
        interval and clipped at the midpoint to the neighbouring tier-1
        anchors — so no read end is ever counted twice.  For a cluster
        outside every candidate this is its whole count; for a cluster
        inside a peak it recovers the reads the coverage caller counted
        into ``cb_positions`` but never put in a candidate (the streaming
        loop drops the ends still in ``data_array`` when a peak closes, so
        ``peak_list`` — and with it every candidate interval — stops
        ~``seq_len`` bp short of the pile's 3' end; on ``+`` those are
        exactly the reads AT the cleavage site);
      - a cluster that ends up with no reads from either source (possible
        only on degenerate geometries) falls back to its own clip reads so
        a tier-1 row is never empty.

    ``stats`` accumulates the mass accounting across flushes (reads taken
    from candidate partitions / from windows / rows on the clip fallback,
    and the tier sizes) so a run can report where every count came from.

    * TIER 2 — coverage PAS that overlap no kept cluster, emitted unchanged
      (their ``cb_dict`` straight from ``get_cb_dict_for_pas``) with their
      clip-read support near the strand-aware 3' base as the score.  A peak
      counts as overlapping when a kept cluster's mode falls in
      ``[start - window, end + window]`` **or** any of that cluster's member
      sites lies within ``+/-window`` of the peak's 3' base — the second
      clause is what makes the score an exact tier tag: at the default
      ``min_umis=1`` no cluster is ever rejected, so every surviving tier-2
      peak scores exactly 0.  ``score == 0`` therefore means coverage_only,
      and pasbed.bed stays plain BED6.  (Raise ``min_umis`` and a tier-2
      peak may score >0 — read that as "has clip evidence, but not enough
      for a call of its own", which is exactly what ``--polya-mode filter``
      should keep.)

    Coordinates and tier tags are a pure function of the clip sites and the
    coverage candidates — the counting change above never moves a PAS.
    Nor does the support unit: a cluster's mode is read-weighted over every
    accepted clip read whatever ``--polya-clip-filter`` is set to, so
    switching the unit can only drop a cluster (when its molecule count
    falls below ``min_umis``), never move one.

    Coordinate spaces: the peak caller tracks each read by ``end1``
    (``start1 + seq_len``), so on the minus strand a read whose transcript
    3' end is the cleavage site ``s`` sits at ``end1 == s + seq_len``.  All
    partitioning and windowing therefore happens in ``end1`` space, with a
    cluster's *anchor* being ``mode`` (``+``) or ``mode + seq_len`` (``-``).

    The caller owns pasnumber assignment and the actual writes, so the same
    class serves the monolithic loop, the 3-stage pipeline writer, and the
    tile workers.
    """

    __slots__ = ("direction", "seed_window", "min_umis", "window",
                 "count_up", "count_down", "stream", "coverage", "stats",
                 "strict", "clip_filter", "_cb_ids", "_cb_names",
                 "geometry", "seq_len", "features")

    def __init__(self, direction: bool, seed_window: int = 25,
                 min_umis: int = 1, window: int = 100,
                 count_window=(-1, 25), clip_filter: str = "none",
                 min_reads: int | None = None, *,
                 geometry: str = "fixed", seq_len: int | None = None,
                 features: bool = False) -> None:
        """
        Args:
            min_umis: Minimum DISTINCT MOLECULES for a clip cluster to be
                emitted as a tier-1 PAS (``--polya-min-umis``).
            clip_filter: ``"none"`` (default) counts every accepted clip
                read's molecule; ``"f3844"`` counts only molecules seen on
                reads passing :func:`clip_read_ok`, which also drops
                clusters whose evidence is entirely secondary /
                supplementary / duplicate / qcfail alignments.
            min_reads: Deprecated alias for *min_umis* (the flag was
                ``--polya-min-reads`` and always gated on molecules).
            geometry: ``--read-geometry`` (see
                :data:`ema.countmatrix.read.READ_GEOMETRIES`).  Anything but
                the v2 ``"fixed"`` breaks the ``end1 - start1 == seq_len``
                invariant, so the stream is told the geometry and the
                configured *seq_len* instead of inferring them.
            seq_len: ``--seq-len``.  Required when *geometry* is not
                ``"fixed"``; ignored (learned from the reads, exactly as v2
                does) when it is.
            features: ``--pas-features on`` (peakAtail-prime).  Adds
                ``clip_positions`` / ``clip_span`` to every support dict
                :meth:`flush` produces.  False (v2) leaves both keys absent,
                and :func:`ema.countmatrix.paswrite.support_write` then
                writes exactly v2's seven columns.
        """
        if min_reads is not None:
            min_umis = min_reads
        if clip_filter not in CLIP_FILTERS:
            raise ValueError(
                f"clip_filter must be one of {CLIP_FILTERS}, got {clip_filter!r}"
            )
        self.direction = direction
        self.seed_window = seed_window
        self.min_umis = min_umis
        self.clip_filter = clip_filter
        self.strict = clip_filter == "f3844"
        self.window = window
        self.count_up, self.count_down = parse_count_window(count_window)
        self.geometry = geometry
        self.seq_len = seq_len
        self.features = bool(features)
        self.stream = self._new_stream()
        self.coverage: list[tuple[int, int, dict, _CandidateSlice | None]] = []
        # cell-barcode intern table for the candidate slices (the stream has
        # its own: in the pipeline the writer's stream is replaced by the
        # finder's at flush time, so the two tables must not be shared)
        self._cb_ids: dict = {}
        self._cb_names: list = []
        self.stats = {
            "tier1": 0, "tier2": 0, "suppressed_candidates": 0,
            "reads_from_candidates": 0, "reads_from_window": 0,
            "reads_tier2": 0, "clip_fallback_rows": 0,
            "reads_from_clip_fallback": 0, "empty_tier1_rows": 0,
        }

    # -- streaming hooks ---------------------------------------------------
    def add_clip(self, site: int, cb, umi=None, end1=None,
                 primary: bool = True, start1: int | None = None) -> None:
        """Record one clip read.

        *end1* is stored so the tier-1 clip fallback can tell which cluster's
        midpoint territory the read belongs to, and it must therefore live in
        the SAME space as the read ends :meth:`add_read` stores.  Pass
        *start1* as well and the seeder maps it through
        :meth:`ClipStream.three_key`; under v2 geometry the mapping is the
        identity, so v2 callers may keep omitting it.
        """
        if start1 is not None and self.geometry != "fixed" and end1 is not None:
            end1 = self.stream.three_key(start1, end1)
        self.stream.add_clip(site, cb, umi, end1, primary)

    def add_read(self, start1: int, end1: int, cb) -> None:
        self.stream.add_read(start1, end1, cb)

    def add_coverage_pas(self, pas_1: int, pas_2: int, cb_dict: dict,
                         peak=None) -> None:
        """Buffer one coverage candidate.  *peak* (the ``Peak`` it was found
        in) is what lets several suppressing clusters split the candidate's
        positional counts; its ``cb_positions`` are compacted to the
        candidate's range immediately (the ``Peak`` is not retained).
        Without it the whole ``cb_dict`` goes to the nearest suppressing
        cluster."""
        view = None
        cb_positions = getattr(peak, "cb_positions", None)
        if cb_positions:
            view = _CandidateSlice(cb_positions, min(pas_1, pas_2),
                                   max(pas_1, pas_2), self._intern)
        self.coverage.append((pas_1, pas_2, cb_dict, view))

    def _intern(self, cb) -> int:
        cid = self._cb_ids.get(cb)
        if cid is None:
            cid = len(self._cb_names)
            self._cb_ids[cb] = cid
            self._cb_names.append(cb)
        return cid

    @property
    def min_reads(self) -> int:
        """Deprecated alias of :attr:`min_umis` (the gate was always in
        molecules; only the flag name and the BED score column lied)."""
        return self.min_umis

    @property
    def accum(self) -> ClipAccumulator:
        return self.stream.accum

    def load_sites(self, sites: dict) -> None:
        """Adopt a pre-built site map (no read ends: window counts are empty)."""
        self.stream.accum.sites = sites

    def load_stream(self, stream: ClipStream) -> None:
        """Adopt a finder-built :class:`ClipStream` (pipeline writer)."""
        self.stream = stream

    # -- emission ------------------------------------------------------------
    def flush(self, support_out: list | None = None
              ) -> list[tuple[int, int, int, dict]]:
        """Emit and reset. Returns ``[(bed_start, bed_end, score, cb_dict)]``
        sorted by coordinate.

        ``score`` (BED column 5) is DISTINCT MOLECULES — the unit
        ``--polya-min-umis`` gates on — not raw clip reads.

        Args:
            support_out: Optional list; when given, one support dict per
                emitted record is appended IN THE SAME ORDER (keys
                ``clip_reads``, ``clip_umis``, ``clip_reads_f3844``,
                ``clip_umis_f3844``, ``window_reads``, ``tier``, and -- only
                when this seeder was built with ``features=True`` --
                ``clip_positions``, ``clip_span``).  This is the sidecar
                (``pas_support.tsv``) source: the raw clip-read count stays
                available without overloading BED6.
        """
        stream = self.stream
        # `count_ends` bisects; outside v2 geometry the arrival order is not
        # guaranteed sorted.  No-op (one bool test) whenever it already is.
        stream.finalize()
        sites = stream.sites
        strict = self.strict
        umi_key = "numis_f3844" if strict else "numis"
        clusters = cluster_clip_sites(sites, self.seed_window)
        kept = [c for c in clusters if c[umi_key] >= self.min_umis]
        n_kept = len(kept)

        stats = self.stats
        seq_len = stream.seq_len or 0
        shift = seq_len if self.direction else 0
        anchors = [c["mode"] + shift for c in kept]      # end1 space, sorted
        counts: list[dict] = [{} for _ in kept]
        records: list[tuple[int, int, int, dict]] = []
        supports: list[dict] = []
        cov_intervals: list[tuple[int, int]] = []

        if self.coverage:
            modes = [c["mode"] for c in kept]  # already sorted
            support = _SupportIndex(sites)
            # Sites belonging to KEPT clusters: a coverage peak drawing its
            # support from one of these is tier-1 territory, so suppressing
            # on them keeps `score == 0` an exact coverage_only tag.
            kept_site_map = {s: sites[s] for c in kept for s in c["sites"]}
            kept_sites = _SupportIndex(kept_site_map)
            owner = {s: i for i, c in enumerate(kept) for s in c["sites"]}
            kept_pos = kept_sites.pos
            w = self.window
            for pas_1, pas_2, cb_dict, view in self.coverage:
                bed_start = min(pas_1, pas_2)
                bed_end = max(pas_1, pas_2)
                cov_intervals.append((bed_start, bed_end))
                three = bed_start if self.direction else bed_end - 1
                # suppressed when a kept cluster's mode falls within the
                # window-expanded interval...
                claim = set(range(bisect_left(modes, bed_start - w),
                                  bisect_right(modes, bed_end + w)))
                # ...or when a kept cluster's own sites are what would be
                # scored as this peak's support
                if kept_sites.reads_within(three, w):
                    a = bisect_left(kept_pos, three - w)
                    b = bisect_right(kept_pos, three + w)
                    claim.update(owner[kept_pos[k]] for k in range(a, b))
                if not claim:
                    _score = support.molecules_within(three, w, strict)
                    records.append((bed_start, bed_end, _score, cb_dict))
                    _row2 = {
                        "clip_reads": support.reads_within(three, w),
                        "clip_umis": support.molecules_within(three, w),
                        "clip_reads_f3844": support.reads_within(three, w, True),
                        "clip_umis_f3844": support.molecules_within(three, w, True),
                        "window_reads": sum(cb_dict.values()),
                        "tier": 2,
                    }
                    if self.features:
                        # A tier-2 row has no cluster of its own, so its clip
                        # geometry is that of the clip positions inside the
                        # SAME +/-window its four clip counts come from.
                        _n_pos, _span = support.sites_within(three, w)
                        _row2["clip_positions"] = _n_pos
                        _row2["clip_span"] = _span
                    supports.append(_row2)
                    stats["reads_tier2"] += sum(cb_dict.values())
                    continue
                stats["suppressed_candidates"] += 1
                stats["reads_from_candidates"] += sum(cb_dict.values())
                self._assign_candidate(sorted(claim), three, cb_dict, view,
                                       anchors, counts)

        # Every cluster also counts the read ends in its cleavage window that
        # belong to no coverage candidate (clipped to its neighbours).
        cov_starts, cov_ends = _merge_intervals(cov_intervals)
        if n_kept and len(stream):
            up = self.count_up if self.count_up >= 0 else seq_len
            down = self.count_down
            for i in range(n_kept):
                a = anchors[i]
                if self.direction:
                    lo, hi = a - down, a + up
                else:
                    lo, hi = a - up, a + down
                lo, hi = self._clip_to_neighbours(i, lo, hi, anchors, n_kept)
                for s, e in _subtract_intervals(lo, hi, cov_starts, cov_ends):
                    stats["reads_from_window"] += stream.count_ends(s, e, counts[i])

        for i, c in enumerate(kept):
            cb_dict = counts[i]
            if not cb_dict:
                cb_dict = self._clip_fallback_counts(
                    i, c, sites, anchors, n_kept, cov_starts, cov_ends
                )
                if cb_dict:
                    stats["clip_fallback_rows"] += 1
                    stats["reads_from_clip_fallback"] += sum(cb_dict.values())
                else:
                    stats["empty_tier1_rows"] += 1
            records.append((c["mode"], c["mode"] + 1, c[umi_key], cb_dict))
            _row1 = {
                "clip_reads": c["nreads"],
                "clip_umis": c["numis"],
                "clip_reads_f3844": c["nreads_f3844"],
                "clip_umis_f3844": c["numis_f3844"],
                "window_reads": sum(cb_dict.values()),
                "tier": 1,
            }
            if self.features:
                # cluster_clip_sites() already carries the member positions,
                # ascending; this is a len() and a subtraction, not a rescan.
                _members = c["sites"]
                _row1["clip_positions"] = len(_members)
                _row1["clip_span"] = _members[-1] - _members[0]
                # Clip-anchored cleavage offset (TASK E): where the poly(A)
                # tails actually pinned the last aligned base, relative to the
                # base this cluster REPORTS, in transcript orientation
                # (positive = downstream of the call).  Read-weighted over the
                # cluster's own members, which are already in hand -- one pass
                # over a list whose median length is 1.  A single-member
                # cluster is exactly 0 by construction.
                if len(_members) == 1:
                    _row1["clip_offset_mean"] = "0.00"
                else:
                    _mode = c["mode"]
                    _num = _den = 0
                    for _p in _members:
                        _w = sites[_p][0]
                        _num += (_mode - _p if self.direction else _p - _mode) * _w
                        _den += _w
                    _row1["clip_offset_mean"] = "%.2f" % (_num / _den) if _den else "0.00"
            supports.append(_row1)
        stats["tier1"] += n_kept
        stats["tier2"] += len(records) - n_kept

        order = sorted(range(len(records)), key=lambda k: (records[k][0], records[k][1]))
        records = [records[k] for k in order]
        if support_out is not None:
            support_out.extend(supports[k] for k in order)
        self.stream = self._new_stream()
        self.coverage = []
        self._cb_ids = {}
        self._cb_names = []
        return records

    def _new_stream(self) -> ClipStream:
        """A stream for this seeder's geometry (v2 keeps the bare default)."""
        if self.geometry == "fixed":
            return ClipStream()
        return ClipStream(seq_len=self.seq_len, geometry=self.geometry,
                          direction=self.direction)

    @staticmethod
    def _clip_to_neighbours(i: int, lo: int, hi: int, anchors: list[int],
                            n_kept: int) -> tuple[int, int]:
        """Clip ``[lo, hi]`` at the midpoints to the neighbouring anchors so
        two clusters never claim the same read end."""
        a = anchors[i]
        if i > 0:
            lo = max(lo, (anchors[i - 1] + a) // 2 + 1)
        if i + 1 < n_kept:
            hi = min(hi, (a + anchors[i + 1]) // 2)
        return lo, hi

    def _clip_fallback_counts(self, i: int, cluster: dict, sites: dict,
                              anchors: list[int], n_kept: int,
                              cov_starts: list[int], cov_ends: list[int]) -> dict:
        """Last-resort per-CB counts for a tier-1 row that got no read end
        from either the candidate partition or its cleavage window.

        The cluster's own clip reads are used, but ONLY those whose read-end
        coordinate falls in this cluster's midpoint territory and in no
        coverage candidate — a clip read whose ``end1`` sits past the
        midpoint is already counted in the neighbouring cluster's window
        (measured at 13:43135248+ on mouse1: 37 clip reads with ``end1``
        beyond the midpoint, counted twice), and one inside a candidate is
        already in that candidate's partition.  Rows can legitimately come
        out empty: the reads at that cleavage site belong to the neighbour.

        The exclusion needs each clip read's ``end1``.  When none was
        recorded — a site map adopted through :meth:`load_sites`, or a
        hand-built accumulator — the cluster's whole clip-read count is
        used, exactly as before this fix.
        """
        a = anchors[i]
        lo = (anchors[i - 1] + a) // 2 + 1 if i > 0 else None
        hi = (a + anchors[i + 1]) // 2 if i + 1 < n_kept else None
        out: dict = {}
        have_ends = False
        for s in cluster["sites"]:
            rec = sites.get(s)
            if rec is None:
                continue
            if rec[4]:
                have_ends = True
            for end1, cbs in rec[4].items():
                if lo is not None and end1 < lo:
                    continue
                if hi is not None and end1 > hi:
                    continue
                if _point_covered(end1, cov_starts, cov_ends):
                    continue
                _merge_counts(out, cbs)
        if not have_ends:
            return dict(cluster["cb_counts"])
        return out

    def _assign_candidate(self, claim: list[int], three: int, cb_dict: dict,
                          view, anchors: list[int], counts: list[dict]) -> None:
        """Give a suppressed candidate's counts to the clusters claiming it."""
        if not cb_dict:
            return
        nearest = min(claim, key=lambda i: (abs(anchors[i] - three), i))
        if len(claim) == 1 or not view:
            _merge_counts(counts[nearest], cb_dict)
            return
        from ema.strategies.utils import partition_peak_region
        # partition_peak_region only consults min/max of the positions
        regions = partition_peak_region([anchors[i] for i in claim],
                                        [view.pos[0], view.pos[-1]])
        names = self._cb_names
        shares = [view.counts(rs, re, names) for rs, re in regions]
        if sum(sum(d.values()) for d in shares) != sum(cb_dict.values()):
            # A non-positional coverage strategy (e.g. ``original`` returns
            # the whole peak's cb_dict for any interval) cannot be split
            # without inventing reads: keep the mass whole, nearest cluster.
            _merge_counts(counts[nearest], cb_dict)
            return
        for i, d in zip(claim, shares):
            _merge_counts(counts[i], d)


def _sample_head(bam, min_clip, min_purity, barcode_tag, max_reads):
    """v2's sampler: scan from record 0 until *max_reads* CB reads are seen.

    On a coordinate-sorted BAM that is the head of the first contig, which is
    why the estimate it returns is not the file's rate (2.2565 % against a
    0.5364 % truth on PBMC 10k v3).  Kept because v2 output must stay
    reachable, not because it is right.
    """
    n_cb = n_clip = 0
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
    return n_cb, n_clip, {"mode": "head", "contigs": 1, "strata": 1}


def _sample_strided(bam, min_clip, min_purity, barcode_tag, max_reads,
                    strata=CLIP_RATE_STRATA):
    """Spread the sample across the whole BAM, weighted by read density.

    The sample is a set of equally spaced GENOMIC WINDOWS from which EVERY
    read is taken, not a set of equally sized read blocks.  That distinction
    is the whole estimator, and it was arrived at by measurement rather than
    by taste:

      A first implementation gave each stratum an equal READ quota and started
      it at an evenly spaced coordinate.  That weights read-sparse coordinate
      regions equally with read-dense ones, which is not what the file's rate
      is, and it showed: on the PBMC chr19+21 slice (whole-file truth
      0.6169 %) it returned 0.9111 % at 25 strata, 0.8489 % at 100, 0.3988 %
      at 400 and 0.1853 % at 2,500 -- a 3x swing driven by nothing but the
      stratum count, which is the signature of a biased estimator rather than
      a noisy one.

    Taking every read inside a coordinate-uniform window instead makes dense
    regions contribute proportionally more reads, exactly as they do in the
    whole file, so the pooled ratio estimates the same quantity a full scan
    would.  A read is counted when its START is inside the window, so a read
    is never counted twice and long reads are not over-sampled at the edges.

    Window width is chosen so the expected yield is ``max_reads``: with
    sampling fraction ``f = max_reads / total_mapped``, a contig of length L
    sampled in S windows uses ``W = f * L / S``.

    Returns ``(n_cb, n_clip, diagnostics)``, or ``(0, 0, None)`` when the BAM
    has no usable index -- the caller then falls back to the head scan and
    says so.
    """
    try:
        stats = bam.get_index_statistics()
    except (ValueError, AttributeError):      # no index / unindexed stream
        return 0, 0, None
    mapped = [(s.contig, s.mapped) for s in stats if s.mapped > 0]
    total = sum(m for _, m in mapped)
    if not mapped or total <= 0:
        return 0, 0, None

    lengths = dict(zip(bam.references, bam.lengths))
    fraction = min(1.0, float(max_reads) / float(total))
    n_cb = n_clip = 0
    n_contigs = n_windows = 0
    for contig, _m in mapped:
        length = lengths.get(contig, 0)
        if length <= 0:
            continue
        # bp to sample on this contig -- ALWAYS the same fraction of it, so a
        # short, very deep contig (MT: 16.5 kb and often a tenth of a 10x
        # library) cannot swamp the sample.  Letting the minimum width win
        # here took 100 % of MT and of every scaffold and returned 36.2 M
        # reads for a 200,000-read budget, at 0.3492 % against a 0.5364 %
        # truth -- measured on the full PBMC BAM, not reasoned about.
        span = max(1, int(fraction * length))
        k = max(1, int(strata))
        width = max(1, min(span, max(CLIP_RATE_MIN_WINDOW, span // k)))
        k = max(1, span // width)
        n_contigs += 1
        for j in range(k):
            lo = int(length * j / k)
            hi = min(length, lo + width)
            if hi <= lo:
                continue
            n_windows += 1
            for read in bam.fetch(contig, lo, hi):
                if read.reference_start < lo or read.reference_start >= hi:
                    continue          # count each read in exactly one window
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if not read.has_tag(barcode_tag):
                    continue
                n_cb += 1
                if clip_site(read, min_clip, min_purity) is not None:
                    n_clip += 1
    return n_cb, n_clip, {"mode": "strided", "contigs": n_contigs,
                          "strata": n_windows}


#: Set by :meth:`ClipRateCounter.report` so a caller in the SAME process can
#: pick the two integers up without widening ``peak_calling``'s return type --
#: which several paths ignore.  ``chrom_parallel._run_job`` reads it straight
#: after ``peak_calling`` returns and puts the counts in its result dict, so the
#: dispatcher can report one exact run-level rate.
_last_clip_counts: dict = {"n_reads": 0, "n_clip": 0}


def take_last_clip_counts() -> tuple[int, int]:
    """Return and clear the last job's ``(accepted reads, qualifying clips)``."""
    n_reads = int(_last_clip_counts.get("n_reads", 0))
    n_clip = int(_last_clip_counts.get("n_clip", 0))
    _last_clip_counts["n_reads"] = 0
    _last_clip_counts["n_clip"] = 0
    return n_reads, n_clip


def _warn_low_clip_rate(label, rate, n_reads, n_clip) -> None:
    log.warning(
        "LOW POLY(A) CLIP RATE (exact, %s): %.4f%% of %d accepted CB reads "
        "carry a qualifying poly(A) soft clip (%d reads). Expected >~0.3%% for "
        "10x cDNA kept untrimmed; PBMC 10k v3 reference on this denominator: "
        "0.5730%% genome-wide, 0.3669%%-0.8179%% per chromosome. The clip "
        "evidence channel looks destroyed (poly(A) trimmed upstream, unusual "
        "chemistry, or aligner hard-clipping): clip_seeded calls and BED "
        "score-column support values from this BAM are NOT reliable.",
        label, rate * 100.0, n_reads, n_clip,
    )


def report_clip_rate_total(n_reads: int, n_clip: int, label: str) -> float:
    """Report the EXACT run-level clip rate, and warn once if it is too low."""
    if n_reads <= 0:
        return 0.0
    rate = n_clip / n_reads
    if rate < LOW_CLIP_RATE:
        _warn_low_clip_rate(label, rate, n_reads, n_clip)
    else:
        log.info(
            "poly(A) clip rate (exact, %s): %.4f%% of %d accepted CB reads "
            "(%d clips) -- counted during peak calling, not sampled",
            label, rate * 100.0, n_reads, n_clip,
        )
    return rate


class ClipRateCounter:
    """Exact poly(A) clip rate over the reads the caller actually accepted.

    ``--clip-rate-sampling pass``.  The peak-calling loop already calls
    :func:`clip_site` on every read that clears ``read_check``; this class adds
    two integers to that, so the rate it reports is the file's, on the
    caller's own denominator, for no extra I/O and no sampling theory.

    Reported per (contig, strand) job rather than per run because that is the
    unit the caller parallelises over and the unit whose log line already
    exists -- and because the per-chromosome spread is real: 0.3669 % (chr21)
    to 0.8179 % (chr19) on PBMC 10k v3.
    """

    __slots__ = ("n_reads", "n_clip")

    def __init__(self) -> None:
        self.n_reads = 0
        self.n_clip = 0

    def add(self, has_clip: bool) -> None:
        self.n_reads += 1
        if has_clip:
            self.n_clip += 1

    @property
    def rate(self) -> float:
        return (self.n_clip / self.n_reads) if self.n_reads else 0.0

    def report(self, label: str, warn: bool = False) -> float:
        """Log the exact rate for one unit of work.

        ``warn`` is False for a per-(contig, strand) job on purpose: the
        per-chromosome rate legitimately ranges 0.3669 % (chr21) to 0.8179 %
        (chr19) on a healthy PBMC 10k v3 library, and chr21's ``+`` strand
        alone measures 0.2878 % -- under the 0.3 % bar.  Warning per job would
        fire on a perfectly good library.  The alarm belongs to the run-level
        total (:func:`report_clip_rate_total`).
        """
        rate = self.rate
        if not self.n_reads:
            return rate
        _last_clip_counts["n_reads"] = self.n_reads
        _last_clip_counts["n_clip"] = self.n_clip
        if warn and rate < LOW_CLIP_RATE:
            _warn_low_clip_rate(label, rate, self.n_reads, self.n_clip)
        else:
            log.info(
                "poly(A) clip rate (exact, %s): %.4f%% of %d accepted CB reads "
                "(%d clips)", label, rate * 100.0, self.n_reads, self.n_clip,
            )
        return rate


def check_clip_rate(bam_path: str, min_clip: int = 6, min_purity: float = 0.8,
                    barcode_tag: str = "CB", max_reads: int | None = None,
                    sampling: str | None = None):
    """Estimate the poly(A) clip rate over *max_reads* CB reads and warn
    loudly when it is below :data:`LOW_CLIP_RATE`.

    R2 mitigation (plan section 5): a pipeline that trims poly(A) before
    alignment destroys the evidence channel entirely; the caller must say so
    at startup rather than silently emitting an unsupported call set.

    **How the sample is drawn matters more than how big it is.**  v2 read the
    first *max_reads* CB reads of the file; on a coordinate-sorted BAM that is
    the head of chr1, and on PBMC 10k v3 it returns **2.2565 %** where the
    whole-file rate on the same denominator is **0.5364 %** -- a ~4x
    over-estimate, in the direction that would MASK a genuinely destroyed
    channel.  ``sampling="strided"`` (the peakAtail-prime default) allocates
    the same budget across every mapped contig in proportion to its mapped
    reads and across :data:`CLIP_RATE_STRATA` evenly spaced strata inside each
    contig.

    Cached per BAM per process, so tile/pipeline workers and repeat strand
    passes never rescan.

    Args:
        bam_path: BAM to sample.
        min_clip: ``--polya-min-clip``.
        min_purity: ``--polya-min-purity``.
        barcode_tag: Cell-barcode SAM tag.
        max_reads: CB-read budget for the sample.  ``None`` (default) means
            v2's 200,000 for ``"head"`` and :data:`CLIP_RATE_STRIDED_READS`
            for ``"strided"``, which needs more reads for the same confidence
            because it samples clustered rare events.
        sampling: One of :data:`CLIP_RATE_SAMPLINGS`.  ``None`` (default)
            reads ``variable_config.clip_rate_sampling``, whose legacy-global
            default is the BRANCH value; pass an explicit value to pin it.

    Returns the estimated rate (float) or ``None`` when the BAM could not be
    sampled.
    """
    key = str(bam_path)
    if key in _clip_rate_warned:
        return _clip_rate_warned[key]

    if sampling is None:
        from ema.config import variable_config as _vc
        sampling = str(getattr(_vc, "clip_rate_sampling", V2_CLIP_RATE_SAMPLING))
    if sampling not in CLIP_RATE_SAMPLINGS:
        raise ValueError(
            "clip_rate_sampling must be one of %r, got %r"
            % (list(CLIP_RATE_SAMPLINGS), sampling)
        )
    if sampling == "pass":
        # Nothing to sample: the peak-calling pass counts the real thing.
        # See ClipRateCounter and peak_calling()'s per-job report.
        return None

    if max_reads is None:
        max_reads = (CLIP_RATE_STRIDED_READS if sampling == "strided"
                     else CLIP_RATE_HEAD_READS)

    try:
        import pysam
        with pysam.AlignmentFile(key, "rb") as bam:
            diag = None
            if sampling == "strided":
                n_cb, n_clip, diag = _sample_strided(
                    bam, min_clip, min_purity, barcode_tag, max_reads)
                if diag is None:
                    log.warning(
                        "poly(A) clip-rate QC: %s has no usable index, so the "
                        "sample cannot be spread across the file; falling back "
                        "to the head-of-file scan, which OVER-ESTIMATES the "
                        "rate on a coordinate-sorted BAM (2.2565%% vs a "
                        "0.5364%% truth on PBMC 10k v3).", key,
                    )
                    n_cb, n_clip, diag = _sample_head(
                        bam, min_clip, min_purity, barcode_tag,
                        min(max_reads, CLIP_RATE_HEAD_READS))
            else:
                n_cb, n_clip, diag = _sample_head(
                    bam, min_clip, min_purity, barcode_tag, max_reads)
        rate = (n_clip / n_cb) if n_cb else 0.0
    except Exception as exc:  # never fail a run over the QC estimate
        log.warning("poly(A) clip-rate estimate failed for %s: %s", key, exc)
        _clip_rate_warned[key] = None
        return None

    _clip_rate_warned[key] = rate
    # v2's log line, byte-for-byte, when v2's sampler was asked for: the
    # compatibility guarantee covers the run journal too, and identity_check.py
    # compares it.  The diagnostics ride only on the new sampler.
    if diag["mode"] == "head":
        if rate < LOW_CLIP_RATE:
            log.warning(
                "LOW POLY(A) CLIP RATE: %.4f%% of the first %d CB reads in %s "
                "carry a poly(A) soft clip (expected >~0.3%% for 10x cDNA kept "
                "untrimmed; PBMC 10k v3 genome-wide reference: 0.54%% on this "
                "denominator). The clip evidence channel looks destroyed "
                "(poly(A) trimmed upstream, unusual chemistry, or aligner "
                "hard-clipping). clip_seeded calls and BED score-column "
                "support values from this BAM are NOT reliable. NOTE: "
                "--clip-rate-sampling head reads only the HEAD of the file, "
                "which on a coordinate-sorted BAM is one contig.",
                rate * 100.0, n_cb, key,
            )
        else:
            log.info("poly(A) clip rate for %s: %.4f%% of %d sampled CB reads",
                     key, rate * 100.0, n_cb)
        return rate
    if rate < LOW_CLIP_RATE:
        log.warning(
            "LOW POLY(A) CLIP RATE: %.4f%% of %d sampled CB reads in %s carry "
            "a poly(A) soft clip (sampling=%s over %d contig(s)/%d strata; "
            "expected >~0.3%% for 10x cDNA kept untrimmed; PBMC 10k v3 "
            "genome-wide reference: 0.54%% on this denominator, 0.57%% on the "
            "caller's accepted-read denominator). The clip evidence channel "
            "looks destroyed (poly(A) trimmed upstream, unusual chemistry, or "
            "aligner hard-clipping). clip_seeded calls and BED score-column "
            "support values from this BAM are NOT reliable.",
            rate * 100.0, n_cb, key, diag["mode"], diag["contigs"],
            diag["strata"],
        )
    else:
        log.info(
            "poly(A) clip rate for %s: %.4f%% of %d sampled CB reads "
            "(sampling=%s, %d contig(s), %d windows). This is a QC ALARM, not "
            "a measurement: clips are rare and clustered at 3' ends, so a "
            "sample lands within roughly 0.85x-1.3x of the file's rate.",
            key, rate * 100.0, n_cb, diag["mode"], diag["contigs"],
            diag["strata"],
        )
    return rate
