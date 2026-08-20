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
from array import array
from bisect import bisect_left, bisect_right
from collections import Counter

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

    ``read_check`` pads/drops reads so that ``end1 - start1 == seq_len`` for
    every accepted read; the BAM is coordinate-sorted, so ``ends`` is
    non-decreasing in arrival order and any ``[lo, hi]`` window is an
    ``O(log n)`` bisect plus a ``Counter`` over the slice.  Memory is 8 bytes
    per accepted read (plus one interned string per distinct barcode), i.e.
    tens of MB for the deepest chromosome of a 200M-read BAM — and it is
    released at every chromosome flush.  ``seq_len`` is learned from the
    first read (the invariant above), so the minus-strand coordinate shift
    below never depends on config plumbing across spawned processes.

    The object is picklable (the 3-stage pipeline hands one per chromosome
    from the finder to the writer); the intern dict is dropped on pickling
    because the receiving side only reads.
    """

    __slots__ = ("accum", "ends", "cbids", "cb_names", "seq_len", "_cb_ids")

    def __init__(self) -> None:
        self.accum = ClipAccumulator()
        self.ends = array("i")
        self.cbids = array("i")
        self.cb_names: list = []
        self.seq_len: int | None = None
        self._cb_ids: dict = {}

    # -- streaming hooks ---------------------------------------------------
    def add_clip(self, site: int, cb, umi=None) -> None:
        self.accum.add(site, cb, umi)

    def add_read(self, start1: int, end1: int, cb) -> None:
        if self.seq_len is None:
            self.seq_len = end1 - start1
        cid = self._cb_ids.get(cb)
        if cid is None:
            cid = len(self.cb_names)
            self._cb_ids[cb] = cid
            self.cb_names.append(cb)
        self.ends.append(end1)
        self.cbids.append(cid)

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
        return (self.accum.sites, self.ends, self.cbids, self.cb_names, self.seq_len)

    def __setstate__(self, state):
        sites, self.ends, self.cbids, self.cb_names, self.seq_len = state
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

    * TIER 1 — clip clusters clearing ``min_reads`` distinct molecules,
      emitted as the 1-bp interval ``[mode, mode+1)`` with the cluster's
      raw clip-read count in the BED score column (a SUPPORT annotation).
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
      ``min_reads=1`` no cluster is ever rejected, so every surviving tier-2
      peak scores exactly 0.  ``score == 0`` therefore means coverage_only,
      and pasbed.bed stays plain BED6.  (Raise ``min_reads`` and a tier-2
      peak may score >0 — read that as "has clip evidence, but not enough
      for a call of its own", which is exactly what ``--polya-mode filter``
      should keep.)

    Coordinates and tier tags are a pure function of the clip sites and the
    coverage candidates — the counting change above never moves a PAS.

    Coordinate spaces: the peak caller tracks each read by ``end1``
    (``start1 + seq_len``), so on the minus strand a read whose transcript
    3' end is the cleavage site ``s`` sits at ``end1 == s + seq_len``.  All
    partitioning and windowing therefore happens in ``end1`` space, with a
    cluster's *anchor* being ``mode`` (``+``) or ``mode + seq_len`` (``-``).

    The caller owns pasnumber assignment and the actual writes, so the same
    class serves the monolithic loop, the 3-stage pipeline writer, and the
    tile workers.
    """

    __slots__ = ("direction", "seed_window", "min_reads", "window",
                 "count_up", "count_down", "stream", "coverage", "stats",
                 "_cb_ids", "_cb_names")

    def __init__(self, direction: bool, seed_window: int = 25,
                 min_reads: int = 1, window: int = 100,
                 count_window=(-1, 25)) -> None:
        self.direction = direction
        self.seed_window = seed_window
        self.min_reads = min_reads
        self.window = window
        self.count_up, self.count_down = parse_count_window(count_window)
        self.stream = ClipStream()
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
        }

    # -- streaming hooks ---------------------------------------------------
    def add_clip(self, site: int, cb, umi=None) -> None:
        self.stream.add_clip(site, cb, umi)

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
    def accum(self) -> ClipAccumulator:
        return self.stream.accum

    def load_sites(self, sites: dict) -> None:
        """Adopt a pre-built site map (no read ends: window counts are empty)."""
        self.stream.accum.sites = sites

    def load_stream(self, stream: ClipStream) -> None:
        """Adopt a finder-built :class:`ClipStream` (pipeline writer)."""
        self.stream = stream

    # -- emission ------------------------------------------------------------
    def flush(self) -> list[tuple[int, int, int, dict]]:
        """Emit and reset. Returns ``[(bed_start, bed_end, score, cb_dict)]``
        sorted by coordinate."""
        stream = self.stream
        sites = stream.sites
        clusters = cluster_clip_sites(sites, self.seed_window)
        kept = [c for c in clusters if c["numis"] >= self.min_reads]
        n_kept = len(kept)

        stats = self.stats
        seq_len = stream.seq_len or 0
        shift = seq_len if self.direction else 0
        anchors = [c["mode"] + shift for c in kept]      # end1 space, sorted
        counts: list[dict] = [{} for _ in kept]
        records: list[tuple[int, int, int, dict]] = []
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
                    records.append((bed_start, bed_end,
                                    support.reads_within(three, w), cb_dict))
                    stats["reads_tier2"] += sum(cb_dict.values())
                    continue
                stats["suppressed_candidates"] += 1
                stats["reads_from_candidates"] += sum(cb_dict.values())
                self._assign_candidate(sorted(claim), three, cb_dict, view,
                                       anchors, counts)

        # Every cluster also counts the read ends in its cleavage window that
        # belong to no coverage candidate (clipped to its neighbours).
        if n_kept and len(stream):
            up = self.count_up if self.count_up >= 0 else seq_len
            down = self.count_down
            cov_starts, cov_ends = _merge_intervals(cov_intervals)
            for i in range(n_kept):
                a = anchors[i]
                if self.direction:
                    lo, hi = a - down, a + up
                else:
                    lo, hi = a - up, a + down
                if i > 0:
                    lo = max(lo, (anchors[i - 1] + a) // 2 + 1)
                if i + 1 < n_kept:
                    hi = min(hi, (a + anchors[i + 1]) // 2)
                for s, e in _subtract_intervals(lo, hi, cov_starts, cov_ends):
                    stats["reads_from_window"] += stream.count_ends(s, e, counts[i])

        for i, c in enumerate(kept):
            cb_dict = counts[i]
            if not cb_dict:
                cb_dict = c["cb_counts"]
                stats["clip_fallback_rows"] += 1
            records.append((c["mode"], c["mode"] + 1, c["nreads"], cb_dict))
        stats["tier1"] += n_kept
        stats["tier2"] += len(records) - n_kept

        records.sort(key=lambda r: (r[0], r[1]))
        self.stream = ClipStream()
        self.coverage = []
        self._cb_ids = {}
        self._cb_names = []
        return records

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
