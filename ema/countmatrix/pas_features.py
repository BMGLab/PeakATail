"""Per-site scoring features for the ``pas_support.tsv`` sidecar (peakAtail-prime).

WHAT THIS IS FOR
----------------
``manuscript/23_algorithm_roadmap.md`` §3 Step 1 and ``PRIME_PLAN.md`` Change 2:
the one measured-positive algorithmic change available to this caller is a
**calibrated per-site score used as a re-ranker inside the existing gates**,
replacing only the ``>=2 clip molecules`` threshold.  That score needs per-site
covariates emitted at call time.  This module computes them; it does not score
anything and it never decides whether a PAS is kept.

WHAT THIS DELIBERATELY DOES **NOT** DO
--------------------------------------
The measurement programme (``results/algo_headroom/VERIFY/VERIFICATION_REPORT.md``,
which supersedes the raw A1-A5 reports wherever they disagree) closed several
doors, and this module stays on the right side of all of them:

* **No second BAM pass.**  Distinct clip cells, top-cell share, end counts at
  +/-5/25/100, pileup sharpness and end-position entropy are worth
  **0.000-0.002 held-out AUC** [V].  They are not emitted and must not be added.
* **No end-position / molecule-end pileup statistics.**  After depth matching, a
  molecule-end pileup at a true missed PAS is as likely as at a random position
  of the same local depth (0.95-1.01x) against 11.36x for the clip channel
  [V §2.1].  The channel carries no information beyond depth.
* **The internal-priming veto stays a HARD GATE.**  ``ip_tool_flag`` is emitted
  here as a *covariate in addition to* the veto, never as a replacement for it:
  every configuration in which a score was allowed to override the veto looked
  spectacular on the atlas and no better on long reads [V §1.4, §7.1].
* **No extra FASTA pass.**  The sequence features ride inside the single pass
  the internal-priming filter already makes over the PAS BED with the genome
  open (``ema.experimental.internal_priming.filter_internal_priming``).  A
  feature that would need its own pass is not emitted; see the module docs in
  ``docs/prime/pas_features.md`` for the list and the reason.

THE ONE COVARIATE THAT MUST BE EXACT
------------------------------------
The verifier found that the tool's **own** ``ip_tool_afrac`` -- the A fraction
over the caller's internal-priming window -- **inverted, alone scores AUC
0.7790** at separating genuine Kinnex long-read termini from internal-priming
decoys, beating the whole A3 model's 0.7655 [V §6(c)].  It is therefore emitted
from the very string the veto tests (``call_internal_priming``'s ``tested``
return value), so it tracks ``--ip-window-left/right`` exactly and cannot drift
away from the rule it summarises.  The ``ip_tool_`` prefix is kept -- rather
than the shorter ``ip_``  -- because that is the name the offline analysis
(``results/algo_headroom/A3_scoring_model/data/FEATURE_TABLE_SCHEMA.tsv``)
joins on, and because it says *whose* window it is.

TRANSCRIPT-RELATIVE COORDINATES
-------------------------------
``r`` is the transcript-relative offset from the cleavage base ``c``:

* ``c`` is the **BED end - 1** on ``'+'`` and the **BED start** on ``'-'`` --
  the same convention :func:`ema.experimental.internal_priming.ip_window` has
  always used (it takes ``pas_pos`` = BED end on ``'+'``, BED start on ``'-'``);
* ``r = 0`` is ``c``; ``r > 0`` is **DOWNSTREAM in transcript orientation**,
  which is where genome-encoded A-stretches cause oligo-dT mispriming, and
  which runs toward **decreasing** genomic coordinate on ``'-'``;
* on ``'-'`` the fetched forward-strand window is reverse-complemented before
  anything is measured -- the strand handling fixed in #96
  (``fix/ip-filter-strand``), applied identically here.

The window fetched per PAS is ``r in [-40, +30]`` (71 nt), the union of the
hexamer window (``-40..-5``) and the downstream-A windows (``+1..+18``,
``+1..+30``).  Nothing emitted here needs a wider one.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right

# --- window geometry -------------------------------------------------------
#: Transcript-relative span fetched per PAS: hexamer window (-40..-5) union the
#: downstream-A windows (+1..+18, +1..+30).  Widening it costs no extra FASTA
#: pass (it is the same fetch), but every column below fits inside it.
FEATURE_R_MIN = -40
FEATURE_R_MAX = 30

#: Value written for a feature that could not be computed (no genome FASTA, a
#: contig missing from it, or a PAS too close to a contig edge for the window).
#: Chosen to match the offline feature tables so a join needs no translation.
NA = "NA"

# --- sequence constants ----------------------------------------------------
#: The two strong canonical poly(A) signals.
HEX_STRONG: tuple[str, ...] = ("AATAAA", "ATTAAA")
#: Beaudoing et al. (2000) canonical hexamer set, in the order the offline
#: feature table uses.  ``hex_any12`` / ``hex_n_types`` count over this set.
HEX12: tuple[str, ...] = (
    "AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AATACA", "CATAAA",
    "GATAAA", "AATATA", "AATAGA", "AAAAAG", "ACTAAA", "AAGAAA",
)
#: Hexamer search window, transcript-relative, inclusive.
HEX_R_LO, HEX_R_HI = -40, -5
#: The Kinnex internal-priming rule's numerator window (``>=12 A of 18``).
KIN_R_LO, KIN_R_HI = 1, 18
#: Kinnex internal-priming threshold: A count in ``+1..+18``.
KIN_A_MIN = 12

_RC_TABLE = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")

#: Sentinel distance for "no same-strand neighbour on this contig", matching
#: the offline feature table's ``d_prev_cand`` / ``d_next_cand`` encoding.
NO_NEIGHBOUR = 1_000_000

# --- column groups ---------------------------------------------------------
# APPEND-ONLY.  These tuples are a file format: a reader may rely on the
# position of every column that already exists, so new features go on the END
# of the last group and nothing is ever reordered or removed.

#: Written by the caller itself, per emitted PAS (no genome, no second pass).
CALL_FEATURE_COLUMNS: tuple[str, ...] = (
    "clip_positions",   # distinct poly(A) clip POSITIONS backing this PAS
    "clip_span",        # bp between the first and last of them (0 if <2)
)

#: Written at the internal-priming seam, from the genome FASTA.
SEQ_FEATURE_COLUMNS: tuple[str, ...] = (
    "seq_ok",           # 1 = the whole r in [-40,+30] window was readable
    "ip_tool_flag",     # the caller's OWN internal-priming call (the veto)
    "ip_tool_afrac",    # A fraction over the caller's IP window
    "ip_tool_arun",     # longest A run over the caller's IP window
    "a_count_d18",      # A count in r +1..+18
    "a_frac_d18",       # = a_count_d18 / 18
    "a_run_d18",        # longest A run in r +1..+18
    "a_frac_d30",       # A fraction in r +1..+30
    "a_run_d30",        # longest A run in r +1..+30
    "kin_ip_flag",      # Kinnex-style indicator: a_count_d18 >= 12
    "hex_strong",       # AATAAA|ATTAAA present in r -40..-5
    "hex_any12",        # any canonical hexamer present in r -40..-5
    "hex_n_types",      # how many distinct canonical hexamers matched
    "hex_best_off",     # r of the last base of the 3'-most hit (0 = none)
    "hex_strong_off",   # same, restricted to AATAAA|ATTAAA (0 = none)
)

#: Written at the same seam, from the PAS BED alone (no genome needed).
CONTEXT_FEATURE_COLUMNS: tuple[str, ...] = (
    "d_prev_cand",      # bp to the previous same-strand candidate
    "d_next_cand",      # bp to the next same-strand candidate
    "n_cand_100",       # OTHER same-strand candidates within +/-100 bp
    "n_cand_500",       # OTHER same-strand candidates within +/-500 bp
    "mol_500_sum",      # sum of BED score over the +/-500 bp neighbourhood
    "is_local_mol_max", # this candidate holds that neighbourhood's max
    "mol_frac_local",   # this candidate's share of mol_500_sum
)

#: Everything the seam appends, in file order.
SEAM_FEATURE_COLUMNS: tuple[str, ...] = SEQ_FEATURE_COLUMNS + CONTEXT_FEATURE_COLUMNS


# ---------------------------------------------------------------------------
# small sequence helpers
# ---------------------------------------------------------------------------
def reverse_complement(seq: str) -> str:
    """Reverse complement (IUPAC letters other than ACGTUN pass through).

    Identical to :func:`ema.experimental.internal_priming.reverse_complement`;
    duplicated as a one-liner rather than imported so this module has no
    dependency on ``ema.experimental``.
    """
    return seq.translate(_RC_TABLE)[::-1]


def longest_run(seq: str, base: str = "A") -> int:
    """Length of the longest run of *base* in *seq*."""
    best = cur = 0
    for ch in seq:
        if ch == base:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 0
    return best


def feature_window(pas_pos: int, strand: str,
                   r_min: int = FEATURE_R_MIN,
                   r_max: int = FEATURE_R_MAX) -> tuple[int, int]:
    """Genomic half-open ``[start, end)`` covering transcript offsets
    ``r_min..r_max`` around one PAS.

    *pas_pos* is the BED ``end`` on ``'+'`` and the BED ``start`` on ``'-'``,
    i.e. exactly what :func:`ema.experimental.internal_priming.ip_window`
    takes, so the two windows are anchored on the same base.  Expressed in
    that function's terms this is ``ip_window(pas_pos, strand,
    window_left=1 - r_min, window_right=r_max)``.  ``start`` is clamped at 0;
    the caller's slice clamps ``end`` at the contig length.
    """
    left = 1 - r_min
    if strand == "+":
        return max(0, pas_pos - left), pas_pos + r_max
    return max(0, pas_pos - r_max), pas_pos + left


def orient_window(window_seq: str, pas_pos: int, strand: str,
                  win_start: int, r_min: int = FEATURE_R_MIN) -> tuple[str, int]:
    """Turn a fetched FORWARD-strand window into transcript orientation.

    Args:
        window_seq: the bases actually returned for ``[win_start, ...)`` --
            possibly shorter than requested if the contig ended.
        pas_pos: BED ``end`` on ``'+'``, BED ``start`` on ``'-'``.
        strand: ``'+'`` or ``'-'``.
        win_start: the (already 0-clamped) genomic start the bases came from.
        r_min: the transcript offset the *requested* window started at.

    Returns:
        ``(tseq, r_first)`` -- the upper-cased transcript-oriented sequence and
        the transcript offset ``r`` of its FIRST base.  ``r_first`` is derived
        from the bases actually returned, so clamping at the contig start
        (``'+'``) or truncation at the contig end (``'-'``) is accounted for
        rather than assumed away.
    """
    seq = window_seq.upper()
    n = len(seq)
    if strand == "+":
        # r of forward index i is (win_start + i) - c, with c = pas_pos - 1.
        return seq, win_start - (pas_pos - 1)
    # '-': transcript index j is forward index n-1-j; c = pas_pos.
    return reverse_complement(seq), pas_pos - win_start - n + 1


def _sub(tseq: str, r_first: int, lo: int, hi: int) -> str | None:
    """``tseq`` over transcript offsets ``lo..hi`` inclusive, or ``None`` if
    the window does not fully cover them."""
    a = lo - r_first
    b = hi - r_first + 1
    if a < 0 or b > len(tseq):
        return None
    return tseq[a:b]


def ip_covariates(tested_seq: str) -> tuple[float, int]:
    """``(A fraction, longest A run)`` over the caller's OWN internal-priming
    window.

    *tested_seq* is the second element of
    :func:`ema.experimental.internal_priming.call_internal_priming`'s return
    value -- already upper-cased and already reverse-complemented on ``'-'``.
    Taking it from there rather than recomputing a window is the point: these
    two numbers then summarise **the exact string the veto tested**, under
    whatever ``--ip-window-left`` / ``--ip-window-right`` the run used.
    """
    if not tested_seq:
        return 0.0, 0
    return tested_seq.count("A") / len(tested_seq), longest_run(tested_seq, "A")


# ---------------------------------------------------------------------------
# sequence features
# ---------------------------------------------------------------------------
def sequence_features(tseq: str, r_first: int) -> dict:
    """The :data:`SEQ_FEATURE_COLUMNS` that come from genomic sequence, minus
    the three ``ip_tool_*`` ones (those come from :func:`ip_covariates`).

    Returns ``{"seq_ok": 0}`` plus :data:`NA` for every value when the window
    does not cover ``r in [FEATURE_R_MIN, FEATURE_R_MAX]`` -- the honest answer
    at a contig edge, and one an offline model can filter on rather than
    silently averaging a truncated denominator.
    """
    hexwin = _sub(tseq, r_first, HEX_R_LO, HEX_R_HI)
    d18 = _sub(tseq, r_first, KIN_R_LO, KIN_R_HI)
    d30 = _sub(tseq, r_first, 1, 30)
    if hexwin is None or d18 is None or d30 is None:
        return _na_sequence_features()

    hits = [h for h in HEX12 if h in hexwin]
    # Offset convention: the transcript offset of the hexamer's LAST base
    # (so always negative, in [-35, -5]); 0 means "no hit".  Among hits the
    # 3'-most -- closest to the cleavage site -- wins.
    best_off = 0
    for h in hits:
        p = hexwin.rfind(h)
        off = HEX_R_LO + p + 5
        if best_off == 0 or off > best_off:
            best_off = off
    strong_off = 0
    for h in HEX_STRONG:
        p = hexwin.rfind(h)
        if p >= 0:
            off = HEX_R_LO + p + 5
            if strong_off == 0 or off > strong_off:
                strong_off = off

    a18 = d18.count("A")
    return {
        "seq_ok": 1,
        "a_count_d18": a18,
        "a_frac_d18": round(a18 / len(d18), 4),
        "a_run_d18": longest_run(d18, "A"),
        "a_frac_d30": round(d30.count("A") / len(d30), 4),
        "a_run_d30": longest_run(d30, "A"),
        "kin_ip_flag": int(a18 >= KIN_A_MIN),
        "hex_strong": int(any(h in HEX_STRONG for h in hits)),
        "hex_any12": int(bool(hits)),
        "hex_n_types": len(hits),
        "hex_best_off": best_off,
        "hex_strong_off": strong_off,
    }


def _na_sequence_features() -> dict:
    out = {c: NA for c in SEQ_FEATURE_COLUMNS}
    out["seq_ok"] = 0
    return out


def na_features() -> dict:
    """Every SEQUENCE column at :data:`NA` with ``seq_ok`` 0 -- the row for a
    PAS whose sequence could not be read at all (contig absent from the
    FASTA).  Distinct from :func:`no_genome_features`, where ``seq_ok`` is
    itself ``NA`` because no FASTA was supplied."""
    out = {c: NA for c in SEQ_FEATURE_COLUMNS}
    out["seq_ok"] = 0
    return out


def no_genome_features() -> dict:
    """Every SEQUENCE column at :data:`NA`, with ``seq_ok`` itself ``NA``.

    Distinguishes "no genome FASTA was supplied, so no sequence feature could
    exist" (``seq_ok = NA``) from "a FASTA was supplied but this site sits too
    close to a contig edge" (``seq_ok = 0``).  Context columns are unaffected:
    they need only the BED.
    """
    return {c: NA for c in SEQ_FEATURE_COLUMNS}


# ---------------------------------------------------------------------------
# local candidate context
# ---------------------------------------------------------------------------
def local_context(rows: list[tuple[str, str, int, int, object]],
                  windows: tuple[int, int] = (100, 500)) -> dict:
    """:func:`context_iter` as a ``{pas_id: {column: value}}`` dict.

    Materialises one small dict per candidate, so it is the readable form for
    tests and library callers and NOT what :class:`FeatureCollector` uses --
    see :meth:`FeatureCollector.finish`.
    """
    return dict(context_iter(rows, windows))


def context_iter(rows: list[tuple[str, str, int, int, object]],
                 windows: tuple[int, int] = (100, 500)):
    """Neighbourhood features for a whole candidate set, one at a time.

    Args:
        rows: ``(chrom, strand, cleavage, molecules, pas_id)`` for **every**
            candidate the seam saw, in any order.  ``cleavage`` is the
            transcript 3' base (BED ``end - 1`` on ``'+'``, BED ``start`` on
            ``'-'``); ``molecules`` is BED column 5.
        windows: the two radii, in bp.

    Yields:
        ``(pas_id, {column: value})`` over :data:`CONTEXT_FEATURE_COLUMNS`,
        in no particular order.  Yielding rather than returning a dict of
        dicts is what keeps the collector's memory flat: a genome-wide run
        has ~650 k candidates and a 7-key dict each would cost ~0.5 GB for
        values that are about to be turned into text anyway.

    The neighbourhood is **same contig, same strand** and is defined over the
    candidate set present at the seam -- which is the caller's full candidate
    set (both tiers, before the internal-priming veto drops anything), the
    same population the offline feature table used.  A ``--polya-mode filter``
    run reaches the seam with the coverage-only tier already gone, so its
    context columns describe the smaller set; that is a property of the run,
    and ``run_config.json`` records it.
    """
    w_small, w_big = windows
    by_key: dict[tuple[str, str], list[tuple[int, int, object]]] = {}
    for chrom, strand, pos, mols, pas_id in rows:
        by_key.setdefault((chrom, strand), []).append((pos, mols, pas_id))

    for group in by_key.values():
        group.sort(key=lambda t: t[0])
        pos = [t[0] for t in group]
        mols = [t[1] for t in group]
        # prefix sums make mol_500_sum O(1) per candidate; the max over the
        # window needs the slice, but a +/-500 bp slice is a handful of rows.
        csum = [0]
        for m in mols:
            csum.append(csum[-1] + m)
        n = len(group)
        for i, (p, m, pas_id) in enumerate(group):
            lo_s = bisect_left(pos, p - w_small)
            hi_s = bisect_right(pos, p + w_small)
            lo_b = bisect_left(pos, p - w_big)
            hi_b = bisect_right(pos, p + w_big)
            mol_sum = csum[hi_b] - csum[lo_b]
            local_max = max(mols[lo_b:hi_b]) if hi_b > lo_b else m
            yield pas_id, {
                "d_prev_cand": (p - pos[i - 1]) if i > 0 else NO_NEIGHBOUR,
                "d_next_cand": (pos[i + 1] - p) if i + 1 < n else NO_NEIGHBOUR,
                # "OTHER candidates": self is inside the slice, so subtract it.
                "n_cand_100": (hi_s - lo_s) - 1,
                "n_cand_500": (hi_b - lo_b) - 1,
                "mol_500_sum": mol_sum,
                "is_local_mol_max": int(m >= local_max),
                "mol_frac_local": round(m / mol_sum, 4) if mol_sum else 0.0,
            }


# ---------------------------------------------------------------------------
# writing
# ---------------------------------------------------------------------------
def format_row(features, columns: tuple[str, ...]) -> str:
    """Tab-joined values for *columns*, :data:`NA` for anything missing.

    *features* may already BE the tab-joined string (what
    :class:`FeatureCollector` stores, so a genome-wide run does not hold a
    24-key dict per candidate); it is then returned unchanged.
    """
    if isinstance(features, str):
        return features
    return "\t".join(str(features.get(c, NA)) for c in columns)


def append_columns(support_path, features: dict,
                   columns: tuple[str, ...] = SEAM_FEATURE_COLUMNS) -> int:
    """Append *columns* to an existing ``pas_support.tsv``, in place.

    APPEND ONLY: every pre-existing column keeps its position and its bytes;
    the header gains the new names and each row gains the new values.  Rows
    whose ``pas_id`` is absent from *features* get :data:`NA` throughout,
    which is what happens to a PAS the seam never saw.

    Returns the number of data rows rewritten.
    """
    from pathlib import Path

    path = Path(support_path)
    tmp = path.with_suffix(path.suffix + ".featuretmp")
    n = 0
    with open(path) as src, open(tmp, "w") as dst:
        header = src.readline().rstrip("\n")
        if not header:
            return 0
        dst.write(header + "\t" + "\t".join(columns) + "\n")
        missing = "\t".join([NA] * len(columns))
        for line in src:
            line = line.rstrip("\n")
            if not line:
                continue
            pas_id = line.split("\t", 1)[0]
            row = features.get(pas_id)
            dst.write(line + "\t"
                      + (missing if row is None else format_row(row, columns))
                      + "\n")
            n += 1
    tmp.replace(path)
    return n


def write_features_tsv(out_path, features: dict,
                       columns: tuple[str, ...] = SEAM_FEATURE_COLUMNS) -> int:
    """Write a standalone ``pas_id`` + *columns* table.

    Used when there is no run-root ``pas_support.tsv`` to extend -- the
    multi-BAM path re-keys its PAS ids in ``merge_pas_beds``, so its per-caller
    sidecars are in a different id space and must not be joined to.
    """
    n = 0
    with open(out_path, "w") as fh:
        fh.write("pas_id\t" + "\t".join(columns) + "\n")
        for pas_id in sorted(features, key=_pas_sort_key):
            fh.write(f"{pas_id}\t{format_row(features[pas_id], columns)}\n")
            n += 1
    return n


def _pas_sort_key(pas_id):
    try:
        return (0, int(pas_id), "")
    except (TypeError, ValueError):
        return (1, 0, str(pas_id))


# ---------------------------------------------------------------------------
# collector
# ---------------------------------------------------------------------------
class FeatureCollector:
    """Accumulates per-PAS features across the pos and neg BEDs.

    One instance is threaded through the internal-priming pass (see
    :func:`ema.experimental.internal_priming.filter_internal_priming`'s
    ``features=`` argument) so the sequence features ride inside the pass that
    already has the genome open -- **there is no second FASTA pass**.  The
    local-context columns need only the BED, so they are computed by
    :meth:`finish` once both strands have been seen.

    The candidate coordinates are kept as a flat list rather than looked up
    again later, so nothing re-reads the BED either.
    """

    __slots__ = ("features", "rows", "n_seq", "n_seq_ok", "_finished")

    #: The sequence block for a candidate with no genome, preformatted once.
    _NO_GENOME_ROW = "\t".join([NA] * len(SEQ_FEATURE_COLUMNS))

    def __init__(self) -> None:
        # {pas_id: tab-joined SEQ_FEATURE_COLUMNS, then SEAM_FEATURE_COLUMNS
        # once finish() has folded the context in}.  A STRING and not a dict:
        # a genome-wide run has ~650 k candidates, and a 24-key dict each
        # costs ~1.5 kB against ~0.15 kB for the text those values are about
        # to become anyway.  Use `parsed()` when you want them by name.
        self.features: dict = {}
        self.rows: list[tuple[str, str, int, int, object]] = []
        self.n_seq = 0        # rows for which sequence features were attempted
        self.n_seq_ok = 0     # ...of which the whole window was readable
        self._finished = False

    def add(self, pas_id, chrom: str, strand: str, cleavage: int,
            molecules: int, feats: dict | None = None) -> None:
        """Record one candidate.  *feats* is the sequence block (or ``None``
        when no genome FASTA was available)."""
        if feats:
            self.n_seq += 1
            if feats.get("seq_ok") == 1:
                self.n_seq_ok += 1
            row = format_row(feats, SEQ_FEATURE_COLUMNS)
        else:
            row = self._NO_GENOME_ROW
        pas_id = str(pas_id)
        self.features[pas_id] = row
        self.rows.append((chrom, strand, cleavage, molecules, pas_id))

    def finish(self) -> dict:
        """Fold the local-context columns in and return
        ``{pas_id: tab-joined SEAM_FEATURE_COLUMNS}``.

        The context is folded in one candidate at a time (:func:`context_iter`)
        and the coordinate list is released, so the collector never holds two
        representations of the same run.  Idempotent.
        """
        if self._finished:
            return self.features
        for pas_id, ctx in context_iter(self.rows):
            self.features[pas_id] += "\t" + format_row(
                ctx, CONTEXT_FEATURE_COLUMNS)
        self.rows = []
        self._finished = True
        return self.features

    def parsed(self, pas_id) -> dict:
        """One candidate's features as ``{column: value}`` (strings).

        The stored form is text; this is the by-name view for tests, logging
        and library callers.  Values are exactly the characters that reach
        ``pas_support.tsv``.
        """
        cols = SEAM_FEATURE_COLUMNS if self._finished else SEQ_FEATURE_COLUMNS
        return dict(zip(cols, self.features[str(pas_id)].split("\t")))

    def stats(self) -> dict:
        return {
            "n_pas": len(self.features),
            "n_sequence_scored": self.n_seq,
            "n_sequence_window_complete": self.n_seq_ok,
        }


def bed_cleavage(start: int, end: int, strand: str) -> int:
    """The transcript 3' base of a BED interval: ``end - 1`` on ``'+'``,
    ``start`` on ``'-'``.  The same base ``r = 0`` refers to."""
    return (end - 1) if strand == "+" else start


def collect_from_bed(bed_path, collector: FeatureCollector) -> int:
    """Add every row of a PAS BED to *collector* with NO sequence features.

    The path taken when ``--pas-features on`` but no ``--genome-fasta`` was
    given: the context columns are still real, the sequence columns are
    :data:`NA`, and ``seq_ok`` is :data:`NA` too (not 0) so a reader can tell
    "no genome was supplied" from "this site is at a contig edge".
    """
    n = 0
    with open(bed_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            start, end, strand = int(parts[1]), int(parts[2]), parts[5]
            try:
                mols = int(parts[4])
            except ValueError:
                mols = 0
            collector.add(parts[3], parts[0], strand,
                          bed_cleavage(start, end, strand), mols, None)
            n += 1
    return n
