from ema.config import variable_config
from ema.countmatrix.cb_encode import encode_cb

#: Accepted values of ``--read-geometry`` (peakAtail-prime).
#:
#: ``fixed``
#:     **v2 / compatibility.**  Discard any read whose REFERENCE span exceeds
#:     ``--seq-len``; rewrite a shorter read's end to ``start + seq_len``.
#:     Every accepted read is exactly ``seq_len`` bp of reference.
#: ``keep``
#:     Stop discarding: accept on QUERY length instead (a spliced alignment is
#:     no longer thrown away for the length of its intron) but keep v2's
#:     fixed-length interval.  This is the ablation arm -- it isolates "stop
#:     discarding" from "stop fabricating the 3' end".
#: ``true``
#:     Stop discarding **and** stop fabricating: the read's true aligned
#:     reference footprint, soft clips excluded (they are not aligned), D
#:     inside the span, N (intron) removed.
READ_GEOMETRIES = ("fixed", "keep", "true")

#: v2's value, and what :data:`ema.cli.config_schema.V2_COMPAT_FLAGS` -- the
#: documented v2-compatibility command line, not a ``--compat`` flag (there is
#: none) -- pins ``--read-geometry`` to.
V2_READ_GEOMETRY = "fixed"

# Fallback sample_id used when a read has no RG tag. main.py sets this to the
# current dataset_id before each peak_calling call so multi-sample BAMs without
# pre-existing RG tags still produce unique CB prefixes per dataset.
_default_sample_id = "default"

# Interned ``"<RG>_<CB>"`` composites, keyed by (RG, CB) -- one string object
# per cell instead of one per read (see the note at the end of read_check).
_composite_cb: dict[tuple[str, str], str] = {}


def set_default_sample_id(sample_id: str) -> None:
    global _default_sample_id
    _default_sample_id = sample_id


def reset_composite_cache() -> None:
    """Drop the interned composite-barcode table (between datasets/tests)."""
    _composite_cb.clear()


def read_check(
    read,
    direction: bool,
    barcode=None,
    barcode_len=None,
    seq_len=None,
    ignore_chro=None,
    *,
    sample_id: str | None = None,
    geometry: str | None = None,
    exclude_flags: int | None = None,
):
    """Validate a BAM read and extract fields needed for peak calling.

    Called once per read in the hot path (~14 M times per run); keep the
    function body as lean as possible.

    Args:
        read: A ``pysam.AlignedSegment`` object.
        direction: Expected strand direction (``True`` = reverse/negative,
            ``False`` = forward/positive).  Reads on the wrong strand are
            discarded.
        barcode: SAM tag name for the cell barcode (default: CB).
        barcode_len: Expected barcode length in nucleotides (default: 16).
        seq_len: Expected/maximum read length in bp.
        ignore_chro: Set of chromosome names to skip entirely.
        sample_id: Explicit fallback sample identifier used when the read
            has no RG tag.  When ``None`` (default), the module-level
            ``_default_sample_id`` is used -- preserving backward
            compatibility for all existing callers.  New code (e.g. parallel
            workers that know their dataset ID) should pass an explicit value
            so the module singleton is not consulted.
        geometry: Read acceptance geometry, one of :data:`READ_GEOMETRIES`.
            ``None`` (default) reads ``variable_config.read_geometry``, whose
            legacy-global default is v2's ``"fixed"`` -- so every existing
            library caller and unit test keeps v2 behaviour.  The BRANCH
            default lives in the CLI schema (``RunConfig.read_geometry``),
            which is what ``ema run`` bridges into the legacy global.
        exclude_flags: SAM flag mask; a read with any of these bits set is
            skipped (``samtools view -F``).  ``None`` (default) reads
            ``variable_config.read_exclude_flags``, whose default is ``0`` =
            v2 = no flag filtering at all on the coverage channel.

    Returns:
        A 5-tuple ``(chro, start, end, strand, cb)`` for valid reads, or
        ``(0, 0, 0, 0, 0)`` for reads that should be skipped.

    Note:
        Under ``geometry="fixed"`` every accepted read satisfies
        ``end - start == seq_len``; under the other two modes it does not, and
        anything downstream that relied on that invariant (see
        :class:`ema.countmatrix.polya.ClipStream`) must be told the geometry
        explicitly rather than learning ``seq_len`` from the first read.
    """
    # Look up defaults from variable_config at CALL time (not at function
    # definition time) so that the value reflects the live YAML/CLI bridge.
    # Previously these were captured as Python default arguments at module
    # import — meaning a None at import time was frozen forever, causing
    # `pysam.AlignedSegment.get_tag(None)` to crash with
    # "expected bytes, NoneType found" when the user did not pass
    # --barcode-tag.
    if barcode is None:
        barcode = variable_config.barcode_tag
    if barcode_len is None:
        barcode_len = variable_config.cb_len
    if seq_len is None:
        seq_len = variable_config.seqlen
    if ignore_chro is None:
        ignore_chro = variable_config.ignore_chro
    if geometry is None:
        geometry = variable_config.read_geometry
    if exclude_flags is None:
        exclude_flags = variable_config.read_exclude_flags

    # ORDER MATTERS FOR SPEED ONLY: every rejection returns the same
    # sentinel, so the checks may run in any order — and the cheapest,
    # most selective ones run first.  The strand test is a single flag bit
    # and discards roughly half of all records before the ~0.7 us
    # `get_tag(CB)` lookup, the reference_end computation (walks the CIGAR)
    # and the composite-barcode build.  Measured on the PBMC chr19+21
    # slice: ~2.2 us/record less than the CB-tag-first order.
    read_strand = read.is_reverse
    if direction != read_strand:
        return 0, 0, 0, 0, 0

    # Optional SAM-flag veto on the COVERAGE channel (--read-exclude-flags).
    # Default 0 == v2: read_check applies no -F filter, so a read aligned to
    # N places contributes N reads of coverage and N matrix counts.  Kept off
    # by default because it is a call-set change, not a bug fix -- see the
    # measurement in results/prime/read_geometry_slice.tsv.
    if exclude_flags and (read.flag & exclude_flags):
        return 0, 0, 0, 0, 0

    # Unmapped (or CIGAR-less) reads have reference_end None; CellRanger BAMs
    # keep unmapped reads with barcodes, so this must be guarded or the run
    # crashes hours in with a TypeError on `read_end - read_start`.
    if read.is_unmapped:
        return 0, 0, 0, 0, 0
    read_chro = read.reference_name
    read_end = read.reference_end
    if read_end is None or read_chro is None:
        return 0, 0, 0, 0, 0

    # chromosomes the user wants to ignore
    if read_chro in ignore_chro:
        return 0, 0, 0, 0, 0

    read_start = read.reference_start

    # --- read acceptance geometry (peakAtail-prime, --read-geometry) --------
    # v2 ("fixed") normalises every accepted read to exactly `seq_len` bp of
    # REFERENCE span: longer reads are DISCARDED and shorter ones have their
    # end rewritten to `read_start + seq_len`.  Measured on the PBMC chr19+21
    # slice (results/prime/read_geometry_census_pbmc_slice.tsv): the discard
    # removes 11,768,752 / 48,647,964 = 24.19 % of valid-CB reads, 98.08 % of
    # them spliced (genome-wide A2 census: 13.74 %, 96 % spliced), and the
    # rewrite fabricates the 3' end of a further 4,752,308 reads (9.77 %) by
    # a mean +13.53 bp.  See ema/countmatrix/read.py's module docstring for
    # what each mode does; the flag is --read-geometry / yaml `read_geometry`.
    span = read_end - read_start
    if geometry == "fixed":
        if span > seq_len:
            return 0, 0, 0, 0, 0
        elif span < seq_len:
            read_end = read_start + seq_len
    else:
        # "keep" and "true" replace v2's REFERENCE-SPAN test with the same test
        # applied to the read's ALIGNED FOOTPRINT -- the reference span minus
        # its introns (CIGAR N).  A spliced alignment is then no longer thrown
        # away for the length of its intron, which is the whole point.
        #
        # Why the footprint and not the query length (`len(SEQ)`), which is the
        # more obvious reading of "a span-vs-query-length test": the query
        # length test is NOT a strict relaxation of v2.  A short alignment with
        # a long terminal soft clip has span <= seq_len (v2 KEEPS it) but
        # query length > seq_len (a query-length test DROPS it) -- exactly the
        # shape of a poly(A) clip read, the evidence this branch exists to
        # keep.  tests/test_read_geometry.py pins the property that matters:
        # every read v2 accepts, "keep" and "true" accept too.  The footprint
        # test has that property by construction (footprint <= span).
        #
        # The two rules are NOT interchangeable on the PBMC slice either, and
        # in the direction that matters for a *quantification* claim: the
        # census (results/prime/read_geometry_census_pbmc_slice.tsv) shows
        # qlen_max 91 and qlen_gt_seqlen 0, so a query-length rule would accept
        # ALL 11,768,752 reads v2 discards, while the footprint rule accepts
        # 11,526,534 of them and still rejects 242,218 whose real ALIGNED
        # footprint exceeds --seq-len (long deletions and the like).  The
        # footprint rule is the conservative one, and that is deliberate.
        #
        # Soft clips are excluded at BOTH ends: they are not aligned to the
        # reference, and the terminal poly(A) clip in particular is the clip
        # detector's business (ema.countmatrix.polya.clip_site), not the
        # coverage geometry's -- counting it would push a clipped read's 3'
        # end PAST its own cleavage site.  D (deletion) advances the reference
        # and stays inside the footprint.  N (intron) does not: a spliced read
        # does not cover its intron, and measured over 7.82 M valid-CB reads
        # of the PBMC chr19+21 slice the introns carry 10.02 Gb -- about 14x
        # the real read mass over the same 105 Mb -- so admitting them would
        # turn the coverage state machine into a gene-body detector.
        #
        # The CIGAR walk is gated on `span > seq_len`.  That gate is EXACT for
        # the acceptance decision (footprint <= span, so span <= seq_len can
        # never be a rejection) and costs nothing on the ~76 % of slice reads
        # that are unspliced.  It is very slightly approximate for `true`'s
        # end coordinate: a read with span <= seq_len keeps any intron it has
        # inside its interval, which needs an intron shorter than
        # seq_len - (aligned length).  Measured: 1 read in 7,824,805 on the
        # slice, carrying 1 bp.
        if span > seq_len:
            nlen = 0
            cigar = read.cigartuples
            if cigar:
                for _op, _n in cigar:
                    if _op == 3:                # BAM_CREF_SKIP (N)
                        nlen += _n
            if span - nlen > seq_len:
                return 0, 0, 0, 0, 0
            if geometry == "true" and nlen:
                read_end = read_start + (span - nlen)
        if geometry == "keep":
            # v2's fixed-length interval, minus the discard: this arm isolates
            # "stop discarding" from "stop fabricating the 3' end".
            if read_end - read_start != seq_len:
                read_end = read_start + seq_len

    try:  # do not calculate reads don't have CB
        cb = read.get_tag(barcode)
        # CellRanger appends a GEM-group suffix to corrected barcodes
        # ("AAAC...GTT-1"); STARsolo does not. Strip a trailing "-<digits>"
        # so stock CellRanger BAMs work — without this, every read fails the
        # length check below and is silently dropped, producing an empty run.
        if len(cb) != barcode_len:
            dash = cb.rfind("-")
            if dash == barcode_len and cb[dash + 1:].isdigit():
                cb = cb[:dash]
            if len(cb) != barcode_len:
                return 0, 0, 0, 0, 0
    except KeyError:
        return 0, 0, 0, 0, 0

    # NOTE: invalid-char filtering is NOT done here — encode_cb is too slow
    # (~3μs × 14M reads = 40+ sec). BarcodeIndex.get_index handles invalid
    # CBs by encoding them to -1 (tuple key (sample_id, -1)); these end up in
    # the same column but are rare (1-letter Ns are <0.1% of CB sequencing data).

    try:
        rg = read.get_tag('RG')
    except (KeyError, ValueError):
        # Use the caller-supplied sample_id if given; otherwise fall back to
        # the module-level singleton (backward compatibility).
        rg = sample_id if sample_id is not None else _default_sample_id
    # The RG is used VERBATIM.  The composite built below is decoded by
    # ema.countmatrix.indexing.split_cb, which splits on the LAST underscore
    # -- the barcode half is a fixed-length ACGTN string and can never hold
    # one, so any RG (underscores included) round-trips exactly.
    #
    # Do NOT sanitise or drop underscore-bearing RGs here.  `ema merge`
    # stamps RG = dataset_id and `samtools merge` derives RG ids from file
    # names, so `sampleA_rep1` / `sampleB_rep1` are ordinary values; mapping
    # them onto a shared fallback (or onto `_`->`-`, which collides
    # `a_b` with `a-b`) puts two samples' cells in one matrix column, and
    # rewriting the prefix also breaks the per-dataset column selector in
    # ema/main.py and ema/reannotate.py, which matches the sample half of the
    # composite against the dataset id exactly.

    # The composite is interned per distinct (RG, CB) pair: it is built once
    # per cell instead of once per read, and every downstream structure that
    # holds it (peak cb_dicts, ClipStream.cb_names, the clip accumulator)
    # then shares ONE string object per cell rather than one per read —
    # tracemalloc measured 102 MB of duplicate composites alive at a single
    # chromosome flush of the PBMC slice.
    key = (rg, cb)
    composite = _composite_cb.get(key)
    if composite is None:
        composite = f"{rg}_{cb}"
        _composite_cb[key] = composite
    return read_chro, read_start, read_end, read_strand, composite
    
if __name__ == "__main__":
    read_check()