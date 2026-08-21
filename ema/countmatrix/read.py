from ema.config import variable_config
from ema.countmatrix.cb_encode import encode_cb

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

    Returns:
        A 5-tuple ``(chro, start, end, strand, cb)`` for valid reads, or
        ``(0, 0, 0, 0, 0)`` for reads that should be skipped.
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
    # this block just get reads that len are standard len
    # if it is > just skip read
    # if it is < change to standard
    # TODO: this block will be fixed
    span = read_end - read_start
    if span > seq_len:
        return 0, 0, 0, 0, 0
    elif span < seq_len:
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