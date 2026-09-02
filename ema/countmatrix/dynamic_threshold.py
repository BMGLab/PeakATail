"""The dynamic-threshold look-back index, and the option that bounds it.

WHAT THIS IS FOR
----------------
Both peak-calling loops (:mod:`ema.countmatrix.peackcalling`, the monolithic
loop every parallel worker runs, and :mod:`ema.countmatrix.peak_pipeline`, the
streaming variant) compute the current peak's right edge as::

    l_end = data_array[-current_threshold]

``data_array`` is the live window of read 3' ends; ``current_threshold`` is the
peak-height threshold.  With ``--dynamic-threshold`` the threshold is
recomputed per read from the local lambda::

    current_threshold = max(floor_threshold, int(local_lambda * lambda_fold_change))

and NOTHING bounds it by ``len(data_array)``.  When the threshold outgrows the
live window the index runs off the front of the list and the run dies with a
bare ``IndexError: list index out of range``.

MEASURED, NOT INFERRED (2026-08-23, manuscript/28 Lead D):
  * ``--dynamic-threshold --lambda-fold-change 2.0`` -- the reference's own
    documented value for "find more PAS" -- aborts the PBMC chr19+21 dev slice
    at ``peackcalling.py`` ``l_end = data_array[-current_threshold]``, raised
    through ``chrom_parallel.py``.  The accuracy verifier reproduced the same
    abort on the GSE104556 mouse1 chr18+19 slice from its own runner, so the
    defect is SPECIES-INDEPENDENT: it is not a property of one library.
  * ``--lambda-fold-change`` is read at exactly two places in the tree, both
    inside ``if dynamic_threshold:``, so without ``--dynamic-threshold`` (off
    by default) neither the flag nor this defect is reachable.  That is why
    every published PeakATail number is unaffected.

THE DEFAULT IS ``off`` AND THAT IS DELIBERATE.  ``off`` is v2 (commit
``9dfdefb``) to the character: the same expression, the same ``IndexError``,
the same abort.  This branch's cardinal rule is that v2 stays reproducible
byte-for-byte, and a fix that changes what a run *emits* would have to clear
``manuscript/28`` section 2 first.  Turning the clamp ``on`` cannot change any
run that does not already crash -- see :func:`resolve_l_end`, which is the
single copy of the rule and is exhaustively tested for exactly that in
``tests/test_prime_dynamic_threshold_clamp.py``.
"""
from __future__ import annotations

import logging

log = logging.getLogger(__name__)

#: Branch default for ``--dynamic-threshold-clamp``.  ``False`` == v2.
DYNAMIC_THRESHOLD_CLAMP_DEFAULT = False

#: What a caller of :func:`resolve_l_end` is told when the index goes out of
#: range anyway.  The message names the flag, because a bare ``IndexError``
#: from inside a spawned worker cost this project a whole sweep arm before
#: anyone read the source.  NOTE: the peak-calling loops route through this
#: function ONLY when ``--dynamic-threshold-clamp on`` was asked for -- the
#: default path keeps v2's literal expression, so a default run's abort stays
#: exactly as bare as v2's (byte-for-byte rule; the flag is named in the
#: CHANGELOG and the issue instead).
OUT_OF_RANGE_HINT = (
    "dynamic threshold (%d) outgrew the live read window (%d ends): "
    "l_end = data_array[-%d] is out of range. This is a v2 (9dfdefb) defect, "
    "reproduced on both dev slices and on both species; it is reachable only "
    "with --dynamic-threshold. Re-run with --dynamic-threshold-clamp on to "
    "bound the index, or drop --dynamic-threshold (the default)."
)


def resolve_l_end(data_array, current_threshold: int, clamp: bool):
    """Return the current peak's right edge, optionally bounding the index.

    Args:
        data_array: the live window of read 3' ends, ascending.  Anything
            supporting ``len()`` and integer indexing (the callers pass a
            ``sortedcontainers.SortedList``).
        current_threshold: the active peak-height threshold, i.e. how far back
            from the newest end to look.
        clamp: ``False`` (the branch default) evaluates v2's expression
            unchanged, so an out-of-range index raises ``IndexError`` exactly
            as v2 does -- only with :data:`OUT_OF_RANGE_HINT` logged first.
            ``True`` bounds the index to the window, which returns the OLDEST
            end still in the window (``data_array[0]``).

    Returns:
        The element v2 would have returned whenever v2 could return one.

    Raises:
        IndexError: when ``clamp`` is False and the index is out of range, and
            when ``data_array`` is empty regardless of ``clamp`` (an empty
            window has no oldest end, so clamping has nothing to clamp TO and
            must not invent one).

    THE IDENTITY GUARANTEE: for every ``0 < current_threshold <=
    len(data_array)`` the two branches return the SAME element.  ``clamp=True``
    can therefore only change a run that would otherwise have aborted.
    """
    n = len(data_array)
    if clamp and n and current_threshold > n:
        return data_array[0]
    try:
        return data_array[-current_threshold]
    except IndexError:
        log.error(OUT_OF_RANGE_HINT, current_threshold, n, current_threshold)
        raise


class DynamicThresholdGuard:
    """Stateful wrapper over :func:`resolve_l_end` that counts clamp hits.

    One per peak-calling job.  ``None`` is used in the loops instead of an
    ``off`` guard so the default path keeps v2's literal expression and pays
    nothing for an option it is not using.
    """

    __slots__ = ("clamp", "n_clamped", "worst_threshold", "worst_len")

    def __init__(self, clamp: bool = DYNAMIC_THRESHOLD_CLAMP_DEFAULT) -> None:
        self.clamp = bool(clamp)
        self.n_clamped = 0
        self.worst_threshold = 0
        self.worst_len = 0

    def l_end(self, data_array, current_threshold: int):
        n = len(data_array)
        if self.clamp and n and current_threshold > n:
            self.n_clamped += 1
            if current_threshold - n > self.worst_threshold - self.worst_len:
                self.worst_threshold, self.worst_len = current_threshold, n
        return resolve_l_end(data_array, current_threshold, self.clamp)

    def report(self, where: str = "") -> None:
        """Log the clamp census, or nothing at all when it never fired."""
        if self.n_clamped:
            log.warning(
                "--dynamic-threshold-clamp bounded the look-back index %d "
                "time(s)%s (worst: threshold %d vs window %d). v2 (9dfdefb) "
                "aborts with IndexError here; this run did not.",
                self.n_clamped, (" in %s" % where) if where else "",
                self.worst_threshold, self.worst_len,
            )
