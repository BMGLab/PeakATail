"""The dynamic-threshold look-back index, and the bound that keeps it legal.

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

and nothing used to bound it by ``len(data_array)``.  When the threshold
outgrew the live window the negative index ran off the front of the list and
the run died with a bare ``IndexError: list index out of range`` raised from
inside a spawned chromosome worker (issue #101).

MEASURED, NOT INFERRED (issue #101):
  * ``--dynamic-threshold --lambda-fold-change 2.0`` -- the parameter
    reference's own documented value for "find more PAS" -- aborted the PBMC
    chr19+21 dev slice at ``l_end = data_array[-current_threshold]``.  The
    same abort reproduced on the GSE104556 mouse1 chr18+19 slice, so the
    defect is SPECIES-INDEPENDENT: not a property of one library's depth.
  * ``--lambda-fold-change`` is read at exactly two places in the tree, both
    inside ``if dynamic_threshold:``, and ``--dynamic-threshold`` is off by
    default, so no published PeakATail number is affected.

THE BOUND IS ON BY DEFAULT.  A crash is not a behaviour worth preserving, and
:func:`resolve_l_end` carries an identity guarantee -- for every in-range
threshold the bounded and unbounded branches return the SAME element -- so
bounding can only change a run that would otherwise have aborted.

ONE IMPLEMENTATION, NOT TWO.  ``peakAtail-prime`` (PR #100) landed its own
copy of this module behind ``--dynamic-threshold-clamp``, defaulting to
``off`` so that v2 (commit ``9dfdefb``) stayed reproducible down to its
``IndexError``.  The merge of ``develop`` into that branch consolidated the
two on THIS one: the bound is on by default and no command line, YAML key or
library keyword can turn it off.

``--dynamic-threshold-clamp`` therefore survives only as an ACCEPTED NO-OP,
kept so that the branch's documented v2-compat command line and every
existing config keep parsing.  Both of its values now mean the same thing
(bounded), and so does the ``dynamic_threshold_clamp=`` keyword the peak
loops still accept.  What the branch's byte-for-byte v2 promise still buys is
unchanged and exact, because of the identity guarantee: every v2 run that
PRODUCED OUTPUT produces the same output here.  The only runs that differ are
the ones v2 aborted, and those emitted nothing to be identical to.

The ``clamp`` argument on :func:`resolve_l_end` and
:class:`DynamicThresholdGuard` is retained for the unit tests that pin the
identity guarantee, which have to be able to evaluate the unbounded branch to
prove the bounded one agrees with it.  No production caller passes it.
"""
from __future__ import annotations

import logging

log = logging.getLogger(__name__)

#: Default for the look-back bound.  ``True`` == bounded == no IndexError.
DYNAMIC_THRESHOLD_CLAMP_DEFAULT = True

#: What a caller is told when the index goes out of range anyway (an empty
#: window, or ``clamp=False``).  The message names the flag, the contig and
#: the remedy, because a bare ``IndexError`` from inside a spawned worker
#: cost this project a whole parameter-sweep arm before anyone read the
#: source.
OUT_OF_RANGE_HINT = (
    "dynamic threshold (%d) outgrew the live read window (%d ends): "
    "l_end = data_array[-%d] is out of range%s. This is issue #101, and it "
    "is reachable only with --dynamic-threshold. The look-back index is "
    "bounded on every production path, so reaching this means the live "
    "window was EMPTY; lower --lambda-fold-change, or drop "
    "--dynamic-threshold (the default)."
)


def validate_floor_threshold(floor_threshold: int) -> int:
    """Reject a ``--floor-threshold`` that cannot mean anything.

    ``current_threshold`` is used as ``data_array[-current_threshold]``, so
    ``0`` silently reads ``data_array[-0] == data_array[0]`` -- the OLDEST end
    in the window rather than the newest-minus-N -- and a negative value
    indexes forward from the front.  Neither crashes and neither is what the
    flag says it does, so both are refused here rather than producing a wrong
    peak edge without a word.

    Args:
        floor_threshold: the value to check.

    Returns:
        The value, unchanged, when it is usable.

    Raises:
        ValueError: when it is below 1.
    """
    value = int(floor_threshold)
    if value < 1:
        raise ValueError(
            "--floor-threshold must be >= 1 (got %d): it is used as the "
            "look-back distance data_array[-floor_threshold], so 0 silently "
            "returns the oldest end in the window instead of the peak edge "
            "and a negative value indexes from the wrong end (issue #101)."
            % value
        )
    return value


def resolve_l_end(data_array, current_threshold: int,
                  clamp: bool = DYNAMIC_THRESHOLD_CLAMP_DEFAULT):
    """Return the current peak's right edge, bounding the index by the window.

    Args:
        data_array: the live window of read 3' ends, ascending.  Anything
            supporting ``len()`` and integer indexing (the callers pass a
            ``sortedcontainers.SortedList``).
        current_threshold: the active peak-height threshold, i.e. how far back
            from the newest end to look.
        clamp: ``True`` (the default) bounds the index to the window, which
            returns the OLDEST end still in it (``data_array[0]``).  ``False``
            evaluates the historical expression unchanged, so an out-of-range
            index raises ``IndexError`` exactly as it always did -- only with
            :data:`OUT_OF_RANGE_HINT` logged first.

    Returns:
        The element the unbounded expression would have returned whenever it
        could return one.

    Raises:
        IndexError: when ``data_array`` is empty (an empty window has no
            oldest end, so there is nothing to bound TO and one must not be
            invented), and when ``clamp`` is False and the index is out of
            range.

    THE IDENTITY GUARANTEE: for every ``0 < current_threshold <=
    len(data_array)`` both branches return the SAME element.  Bounding can
    therefore only change a run that would otherwise have aborted.
    """
    n = len(data_array)
    if clamp and n and current_threshold > n:
        return data_array[0]
    try:
        return data_array[-current_threshold]
    except IndexError:
        log.error(OUT_OF_RANGE_HINT, current_threshold, n, current_threshold,
                  " (the live window is empty)" if not n else "")
        raise


class DynamicThresholdGuard:
    """Stateful wrapper over :func:`resolve_l_end` that counts bound hits.

    One per peak-calling job.  The loops build one only on the dynamic path,
    so a run without ``--dynamic-threshold`` keeps the literal expression and
    pays nothing for a bound it cannot need.
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
        """Log the census of bound hits, or nothing at all when it never fired."""
        if self.n_clamped:
            log.warning(
                "--dynamic-threshold: the look-back index was bounded by the "
                "live read window %d time(s)%s (worst: threshold %d vs window "
                "%d). Before issue #101 this aborted the run with a bare "
                "IndexError; this run continued.",
                self.n_clamped, (" in %s" % where) if where else "",
                self.worst_threshold, self.worst_len,
            )
