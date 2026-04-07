from scipy.stats import poisson
from typing import Optional
import numpy as np

def poisson_pvalue(observed: int, expected_lambda: float) -> float:
    """Compute Poisson p-value for peak significance.

    Tests whether the observed coverage is significantly higher than
    expected under a Poisson model with rate = expected_lambda.

    Args:
        observed: Observed read count (peak height)
        expected_lambda: Expected background rate (local lambda)

    Returns:
        p-value: P(X >= observed | Poisson(expected_lambda))
    """
    if expected_lambda <= 0:
        expected_lambda = 1.0
    if observed <= 0:
        return 1.0

    return float(poisson.sf(observed - 1, expected_lambda))


def nb_pvalue(observed: int, mu: float, dispersion: float) -> float:
    """Compute Negative Binomial p-value for peak significance.

    NB is more appropriate than Poisson for scRNA-seq data where
    variance > mean (overdispersion).

    Args:
        observed: Observed read count
        mu: Expected mean (from local lambda)
        dispersion: Dispersion parameter r (smaller = more overdispersed)

    Returns:
        p-value: P(X >= observed | NB(mu, dispersion))
    """
    from scipy.stats import nbinom

    if mu <= 0 or dispersion <= 0:
        return 1.0
    if observed <= 0:
        return 1.0

    # scipy parameterization: n=dispersion, p=dispersion/(dispersion+mu)
    n = dispersion
    p = dispersion / (dispersion + mu)

    return float(nbinom.sf(observed - 1, n, p))


def estimate_dispersion(counts: list) -> float:
    """Estimate NB dispersion parameter from a list of per-cell counts.

    Method of moments: r = mean^2 / (variance - mean)
    If variance <= mean (no overdispersion), returns a large r (Poisson-like).

    Args:
        counts: List of per-cell read counts for a PAS

    Returns:
        Estimated dispersion parameter r
    """
    if len(counts) < 2:
        return 100.0  # default: near-Poisson

    mean_val = np.mean(counts)
    var_val = np.var(counts, ddof=1)

    if var_val <= mean_val or mean_val <= 0:
        return 100.0  # no overdispersion detected

    r = mean_val ** 2 / (var_val - mean_val)
    return max(r, 0.01)  # minimum to avoid numerical issues
