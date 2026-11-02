"""Behavior-preserving self-check for the Bayes-factor scoring.

These reference values were computed from the reference analysis prototype
(``bin/classify_genotypes_bayes_factor.py``) with the default weights
(hom_ref=0.72, het=0.19, hom_alt=0.09). Any change to the math that alters the
scores will fail these assertions.

The module under test imports pysam lazily inside ``main``, so importing the
scoring functions here does not require pysam to be installed.

Run:
    python -m pytest test_clean_genotype_log10bf.py
    # or, without pytest:
    python test_clean_genotype_log10bf.py
"""

import math
import os
import sys

# The module under test lives one directory up; add it to the path before the
# import so this self-check can run in place without installation.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# pylint: disable=wrong-import-position
from classify_genotypes_bayes_factor import (  # noqa: E402
    clean_genotype_log10bf,
    log_dbetabinom,
)

_WEIGHTS = {"hom_ref": 0.72, "het": 0.19, "hom_alt": 0.09}

# Default "skewed" threshold on the Jeffreys scale, matching the workflow
# default bf_threshold.
_BF_THRESHOLD = -0.5

# (ref_ad, alt_ad) -> expected CLEAN_LOG10BF with the default weights.
_REFERENCE = {
    (30, 0): 1.2337396787,
    (15, 15): -0.1273459431,
    (0, 30): 0.3306497205,
    (28, 2): -0.0377409041,
    (20, 10): -0.6717329438,
    (25, 5): -1.9039121919,
}


def test_reference_scores() -> None:
    """Each reference count pair reproduces the source CLEAN_LOG10BF value."""
    for (ref_ad, alt_ad), expected in _REFERENCE.items():
        got = clean_genotype_log10bf(ref_ad, alt_ad, _WEIGHTS)
        assert math.isclose(got, expected, rel_tol=0, abs_tol=1e-9), (
            f"ref={ref_ad} alt={alt_ad}: expected {expected}, got {got}"
        )


def test_clean_hom_ref_scores_positive() -> None:
    """A clean hom_ref site (all ref reads) fits a clean genotype well."""
    assert clean_genotype_log10bf(30, 0, _WEIGHTS) > 0


def test_intermediate_af_scores_negative() -> None:
    """A ~1/3 alt fraction fits no clean genotype and is flagged skewed."""
    assert clean_genotype_log10bf(20, 10, _WEIGHTS) < _BF_THRESHOLD


def test_log_dbetabinom_matches_uniform() -> None:
    """Beta(1, 1) Beta-Binomial reduces to a discrete uniform over k.

    P(k | n) = 1 / (n + 1) for all k in [0, n].
    """
    n = 10
    for k in range(n + 1):
        assert math.isclose(
            math.exp(log_dbetabinom(k, n, 1, 1)), 1.0 / (n + 1), abs_tol=1e-12
        )


if __name__ == "__main__":
    test_reference_scores()
    test_clean_hom_ref_scores_positive()
    test_intermediate_af_scores_negative()
    test_log_dbetabinom_matches_uniform()
    print("All self-checks passed.")
