#!/usr/bin/env python3
"""Score normal-sample SNP sites by a clean-genotype Bayes factor.

Score each site in a per-sample filtered normal VCF by how well its allele
counts fit *any* clean germline genotype, using a Bayes factor.

Rather than calling a genotype, this computes a model-check statistic: the
log10 Bayes factor comparing a "clean genotype" hypothesis against a diffuse
alternative.

  clean: a weighted mixture of three Beta-Binomial genotype models
           hom_ref ~ Beta(1, 99)
           het     ~ Beta(50, 50)
           hom_alt ~ Beta(99, 1)
         with mixture weights (genotype priors) w_hom_ref, w_het, w_hom_alt.
         Each shape has concentration alpha + beta = 100, so the three
         genotypes are equally confident; hom_ref/hom_alt sit at a realistic
         ~1% alt/ref error rate rather than at exactly 0 or 1.
  alt:   a diffuse Beta(1, 1) (uniform allele-fraction) model, i.e. the
         marginal likelihood of the counts when no clean genotype is assumed.

A LARGE positive value means the counts fit a clean genotype well; a NEGATIVE
value means the diffuse alternative explains the data better, i.e. the site
does not look like a clean germline genotype (e.g. mapping/allele-balance
skew).  Downstream steps threshold CLEAN_LOG10BF (e.g. < -0.5 on the Jeffreys
scale) and aggregate across normals to build the SNP blocklist.

The default genotype weights are the mean Hardy-Weinberg priors over common
gnomAD SNPs (AF >= 0.01), as computed by ``calculate_genotype_priors.py``:
  hom_ref = 0.72, het = 0.19, hom_alt = 0.09
They are hard-coded defaults but adjustable via CLI flags.

Origin: ported from an internal analysis prototype
(``bin/classify_genotypes_bayes_factor.py``) into the XOOS Copy Number Caller
SNP blocklist workflow (``copy_number_caller/scripts/snp_blocklist/``).

Usage:
    python classify_genotypes_bayes_factor.py \\
        --input  <normal_sample>.normal_dp15.vcf.gz \\
        --output <normal_sample>.normal_dp15.genotypes.vcf.gz
"""

import argparse
import math
import sys
from typing import Optional

# Beta-Binomial shape parameters for each clean genotype model.  Each shape
# has concentration alpha + beta = 100, so the three genotypes are equally
# confident:
#   hom_ref ~ Beta(1, 99)   -> mean alt fraction 0.01 (a realistic ~1% WGS
#                              error rate), NOT 0; a tighter shape such as
#                              Beta(1, 1000) demands near-zero alt reads and
#                              falsely flags clean hom_ref sites that carry a
#                              couple of noise reads.
#   het     ~ Beta(50, 50)  -> mean alt fraction 0.50
#   hom_alt ~ Beta(99, 1)   -> mean alt fraction 0.99
GENOTYPE_SHAPES = {
    "hom_ref": (1, 99),
    "het": (50, 50),
    "hom_alt": (99, 1),
}


def log_dbetabinom(k: int, n: int, alpha: float, beta: float) -> float:
    """Return the log Beta-Binomial PMF: log P(k | n, alpha, beta)."""
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
        + math.lgamma(k + alpha)
        + math.lgamma(n - k + beta)
        - math.lgamma(n + alpha + beta)
        - math.lgamma(alpha)
        - math.lgamma(beta)
        + math.lgamma(alpha + beta)
    )


def clean_genotype_log10bf(
    ref_ad: int, alt_ad: int, weights: dict[str, float]
) -> float:
    """Return the clean-genotype log10 Bayes factor for the allele counts.

    Compares the clean-genotype Beta-Binomial mixture against a diffuse
    Beta(1, 1) alternative.  ``weights`` maps each genotype name to its mixture
    weight (need not be normalized; only relative values matter for the Bayes
    factor).
    """
    n = ref_ad + alt_ad
    k = alt_ad

    # log marginal likelihood under the clean-genotype mixture, computed in a
    # numerically stable way via the log-sum-exp trick.
    log_terms = [
        math.log(weights[g]) + log_dbetabinom(k, n, alpha, beta)
        for g, (alpha, beta) in GENOTYPE_SHAPES.items()
    ]
    m = max(log_terms)
    log_clean = m + math.log(sum(math.exp(t - m) for t in log_terms))

    # log marginal likelihood under the diffuse alternative: Beta(1, 1).
    log_alt = log_dbetabinom(k, n, 1, 1)

    return (log_clean - log_alt) / math.log(10)


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Score normal-sample VCF sites by the log10 Bayes factor of a "
            "clean-genotype mixture vs a diffuse Beta(1, 1) alternative, and "
            "annotate the result as FORMAT field CLEAN_LOG10BF."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input single-sample VCF/VCF.gz (normal sample)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output bgzipped VCF.gz with CLEAN_LOG10BF FORMAT field",
    )
    parser.add_argument(
        "--weight-hom-ref",
        type=float,
        default=0.72,
        help=(
            "Mixture weight (genotype prior) for the hom_ref model "
            "[default: 0.72, mean HWE prior over common gnomAD SNPs]"
        ),
    )
    parser.add_argument(
        "--weight-het",
        type=float,
        default=0.19,
        help=(
            "Mixture weight (genotype prior) for the het model "
            "[default: 0.19, mean HWE prior over common gnomAD SNPs]"
        ),
    )
    parser.add_argument(
        "--weight-hom-alt",
        type=float,
        default=0.09,
        help=(
            "Mixture weight (genotype prior) for the hom_alt model "
            "[default: 0.09, mean HWE prior over common gnomAD SNPs]"
        ),
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    """Annotate every site of the input VCF with CLEAN_LOG10BF."""
    args = _parse_args(argv)

    weights = {
        "hom_ref": args.weight_hom_ref,
        "het": args.weight_het,
        "hom_alt": args.weight_hom_alt,
    }
    if min(weights.values()) <= 0:
        sys.exit("All genotype weights must be positive.")

    # Imported lazily so the scoring functions above can be imported and tested
    # without pysam installed (see tests/test_clean_genotype_log10bf.py).
    import pysam  # pylint: disable=import-outside-toplevel

    vcf_in = pysam.VariantFile(args.input, "r")  # pylint: disable=no-member

    # Add CLEAN_LOG10BF FORMAT field to the header before opening the output.
    vcf_in.header.formats.add(
        "CLEAN_LOG10BF",
        1,
        "Float",
        "log10 Bayes factor of a clean-genotype Beta-Binomial mixture "
        "(hom_ref/het/hom_alt) vs a diffuse Beta(1,1) alternative; negative "
        "values indicate the site does not fit a clean germline genotype",
    )

    samples = list(vcf_in.header.samples)
    if len(samples) != 1:
        sys.exit(
            f"Input VCF must contain exactly 1 sample (the normal); found {len(samples)}."
        )
    sample_name = samples[0]
    vcf_out = pysam.VariantFile(  # pylint: disable=no-member
        args.output, "wz", header=vcf_in.header
    )

    n_total = 0
    n_scored = 0

    for rec in vcf_in:
        ad = rec.samples[sample_name].get("AD")
        if ad is None or len(ad) < 2 or ad[0] is None or ad[1] is None:
            log10bf = None
        else:
            log10bf = clean_genotype_log10bf(int(ad[0]), int(ad[1]), weights)
            n_scored += 1

        rec.samples[sample_name]["CLEAN_LOG10BF"] = log10bf
        vcf_out.write(rec)
        n_total += 1

    vcf_in.close()
    vcf_out.close()

    print(
        f"Scored {n_scored} of {n_total} variants with CLEAN_LOG10BF "
        f"(weights hom_ref={weights['hom_ref']}, het={weights['het']}, "
        f"hom_alt={weights['hom_alt']})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
