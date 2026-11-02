#!/usr/bin/env python3
"""Compute per-site Hardy-Weinberg genotype priors from population AF.

Reference/analysis tool: it is NOT on the runtime SNP blocklist path.  It
derives the default genotype mixture weights used by
``classify_genotypes_bayes_factor.py``; running it over the SNP sites of
interest and averaging the ``prior_*`` columns is how those defaults
(hom_ref=0.72, het=0.19, hom_alt=0.09) were chosen.  It is kept here for
provenance and to allow re-deriving the weights against a different population
VCF (e.g. a new gnomAD version or a different AF cutoff).

For a biallelic site with population alt-allele frequency p, the
Hardy-Weinberg expected genotype frequencies are:
  hom_ref = (1 - p)^2
  het     = 2 * p * (1 - p)
  hom_alt = p^2

Origin: ported from an internal analysis prototype
(``bin/calculate_genotype_priors.py``) into the XOOS Copy Number Caller SNP
blocklist workflow (``copy_number_caller/scripts/snp_blocklist/``).

Usage:
    python calculate_genotype_priors.py \\
        --vcf af-only-gnomad.hg38.vcf.gz \\
        --output genotype_priors.tsv \\
        --min-af 0.01 \\
        --snps-only

Output columns (TSV):
    contig, pos, ref, alt, af, prior_hom_ref, prior_het, prior_hom_alt
"""

import argparse
import sys
from typing import Optional


def hwe_priors(p: float) -> tuple[float, float, float]:
    """Return (hom_ref, het, hom_alt) Hardy-Weinberg priors for alt freq p."""
    q = 1.0 - p
    return q * q, 2.0 * p * q, p * p


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-site Hardy-Weinberg genotype priors from population "
            "allele frequencies (gnomAD AF) and write them to a TSV."
        )
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Input VCF/VCF.gz carrying an INFO/AF field (e.g. gnomAD)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV of per-site genotype priors",
    )
    parser.add_argument(
        "--min-af",
        type=float,
        default=0.0,
        help="Only emit sites with alt AF >= this value [default: 0.0]",
    )
    parser.add_argument(
        "--snps-only",
        action="store_true",
        help="Restrict to biallelic SNP sites (single-base ref and alt)",
    )
    return parser.parse_args(argv)


def _format_row(rec, snps_only: bool, min_af: float) -> Optional[str]:
    """Return the TSV row for a VCF record, or None if it should be skipped.

    Skips non-biallelic sites, indels (when ``snps_only``), sites without an
    ``INFO/AF``, and sites whose alt frequency is out of ``(0, 1]`` or below
    ``min_af``.
    """
    # Only biallelic sites have a single, unambiguous alt allele.
    if len(rec.alts or ()) != 1:
        return None
    alt = rec.alts[0]
    ref = rec.ref

    if snps_only and (len(ref) != 1 or len(alt) != 1):
        return None

    af = rec.info.get("AF")
    # INFO/AF is Number=A; for a biallelic record it is a 1-tuple.
    if isinstance(af, tuple):
        af = af[0] if af else None
    if af is None:
        return None

    p = float(af)
    if p < min_af or p <= 0.0 or p > 1.0:
        return None

    hom_ref, het, hom_alt = hwe_priors(p)
    return (
        f"{rec.contig}\t{rec.pos}\t{ref}\t{alt}\t{p:.6g}\t"
        f"{hom_ref:.6g}\t{het:.6g}\t{hom_alt:.6g}\n"
    )


def main(argv: Optional[list[str]] = None) -> None:
    """Write per-site Hardy-Weinberg genotype priors to a TSV."""
    args = _parse_args(argv)

    # Imported lazily so hwe_priors can be imported and tested without pysam.
    import pysam  # pylint: disable=import-outside-toplevel

    vcf_in = pysam.VariantFile(args.vcf, "r")  # pylint: disable=no-member

    n_total = 0
    n_written = 0

    with open(args.output, "w", encoding="utf-8") as out:
        header = "\t".join(
            [
                "contig",
                "pos",
                "ref",
                "alt",
                "af",
                "prior_hom_ref",
                "prior_het",
                "prior_hom_alt",
            ]
        )
        out.write(header + "\n")

        for rec in vcf_in:
            n_total += 1
            row = _format_row(rec, args.snps_only, args.min_af)
            if row is not None:
                out.write(row)
                n_written += 1

    vcf_in.close()

    print(
        f"Read {n_total} records; wrote {n_written} per-site priors to {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
