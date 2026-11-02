# File formats

This page is the reference for the output and custom input formats that the Tumor Fraction Estimator produces and consumes.
Usage details and the command-line option reference for [`tumor-fraction`](usage-tumor-fraction.md) and [`contamination`](usage-contamination.md) are on their command pages.
Filter definitions referenced below are on the [Appendix](appendix.md) page.

The `tumor-fraction` and `contamination` subcommands write the same set of files.
Both produce a results file and *metrics/per_site_counts.tsv*; *metrics/read_features.tsv* is written only when `--write-read-features` is set.
The results file differs by subcommand: `tumor-fraction` writes *tumor_fraction_estimation_results.tsv* with a `tumor_fraction` metric, while `contamination` writes *contamination_estimation_results.tsv* with a `contamination_fraction` metric.
Both TSV files include a metadata header with the command-line invocation.

## Results file

Applies to: `tumor-fraction` (*tumor_fraction_estimation_results.tsv*) and `contamination` (*contamination_estimation_results.tsv*).

This file contains one row per metric in a two-column layout (`metric_name`, `value`).

| Metric name | Definition |
|:-------------|:------------|
| tumor_fraction | Estimated tumor fraction: `2 * mean_vaf`. For `contamination`, this metric is labeled `contamination_fraction`. |
| mean_vaf | Mean variant allele frequency across all passing query variant sites. |
| error_rate | Estimated sequencing error rate from the model. |
| total_alt | Total alternate allele read count across passing sites. |
| total_ref | Total reference allele read count across passing sites. |
| total_depth | Total depth across passing sites. |
| total_non_query_alts | Total count of non-reference, non-query alternate alleles across passing sites. |
| total_adjusted_alt | Alternate count adjusted for non-query-alt noise: `total_alt - total_non_query_alts / 2`. |
| mean_vaf_non_query_alts | Mean VAF of non-query alternate alleles. |
| total_variants | Total number of query variant records parsed from the VCF. |
| variants_detected | Number of query variants with at least one alternate read. |
| filtered_variants_* | Count of variants removed by each [variant-level filter](appendix.md#variant-level-filters) and certain [pileup filters](appendix.md#pileup-filters) (`LowCoverage`, `HighCoverage`, `HighObservedAltVAF`, `HighObservedNonQueryAltVAF`, `HighAltCoverage`). One row per filter type. |
| passing_variants | Number of variants that passed all variant-level filters. |
| passing_variants_detected | Number of passing variants with at least one alternate read. |
| passing_variants_with_coverage | Number of passing variants with non-zero depth. |
| filtered_reads_* | Count of reads removed by each [read-level filter](appendix.md#read-level-filters), [base-level filter](appendix.md#base-level-filters), and [pileup filter](appendix.md#pileup-filters) excluding pileup filters that are represented by a `filtered_variants_*` column. One row per filter type. |
| filtered_low_quality_alts_* | Count of alt reads removed by each [low-quality alt filter](appendix.md#low-quality-alt-filters) (`HighLowQualityAltPercentage`, `HighFilteredAltPercentage`). One row per filter type. |

{% hint style="info" %}
The reporting limit for contamination estimation is 0.0001 (0.01%). Contamination fractions below this threshold cannot be reliably measured with enough precision are reported as `<0.0001` instead of a numeric value.
{% endhint %}

When `--enable-noise-estimation` is set, additional noise metrics are appended:

| Metric name | Definition |
|:-------------|:------------|
| noise_fraction | Estimated noise fraction: `2 * noise_mean_vaf`. |
| noise_mean_vaf | Mean VAF at noise probe sites. |
| noise_total_alt | Total alternate count at noise probe sites. |
| noise_total_ref | Total reference count at noise probe sites. |
| noise_total_depth | Total depth at noise probe sites. |
| noise_total_non_query_alts | Total non-query-alt count at noise probe sites. |
| noise_total_adjusted_alt | Adjusted alternate count at noise probe sites. |
| noise_mean_vaf_non_query_alts | Mean non-query-alt VAF at noise probe sites. |
| noise_probe_sites_total | Total number of noise probe sites generated. |
| noise_probe_sites_detected | Number of noise probe sites with at least one alternate read. |
| passing_noise_probe_sites | Number of noise probe sites that passed all filters. |
| passing_noise_probe_sites_detected | Number of passing noise probe sites with at least one alternate read. |
| passing_noise_probe_sites_with_coverage | Number of passing noise probe sites with non-zero depth. |

Back to [`tumor-fraction` output](usage-tumor-fraction.md#output) and [`contamination` output](usage-contamination.md#output).

## Counts file

Applies to: `tumor-fraction` and `contamination` (*metrics/per_site_counts.tsv*).

This file contains one row per variant record (query variant and noise probe) with the following columns:

| Column | Definition |
|:--------|:------------|
| chrom | Chromosome. |
| pos | 1-based genomic position. |
| ref | Reference allele. |
| alt | Alternate allele. |
| var_type | `QueryVariant` or `NoiseProbe`. |
| ref_count | Number of reads supporting the reference allele. |
| alt_count | Number of reads supporting the alternate allele. |
| dp | Total depth at this position. |
| non_query_alt_count | Number of reads supporting alleles other than ref or the query alt. |
| seq_context | Sequence context annotation (for example, `HighConf`, `IncludedRegions`). |
| sequence | 7 bp sequence context around the variant (query variant and noise probe records). |
| filter | [Variant-level filter](appendix.md#variant-level-filters) status and certain [pileup filters](appendix.md#pileup-filters) (`LowCoverage`, `HighCoverage`, `HighObservedAltVAF`, `HighObservedNonQueryAltVAF`, `HighAltCoverage`). |
| alt_* | Per-filter type alt allele counts for each type of [read-level filter](appendix.md#read-level-filters), [base-level filter](appendix.md#base-level-filters), and [pileup filter](appendix.md#pileup-filters) excluding pileup filter types that are included in the `filter` column. |
| low_quality_alt_* | Per-filter type alt allele counts removed by each [low-quality alt filter](appendix.md#low-quality-alt-filters). |

Back to [`tumor-fraction` output](usage-tumor-fraction.md#output) and [`contamination` output](usage-contamination.md#output).

## Read features file

Applies to: `tumor-fraction` and `contamination`, only when `--write-read-features` is set.

The *metrics/read_features.tsv* file contains a per-read feature table used for development and debugging.
It is written to the *metrics/* subdirectory alongside the counts file.

Back to [`tumor-fraction` output](usage-tumor-fraction.md#output) and [`contamination` output](usage-contamination.md#output).

## Systematic noise TSV

Applies to: `tumor-fraction` and `contamination`, supplied with `--systematic-noise-tsv`.

The systematic noise file is a pre-computed list of genomic positions with known recurrent technical artifacts.
Variants matching an entry in this file are excluded from the tumor fraction calculation with the `SystematicNoise` [variant-level filter](appendix.md#variant-level-filters).

The file is tab-delimited with the following columns:

```text
chrom  pos  ref  alt
```

Rules:

- The `pos` column must use 1-based VCF coordinates.
- Each row is one artifact position defined by its chromosome, position, reference allele, and alternate allele.
- Blank lines, lines beginning with `#`, and extra trailing columns are ignored.

Accepted formats: plain-text, gzip, and bgzip. The *compression* is detected from file content, not the filename, so any readable file is accepted regardless of its extension. For genome-wide noise tracks, bgzip the file and build a tabix/CSI index (`.tbi`/`.csi`, e.g. `tabix -s1 -b2 -e2`); TFE then reads only the indexed blocks overlapping candidate variant positions instead of streaming the whole file, so memory and I/O scale with the candidate count rather than the noise-file size.

Back to [`tumor-fraction` input](usage-tumor-fraction.md#input) and [`contamination` input](usage-contamination.md#input).
