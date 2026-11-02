# Metrics

This page is the reference for the summary metrics that Copy Number Caller writes to *metrics/metrics.tsv*.
The `germline-wgs`, `tumor-normal-wgs`, and `tumor-only-wgs` commands emit this file alongside their callsets; the metrics reported depend on whether the command runs a germline or a somatic workflow.
The `predict-somatic-cnv` command does not write a metrics file.

Each section lists the commands it **applies to** (linked to their usage pages) and, at the end, links back to the command output where the file is produced.

*metrics/metrics.tsv* is a tab-delimited file with a two-column header:

```text
Metric  Value
```

## Coverage-check metrics

Applies to: [`germline-wgs`](usage-germline-wgs.md), [`tumor-normal-wgs`](usage-tumor-normal-wgs.md), and
[`tumor-only-wgs`](usage-tumor-only-wgs.md).

These metrics are written only when a workflow exits early because coverage is below the required threshold. They are
not present in metrics files from successful runs.

| Metric name | Definition |
|:---------|:-----------------------|
| Median coverage | Median coverage across the analysis windows. A non-negative value, formatted to two decimal places. |
| Percentage of low coverage windows | Percentage of analysis windows with coverage below the configured threshold. A value in [0, 100], formatted to two decimal places. |

Back to [`germline-wgs` output](usage-germline-wgs.md#output), [`tumor-normal-wgs` output](usage-tumor-normal-wgs.md#output), and
[`tumor-only-wgs` output](usage-tumor-only-wgs.md#output).

## Germline metrics

Applies to: [`germline-wgs`](usage-germline-wgs.md).

| Metric name | Definition |
|:---------|:-----------------------|
| Number of targets   | Total number of targets used in the analysis. An integer >= 0. |
| Sex   | Predicted or specified sex of the sample. Possible values are male, female, or unknown. |
| Number of segments | Number of predicted segments. An integer >= 0. |
| Number of gain segments   | Number of predicted gain segments. An integer >= 0. |
| Number of loss segments   | Number of predicted loss segments. An integer >= 0. |
| Number of passing gain segments  | Number of predicted gain segments that passed filtering. An integer >= 0. |
| Number of passing loss segments  | Number of predicted loss segments that passed filtering. An integer >= 0. |

Back to [`germline-wgs` output](usage-germline-wgs.md#output).

## Somatic metrics

Applies to: [`tumor-normal-wgs`](usage-tumor-normal-wgs.md), [`tumor-only-wgs`](usage-tumor-only-wgs.md).

| Metric name | Definition |
|:---------|:-----------------------|
| Number of targets   | Total number of targets used in the analysis. An integer >= 0. |
| Number of SNPs | Total number of SNPs used in the analysis. An integer >= 0. |
| Sex   | Predicted or specified sex of the sample. Possible values are male, female, or unknown. |
| Tumor Purity | Predicted tumor purity as a fraction in [0, 1], formatted to two decimal places, or `NA` when unavailable. |
| Tumor Ploidy | Predicted tumor ploidy. Estimates may be fractional (the search grid uses 0.1 increments); values are formatted to two decimal places, or `NA` when unavailable. |
| Number of segments | Number of predicted segments. An integer >= 0. |
| Number of somatic gain (non-LOH) segments   | Number of predicted somatic gain (non-LOH) segments. An integer >= 0. |
| Number of somatic loss segments   | Number of predicted somatic loss segments. An integer >= 0. |
| Number of somatic copy-neutral LOH segments   | Number of predicted somatic copy-neutral LOH segments. An integer >= 0. |
| Number of somatic gain LOH segments   | Number of predicted gain LOH segments. An integer >= 0. |
| Number of passing somatic gain (non-LOH) segments   | Number of predicted somatic gain (non-LOH) segments that passed filtering. An integer >= 0. |
| Number of passing somatic loss segments   | Number of predicted somatic loss segments that passed filtering. An integer >= 0. |
| Number of passing somatic copy-neutral LOH segments   | Number of predicted somatic copy-neutral LOH segments that passed filtering. An integer >= 0. |
| Number of passing somatic gain LOH segments   | Number of predicted gain LOH segments that passed filtering. An integer >= 0. |

### Tumor-normal only metrics

The following metrics are present only in [`tumor-normal-wgs`](usage-tumor-normal-wgs.md) output; they are not written by [`tumor-only-wgs`](usage-tumor-only-wgs.md), which has no matched normal VCF.

These rows are written only when the tumor-normal workflow exits early with a VCF check: when no variants remain after
filtering, or when a somatic-flagged VCF has too low a heterozygous germline SNP fraction. They are not included in
successful calling runs or low-coverage early-exit metrics.

| Metric name | Definition |
|:---------|:-----------------------|
| Number of heterozygous germline SNPs | Number of biallelic SNP records retained as heterozygous germline variants after the normal-sample depth and BAF filters and requiring nonzero tumor depth. An integer >= 0. |
| Number of somatic variants | Number of biallelic SNP records marked `PASS` in the input VCF; these are counted regardless of depth. An integer >= 0. |
| Het SNP fraction | Number of heterozygous germline SNPs divided by the sum of heterozygous germline SNPs and somatic variants. A value in [0, 1], or `NA` when the denominator is zero. |

Back to [`tumor-normal-wgs` output](usage-tumor-normal-wgs.md#output) and
[`tumor-only-wgs` output](usage-tumor-only-wgs.md#output).
