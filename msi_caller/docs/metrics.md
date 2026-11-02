# Metrics

This page describes the run-level quality-control metrics that the MSI Caller produces in *msi_statistics.tsv*.
The JSON report and per-site VCF calls are described on the [File formats](file-formats.md) page.

## TSV metrics

Applies to: `tumor-only`, `tumor-normal`.

The TSV begins with a `##RocheCommandLine=<...>` provenance line.
When the tumor BAM contains a read group, the following `##Sample=<...>` line records nonempty `SM`, `LB`, `PL`, and `PM` values from its first `@RG` entry as `Name`, `LibraryId`, `Platform`, and `PlatformModel`.
The `##Sample` line is omitted when the BAM has no read group.

The TSV output contains `metric` and `value` columns with the following rows:

| Metric name | Definition |
|:------------|:-----------|
| mode | Analysis mode: `TUMOR_ONLY` or `TUMOR_NORMAL`. |
| total_microsatellite_sites | Number of evaluated sites. A nonnegative integer. |
| total_tumor_usable_sites | Number of usable tumor sites. In tumor-normal mode, a site is counted only when both tumor and matched normal pass the filters. A nonnegative integer. |
| total_normal_usable_sites | Number of usable matched-normal sites in tumor-normal mode. This is the same paired-site count as `total_tumor_usable_sites`. A nonnegative integer, or `NA` in tumor-only mode. |

Back to [`tumor-only` output](usage-tumor-only.md#output) or [`tumor-normal` output](usage-tumor-normal.md#output).
