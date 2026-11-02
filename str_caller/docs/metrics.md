# Metrics

This page describes the run-level metric file that the STR Caller produces.
Running the caller and its CLI options are on the [STR caller](usage-str-caller.md) page; custom non-metric format schemas are on the [File formats](file-formats.md) page.

## Summary statistics

Applies to: `str_caller`.

*metrics.tsv* is the primary metric file containing run-level summary statistics on STR calling.
It is written to the `metrics/` subdirectory of the output directory (`<output-dir>/metrics/metrics.tsv`).
Metrics that do not apply to a given run have a value and percentage of `NA`.

A `##RocheCommandLine` metadata line is prepended to the file, recording the command used to produce it.
The first data line is the header:

```text
metric_name value percentage denominator
```

The percentage is the value divided by the value of the `denominator` metric, multiplied by 100.
When the denominator metric is `0`, the percentage is `NA`.

| Metric name | Denominator | Definition                                                        |
|:------------|:------------|:------------------------------------------------------------------|
| num_sites | NA | Number of sites in the STR catalog.                               |
| num_low_depth | num_sites | Number of sites that failed low depth quality control check.      |
| num_low_confidence | num_sites | Number of sites that failed low confidence quality control check. |

Back to [Output](usage-str-caller.md#output).
