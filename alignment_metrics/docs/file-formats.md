# File formats

This page is the reference for the tool-specific formats that Alignment Metrics consumes and produces.
For Alignment Metrics, the primary deliverable is the metrics output itself,
so this page is intentionally narrower than the Metrics page and only documents the supporting custom formats around that metrics workflow.
Usage details for [`all`](usage-all.md),
[`coverage`](usage-coverage.md),
[`accuracy`](usage-accuracy.md),
[`read`](usage-read.md),
and [`merge_and_label_bed`](usage-merge-and-label-bed.md) are on their command pages.
Metric file schemas are on the [Metrics](metrics.md) page.

## Labeled BED for `--stratify-by-label`

Applies to: `all`, `coverage`, `accuracy`, `read`, and `merge_and_label_bed`.

The labeled BED format is a standard BED interval with one additional required field.
Columns 1 to 3 are the genomic interval,
and column 4 is a comma-separated list of one or more region labels.

```text
chrom  start  end  label_a,label_b
```

Example:

```text
chr1 100 200 autosome
chr1 200 300 autosome,high_confidence
chr1 300 400 high_confidence
```

Rules:

- Every interval must have a column-4 label list when `--stratify-by-label` is enabled.
- Labels are comma-separated, with no whitespace between labels.
- At most 16 distinct labels are supported in one file.
- A label must be safe for use as a directory name, so it must not contain path separators, control characters, or whitespace.
- The output of `merge_and_label_bed` is already normalized into sorted, disjoint intervals and is the recommended way to build this file.

When used with `alignment_metrics ... --stratify-by-label`,
the Alignment Metrics tool attributes metrics independently to each label and writes one output subtree per label.

Back to [`merge_and_label_bed` output](usage-merge-and-label-bed.md#output),
[`all` input](usage-all.md#input),
[`coverage` input](usage-coverage.md#input),
[`accuracy` input](usage-accuracy.md#input),
and [`read` input](usage-read.md#input).
