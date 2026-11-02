# Metrics

This page is the reference for the metric files that Alignment Metrics produces.
For Alignment Metrics, those metric files are the primary deliverable,
unlike many other modules where metrics are only supporting outputs.
Usage details for [`all`](usage-all.md),
[`coverage`](usage-coverage.md),
[`accuracy`](usage-accuracy.md),
[`gcbias`](usage-gcbias.md),
and [`read`](usage-read.md) are on their command pages.
Custom non-metric format schemas are on the [File formats](file-formats.md) page.

## Coverage metrics

Applies to: `coverage` and `all`.

### `coverage_histograms.tsv`

This file contains a histogram of coverage depths.
The homopolymer-only variant `coverage_histograms_hp.tsv` uses the same schema and is emitted when `--enable-homopolymer-metrics` is enabled.

All dataset types use these columns:

```text
coverage    position_count    post_filter_position_count
```

`post_filter_position_count` counts positions that pass the configured base-call quality and base-type filters.

| Metric name | Definition |
|:-------------|:------------|
| `coverage` | Number of reads covering a position. |
| `position_count` | Number of positions observed at the given coverage depth. |
| `post_filter_position_count` | Number of positions at the given depth that pass the configured base-level filters. |

### `coverage_stats.tsv`

This file contains coverage summary statistics.
The homopolymer-only variant `coverage_stats_hp.tsv` uses the same schema and is emitted when `--enable-homopolymer-metrics` is enabled.

All dataset types use these columns:

```text
metric_name    any_coverage    post_filter_coverage
```

| Metric name | Definition |
|:-------------|:------------|
| `median` | Median coverage for the corresponding column. |
| `mean` | Mean coverage for the corresponding column. |
| `min` | Minimum coverage for the corresponding column. |
| `max` | Maximum coverage for the corresponding column. |
| `percentile_X` | Coverage at the Xth percentile for the corresponding column. |
| `ratio_X_to_Y_percentile` | Ratio of the Xth percentile to the Yth percentile for the corresponding column. |

Notes:

- `post_filter_coverage` is coverage after base-call quality and base-type filtering.
- If no BED file is provided, zero-coverage positions are excluded from the mean, median, percentile, and ratio calculations.
- If every position has zero coverage, the summary reports `NA` when no BED file is provided, or `0` when a BED file is provided.

### `coverage_distribution_summary.tsv`

This file summarizes counts and percentages of positions that meet coverage thresholds.
The homopolymer-only variant `coverage_distribution_summary_hp.tsv` uses the same schema and is emitted when `--enable-homopolymer-metrics` is enabled.

All dataset types use these columns:

```text
metric_name    any_coverage    post_filter_coverage
```

| Metric name | Definition | Denominator |
|:-------------|:------------|:-------------|
| `total_positions` | Total number of reference positions considered. | N/A |
| `positions_with_no_coverage` | Total number of positions with 0 coverage in the corresponding column. | N/A |
| `percentage_of_positions_with_no_coverage` | Percentage of positions with 0 coverage in the corresponding column. | `total_positions` |
| `positions_with_coverage` | Total number of positions with at least 1x coverage in the corresponding column. | N/A |
| `percentage_of_positions_with_coverage` | Percentage of positions with at least 1x coverage in the corresponding column. | `total_positions` |
| `positions_with_at_least_Xx_coverage` | Total number of positions with at least X coverage in the corresponding column. | N/A |
| `percentage_of_positions_with_at_least_Xx_coverage` | Percentage of positions with at least X coverage in the corresponding column. | `total_positions` |

Constraint:

- `total_positions = positions_with_no_coverage + positions_with_coverage`

### `mean_coverage_histogram.tsv`

This file contains a histogram of mean coverage values across target regions.

```text
mean_coverage    target_region_count
```

| Metric name | Definition |
|:-------------|:------------|
| `mean_coverage` | Mean coverage across a target region. |
| `target_region_count` | Number of target regions with the given mean coverage. |

### `coverage_uniformity_summary.tsv`

This file contains targeted-enrichment coverage uniformity summaries.

```text
metric_name    value
```

| Metric name | Definition |
|:-------------|:------------|
| `fold_80_base_penalty` | Mean target coverage divided by the 20th percentile coverage. Reports `NA` if the 20th percentile coverage is zero. |
| `pct_bases_at_0.5x_to_2x_mean_coverage` | Percentage of targeted positions observed within 50% to 200% of the mean target coverage. |
| `mean_coverage` | Mean coverage across all target regions, including positions with no coverage. |
| `median_coverage` | Median coverage across all target regions, including positions with no coverage. |
| `target_region_count` | Total number of target regions specified in the BED file. |
| `target_regions_with_no_coverage` | Number of target regions with zero mean coverage. |

Back to [`coverage` output](usage-coverage.md#output) and [`all` output](usage-all.md#output).

## Accuracy metrics

Applies to: `accuracy` and `all`.

### `base_level_accuracy_summary.tsv`

This file reports counts, percentages, and Phred-scaled error scores for substitutions and indels.

```text
type    count    denominator    phred    percentage
```

`type` is one of `total_aligned_bases`, `substitutions`, `insertion_events`, `inserted_bases`, `deletion_events`, `deleted_bases`, `indel_events`, `indel_bases`, `all_errors_with_indel_events`, or `all_errors_with_indel_bases`.

| Metric name | Definition |
|:-------------|:------------|
| `total_aligned_bases` | Total number of aligned bases. |
| `substitutions` | Total number of substitution errors. |
| `insertion_events` | Total number of insertion events, counting each inserted stretch as one event. |
| `inserted_bases` | Total number of inserted bases. |
| `deletion_events` | Total number of deletion events, counting each deleted stretch as one event. |
| `deleted_bases` | Total number of deleted bases. |
| `indel_events` | Sum of insertion events and deletion events. |
| `indel_bases` | Sum of inserted bases and deleted bases. |
| `all_errors_with_indel_events` | Sum of substitutions and indel events. |
| `all_errors_with_indel_bases` | Sum of substitutions and indel bases. |
| `total_positions` | Total number of reference positions considered. |
| `positions_skipped_due_to_low_depth` | Number of positions skipped because depth is below `--min-depth`. |
| `positions_skipped_due_to_high_alt_allele_fraction` | Number of positions skipped because the alt allele fraction exceeds `--max-alt-allele-fraction`. |

Notes:

- `substitutions` appears twice with different denominators because the Alignment Metrics tool reports substitution rates once against base-plus-indel-events and once against base-plus-indel-bases.
- The Phred-scaled error score is `min(93, -10 * log10(errors / (total bases + indels)))`.

Constraints:

- `all_errors_with_indel_events = substitutions + indel_events`
- `all_errors_with_indel_bases = substitutions + indel_bases`
- `indel_events = insertion_events + deletion_events`
- `indel_bases = inserted_bases + deleted_bases`

### `errors_by_substitution_type.tsv`

This file breaks substitutions down by substitution class and, when available, by strand, read type, or cluster type.

Possible column sets include:

```text
substitution_type    substitutions_total    substitutions_forward_strand    substitutions_reverse_strand
```

```text
substitution_type    substitutions_total    substitutions_full_read    substitutions_partial_read    substitutions_forward_strand    substitutions_reverse_strand
```

```text
substitution_type    substitutions_total    substitutions_mixed_strand_cluster    substitutions_forward_cluster    substitutions_reverse_cluster    substitutions_full_cluster    substitutions_mixed_full_and_partial_cluster    substitutions_partial_cluster    substitutions_forward_strand    substitutions_reverse_strand
```

The rows are the 12 substitution classes from `A->C` through `T->G`.

| Metric name | Applicable scenarios | Definition |
|:-------------|:----------------------|:------------|
| `substitution_type` | All | Substitution class. |
| `substitutions_total` | All | Number of substitution errors of the given type. |
| `substitutions_forward_strand` | All | Number of substitution errors of the given type on forward reads. |
| `substitutions_reverse_strand` | All | Number of substitution errors of the given type on reverse reads. |
| `substitutions_full_read` | Datasets with read-type information | Number of substitution errors of the given type in full reads. |
| `substitutions_partial_read` | Datasets with read-type information | Number of substitution errors of the given type in partial reads. |
| `substitutions_mixed_strand_cluster` | Consensus datasets with cluster metadata | Number of substitution errors of the given type in mixed-strand cluster reads. |
| `substitutions_forward_cluster` | Consensus datasets with cluster metadata | Number of substitution errors of the given type in forward-cluster reads. |
| `substitutions_reverse_cluster` | Consensus datasets with cluster metadata | Number of substitution errors of the given type in reverse-cluster reads. |
| `substitutions_full_cluster` | Consensus datasets with cluster metadata | Number of substitution errors of the given type in full-cluster reads. |
| `substitutions_mixed_full_and_partial_cluster` | Consensus datasets with cluster metadata | Number of substitution errors of the given type in mixed full-and-partial cluster reads. |
| `substitutions_partial_cluster` | Consensus datasets with cluster metadata | Number of substitution errors of the given type in partial-cluster reads. |

Constraints for each substitution type:

- `substitutions_total = substitutions_forward_strand + substitutions_reverse_strand`
- `substitutions_total = substitutions_full_read + substitutions_partial_read` when read-type information is available
- `substitutions_total = substitutions_mixed_strand_cluster + substitutions_forward_cluster + substitutions_reverse_cluster` when cluster information is available
- `substitutions_total = substitutions_full_cluster + substitutions_mixed_full_and_partial_cluster + substitutions_partial_cluster` when cluster information is available

### `errors_by_read_type.tsv`

This file reports error counts by strand and read type.

```text
type    count    denominator    percentage
```

| Metric name | Definition | Denominator for percentage |
|:-------------|:------------|:----------------------------|
| `substitutions_total` | Total number of substitution errors. | N/A |
| `substitutions_forward_strand` | Number of substitution errors from forward reads. | `substitutions_total` |
| `substitutions_reverse_strand` | Number of substitution errors from reverse reads. | `substitutions_total` |
| `substitutions_full_read` | Number of substitution errors from full reads. | `substitutions_total` |
| `substitutions_partial_read` | Number of substitution errors from partial reads. | `substitutions_total` |
| `insertions_total` | Total number of insertion errors. | N/A |
| `insertions_forward_strand` | Number of insertion errors from forward reads. | `insertions_total` |
| `insertions_reverse_strand` | Number of insertion errors from reverse reads. | `insertions_total` |
| `insertions_full_read` | Number of insertion errors from full reads. | `insertions_total` |
| `insertions_partial_read` | Number of insertion errors from partial reads. | `insertions_total` |
| `deletions_total` | Total number of deletion errors. | N/A |
| `deletions_forward_strand` | Number of deletion errors from forward reads. | `deletions_total` |
| `deletions_reverse_strand` | Number of deletion errors from reverse reads. | `deletions_total` |
| `deletions_full_read` | Number of deletion errors from full reads. | `deletions_total` |
| `deletions_partial_read` | Number of deletion errors from partial reads. | `deletions_total` |

### `errors_by_cluster_size.tsv`

This file reports error counts by consensus-cluster size and cluster composition.

| Metric name | Definition |
|:-------------|:------------|
| `cluster_size` | Number of reads in the cluster used to generate the consensus read. |
| `total_bases` | Number of bases in consensus reads of the given cluster size. |
| `total_bases_mixed_strand_cluster` | Number of bases in mixed-strand cluster reads of the given cluster size. |
| `total_bases_forward_cluster` | Number of bases in forward-cluster reads of the given cluster size. |
| `total_bases_reverse_cluster` | Number of bases in reverse-cluster reads of the given cluster size. |
| `total_bases_full_cluster` | Number of bases in full-cluster reads of the given cluster size. |
| `total_bases_mixed_full_and_partial_cluster` | Number of bases in mixed full-and-partial cluster reads of the given cluster size. |
| `total_bases_partial_cluster` | Number of bases in partial-cluster reads of the given cluster size. |
| `overall_phred` | Overall Phred-scaled error score for consensus reads of the given cluster size. |
| `substitutions_total` | Number of substitution errors in consensus reads of the given cluster size. |
| `substitution_phred` | Phred-scaled error score for substitution errors in consensus reads of the given cluster size. |
| `substitutions_mixed_strand_cluster` | Number of substitution errors in mixed-strand cluster reads of the given cluster size. |
| `substitutions_forward_cluster` | Number of substitution errors in forward-cluster reads of the given cluster size. |
| `substitutions_reverse_cluster` | Number of substitution errors in reverse-cluster reads of the given cluster size. |
| `substitutions_full_cluster` | Number of substitution errors in full-cluster reads of the given cluster size. |
| `substitutions_mixed_full_and_partial_cluster` | Number of substitution errors in mixed full-and-partial cluster reads of the given cluster size. |
| `substitutions_partial_cluster` | Number of substitution errors in partial-cluster reads of the given cluster size. |
| `insertions_total` | Number of insertion errors in consensus reads of the given cluster size. |
| `insertion_phred` | Phred-scaled error score for insertion errors in consensus reads of the given cluster size. |
| `insertions_mixed_strand_cluster` | Number of insertion errors in mixed-strand cluster reads of the given cluster size. |
| `insertions_forward_cluster` | Number of insertion errors in forward-cluster reads of the given cluster size. |
| `insertions_reverse_cluster` | Number of insertion errors in reverse-cluster reads of the given cluster size. |
| `insertions_full_cluster` | Number of insertion errors in full-cluster reads of the given cluster size. |
| `insertions_mixed_full_and_partial_cluster` | Number of insertion errors in mixed full-and-partial cluster reads of the given cluster size. |
| `insertions_partial_cluster` | Number of insertion errors in partial-cluster reads of the given cluster size. |
| `deletions_total` | Number of deletion errors in consensus reads of the given cluster size. |
| `deletion_phred` | Phred-scaled error score for deletion errors in consensus reads of the given cluster size. |
| `deletions_mixed_strand_cluster` | Number of deletion errors in mixed-strand cluster reads of the given cluster size. |
| `deletions_forward_cluster` | Number of deletion errors in forward-cluster reads of the given cluster size. |
| `deletions_reverse_cluster` | Number of deletion errors in reverse-cluster reads of the given cluster size. |
| `deletions_full_cluster` | Number of deletion errors in full-cluster reads of the given cluster size. |
| `deletions_mixed_full_and_partial_cluster` | Number of deletion errors in mixed full-and-partial cluster reads of the given cluster size. |
| `deletions_partial_cluster` | Number of deletion errors in partial-cluster reads of the given cluster size. |

### `hp_errors.tsv`

This file summarizes homopolymer support and homopolymer error rates.

| Metric name | Definition | Denominator |
|:-------------|:------------|:-------------|
| `hp_base` | Homopolymer base, or `any` for all homopolymers of a given length. | N/A |
| `hp_length` | Homopolymer length. | N/A |
| `hp_count` | Number of homopolymers of the given type in the reference. | N/A |
| `total_reads` | Number of reads passing filters that cover the first base of the homopolymer. | N/A |
| `spanning_reads` | Number of reads passing filters that span the homopolymer with flanking context. | N/A |
| `percentage_spanning` | Percentage of reads that span the homopolymer. | `total_reads` |
| `mean_spanning_coverage` | `spanning_reads / hp_count`. | N/A |
| `discordant_reads` | Number of spanning reads with at least one discordant duplex base or zero-quality base in the homopolymer. | N/A |
| `percentage_discordant` | Percentage of discordant reads. | `total_reads` |
| `low_quality_reads` | Number of spanning reads with at least one base that fails the base-call quality filter. | N/A |
| `percentage_low_quality` | Percentage of low-quality reads. | `total_reads` |
| `effective_reads` | Number of spanning reads whose homopolymer-aligned bases pass the quality and base-type filters and whose insertions are homogeneous unless heterogeneous insertions are allowed. | N/A |
| `percentage_effective` | Percentage of effective reads. | `total_reads` |
| `mean_effective_coverage` | `effective_reads / hp_count`. | N/A |
| `effective_reads_with_insertion` | Number of effective reads with at least one insertion in the homopolymer. | N/A |
| `effective_reads_with_deletion` | Number of effective reads with at least one deletion in the homopolymer. | N/A |
| `effective_reads_with_substitution` | Number of effective reads with at least one substitution in the homopolymer. | N/A |
| `insertion_error_rate` | Percentage of effective reads with at least one insertion in the homopolymer. | `effective_reads` |
| `deletion_error_rate` | Percentage of effective reads with at least one deletion in the homopolymer. | `effective_reads` |
| `phred` | Phred-scale homopolymer accuracy considering insertions and deletions. | `effective_reads` |

### `qscore_stats.tsv`

This file compares predicted base-call quality scores with empirical error rates.

| Metric name | Definition |
|:-------------|:------------|
| `q_score` | Predicted base-call quality score. |
| `matches` | Number of bases with that predicted quality that match the reference. |
| `mismatches` | Number of bases with that predicted quality that do not match the reference. |
| `empirical_phred_score` | Empirical Phred-scaled error score derived from the observed matches and mismatches. |
| `forward_cluster_matches` | Matches from forward-cluster consensus reads at that predicted quality. |
| `forward_cluster_mismatches` | Mismatches from forward-cluster consensus reads at that predicted quality. |
| `forward_cluster_empirical_phred_score` | Empirical Phred-scaled error score for forward-cluster consensus reads. |
| `reverse_cluster_matches` | Matches from reverse-cluster consensus reads at that predicted quality. |
| `reverse_cluster_mismatches` | Mismatches from reverse-cluster consensus reads at that predicted quality. |
| `reverse_cluster_empirical_phred_score` | Empirical Phred-scaled error score for reverse-cluster consensus reads. |
| `mixed_strand_cluster_matches` | Matches from mixed-strand cluster reads at that predicted quality. |
| `mixed_strand_cluster_mismatches` | Mismatches from mixed-strand cluster reads at that predicted quality. |
| `mixed_strand_cluster_empirical_phred_score` | Empirical Phred-scaled error score for mixed-strand cluster reads. |

Back to [`accuracy` output](usage-accuracy.md#output) and [`all` output](usage-all.md#output).

## Read metrics

Applies to: `read` and `all`.

### `read_metrics_summary.tsv`

This file summarizes read counts, alignment counts, and base-count summaries.

```text
metric_name    count    denominator    percentage
```

| Metric name | Definition | Denominator for percentage |
|:-------------|:------------|:----------------------------|
| `alignments_plus_unmapped_reads` | Total number of alignment records plus unmapped reads. | N/A |
| `alignment_records` | Number of primary, secondary, and supplementary alignment records. | `alignments_plus_unmapped_reads` |
| `secondary_alignments` | Number of secondary alignment records. | `alignments_plus_unmapped_reads` |
| `supplementary_alignments` | Number of supplementary alignment records. | `alignments_plus_unmapped_reads` |
| `total_reads` | Total number of reads. | N/A |
| `unmapped_reads` | Total number of unmapped reads, including unplaced unmapped reads. | `total_reads` |
| `mapped_reads` | Total number of mapped reads. | `total_reads` |
| `forward_reads` | Number of forward mapped reads. | `mapped_reads` |
| `reverse_reads` | Number of reverse mapped reads. | `mapped_reads` |
| `duplicate_reads` | Number of duplicate reads. | `mapped_reads` |
| `mapq_zero_reads` | Number of mapped reads with mapping quality 0. | `mapped_reads` |
| `reads_passing_filter` | Number of reads that pass the mapping-quality and flag filters. | `total_reads` |
| `forward_reads_passing_filter` | Number of forward mapped reads that pass the filters. | `reads_passing_filter` |
| `reverse_reads_passing_filter` | Number of reverse mapped reads that pass the filters. | `reads_passing_filter` |
| `total_bases` | Total number of bases in all reads, including soft clips, inserted bases, and bases from unmapped reads. | N/A |
| `unmapped_bases` | Number of bases in unmapped reads. | `total_bases` |
| `aligned_bases` | Number of aligned bases in the primary and supplementary alignments of each read. | `total_bases` |
| `unique_aligned_bases` | Number of distinct primary read bases that are aligned at least once. | `total_bases` |
| `soft_clipped_bases` | Number of soft-clipped bases in the primary alignment of mapped reads that pass the filters. | `total_bases` |
| `unmapped_soft_clipped_bases` | Number of soft-clipped bases not mapped by any supplementary alignment. | `total_bases` |
| `gc_bases` | Number of G or C bases in all reads. | `total_bases` |
| `total_bases_passing_filter` | Number of bases in reads that pass the mapping-quality and flag filters. | `total_bases` |
| `aligned_bases_passing_filter` | Number of aligned bases in reads that pass the filters. | `total_bases_passing_filter` |
| `unique_aligned_bases_passing_filter` | Number of distinct aligned primary read bases in reads that pass the filters. | `total_bases_passing_filter` |
| `soft_clipped_bases_passing_filter` | Number of soft-clipped bases in reads that pass the filters. | `total_bases_passing_filter` |
| `unmapped_soft_clipped_bases_passing_filter` | Number of soft-clipped bases not mapped by any supplementary alignment in reads that pass the filters. | `total_bases_passing_filter` |
| `gc_bases_passing_filter` | Number of G or C bases in reads that pass the filters. | `total_bases_passing_filter` |

When read-type information is present, the file also includes `full_length_reads`, `partial_length_reads`, `full_length_reads_passing_filter`, and `partial_length_reads_passing_filter`.
When consensus-cluster metadata is present, the file also includes cluster-composition read counts and post-filter cluster-composition read counts.

Constraints:

- `alignments_plus_unmapped_reads = total_reads + supplementary_alignments + secondary_alignments`
- `alignment_records = alignments_plus_unmapped_reads - unmapped_reads`
- `total_reads = mapped_reads + unmapped_reads`
- `mapped_reads = forward_reads + reverse_reads`

The file aligns with `samtools flagstat` for several key counts:

| `read_metrics_summary.tsv` metric | `samtools flagstat` metric |
|:-----------------------------------|:----------------------------|
| `total_reads` | `primary` |
| `mapped_reads` | `primary mapped` |
| `secondary_alignments` | `secondary` |
| `supplementary_alignments` | `supplementary` |
| `duplicate_reads` | `duplicates` |

### `read_metrics_te.tsv`

This file contains targeted-enrichment read summaries.

```text
metric_name    count    denominator    percentage
```

| Metric name | Definition | Denominator for percentage |
|:-------------|:------------|:----------------------------|
| `total_mapped_alignments_in_bam` | Total number of mapped alignments in the BAM, including primary, secondary, and supplementary alignments. | N/A |
| `on_target_alignments` | Number of alignments that intersect at least one BED region. | `total_mapped_alignments_in_bam` |
| `on_target_reads_passing_filter` | Number of reads that intersect at least one BED region and pass the read filters. | `total_mapped_alignments_in_bam` |

Note:

- `on_target_alignments` equals `mapped_reads + secondary_alignments + supplementary_alignments` from `read_metrics_summary.tsv`.

### `read_length_histograms.tsv`

This file contains read-length histograms including soft-clipped bases.

```text
read_length    mapped_read_count    mapped_partial_read_count    mapped_full_read_count    post_filter_read_count    post_filter_partial_read_count    post_filter_full_read_count
```

| Metric name | Definition |
|:-------------|:------------|
| `read_length` | Read length including soft-clipped bases. |
| `mapped_read_count` | Number of mapped reads of the given length. |
| `mapped_partial_read_count` | Number of mapped partial reads of the given length. |
| `mapped_full_read_count` | Number of mapped full reads of the given length. |
| `post_filter_read_count` | Number of mapped reads of the given length that also pass the mapping-quality and flag filters. |
| `post_filter_partial_read_count` | Number of mapped partial reads of the given length that also pass the filters. |
| `post_filter_full_read_count` | Number of mapped full reads of the given length that also pass the filters. |

Constraints:

- `mapped_read_count = mapped_partial_read_count + mapped_full_read_count` when read-type information is available
- `post_filter_read_count = post_filter_partial_read_count + post_filter_full_read_count` when read-type information is available

### `read_length_excluding_soft_clips_histograms.tsv`

This file matches `read_length_histograms.tsv`, but uses `read_length_excluding_soft_clipped_bases` as the first column and computes the histogram over read lengths that exclude soft-clipped bases.

### `read_length_summary.tsv`

This file contains summary statistics for read length including soft-clipped bases.

```text
metric_name    mapped_read_length    mapped_partial_read_length    mapped_full_read_length    post_filter_read_length    post_filter_partial_read_length    post_filter_full_read_length
```

| Metric name | Definition |
|:-------------|:------------|
| `mean` | Mean read length. |
| `median` | Median read length. |
| `stddev` | Standard deviation of read length. |
| `min` | Minimum read length. |
| `max` | Maximum read length. |
| `percentile_X` | Read length at the Xth percentile. |

### `read_length_excluding_soft_clips_summary.tsv`

This file matches `read_length_summary.tsv`, but computes the same statistics over read lengths that exclude soft-clipped bases.

### `read_counts_by_cluster_size.tsv`

This file contains post-filter mapped-read counts by consensus-cluster size and cluster composition.

| Metric name | Definition |
|:-------------|:------------|
| `cluster_size` | Number of reads in the cluster used to generate the consensus read. |
| `mixed_strand_cluster_reads_passing_filter` | Number of mixed-strand cluster reads that pass the filters. |
| `forward_cluster_reads_passing_filter` | Number of forward-cluster reads that pass the filters. |
| `reverse_cluster_reads_passing_filter` | Number of reverse-cluster reads that pass the filters. |
| `full_cluster_reads_passing_filter` | Number of full-cluster reads that pass the filters. |
| `mixed_full_and_partial_cluster_reads_passing_filter` | Number of mixed full-and-partial cluster reads that pass the filters. |
| `partial_cluster_reads_passing_filter` | Number of partial-cluster reads that pass the filters. |

Back to [`read` output](usage-read.md#output) and [`all` output](usage-all.md#output).

## GC bias metrics

Applies to: `gcbias`, and `coverage` or `all` when `--enable-gc-bias-metrics` is set.

### `gc_bias.tsv`

This file contains per-GC-percentage-bin metrics that quantify the relationship between GC content and read coverage.

```text
gc_percentage    window_count    read_count    post_filter_read_count    normalized_coverage    normalized_post_filter_coverage    error_bar    post_filter_error_bar
```

| Metric name | Definition |
|:-------------|:------------|
| `gc_percentage` | GC content bin (0 to 100). |
| `window_count` | Number of reference windows at this GC percentage. |
| `read_count` | Number of reads whose start position maps to a window at this GC percentage. |
| `post_filter_read_count` | Number of reads passing read-level filters whose start position maps to a window at this GC percentage. |
| `normalized_coverage` | (reads_in_bin / windows_in_bin) / mean_reads_per_window. A value of 1.0 indicates no bias at this GC level. |
| `normalized_post_filter_coverage` | Normalized coverage computed from post-filter read counts. |
| `error_bar` | Error bar width: sqrt(reads_in_bin) / windows_in_bin / mean_reads_per_window. |
| `post_filter_error_bar` | Error bar width computed from post-filter read counts. |

### `gc_bias_summary.tsv`

This file contains summary-level GC bias metrics.

```text
metric_name    value
```

| Metric name | Definition |
|:-------------|:------------|
| `total_reads` | Total reads considered for GC bias analysis. |
| `total_windows` | Total reference windows considered. |
| `at_dropout` | AT dropout: sum of positive per-bin dropout values for GC bins 0 through 50. Higher values indicate under-representation of AT-rich regions. |
| `gc_dropout` | GC dropout: sum of positive per-bin dropout values for GC bins 51 through 100. Higher values indicate under-representation of GC-rich regions. |
| `first_quintile_normalized_coverage` | Normalized coverage for GC bins 0-19. |
| `second_quintile_normalized_coverage` | Normalized coverage for GC bins 20-39. |
| `third_quintile_normalized_coverage` | Normalized coverage for GC bins 40-59. |
| `fourth_quintile_normalized_coverage` | Normalized coverage for GC bins 60-79. |
| `fifth_quintile_normalized_coverage` | Normalized coverage for GC bins 80-100. |
| `post_filter_total_reads` | Post-filter total reads. |
| `post_filter_at_dropout` | AT dropout computed from post-filter read counts. |
| `post_filter_gc_dropout` | GC dropout computed from post-filter read counts. |
| `post_filter_*_quintile_normalized_coverage` | Quintile normalized coverages computed from post-filter read counts. |

Back to [`gcbias` output](usage-gcbias.md#output), [`coverage` output](usage-coverage.md#output), and [`all` output](usage-all.md#output).
