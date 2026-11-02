# Metrics

Demux writes summary metrics to a `metrics/` directory under `--output-dir`. For
how these metric files are produced and where they appear in the workflow, see
[Simplex demultiplexing](usage-simplex.md#output) and
[Duplex demultiplexing](usage-duplex.md#output). Format details for non-metric
inputs and outputs live on the [File formats](file-formats.md) page. The metric
set depends on the pipeline: [simplex](usage-simplex.md) and
[duplex](usage-duplex.md) share some metrics but each has its own run- and
sample-level schema and its own set of distribution files. Metrics common to
both pipelines are described once in [Adapter code metrics](#adapter-code-metrics).

## Simplex metrics

Applies to: `demux` simplex.

### Simplex run metrics

*run_stats.tsv* contains run-level summary statistics. The first line of the TSV
file is the following header:

```text
metric_name    value   percentage   denominator
```

The percentage is calculated by taking the value, dividing it by the value of the
`denominator` metric, and multiplying by 100. The `denominator` column names that
base metric. Both `percentage` and `denominator` are `NA` for rows that are not
percentages (for example `mean_passing_read_length`). The `percentage` is
also `NA` when the denominator is zero.

| Metric name | Description | Percentage denominator | Relationship |
|:------------|:------------|:-----------------------|:-------------|
| total_reads | Number of reads parsed | total_reads | `preassignment_passing_reads + failed_reads` |
| preassignment_passing_reads | Number of reads passing filter before assignment | total_reads | `assigned_reads + unassigned_reads` |
| failed_reads | Number of reads filtered before adapter processing | total_reads | `too_short_reads + too_long_reads` |
| too_short_reads | Number of reads filtered by min read length | total_reads | |
| too_long_reads | Number of reads too long to fit in the fixed-size record buffer | total_reads | |
| assigned_reads | Number of reads assigned to a SID | total_reads | `passing_reads + failed_assigned_reads` |
| unassigned_reads | Number of reads that could not be assigned to a SID | total_reads | |
| passing_reads | Number of SID-assigned reads that pass all post-assignment filters and are written to the SID-assigned output files | total_reads | `full_simplex_reads + partial_simplex_reads` |
| failed_assigned_reads | Number of reads that have valid adapters/SID but subsequently failed | total_reads | `too_short_trimmed_reads + sid_discordant_discarded_reads + full_only_discarded_reads + umi_failed_reads` |
| too_short_trimmed_reads | Number of assigned reads filtered by min length after trim | total_reads | |
| sid_discordant_discarded_reads | Number of discordant SID reads discarded by `--discordant-sid-mode` | total_reads | |
| full_only_discarded_reads | Number of partial reads discarded because `--read-length-mode full-only` is set | total_reads | |
| umi_failed_reads | Number of UMI-simplex (YSU) reads assigned by SID but discarded because the UMI detection requirement was not met | total_reads | |
| full_simplex_reads | Number of passing reads with a landmark (SID or UMI) detected on both ends | total_reads | counted for `passing_reads` only |
| partial_simplex_reads | Number of passing reads missing a landmark on at least one end | total_reads | counted for `passing_reads` only |
| both_sid_detected_reads | Number of reads with SID detected on both ends | total_reads | `sid_concordant_reads + sid_discordant_reads` |
| sid_concordant_reads | Number of reads with SIDs detected on both ends and in agreement | both_sid_detected_reads | |
| sid_discordant_reads | Number of reads with SIDs detected on both ends but in disagreement | both_sid_detected_reads | |
| perfect_index_reads | Number of reads with concordant, perfect (zero edit distance) SID matches on both ends | both_sid_detected_reads | |
| raw_bases | Raw (pre-trim) bases across all `total_reads` | raw_bases | `total_bases + trimmed_bases + unassigned_bases + failed_assigned_bases + filtered_bases` |
| total_bases | Bases in `passing_reads` after trimming | raw_bases | |
| trimmed_bases | Bases removed by trimming from `passing_reads` | raw_bases | |
| filtered_bases | Raw bases of `filtered_reads` | raw_bases | |
| unassigned_bases | Raw bases of `unassigned_reads` | raw_bases | |
| failed_assigned_bases | Raw bases of `failed_assigned_reads` | raw_bases | |
| mean_passing_read_length | Mean length of `passing_reads` | NA | `total_bases / passing_reads` |
| num_expected_sids | Number of expected SIDs in sample sheet | num_expected_sids | |
| num_sids | Number of unique SIDs detected in run | num_expected_sids | |

The [adapter code metrics](#adapter-code-metrics) (16 combinations) are also
reported at the run level.

### Simplex sample metrics

*sample_stats.tsv* contains per-sample statistics. Unlike `run_stats.tsv`, this
file is transposed: it has a leading `metric` column, followed by one column per
SID and an `Unassigned` column, so each row holds that metric's values across the
samples.

| Name | Description |
|:-----|:------------|
| index_sequence | Canonical adapter-bundle SID without a sample sheet; otherwise, the supplied `sample_sid`.<sup>*</sup> |
| assigned_reads | Number of reads assigned to this SID |
| assigned_reads_percentage | Percentage of reads assigned to this SID. Denominator: run-level `total_reads`. |
| full_simplex_reads | Number of this SID's passing reads with a landmark (SID or UMI) detected on both ends |
| full_simplex_reads_percentage | Percentage of full simplex reads for this SID. Denominator: sample-level `assigned_reads`. |
| partial_simplex_reads | Number of this SID's passing reads missing a landmark on at least one end |
| partial_simplex_reads_percentage | Percentage of partial simplex reads for this SID. Denominator: sample-level `assigned_reads`. |
| both_sid_detected_reads | Number of reads with SID detected on both ends |
| both_sid_detected_reads_percentage | Percentage of reads with SID detected on both ends. Denominator: sample-level `assigned_reads`. |
| sid_concordant_reads | Number of reads with SIDs detected on both ends and in agreement |
| sid_concordant_reads_percentage | Percentage of reads with concordant SIDs. Denominator: sample-level `both_sid_detected_reads`. |
| sid_discordant_reads | Number of reads with discordant SID |
| sid_discordant_reads_percentage | Percentage of reads with discordant SIDs. Denominator: sample-level `both_sid_detected_reads`. |
| perfect_index_reads | Number of reads with concordant, perfect (zero edit distance) SID matches on both ends |
| perfect_index_reads_percentage | Percentage of reads with concordant, perfect SID matches on both ends. Denominator: sample-level `both_sid_detected_reads`. |
| umi_failed_reads | Number of UMI-simplex (YSU) reads assigned by SID but discarded because the UMI detection requirement was not met |
| umi_failed_reads_percentage | Percentage of UMI-simplex reads that failed the UMI detection requirement. Denominator: sample-level `assigned_reads`. |
| passing_reads | Number of passing reads assigned to this SID |
| passing_reads_percentage | Percentage of passing reads assigned to this SID. Denominator: sample-level `assigned_reads`. |
| passing_reads_percentage_all_samples | Percentage of passing reads for this SID across configured samples. Denominator: sum of sample-level `passing_reads` across configured samples. |
| total_bases | Insert (post-trim) bases of this SID's passing reads |
| total_bases_percentage_all_samples | Percentage of total bases for this SID across configured samples. Denominator: sum of sample-level `total_bases` across configured samples. |
| trimmed_bases | Bases removed by trimming from this SID's passing reads |
| trimmed_bases_percentage | Percentage of bases removed by trimming from this SID's passing reads. Denominator: sample-level `total_bases + trimmed_bases`. |
| mean_passing_read_length | Mean passing read length for this SID (`total_bases / passing_reads`) |

Note that the denominator for a sample-level percentage can differ from that of its
corresponding run-level percentage metric.

<sup>*</sup> For legacy sample-sheet support, Demux reports canonical SID sequences when every supplied SID is uniformly reverse-complemented or uses the same recognized adapter-structure extension; Demux logs a warning for either conversion.

The [adapter code metrics](#adapter-code-metrics) (16 combinations) are also
reported for each sample, so the file has 41 (13 base rows, 12 percentage rows,
and 16 adapter-code rows).

### Sample assignment metrics

*sample_assignment_metrics.tsv* reports SID collision statistics (7 columns).

| Name | Description |
|:-----|:------------|
| expected_sid | Expected SID name |
| collision_sid | Collided SID name |
| 5'collision | Number of 5′ SID collisions |
| 3'collision | Number of 3′ SID collisions |
| total_collisions | Number of SID collisions; sum of `5'collision` and `3'collision` |
| both_sid_detected_reads | Number of reads with both 5′ and 3′ SIDs detected |
| percent_sid_collision | Percent of SID collisions; `total_collisions` divided by `both_sid_detected_reads` |

### Simplex read length distribution metrics

Simplex emits run-level distributions in `metrics/read_length_distr/`:

1. `total_read_length_distr.tsv` — all raw input reads.
2. `raw_failed_read_length_distr.tsv` — reads that would be written to `raw_failed` when `--write-failed-reads` is enabled.
3. `total_untrimmed_simplex_passing_read_length_distr.tsv` — raw lengths of passing reads.
4. `total_trimmed_simplex_passing_read_length_distr.tsv` — insert lengths of passing reads after trimming.
5. `total_endadapter_position_distr.tsv` — 1-based 3′ adapter-start positions for all passing reads; the `0` bin represents reads without a 3′ adapter.

The two passing-read run-level files have `read_length`, `passing_reads`,
`full_passing_reads`, and `partial_passing_reads` columns. Each bin satisfies
`passing_reads = full_passing_reads + partial_passing_reads`.

Per-sample distributions are written in `metrics/read_length_distr/<sample-name>_read_length_distr/`,
where `<sample-name>` matches the Demux output directory. Passing-read distributions use the same three count
columns as their run-level counterparts. Each sample additionally has an
`sample_endadapter_position_distr.tsv` file.

Lengths exceeding `--max-read-length-bin` are grouped into a single
bucket (e.g., `1000+`).

{% hint style="info" %}
When `--adapter-design-name SIMPLEX-10X` is used for 10x scRNA libraries, only the
5′ SID is detected and trimmed. All assigned reads are therefore reported as
`partial_simplex_reads` with adapter code `s5xXxXxX`, and metrics that depend on
3′ SID detection (`full_simplex_reads`, `both_sid_detected_reads`,
`sid_concordant_reads`, `sid_discordant_reads`, `perfect_index_reads`) report a
value of zero.
{% endhint %}

### Adapter code metrics

Applies to: `demux` simplex.

Adapter metrics use an eight-letter code combining `s5`/`s3`/`u5`/`u3` (found) and
`xX` (not found) for 5′ SID, 3′ SID, 5′ UMI, and 3′ UMI respectively — 16
combinations total. They are reported in the simplex
[run](#simplex-run-metrics) and [sample](#simplex-sample-metrics) metric files.

| Name | 5′ SID | 3′ SID | 5′ UMI | 3′ UMI |
|:-----|:-------|:-------|:-------|:-------|
| s5s3u5u3 | yes | yes | yes | yes |
| s5s3u5xX | yes | yes | yes | no |
| s5s3xXu3 | yes | yes | no | yes |
| s5s3xXxX | yes | yes | no | no |
| s5xXu5u3 | yes | no | yes | yes |
| s5xXu5xX | yes | no | yes | no |
| s5xXxXu3 | yes | no | no | yes |
| s5xXxXxX | yes | no | no | no |
| xXs3u5u3 | no | yes | yes | yes |
| xXs3u5xX | no | yes | yes | no |
| xXs3xXu3 | no | yes | no | yes |
| xXs3xXxX | no | yes | no | no |
| xXxXu5u3 | no | no | yes | yes |
| xXxXu5xX | no | no | yes | no |
| xXxXxXu3 | no | no | no | yes |
| xXxXxXxX | no | no | no | no |

Back to [Simplex demultiplexing output](usage-simplex.md#output).

## Duplex metrics

Applies to: `demux` duplex.

### Duplex run metrics

*run_stats.tsv* contains run-level summary statistics. The first line of the TSV
file is the following header:

```text
metric_name    value   percentage   denominator
```

The percentage is calculated by taking the value, dividing it by the value of the
`denominator` metric, and multiplying by 100. The `denominator` column names that
base metric. Both `percentage` and `denominator` are `NA` for rows that are not
percentages (for example `mean_passing_read_length`). The `percentage` is
also `NA` when the denominator is zero.

| Metric name | Description | Percentage denominator | Relationship |
|:------------|:------------|:-----------------------|:-------------|
| total_reads | Raw reads processed. | total_reads | `assigned_reads + unassigned_reads` |
| assigned_reads | Reads with a hairpin that contains valid SIDs and are not filtered by raw read length. | total_reads | `passing_reads + failed_assigned_reads`, `found_by_string_compare + found_by_global_symmetry + found_by_local_symmetry` |
| passing_reads | Assigned reads that pass all filters such as trimmed length, error rate and buffer size. | total_reads | `full_duplex_reads + partial_duplex_reads` |
| full_duplex_reads | Passing reads without a simplex region. | total_reads | |
| partial_duplex_reads | Passing reads with a simplex region. | total_reads | |
| unassigned_reads | Reads filtered before or during SID assignment. | total_reads | `no_hairpin_reads + too_long_reads + too_short_reads` |
| no_hairpin_reads | Unassigned reads without a hairpin. | total_reads | |
| too_long_reads | Unassigned reads too long to fit in the fixed-size record buffer. | total_reads | |
| too_short_reads | Unassigned reads too short according to length filter (default 50 bp). | total_reads | |
| failed_assigned_reads | Assigned reads that failed subsequent filters. | total_reads | `too_many_errors_reads + too_short_trimmed_reads + too_long_consensus_reads + failed_hairpin_stem_trim_reads` |
| too_many_errors_reads | Failed assigned reads with too many errors in duplex region (default 10%). | total_reads | |
| too_short_trimmed_reads | Failed assigned reads filtered after trimming (default 50 bp). | total_reads | |
| too_long_consensus_reads | Failed assigned reads with consensus sequence too long for output buffer (currently 1344). | total_reads | |
| failed_hairpin_stem_trim_reads | Assigned read that failed to trim on the hairpin side. | total_reads | |
| found_by_string_compare | Hairpins found via string comparison algorithm. | assigned_reads | |
| found_by_global_symmetry | Hairpins found via global symmetry algorithm. | assigned_reads | |
| found_by_local_symmetry | Hairpins found via local symmetry algorithm. | assigned_reads | |
| longer_r2_reads | Reads where raw R1 is shorter than raw R2. | passing_reads | |
| longer_r2_full_duplex_reads | Full duplex reads where raw R1 is shorter than raw R2. | passing_reads | |
| both_umi_reads | Reads where both UMIs (relative to R1) were found. | passing_reads | |
| 5p_umi_reads | Reads where only the 5p UMI (relative to R1) was found. | passing_reads | |
| 3p_umi_reads | Reads where only the 3p UMI (relative to R1) was found. | passing_reads | |
| no_endadapter_reads | Passing reads without a detectable endadapter after consensus (left untrimmed). | passing_reads | |
| strand_fw_reads | Forward strand classified reads (`--strand-index` mode only). | passing_reads | |
| strand_rv_reads | Reverse strand classified reads (`--strand-index` mode only). | passing_reads | |
| strand_fw_sig_reads | Forward strand significantly classified reads (`--strand-index` mode only). | passing_reads | |
| strand_rv_sig_reads | Reverse strand significantly classified reads (`--strand-index` mode only). | passing_reads | |
| total_bases | Total bases in all passing reads. Its percentage of `raw_bases` shows base conversion efficiency to passing output. | raw_bases | `duplex_bases + non_duplex_bases` |
| non_duplex_bases | Bases in non-duplex regions in passing reads. | total_bases | |
| duplex_bases | Duplex bases that are concordant, ambiguous, or discordant. | total_bases | `concordant_duplex_bases + ambiguous_duplex_bases + discordant_duplex_bases` |
| concordant_duplex_bases | Bases matching duplex pair in passing reads. | total_bases | |
| discordant_duplex_bases | Bases not matching in passing reads. | total_bases | |
| ambiguous_duplex_bases | Bases at topologically ambiguous positions in passing reads. Reported only when the ambiguous topology quality mode is enabled. | total_bases | |
| raw_bases | Raw bases in all input reads. | raw_bases | `unassigned_bases + failed_assigned_bases + trimmed_bases + total_bases` |
| unassigned_bases | Raw bases in reads that were determined unassigned. | raw_bases | |
| failed_assigned_bases | Raw bases in reads that were assigned but failed subsequent filters. | raw_bases | |
| trimmed_bases | Bases removed from passing reads during adapter trimming. | raw_bases | |
| num_expected_sids | Number of expected SIDs in sample sheet. | num_expected_sids | |
| num_sids | Number of unique SIDs detected in run. | num_expected_sids | |
| multi_candidate_alignments | Passing reads where multiple equal-edit-distance alignment candidates existed. | passing_reads | |
| mean_passing_read_length | Mean read length of passing reads. | NA | `total_bases / passing_reads` |

### Duplex sample metrics

*sample_stats.tsv* contains per-sample statistics in a transposed layout.
The leading `metric` column identifies the metric, followed by one column per SID.

Except for the following differences, Duplex sample metrics use the same metric names and definitions as `run_stats.tsv`:

| Name | Description |
|:-----|:------------|
| index_sequence | Canonical adapter-bundle SID without a sample sheet; otherwise, the supplied `sample_sid`.<sup>*</sup> |

<sup>*</sup> For legacy sample-sheet support, Demux reports canonical SID sequences when every supplied SID is uniformly reverse-complemented or uses the same recognized adapter-structure extension; Demux logs a warning for either conversion.

The following run-level metrics are **not** reported per sample:

- `total_reads`
- `unassigned_reads`
- `no_hairpin_reads`
- `too_long_reads`
- `too_short_reads`
- `raw_bases`
- `unassigned_bases`
- `num_expected_sids`
- `num_sids`

In addition to the metrics inherited from `run_stats.tsv`, Duplex sample metrics include the following percentage metrics immediately after their associated count/base rows. Their denominators can differ from the corresponding run-level percentage metrics.

| Percentage metric | Percentage denominator |
|:------------------|:------------|
| assigned_reads_percentage | Run-level `total_reads` |
| passing_reads_percentage | Sample-level `assigned_reads` |
| passing_reads_percentage_all_samples | Sum of sample-level `passing_reads` across configured samples |
| full_duplex_reads_percentage | Sample-level `assigned_reads` |
| partial_duplex_reads_percentage | Sample-level `assigned_reads` |
| failed_assigned_reads_percentage | Sample-level `assigned_reads` |
| both_umi_reads_percentage | Sample-level `passing_reads` |
| 5p_umi_reads_percentage | Sample-level `passing_reads` |
| 3p_umi_reads_percentage | Sample-level `passing_reads` |
| no_endadapter_reads_percentage | Sample-level `passing_reads` |
| total_bases_percentage_all_samples | Sum of sample-level `total_bases` across configured samples |
| concordant_duplex_bases_percentage | Sample-level `total_bases` |
| discordant_duplex_bases_percentage | Sample-level `total_bases` |
| ambiguous_duplex_bases_percentage | Sample-level `total_bases` |
| duplex_bases_percentage | Sample-level `total_bases` |
| non_duplex_bases_percentage | Sample-level `total_bases` |

### Duplex length distribution metrics

Demux outputs read and region size distributions in TSV format. Each filename
corresponds to a metric above.

Pre-SID-assignment distribution files:

1. `total_read_length_distr.tsv`
2. `raw_failed_read_length_distr.tsv`
3. `no_hairpin_read_length_distr.tsv`
4. `total_untrimmed_duplex_passing_read_length_distr.tsv`
5. `total_trimmed_duplex_passing_read_length_distr.tsv`
6. `total_endadapter_position_distr.tsv`

The passing files contain `read_length`, `passing_reads`, `full_passing_reads`,
and `partial_passing_reads`; each bin satisfies `passing_reads =
full_passing_reads + partial_passing_reads`.

Per-sample files are written in `metrics/read_length_distr/<sample-name>_read_length_distr/`:

1. `sample_untrimmed_duplex_passing_read_length_distr.tsv`
2. `sample_trimmed_duplex_passing_read_length_distr.tsv`
3. `sample_endadapter_position_distr.tsv`

Back to [Duplex demultiplexing output](usage-duplex.md#output).
