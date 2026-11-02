# Troubleshooting

## Excessive memory usage

**Symptom:** The run fails with `std::bad_alloc` or another out-of-memory failure.

**Cause:** Alignment Metrics processes metrics by region.
If the effective region size is too large for the available memory,
the per-thread working set can exceed memory limits.

**Resolution:** Reduce `--region-size`, reduce `--threads`, or provide more memory.

## Excessive memory usage with target-enrichment metrics

**Symptom:** A run with `--enable-target-enrichment-metrics` uses much more memory than expected.

**Cause:** Target-enrichment metrics disable region partitioning so the Alignment Metrics tool can preserve correct target-region counting and mean-coverage calculations.
Large BED intervals therefore stay large.

**Resolution:** Confirm the BED file is a true target-panel BED with small intervals.
If the intervals are genuinely large,
reduce `--threads`,
increase memory,
or split the BED intervals only if that tradeoff is acceptable for the analysis.

## Unexpectedly high post-filter coverage in duplex outputs

**Symptom:** `post_filter_coverage` is higher than expected in duplex outputs.

**Cause:** Post-filter coverage means the bases passed the configured base-call quality and base-type filters.
It does not imply the bases are strictly concordant duplex bases.

**Resolution:** Review `--min-base-type` and `--min-base-quality` for the Alignment Metrics tool.
For duplex BAMs, `--min-base-type concordant` is usually the intended strict setting.

## Post-filter coverage is 0 or `NA`

**Symptom:** Post-filter coverage outputs are 0 or `NA`.

**Cause:** Simplex datasets,
consensus datasets,
and BAMs without YC tags treat all bases as simplex bases.
If `--min-base-type` is set to `concordant` in that situation,
no bases pass the filter.

**Resolution:** Set `--min-base-type simplex` for simplex,
consensus,
or YC-tag-free BAMs.

## Malformed YC tag

**Symptom:** The run fails while decoding a YC tag from a BAM record.

**Cause:** The BAM contains a malformed YC tag.

**Resolution:** Inspect the error message for the offending read name and fix the upstream BAM if possible.
If base-type decoding is not required for the analysis,
rerun with `--skip-base-type-decoding`.

## Error reading reference sequence for region

**Symptom:** The run fails with an error indicating that a reference sequence could not be retrieved for a specific region.

**Cause:** The reference FASTA provided to the Alignment Metrics tool does not contain the requested region, either because the contig is missing or the contig does not fully cover the requested region.

**Resolution:** Ensure that the reference FASTA provided to the Alignment Metrics tool is the same one used during the alignment that produced the BAM file and that it contains all contigs present in the BAM header with the correct lengths.

## Contig listed in BAM header but has no mapped reads and is missing from the reference FASTA

**Symptom:** The run logs a warning indicating that a contig listed in the BAM header has no mapped reads and is missing from the reference FASTA.

**Cause:** The contig is listed in the BAM header but has no mapped reads and is missing from the reference FASTA. Even if a contig has no mapped reads, the Alignment Metrics tool still needs to access the reference sequence for coverage metrics (for counting positions with zero coverage) and GC bias metrics calculation (for counting GC windows).

**Resolution:** Ensure that the reference FASTA provided to the Alignment Metrics tool is the same one used during the alignment that produced the BAM file and that it contains all contigs present in the BAM header with the correct lengths. If not calculating GC bias metrics, this warning does not cause the run to fail because contig length can be approximated from the BAM header even if the reference sequence is missing. However, if GC bias metrics are enabled, the run will fail because the Alignment Metrics tool cannot collect GC windows and calculate GC content for the missing contig. It is therefore important to provide a complete and correct reference FASTA instead of ignoring the warning message.
