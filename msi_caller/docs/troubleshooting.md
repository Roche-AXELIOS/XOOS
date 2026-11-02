# Troubleshooting

## Output files already exist

**Symptom:** The command stops with `Output file(s) already exist (use --overwrite to replace)`.

**Cause:** At least one resolved VCF, JSON, or TSV output path already exists, and `--overwrite` is not set.

**Resolution:** Select a different `--output-dir`, prefix, or explicit output filename, or add `--overwrite` to replace all conflicting outputs.

## Insufficient coverage

**Symptom:** The aggregate status or score is `INSUFFICIENT_COVERAGE`, and some VCF records have the `INSUFFICIENT_COVERAGE` filter.

**Cause:** A site is unusable when its tumor distribution, or its matched-normal distribution in `tumor-normal` mode, fails the configured read filters.
The aggregate score is unavailable when fewer than `--min-sites` sites remain usable.

**Resolution:** Confirm that the BAM covers the configured sites, inspect the read-filter settings, and verify that the background model or sites file matches the analyzed data.
Change `--min-depth` or `--min-sites` only when the resulting limits are appropriate for the analysis.

## Aggregate threshold is unspecified

**Symptom:** The aggregate status is `THRESHOLD_UNSPECIFIED` even though the VCF contains stable and unstable site calls.

**Cause:** The command can classify individual sites without `--min-msi-classification-fraction`, but it needs this option to convert the MSI score into an aggregate `MSI-H` or `MSS` status.

**Resolution:** Set `--min-msi-classification-fraction 0.3` when an aggregate classification is required. Use another calibrated value only when appropriate for the analysis.

## Threshold column does not match the metric

**Symptom:** The command reports that the sites file is missing the threshold column required by `--microsatellite-sites-thresholds-tsv`.

**Cause:** The threshold-bearing sites file lacks a column whose name matches `--distribution-comparison-metric`.
For example, `--distribution-comparison-metric js-distance` requires a `js-distance` column.

**Resolution:** Select the metric that matches the threshold column in the sites file.

## BAM or BAM index cannot be read

**Symptom:** Input validation reports that a BAM file does not exist, is not readable, or is not indexed.

**Cause:** The BAM path is incorrect, the BAM is unavailable, or its adjacent BAI index is missing or unreadable.

**Resolution:** Provide a readable, coordinate-sorted BAM and its BAI index.
For *sample.bam*, place the index at *sample.bam.bai* or another index path supported by the BAM reader.

## Simplex input is rejected

**Symptom:** The command reports that SBX-S input is not supported.

**Cause:** MSI Caller requires a BAM identified as duplex data.

**Resolution:** Run the command with a supported duplex BAM.
