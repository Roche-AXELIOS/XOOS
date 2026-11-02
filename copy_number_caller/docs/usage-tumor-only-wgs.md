# Call tumor-only somatic copy numbers

## How it works

The `tumor-only-wgs` subcommand characterizes somatic CNVs (aka. somatic copy number variations) from a tumor sample when no matched germline normal sample is available.
Like the tumor-normal workflow, it is optimized for oncology research, but it derives the reference baseline and germline heterozygous sites from the tumor sample itself rather than from a matched normal.

Without a matched normal, the two signals used for calling are obtained differently.
The log-ratio signal is computed by self-normalizing the tumor coverage (each interval against the median autosomal coverage) instead of dividing tumor coverage by normal coverage.
The B-allele fraction (BAF) signal is extracted at likely germline heterozygous sites, which are identified by intersecting the single-sample tumor VCF with a population SNP panel (for example, gnomAD) supplied through `--population-vcf`.

If the tumor BAM has a median coverage below 10x, or if the input VCF has insufficient data quality (too few germline SNPs or zero variants after filtering), the module exits with code 0 and writes structurally valid but data-empty output files. This avoids numerical errors during GC correction on very low or zero-coverage samples and prevents unreliable calls from poor-quality VCF input.

Expand the details below for the step-by-step design of this command.

<details>
<summary>Design details.</summary>

1. **Gather input files**, which include the aligned reads from the tumor sample in BAM format, a single-sample tumor VCF, a population SNP VCF (for germline heterozygous site identification), as well as genome resource files (for example, a FASTA reference and a mappability file).
2. **Calculate coverage** of the aligned reads across fixed-width intervals over the whole genome.
3. **Correct and normalize coverage** to mitigate bias introduced by regions with high or low GC content, as well as to mitigate outliers and noise. Without a matched normal, the caller normalizes data by self-normalizing the tumor coverage against the median autosomal coverage.
4. **Segment the genome** into stretches of adjacent intervals with similar coverage. The first round segments by the self-normalized log-ratios, and the second round further segments by the B-allele fractions at likely germline heterozygous sites.
5. **Call copy number** of each segment using a likelihood model. Purity and ploidy are estimated from the SCNV grid search.
6. **Write** output summary files and metrics.

</details>

## Example command

```bash
copy_number_caller tumor-only-wgs \
        --tumor-bam /path/to/tumor_bam \
        --tumor-sample-name tumor_sample_name \
        --vcf /path/to/tumor/variant/calls \
        --population-vcf /path/to/population_vcf \
        --reference /path/to/reference \
        --threads 16 \
        --output-dir /path/to/output
```

For more details on available CLI options, see the [CLI options](#cli-options) section below.

## Input

| Input | Required | Notes |
|:---------|:-----------------------|:-----------------------|
| Tumor BAM | Yes | Position-sorted, indexed tumor BAM (`--tumor-bam`). |
| Tumor VCF | Yes | Bgzip-compressed, tabix-indexed, position-sorted single-sample tumor VCF of variant sites for germline heterozygous site identification, [BAF](file-formats.md#bafs) extraction, and BAF visualization (`--vcf`). Records must be sorted by position within each contig and kept in contiguous contig blocks (contig block order can be arbitrary). The `.tbi` index must be next to the VCF. |
| Population SNP VCF | Yes | Bgzip-compressed population SNP VCF and its co-located `.tbi` index (`--population-vcf`), used to identify likely germline heterozygous sites. |
| Reference FASTA | Yes | Indexed reference genome (`--reference`). |
| Mappability file | No | [Genome mappability file](appendix.md#genome-mappability) (`--mappability-bigwig`). When omitted, the caller selects the bundled resource from the library type detected from the tumor BAM. |
| Seed segments | No | [Seed segment file](appendix.md#seed-segments) (`--seed-segments`). When omitted, uses `/resources/hg38_seed.seg`; an explicit value overrides it. |
| Tumor sample name | Yes | Tumor sample name as it appears in the input VCF (`--tumor-sample-name`). |
| Blocklist BED | No | 0-based regions to exclude before generating intervals (`--blocklist-bed`). When omitted, uses the [hg38 tumor-only WGS blocklist](appendix.md#hg38-tumor-only-wgs-blocklist-bed-file); an explicit value overrides it. |
| SNP blocklist BED | No | Regions to exclude from BAF extraction (`--snp-blocklist-bed`). When omitted, the bundled hg38 resource at `/resources/hg38_snp_blocklist.bed` is selected automatically. |

The recommended way to generate the single-sample tumour VCF is to use Mutect2 in tumor-only mode (see [Tips for analysis](#tips-for-analysis) for more details).
Germline heterozygous calls are identified from this single-sample tumor VCF by intersecting it with the population panel provided through `--population-vcf`.
Only biallelic SNPs with a matching population-panel record at the same position and with the same REF and ALT alleles are considered. The matching record must have an `INFO/AF` value at or above `--min-population-af` (default: `0.01`).
For a site to be used, the `AD` (Allele Depth) tag must be present in the tumor `FORMAT` field.
The `AD` tag describes a vector of two values, corresponding to the read depths over the reference and alternative alleles, respectively.
The B-allele fraction (BAF) is derived by dividing the number of alternate reads (alt B reads) by the total number of reads (ref + alt reads).
The caller keeps only sites with total depth at or above `--min-tumor-depth` (default: `10`) and a BAF within the inclusive range set by `--min-tumor-baf`/`--max-tumor-baf` (default: `0.1 <= BAF <= 0.9`).
Sites overlapping the regions in `--snp-blocklist-bed` (for example, known problematic or recurrent-artifact loci) are excluded from BAF extraction.

The recommended population VCF is the GATK gnomAD germline resource.
Download both the VCF and its `.tbi` index so the index is present alongside the VCF:

```bash
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz .
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi .
```

## Output

Files are written to the current working directory by default.
The output folder can be specified with `--output-dir`.
A custom prefix (for example a sample ID) can be prepended to every output filename with `--output-prefix`; the prefix and the default name are joined with an underscore, and the subfolder is preserved (for example, `--output-prefix TO` produces `TO_somatic_copy_number_callset.vcf.gz` and `metrics/TO_metrics.tsv`). When `--output-prefix` is omitted, the default filenames below are used.
If any output file already exists, the command exits with an error unless `--overwrite` is specified.

Output files are organized into subfolders:

```text
<output-dir>/
├── somatic_copy_number_callset.seg
├── somatic_copy_number_callset.vcf.gz
├── somatic_copy_number_callset.vcf.gz.tbi
├── metrics/       # QC and statistics
├── tracks/        # visualization (BigWig, IGV session)
└── supplemental/  # purity/ploidy search solutions, raw tracks, coverage, intermediate segments
```

Principal output:

| Output | Description |
|:---------|:-----------------------|
| *somatic_copy_number_callset.seg*  | Copy number estimations for each segment along with their likelihoods in [SEG](file-formats.md#segments) format.          |
| *somatic_copy_number_callset.vcf.gz* | Copy number estimations for each segment in [VCF](file-formats.md#vcf) format. |
| *somatic_copy_number_callset.vcf.gz.tbi* | Tabix index for the compressed VCF. |
| *supplemental/purity_ploidy_grid.tsv* | Ranked purity ploidy solutions generated by the caller and their likelihoods in [purity ploidy grid](file-formats.md#purity-ploidy-grid) format. Not generated if `--tumor-purity` and `--tumor-ploidy` provided. |
| *tracks/log_ratios.bw*  | Output file for the self-normalized tumor log-ratios across the genome in BigWig format. |
| *tracks/bafs.bw*  |  Output file for the tumor BAF values at putative germline heterozygous SNPs across the genome in BigWig format.   |
| *tracks/igv_visualization.xml*  | Pre-configured IGV session file to visualize log ratios and BAF values.           |
| *metrics/metrics.tsv* | Summary metrics for copy number calls. |

Additional output:

| Output | Description |
|:---------|:-----------------------|
| *supplemental/augmented_intervals.bed*  | Output file for the whole genome intervals that were automatically generated in [intervals](file-formats.md#intervals) BED format. |
| *supplemental/tumor_corrected_coverage.bed*  | Output file for the GC-corrected tumor sample coverage in [coverage](file-formats.md#coverage) BED format. |
| *supplemental/tumor_coverage.bed*  | Output file for the raw tumor sample coverage in [coverage](file-formats.md#coverage) BED format. |
| *supplemental/log_ratio_segments.seg*  | Output file for the log ratios in [SEG](file-formats.md#segments) format. |
| *supplemental/log_ratios.bed*  | Output file for the log ratios in [log ratios](file-formats.md#log-ratios) BED format. |
| *supplemental/mapping_quality.bed*  | Output file for the average mapping quality per interval in [observations](file-formats.md#observations) BED format. |
| *supplemental/bafs.bed*  | Output file for the BAF values in [BAFs](file-formats.md#bafs) BED format. |

See [File formats](file-formats.md) for the schema of each output and [Metrics](metrics.md#somatic-metrics) for the contents of *metrics/metrics.tsv*. In tumor-only mode the `LogR` column holds the log2 of the self-normalized read depth ratios (each interval against the median autosomal coverage) rather than the tumor/normal ratio.

Notes:

- The purity/ploidy solution the caller uses to call allele-specific variants is printed to each line of the SEG file.
- The predicted sex of the individual is printed to the header of the VCF file and to each line of the SEG file.
  `--sex auto-detect` predicts from tumor coverage and is the default.
  Use `male` or `female` to override prediction, or `skip-prediction` to leave the resolved output sex as `unknown`.
- When the sex is `unknown`, the caller excludes allosome segments (chromosomes X and Y) from copy-number calling because it cannot determine whether those segments are haploid or diploid. Autosomal segments remain eligible for calling.
- The sample ID in the output is the `--tumor-sample-name` value.
- Absent values in BED, SEG, and TSV outputs are written as `NA`.

The tumor-only *metrics/metrics.tsv* uses the [somatic metrics](metrics.md#somatic-metrics); the three [tumor-normal-only rows](metrics.md#tumor-normal-only-metrics) are not present because there is no matched normal VCF.

## Tips for analysis

The following recommendations apply to tumor-only analysis and are documented in the appendix and troubleshooting guide:

- Generating variant calls — See the [Small Variant Caller — Tumor-only whole-genome](../../small_variant_caller/docs/generating-candidate-calls.md#tumor-only-whole-genome) section for the recommended Mutect2 command and parameters.
- [Excluding problematic regions](appendix.md#excluding-problematic-regions) — masking paralogous and other error-prone regions with the default `/resources/hg38_tumor_only_wgs_blocklist.bed` to protect BAF-based calling. Pass `--blocklist-bed` only to override that bundled resource (see [hg38 tumor-only WGS blocklist BED file](appendix.md#hg38-tumor-only-wgs-blocklist-bed-file)).
- [Preventing oversegmentation](troubleshooting.md#preventing-oversegmentation) — tuning `--segment-undoing` for noisy samples such as FFPE.
- [Dealing with low SCNV signal tumors](troubleshooting.md#low-scnv-signal-tumors) — inspect the log-ratio and BAF channels; calls from tumors with predicted purity below 30% should generally be considered unreliable. The somatic flagged VCF approach described in the troubleshooting guide applies to tumor-normal mode only.
- [Predicting SCNV in whole genome doubling tumors](appendix.md#whole-genome-doubling) — investigating alternative purity/ploidy solutions and using `predict-somatic-cnv` to regenerate calls.

## CLI options

This command also accepts the common `-h, --help`, `-v, --version`, and `-l, --log-level` (`debug`, `info`, `warn`, `error`; default `info`) options.
Parameters in **bold** are required.

### Input options

| Parameter | Description | Value(s) |
|:----|:----|:----|
| **--tumor-bam** | Path to the input tumor BAM file. The BAM file must be coordinate-sorted and indexed. | File path |
| **--vcf** | Path to the bgzip-compressed, position-sorted, single-sample tumor VCF file with variant sites. The VCF must have a co-located `.tbi` index. | File path |
| **--population-vcf** | Path to the bgzip-compressed population SNP VCF (e.g., gnomAD). Used to identify likely germline heterozygous sites. The VCF must have a co-located `.tbi` index. | File path |
| --seed-segments | Path to pre-generated segments file from which to start segmentation. The seed segments file must be in SEG format. When omitted, uses `/resources/hg38_seed.seg`; an explicit value overrides it. | File path |
| **--reference** | Path to the input reference FASTA file. The FASTA file must be indexed with a .fai file. | File path |
| --mappability-bigwig | Path to an input BigWig file describing mappability across the genome. When omitted, the library type is detected from `--tumor-bam`: SBX-S uses `/resources/hg38_500mer_mappability.bw`; SBX-D, SBX-FAST, and unresolved library types use `/resources/hg38_mappability.bw`. | File path; [default: library dependent] |
| --snp-blocklist-bed | Path to a BED file of 0-based regions whose overlapping SNP sites are excluded from BAF extraction (e.g., known problematic or recurrent-artifact loci). When omitted, uses `/resources/hg38_snp_blocklist.bed`. | File path |
| --blocklist-bed | Path to a BED file of 0-based regions to exclude from the analysis before generating intervals. When omitted, uses `/resources/hg38_tumor_only_wgs_blocklist.bed`; an explicit value overrides it. | File path |

### Output options

| Parameter | Description | Value(s) |
|:----|:----|:----|
| --output-dir | Path to output directory. | Directory path [default: .] |
| --output-prefix | Prefix prepended to all output filenames (for example a sample ID), joined with an underscore. If unset, default filenames are used. | Text [default: none] |
| --overwrite | Overwrite output files if they already exist. By default, the command exits with an error if an output file already exists. | Flag |

### Region options

| Parameter | Description | Value(s) |
|:----|:----|:----|
| --interval-size | The size (in base pairs) for each target region of the genome for which log ratios are calculated. | Uint; [default: 1000] |
| --min-interval-mappability | Minimum allowable mappability score for an interval. Intervals with a score less than this value are excluded to ensure high-quality read alignment. Mappability scores come from `--mappability-bigwig`. This often excludes low complexity regions (i.e. segdups). | Float; [default: 0.6] |

### Sample metadata options

| Parameter | Description | Value(s) |
|:----|:----|:----|
| **--tumor-sample-name** | Name of tumor sample in the input VCF header and used as the sample ID in the output. | Text |
| --sex | How to obtain the sample sex. `auto-detect` predicts from tumor coverage; `skip-prediction` records the resolved sex as `unknown`. | Enum: {auto-detect, male, female, skip-prediction}; [default: auto-detect] |
| --tumor-purity | Tumor purity. Use if the tumor purity of the sample is known prior to the analysis. Requires: --tumor-ploidy. | Float; Range: [0.15 - 1] |
| --tumor-ploidy | Tumor ploidy. Use if the tumor ploidy of the sample is known prior to the analysis. Requires: --tumor-purity. | Float:positive |

### Read filtering options

| Parameter | Description | Value(s) |
|:----|:----|:----|
| --exclude-flags | Exclude alignments if any bits in their FLAG field match the specified integer. Refer to <https://broadinstitute.github.io/picard/explain-flags.html> for flag generation. | Text; [default: 1796] |
| --min-tumor-coverage | Minimum allowable coverage in the tumor sample for an interval to be kept. Intervals below this value are excluded. | Int; [default: 0] |

### Segmentation options

| Parameter | Description | Value(s) |
|:----|:----|:----|
| --segment-undoing | Strategy for undoing log-ratio segments (re-merging segments with similar median log-ratio values). `iterative` relaxes `--initial-segment-undoing-factor` progressively each round (by an increasing multiple of `--segment-undoing-increment`) until `--max-segments-for-undoing` is reached; `non-iterative` performs a single merge pass at `--initial-segment-undoing-factor`; `disabled` skips undoing. | {disabled, iterative, non-iterative}; [default: iterative] |
| --initial-segment-undoing-factor | The initial standard deviation multiplier for segment undoing. Adjacent segments with a difference less than or equal to the threshold determined by this multiplier are merged. | Float:positive; [default: 3] |
| --max-segments-for-undoing | Target segment count for the segment undoing procedure (`iterative` only). The process stops once the total number of segments is less than or equal to this value. To skip segment undoing, set `--segment-undoing disabled` instead. | Integer:positive; [default: 500] |
| --hierarchical-pruning-height | Clustering parameter for hierarchical pruning. | Float; [default: 0.25] |
| --baf-distance-weight | Weight applied to the mirrored-BAF axis of the distance used for hierarchical pruning, applied only when both segments in a pair are copy-neutral. Values greater than 1 increase the separation of allelic-imbalance segments (for example, copy-neutral LOH) from balanced segments, making them less likely to be merged. | Float:positive; [default: 2] |

### Purity/ploidy search options

| Parameter | Description | Value(s) |
|:----|:----|:----|
| --baf-likelihood-weight-for-purity-ploidy-search-rerun | Weight applied to the BAF term of the joint likelihood during the [WGD uniform-prior re-run](appendix.md#whole-genome-doubling) only. The initial [diploid-bias search](appendix.md#whole-genome-doubling) is unaffected. Values greater than 1 give more weight to BAF than to log ratio. | Float:non-negative; [default: 1.1] |

### VCF parsing options

| Parameter | Description | Value(s) |
|:----|:----|:----|
| --min-population-af | Minimum population allele frequency to consider a site as a likely germline het. | Float; Range: [0 - 1]; [default: 0.01] |
| --min-tumor-depth | Minimum total depth (ref + alt) at a site in the tumor VCF. | Int:positive; [default: 10] |
| --min-tumor-baf | Minimum BAF to include in the output. | Float; Range: [0 - 1]; [default: 0.1] |
| --max-tumor-baf | Maximum BAF to include in the output. | Float; Range: [0 - 1]; [default: 0.9] |

### Performance options

| Parameter | Description | Value(s) |
|:----|:----|:----|
| --threads | Number of threads. | Int; [default: 1] |
