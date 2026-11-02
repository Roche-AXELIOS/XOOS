# STR caller

The STR Caller genotypes and detects repeat expansions in short tandem repeats (STRs) from Roche SBX sequencing data.
Input: **BAM** + **reference FASTA** → Output: **VCF**, **JSON**, **metrics**, and optional realigned and visualization outputs.

## How it works

<details>

<summary>Algorithm detail</summary>

The caller processes each STR catalog entry as an independent analysis task, parallelized across threads.
Read recruiters and references are pooled per thread; sample sex and name are inferred from the BAM when not provided.

The workflow has four layers — input, recruitment, analysis, and output.
For each catalog entry, an analysis task runs the following steps.

### 1. Create the local alignment

The alignment blueprint describes the repeat structure — each repeat unit and its flanking regions — as a linear vector.
Each genomic region in the blueprint becomes a node in the alignment data structure.
Repeats are represented as cycles: the repeat node has a self-loop edge.
Edges connect each feature to the next feature and to each downstream skippable feature.

### 2. Recruit reads

Reads overlapping the repeat region plus padding are recruited from the BAM.
When the catalog entry specifies off-target regions, reads are also recruited from those regions, along with unmapped reads whose weighted purity score against the repeat unit is high enough.

### 3. Generate the optimal local realignment

Recruited reads are aligned against the local realignment data structure.
Because the repeat node can be visited multiple times, each visit records its own alignment metrics.
The caller intersects each alignment path with the next to choose the best alignment for a read, falling back to the first path when no common path exists.

Key concepts:

- **LF** = Left Flank node, **R** = Repeat node, **RF** = Right Flank node
- Repeat nodes can be visited multiple times (for example, R → R → R for three repeat units)
- Each path represents a different hypothesis about the repeat count
- The optimal alignment is chosen by scoring and intersecting paths

### 4. Analyze the repeat site

Reads are classified as spanning, flanking, or in-repeat, and an alignment score is recorded for each plausible motif repeat count.
The caller computes the maximum-likelihood genotype and confidence interval that best explains the repeat counts observed across the aligned reads.

### 5. Experimental features: Variant calling and visualization

Use `--call-str-region-variants` to call SNVs and indels from the projected haplotype pileups within and around each STR locus.
The caller derives the projected haplotype contigs from the STR genotype and calls variants from reads anchored to a flank on at least one side.
It writes the individual calls to *str-region-variants.vcf.gz* and marks repeat-VCF records with `LocusSequenceVariant` when at least one SNV is called for that catalog locus.
The individual calls provide the sequence-level evidence to review when a locus has an interruption or other sequence variant.

Use `--add-pileup-visualization` to write an SVG pileup for visual review.
Use `--add-outputs-for-igv` to write the haplotype FASTA and projected-read BAM for interactive IGV review.
For an interactive review, load *str_haplotypes.fa* as the IGV genome and then load *str_projected_reads.bam*.
See [IGV visualization files](file-formats.md#igv-visualization-files) and [STR-region variant VCF output](file-formats.md#str-region-variant-vcf-output) for the output formats.

</details>

## Example command

### Minimal command

The Docker image includes the default hg38 catalog.

```shell
str_caller \
  --bam sample.bam \
  --reference hg38.fa
```

## Input

| Input | Required | Notes |
|:------|:---------|:------|
| BAM file | Yes | Passed with `--bam`. Coordinate-sorted, with an accompanying index (*.bai*) in the same directory. |
| Reference genome | Yes | Passed with `--reference`. Reference FASTA with a *.fai* index. |

> **Note:** SBX-S (simplex) input is not yet supported.
> If the input BAM is detected as SBX-S, or cannot be confirmed as duplex, the STR caller logs a warning and writes empty output files instead of calling repeats.
> The STR caller is designed for duplex (SBX-D / SBX-FAST) input.

## Output

Output is written to the directory given by `--output-dir`.

| Output | Description |
|:-------|:------------|
| *output.vcf.gz* (and *output.vcf.gz.tbi*) | One VCF record per repeat with location and genotype. Custom symbolic alleles and fields are described in [STR VCF output](file-formats.md#str-vcf-output). |
| *output.json* | Sample parameters and per-locus analysis results. Described in [STR JSON output](file-formats.md#str-json-output). |
| *metrics/metrics.tsv* | Run-level summary metrics. Described in [Metrics](metrics.md#summary-statistics). |
| *output.bam* | Optional BAM of realigned reads around the analyzed repeat regions. Enable with `--output-bam`. |
| *pileups/<locus>.svg* | Optional experimental feature to generate haplotype-aware read-pileup visualization for each selected locus. Enable with `--add-pileup-visualization`. Described in [Pileup visualization SVG output](file-formats.md#pileup-visualization-svg-output). |
| *pileups/str_haplotypes.fa*, *str_haplotypes.fa.fai*, *str_projected_reads.bam*, and *str_projected_reads.bam.bai* | Optional experimental feature that generates a projected-haplotype reference and coordinate-sorted projected reads for IGV inspection. Enable with `--add-outputs-for-igv`. |
| *str-region-variants.vcf.gz* (and *.tbi*) | Optional experimental feature that calls small variants inside and around STR loci. Enable with `--call-str-region-variants`. Described in [STR-region variant VCF output](file-formats.md#str-region-variant-vcf-output). |

The caller writes a pileup SVG for every locus with visualization data when `--add-pileup-visualization` is set.
Use `--pileup-visualization-loci` to limit SVG files to selected catalog loci without changing the STR calls or other outputs.
`--add-outputs-for-igv` writes the FASTA, BAM, and indexes without writing SVG files unless `--add-pileup-visualization` is also set.
Note: The experimental features for detecting and visualizing interruptions/variants, `--add-pileup-visualization`, `--add-outputs-for-igv`, and `--call-str-region-variants`, have gone through limited testing, and results can be less accurate for untested data.
See [Experimental features and limitations](appendix.md#experimental-features-and-limitations) for interpretation guidance.

## Tips for analysis

Expand the sections below for guidance on interpreting output and tuning the run.

<details>

<summary>Assessing call reliability</summary>

Review the genotype confidence intervals and supporting evidence (`REPCN`/`REPCI`, and the read-support fields `LC`, `ADSP`, `ADFL`, `ADIR`) to judge each call.
Exclude loci flagged `LowDepth` (fewer than 10 supporting reads, not currently configurable) to avoid unreliable calls, and compare `<STRn>` symbolic alleles to published repeat-count thresholds for research interpretation.

</details>

<details>

<summary>Interpreting symbolic alleles</summary>

A `<STRn>` allele denotes `n` repeat units. Use the `REF` and `RU` INFO fields to convert repeat counts to sequence length and to compare against the reference.
For repeats expanding beyond read length, rely on the reported confidence intervals rather than the point estimate.

</details>

<details>

<summary>Cross-referencing and batching</summary>

Map `VariantId` or `LocusId` back to the catalog to interpret biological relevance, and inspect read alignments in a genome browser at loci of interest.
Aggregating results across samples helps identify outliers and population-level patterns.

</details>

## CLI options

Required parameters are shown in **bold**.
This command also accepts the common `-h, --help`, `-v, --version`, and `-l, --log-level` options.

### Input options

| Parameter | Description | Value(s) |
|:----------|:-------------|:---------|
| **--bam** | Path to the BAM file and associated index for the sample to be analyzed. The BAM file must be coordinate-sorted and indexed. | Path to a sorted and indexed BAM file and associated `*.bai` file. |
| **--reference** | Path to the reference FASTA file used to compare the reads against. The FASTA file must be indexed with a `.fai` file. | Path to an indexed reference FASTA file |
| --legacy-duplex-bam | If set, duplex reads in the BAM file are encoded based on raw read quality scores for legacy duplex and the YC tag is ignored and thus R1 and R2 sequences are also not generated. | Flag (no value; i.e., false) |

### Output options

| Parameter | Description | Value(s) |
|:----------|:-------------|:---------|
| --output-dir | Path to the output directory where the metrics will be written to. | Path to a directory. [default: `.`] |
| --output-json | If set, a JSON file summarizing analysis results will be written to the output directory (e.g., <output_dir>/output.json). | Filename under output directory. [default: `output.json`] |
| --output-metrics | Custom filename for the metrics TSV file (written to `<output-dir>/metrics/`). | Filename under the metrics subdirectory. [default: `metrics.tsv`] |
| --output-vcf | Sets a custom VCF file name summarizing analysis results (e.g., <output_dir>/output.vcf). | Filename under output directory. [default: `output.vcf.gz` and `output.vcf.gz.tbi`] |
| --output-bam | If set, a realigned BAM file around the repeat regions analyzed will be written to the output directory (e.g., <output_dir>/output.bam). | Filename under output directory. [default: `none`] |
| --overwrite | Overwrite existing VCF/JSON output files instead of failing. | Flag (no value) [default: `false`] |
| --skip-sex-output | If set, the sex output is excluded from all output files. | Flag (no value) [default: `false`] |
| --write-reads-as-unmapped | If set, reads in output bam will be unmapped. | Flag (no value) [default: `false`] |
| --add-pileup-visualization | Writes one haplotype-aware read-pileup SVG per locus under `<output-dir>/pileups/`. See [Pileup visualization SVG output](file-formats.md#pileup-visualization-svg-output). | Flag (no value) [default: `false`] |
| --pileup-visualization-loci | Limits pileup SVG output to comma-separated exact locus IDs from the STR catalog. Requires `--add-pileup-visualization`; it does not change STR calling or other outputs. Whitespace is ignored, duplicates are ignored, and every requested ID must occur in the catalog. | Comma-separated locus IDs, for example, `ATXN3,FMR1` |
| --add-outputs-for-igv | Writes *str_haplotypes.fa*, its *.fai* index, *str_projected_reads.bam*, and its *.bai* index under `<output-dir>/pileups/` for IGV inspection. | Flag (no value) [default: `false`] |
| --call-str-region-variants | Experimental: Writes `<output-dir>/str-region-variants.vcf.gz` and its tabix index. The main STR VCF also includes the `LocusSequenceVariant` FORMAT field. See [STR-region variant VCF output](file-formats.md#str-region-variant-vcf-output). | Flag (no value) [default: `false`] |

When `--add-outputs-for-igv` is set, the caller writes `str_haplotypes.fa`, `str_haplotypes.fa.fai`, `str_projected_reads.bam`, and `str_projected_reads.bam.bai` under `<output-dir>/pileups/`.
When `--call-str-region-variants` is set, the caller writes the experimental `<output-dir>/str-region-variants.vcf.gz` file and its tabix index, with SNV and indel calls for the STR region and flanks.
For each STR locus with an SNV, the main repeat-count VCF includes the `LocusSequenceVariant` FORMAT field.

### Sample metadata options

| Parameter | Description | Value(s) |
|:----------|:-------------|:---------|
| --sample-name | Name of the BAM sample.  If not provided it will be inferred from the bam file or marked as unknown. | String sample name description |
| --sex | How to obtain the BAM sample sex for genotype priors. `auto-detect` infers it from the BAM and is the default; `male` or `female` overrides prediction; `skip-prediction` leaves the resolved output sex as `unknown`. | `auto-detect`, `male`, `female`, or `skip-prediction` [default: `auto-detect`] |

### Read filtering options

| Parameter | Description | Value(s) |
|:----------|:-------------|:---------|
| --region-padding | Padding for pulling additional bases on either side of the reference region as well as adjusting the flanking region when building the local realignment.  NOTE: using the non-default padding is not recommended as it can affect STR calling accuracy at certain sites. | Integer between 0 and 100000 [default: `800`] |
| --min-base-type | Minimum base type of duplex reads supporting a base for it to be considered in the analysis (only applies to duplex reads with YC tags). Bases decoded below this threshold are treated as low quality. | `concordant`, `simplex`, `discordant` [default: `concordant`] |
| --min-weighted-purity-score | Weighted purity score cutoff for off-target read recruitment; reads scoring below this threshold are excluded from recruitment | Float between 0 and 1 [default: `0.9`] |
| --skip-r1-r2-decoding | If set, skip generating additional R1 and R2 reads from the YC tag. This can decrease memory usage but can also reduce accuracy in certain scenarios. | Flag (no value; i.e., false) |
| --include-discordant-bases | If set, include discordant bases in the consensus read (Note: by default discordant bases are converted to N in the consensus read; this flag only applies to consensus reads and does not affect decoded R1/R2 reads). | Flag (no value; i.e., false) |

### Genotyping options

| Parameter | Description                                                                                                                                                                                                                                                                       | Value(s) |
|:----------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------|
| --aligner | Aligner used to produce the input BAM. This affects how reads are interpreted and processed during analysis. When omitted, or set to `unknown`, the aligner is inferred from the BAM `@PG` header lines, defaulting to `bwa` if no known aligner is found. | One of `unknown`, `bwa`, `giraffe` [default: `unknown` (auto-detected from the BAM `@PG` header)] |
| --min-spanning-reads-hom-fraction | For RepeatFXN variant type, sites with a spanning read fraction below this threshold are used to identify HOM expansions.                                                                                                                                                         | Float between 0 and 1 [default: `0.25`] |
| --max-spanning-reads-low-confidence-fraction | For RepeatFXN variant type, sites with a spanning read fraction at or above --min-spanning-reads-hom-fraction but below this threshold are treated as low confidence (i.e possible expansion). For RareRepeatFXN this threshold is used to determine off-target read recruitment. | Float between 0 and 1 [default: `0.35`] |

### Performance options

| Parameter | Description | Value(s) |
|:----------|:-------------|:---------|
| --threads | Number of threads to use for parallel processing.  0 for all available threads. | Integer greater or equal to 0. [default: `1`] |
