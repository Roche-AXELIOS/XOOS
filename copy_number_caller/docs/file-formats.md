# File formats

This page is the reference for the tool-specific formats that Copy Number Caller produces.
For Copy Number Caller, the primary deliverables are the segment ([SEG](#segments)) and [VCF](#vcf) callsets; the remaining formats are intermediate coverage, interval, log-ratio, and BAF tracks plus purity/ploidy analysis outputs.

Coordinate conventions differ by file type: BED-based outputs use 0-based, half-open coordinates, while SEG outputs use 1-based, inclusive coordinates. Each section below states the convention it uses.
Absent values in BED, SEG, and TSV outputs are written as `NA`.

Each section lists the commands it **applies to** (linked to their usage pages) and, at the end, links back to the command output where the format is produced.

## Coverage

Applies to: [`germline-wgs`](usage-germline-wgs.md), [`tumor-normal-wgs`](usage-tumor-normal-wgs.md), [`tumor-only-wgs`](usage-tumor-only-wgs.md).

Coverage data is provided in BED format (tab-delimited, with header, 0-based half-open coordinates):

Raw coverage files (`coverage.bed`) report whole-number `Counts` and `TotalCoverage` values.
GC-corrected coverage files (`corrected_coverage.bed`, `normal_corrected_coverage.bed`, `tumor_corrected_coverage.bed`) report the same two columns after per-region GC-bias scaling, so their values are fractional and are not rounded to whole numbers.

| Column | Definition |
|:---------|:-----------------------|
| Contig   | The chromosome or contig name. |
| Start   | The 0-based start position of the interval (inclusive). |
| End   | The end position of the interval (exclusive). |
| Counts   | The read counts for the interval: whole numbers in `coverage.bed`, fractional in GC-corrected coverage files. Uses `NA` when a value is absent. |
| TotalCoverage   | The normalized or raw coverage value: whole numbers in `coverage.bed`, fractional in GC-corrected coverage files. Uses `NA` when a value is absent. |
| OnTarget   | Indicates if the region is within the intended target set. `1` is true and `0` is false. |

File example:

```text
#Contig Start End Counts TotalCoverage OnTarget
chr1    0       500     540     5400    1
chr1    500     1000    469     4690    1
chrX    0       500     500     5000    1
chrX    500     1000    500     5000    1
chrY    0       500     500     5000    1
chrY    500     1000    500     5000    1
```

Back to [`germline-wgs` output](usage-germline-wgs.md#output),
[`tumor-normal-wgs` output](usage-tumor-normal-wgs.md#output), and
[`tumor-only-wgs` output](usage-tumor-only-wgs.md#output).

## Intervals

Applies to: [`germline-wgs`](usage-germline-wgs.md), [`tumor-normal-wgs`](usage-tumor-normal-wgs.md), [`tumor-only-wgs`](usage-tumor-only-wgs.md).

Augmented intervals data is provided in BED format (tab-delimited, with header, 0-based half-open coordinates):

| Column | Definition |
|:---------|:-----------------------|
| Contig   | The chromosome or contig name. |
| Start   | The 0-based start position of the interval (inclusive). |
| End   | The end position of the interval (exclusive). |
| GCBias   | The fraction of G/C bases in the interval (range `0.000-1.000`). |
| Mappability   | The uniqueness score of the interval (range `0.000-1.000`). |
| OnTarget   | Indicates if the region is within the intended target set. |

File example:

```text
#Contig Start End GCBias Mappability OnTarget
chr1    0       500     0.510000        1.000000        true
chr1    500     1000    0.476000        1.000000        true
chrX    0       500     0.508000        1.000000        true
chrX    500     1000    0.498000        1.000000        true
chrY    0       500     0.490000        1.000000        true
```

Back to [`germline-wgs` output](usage-germline-wgs.md#output),
[`tumor-normal-wgs` output](usage-tumor-normal-wgs.md#output), and
[`tumor-only-wgs` output](usage-tumor-only-wgs.md#output).

## Observations

Applies to: [`germline-wgs`](usage-germline-wgs.md), [`tumor-normal-wgs`](usage-tumor-normal-wgs.md), [`tumor-only-wgs`](usage-tumor-only-wgs.md).

Mapping quality data is provided in BED format (tab-delimited, with header, 0-based half-open coordinates):

| Column | Definition |
|:---------|:-----------------------|
| Contig   | The chromosome or contig name. |
| Start   | The 0-based start position of the interval (inclusive). |
| End   | The end position of the interval (exclusive). |
| MeanMapQ   | The mean of the average mapping quality (MAPQ) scores of all targets contained within the segment. Uses `NA` when a value is absent. |

File example:

```text
#Contig Start End MeanMapQ
chr1    0       500     60.000000
chr1    500     1000    60.000000
chrX    0       500     60.000000
chrX    500     1000    60.000000
chrY    0       500     60.000000
chrY    500     1000    60.000000
```

Back to [`germline-wgs` output](usage-germline-wgs.md#output),
[`tumor-normal-wgs` output](usage-tumor-normal-wgs.md#output), and
[`tumor-only-wgs` output](usage-tumor-only-wgs.md#output).

## Log ratios

Applies to: [`germline-wgs`](usage-germline-wgs.md), [`tumor-normal-wgs`](usage-tumor-normal-wgs.md), [`tumor-only-wgs`](usage-tumor-only-wgs.md), and as input to [`predict-somatic-cnv`](usage-predict-somatic-cnv.md).

Log ratios are provided in BED format (tab-delimited, with header, 0-based half-open coordinates). The same values are also emitted in BigWig format (`tracks/log_ratios.bw`) for visualization.

| Column | Definition |
|:---------|:-----------------------|
| Contig   | The chromosome or contig name. |
| Start   | The 0-based start position of the interval (inclusive). |
| End   | The end position of the interval (exclusive). |
| LogR   | The Log2 of the read depth ratio for the interval. |

The `LogR` value is computed differently per command:

- `germline-wgs`: Log2 of the read depth ratio between the interval and the median autosome interval read depth.
- `tumor-normal-wgs`: Log2 of the tumor/normal read depth ratio.
- `tumor-only-wgs`: Log2 of the self-normalized read depth ratio (each interval against the median autosomal coverage).

File example:

```text
#Contig Start End LogR
chr1    0       10      0.000000
chr1    10      20      -0.016639
chr1    20      30      0.008143
chr1    30      40      -0.016639
...
```

Back to [`germline-wgs` output](usage-germline-wgs.md#output),
[`tumor-normal-wgs` output](usage-tumor-normal-wgs.md#output),
[`tumor-only-wgs` output](usage-tumor-only-wgs.md#output), and
[`predict-somatic-cnv` input](usage-predict-somatic-cnv.md#input).

## BAFs

Applies to: [`germline-wgs`](usage-germline-wgs.md) (only when `--vcf` is provided), [`tumor-normal-wgs`](usage-tumor-normal-wgs.md), [`tumor-only-wgs`](usage-tumor-only-wgs.md), and as input to [`predict-somatic-cnv`](usage-predict-somatic-cnv.md).

B-allele fractions (BAFs) are provided in BED format (tab-delimited, with header, 0-based half-open coordinates). The BAF values are also emitted in BigWig format (`tracks/bafs.bw`) for visualization.

| Column | Definition |
|:---------|:-----------------------|
| Contig   | The chromosome or contig name in which the current record is located. |
| Start   | The 0-based start position (inclusive). |
| End   | The end position (exclusive). For single-base variants, End = Start + 1. |
| BAF   | The B-allele fraction, calculated as `Alt_AD / (Ref_AD + Alt_AD)` in the inclusive range [0, 1]. |
| Ref_AD   | The read depth over the reference allele. |
| Alt_AD   | The read depth over the alternate allele. |

File example:

```text
#Contig Start End BAF Ref_AD Alt_AD
chr1    0       1       0.455556    49      41
chr1    13      14      0.488889    46      44
chr1    21      22      0.533333    42      48
chr1    39      40      0.455556    49      41
...
```

Back to [`germline-wgs` output](usage-germline-wgs.md#output),
[`tumor-normal-wgs` output](usage-tumor-normal-wgs.md#output),
[`tumor-only-wgs` output](usage-tumor-only-wgs.md#output), and
[`predict-somatic-cnv` input](usage-predict-somatic-cnv.md#input).

## Segments

Applies to: [`germline-wgs`](usage-germline-wgs.md), [`tumor-normal-wgs`](usage-tumor-normal-wgs.md), [`tumor-only-wgs`](usage-tumor-only-wgs.md), and [`predict-somatic-cnv`](usage-predict-somatic-cnv.md) (both input and output).

All `.seg` files follow the standard SEG format. Coordinates are 1-based, inclusive.

The first few lines of the TSV file contain a `##RocheCommandLine=<...>` metadata header followed by the column description header:

```text
Id      Contig  Start   End
```

| Column | Definition |
|:---------|:-----------------------|
| Id   | Sample ID.  |
| Contig  | The chromosome or contig name.  |
| Start  | 1-based start position (inclusive). |
| End  |  1-based end position (inclusive). |

Additional columns can also be provided:

| Column | Application | Definition |
|:---------|:-----------------------|:-----------------------|
| AvgMeanMapq | Germline | Mean of the average mapping quality (MAPQ) scores of all targets contained within the segment. High values indicate unique genomic regions, while low values suggest repetitive regions where copy number calls can be less reliable. |
| ExpectedTotalCopyNumber | All | Estimated raw total copy number from the segment's mean log-ratio. |
| NumLogRatio | All | Number of log-ratio observations contained in the segment. |
| NumSnp | Somatic | Number of SNP BAF observations contained in the segment. |
| TotalCopyNumber | All | Estimated final total copy number for the segment. |
| MajorCopyNumber | Somatic | Estimated final copy number of the major allele. |
| MinorCopyNumber | Somatic | Estimated final copy number of the minor allele. |
| LogRatioLikelihood | All | Likelihood of observed log-ratios in the segment given purity and ploidy, if applicable, and total copy number. |
| BAlleleFrequencyLikelihood | Somatic | Likelihood of observed BAF values in the segment given purity, total copy number, and minor copy number. |
| JointLikelihood | Somatic | Combined statistical likelihood (sum of LogRatioLikelihood and BAlleleFrequencyLikelihood). |
| Purity | Somatic | Estimated tumor purity used for likelihood calculations. |
| Ploidy | Somatic | Estimated tumor ploidy used for likelihood calculations. |
| Sex | All | Predicted biological sex of the sample. One of `male`, `female`, or `unknown`. |
| Filter | Germline | One of three possible filter flags, lowMmMapq, nonIntegerTcn, or cnvLength. See below for more details. |
| IsInAllosome | All | True if the segment is on a sex chromosome (e.g., X or Y); false otherwise. |
| IsInPseudoAutosomalRegion | All | True if the segment is within a pseudo-autosomal region (PAR); false otherwise. |
| MirroredBaf | Somatic | Mirrored B-Allele Fraction for the segment (between 0.0 and 0.5 inclusive, or `NA`). |
| MeanDh | Somatic | Mean Degree of Heterozygosity (DH) derived from BAF observations in the segment. |
| IsGermlineCnv | Tumor-only (`tumor-only-wgs` and `predict-somatic-cnv --tumor-only`) | True when the segment fits the germline log-ratio likelihood model better than the somatic model; otherwise false. |
| MeanLogRatio | All | Mean of all log-ratio observations in the segment. Note: If present, this will always appear as the last column. |

### Segment filters

Segment filters are only included in germline application output (i.e. `germline-wgs`).

| Filter | Definition |
|:---------|:-----------------------|
| lowMmMapq | Low mapping quality. The segment's AvgMeanMapq falls below the required threshold, indicating low confidence. |
| nonIntegerTcn | Non-integer copy number. For SBX-S, the accepted expected total copy number is 0.8–1.2 for heterozygous deletions. Duplications require at least 2.65 in diploid regions or 1.65 in haploid regions, such as a male non-PAR allosome. Other or unresolved library types use 0.7–1.3 for heterozygous deletions and do not apply the duplication threshold. Homozygous deletions retain a maximum of 0.3 for all library types. |
| cnvLength | Small variant length. The segment CNV length is below the minimum detection threshold (<2000bp). |

Back to [`germline-wgs` output](usage-germline-wgs.md#output),
[`tumor-normal-wgs` output](usage-tumor-normal-wgs.md#output),
[`tumor-only-wgs` output](usage-tumor-only-wgs.md#output), and
[`predict-somatic-cnv` output](usage-predict-somatic-cnv.md#output).

## Purity ploidy grid

Applies to: [`tumor-normal-wgs`](usage-tumor-normal-wgs.md), [`tumor-only-wgs`](usage-tumor-only-wgs.md). Not generated when `--tumor-purity` and `--tumor-ploidy` are both provided.

The purity ploidy grid file describes the purity-ploidy solutions considered for a tumor sample and their associated likelihoods.
The file is in a tab-delimited format.

The first line of the TSV file is the following header:

```text
Purity  Ploidy  TotalJointLikelihood    TotalLogRatioLikelihood TotalBAlleleFrequencyLikelihood
```

`TotalJointLikelihood` is the sum of `TotalLogRatioLikelihood` and `TotalBAlleleFrequencyLikelihood` (if they are not both `NA`).

| Column | Definition |
|:---------|:-----------------------|
| Purity   | The estimated fraction of the sample composed of tumor cells (range `0.15-1.00`, inclusive). |
| Ploidy   | The estimated average copy number across the tumor genome. |
| TotalJointLikelihood   | The combined statistical confidence score for a specific purity/ploidy pair, calculated as the sum of Log-Ratio and BAF likelihoods. |
| TotalLogRatioLikelihood   | The likelihood score based on how well the purity/ploidy model fits the observed total read depth (copy number intensity) data. |
| TotalBAlleleFrequencyLikelihood   | The likelihood score based on how well the purity/ploidy model fits the observed allelic imbalance (minor allele proportions). |

File example:

```text
Purity  Ploidy  TotalJointLikelihood    TotalLogRatioLikelihood TotalBAlleleFrequencyLikelihood
0.15    1.4     -38.938400      -42.861200      3.922820
0.15    1.5     -38.567700      -42.751800      4.184150
0.15    1.6     -34.050500      -38.747500      4.696950
...
```

The best solution is typically the one with the highest `TotalJointLikelihood`.

Back to [`tumor-normal-wgs` output](usage-tumor-normal-wgs.md#output) and
[`tumor-only-wgs` output](usage-tumor-only-wgs.md#output).

## Segment likelihood

Applies to: [`tumor-normal-wgs`](usage-tumor-normal-wgs.md), [`tumor-only-wgs`](usage-tumor-only-wgs.md). Written only when the hidden `--write-segment-likelihood` flag is set.

The segment likelihood files are lossless purity/ploidy analysis outputs emitted from the search grid.
One file is written per search pass:

- *supplemental/segment_likelihood.diploid_bias.tsv* — the [diploid-bias prior pass](troubleshooting.md#purity-and-ploidy-ambiguity-in-whole-genome-doubled-tumors).
- *supplemental/segment_likelihood.uniform.tsv* — the [uniform prior pass](troubleshooting.md#purity-and-ploidy-ambiguity-in-whole-genome-doubled-tumors), written only when the [diploid-bias solution](troubleshooting.md#purity-and-ploidy-ambiguity-in-whole-genome-doubled-tumors) is rejected by whole-genome-doubling (WGD) QC and the search is re-run.

Each file contains one row per contributing segment per purity/ploidy grid cell, so the per-segment
contributions that are aggregated into the `supplemental/purity_ploidy_grid.tsv` totals can be inspected directly.
The grid axes are purity from `0.15` to `0.95` (step `0.01`) and ploidy from `1.4` to `6.0` (step `0.1`).
The file is in a tab-delimited format with metadata lines (prefixed with `##`) describing the run, including
`##TotalCopyNumberPrior=<prior>` identifying the pass.

The analytical layout depends on the run mode. Tumor/normal runs emit somatic columns only.
In tumor-only runs, additional germline/classification/contribution columns are included; they are populated only on the [diploid-bias pass](appendix.md#whole-genome-doubling) and are `NA` on the [uniform pass](appendix.md#whole-genome-doubling).
Allosome segments are omitted in tumor-only mode (matching the purity/ploidy search).

| Column | Definition |
|:---------|:-----------------------|
| Contig   | The contig of the segment. |
| Start   | The 0-based start coordinate of the segment. |
| End   | The end coordinate of the segment. |
| Purity   | The purity grid value for this row. |
| Ploidy   | The ploidy grid value for this row. |
| SomaticLogRLikelihood   | The segment's log-ratio likelihood under the somatic model at this purity/ploidy. |
| SomaticBafLikelihood   | The segment's BAF likelihood under the somatic model at this purity/ploidy. `NA` when BAF is not used. |
| SomaticJointLikelihood   | The combined somatic likelihood (sum of the log-ratio and BAF likelihoods). Equals `SomaticLogRLikelihood` when BAF is not used (i.e. when `SomaticBafLikelihood` is `NA`), since only the log-ratio is available to the search. |
| SomaticTotalCopyNumber   | The maximum-likelihood somatic total copy number for the segment at this purity/ploidy. |
| SomaticMinorCopyNumber   | The maximum-likelihood somatic minor copy number for the segment at this purity/ploidy. `NA` when BAF is not used. |
| GermlineLogRLikelihood   | Present only in tumor-only analysis output; `NA` on the uniform pass. |
| GermlineTotalCopyNumber   | Present only in tumor-only analysis output; `NA` on the uniform pass. |
| Classification   | Present only in tumor-only analysis output; `NA` on the uniform pass. |
| ContributionLikelihood   | Present only in tumor-only analysis output; `NA` on the uniform pass. |

Absent values are reported as `NA`.

Back to [`tumor-normal-wgs` output](usage-tumor-normal-wgs.md#output) and
[`tumor-only-wgs` output](usage-tumor-only-wgs.md#output).

## VCF

Applies to: [`germline-wgs`](usage-germline-wgs.md), [`tumor-normal-wgs`](usage-tumor-normal-wgs.md), [`tumor-only-wgs`](usage-tumor-only-wgs.md), and [`predict-somatic-cnv`](usage-predict-somatic-cnv.md).

The VCF output follows the standard [VCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf), with a number of additional or modified fields.
Use the tables below to interpret the variant calls and metadata.

### Header

The predicted sex of the sample is recorded in the header.

- Format: \#\#Sex=\<male|female|unknown\>
- Values: `female` (Female/XX), `male` (Male/XY), `unknown` (Unknown).

### ALT

| ID | Description |
|:----|:----|
| CNV | Copy number variant region |
| DEL | Deletion (lower copy number relative to the reference or a deletion breakpoint) |
| DUP | Duplication (elevated copy number relative to the reference or a tandem duplication breakpoint) |

### INFO

| ID | Description |
|:----|:----|
| REFLEN | Number of REF positions included in the record |
| SVLEN | Difference in length between REF and ALT alleles |
| SVTYPE | Type of structural variant |
| END | End position of the variant |
| GERMLINE | Tumor-only flag, emitted when the segment fits the germline copy-number model better than the somatic model. |

### FILTER

| ID | Applications | Description |
|:----|:----|:----|
| cnvLength | Germline | CNV with length below 2000 bp |
| lowMmMapq | Germline | Low mean of mean MAPQ for overlapping targets |
| nonIntegerTcn | Germline | Copy-number call where observed LR deviates beyond the library-type-specific expected total copy number range |

### FORMAT

| ID | Applications | Description |
|:----|:----|:----|
| GT | All | Genotype |
| CN | All | Estimated total copy number |
| LR | All | Log2 of the Mean Read Depth Ratio (MRDR) |
| MRDR | All | Mean of the read depth ratios of overlapping targets |
| NUMTARGETS | All | Number of overlapping targets |
| MMMAPQ | All | Mean of mean MAPQ of overlapping targets |
| MCN | Somatic | Estimated minor copy number |
| MBAF | Somatic | Mean mirrored B-allele fraction of overlapping SNPs. SNPs with BAF > 0.5 are mirrored (values range 0.5 to 1.0). |
| NUMSNPS | Somatic | Number of overlapping SNPs |

Back to [`germline-wgs` output](usage-germline-wgs.md#output),
[`tumor-normal-wgs` output](usage-tumor-normal-wgs.md#output),
[`tumor-only-wgs` output](usage-tumor-only-wgs.md#output), and
[`predict-somatic-cnv` output](usage-predict-somatic-cnv.md#output).
