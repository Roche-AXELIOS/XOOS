# Measuring Error Rate for SBX Duplex Data

## Overview

This tutorial demonstrates how to measure empirical error rates in SBX duplex sequencing data using the `best` tool.
The workflow compares sequencing results against NIST Genome in a Bottle truth data to calculate accurate quality scores and error metrics.

{% stepper %}
{% step %}

## Prerequisites

- [samtools](https://github.com/samtools/samtools) installed
- [BWA-MEM](https://github.com/lh3/bwa) installed
- `best` installed. **Note:** Pre-compiled versions of `best` are only available for x86_64 Linux, for other operating systems or architectures please see [the best README](https://github.com/google/best/?tab=readme-ov-file#installing).
- The HG001 BAM downloaded and accessible, data can be downloaded from the dataset titled [091025 Webinar GIAB BAMs Giraffe](https://web.sbxdata.kamino.platform.navify.com/files/)

## Download Additional Prerequisites

**Reference Data for HG001 (NA12878) sample:**

- `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz` - GRCh38 human reference genome for read alignment
- `HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` - NIST Genome in a Bottle truth set containing known variants for the NA12878 reference sample
- `HG001_GRCh38_1_22_v4.2.1_benchmark.bed` - High-confidence regions where the truth set variants are reliable for benchmarking

These files enable comparison of sequencing results against established truth data to accurately measure error rates.

```shell
curl -OL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
gzip -d GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz

curl -OL https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl -OL https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
```

{% endstep %}
{% step %}

## Convert BAM to FASTQ and Align

This step converts the BAM file back to FASTQ before aligning by:

1. Extracting reads from chromosome 10 only (for faster processing)
2. Converting BAM to FASTQ format while preserving YC tags (duplex consensus information)
3. Aligning reads to the reference genome using BWA-MEM with `-C` flag to preserve FASTQ comments
4. Filtering alignments with mapping quality ≥4 to remove poor quality alignments
5. Sorting the resulting BAM file and creating an index for downstream analysis

```shell
samtools view -h ${bam} chr10 \
  | samtools bam2fq -@ ${threads} -T YC \
  | bwa mem -C -t ${threads} GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna - \
  | samtools view  -q 4 \
  | samtools sort -@ ${threads} --write-index -o HG001.bwa.bam##idx##HG001.bwa.bam.bai
```

{% endstep %}
{% step %}

## Run best

This step calculates empirical quality scores and error rates by:

1. Creating a refined high-confidence region BED file by subtracting known variant positions from the benchmark regions (avoiding areas where differences are expected)
2. Running the `best` tool to analyze the aligned BAM file against the reference genome within high-confidence intervals
3. Generating quality score statistics comparing reported quality scores with empirical error rates
4. Outputting results showing both overall alignment statistics and high-confidence region-specific metrics

```shell
bedtools subtract -a HG001_GRCh38_1_22_v4.2.1_benchmark.bed -b HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  | awk '{ print $0 "\thigh_confidence" }' \
  > HG001_high_confidence.bed

best --bed-intervals HG001_high_confidence.bed -t ${threads} HG001.bwa.bam GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna HG001

cat HG001.summary_qual_score_stats.csv
feature,qual_score,empirical_qv
all_alignments,5,5.75
all_alignments,22,23.04
all_alignments,39,28.10
high_confidence,5,5.75
high_confidence,22,24.56
high_confidence,39,39.84
```

{% endstep %}
{% endstepper %}
