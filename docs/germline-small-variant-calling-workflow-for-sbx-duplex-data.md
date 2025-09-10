# Germline Small Variant Calling Workflow for SBX Duplex Data

## Overview

This tutorial reproduces the steps demonstrated in the ["Germline Small Variant Calling Workflow for SBX Duplex Data" webinar](http://sequencing.roche.com/sbxdata), covering three key analysis workflows for germline genome sequencing data.
Each section builds upon Docker-based tools to ensure reproducible results.

{% stepper %}
{% step %}

## Prerequisites

- [Docker](../overview/docker-guide) installed
- [samtools](https://github.com/samtools/samtools) installed
- [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools) installed
- GATK Docker image [broadinstitute/gatk-nightly:2025-08-19-4.6.2.0-17-g2a1f41bf3-NIGHTLY-SNAPSHOT](https://hub.docker.com/r/broadinstitute/gatk-nightly/tags?name=2025-08-19-4.6.2.0-17) accessible
- XOOS Small Variant Caller Docker image accessible
- The HG001 BAM downloaded and accessible, data can be downloaded from the dataset titled [091025 Webinar GIAB BAMs Giraffe](https://web.sbxdata.kamino.platform.navify.com/files/)

{% hint style="info" %}
**Some commands** do use Docker, but to keep the commands short and readable the boilerplate for the actual `docker run` is removed, please see [the Docker guide](../overview/docker-guide) for more information about how to use Docker for data analysis.
{% endhint %}

{% endstep %}
{% step %}

## Downsample BAM File

The BAM files shared during the webinar are full coverage, but the current pre-trained multi-sample models were trained on HG002-HG007 30x data with HG001 left out for evaluation. Therefore, we will analyze the HG001 BAM.
Additionally, when computing coverage, we consider only the concordant duplex bases.
We have provided the correct subsampling ratios to achieve 30x for each shared BAM.

| Sample | Downsampling Ratio |
|--------|--------------------|
| HG001  | 0.69               |
| HG002  | 0.63               |
| HG003  | 0.69               |
| HG004  | 0.68               |
| HG005  | 0.63               |
| HG006  | 0.68               |
| HG007  | 0.68               |

The following samtools command will subsample the HG001 BAM to 30x.

```shell
samtools view \
  -@ ${threads} \
  -s 0.69 \
  --subsample-seed 1234 \
  -b \
  --write-index \
  -o HG001.30x.bam##idx##HG001.30x.bam.bai \
  HG001.bam
```

{% endstep %}
{% step %}

## Call Variants with GATK HaplotypeCaller

This step performs variant calling using GATK HaplotypeCaller with optimized parameters for duplex sequencing data.
The following must be executed using the [broadinstitute/gatk-nightly:2025-08-19-4.6.2.0-17-g2a1f41bf3-NIGHTLY-SNAPSHOT](https://hub.docker.com/r/broadinstitute/gatk-nightly/tags?name=2025-08-19-4.6.2.0-17) version of the GATK Docker image.
From our experience, GATK HaplotypeCaller is not able to utilize more than 4 cores, because of this we recommend that in practice HaplotypeCaller is run in parallel on individual chromosomes to improve the turn-around time.

```shell
curl -OL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
gzip -d GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
curl -OL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai

gatk CreateSequenceDictionary -R GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

gatk HaplotypeCaller \
  -I HG001.30x.bam \
  -R GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
  -O HG001.30x.vcf.gz \
  -OVI \
  -bamout HG001.30x.bamout.bam \
  -OBI \
  -A AssemblyComplexity \
  -A TandemRepeat \
  -RF MappingQualityReadFilter \
  --activeregion-alt-multiplier 5 \
  --adaptive-pruning true \
  --enable-dynamic-read-disqualification-for-genotyping true \
  --mapping-quality-threshold-for-genotyping 1 \
  --minimum-mapping-quality 1 \
  --min-base-quality-score 6 \
  --native-pair-hmm-threads 4 \
  --smith-waterman FASTEST_AVAILABLE
```

{% endstep %}
{% step %}

## Filter Variants

This step runs the Roche Small Variant Caller on the GATK VCF leveraging additional information from the GATK BAM and gnomAD population allele frequency database.
Please see [the Roche Small Variant Caller Documentation](../small_variant_caller) for more detailed information.

```shell
curl -OL https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
curl -OL https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

filter_variants \
  --bam-input HG001.30x.bamout.bam \
  --vcf-input HG001.30x.vcf.gz \
  --pop-af-vcf af-only-gnomad.hg38.vcf.gz \
  --workflow germline \
  --threads ${threads} \
  --genome GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
  --model /resources/model-germline-sbxd-multisample-snv.txt.gz /resources/model-germline-sbxd-multisample-indel.txt.gz \
  --vcf-output HG001.30x.filtered.vcf.gz
```

{% endstep %}
{% step %}

## Evaluate Variants

This step compares the filtered variants against truth data to assess accuracy and performance metrics.
We use the `vcfeval` functionality of the RTG Tools suite, for more information please see [the RTG Tools README](https://github.com/RealTimeGenomics/rtg-tools/blob/master/README.md).

```shell
curl -OL https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl -OL https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
curl -OL https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed

rtg RTG_JAVA_OPTS="-Xmx2G" format -o GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.sdf GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

rtg RTG_JAVA_OPTS="-Xmx2G" vcfeval \
  --threads ${threads} \
  --template GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.sdf \
  --evaluation-regions HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
  --baseline HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --calls HG001.30x.vcf.gz \
  --decompose \
  --ref-overlap \
  --output HG001-vcfeval
```

For detailed documentation on the output of `vcfeval` refer to [the RTG Tools documentation](https://realtimegenomics.github.io/rtg-tools/rtg_command_reference.html#vcfeval-outputs).
To simply print the F1 scores for SNVs and indels, execute the following.

```shell
zcat HG001-vcfeval/snp_roc.tsv.gz | tail -1 | awk '{ print $8 }'

zcat HG001-vcfeval/non_snp_roc.tsv.gz | tail -1 | awk '{ print $8 }'
```

{% endstep %}
{% endstepper %}

### Final Output Files

| Filename | Description |
|----------|-------------|
| `HG001.30x.filtered.vcf.gz` | The VCF produced by the Small Variant Caller, [detailed documentation](../small_variant_caller#output-files-for-germline-variant-filtering) |
| `HG001-vcfeval` | The output of RTG Tools `vcfeval`, [detailed documentation](https://realtimegenomics.github.io/rtg-tools/rtg_command_reference.html#vcfeval-outputs) |
