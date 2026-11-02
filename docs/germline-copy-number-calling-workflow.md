# Germline Copy Number Calling Workflow for SBX Duplex

## Overview

This tutorial provides the steps to run the XOOS copy number caller on HG002 SBX data and how to evaluate the calls.

{% stepper %}
{% step %}

## Prerequisites

- [Docker](docker-guide.md) installed
- [CrossMap (0.7.0)](https://github.com/liguowang/CrossMap) installed
  - Newer versions of CrossMap have an issue where `INFO/END` is not lifted over correctly.
    It works in version 0.7.0.
- [bcftools](https://samtools.github.io/bcftools/bcftools.html) installed
- [Bgzip](https://www.htslib.org/doc/bgzip.html) installed
- [Tabix](https://www.htslib.org/doc/tabix.html) installed
- [witty.er](https://github.com/Illumina/witty.er) installed
- [jq](https://jqlang.org/) installed (optional)
  - Only used for formatting the output into a tabular view.
- XOOS Copy Number Caller Docker image accessible
- The HG002 BAM and BAI downloaded and accessible:
  - These files can be retrieved from <https://web.sbxdata.kamino.platform.navify.com/files>
- hg38 FASTA and FAI files. For instance:
  - `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz`
  - `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai`

{% hint style="info" %}
**Some commands** do use Docker, but to keep the commands short and readable the boilerplate for the actual `docker run` is removed, please see [the Docker guide](docker-guide.md) for more information about how to use Docker for data analysis.
{% endhint %}

{% endstep %}
{% step %}

## Call Germline Copy Number Variants

This step runs the XOOS copy number caller.

```shell
copy_number_caller germline-wgs \
        --bam HG002.bam \
        --reference GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
        --mappability-bigwig /resources/hg38_mappability.bw \
        --seed-segments /resources/hg38_seed.seg \
        --threads 16 \
        --output-dir cnvs/HG002
```

Note that `/resources/hg38_mappability.bw` and `/resources/hg38_seed.seg` are resource files in the Docker image.

In the `--output-dir`, you will find several output files with germline CNV calls being found in `germline_copy_number_callset.vcf.gz`

{% endstep %}
{% step %}

## Evaluate Variants

To evaluate the germline CNV calls, we use [witty.er](https://github.com/Illumina/witty.er).

The ground-truth CNV calls can be retrieved from [Zook et al. 2020. Nature Biotechnology](https://www.nature.com/articles/s41587-020-0538-8).

```bash
mkdir -p ground_truth

# SV calls
wget \
        --output-document "ground_truth/HG002_SVs_Tier1_v0.6.vcf.gz" \
        https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz

# The VCF has 3 duplicated INFO header lines.
# We remove them using awk and then bgzip and tabix it
zcat ground_truth/HG002_SVs_Tier1_v0.6.vcf.gz \
        | awk '!a[$0]++' \
        | bgzip --stdout \
        > ground_truth/HG002_SVs_Tier1_v0.6.no_dup_lines.vcf.gz

tabix ground_truth/HG002_SVs_Tier1_v0.6.no_dup_lines.vcf.gz

# High confidence CNV regions
wget \
        --output-document "ground_truth/HG002_SVs_Tier1_v0.6.bed" \
        ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed
```

These ground-truth files are in hg19. We lift them over to hg38 for evaluation:

```bash
# Download the chain file to enable liftover
wget \
        --output-document ground_truth/hg19ToHg38.over.chain.gz \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# Liftover high confidence BED to hg38
CrossMap bed \
        --chromid l \
        --unmap-file ground_truth/HG002_SVs_Tier1_v0.6.hg38.bed.unmap \
        ground_truth/hg19ToHg38.over.chain.gz \
        ground_truth/HG002_SVs_Tier1_v0.6.bed \
        ground_truth/HG002_SVs_Tier1_v0.6.hg38.bed

# Liftover SV VCF to hg38
CrossMap vcf \
        --chromid l \
        ground_truth/hg19ToHg38.over.chain.gz \
        ground_truth/HG002_SVs_Tier1_v0.6.no_dup_lines.vcf.gz \
        GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
        ground_truth/HG002_SVs_Tier1_v0.6.no_dup_lines.hg38.vcf

# Lifted over VCF requires sorting
bcftools sort \
        --output-type z \
        --output ground_truth/HG002_SVs_Tier1_v0.6.no_dup_lines.hg38.vcf.gz \
        ground_truth/HG002_SVs_Tier1_v0.6.no_dup_lines.hg38.vcf

tabix -p vcf ground_truth/HG002_SVs_Tier1_v0.6.no_dup_lines.hg38.vcf.gz
```

This will produce

1. `ground_truth/HG002_SVs_Tier1_v0.6.hg38.bed`: Genomic regions to evaluate on
2. `ground_truth/HG002_SVs_Tier1_v0.6.no_dup_lines.hg38.vcf.gz:` VCF that will act as the ground-truth to evaluate against.

Now run witty.er to evaluate the calls:

```bash
Wittyer \
        -i cnvs/HG002/germline_copy_number_callset.vcf.gz \
        -t ground_truth/HG002_SVs_Tier1_v0.6.no_dup_lines.hg38.vcf.gz \
        --configFile config_wittyer_CNV.json \
        -em CrossTypeAndSimpleCounting \
        -o cnvs/HG002/wittyer
```

For the `config_wittyer_CNV.json` file (based on <https://github.com/srbehera/DRAGEN_Analysis/blob/main/config_wittyer_CNV.json>), you can use:

```json
[
  {
    "variantType": "CopyNumberReference",
    "binSizes": "!1,1000,5000,10000,20000,50000",
    "bpDistance": 10000,
    "percentDistance": 0.25,
    "includedFilters": "PASS",
    "excludedFilters": "",
    "includeBed": "ground_truth/HG002_SVs_Tier1_v0.6.hg38.bed"
  },
  {
    "variantType": "Duplication",
    "binSizes": "!1,1000,5000,10000,20000,50000",
    "bpDistance": 10000,
    "percentDistance": 0.25,
    "includedFilters": "PASS",
    "excludedFilters": "",
    "includeBed": "ground_truth/HG002_SVs_Tier1_v0.6.hg38.bed"
  },
  {
    "variantType": "Deletion",
    "binSizes": "!1,1000,5000,10000,20000,50000",
    "bpDistance": 10000,
    "percentDistance": 0.25,
    "includedFilters": "PASS",
    "excludedFilters": "",
    "includeBed": "ground_truth/HG002_SVs_Tier1_v0.6.hg38.bed"
  }
]
```

Then you can use `jq` to parse the wittyer output to form a nice table summarizing the F1, recall, and precision statistics of deletions:

```bash
jq -r '
  [
    .PerSampleStats[0].DetailedStats[]
    | select(.VariantType == "CopyNumberLoss" or .VariantType == "Deletion")
    | .PerBinStats[]
    | . as $bin
    | (.Stats[] | select(.StatsType == "Event"))
    | {
        bin:      $bin.Bin,
        truth_tp: .TruthTpCount,
        truth_fn: .TruthFnCount,
        query_tp: .QueryTpCount,
        query_fp: .QueryFpCount
      }
  ]
  | group_by(.bin)
  | map({
      bin:      .[0].bin,
      truth_tp: (map(.truth_tp) | add),
      truth_fn: (map(.truth_fn) | add),
      query_tp: (map(.query_tp) | add),
      query_fp: (map(.query_fp) | add)
    })
  | sort_by(
      if .bin == "50000+" then 999999999
      else (.bin | ltrimstr("[") | split(",")[0] | tonumber)
      end
    )[]
  | . as $r
  | [
      $r.bin,
      $r.truth_tp, $r.truth_fn,
      $r.query_tp, $r.query_fp,
      (if ($r.truth_tp + $r.truth_fn) > 0
       then ($r.truth_tp / ($r.truth_tp + $r.truth_fn) * 1000 | round / 1000)
       else "NaN" end),
      (if ($r.query_tp + $r.query_fp) > 0
       then ($r.query_tp / ($r.query_tp + $r.query_fp) * 1000 | round / 1000)
       else "NaN" end)
    ]
  | . as $row
  | if (($row[5] | type) == "number") and (($row[6] | type) == "number")
    then $row + [(2 * $row[5] * $row[6] / ($row[5] + $row[6]) * 1000 | round / 1000)]
    else $row + ["NaN"] end
  | @tsv
' cnvs/HG002/wittyer/Wittyer.Stats.json \
  | (echo -e "Bin\tTruthTP\tTruthFN\tQueryTP\tQueryFP\tRecall\tPrecision\tF1"; cat) \
  | column -t -s $'\t'
```

This produces a result like:

```bash
[1000, 5000)    155      233      166      12       0.399   0.933      0.559
[5000, 10000)   46       58       38       1        0.442   0.974      0.608
[10000, 20000)  16       2        18       0        0.889   1          0.941
[20000, 50000)  9        1        5        0        0.9     1          0.947
50000+          3        0        3        0        1       1          1
```

This performance was generated using XOOS 1.2.

{% endstep %}
{% endstepper %}
