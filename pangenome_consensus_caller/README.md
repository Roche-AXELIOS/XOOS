# pangenome_consensus_caller README

## Overview

During duplex consensus, if R1 and R2 bases disagree, the R1 base is typically chosen.
However, pangenome alignment of the consensus may reveal that selecting the R2 base results in a better alignment.
`pangenome_consensus_caller` is designed to identify these cases and update the consensus to use the R2 base.
To ensure compatibility, input BAM files must be aligned with [Giraffe](https://github.com/vgteam/vg) using the `--add-graph-aln` parameter.

{% hint style="warning" %}
**WARNING:** This module has undergone limited testing. It is intended for use only with the 7 BAM files provided with the Early Access webinar release.
{% endhint %}

## Example Usage

```bash
curl -OL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai

curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz \
  | zcat \
  | grep Simple_repeat \
  | cut -f 6,7,8 \
  | awk 'NR==FNR{valid[$1]=1; next} $1 in valid' <(cut -f1 GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai) - \
  | sort -k1,1 -k2,2n -k3,3n \
  | bedtools merge \
  | bedtools complement -g <(cut -f1,2 GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai | sort -k1,1) -i - \
    > no_repeat.bed

pangenome_consensus_caller -b no_repeat.bed -t ${thread_count} < ${input_bam} > ${output_bam}
samtools sort -@ ${thread_count} --write-index -o ${output_sorted_bam}##idx##${output_sorted_bam}.bai ${output_bam}
```

**Note:** If multi-threading is enabled, the output BAM will be unsorted and must be sorted as shown above.

| Parameter | Description                                                                      | Value(s) |
|-----------|----------------------------------------------------------------------------------|----------|
| -b        | A bed file containing the complement of the UCSC Browser's simple repeat track.  |          |
| -t        | The number of threads.                                                           |          |

## Troubleshooting

| Issue                       | Description and Mitigation                                                                             |
|-----------------------------|--------------------------------------------------------------------------------------------------------|
| The output BAM is unsorted. | Sort the output BAM, or run without multi-threading. This issue will be addressed in a future version. |
