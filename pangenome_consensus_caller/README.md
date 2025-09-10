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
pangenome_consensus_caller -b ${no_repeat_bed} -t ${thread_count} < ${input_bam} > ${output_bam}
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
