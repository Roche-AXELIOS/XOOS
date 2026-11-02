# File formats

This page describes the tool-specific output format of the Pan-Genome Consensus Caller.
Usage details are on the [pangenome_consensus_caller](usage-pangenome-consensus-caller.md) page.

## Output BAM

Applies to: `pangenome_consensus_caller`.

The output is a BAM file with the same records as the input, with duplex consensus bases and base qualities corrected where the pangenome graph alignment is informative.
For each read, the following fields may change relative to the input:

| Field | Change |
|:-------|:--------|
| SEQ | A discordant consensus base is switched to the base supported by the graph alignment. |
| QUAL | Base qualities at corrected discordant positions are raised to the adjusted base quality (`-q`, default 22). |
| CIGAR | Bases clipped because the graph supports the shorter duplex read are removed, adjusting the CIGAR. |
| POS | The alignment start position is adjusted when leading operations are removed. |
| YC | The duplex-structure tag is re-encoded to stay consistent with the corrected sequence. |

The `GR` tag written by Giraffe is read as input.
Reads that lack a `YC` tag, lack a `GR` tag, or carry a `GR` value of `*` are written unchanged.

When multi-threading is enabled (`-t` greater than 1), records are written in a non-deterministic order and the output BAM is unsorted.
Sort the output before indexing:

```bash
samtools sort -@ ${thread_count} --write-index \
  -o ${output_sorted_bam}##idx##${output_sorted_bam}.bai ${output_bam}
```

Back to [pangenome_consensus_caller output](usage-pangenome-consensus-caller.md#output).
