# Troubleshooting

## Output BAM is unsorted

**Symptom:** The output BAM is not coordinate-sorted, and downstream tools that require a sorted, indexed BAM fail.

**Cause:** When multi-threading is enabled (`-t` greater than 1), records are written in the order worker threads complete, which is not coordinate order.

**Resolution:** Sort the output with `samtools sort`, or run single-threaded to preserve input order:

```bash
samtools sort -@ ${thread_count} --write-index \
  -o ${output_sorted_bam}##idx##${output_sorted_bam}.bai ${output_bam}
```

## Reads are not corrected

**Symptom:** The output BAM is identical to the input, or specific reads are unchanged.

**Cause:** The Pan-Genome Consensus Caller only corrects reads that carry both a `YC` duplex tag and a valid `GR` graph-alignment tag.
Reads without a `YC` tag, without a `GR` tag, or with a `GR` value of `*` are passed through unchanged.
Positions inside the excluded-regions BED (`-b`) and positions in soft-clips of the graph alignment are also left unchanged.

**Resolution:** Confirm the input BAM was aligned with [Giraffe](https://github.com/vgteam/vg) using the `--add-graph-aln` parameter so that the `GR` tag is present.
Check that the regions of interest are not covered by the excluded-regions BED.

If you have issues or concerns not documented here, please create an issue on [GitHub](https://github.com/Roche-AXELIOS/XOOS/issues) or reach out to your local Roche affiliate.
