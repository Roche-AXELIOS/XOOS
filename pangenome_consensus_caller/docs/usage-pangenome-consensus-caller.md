# Consensus correction

## How it works

At each discordant position that falls outside the excluded regions, the module compares the R1 and R2 bases against the base the pangenome graph aligns to that position:

- If the current consensus base matches the graph, its base quality is raised to the adjusted base quality (`-q`), reflecting the extra support.
- If the graph supports the discordant (R2) base instead, the consensus base is switched to it and its base quality is raised to the adjusted base quality.
- If the graph supports the shorter of the two duplex reads at an indel, the base is clipped from the output.

Positions inside the regions listed in the excluded-regions BED (`-b`), and positions that fall in soft-clips of the graph alignment, are treated as uninformative and are left unchanged.
Reads without both a `YC` tag and a valid `GR` tag are passed through unchanged.

To use the module, input BAM files must be aligned with [Giraffe](https://github.com/vgteam/vg) using the `--add-graph-aln` parameter, which records the graph alignment in each read's `GR` tag.

{% hint style="warning" %}
**Warning:** This module has undergone limited testing.
It is intended for use only with the 7 BAM files provided with the Early Access webinar release.
{% endhint %}

## Example command

The excluded-regions BED restricts correction to non-repetitive regions.
The commands below build a BED of the complement of the UCSC Browser simple-repeat track, then run the module and sort the output:

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

{% hint style="info" %}
**Note:** If multi-threading is enabled, the output BAM is unsorted and must be sorted as shown above.
{% endhint %}

## Input

The module reads a BAM file from a positional argument or from standard input.

| Input | Required | Notes |
|:---|:---|:---|
| BAM | Yes | Duplex reads aligned with Giraffe using `--add-graph-aln`, carrying `YC` and `GR` tags. Supplied as a positional argument or on standard input. |
| BED | No | Excluded regions (`-b`); positions inside these regions are not corrected. |

## Output

The module writes a BAM file to standard output.

| Output | Description |
|:---|:---|
| BAM | Corrected reads written to standard output. Discordant consensus bases and their base qualities are corrected where the graph alignment is informative. For the fields that may change, see [Output BAM](file-formats.md#output-bam). |

When multi-threading is enabled, the output is unsorted; sort it as shown in [Example command](#example-command).

## Tips for analysis

- Use the excluded-regions BED (`-b`) to skip repetitive regions, where graph alignments are less reliable.
- Run single-threaded to keep the output in input order when a sorted output is required and re-sorting is inconvenient.

## CLI options

The Pan-Genome Consensus Caller reads a BAM from a positional argument or standard input and writes a BAM to standard output.
This command also accepts the common `-h, --help` option.

### Input options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| \<bam\> | Input BAM of Giraffe-aligned duplex reads, given as a positional argument. Omit to read from standard input. | File path [default: standard input] |
| -b, --bed | BED file of regions to exclude from pangenome-based correction. | File path [default: none] |

### Consensus options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| -q, --adjusted-qual | Base quality assigned to a base whose correction is supported by the graph alignment. | Integer [default: `22`] |
| -i, --bdry-ins-excl | Treat insertions on an excluded-region boundary as excluded. | No value (omit to set false; include to set true) |

### Performance options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| -t, --threads | Number of worker threads. Values greater than 1 produce an unsorted output. | Integer [default: `1`] |
