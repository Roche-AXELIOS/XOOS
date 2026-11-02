# Usage

The full process to recover mitochondrial reads that have been misaligned to NuMT sequences is implemented as part of the [xoosnf pipeline](../../docs/xoosnf/README.md), with the specific workflow within the pipeline documented [here](../../docs/xoosnf/mitochondrial-bam-post-processing.md). See this documentation section for an overview of how the mitochondrial BAM post-processing steps fit into the variant-calling workflow.

## Running the end-to-end workflow

The docker image for the module contains a helper script to run the steps using  `process-mitochondrial-bam` and `vg`, that will take in an input BAM and produce a realigned chrM BAM with updated chrM read alignments.

Input: **BAM file and pangenome indexes** → Output: **sorted and indexed mitochondrial BAM file**.

```shell
run-process-mitochondrial-bam.sh \
  --orig-bam sample.bam \
  --gbz graph.gbz \
  --min graph.min \
  --dist graph.dist \
  --zipcodes graph.zipcodes \
  --out-bam sample.chrM.bam
```

To run the helper script, go to [Running the end-to-end workflow](usage-run.md).

## Recruiting reads

Select mitochondrial and NuMT-associated reads for graph realignment.

Input: **NR-tagged BAM file and NuMT resources** → Output: **read names, merged regions, and FASTQ file**.

```shell
process-mitochondrial-bam recruit \
  --orig-bam sample.bam \
  --nr-positions resources/numt_catalog.grch38.nr.tsv \
  --ucsc-bed resources/nuMtSeq_hg38.ucsc.bed \
  --out-fastq recruited.fastq
```

To run the recruitment stage, go to [Recruiting reads](usage-recruit.md).

## Assigning mitochondrial reads

Estimate mitochondrial origin probabilities and retain mitochondrial graph placements.

Input: **multimapping GAF and original BAM files** → Output: **mitochondrial GAF file with updated mapping qualities**.

```shell
process-mitochondrial-bam call \
  --gaf all.gaf \
  --paths-gaf resources/chrm_paths.grch38.gaf \
  --orig-bam sample.bam \
  --control-regions resources/control_regions.grch38.bed \
  --out-gaf chrM.gaf
```

To run the assignment stage, go to [Assigning mitochondrial reads](usage-call.md).
