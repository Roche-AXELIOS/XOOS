# process-mitochondrial-bam troubleshooting

## Original BAM index is missing

**Symptom:** The command reports that the original BAM file is not indexed or cannot open its index.

**Cause:** The input BAM file has no adjacent BAI file, or its index is stale.

**Resolution:** Coordinate-sort the BAM file and create a new index with `samtools index`.

## Mitochondrial path is not found

**Symptom:** `process-mitochondrial-bam call` cannot find the configured mitochondrial path or reports an empty mitochondrial node set.

**Cause:** The paths GAF file does not contain the value passed to `--chrm-name`.

**Resolution:** Use a paths GAF file generated from the same graph and mitochondrial path name as the alignment.

## No reads remain in the output

**Symptom:** The output GAF or BAM file contains no reads.

**Cause:** No candidate read has a mitochondrial placement that passes the assignment criteria.

**Resolution:** Confirm that recruitment produced FASTQ records, the graph contains the mitochondrial path, the mitochondrial path matches the GBZ used for realignment, and `vg giraffe` emitted multiple candidate placements.

## Graph resources do not match

**Symptom:** Graph alignment fails, or mitochondrial placements are not recognized.

**Cause:** The GBZ, minimizer, distance, zipcode, and paths GAF files originate from different graph builds.

**Resolution:** Use companion indexes and regenerate the paths GAF file from the selected GBZ.
