# process-mitochondrial-bam appendix

## Concepts and terminology

**Mitochondrial chromosome (`chrM`)**
: The mitochondrial reference contig and pangenome graph path used by the current workflow.

**Nuclear mitochondrial DNA segment (NuMT)**
: Mitochondrial-derived sequence located in the nuclear genome.
Such regions can attract reads that originated from mitochondrial DNA.

**NR tag**
: A BAM auxiliary tag that records off-reference pangenome graph-placement information after projecting from graph alignments to BAM.

**Graph Alignment Format (GAF)**
: A text format that represents read alignments to pangenome graph paths.

**Abundance expectation-maximization model**
: The 2-class model that estimates mitochondrial and non-mitochondrial abundance and produces per-read posterior probabilities.

## Assignment model

`process-mitochondrial-bam call` reads each placement score from the GAF `AS:i:` tag.
It initializes mitochondrial and non-mitochondrial abundance from coverage and adjusts class abundance for the effective sequence lengths.

Reference NuMT intervals always contribute to the non-mitochondrial length.
Polymorphic NuMT insertion lengths contribute when the original BAM contains supporting reads near the corresponding anchor.

The command retains a read when its posterior mitochondrial probability exceeds the threshold or its strongest mitochondrial placement is at least as strong as its strongest nuclear placement.
It converts the posterior probability to a capped Phred mapping quality before writing the retained GAF lines.

## Packaged resources

The module includes GRCh38 reference NuMT intervals, polymorphic NuMT anchors, autosomal control regions, and a precomputed mitochondrial paths GAF file.
See the [resource manifest](../resources/MANIFEST.md) for checksums and provenance.

## Performance and cost notes

Recruitment scans the original BAM file and performs indexed interval queries.
Memory use during assignment increases with the number of reads and candidate GAF placements.
The end-to-end runtime is typically dominated by loading the graph and running `vg giraffe`.
