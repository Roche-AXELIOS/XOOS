# Usage

The STR Caller exposes a single `str_caller` command. This page links to its
dedicated Usage page, where you will find how it works, an example command,
input and output details, tips, and the full CLI options.

## STR caller

Genotypes and detects repeat expansions in short tandem repeats (STRs) from
Roche SBX sequencing data. Input: **BAM** + **reference FASTA** → Output:
**VCF**, **JSON**, **metrics**, and an optional realigned **BAM**.

```shell
str_caller --bam sample.bam --reference hg38.fa --output-dir results
```

The caller uses the catalog packaged in the Docker image.

For how it works, inputs, outputs, tips, and the full CLI options, go to the
[STR caller](usage-str-caller.md) page.
