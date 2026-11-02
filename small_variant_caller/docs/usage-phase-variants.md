# Phase variants

`phase_variants` uses read-level evidence from an aligned BAM to determine
whether nearby variants occur on the same or opposite haplotypes.
It annotates phased records in a VCF and can merge phased variants into MNV or
complex records.

**NOTE**: `phase_variants` currently only supports phasing somatic WGS samples with either the `tumor-normal-wgs` or `tumor-only-wgs` workflows.

## How it works

The command performs the following steps:

- Bi-allelic variant records from the input VCF are read in with nearby variants using `--max-phase-distance`.
- Reads that span both variants are counted and classified as cis or
  trans.
- Phase blocks are then created from pairwise results that meet the read-count and
  cis-fraction thresholds.
- Phased records are annotated with `PGT`, `PID`, and `PS` FORMAT fields, with merged MNV or complex records optionally emitted. If MNV or complex records are emitted their component records are either marked or removed according to the merge options.

The command processes only bi-allelic records.
It skips filtered records unless `--include-filtered` or `--rescue-filtered` is
enabled.

## Example command

```bash
phase_variants \
    tumor-normal-wgs \
    --bam-input MUTECT2_BAMOUT.bam \
    --vcf-input SVC_OUTPUT.vcf.gz \
    --output-dir phased \
    --vcf-output PHASED_OUTPUT.vcf.gz \
    --tumor-sample-name TUMOR_READ_GROUP
```

For `tumor-normal-wgs`, `--tumor-sample-name` is the tumor read-group
identifier in the BAM; it is not the tumor sample name from the VCF. Use this
option with a multi-sample VCF to select the tumor reads.

For `tumor-only-wgs`, use a single-sample VCF and all tagged BAM reads are
used. This subcommand does not accept `--tumor-sample-name`.

## Input

| Input | Required | Notes                                                                                                                                                                 |
|:---|:---|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Indexed BAM | Yes | Passed with `--bam-input`. A mutect2 output BAM file. The BAM must contain read-group tags. When `--tumor-sample-name` is set, it must match the tumor read-group ID. |
| Candidate VCF | Yes | Passed with `--vcf-input`. A SVC filtered VCF file. The command reads bi-allelic candidate records and preserves the original records in the output.                  |
| Reference FASTA | No | Passed with `--genome`. The indexed reference provides the backbone when merging phased variants across gaps. The input BAM provides a fallback. |

## Output

| Output | Description |
|:---|:---|
| Phased VCF | Written to `<output-dir>/<vcf-output>`. The default is `./output.phased.vcf.gz`. The output contains the input records with phasing annotations and any merged records. |
| VCF index | A Tabix index is written alongside compressed VCF output as `<vcf-output>.tbi`. |

The phasing fields and merge-related tags are described in [Variant phaser
VCF fields](file-formats.md#variant-phaser-vcf-fields).

## Tips for analysis

- A Mutect2 BAM output file and the corresponding SVC filtered somatic VCF are required so that the read alignments and candidate records describe the same analysis.
- The subcommands `tumor-normal-wgs` and `tumor-only-wgs` are intended to be used alongside their `filter_variants` counterparts. For phasing with the `tumor-normal-wgs` subcommand specify `--tumor-sample-name` as the tumor read-group ID in the BAM. For `tumor-only-wgs` ensure all tagged BAM reads belong to the same tumor sample.
- The default maximum phase distance is 10 bp, and the default minimum base
  quality is 23, matching the SVC `tumor-normal-wgs` setting.
- MNV and complex merging, and retention of component records, are enabled by
  default.
- Use `--include-filtered` to include filtered variants in phasing.
  Use `--rescue-filtered` when a filtered variant can be retained after it
  phases with a PASS variant.
- Use `--skip-ins-snv` to skip insertion and SNV pairs when those pairs produce
  unwanted phase relationships.

## CLI options

Required parameters are shown in **bold**.

This command also accepts the common `-h, --help`, `-v, --version`, and
`-l, --log-level` options.

### Workflow subcommands

| Subcommand | Description |
|:---|:---|
| **tumor-normal-wgs** | Phase a tumor-normal VCF. Requires `--tumor-sample-name` to identify tumor reads in the BAM. |
| **tumor-only-wgs** | Phase a single-sample tumor-only VCF. The command uses all tagged BAM reads and does not accept `--tumor-sample-name`. |

### Input options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--bam-input** | Path to the indexed Mutect2 BAM output. | Path to an indexed BAM file |
| **--vcf-input** | Path to the candidate VCF. | Path to a non-empty VCF file |
| --genome | Optional indexed reference FASTA used as the merge backbone. The input BAM provides a fallback. | Path to an indexed FASTA file |

### Output options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --output-dir | Directory for output files. | Directory path [default: `.`] |
| --vcf-output | Phased VCF filename, relative to `--output-dir`. | Filename [default: `output.phased.vcf.gz`] |

### Sample metadata options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --tumor-sample-name | Tumor read-group identifier in the BAM. Available and required for the `tumor-normal-wgs` subcommand. | Read-group ID |

### Phasing options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --max-phase-distance | Maximum distance between variants considered for phasing. | Integer >= 0 [default: `10`] |
| --min-shared-reads | Minimum number of reads spanning both variant positions. | Integer >= 0 [default: `2`] |
| --min-cis-reads | Minimum number of reads supporting the cis configuration. | Integer >= 0 [default: `2`] |
| --min-cis-fraction | Minimum fraction of informative reads supporting cis. | Float from `0` to `1` [default: `0.8`] |
| --min-base-quality | Minimum base quality at the variant position for a read to count. | Integer from `0` to `93` [default: `23`] |

### Annotation and merging options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --report-mnv | Add the `MNVTAG` INFO field to phased records. Merging also enables this annotation on the affected records. | Flag [default: `false`] |
| --merge-mnv | Merge phased SNV-only blocks into MNV records. | Flag [default: `true`] |
| --merge-complex | Merge phased blocks that include indels into complex records. | Flag [default: `true`] |
| --keep-components | Keep the original component records alongside merged records. | Flag [default: `true`] |

### Filtering options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --include-filtered | Include filtered variants in the phasing analysis. | Flag [default: `false`] |
| --rescue-filtered | Phase PASS and non-PASS pairs and rescue non-PASS variants when they phase with a PASS variant. | Flag [default: `false`] |
| --skip-ins-snv | Skip insertion and SNV pairs. | Flag [default: `false`] |

### Performance options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --threads | Number of worker threads. | Integer >= 0, where `0` uses available hardware threads [default: `1`] |
| --warn-as-error | Treat warning messages as errors. | Flag [default: `false`] |
