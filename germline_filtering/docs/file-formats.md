# File formats

This page is the reference for the input and output files that Germline Variant Filtering reads and produces.
Running each subcommand, with its full CLI options, is covered on the [Train a model](usage-train.md) and [Apply a model to a sample](usage-predict.md) pages.
The per-sample QC metrics file is documented on the [Metrics](metrics.md) page.

The tool does **not** generate its input VCF.
You must produce the input with your own variant-calling workflow and ensure it contains the header lines and fields described below. The `predict` subcommand generates the tagged output VCF documented below.

## JSON configuration file

Applies to: `train`.

The `train` subcommand reads a JSON configuration file that lists the samples and their VCF/SEG inputs.
The `predict` subcommand does not use JSON; it takes `--vcf` and `--seg` directly.

```json
{
  "samples": [
    {
      "sample": "tumor_sample_1",
      "vcf": "path/to/sample1.vcf",
      "seg": "path/to/sample1.seg"
    },
    {
      "sample": "tumor_sample_2",
      "vcf": "path/to/sample2.vcf",
      "seg": "path/to/sample2.seg"
    }
  ]
}
```

| Field | Definition |
|:------|:-----------|
| samples | Array of sample objects. |
| sample | Unique sample name. |
| vcf | Path to the VCF file with variant calls. |
| seg | Path to the copy-number SEG file. |

See *tests/resources/template.json* for a working example.

## VCF file

Applies to: `train`, `predict`.

The same VCF format is used for both subcommands.
The only difference is that training requires ground-truth germline/somatic labels in the `FILTER` column, whereas prediction does not.

### Required header lines

| Header line | Definition |
|:------------|:-----------|
| ##tumor_sample=\<sample_id\> | Identifies which sample column holds the tumor `DP`/`AF` values. `<sample_id>` must exactly match one of the sample names in the `#CHROM` line. **Mandatory for multi-sample VCFs**: if it is missing or matches no sample column, the tool fails with an error. For a single-sample VCF, a missing or non-matching value logs a warning and the only sample column is used as the tumor sample. |
| ##INFO=\<ID=POPAF,...\> | Declares the population allele-frequency INFO field. |
| ##FORMAT=\<ID=DP,...\> | Declares the per-sample read-depth field. |
| ##FORMAT=\<ID=AF,...\> | Declares the per-sample allele-fraction field. |
| ##FILTER=\<ID=germline,...\> | Declares the germline label filter used for training. `##FILTER=<ID=normal_artifact,...>` is optional. |

### Required fields per record

| Field | Location | Type | Definition |
|:------|:---------|:-----|:-----------|
| POPAF | `INFO` | Optional float | Population allele frequency (for example, from gnomAD). An absent value is NaN-encoded; an explicit `0` remains a real observation. |
| DP | `FORMAT` (tumor sample) | Integer | Read depth at the variant position in the tumor sample. |
| AF | `FORMAT` (tumor sample) | Float | Variant allele fraction in the tumor sample. |

`DP` and `AF` are read from the sample column named in `##tumor_sample` (or from the only sample column of a single-sample VCF).
`POPAF` is read from the shared `INFO` column.

### FILTER conventions (labels)

The `FILTER` column carries the ground-truth label used for training and for optional predict metrics.
The tool recognizes exactly these combinations:

| FILTER value | Meaning | Label |
|:-------------|:--------|:------|
| PASS | Somatic variant | `0` (somatic) |
| germline | Germline variant | `1` (germline) |
| germline;normal_artifact | Germline variant | `1` (germline) |
| anything else | Not a usable label | Record skipped (counted as "filtered by other") |

For prediction on unlabeled tumor-only (TO) data, records should carry `FILTER=PASS` so they are scored.
Records with any other filter are not scored, but they can still pass through to the tagged VCF with `--emit-all`.

### Filtering behavior applied by the tool

- Both SNVs and indels (`len(REF) != len(ALT)`) are scored.
  The `is_indel` column is physical-matrix metadata that records the variant type and routes each record to its
  matching per-type model; it is not a model input.
  The number of indels flowing to scoring is reported in the QC metrics (`features_indel_scored`).
- Phased compound records (`PASS`, `len(REF) > 1` and `len(ALT) > 1`) are scored too, but `is_indel` is decided
  from the first variant in `INFO/MNVTAG` rather than the compound record's own REF/ALT length; all other
  `INFO/MNVTAG` entries are ignored. The variant ID and remaining features are unaffected. A phased compound
  with a missing or invalid `INFO/MNVTAG` is skipped and counted in
  `vcf_phased_variants_skipped_invalid_mnvtag`.
- Sex chromosomes (`chrX`, `chrY`) are excluded.
- Variants outside all SEG segments (for example, subtelomeric or inter-segment) use a diploid copy-number
  fallback (total copy number 2, minor allele count 1) during feature generation, so they are not dropped before
  scoring. In `predict`, these records are nevertheless emitted with
  `FILTER=nocall_segment_missing` and `INFO/GERMLINE_STATUS=nocall_segment_missing`; they are excluded from both
  the germline and somatic sets. The count is reported in the QC metrics (`features_cnv_gap_fallback`).
- Variants with tumor `DP` below `--min-tumor-depth` (default 20) are dropped.
- Duplicate variant IDs within a sample are de-duplicated (last record wins).

### Minimal example

```text
##fileformat=VCFv4.2
##tumor_sample=003T
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=germline,Description="Germline variant">
##INFO=<ID=POPAF,Number=1,Type=Float,Description="Population allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 003T 003N
chr1 75206 . G A 105.64 germline POPAF=0.233 GT:AD:AF:DP 0/1:50,8:0.5:59 0/0:0,0:0:59
chr1 930188 . C T 220.10 PASS POPAF=1e-06 GT:AD:AF:DP 0/1:40,20:0.33:60 0/0:0,0:0:60
```

Here `003T` is the tumor sample declared in `##tumor_sample`.
The first record is a labeled germline variant (training label 1); the second is a labeled somatic variant (training label 0).

**Reference genome:** GRCh38 (based on test data).
The VCF and SEG files must all use the same build.

### How to create these VCF files

The training/prediction VCF is expected to be the output of a tumor-only somatic caller with a few annotations layered on top.
A typical recipe (tools are examples; use your organization's validated pipeline):

1. **Call variants** on the tumor BAM with a somatic caller (for example, Mutect2 in tumor-only mode, or your in-house caller).
   This produces raw SNV/indel calls with per-sample `FORMAT` fields.
2. **Ensure `DP` and `AF` are present** in the tumor sample's `FORMAT` column.
   Most callers emit `DP` directly; if your caller emits allelic depths (`AD`) but not `AF`, compute `AF = alt_AD / DP` and add it as a `FORMAT/AF` field.
3. **Annotate population allele frequency** as `INFO/POPAF`.
   Intersect the calls against a population database such as gnomAD (for example, with `bcftools annotate -a gnomad.vcf.gz -c INFO/POPAF`) so records carry `POPAF`.
   `POPAF` is used as the raw fraction; a value of `0` (variant absent from the database) is a real, informative observation and is kept as such. A record whose `POPAF` annotation is entirely absent is still featurized and scored — its `popaf` feature is marked missing (`popaf_is_missing`) and the model routes on the remaining features — but annotating `POPAF` is strongly recommended, since population frequency is one of the most discriminative signals.
4. **Set the `##tumor_sample` header line** to the tumor sample name so the tool knows which sample column to read (for example, `bcftools reheader` or `sed` in the header block).
5. **Assign FILTER labels.** This is the step that distinguishes training input from prediction input:
   - *For training*, label each record's `FILTER` as `germline` (or `germline;normal_artifact`) for known germline variants and `PASS` for known somatic variants.
     Labels typically come from a matched-normal comparison, a germline database, or a curated truth set.
     Only labeled records are used; unlabeled records are skipped.
   - *For prediction*, leave callable variants as `FILTER=PASS`.
     No germline labels are required.
6. **(Optional) bgzip and tabix** the VCF.
   Plain *.vcf* and *.vcf.gz* inputs are both accepted.

The corresponding per-sample copy-number SEG file (with purity/ploidy) is produced separately by the copy-number caller; see the SEG file requirements below.

## SEG file

Applies to: `train`, `predict`.

Tab-delimited copy-number segmentation file with a header row.

**Required columns:**

| Column | Definition |
|:-------|:-----------|
| Contig | Chromosome (for example, chr1, chr2). |
| Start | Segment start position (1-based). |
| End | Segment end position (inclusive). |
| NumSnp | Number of SNPs supporting the segment's allelic (BAF) call. `0` marks the segment as having no allelic evidence. |
| TotalCopyNumber | Total copy number for the segment. |
| MinorCopyNumber | Minor allele copy number. `0` means loss of heterozygosity. |
| MajorCopyNumber | Major allele copy number. |
| Purity | Tumor purity estimate (0 to 1). |
| Ploidy | Tumor ploidy estimate. |
| MirroredBaf | Mirrored B-allele frequency for the segment. |
| MeanDh | Mean decrease of heterozygosity for the segment. |
| IsInAllosome | Whether the segment lies on a sex chromosome. Allosomes are excluded from the genome-wide LOH burden. |

**Example of required columns:**

```tsv
Contig Start End NumSnp TotalCopyNumber MinorCopyNumber MajorCopyNumber Purity Ploidy MirroredBaf MeanDh IsInAllosome
chr1 1 50000000 1200 2 1 1 0.75 2.1 0.499 0.09 false
chr1 50000001 100000000 800 3 1 2 0.75 2.1 0.663 0.21 false
```

`NumSnp`, `MirroredBaf`, `MeanDh`, and the two allele-specific copy numbers feed the [LOH features](appendix.md#loh-features).
The copy-number caller writes `nan` in the BAF-derived columns for segments with `NumSnp` of `0`; the tool detects that and marks those variants `loh_is_unknown` rather than treating the resulting `0` minor copy number as loss of heterozygosity.

The SEG file may contain more columns.
The tool assumes that the SEG file contains genome-wide segmentation data.
Ensure that the SEG file is properly formatted and corresponds to the same reference genome as the VCF.

## Model files

Applies to: `train` (writes), `predict` (reads).
Default `train` paths: `<output-prefix>_snv_model.txt`, `<output-prefix>_indel_model.txt`.

Model files are LightGBM text files: one model trained on SNVs and one trained on indels. Run `train` with a
labeled cohort to create both files. Then pass them to `predict` with `--snv-model` and `--indel-model`; `predict`
routes each variant to the model matching its type.

When model paths are omitted, `predict` uses the SNV and indel models shipped with the tool. A model file may be
plain text or gzip-compressed (*.gz*).

## Prediction CSV

Applies to: `predict`.
Default path: `<output-prefix>_prediction.csv`.

Contains all feature columns plus:

| Column | Definition |
|:-------|:-----------|
| prediction | Predicted probability of being germline (0 to 1). |
| original_label | Ground-truth label when available. |

**Classification rule:** a fixed candidate-tier cutoff (`--min-probability-score`, default `0.96`) is applied to candidate variants and a separate, lower fixed cutoff to population-known variants.
Two out-of-model overrides can additionally reclassify a model-somatic call as germline; see [Classification overrides](appendix.md#classification-overrides) for the AF override, the POPAF override, and the tumor-fraction-scaled effective cutoff.

## Tagged VCF

Applies to: `predict`.
Default path: `<output-prefix>.vcf.gz` (with a `.vcf.gz.tbi` tabix index).

The single-sample tagged VCF, written alongside the prediction CSV.
The output is **BGZF-compressed** (*.vcf.gz*) and a **tabix index** (*.vcf.gz.tbi*) is written next to it so downstream tools can random-access it.

**Annotations added to every emitted record:**

| Field | Type | Definition |
|:------|:-----|:-----------|
| INFO/GERMLINE_STATUS | String | The final per-variant status token (see the status table below). Present on **every** emitted record. |
| INFO/GERMLINE_PRED | String | The model's germline probability. Present **only** on `germline` and `somatic` (PASS) records; the key is omitted entirely for every other status (no `NA` sentinel). The field is named to avoid confusion with any upstream `INFO/PRED` (small-variant-caller real-call probability) present in the input. |

**Status vocabulary.** `INFO/GERMLINE_STATUS` carries one of these VCF-safe tokens, and (for the statuses that have one) the matching `FILTER` value is set. A no-call or ambiguous status excludes the record from **both** the germline and the somatic (PASS) sets, so `bcftools view -f PASS` drops it.

| GERMLINE_STATUS | FILTER value | Meaning |
|:----------------|:-------------|:--------|
| germline | germline | Model-germline record (`GERMLINE_PRED` at or above the applicable resolved cutoff). |
| somatic | PASS | Somatic record (input `PASS` preserved). |
| high_af | high_af | Model kept somatic but the tumor AF triggers the AF override. Standalone filter (no `germline`). |
| common_germline | common_germline | Model and AF override kept somatic but `INFO/POPAF` triggers the POPAF override. Standalone filter. |
| nocall_cyp2d6 | nocall_cyp2d6 | Variant in the CYP2D6/CYP2D7 locus (paralog/pseudogene mismapping). |
| nocall_hla | nocall_hla | Variant in the hyperpolymorphic MHC/HLA region. |
| nocall_purity | nocall_purity | Sample purity above 0.95, where the tumor/normal AF geometry is degenerate. |
| nocall_segment_missing | nocall_segment_missing | No copy-number segment covers the variant (diploid fallback used); applies to scored and unscored records when copy-number input quality is supplied. |
| ambiguous_cna | ambiguous_cna | The covered segment's copy number is untrustworthy (no BAF evidence, or CN=0); applies to scored and unscored records when copy-number input quality is supplied. |
| nocall_coverage | nocall_coverage | Unscored, not covered by a copy-number input-quality status, and tumor `DP` is below `--min-tumor-depth`. |
| unknown | *(unchanged)* | Unscored and not caught above (sex chromosome / NaN features); the input `FILTER` is left unchanged. |

The `high_af` and `common_germline` FILTER values identify variants reclassified by the corresponding
override. The AF override is evaluated before the POPAF override.

The output VCF adds header definitions for `GERMLINE_STATUS`, `GERMLINE_PRED`, and every FILTER value in the
table except `unknown`. Input FILTER header definitions are preserved.

**Provenance metadata.** The output header also records how the file was produced:

```text
##RocheCommandLine=<ID=germline_filtering,Version="<version>",CommandLine="<rendered command>">
```

The rendered command line includes options resolved from defaults (for example `--snv-model`) as well as
those passed explicitly, so it names the inputs the run actually used. Any `##RocheCommandLine` lines
already present on the input VCF are **preserved** and the germline-filtering line is appended after them,
so a file carries the full provenance chain of every tool that wrote it. The header is identical whether
`predict` runs single-threaded or in parallel (`--threads`).

`Version` is taken from the build's `VERSION` value and reads `unknown` for a build configured without it.

Use `INFO/GERMLINE_STATUS` to read the final classification. See
[Classification overrides](appendix.md#classification-overrides) for override cutoffs and
[No-call statuses](appendix.md#no-call-statuses) for no-call conditions.

## Emitted records

The tool scores `PASS` SNVs and indels. PASS records that cannot be scored, such as sex-chromosome variants or
records below `--min-tumor-depth`, remain in the output with an appropriate no-call or `unknown` status.

`--emit-all` controls whether records that were already filtered in the input are retained:

| Option | Records written |
|:-------|:----------------|
| `--emit-all` (default) | Every input record. |
| `--no-emit-all` | Only records with an input `FILTER=PASS`. |

This VCF is the recommended input for downstream somatic-only filtering.
