# SNP panel-of-normals blocklist workflow

Build a panel-of-normals (PoN) SNP blocklist for tumor-only WGS BAF filtering.
The blocklist marks germline SNP loci whose normal allele balance is unreliable
across a panel of normals (for example mapping artefacts or systematic
allele-balance skew), so the Copy Number Caller can drop them from its
B-allele-frequency signal.

Each site is scored by a **clean-genotype log10 Bayes factor** (`CLEAN_LOG10BF`)
rather than a hard genotype call. This recovers allele-balance-skewed regions
(for example `chrX:118,283,424-118,351,264`) that an earlier 4-state
hard-classifier blocklist missed.

Ported from an internal analysis prototype
(`libs/snakemake/10_improved-snp-pon.smk` and `bin/`).

## How the score works

For each site with normal allele counts `(ref_ad, alt_ad)`, the workflow
computes:

```text
CLEAN_LOG10BF = log10( P(counts | clean-genotype mixture) / P(counts | Beta(1,1)) )
```

- **Clean-genotype mixture** — a weighted mixture of three Beta-Binomial
  genotype models: `hom_ref ~ Beta(1, 99)`, `het ~ Beta(50, 50)`,
  `hom_alt ~ Beta(99, 1)`. Each shape has concentration `alpha + beta = 100`, so
  the three genotypes are equally confident and `hom_ref`/`hom_alt` sit at a
  realistic ~1% error rate rather than exactly 0 or 1.
- **Alternative** — a diffuse `Beta(1, 1)` (uniform allele-fraction) model.

A strongly **negative** `CLEAN_LOG10BF` means the diffuse alternative explains
the counts better, i.e. the site does not look like a clean germline genotype
and is a blocklist candidate.

## Pipeline

```text
raw somatic (mutect2) T-N VCFs           [prefilter=True only]
  └─ sort_vcf                            bcftools sort + index
       └─ filter_normal_snps            biallelic SNPs, normal DP > 15,
       │                                gnomAD-overlap, keep normal only
       ▼
per-normal filtered SNP VCFs (*.normal_dp15.vcf.gz)
  └─ classify_genotypes_bayes_factor    annotate CLEAN_LOG10BF (per sample)
       └─ slim_normal_genotypes         keep only CLEAN_LOG10BF, write BCF
            └─ merge_normal_genotypes   bcftools merge -> multi-sample VCF
                 └─ extract_blocklist   threshold + aggregate -> snp_blocklist.bed
```

The two `prefilter` steps run only when `prefilter=True`; otherwise the
workflow starts from already-filtered per-normal VCFs.

## Inputs

The workflow accepts inputs at **either of two boundaries**, chosen with the
`prefilter` flag:

- **`prefilter=False` (default)** — inputs are **per-normal filtered SNP VCFs**
  (`*.normal_dp15.vcf.gz`), one bgzipped VCF per normal sample. They are scored
  directly. Use this when the filtered VCFs already exist.
- **`prefilter=True`** — inputs are **raw somatic (mutect2) tumour-normal
  VCFs**. The workflow first runs the pre-filter (`sort_vcf` +
  `filter_normal_snps`) to derive the per-normal `normal_dp15` VCFs, then scores
  them. Use this to reproduce `snp_blocklist.bed` from raw VCFs, or to build a
  fresh blocklist from a new set of raw T-N VCFs. Requires `population_vcf`.

### Supplying inputs

Independently of `prefilter`, the sample set is supplied one of two ways
(exactly one is required):

1. **Directory + glob (ad hoc)** — set `input_dir` to a directory of VCFs;
   `input_glob` selects them and the sample name is derived from each filename
   (the VCF suffix and the `.normal_dp15` tag are stripped, so
   `DLS_7.normal_dp15.vcf.gz` becomes sample `DLS_7`).
2. **Sample sheet (reproducible)** — set `sample_sheet` to a TSV with a
   `sample<TAB>vcf` header and an optional third `normal_sample` column. This
   names each sample and its VCF path explicitly, pinning a fixed, committed
   panel independent of on-disk filenames. In `prefilter` mode, `normal_sample`
   names the normal column in the raw T-N VCF; when omitted it defaults to
   `sample + normal_sample_suffix` (default suffix `_N`).

### Reproducing the reference blocklist

A committed config and sample sheet rebuild the exact reference panel from raw
VCFs:

- `config/blocklist.yaml` — `prefilter: true` plus the reference thresholds and
  classifier weights.
- `config/samples.tsv` — the 50 DLS normals (`DLS_1`..`DLS_50`) with
  **relative** `../inputs/vcfs/<file>` paths to their raw mutect2 T-N VCFs and
  normal-sample columns. Sheet paths are resolved relative to the sheet's own
  directory (`config/`), so `../inputs/vcfs/` lands at the workflow-root
  `inputs/vcfs/` that `stage_inputs.sh` populates.

All inputs are staged from durable storage — no site-specific filesystem path is
required. `scripts/stage_inputs.sh` downloads the 50 raw T-N VCFs (~144 GiB)
from the S3 data catalog and the public af-only gnomAD VCF (~3 GB) from GCS:

```bash
scripts/stage_inputs.sh          # syncs VCFs from S3, fetches gnomAD
snakemake -s snp_blocklist.smk \
    --configfile config/blocklist.yaml \
    --use-singularity --cores 8
```

The raw T-N VCFs are stored at:

```text
s3://prod-data-catalog-20221116223806161200000001/cnc_resources/snp_blocklist_bed/inputs/vcfs/
```

Reading them requires AWS credentials with `GetObject` on the data-catalog
bucket (the bucket enforces SSE-KMS). To reproduce the blocklist in a different
environment, either run `stage_inputs.sh` (which syncs into `inputs/vcfs/`) or
provide your own `sample_sheet` pointing at locally staged VCFs.

### Pre-filter provenance

In `prefilter` mode the two upstream rules reproduce the pre-filter
from `02_bafs.smk`:

```text
{sample}_mutect2.vcf.gz               (external: Mutect2 tumor/normal VCF)
  └─ sort_vcf (bcftools sort)
       └─ filter_normal_snps (bcftools view)
            • biallelic SNPs only
            • normal-sample depth > 15  (FMT/DP)
            • restricted to gnomAD common sites (--targets-file gnomAD AF VCF)
            • keep ONLY the normal sample
            └─ {case}.normal_dp15.vcf.gz    <-- scored by this workflow
```

### Input contract

Each input VCF must satisfy:

- One bgzipped, tabix-indexed VCF per normal sample; a **single sample** per VCF
  (the normal).
- **Biallelic SNP sites only.**
- An `AD` FORMAT field per site holding `[ref_ad, alt_ad]`. This is what
  `classify_genotypes_bayes_factor.py` reads to compute `CLEAN_LOG10BF`; sites
  with a missing or incomplete `AD` are scored as missing and do not count as
  informative in `extract_blocklist`.
- Recommended (matching the prototype, not enforced here): normal depth >= 15
  and restriction to gnomAD common SNPs (AF >= 0.01).

In `prefilter=True` mode this contract applies to the workflow's
`filter_normal_snps` output rather than to the raw inputs; the raw inputs are
Mutect2 tumour-normal VCFs and only need the normal-sample `AD`/`DP` fields.

Sample names come from the sample sheet when one is given, otherwise they are
derived from the input filenames (the VCF suffix and the `.normal_dp15` tag are
stripped, so `DLS_7.normal_dp15.vcf.gz` becomes sample `DLS_7`).

## Output

`{output_dir}/extract_blocklist/snp_blocklist.bed` — a 3-column BED
(`contig`, `start` 0-based, `end`) of blocklisted SNP sites.

A site is written when at least `min_informative_samples` normals have a
non-missing `CLEAN_LOG10BF` **and** at least `min_skewed_fraction` of those
informative normals are skewed (`CLEAN_LOG10BF < bf_threshold`).

## Configuration

Override with `--config key=value`:

| Key | Default | Description |
|---|---|---|
| `input_dir` | `None` | Directory containing input VCFs. Provide exactly one of `input_dir` or `sample_sheet`. |
| `input_glob` | `*.normal_dp15.vcf.gz` | Glob to match VCFs within `input_dir` |
| `sample_sheet` | `None` | TSV (`sample<TAB>vcf[<TAB>normal_sample]` header) of samples and VCF paths. Provide exactly one of `input_dir` or `sample_sheet`. |
| `output_dir` | `outputs` | Output directory |
| `prefilter` | `False` | Treat inputs as raw T-N VCFs and run the `normal_dp15` pre-filter first |
| `population_vcf` | *(required if `prefilter`)* | Population AF VCF (e.g. af-only gnomAD) for the pre-filter |
| `min_normal_depth` | `15` | Minimum normal-sample depth in the pre-filter |
| `normal_sample_suffix` | `_N` | Normal column name = sample + this suffix (prefilter, no explicit `normal_sample`) |
| `bf_threshold` | `-0.5` | `CLEAN_LOG10BF` below this is "skewed" |
| `min_skewed_fraction` | `0.25` | Fraction of informative normals that must be skewed |
| `min_informative_samples` | `10` | Minimum informative normals per site |
| `weight_hom_ref` | `0.72` | `hom_ref` mixture weight |
| `weight_het` | `0.19` | `het` mixture weight |
| `weight_hom_alt` | `0.09` | `hom_alt` mixture weight |
| `classifier` | *(alongside workflow)* | Path to `classify_genotypes_bayes_factor.py` |
| `python` | `python3` | Python interpreter used to run the classifier (must have pysam) |

The default genotype weights are the mean Hardy-Weinberg priors over common
gnomAD SNPs (AF >= 0.01); see [Reference tools](#reference-tools).

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/)
- Python 3.11+ with [pysam](https://pysam.readthedocs.io/) (for the classifier)
- `bcftools` 1.20 and `tabix` (the bcftools rules pin the
  `quay.io/biocontainers/bcftools:1.20` container; run with
  `--use-singularity` / `--use-conda`, or provide them on `PATH`).

## Usage

Score already-filtered per-normal VCFs from a directory (ad hoc):

```bash
snakemake -s snp_blocklist.smk \
    --config input_dir=/path/to/normal_vcfs output_dir=results \
    --cores 8
```

Build a fresh blocklist from raw T-N VCFs listed in a sample sheet:

```bash
snakemake -s snp_blocklist.smk \
    --config sample_sheet=my_panel.tsv prefilter=true \
             population_vcf=/path/to/af-only-gnomad.hg38.vcf.gz \
             output_dir=results \
    --cores 8
```

Reproduce the committed reference panel (see
[Reproducing the reference blocklist](#reproducing-the-reference-blocklist)):

```bash
snakemake -s snp_blocklist.smk --configfile config/blocklist.yaml --cores 8
```

Dry run to inspect the plan without executing:

```bash
snakemake -s snp_blocklist.smk --config input_dir=/path/to/normal_vcfs -n
```

## Reference tools

### `calculate_genotype_priors.py`

A reference/analysis tool that is **not on the runtime blocklist path** (it is
not part of `rule all`). It computes per-site Hardy-Weinberg genotype priors
from population allele frequencies (gnomAD `INFO/AF`):

```text
hom_ref = (1 - p)^2   het = 2p(1 - p)   hom_alt = p^2
```

Averaging the `prior_*` columns over common gnomAD SNPs (AF >= 0.01) is how the
classifier's default mixture weights (`hom_ref=0.72`, `het=0.19`,
`hom_alt=0.09`) were chosen. It is kept here for provenance and to allow
re-deriving the weights against a different population VCF or AF cutoff.

```bash
python calculate_genotype_priors.py \
    --vcf af-only-gnomad.hg38.vcf.gz \
    --output genotype_priors.tsv \
    --min-af 0.01 \
    --snps-only
```

## Tests

A behavior-preserving self-check pins `CLEAN_LOG10BF` values against the
prototype:

```bash
python -m pytest tests/test_clean_genotype_log10bf.py
# or, without pytest:
python tests/test_clean_genotype_log10bf.py
```

## Consuming the blocklist

The `snp_blocklist.bed` produced here is a **SNP-site** blocklist consumed by
the copy-number caller's tumor-only BAF filtering via the
`--snp-blocklist-bed` option (loaded by `SnpBlocklist` in
`src/tumor-only/snp-blocklist.cpp`). It is distinct from the region blocklist
passed to `--blocklist-bed` (shipped as `resources/hg38_blocklist.bed`).

This SNP-site blocklist can be generated with this workflow, and a prebuilt hg38
SNP blocklist is also available as a DVC-tracked resource at
`copy_number_caller/resources/hg38_snp_blocklist.bed` (materialize it with `dvc pull`).
To use a blocklist, pass the resulting BED to `--snp-blocklist-bed`. This README,
together with `config/blocklist.yaml` and `config/samples.tsv`, records how the
reference resource is produced.
