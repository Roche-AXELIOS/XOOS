# Appendix

## Concepts and terminology

| Term | Definition |
|:-----|:-----------|
| Germline variant | An inherited variant present in an individual's normal cells. |
| Somatic mutation | A tumor-specific mutation not present in the germline. |
| Tumor-only (TO) | Sequencing of the tumor without a matched normal sample. |
| POPAF | Population allele frequency (`INFO/POPAF`), for example from gnomAD. |
| SEG file | Copy-number segmentation file with purity and ploidy. |
| Purity | Fraction of tumor cells in the sample. |
| Ploidy | Average genome copy number of the tumor. |
| AF | Allele fraction — the fraction of reads supporting the variant. |
| GBDT | Gradient Boosted Decision Trees, the LightGBM model family used. |
| POPAF-blind | Optional (`--popaf-blind`) training mode in which the model does not learn from POPAF; population frequency is applied as an outside rule. Off by default. |
| LOH | Loss of heterozygosity — a genomic segment that has lost one parental allele, so its minor copy number is 0. |
| BAF | B-allele frequency — the allelic-ratio signal the copy-number caller infers minor copy number from. Absent on segments with no overlapping SNPs. |

## Feature engineering

The tool generates 31 features per variant:

| # | Feature | Definition |
|:--|:--------|:-----------|
| 1 | `count` | Total number of variants in the sample. |
| 2 | `af` | Observed variant allele fraction from the VCF. |
| 3 | `popaf` | Population allele frequency (`INFO/POPAF`, for example gnomAD), used as the raw fraction in `[0, 1]` with no transform (lower values = rarer variant). A raw `popaf == 0` (variant absent from the population database, the rarest / most private class) is kept as a real `0.0`; only an entirely absent POPAF annotation is NaN-encoded (see `popaf_is_missing`). Using the raw fraction keeps the model, the input VCF, and the outside POPAF override on identical units, so there is no train/predict encoding to keep in sync. |
| 4 | `dp` | Sequencing depth at the variant position. |
| 5 | `common_variant_cluster` | Number of *other* common variants — tumor AF ≥ 0.85 — within a ±2 Mbp window (4 Mbp total) centred on the variant. Binary search brackets the window on sorted positions, then the bracketed range is scanned for the AF cut. The variant of interest is matched by id and excluded, so co-located records (multi-allelic sites, an SNV and an indel at one coordinate) still count as neighbours. Restricting the count to near-homozygous calls turns a raw variant density into a local germline-haplotype density: linked germline variants share the same high-AF regime, so a somatic variant sitting inside such a block is measured against that linkage rather than against how variant-dense the region happens to be. |
| 6 | `purity` | Tumor cell fraction estimate (from SEG). |
| 7 | `ploidy` | Average genome copy number (from SEG). |
| 8 | `expected_copy_number` | Expected total copy number at the variant locus (SEG `ExpectedTotalCopyNumber`), real-valued and subclonal-aware. |
| 9 | `minor_allele_count` | Copy number of the minor allele. |
| 10 | `minimum_expected_af_dist_dw` | Minimum distance between observed AF and the expected AF for the major allele `((1-purity) + purity*major_CN) / dna`, the minor allele `((1-purity) + purity*minor_CN) / dna`, or homozygous `1.0`, where `dna = 2*(1-purity) + purity*total_CN`. See [depth-weighted expected allele fractions](#depth-weighted-expected-allele-fractions). |
| 11 | `is_high_cn` | Binary indicator for total copy number > 5. |
| 12 | `is_indel` | Binary indicator for `len(REF) != len(ALT)`. For a phased compound record (`PASS`, `len(REF) > 1` and `len(ALT) > 1`), this is instead derived from the first variant listed in `INFO/MNVTAG` (e.g. `chr1:3574057C>T` in `MNVTAG=chr1:3574057C>T,chr1:3574070C>A`), so the compound is routed by its first component's type rather than its own combined REF/ALT length. A phased compound without a usable first `MNVTAG` variant is skipped. |
| 13 | `popaf_is_missing` | Physical matrix bookkeeping indicating that the POPAF annotation was absent from the record entirely; paired with the NaN-encoded `popaf`. This column is not used by the model. Distinct from a raw `popaf == 0`, which is a real observation (variant absent from the population database) and is *not* flagged missing. |
| 14 | `is_loh` | Binary indicator that the SEG segment covering the variant has lost one parental allele (minor copy number `0`). NaN when no allelic evidence is available; see [LOH features](#loh-features). |
| 15 | `allelic_imbalance` | `(major_CN - minor_CN) / total_CN`, clamped to `[0, 1]`. The continuous generalization of `is_loh`: captures imbalance short of full LOH (for example 3+1), which shifts the germline allele-fraction peak without eliminating it. NaN under the same conditions as `is_loh`. |
| 16 | `mirrored_baf` | Segment `MirroredBaf` — the direct B-allele-frequency observation. NaN when no allelic evidence is available. |
| 17 | `mean_dh` | Segment `MeanDh` (mean decrease of heterozygosity). NaN when no allelic evidence is available. |
| 18 | `sample_loh_fraction` | Genome-wide autosomal LOH burden for the sample: LOH segment length divided by assessable autosomal segment length. Constant across the sample's variants, like `purity` and `ploidy`. |
| 19 | `loh_is_unknown` | Binary indicator that features 14 to 17 are NaN-encoded because no allelic evidence was available at the variant's locus; paired with those NaN values. This column is not used by the model. Distinct from `is_loh == 0`, which is a real observation (the segment retains both alleles). |
| 20 | `expected_af_het_major_dw` | Purity- and copy-number-adjusted allele fraction a germline heterozygous variant is expected at, `((1-purity) + purity*major_CN) / (2*(1-purity) + purity*total_CN)` (the same major-allele expectation used by `minimum_expected_af_dist_dw`). Given to the model directly so it can read the observed AF relative to the shifted germline mode instead of a fixed ~0.5. NaN when total copy number is `0`. See [purity-normalized allele fraction](#purity-normalized-allele-fraction). |
| 21 | `af_germline_ratio_dw` | Observed `af` divided by `expected_af_het_major_dw`; ≈ `1.0` for a germline variant at any purity, which removes the purity-driven shift of the raw `af` feature. NaN when the expectation is undefined (total copy number `0`) or ≈ `0` (the degenerate full-purity, zero-major-copy corner). |
| 22 | `af_zform_dw` | `(af - expected_af_het_major_dw) / expected_af_het_major_dw`; centered at `0` for a germline variant. Equal to `af_germline_ratio_dw - 1`, so a tree splits on it identically; kept for a 0-centered reading of the same signal. NaN under the same conditions as `af_germline_ratio_dw`. |
| 23 | `expected_af_somatic` | Allele fraction a clonal somatic variant on a single tumor copy is expected at, `purity / (2*(1-purity) + purity*total_CN)`. The germline expectations above say where a germline variant should land but leave the somatic mode implicit; in an aneuploid genome the two move independently, so naming the somatic reference point gives the model both sides of the decision. NaN when total copy number is `0`. |
| 24 | `af_somatic_dist` | `abs(af - expected_af_somatic)`, the somatic-side analogue of `minimum_expected_af_dist_dw`. Supplied directly because a tree cannot cheaply synthesize an absolute difference by splitting on `af` and the expectation separately. NaN under the same condition as `expected_af_somatic`. |
| 25 | `z_germline_major_dw` | Observed `af` minus the germline major-allele expectation, divided by the binomial sampling standard deviation at that expectation. The features above measure the same comparison in raw allele-fraction units, where a given distance means very different things at 30x and 300x; dividing by the sampling noise puts every locus on one scale. NaN when the expectations are undefined (expected total copy number `0`) or `dp` is `0`. |
| 26 | `z_germline_minor_dw` | The same deviation from the germline minor-allele expectation. Under LOH the minor allele is the lost one, and this is the mode a germline variant whose allele was deleted falls onto. |
| 27 | `z_somatic_single_dw` | The same deviation from the clonal single-copy somatic expectation. |
| 28 | `min_abs_z_germline_dw` | Distance to the closest germline hypothesis (major, minor, or homozygous), in sampling standard deviations. |
| 29 | `min_abs_z_somatic_dw` | Distance to the closest somatic hypothesis (minor, major, or single copy), in sampling standard deviations. |
| 30 | `z_margin_dw` | `min_abs_z_germline_dw - min_abs_z_somatic_dw`; positive when the observation sits closer to a somatic hypothesis than to any germline one. Supplied directly because a tree cannot cheaply synthesize a difference of two absolute values by splitting on the two columns separately. |
| 31 | `af_separability_dw` | Closest approach between the germline and somatic hypothesis sets, in sampling standard deviations. Depends only on purity, copy number and depth, never on the observation, so it says whether the allele fraction can decide this variant at all. See [Purity and separability](#purity-and-separability). |

### Model features

The tool calculates all 31 features above and writes them to the prediction CSV. The model trains and predicts
with 23 features; these eight generated columns are retained for output and bookkeeping but are not model inputs:

- `count`
- `minor_allele_count`
- `is_high_cn`
- `is_indel`
- `popaf_is_missing`
- `sample_loh_fraction`
- `loh_is_unknown`
- `af_zform_dw`

By default, the remaining 23 features, including `popaf`, are used for both model training and prediction.
Under the optional `--popaf-blind` training mode (see [POPAF-blind training](#popaf-blind-training)), `popaf`
remains a model input but is assigned a LightGBM `feature_contri` weight of `0`, so the booster does not split on it.

Variants outside all SEG segments use a diploid copy-number fallback (total copy number 2, minor allele count 1) so the copy-number-derived features stay valid instead of producing NaN and dropping the variant before scoring. The fallback is an assumption rather than an observation, so it does **not** extend to the LOH features: those variants are marked `loh_is_unknown`. In `predict`, they are emitted with `FILTER=nocall_segment_missing` and excluded from both the germline and somatic sets.

## LOH features

A segment is loss of heterozygosity (LOH) when its **minor copy number is 0** — the same definition the copy-number caller uses for its genome-wide LOH burden.

### Why LOH matters here

In an LOH segment a germline heterozygous variant has lost one allele, so its tumor allele fraction is pushed off the roughly 0.5 germline peak toward 0 or 1 — onto the somatic mode. The LOH features provide this allele-specific context, and `sample_loh_fraction` supplies a genome-wide measure of how much of the sample's autosomal genome is affected.

### Missing allelic evidence

The copy-number caller emits `nan` for the BAF-derived SEG columns on segments with no overlapping SNPs, and those cells parse to `0.0`. A minor copy number of `0` alone therefore cannot distinguish real LOH from absent evidence, so the tool tracks whether a segment carries allelic evidence at all.

When it does not — or when the variant falls in a segmentation gap — `is_loh`, `allelic_imbalance`, `mirrored_baf`, and `mean_dh` are NaN-encoded and `loh_is_unknown` is set. LightGBM routes the NaN natively. Both `loh_is_unknown` and `popaf_is_missing` remain as physical matrix bookkeeping but are masked from model training and prediction. Those variants are still scored on their remaining features rather than dropped.

Segments with no allelic evidence are excluded from **both** sides of the `sample_loh_fraction` ratio, so thin SNP support lowers confidence rather than fabricating LOH. Allosomes are excluded as well, because their log ratios come from an autosome-anchored normalization. The `seg_segments_without_baf` QC metric reports how much of the segmentation was unassessable.

## Purity-normalized allele fraction

The raw `af` feature shifts with tumor purity and local copy number: the allele fraction of a germline heterozygous variant is not fixed at ~0.5 but at the depth-weighted expectation below, which moves as purity and the copy-number state change. A booster trained on `af` therefore has to learn allele-fraction thresholds separately for each purity regime, and generalizes poorly to samples whose purity sits outside the training cohort — the observed cross-sample failure mode on high-purity, aneuploid tumors.

`expected_af_het_major_dw`, `af_germline_ratio_dw`, and `af_zform_dw` recast the observed AF against that germline major-allele expectation so a germline variant lands near a purity-independent fixed point (`af_germline_ratio_dw ≈ 1`, `af_zform_dw ≈ 0`) regardless of purity, while a low-AF somatic variant falls below it. This gives the model a purity-invariant view of the same allele-fraction signal alongside the raw `af` and `minimum_expected_af_dist_dw` features, which are retained.

The ratio and z-form are NaN-encoded (and exempted from the NaN row-drop, like the LOH block and `popaf`) in the two degenerate corners where the expectation is undefined or ~0: total copy number `0`, or full purity with zero major copies. `expected_af_het_major_dw` itself is NaN only when total copy number is `0`.

## Depth-weighted expected allele fractions

Reads are sampled in proportion to DNA, not to cell count. At a locus, the `1 - purity` normal cells contribute two copies each and the `purity` tumor cells contribute `total_CN` each, so every expected allele fraction shares the denominator:

```text
dna = 2*(1-purity) + purity*total_CN
```

A germline variant is present in the normal cells too, which are heterozygous and so contribute `1 - purity` alt copies; a somatic variant is absent from them and they contribute to the denominator only:

```text
expected germline AF (allele_count tumor copies) = ((1-purity) + purity*allele_count) / dna
expected somatic  AF (allele_count tumor copies) = (purity*allele_count) / dna
```

Substituting `total_CN = 2` reduces the germline form to `(1-purity)*0.5 + (allele_count/total_CN)*purity`, so a cell-weighted expectation is correct at diploid loci and drifts everywhere else — increasingly so as the locus departs from `total_CN = 2`. Deletions and amplifications are exactly where the germline and somatic modes are hardest to separate, so the expectation has to hold there.

The depth-weighted form is what the copy-number caller's own `MirroredBaf` column reports: in a clean LOH segment that column *is* the observed germline minor-allele fraction, and it tracks the expression above rather than the cell-weighted one.

The `_dw` suffix on features 10 and 20 to 22 marks this convention. It also acts as the model-compatibility gate: LightGBM validates feature names but not values, so a booster trained against a cell-weighted expectation fails loudly at load time instead of silently scoring against a different feature.

## Classification model

**Algorithm:** Gradient Boosted Decision Trees (GBDT), one model per variant type (SNV, indel).

**Classification threshold:** predict applies one fixed candidate-tier cutoff, `0.96` by default, to every sample. Select a different cutoff with `predict --min-probability-score <FLOAT>`.

**Model handling:** unbalanced datasets are used as-is without up- or down-sampling; `--balance` can enable LightGBM's class balancing.

## Default threshold

The model ranks variants well but its probability is overconfident, so a fixed `0.5` cutoff sits below the point that actually separates the classes and calls too many somatic variants germline.

`--min-probability-score` defaults to a single fixed candidate-tier cutoff of `0.96`, chosen by directly sweeping candidate cutoffs against each validation sample's PASS-filtered prediction and truth call sets and selecting the value with the best average F1 across the cohort. This avoids the extrapolation risk of a per-sample estimate fit on limited data, at the cost of not adapting to any individual sample's characteristics. Each validation sample's F1 at the selected cutoff is close to that sample's own best-achievable F1, indicating the fixed value is not dominated by any single sample. As more validation data becomes available, this fixed value (or the decision to move to a per-sample estimate) should be revisited.

## Tiered germline cutoffs

Variants common in the population are germline with very high probability and vastly outnumber the rest, so they pose a different decision problem from the rare and unannotated variants where the origin is genuinely in doubt.
A single cutoff has to serve both, and the compromise is poor in both directions: loose enough to keep rare somatic calls, it lets common polymorphisms into the somatic set; tight enough to exclude those, it discards rare somatic calls.

Tiering therefore uses two independent values, on two different axes, that are easy to conflate because both are numbers between 0 and 1:

- `--min-popaf` (default raw `INFO/POPAF = 0.001`) decides **which tier** scores a variant: a variant with population allele frequency at or above this value is routed to the population-known tier; everything else, including a missing `POPAF` annotation, uses the candidate tier.
- The population-known tier's own cutoff is a fixed, separate **model probability** of `0.05` -- not currently a CLI option -- applied only after a variant has already been routed there. This lower cutoff calls germline more readily, so those variants must clear a much higher bar of evidence before being released into the somatic set.

Basically: `--min-popaf` answers "is this variant common enough in the population database to get the lenient tier?" (a population-frequency question), while the fixed `0.05` cutoff answers "given the model's own confidence score, is it germline?" (a model-probability question) for variants already in that tier. Raising or lowering one has no bearing on the other.

Developers can inspect `predict --help-advanced` to change the candidate-tier cutoff (`--min-probability-score`) or the population-known POPAF boundary (`--min-popaf`).

This pairs with the `train` option `--train-max-popaf`, which defaults to `0.00001` and removes rows at or
above that population frequency from the *training* set so the model fits the harder population it is
actually asked to decide.
Those variants are still scored at predict time; they are simply judged against their own cutoff.

## Purity and separability

In tumor-only data the only thing distinguishing a germline heterozygote from a somatic mutation is the normal fraction: normal cells contribute alt reads at a germline site and none at a somatic site.
That contribution is proportional to `(1 - purity)`, so the evidence is spent as purity rises.

In balanced diploid the germline mode sits at `0.5` regardless of purity while a clonal somatic variant sits at `purity/2`, so the gap is exactly `(1 - purity)/2` and vanishes in a pure tumor.
Under copy-neutral LOH the confusable pair is instead the germline variant on the lost allele versus the somatic variant, whose gap is `abs(2*purity - 1)/2` -- widest at the extremes and zero at 50% purity.
The two regimes therefore fail at opposite ends of the purity range.

`af_separability_dw` reports that gap in units of sampling noise, which also captures the depth dependence: holding the separation constant as purity rises requires depth to grow as `1/(1 - purity)^2`.
Giving the model the quantity directly lets it learn when to stop trusting allele-fraction geometry and lean on the population-frequency and neighbourhood features instead.

## Classification overrides

After the model scores a variant, two deterministic overrides can reclassify a model-**somatic** call as germline.
Both are disabled by default (their cutoffs default to `1.0`).
The AF override is evaluated first, so a record tripping both is tagged `high_af`.
Each override sets a **standalone** FILTER value (`high_af` or `common_germline`), not the `germline;high_af` / `germline;common_germline` pair used previously: the record is no longer part of the `FILTER=germline` model-call set, and it carries no `GERMLINE_PRED`.

### AF override

A variant the model keeps somatic is reclassified as germline (`FILTER=high_af`) when its tumor allele fraction is at or above the effective cutoff:

```text
effective_cutoff = min(--max-somatic-af, scale * purity)
```

where `purity` comes from the SEG file and `scale` is variant-type specific: `1.0` for SNVs and `1.25` for indels.
The indel distribution has slightly higher AF, so the wider indel cutoff avoids over-filtering true somatic
indels while the tighter SNV cutoff keeps SNV precision high.
The default `--max-somatic-af 1.0` disables the override; set a value below `1.0` to enable it.

See [Per-type AF-override tuning](#per-type-af-override-tuning) for the sweep that motivated the two multipliers.

### POPAF override

A variant the model and AF override both keep somatic is reclassified as germline (`FILTER=common_germline`) when its raw `INFO/POPAF` is at or above the cutoff — a common polymorphism is overwhelmingly germline.
This is an optional hard population-frequency rule, off by default. The model already learns from population frequency; this override lets teams stack an explicit, tunable cutoff on top of the learned score (and is the sole population-frequency stage when the model is trained with `--popaf-blind`).

SNVs use `--common-germline-popaf` and indels use `--indel-common-germline-popaf`, applied independently.
Both default to `1.0` (disabled); set below `1.0` (for example `0.01`, the conventional 1% common-variant line) to enable.
POPAF-missing variants are never tagged and fall through to the model.

## No-call statuses

Beyond the germline/somatic call, the tagged VCF can abstain on a locus it cannot trust, marking it with a no-call or ambiguous status.
A no-call/ambiguous record is excluded from **both** the germline and the somatic (PASS) sets, so `bcftools view -f PASS` drops it.
The status is written to `INFO/GERMLINE_STATUS` on every record and (for statuses that have one) to the `FILTER` column; see [Tagged VCF](file-formats.md#tagged-vcf) for the full token table.

**Precedence** (first match wins, per record):

1. **`nocall_purity`** — the sample's SEG purity is above `0.95`. At near-pure tumor the germline/somatic allele-fraction geometry the model relies on collapses, so every variant in the sample is no-called.
2. **`nocall_cyp2d6` / `nocall_hla`** — the variant falls in the CYP2D6/CYP2D7 locus or the MHC/HLA region. These are decided from the variant's chromosome and position against hardcoded GRCh38, `chr`-prefixed windows (`chr22:42,120,000-42,155,000` and `chr6:29,600,000-33,500,000`).
3. **Input-quality overlay:** `nocall_segment_missing` when no copy-number segment covers the variant (the diploid
   fallback was used), or `ambiguous_cna` when a covered segment's copy number is untrustworthy (no BAF
   evidence, or a homozygous-deletion total copy number of 0); when supplied, this overlay applies to scored
   and unscored records.
4. **Unscored records:** `nocall_coverage` when tumor `DP` is below `--min-tumor-depth`, otherwise `unknown`.
5. **Scored records:** the model call plus the AF/POPAF overrides above.

**Caveats.**

- The region windows are **hardcoded coordinates** for a specific reference build. Confirm them against the deployment reference before relying on the region no-calls; the CYP2D6 window deliberately spans the CYP2D7 paralog and the HLA window covers the classical MHC.
- `ambiguous_cna` indicates that the copy-number SEG input cannot provide reliable allele-specific context at a
  covered locus. The status is assigned when the segment has no allelic (BAF) evidence or a total copy number of
  `0`.
- Because a no-call overrides the model decision, the `tagged_vcf_*` counts (final output) can differ from the `predictions_*` counts (model decisions) for the affected records. This is expected, not a bug.

## Training internals

Training builds a LightGBM classifier from labeled multi-sample data end to end from a single command.
Each sample's VCF and SEG file are parsed into the 31-feature vector, labeled from the `FILTER` value (0 = somatic, 1 = germline), and assembled into a deduplicated feature matrix.
The matrix is split by variant type, and one model is trained per type: rows are shuffled with a `--seed`-seeded RNG, split 60% train / 40% validation, and boosted up to `--num-iterations` (default 1000) with early stopping (`--early-stopping-rounds`, default 50).
Records below `--min-tumor-depth`, records without a recognized `FILTER` label, sex-chromosome variants, and rows with any non-exempt `NaN` cell are removed before boosting.
The dataset is used unbalanced by default.
The SNV booster is written to *\<output-prefix>_snv_model.txt* and the indel booster to *\<output-prefix>_indel_model.txt*.

### Hyperparameters

Fixed in the source (`Model::Train` config string):

| Parameter | Value |
|:----------|:------|
| objective | `binary` (binary cross-entropy) |
| boosting | `gbdt` |
| learning_rate | `0.01` |
| num_leaves | `10` |
| feature_fraction | `0.7` |
| bagging_fraction | `0.7` |
| bagging_freq | `5` |
| min_data_in_leaf | `20` |
| metric | `binary_logloss` |
| boost_from_average | `true` |
| max_bin | `255` |
| num_threads | `12` |
| train/validation split | `60/40` |

## POPAF-blind training

By default the model learns from the raw population-frequency feature `popaf` like any other model feature; `popaf_is_missing` remains in the physical matrix but is excluded from LightGBM in both modes.
The `--popaf-blind` training flag reproduces the earlier design, in which the model was trained *blind* to `popaf`.

### What `--popaf-blind` does

The `popaf` feature is assigned a LightGBM `feature_contri` weight of `0`.
It remains in the matrix — so the model input schema is unchanged and the same feature vector is used at predict time — but the booster is prevented from splitting on it. `popaf_is_missing` is excluded from the model feature set in both modes.
The model then learns germline vs. somatic from the remaining features (allele fraction, depth, copy number, and so on) only, and population frequency is applied outside the model as the deterministic [POPAF override](#popaf-override).

## Per-type AF-override tuning

The AF override is disabled by default (`--max-somatic-af 1.0`). Set `--max-somatic-af` below `1.0` to enable it.
When enabled, it uses per-variant-type multipliers (`scale` = 1.0 for SNVs, 1.25 for indels) rather than a single
global scale.
The indel distribution has slightly higher tumor allele fraction than the SNV distribution, so a wider indel
cutoff avoids reclassifying true somatic indels as germline, while the tighter SNV cutoff keeps SNV precision high.
A single global scale cannot sit at both variant types' optima: on a held-out tumor-only sweep SNV F1 peaks around scale 1.0 to 1.1 and indel F1 peaks around 1.2 to 1.3, and the per-type multipliers give the best overall F1.

## Performance and cost notes

- **Prediction is multithreaded.** `predict --threads` partitions the input VCF into genomic regions scored in parallel and concatenated in order; decompressed VCF content is identical regardless of thread count. Use `--threads 0` to use all available hardware threads on large whole-genome inputs.
- **Training is single-invocation and CPU-bound.** LightGBM boosting uses a fixed internal thread count (see [Hyperparameters](#hyperparameters)); training time scales with the number of labeled variants and iterations (`--num-iterations`, with early stopping).
- **Memory** scales with the number of variants held in the feature matrix. Budget at least 16 GiB for whole-genome cohorts.
