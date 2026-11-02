<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to the XOOS Analysis Pipeline will be documented in this file.

The XOOS Analysis Pipeline consists of two components: the `xoosnf` Nextflow pipeline and the `pipeline_launcher` CLI (`xoos`).

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

- Included `filter_variants` small-variant caller QC metrics from germline (pretrained/retrained), tumor-normal, tumor-only, and mitochondrial workflows in the merged metrics report and optional metrics summary.
- Added `enable_target_enrichment_metrics` and enabled it in the `te_simplex_benchmark`, `te_simplex_umi_benchmark`, `te_duplex_benchmark`, and `te_duplex_umi_benchmark` profiles.
- Added variant phasing for somatic tumor-normal and tumor-only workflows before downstream evaluation, copy-number calling, and germline filtering.
- Added `te_duplex`, `te_duplex_benchmark`, `te_duplex_umi`, and `te_duplex_umi_benchmark` Target Enrichment duplex profiles; UMI profiles enable UMI clustering and benchmark profiles enable HP, GC-bias, and Hybrid Selection metrics.
- Added experimental MSI calling for tumor/normal pairs and tumor-only workflows via `enable_msi_caller`, with `msi_sites_thresholds_tsv` and `msi_background_threshold_json` resource overrides.
- Added `alignment_metrics_additional_regions` and `alignment_metrics_additional_region_beds` for label-stratified alignment metrics over additional BED regions.
- Added somatic tumor-only WGS variant calling (`variant_calling_mode = somatic_to`), including small-variant filtering, tumor-only copy-number calling, optional germline filtering, and the `somatic_wgs_to_duplex` and `somatic_wgs_to_duplex_benchmark` profiles. Added the associated `svc_model_somatic_tumor_only`, `svc_somatic_to_{truth_vcf,eval_bed}_dir`, `cnc_somatic_to_{blocklist_bed,snp_blocklist_bed}`, and `svc_somatic_to_panel_of_normals_{vcf,tbi}` parameters.
- Added indexed panel-of-normals VCF support through `svc_panel_of_normals_vcf` and `svc_panel_of_normals_tbi`.
- Mitochondrial BAM post-processing now consumes Giraffe graph components, supporting static and personalized pangenome graphs.
- Added `dedup_aligner` for selecting the post-read-collapser re-alignment aligner in the consensus deduplication path.
- Added the `germline_wgs_simplex_benchmark` profile for benchmarking germline WGS simplex variant calling.
- Added Somalier-based tumor/normal pair verification with `enable_somatic_tn_pair_verification`, `somatic_tn_pair_mismatch_action`, `somalier_sites_vcf`, and `somatic_tn_pair_concordance_cutoff`.
- Added pre-aligned CRAM input support through the `cram` `file_type`, with reference FASTA/FAI conversion requirements.
- Added germline mitochondrial variant calling through `enable_mitochondrial_variant_caller`, including heteroplasmy detection down to AF >= 1%.
- Added `enable_downsampled_bwa_metrics` and its aligner, read-count, and seed controls.
- Added `cnv_blocklist_bed` support for germline copy-number calling.
- Added guidance for pinning XOOS container images and the `xoosnf` source, running without online lookups, and verifying versions in the MultiQC report.
- Documented the simplex Germline WGS application and its `germline_wgs_simplex` and `germline_wgs_simplex_benchmark` profiles:
  - Use `pb_giraffe` with pangenome personalization by default.
  - Subsample BAMs to 40× before benchmark processing in `germline_wgs_simplex_benchmark`.
  - Disable STR calling because `str_caller` does not support simplex data.
  - Retain partial reads and use XOOS read-collapser MARKDUP instead of Parabricks integrated duplicate marking.
  - Document the `large` and `small` bundled filtering-model families, custom paired SNV/indel model files, fixed `small`-family model retraining, and simplex run-type handling.
- Added `--small_variant_caller_germline_model_type` to select the bundled germline SVC filtering-model family (`large` by default or `small`).
- Added the SBX-D Somatic Tumor-Only WGS application with the `somatic_wgs_to_duplex` profile, tumor-only samplesheet configuration, and downstream small-variant, copy-number, and germline-filtering workflow documentation.
- Added BAM post-processing module (`bam_post_process`) for mitochondrial variant calling and documented the related options:
  - `--enable_mitochondrial_variant_caller` enables germline mitochondrial (chrM) calling.
  - `--enable_giraffe_mitochondrial_bam_post_process` recovers mitochondrial reads from Giraffe-aligned BAMs before calling. When omitted, it is enabled automatically for `pb_giraffe` and `vg_giraffe` when mitochondrial calling is enabled.
  - Added the Giraffe-specific tuning parameters `--giraffe_mitochondrial_bam_post_process_anchor_pad`, `--giraffe_mitochondrial_bam_post_process_window`, `--giraffe_mitochondrial_bam_post_process_chrm_name`, and `--giraffe_mitochondrial_bam_post_process_multimap`, plus documentation for the mitochondrial VCF output location.
  - Documented that the configured mitochondrial contig name defaults to `chrM` and is used consistently by BAM post-processing and mitochondrial variant calling; pre-aligned BAM/CRAM inputs use the filtered Giraffe graph for graph association.

### Changed

- Enabled tumor-fraction estimation by default for the `cfdna_wgs_duplex` and `cfdna_wgs_duplex_benchmark` profiles when a `reporter_vcf` is provided.
- Hybrid Selection metrics now passes `alignment_metrics_min_mapping_quality` to Picard `CollectHsMetrics` and publishes summary, per-target, and per-base coverage outputs.
- Demux 1.1 is no longer supported; the pipeline now uses the current Demux CLI options and canonical metrics headers.
- Added explicit `read_collapser_cluster_by_umi` handling and the `read_collapser_wiggle_room` parameter.
- Germline copy-number calling now maps pipeline run types to the caller's `--library-type` option and fails clearly for unsupported types.
- Copy-number caller invocations preserve the `metrics/`, `tracks/`, and `supplemental/` hierarchy for report and MultiQC consumers.
- Enabled cross-sample contamination estimation and pangenome-personalized `pb_giraffe` defaults in the germline simplex profiles; the benchmark profile subsamples BAMs to 40× and simplex profiles disable STR calling.
- Germline SVC filtering now selects the bundled `large` or `small` model family through `small_variant_caller_germline_model_type`, and retraining uses the unified `germline-wgs` workflow with explicit library type.
- MSI caller invocations now use the caller's embedded default catalogs and thresholds when no overrides are supplied; `msi_sites_thresholds_tsv` and `msi_background_threshold_json` remain available for explicit resource overrides, while the pipeline no longer passes `--min-depth` or `--min-sites`.
- Somatic tumor-normal small-variant calling now defaults to automatic sample-type detection for bundled BWA models.
- Alignment metric keys now include the aligner and use a hierarchical `<stage>/<aligner>/<region>` layout.
- Split FASTQ subsampling into deterministic `seqkit_head` and random `seqkit_random` variants.
- Enabled high-confidence alignment metrics in the cfDNA duplex benchmark profile and renamed the experimental simplex profiles.
- Simplex germline profiles now use the simplex sequencing protocol, retain partial reads, and use XOOS read-collapser MARKDUP.
- Enabled Hybrid Selection metrics in the Target Enrichment simplex benchmark profiles, with `target_bed` fallback when `capture_bed` is absent.
- Updated copy-number caller output naming, read-collapser CLI/output paths (`--bam-input`/`--bed-input` to `--bam`/`--bed`, and `--prefix` to `--output-prefix`), Parabricks read-length handling, and STR caller flags for cross-module consistency.
- Removed the explicit `enable_subsample_bam` setting from the WGS simplex benchmark profile so it inherits the global default.
- `max_parabricks_alignment_read_length` now sets Parabricks `--max-read-length` and `--max-read-length-on-gpu` and defaults to `1024`.
- Extended samplesheet `file_type` to accept CRAM and added `VERIFY_CRAM` validation for reference-checksum (`M5`) tags.
- Parabricks minimap2 now streams FASTQs through a named pipe, removing the command-line length ceiling for many FASTQ shards.
- Alignment metrics now uses the Giraffe linearized index for Giraffe aligners.
- Duplicate marking now uses XOOS read-collapser for `dedup_strategy = markdup`; integrated aligner duplicate marking is disabled.
- Enforced compatibility between `dedup_strategy` and externally deduplicated BAM/CRAM inputs and skipped those inputs in pre-dedup metrics.
- Copy Number Caller WGS invocations now allow the caller to select its bundled mappability track when no override is configured.
- Removed the obsolete `--use_multisample_model` parameter.

### Dependencies

- Updated the Parabricks container to `4.7.1-1`.
- Added `somalier` `0.2.19` for tumor/normal pair verification.
- Updated seqkit to `2.13.0`.

### Fixed

- Normalized somatic tumor-only vcfeval metric names to `vcf_eval/small_variant_pretrained/*`.
- Preserved caller output directories when replacing somatic CNC and MSI VCFs.
- Corrected invalid Q93 simplex pre-dedup accuracy by setting the simplex pre-dedup quality threshold to `0` and reporting `NA` when no bases are evaluated.
- Parabricks Giraffe now emits multimapping alignments needed by copy-number secondary-alignment rescue.
- Samplesheet schema-validation errors are now clean, eager, and reported once.
- Parabricks aligner FIFO writers are cleaned up on early exit, preventing hangs.
- Parabricks minimap2 now emits one sample-level read group with a populated `SM` tag.
- All aligners fail fast when no FASTQ inputs are found.
- Read groups now include `PM=AXELIOS_1` and a structured library name across aligners.
- HaplotypeCaller and Mutect2 indexes are DAG-tracked through producing channels; GATK realigned BAMs receive an explicit indexing step.
- Germline SVC retraining skips samples without matching benchmark resources instead of aborting the pipeline.
- Added regression coverage for nested `metrics/metrics.tsv` discovery and fixed the pipeline completion handler's nested `workflow` binding by capturing the enclosing run success flag.
- Germline copy-number results now use BAM `SM` names consistently and reheader CNC VCFs to samplesheet names.

### Known issues

- The `enable_bam_reheader` param may misalign the ID value of the RG in the BAM header from the RG tag in the reads causing downstream errors and empty variant call sets. Recommended to not use.
- When `enable_giraffe_mitochondrial_bam_post_process` is enabled for input BAMs whose Giraffe alignments were generated with a personalized pangenome, the pipeline passes incorrect graph resources to BAM post-processing and the stage fails.
- Variant Phasing for `somatic tumor-only WGS` and `somatic tumor-normal WGS` requires Mutect BAMs to be merged back together which may increase runtime.
- The `small_variant_caller` for germline mitochondrial variant calling can take in a blocklist BED file. This file is currently not passed by the pipeline and may result in drop in variant calling performance for germline mitochondrial variant calling.

## [1.1.0]

### Added

- Added SBX-D Germline WGS, SBX-D Somatic Tumor/Normal WGS, SBX-D cfDNA WGS.

<!-- Version comparison links -->
[1.2.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.1.0...1.2.0
[1.1.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.1.0
