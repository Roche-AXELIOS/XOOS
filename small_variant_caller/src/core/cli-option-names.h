#pragma once

namespace xoos::svc::cli_opt_name {
// common CLI option names applicable to multiple applications
constexpr auto* kWarnAsError = "--warn-as-error";
constexpr auto* kThreads = "--threads";
constexpr auto* kConfig = "--config";
constexpr auto* kOutputFile = "--output-file";
constexpr auto* kTargetRegions = "--target-regions";
constexpr auto* kGenome = "--genome";
constexpr auto* kVcfInput = "--vcf-input";
constexpr auto* kLeftPad = "--left-pad";
constexpr auto* kRightPad = "--right-pad";
constexpr auto* kCollapseDistance = "--collapse-distance";
constexpr auto* kNormalizeFeatures = "--normalize-features";
constexpr auto* kBamInput = "--bam-input";
constexpr auto* kTumorSampleName = "--tumor-sample-name";

// CLI option names for compute_bam_features
constexpr auto* kSkipVariantsVcf = "--skip-variants-vcf";
constexpr auto* kMaxVariantsPerRead = "--max-variants-per-read";
constexpr auto* kMaxVariantsPerReadNormalized = "--max-variants-per-read-normalized";
constexpr auto* kMinMapq = "--min-mapq";
constexpr auto* kMinBq = "--min-bq";
constexpr auto* kMinDist = "--min-dist";
constexpr auto* kMinFamilySize = "--min-family-size";
constexpr auto* kFilterHomopolymer = "--filter-homopolymer";
constexpr auto* kMinHomopolymerLength = "--min-homopolymer-length";
constexpr auto* kSequencingProtocol = "--sequencing-protocol";
constexpr auto* kDecodeYc = "--decode-yc";
constexpr auto* kMinBaseType = "--min-base-type";
constexpr auto* kMaxRegionSizePerThread = "--max-region-size-per-thread";

// CLI option names for compute_vcf_features
constexpr auto* kInterestRegions = "--interest-regions";
constexpr auto* kPopAfVcf = "--pop-af-vcf";
constexpr auto* kPonVcf = "--pon-vcf";
constexpr auto* kCatalogVcf = "--catalog-vcf";
constexpr auto* kOutputBed = "--output-bed";

// CLI option names for train_model
constexpr auto* kPosBamFeatures = "--positive-bam-features";
constexpr auto* kOutputTrainingDataTsv = "--output-training-data-tsv";
constexpr auto* kPosVcfFeatures = "--positive-vcf-features";
constexpr auto* kNegBamFeatures = "--negative-bam-features";
constexpr auto* kNegVcfFeatures = "--negative-vcf-features";
constexpr auto* kTruthVcfs = "--truth-vcfs";
constexpr auto* kSomaticTruthVcfs = "--somatic-truth-vcfs";
constexpr auto* kGermlineTruthVcfs = "--germline-truth-vcfs";
constexpr auto* kSomaticPositiveRatio = "--somatic-positive-ratio";
constexpr auto* kGermlinePositiveRatio = "--germline-positive-ratio";
constexpr auto* kNegativeRatio = "--negative-ratio";
constexpr auto* kMinSupport = "--min-support";
// Shared BED blocklist option name used by train_model and filter_variants.
constexpr auto* kBlocklistBed = "--blocklist-bed";
constexpr auto* kMaxWeightedScore = "--max-weighted-score";
constexpr auto* kIterations = "--iterations";

// germline-wgs specific CLI option names for train_model
constexpr auto* kSnvIterations = "--snv-iterations";
constexpr auto* kIndelIterations = "--indel-iterations";
constexpr auto* kSnvModelOutput = "--snv-output-file";
constexpr auto* kIndelModelOutput = "--indel-output-file";
constexpr auto* kSnvOutputTrainingDataTsv = "--snv-output-training-data-tsv";
constexpr auto* kIndelOutputTrainingDataTsv = "--indel-output-training-data-tsv";

// CLI option names for filter_variants
constexpr auto* kModel = "--model";
constexpr auto* kOutputDir = "--output-dir";
constexpr auto* kVcfOutput = "--vcf-output";
constexpr auto* kOutputMetrics = "--output-metrics";
constexpr auto* kMaxBamRegionSizePerThread = "--max-bam-region-size-per-thread";
constexpr auto* kMaxVcfRegionSizePerThread = "--max-vcf-region-size-per-thread";
constexpr auto* kSdChrName = "--sd-chr-name";
constexpr auto* kParBedX = "--par-bed-x";
constexpr auto* kParBedY = "--par-bed-y";
constexpr auto* kOutputShapValueTsv = "--output-shap-value-tsv";
constexpr auto* kOutputSnvShapValueTsv = "--output-snv-shap-value-tsv";
constexpr auto* kOutputIndelShapValueTsv = "--output-indel-shap-value-tsv";

// germline-wgs specific CLI option names shared by train_model and filter_variants
constexpr auto* kLibraryType = "--library-type";
constexpr auto* kModelType = "--model-type";

// germline-wgs specific CLI option names for filter_variants
constexpr auto* kSnvModel = "--snv-model";
constexpr auto* kIndelModel = "--indel-model";

// tumor-only-te specific CLI option names for filter_variants
constexpr auto* kBlocklist = "--blocklist";
constexpr auto* kHotspotVcf = "--hotspot-vcf";
constexpr auto* kForcecallBed = "--forcecall-bed";
constexpr auto* kMinAltCounts = "--min-alt-counts";
constexpr auto* kMinAf = "--min-af";
constexpr auto* kMinPhasedAf = "--min-phased-af";
constexpr auto* kMaxPhasedAf = "--max-phased-af";
constexpr auto* kMinWeightedCounts = "--min-weighted-counts";
constexpr auto* kHotspotMinWeightedCounts = "--hotspot-min-weighted-counts";
constexpr auto* kMinMlScore = "--min-ml-score";
constexpr auto* kHotspotMinMlScore = "--hotspot-min-ml-score";
constexpr auto* kPhased = "--phased";

// tumor-normal-wgs specific CLI option names for filter_variants
constexpr auto* kSampleType = "--sample-type";
constexpr auto* kAligner = "--aligner";
constexpr auto* kSnvMinMlScore = "--snv-min-ml-score";
constexpr auto* kIndelMinMlScore = "--indel-min-ml-score";
constexpr auto* kMinTumorSupport = "--min-tumor-support";
constexpr auto* kMaxNormalSupport = "--max-normal-support";
constexpr auto* kMinTumorAf = "--min-tumor-af";
constexpr auto* kMinDpRatio = "--min-dp-ratio";
constexpr auto* kMaxIndelSize = "--max-indel-size";
constexpr auto* kDuplexLowbq = "--duplex-lowbq";
// Relaxed thresholds applied to variants present in the catalog VCF (--catalog-vcf)
constexpr auto* kCatalogMinTumorAf = "--catalog-min-tumor-af";
constexpr auto* kCatalogSnvMinMlScore = "--catalog-snv-min-ml-score";
constexpr auto* kCatalogIndelMinMlScore = "--catalog-indel-min-ml-score";

// germline-mitochondria specific CLI option names for filter_variants
constexpr auto* kMinDp = "--min-dp";
constexpr auto* kMinHomAltAf = "--min-hom-alt-af";

// CLI option names for generate_pon
constexpr auto* kFeatureList = "--feature-list";
constexpr auto* kMinDuplexAf = "--min-duplex-af";
constexpr auto* kMaxDuplexAf = "--max-duplex-af";
constexpr auto* kMinMapqMean = "--min-mapq-mean";

// CLI option names for phase_variants
constexpr auto* kMaxPhaseDistance = "--max-phase-distance";
constexpr auto* kMinSharedReads = "--min-shared-reads";
constexpr auto* kMinCisReads = "--min-cis-reads";
constexpr auto* kMinCisFraction = "--min-cis-fraction";
constexpr auto* kReportMnv = "--report-mnv";
constexpr auto* kMergeMnv = "--merge-mnv";
constexpr auto* kMergeComplex = "--merge-complex";
constexpr auto* kKeepComponents = "--keep-components";
constexpr auto* kIncludeFiltered = "--include-filtered";
constexpr auto* kRescueFiltered = "--rescue-filtered";
constexpr auto* kSkipInsSnv = "--skip-ins-snv";
// Suffixed with `Opt` to avoid shadowing the xoos::svc::kMinBaseQuality config default (cpp:M23_029).
constexpr auto* kMinBaseQualityOpt = "--min-base-quality";

// CLI option names for detect_sample_type
constexpr auto* kTargetBed = "--target-bed";
constexpr auto* kMaxAf = "--max-af";
constexpr auto* kMaxFfpePhred = "--max-ffpe-phred";
constexpr auto* kSampleTypeMaxAf = "--sample-type-max-af";
constexpr auto* kSampleTypeMaxFfpePhred = "--sample-type-max-ffpe-phred";

// CLI option names for features_to_vcf
constexpr auto* kInputFile = "--input-file";
constexpr auto* kSampleName = "--sample-name";

}  // namespace xoos::svc::cli_opt_name
