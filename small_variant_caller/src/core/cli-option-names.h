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
constexpr auto* kOutputBed = "--output-bed";

// CLI option names for train_model
constexpr auto* kPosBamFeatures = "--positive-bam-features";
constexpr auto* kOutputTrainingDataTsv = "--output-training-data-tsv";
constexpr auto* kPosVcfFeatures = "--positive-vcf-features";
constexpr auto* kNegBamFeatures = "--negative-bam-features";
constexpr auto* kNegVcfFeatures = "--negative-vcf-features";
constexpr auto* kTruthVcfs = "--truth-vcfs";
constexpr auto* kBlocklistBed = "--blocklist-bed";
constexpr auto* kMaxWeightedScore = "--max-weighted-score";
constexpr auto* kIterations = "--iterations";

// germline and germline-multi-sample specific CLI option names for train_model
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
constexpr auto* kMaxBamRegionSizePerThread = "--max-bam-region-size-per-thread";
constexpr auto* kMaxVcfRegionSizePerThread = "--max-vcf-region-size-per-thread";
constexpr auto* kOutputVcfBufferSize = "--output-vcf-buffer-size";
constexpr auto* kSdChrName = "--sd-chr-name";
constexpr auto* kParBedX = "--par-bed-x";
constexpr auto* kParBedY = "--par-bed-y";

// germline and germline-multi-sample specific CLI option names for filter_variants
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
constexpr auto* kSnvMinMlScore = "--snv-min-ml-score";
constexpr auto* kIndelMinMlScore = "--indel-min-ml-score";
constexpr auto* kMinTumorSupport = "--min-tumor-support";
}  // namespace xoos::svc::cli_opt_name
