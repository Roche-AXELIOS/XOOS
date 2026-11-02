#pragma once

namespace xoos::cnc {

// Output subfolders.
//   - root         : primary deliverables (callset VCF/SEG) and the ranked purity/ploidy grid
//   - metrics/     : QC and statistics (metrics, chosen purity/ploidy, likelihood analyses)
//   - tracks/      : visualization (BigWig tracks, IGV session)
//   - supplemental/: raw per-interval tracks, coverage, intermediate segments, debug artifacts
constexpr auto* kMetricsSubdir = "metrics";
constexpr auto* kTracksSubdir = "tracks";
constexpr auto* kSupplementalSubdir = "supplemental";

// likelihood output — primary deliverables, kept in the root
constexpr auto* kDefaultSomaticCNCallsetVcfOutput = "somatic_copy_number_callset.vcf.gz";
constexpr auto* kDefaultGermlineCNCallsetVcfOutput = "germline_copy_number_callset.vcf.gz";
constexpr auto* kDefaultSomaticCNCallsetSegOutput = "somatic_copy_number_callset.seg";
constexpr auto* kDefaultGermlineCNCallsetSegOutput = "germline_copy_number_callset.seg";

// visualization tracks
constexpr auto* kDefaultBafBwOutput = "tracks/bafs.bw";
constexpr auto* kDefaultLogRsBwOutput = "tracks/log_ratios.bw";
constexpr auto* kDefaultIgvXmlOutput = "tracks/igv_visualization.xml";

// raw per-interval tracks and intermediate segments — supplemental to the callset
constexpr auto* kDefaultBafOutput = "supplemental/bafs.bed";
constexpr auto* kDefaultBafSegOutput = "supplemental/baf_segments.seg";
constexpr auto* kDefaultLogRsOutput = "supplemental/log_ratios.bed";
constexpr auto* kDefaultLogRsSegOutput = "supplemental/log_ratio_segments.seg";

// coverage output
constexpr auto* kDefaultCoverageOutput = "supplemental/coverage.bed";
constexpr auto* kDefaultCorrectedCoverageOutput = "supplemental/corrected_coverage.bed";
constexpr auto* kDefaultNormalCoverageOutput = "supplemental/normal_coverage.bed";
constexpr auto* kDefaultNormalCorrectedCoverageOutput = "supplemental/normal_corrected_coverage.bed";
constexpr auto* kDefaultTumorCoverageOutput = "supplemental/tumor_coverage.bed";
constexpr auto* kDefaultTumorCorrectedCoverageOutput = "supplemental/tumor_corrected_coverage.bed";

// augmented intervals output
constexpr auto* kDefaultAugmentedIntervalsOutput = "supplemental/augmented_intervals.bed";

// purity/ploidy search solutions and per-segment likelihood analyses
constexpr auto* kDefaultPurityPloidyGridOutput = "supplemental/purity_ploidy_grid.tsv";
// chosen purity/ploidy solution — QC/statistics
constexpr auto* kDefaultPurityPloidyOutput = "metrics/purity_ploidy.tsv";
// Opt-in per-segment likelihood diagnostics, written per purity/ploidy search pass (see --write-segment-likelihood).
constexpr auto* kDefaultSegmentLikelihoodDiploidBiasOutput = "supplemental/segment_likelihood.diploid_bias.tsv";
constexpr auto* kDefaultSegmentLikelihoodUniformOutput = "supplemental/segment_likelihood.uniform.tsv";
constexpr auto* kDefaultMapQOutput = "supplemental/mapping_quality.bed";
constexpr auto* kDefaultMetricsOutput = "metrics/metrics.tsv";
constexpr auto* kDefaultTaskflowGraphOutput = "supplemental/taskflow.dot";

// seg to vcf output — intermediate, supplemental to the callset
constexpr auto* kDefaultVcfOutput = "supplemental/segments.vcf";

// segmentation output — intermediate, supplemental to the callset
constexpr auto* kDefaultSegmentationSegOutput = "supplemental/segments.seg";

}  // namespace xoos::cnc
