#pragma once

#include <string>

namespace xoos::alignment_metrics {

// This file contains the names of the metrics produced by the alignment_metrics module.

static const std::string kNameMetricName = "metric_name";
static const std::string kNameValue = "value";

// Errors by Cluster Size
static const std::string kNameClusterSize = "cluster_size";
static const std::string kNameTotalBases = "total_bases";
static const std::string kNameTotalAlignedBases = "total_aligned_bases";
static const std::string kNameOverallPhred = "overall_phred";
static const std::string kNameTotalBasesMixedStrandCluster = "total_bases_mixed_strand_cluster";
static const std::string kNamePctTotalBasesMixedStrandCluster = "pct_total_bases_mixed_strand_cluster";
static const std::string kNameTotalBasesForwardCluster = "total_bases_forward_cluster";
static const std::string kNamePctTotalBasesForwardCluster = "pct_total_bases_forward_cluster";
static const std::string kNameTotalBasesReverseCluster = "total_bases_reverse_cluster";
static const std::string kNamePctTotalBasesReverseCluster = "pct_total_bases_reverse_cluster";
static const std::string kNameTotalBasesFullCluster = "total_bases_full_cluster";
static const std::string kNamePctTotalBasesFullCluster = "pct_total_bases_full_cluster";
static const std::string kNameTotalBasesMixedFullAndPartialCluster = "total_bases_mixed_full_and_partial_cluster";
static const std::string kNamePctTotalBasesMixedFullAndPartialCluster =
    "pct_total_bases_mixed_full_and_partial_cluster";
static const std::string kNameTotalBasesPartialCluster = "total_bases_partial_cluster";
static const std::string kNamePctTotalBasesPartialCluster = "pct_total_bases_partial_cluster";
static const std::string kNameSubstitutionsTotal = "substitutions_total";
static const std::string kNameSubstitutionsPhred = "substitutions_phred";
static const std::string kNameSubstitutionsMixedStrandCluster = "substitutions_mixed_strand_cluster";
static const std::string kNamePctSubstitutionsMixedStrandCluster = "pct_substitutions_mixed_strand_cluster";
static const std::string kNameSubstitutionsForwardCluster = "substitutions_forward_cluster";
static const std::string kNamePctSubstitutionsForwardCluster = "pct_substitutions_forward_cluster";
static const std::string kNameSubstitutionsReverseCluster = "substitutions_reverse_cluster";
static const std::string kNamePctSubstitutionsReverseCluster = "pct_substitutions_reverse_cluster";
static const std::string kNameSubstitutionsFullCluster = "substitutions_full_cluster";
static const std::string kNamePctSubstitutionsFullCluster = "pct_substitutions_full_cluster";
static const std::string kNameSubstitutionsMixedFullAndPartialCluster = "substitutions_mixed_full_and_partial_cluster";
static const std::string kNamePctSubstitutionsMixedFullAndPartialCluster =
    "pct_substitutions_mixed_full_and_partial_cluster";
static const std::string kNameSubstitutionsPartialCluster = "substitutions_partial_cluster";
static const std::string kNamePctSubstitutionsPartialCluster = "pct_substitutions_partial_cluster";
static const std::string kNameInsertionsTotal = "insertions_total";
static const std::string kNameInsertionsPhred = "insertions_phred";
static const std::string kNameInsertionsMixedStrandCluster = "insertions_mixed_strand_cluster";
static const std::string kNamePctInsertionsMixedStrandCluster = "pct_insertions_mixed_strand_cluster";
static const std::string kNameInsertionsForwardCluster = "insertions_forward_cluster";
static const std::string kNamePctInsertionsForwardCluster = "pct_insertions_forward_cluster";
static const std::string kNameInsertionsReverseCluster = "insertions_reverse_cluster";
static const std::string kNamePctInsertionsReverseCluster = "pct_insertions_reverse_cluster";
static const std::string kNameInsertionsFullCluster = "insertions_full_cluster";
static const std::string kNamePctInsertionsFullCluster = "pct_insertions_full_cluster";
static const std::string kNameInsertionsMixedFullAndPartialCluster = "insertions_mixed_full_and_partial_cluster";
static const std::string kNamePctInsertionsMixedFullAndPartialCluster = "pct_insertions_mixed_full_and_partial_cluster";
static const std::string kNameInsertionsPartialCluster = "insertions_partial_cluster";
static const std::string kNamePctInsertionsPartialCluster = "pct_insertions_partial_cluster";
static const std::string kNameDeletionsTotal = "deletions_total";
static const std::string kNameDeletionsPhred = "deletions_phred";
static const std::string kNameDeletionsMixedStrandCluster = "deletions_mixed_strand_cluster";
static const std::string kNamePctDeletionsMixedStrandCluster = "pct_deletions_mixed_strand_cluster";
static const std::string kNameDeletionsForwardCluster = "deletions_forward_cluster";
static const std::string kNamePctDeletionsForwardCluster = "pct_deletions_forward_cluster";
static const std::string kNameDeletionsReverseCluster = "deletions_reverse_cluster";
static const std::string kNamePctDeletionsReverseCluster = "pct_deletions_reverse_cluster";
static const std::string kNameDeletionsFullCluster = "deletions_full_cluster";
static const std::string kNamePctDeletionsFullCluster = "pct_deletions_full_cluster";
static const std::string kNameDeletionsMixedFullAndPartialCluster = "deletions_mixed_full_and_partial_cluster";
static const std::string kNamePctDeletionsMixedFullAndPartialCluster = "pct_deletions_mixed_full_and_partial_cluster";
static const std::string kNameDeletionsPartialCluster = "deletions_partial_cluster";
static const std::string kNamePctDeletionsPartialCluster = "pct_deletions_partial_cluster";

// Errors by Substitution Type
static const std::string kNameSubstitutionType = "substitution_type";
static const std::string kNameSubstitutionTypeAny = "any";

// Errors by Read Type
static const std::string kNameDenominator = "denominator";
static const std::string kNamePercentage = "percentage";

// Base-level accuracy summary
static const std::string kNameType = "type";
static const std::string kNameCount = "count";
static const std::string kNamePhred = "phred";
static const std::string kNameSubstitutions = "substitutions";
static const std::string kNameInsertions = "insertions";
static const std::string kNameInsertionEvents = "insertion_events";
static const std::string kNameDeletionEvents = "deletion_events";
static const std::string kNameIndelEvents = "indel_events";
static const std::string kNameInsertedBases = "inserted_bases";
static const std::string kNameDeletedBases = "deleted_bases";
static const std::string kNameIndelBases = "indel_bases";
static const std::string kNameDeletions = "deletions";
static const std::string kNameIndel = "indel";
static const std::string kNameAllErrors = "all_errors";
static const std::string kNameAllErrorsWithIndelEvents = "all_errors_with_indel_events";
static const std::string kNameAllErrorsWithIndelBases = "all_errors_with_indel_bases";
static const std::string kNamePositionsSkippedDueToLowDepth = "positions_skipped_due_to_low_depth";
static const std::string kNamePositionsSkippedDueToHighAltAlleleFraction =
    "positions_skipped_due_to_high_alt_allele_fraction";

// Errors by Read Strand
static const std::string kNameTotalBasesForwardStrand = "total_bases_forward_strand";
static const std::string kNameTotalBasesReverseStrand = "total_bases_reverse_strand";
static const std::string kNameSubstitutionsForwardStrand = "substitutions_forward_strand";
static const std::string kNameSubstitutionsReverseStrand = "substitutions_reverse_strand";
static const std::string kNameInsertionsForwardStrand = "insertions_forward_strand";
static const std::string kNameInsertionsReverseStrand = "insertions_reverse_strand";
static const std::string kNameDeletionsForwardStrand = "deletions_forward_strand";
static const std::string kNameDeletionsReverseStrand = "deletions_reverse_strand";

// Errors by Read Type
static const std::string kNameTotalBasesFullRead = "total_bases_full_read";
static const std::string kNameTotalBasesPartialRead = "total_bases_partial_read";
static const std::string kNameSubstitutionsFullRead = "substitutions_full_read";
static const std::string kNameSubstitutionsPartialRead = "substitutions_partial_read";
static const std::string kNameInsertionsFullRead = "insertions_full_read";
static const std::string kNameInsertionsPartialRead = "insertions_partial_read";
static const std::string kNameDeletionsFullRead = "deletions_full_read";
static const std::string kNameDeletionsPartialRead = "deletions_partial_read";

// Coverage Histogram
static const std::string kNameCoverage = "coverage";
static const std::string kNameConcordantDuplexCount = "concordant_duplex_position_count";
static const std::string kNameAnyCoverage = "any_coverage";
static const std::string kNamePostFilterPositionCount = "post_filter_position_count";
static const std::string kNamePostFilterCoverage = "post_filter_coverage";
static const std::string kNameConcordantDuplexCoverage = "concordant_duplex_coverage";
static const std::string kNamePositionCount = "position_count";
static const std::string kNameMedian = "median";
static const std::string kNameMean = "mean";
static const std::string kNameMin = "min";
static const std::string kNameMax = "max";
static constexpr std::string_view kNamePercentileX = "percentile_{}";
static constexpr std::string_view kNameRatioXtoYPercentile = "ratio_{}_to_{}_percentile";
static const std::string kNameTotalPositions = "total_positions";
static const std::string kNamePositionsWithCoverage = "positions_with_coverage";
static const std::string kNamePercentageOfPositionsWithCoverage = "percentage_of_positions_with_coverage";
static const std::string kNamePositionsWithNoCoverage = "positions_with_no_coverage";
static const std::string kNamePercentageOfPositionsWithNoCoverage = "percentage_of_positions_with_no_coverage";
static constexpr std::string_view kNamePositionsAtLeastXCoverage = "positions_with_at_least_{}x_coverage";
static constexpr std::string_view kNamePercentageOfPositionsAtLeastXCoverage =
    "percentage_of_positions_with_at_least_{}x_coverage";

// Read Count by Cluster Size
static const std::string kNamePctMixedStrandClusterReads = "pct_mixed_strand_cluster_reads";
static const std::string kNamePctForwardClusterReads = "pct_forward_cluster_reads";
static const std::string kNamePctReverseClusterReads = "pct_reverse_cluster_reads";
static const std::string kNamePctFullClusterReads = "pct_full_cluster_reads";
static const std::string kNamePctMixedFullAndPartialClusterReads = "pct_mixed_full_and_partial_cluster_reads";
static const std::string kNamePctPartialClusterReads = "pct_partial_cluster_reads";

// Read-level Metrics Summary
static const std::string kNameTotalAlignmentsIncludingUnmapped = "alignments_plus_unmapped_reads";
static const std::string kNameTotalReads = "total_reads";
static const std::string kNameMappedReads = "mapped_reads";
static const std::string kNameUnmappedReads = "unmapped_reads";
static const std::string kNameForwardReads = "forward_reads";
static const std::string kNameReverseReads = "reverse_reads";
static const std::string kNameFullLengthReads = "full_length_reads";
static const std::string kNamePartialLengthReads = "partial_length_reads";
static const std::string kNameNoUmiReads = "mapped_no_umi_reads";
static const std::string kNameMixedStrandClusterReads = "mixed_strand_cluster_reads";
static const std::string kNameForwardClusterReads = "forward_cluster_reads";
static const std::string kNameReverseClusterReads = "reverse_cluster_reads";
static const std::string kNameFullClusterReads = "full_cluster_reads";
static const std::string kNameMixedFullAndPartialClusterReads = "mixed_full_and_partial_cluster_reads";
static const std::string kNamePartialClusterReads = "partial_cluster_reads";
static const std::string kNameSupplementaryAlignments = "supplementary_alignments";
static const std::string kNameSecondaryAlignments = "secondary_alignments";
static const std::string kNameDuplicateReads = "duplicate_reads";
static const std::string kNameMapQZeroReads = "mapq_zero_reads";
static const std::string kNameReadsPassingFilter = "reads_passing_filter";
static const std::string kNameForwardReadsPassingFilter = "forward_reads_passing_filter";
static const std::string kNameReverseReadsPassingFilter = "reverse_reads_passing_filter";
static const std::string kNameFullLengthReadsPassingFilter = "full_length_reads_passing_filter";
static const std::string kNamePartialLengthReadsPassingFilter = "partial_length_reads_passing_filter";
static const std::string kNameNoUmiReadsPassingFilter = "no_umi_reads_passing_filter";
static const std::string kNameMixedStrandClusterReadsPassingFilter = "mixed_strand_cluster_reads_passing_filter";
static const std::string kNameForwardClusterReadsPassingFilter = "forward_cluster_reads_passing_filter";
static const std::string kNameReverseClusterReadsPassingFilter = "reverse_cluster_reads_passing_filter";
static const std::string kNameFullClusterReadsPassingFilter = "full_cluster_reads_passing_filter";
static const std::string kNameMixedFullAndPartialClusterReadsPassingFilter =
    "mixed_full_and_partial_cluster_reads_passing_filter";
static const std::string kNamePartialClusterReadsPassingFilter = "partial_cluster_reads_passing_filter";
static const std::string kNameAlignedBases = "aligned_bases";
static const std::string kNameUnmappedBases = "unmapped_bases";
static const std::string kNameGCBases = "gc_bases";
static const std::string kNameSoftClippedBases = "soft_clipped_bases";
static const std::string kNameTotalBasesPassingFilter = "total_bases_passing_filter";
static const std::string kNameAlignedBasesPassingFilter = "aligned_bases_passing_filter";
static const std::string kNameGCBasesPassingFilter = "gc_bases_passing_filter";
static const std::string kNameSoftClippedBasesPassingFilter = "soft_clipped_bases_passing_filter";

// Read Metrics for TE
static const std::string kNameTotalMappedAlignmentsInBam = "total_mapped_alignments_in_bam";
static const std::string kNameOnTargetAlignments = "on_target_alignments";
static const std::string kNameOnTargetReadsPassingFilter = "on_target_reads_passing_filter";

// Read Length Histogram
static const std::string kNameReadLengthExcludingSoftClips = "read_length_excluding_soft_clipped_bases";
static const std::string kNameReadLength = "read_length";
static const std::string kNamePostFilterPartialReadCount = "post_filter_partial_read_count";
static const std::string kNamePostFilterFullReadCount = "post_filter_full_read_count";
static const std::string kNamePostFilterAllReadCount = "post_filter_read_count";
static const std::string kNamePostFilterAllReadsLength = "post_filter_read_length";
static const std::string kNamePostFilterPartialReadsLength = "post_filter_partial_read_length";
static const std::string kNamePostFilterFullReadsLength = "post_filter_full_read_length";
static const std::string kNamePostFilterAllReadsLengthExcludingSoftClips =
    "post_filter_read_length_excluding_soft_clipped_bases";
static const std::string kNamePostFilterPartialReadsLengthExcludingSoftClips =
    "post_filter_partial_read_length_excluding_soft_clipped_bases";
static const std::string kNamePostFilterFullReadsLengthExcludingSoftClips =
    "post_filter_full_read_length_excluding_soft_clipped_bases";

// HP Accuracy
static const std::string kNameHpBase = "hp_base";
static const std::string kNameHpLength = "hp_length";
static const std::string kNameHpCount = "hp_count";
static const std::string kNameHpTotalReads = "total_reads";
static const std::string kNameHpSpanningReads = "spanning_reads";
static const std::string kNameHpPercentageSpanning = "percentage_spanning";
static const std::string kNameHpMeanSpanningCoverage = "mean_spanning_coverage";
static const std::string kNameHpDiscordantReads = "discordant_reads";
static const std::string kNameHpPercentageDiscordant = "percentage_discordant";
static const std::string kNameHpLowQualityReads = "low_quality_reads";
static const std::string kNameHpPercentageLowQuality = "percentage_low_quality";
static const std::string kNameHpEffectiveReads = "effective_reads";
static const std::string kNameHpPercentageEffective = "percentage_effective";
static const std::string kNameHpMeanEffectiveCoverage = "mean_effective_coverage";
static const std::string kNameHpReadsWithInsertion = "effective_reads_with_insertion";
static const std::string kNameHpReadsWithDeletion = "effective_reads_with_deletion";
static const std::string kNameHpReadsWithSubstitution = "effective_reads_with_substitution";
static const std::string kNameHpInsertionErrorRate = "insertion_error_rate";
static const std::string kNameHpDeletionErrorRate = "deletion_error_rate";
static const std::string kNameHpIndelErrorRate = "indel_error_rate";
static const std::string kNameHpAccuracy = "accuracy";

// Coverage uniformity metrics
static const std::string kNameFold80BasePenalty = "fold_80_base_penalty";
static const std::string kNamePctBasesAt05xTo2xMeanCoverage = "pct_bases_at_0.5x_to_2x_mean_coverage";
static const std::string kNameMeanCoverage = "mean_coverage";
static const std::string kNameMedianCoverage = "median_coverage";
static const std::string kNameTargetRegionCount = "target_region_count";
static const std::string kNameTargetRegionsWithNoCoverage = "target_regions_with_no_coverage";

}  // namespace xoos::alignment_metrics
