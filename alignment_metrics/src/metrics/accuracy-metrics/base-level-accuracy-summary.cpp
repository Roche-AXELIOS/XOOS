#include "metrics/accuracy-metrics/base-level-accuracy-summary.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <fmt/format.h>

#include <csv.hpp>

#include <xoos/types/float.h>
#include <xoos/util/math.h>

#include "metrics/metrics-names.h"
#include "util/format-util.h"

namespace xoos::alignment_metrics {

std::vector<std::string> BaseLevelAccuracySummary::GetHeaders() {
  return {kNameType, kNameCount, kNameDenominator, kNamePhred, kNamePercentage};
}

void BaseLevelAccuracySummary::WriteTsv(const fs::path& output_path, const io::Comments& comments) const {
  // The output of base-level-accuracy-summary is a tsv file with the following columns:
  // - Error type
  // - Count of errors of this type
  // - Phred score for this error type (calculated as -10 * log10(error rate))
  // The rows are total_bases, substitutions, insertions, deletions, indels, and all_errors.
  constexpr u8 kDefaultPhredPrecision = 2;
  constexpr u8 kDefaultPercentagePrecision = 5;

  const auto phred_denominator_with_indel_events = total_bases + insertions + deletions;
  const auto phred_denominator_with_indel_bases = total_bases + inserted_bases + deleted_bases;

  const std::string phred_denominator_name_with_indel_events =
      fmt::format("{} + {}", kNameTotalAlignedBases, kNameIndelEvents);
  const std::string phred_denominator_name_with_indel_bases =
      fmt::format("{} + {}", kNameTotalAlignedBases, kNameIndelBases);

  std::ofstream ofstream(output_path);
  auto writer = csv::make_tsv_writer_buffered(ofstream);
  io::WriteTsvComments(writer, comments);
  writer << GetHeaders();
  writer << std::make_tuple(kNameTotalAlignedBases, total_bases, kNotApplicable, kNotApplicable, kNotApplicable);
  // Same substitution count with different denominators used for phred and percentage calculations
  writer << std::make_tuple(
      kNameSubstitutions,
      substitutions,
      phred_denominator_name_with_indel_events,
      ToStringWithPrecision(math::ErrorRateToPhred(static_cast<f64>(substitutions) /
                                                   static_cast<f64>(phred_denominator_with_indel_events)),
                            kDefaultPhredPrecision),
      ToPercentageWithPrecision(substitutions, phred_denominator_with_indel_events, kDefaultPercentagePrecision));
  writer << std::make_tuple(
      kNameSubstitutions,
      substitutions,
      phred_denominator_name_with_indel_bases,
      ToStringWithPrecision(math::ErrorRateToPhred(static_cast<f64>(substitutions) /
                                                   static_cast<f64>(phred_denominator_with_indel_bases)),
                            kDefaultPhredPrecision),
      ToPercentageWithPrecision(substitutions, phred_denominator_with_indel_bases, kDefaultPercentagePrecision));
  writer << std::make_tuple(
      kNameInsertionEvents,
      insertions,
      phred_denominator_name_with_indel_events,
      ToStringWithPrecision(
          math::ErrorRateToPhred(static_cast<f64>(insertions) / static_cast<f64>(phred_denominator_with_indel_events)),
          kDefaultPhredPrecision),
      ToPercentageWithPrecision(insertions, phred_denominator_with_indel_events, kDefaultPercentagePrecision));
  writer << std::make_tuple(
      kNameInsertedBases,
      inserted_bases,
      phred_denominator_name_with_indel_bases,
      ToStringWithPrecision(math::ErrorRateToPhred(static_cast<f64>(inserted_bases) /
                                                   static_cast<f64>(phred_denominator_with_indel_bases)),
                            kDefaultPhredPrecision),
      ToPercentageWithPrecision(inserted_bases, phred_denominator_with_indel_bases, kDefaultPercentagePrecision));
  writer << std::make_tuple(
      kNameDeletionEvents,
      deletions,
      phred_denominator_name_with_indel_events,
      ToStringWithPrecision(
          math::ErrorRateToPhred(static_cast<f64>(deletions) / static_cast<f64>(phred_denominator_with_indel_events)),
          kDefaultPhredPrecision),
      ToPercentageWithPrecision(deletions, phred_denominator_with_indel_events, kDefaultPercentagePrecision));
  writer << std::make_tuple(
      kNameDeletedBases,
      deleted_bases,
      phred_denominator_name_with_indel_bases,
      ToStringWithPrecision(math::ErrorRateToPhred(static_cast<f64>(deleted_bases) /
                                                   static_cast<f64>(phred_denominator_with_indel_bases)),
                            kDefaultPhredPrecision),
      ToPercentageWithPrecision(deleted_bases, phred_denominator_with_indel_bases, kDefaultPercentagePrecision));
  writer << std::make_tuple(
      kNameIndelEvents,
      insertions + deletions,
      phred_denominator_name_with_indel_events,
      ToStringWithPrecision(math::ErrorRateToPhred(static_cast<f64>(insertions + deletions) /
                                                   static_cast<f64>(phred_denominator_with_indel_events)),
                            kDefaultPhredPrecision),
      ToPercentageWithPrecision(
          insertions + deletions, phred_denominator_with_indel_events, kDefaultPercentagePrecision));
  writer << std::make_tuple(
      kNameIndelBases,
      inserted_bases + deleted_bases,
      phred_denominator_name_with_indel_bases,
      ToStringWithPrecision(math::ErrorRateToPhred(static_cast<f64>(inserted_bases + deleted_bases) /
                                                   static_cast<f64>(phred_denominator_with_indel_bases)),
                            kDefaultPhredPrecision),
      ToPercentageWithPrecision(
          inserted_bases + deleted_bases, phred_denominator_with_indel_bases, kDefaultPercentagePrecision));
  writer << std::make_tuple(
      kNameAllErrorsWithIndelEvents,
      substitutions + insertions + deletions,
      phred_denominator_name_with_indel_events,
      ToStringWithPrecision(math::ErrorRateToPhred(static_cast<f64>(substitutions + insertions + deletions) /
                                                   static_cast<f64>(phred_denominator_with_indel_events)),
                            kDefaultPhredPrecision),
      ToPercentageWithPrecision(
          substitutions + insertions + deletions, phred_denominator_with_indel_events, kDefaultPercentagePrecision));
  writer << std::make_tuple(
      kNameAllErrorsWithIndelBases,
      substitutions + inserted_bases + deleted_bases,
      phred_denominator_name_with_indel_bases,
      ToStringWithPrecision(math::ErrorRateToPhred(static_cast<f64>(substitutions + inserted_bases + deleted_bases) /
                                                   static_cast<f64>(phred_denominator_with_indel_bases)),
                            kDefaultPhredPrecision),
      ToPercentageWithPrecision(substitutions + inserted_bases + deleted_bases,
                                phred_denominator_with_indel_bases,
                                kDefaultPercentagePrecision));

  // Output total positions and skipped positions due to low depth or high alt allele fraction
  writer << std::make_tuple(kNameTotalPositions, total_positions, kNotApplicable, kNotApplicable, kNotApplicable);
  writer << std::make_tuple(
      kNamePositionsSkippedDueToLowDepth,
      positions_skipped_due_to_low_depth,
      kNameTotalPositions,
      kNotApplicable,
      ToPercentageWithPrecision(positions_skipped_due_to_low_depth, total_positions, kDefaultPercentagePrecision));
  writer << std::make_tuple(
      kNamePositionsSkippedDueToHighAltAlleleFraction,
      positions_skipped_due_to_high_alt_allele_fraction,
      kNameTotalPositions,
      kNotApplicable,
      ToPercentageWithPrecision(
          positions_skipped_due_to_high_alt_allele_fraction, total_positions, kDefaultPercentagePrecision));
}

BaseLevelAccuracySummary& BaseLevelAccuracySummary::operator+=(const BaseLevelAccuracySummary& other) {
  total_bases += other.total_bases;
  substitutions += other.substitutions;
  insertions += other.insertions;
  inserted_bases += other.inserted_bases;
  deletions += other.deletions;
  deleted_bases += other.deleted_bases;
  total_positions += other.total_positions;
  positions_skipped_due_to_low_depth += other.positions_skipped_due_to_low_depth;
  positions_skipped_due_to_high_alt_allele_fraction += other.positions_skipped_due_to_high_alt_allele_fraction;
  return *this;
}

BaseLevelAccuracySummary& BaseLevelAccuracySummary::operator+=(
    const CountsByErrorTypeAndReadProperty<u32>& combined_stats) {
  total_bases += combined_stats.totals.total;
  substitutions += combined_stats.substitutions.total;
  insertions += combined_stats.insertions.total;
  deletions += combined_stats.deletions.total;
  return *this;
}

}  // namespace xoos::alignment_metrics
