#include "metrics/accuracy-metrics/errors-by-read-type.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <csv.hpp>

#include "metrics/metrics-names.h"
#include "util/format-util.h"

namespace xoos::alignment_metrics {

std::vector<std::string> ErrorsByReadType::GetHeaders() {
  return {kNameType, kNameCount, kNameDenominator, kNamePercentage};
}

void ErrorsByReadType::WriteTsv(const fs::path& output_path, const io::Comments& comments) const {
  std::ofstream ofstream(output_path);
  auto writer = csv::make_tsv_writer_buffered(ofstream);
  io::WriteTsvComments(writer, comments);
  writer << GetHeaders();
  writer << std::make_tuple(kNameSubstitutionsTotal, substitutions_total, kNotApplicable, kNotApplicable);
  if (dataset_metadata.has_read_type_info) {
    writer << std::make_tuple(
        kNameSubstitutionsFullRead,
        substitutions_full_read,
        kNameSubstitutionsTotal,
        ToPercentageWithPrecision(static_cast<f64>(substitutions_full_read), static_cast<f64>(substitutions_total), 2));
    writer << std::make_tuple(
        kNameSubstitutionsPartialRead,
        substitutions_partial_read,
        kNameSubstitutionsTotal,
        ToPercentageWithPrecision(
            static_cast<f64>(substitutions_partial_read), static_cast<f64>(substitutions_total), 2));
  }
  if (dataset_metadata.has_strand_info) {
    writer << std::make_tuple(
        kNameSubstitutionsForwardStrand,
        substitutions_forward_read,
        kNameSubstitutionsTotal,
        ToPercentageWithPrecision(
            static_cast<f64>(substitutions_forward_read), static_cast<f64>(substitutions_total), 2));
    writer << std::make_tuple(
        kNameSubstitutionsReverseStrand,
        substitutions_reverse_read,
        kNameSubstitutionsTotal,
        ToPercentageWithPrecision(
            static_cast<f64>(substitutions_reverse_read), static_cast<f64>(substitutions_total), 2));
  }
  writer << std::make_tuple(kNameInsertionsTotal, insertions_total, kNotApplicable, kNotApplicable);
  if (dataset_metadata.has_read_type_info) {
    writer << std::make_tuple(
        kNameInsertionsFullRead,
        insertions_full_read,
        kNameInsertionsTotal,
        ToPercentageWithPrecision(static_cast<f64>(insertions_full_read), static_cast<f64>(insertions_total), 2));
    writer << std::make_tuple(
        kNameInsertionsPartialRead,
        insertions_partial_read,
        kNameInsertionsTotal,
        ToPercentageWithPrecision(static_cast<f64>(insertions_partial_read), static_cast<f64>(insertions_total), 2));
  }
  if (dataset_metadata.has_strand_info) {
    writer << std::make_tuple(
        kNameInsertionsForwardStrand,
        insertions_forward_read,
        kNameInsertionsTotal,
        ToPercentageWithPrecision(static_cast<f64>(insertions_forward_read), static_cast<f64>(insertions_total), 2));
    writer << std::make_tuple(
        kNameInsertionsReverseStrand,
        insertions_reverse_read,
        kNameInsertionsTotal,
        ToPercentageWithPrecision(static_cast<f64>(insertions_reverse_read), static_cast<f64>(insertions_total), 2));
  }
  writer << std::make_tuple(kNameDeletionsTotal, deletions_total, kNotApplicable, kNotApplicable);
  if (dataset_metadata.has_read_type_info) {
    writer << std::make_tuple(
        kNameDeletionsFullRead,
        deletions_full_read,
        kNameDeletionsTotal,
        ToPercentageWithPrecision(static_cast<f64>(deletions_full_read), static_cast<f64>(deletions_total), 2));
    writer << std::make_tuple(
        kNameDeletionsPartialRead,
        deletions_partial_read,
        kNameDeletionsTotal,
        ToPercentageWithPrecision(static_cast<f64>(deletions_partial_read), static_cast<f64>(deletions_total), 2));
  }
  if (dataset_metadata.has_strand_info) {
    writer << std::make_tuple(
        kNameDeletionsForwardStrand,
        deletions_forward_read,
        kNameDeletionsTotal,
        ToPercentageWithPrecision(static_cast<f64>(deletions_forward_read), static_cast<f64>(deletions_total), 2));
    writer << std::make_tuple(
        kNameDeletionsReverseStrand,
        deletions_reverse_read,
        kNameDeletionsTotal,
        ToPercentageWithPrecision(static_cast<f64>(deletions_reverse_read), static_cast<f64>(deletions_total), 2));
  }
}

ErrorsByReadType& ErrorsByReadType::operator+=(const ErrorsByReadType& other) {
  substitutions_total += other.substitutions_total;
  substitutions_forward_read += other.substitutions_forward_read;
  substitutions_reverse_read += other.substitutions_reverse_read;
  substitutions_full_read += other.substitutions_full_read;
  substitutions_partial_read += other.substitutions_partial_read;
  insertions_total += other.insertions_total;
  insertions_forward_read += other.insertions_forward_read;
  insertions_reverse_read += other.insertions_reverse_read;
  insertions_full_read += other.insertions_full_read;
  insertions_partial_read += other.insertions_partial_read;
  deletions_total += other.deletions_total;
  deletions_forward_read += other.deletions_forward_read;
  deletions_reverse_read += other.deletions_reverse_read;
  deletions_full_read += other.deletions_full_read;
  deletions_partial_read += other.deletions_partial_read;

  return *this;
}

ErrorsByReadType& ErrorsByReadType::operator+=(const CountsByErrorTypeAndReadProperty<u32>& combined_stats) {
  substitutions_total += combined_stats.substitutions.total;
  substitutions_forward_read += combined_stats.substitutions.forward_read;
  substitutions_reverse_read += combined_stats.substitutions.reverse_read;
  substitutions_full_read += combined_stats.substitutions.full_read;
  substitutions_partial_read += combined_stats.substitutions.partial_read;
  insertions_total += combined_stats.insertions.total;
  insertions_forward_read += combined_stats.insertions.forward_read;
  insertions_reverse_read += combined_stats.insertions.reverse_read;
  insertions_full_read += combined_stats.insertions.full_read;
  insertions_partial_read += combined_stats.insertions.partial_read;
  deletions_total += combined_stats.deletions.total;
  deletions_forward_read += combined_stats.deletions.forward_read;
  deletions_reverse_read += combined_stats.deletions.reverse_read;
  deletions_full_read += combined_stats.deletions.full_read;
  deletions_partial_read += combined_stats.deletions.partial_read;
  return *this;
}

}  // namespace xoos::alignment_metrics
