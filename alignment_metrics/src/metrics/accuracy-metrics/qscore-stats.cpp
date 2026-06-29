#include "metrics/accuracy-metrics/qscore-stats.h"

#include <string>

#include <csv.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/math.h>

#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

QscoreStats::QscoreStats(const DatasetMetadata& dataset_metadata)
    : mismatches_by_qscore(kMaxQscoreValue + 1), dataset_metadata(dataset_metadata) {
}

vec<std::string> QscoreStats::GetHeaders() const {
  if (dataset_metadata.has_cluster_info) {
    return {"q_score",
            "total_bases",
            "mismatches",
            "empirical_phred_score",
            "forward_cluster_bases",
            "forward_cluster_mismatches",
            "forward_cluster_empirical_phred_score",
            "reverse_cluster_bases",
            "reverse_cluster_mismatches",
            "reverse_cluster_empirical_phred_score",
            "mixed_strand_cluster_bases",
            "mixed_strand_cluster_mismatches",
            "mixed_strand_cluster_empirical_phred_score"};
  }
  return {"q_score", "total_bases", "mismatches", "empirical_phred_score"};
}

void QscoreStats::WriteTsv(const fs::path& output_path, const io::Comments& comments) const {
  std::ofstream ofstream(output_path);
  auto writer = csv::make_tsv_writer_buffered(ofstream);
  io::WriteTsvComments(writer, comments);
  writer << GetHeaders();
  for (size_t qscore = 0; qscore < mismatches_by_qscore.size(); ++qscore) {
    if (mismatches_by_qscore[qscore].total_bases == 0) {
      continue;
    }
    const f64 error_rate = static_cast<f64>(mismatches_by_qscore[qscore].total_mismatches) /
                           static_cast<f64>(mismatches_by_qscore[qscore].total_bases);
    if (dataset_metadata.has_cluster_info) {
      const f64 forward_cluster_error_rate = static_cast<f64>(mismatches_by_qscore[qscore].forward_cluster_mismatches) /
                                             static_cast<f64>(mismatches_by_qscore[qscore].forward_cluster_bases);
      const f64 reverse_cluster_error_rate = static_cast<f64>(mismatches_by_qscore[qscore].reverse_cluster_mismatches) /
                                             static_cast<f64>(mismatches_by_qscore[qscore].reverse_cluster_bases);
      const f64 mixed_strand_cluster_error_rate =
          static_cast<f64>(mismatches_by_qscore[qscore].mixed_strand_cluster_mismatches) /
          static_cast<f64>(mismatches_by_qscore[qscore].mixed_strand_cluster_bases);
      writer << std::make_tuple(qscore,
                                mismatches_by_qscore[qscore].total_bases,
                                mismatches_by_qscore[qscore].total_mismatches,
                                math::ErrorRateToPhred(error_rate),
                                mismatches_by_qscore[qscore].forward_cluster_bases,
                                mismatches_by_qscore[qscore].forward_cluster_mismatches,
                                math::ErrorRateToPhred(forward_cluster_error_rate),
                                mismatches_by_qscore[qscore].reverse_cluster_bases,
                                mismatches_by_qscore[qscore].reverse_cluster_mismatches,
                                math::ErrorRateToPhred(reverse_cluster_error_rate),
                                mismatches_by_qscore[qscore].mixed_strand_cluster_bases,
                                mismatches_by_qscore[qscore].mixed_strand_cluster_mismatches,
                                math::ErrorRateToPhred(mixed_strand_cluster_error_rate));
    } else {
      writer << std::make_tuple(qscore,
                                mismatches_by_qscore[qscore].total_bases,
                                mismatches_by_qscore[qscore].total_mismatches,
                                math::ErrorRateToPhred(error_rate));
    }
  }
}

QscoreStats& QscoreStats::operator+=(const QscoreStats& other) {
  for (size_t i = 0; i < mismatches_by_qscore.size(); ++i) {
    mismatches_by_qscore[i].total_bases += other.mismatches_by_qscore[i].total_bases;
    mismatches_by_qscore[i].total_mismatches += other.mismatches_by_qscore[i].total_mismatches;
    mismatches_by_qscore[i].forward_cluster_bases += other.mismatches_by_qscore[i].forward_cluster_bases;
    mismatches_by_qscore[i].forward_cluster_mismatches += other.mismatches_by_qscore[i].forward_cluster_mismatches;
    mismatches_by_qscore[i].reverse_cluster_bases += other.mismatches_by_qscore[i].reverse_cluster_bases;
    mismatches_by_qscore[i].reverse_cluster_mismatches += other.mismatches_by_qscore[i].reverse_cluster_mismatches;
    mismatches_by_qscore[i].mixed_strand_cluster_bases += other.mismatches_by_qscore[i].mixed_strand_cluster_bases;
    mismatches_by_qscore[i].mixed_strand_cluster_mismatches +=
        other.mismatches_by_qscore[i].mixed_strand_cluster_mismatches;
  }
  return *this;
}

QscoreStats& QscoreStats::operator+=(const vec<MismatchCounts<u64>>& other) {
  for (size_t i = 0; i < mismatches_by_qscore.size(); ++i) {
    mismatches_by_qscore[i].total_bases += other[i].total_bases;
    mismatches_by_qscore[i].total_mismatches += other[i].total_mismatches;
    mismatches_by_qscore[i].forward_cluster_bases += other[i].forward_cluster_bases;
    mismatches_by_qscore[i].forward_cluster_mismatches += other[i].forward_cluster_mismatches;
    mismatches_by_qscore[i].reverse_cluster_bases += other[i].reverse_cluster_bases;
    mismatches_by_qscore[i].reverse_cluster_mismatches += other[i].reverse_cluster_mismatches;
    mismatches_by_qscore[i].mixed_strand_cluster_bases += other[i].mixed_strand_cluster_bases;
    mismatches_by_qscore[i].mixed_strand_cluster_mismatches += other[i].mixed_strand_cluster_mismatches;
  }
  return *this;
}

QscoreStats& QscoreStats::operator+=(const ankerl::unordered_dense::map<u8, MismatchCounts<u32>>& other) {
  for (const auto& [qscore, counts] : other) {
    if (qscore < mismatches_by_qscore.size()) {
      mismatches_by_qscore[qscore].total_bases += counts.total_bases;
      mismatches_by_qscore[qscore].total_mismatches += counts.total_mismatches;
      mismatches_by_qscore[qscore].forward_cluster_bases += counts.forward_cluster_bases;
      mismatches_by_qscore[qscore].forward_cluster_mismatches += counts.forward_cluster_mismatches;
      mismatches_by_qscore[qscore].reverse_cluster_bases += counts.reverse_cluster_bases;
      mismatches_by_qscore[qscore].reverse_cluster_mismatches += counts.reverse_cluster_mismatches;
      mismatches_by_qscore[qscore].mixed_strand_cluster_bases += counts.mixed_strand_cluster_bases;
      mismatches_by_qscore[qscore].mixed_strand_cluster_mismatches += counts.mixed_strand_cluster_mismatches;
    }
  }
  return *this;
}

}  // namespace xoos::alignment_metrics
