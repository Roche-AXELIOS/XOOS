#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <htslib/hts.h>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

namespace xoos::alignment_metrics {

/**
 * @brief Metrics for evaluating target enrichment sequencing efficiency.
 *
 * Tracks alignment statistics specific to targeted enrichment sequence runs.
 * Measures how effectively the enrichment process worked by counting alignments that fall within the target regions
 * versus those that fall outside (off-target).
 */
struct TargetEnrichmentReadMetrics {
  u64 total_mapped_alignments_in_bam = 0;
  u64 on_target_alignments = 0;
  u64 on_target_reads_passing_filter = 0;

  TargetEnrichmentReadMetrics() = default;
  explicit TargetEnrichmentReadMetrics(const fs::path& bam_path);
  TargetEnrichmentReadMetrics& operator+=(const TargetEnrichmentReadMetrics& obj);
  static std::vector<std::string> GetHeaders();
  void WriteTsv(const fs::path& output_path, const io::Comments& comments) const;
};

/**
 * @brief Get the total number of mapped alignments in a BAM file.
 *
 * @param options The alignment metrics options containing the BAM file path.
 * @return The total number of mapped alignments.
 */
u64 GetTotalMappedAlignmentsInBam(const fs::path& bam_path);

}  // namespace xoos::alignment_metrics
