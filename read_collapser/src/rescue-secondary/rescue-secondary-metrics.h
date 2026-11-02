#pragma once

#include <optional>

#include <xoos/io/bam/read-group.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

namespace xoos::read_collapser::rescue_secondary {

constexpr std::string_view kDefaultMetricsFilename = "rescue_secondary_metrics.tsv";

struct RescueSecondaryMetrics {
  u64 num_records{};
  u64 num_primary_records{};
  u64 num_secondary_records{};
  u64 num_supplementary_records{};
  u64 num_rescued_secondary_records{};

  // First @RG read group from the input header, used to annotate the metrics TSV with sample info. Captured during
  // processing because the input may be a stream (stdin) that cannot be reopened. std::nullopt when the header has no
  // @RG line.
  std::optional<io::ReadGroup> read_group;

  void WriteToTsv(const fs::path& output, const io::TsvMetadata& tsv_metadata) const;
};

}  // namespace xoos::read_collapser::rescue_secondary
