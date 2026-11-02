#include "rescue-secondary/rescue-secondary-metrics.h"

#include <fstream>
#include <string>

#include <csv.hpp>

#include <xoos/io/metadata-util.h>

#include "metrics/metrics-format-util.h"

namespace xoos::read_collapser::rescue_secondary {

void RescueSecondaryMetrics::WriteToTsv(const fs::path& output, const io::TsvMetadata& tsv_metadata) const {
  auto out = std::ofstream{output};
  auto writer = csv::make_tsv_writer(out).set_auto_flush(false);
  io::WriteTsvMetadata(out, tsv_metadata);

  writer << vec<std::string>{"metric_name", "value", "percentage", "denominator"};
  writer << FormatRow("num_records", num_records, 0, kNA);
  writer << FormatRow("num_primary_records", num_primary_records, num_records, "num_records");
  writer << FormatRow("num_secondary_records", num_secondary_records, num_records, "num_records");
  writer << FormatRow("num_supplementary_records", num_supplementary_records, num_records, "num_records");
  writer << FormatRow("num_rescued_secondary_records", num_rescued_secondary_records, num_records, "num_records");
}

}  // namespace xoos::read_collapser::rescue_secondary
