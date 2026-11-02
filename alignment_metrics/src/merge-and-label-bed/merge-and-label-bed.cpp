#include "merge-and-label-bed/merge-and-label-bed.h"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <utility>

#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/io/fast-bed-parser.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>

#include "core/interval.h"
#include "core/region-labels-set.h"
#include "io/bed-reader.h"
#include "merge-and-label-bed/merge-and-label-bed-options.h"
#include "util/labeled-region-util.h"

namespace xoos::alignment_metrics::merge_and_label_bed {

namespace {

RegionsByChromosome ReadLabeledRegions(const vec<BedAndLabelPair>& inputs, RegionLabelsSet& labels) {
  RegionsByChromosome regions;
  for (const auto& [bed_path, label] : inputs) {
    const u16 mask = labels.AddLabel(label);
    io::ForEachBedRegionInFileFast(bed_path, [&regions, mask](const io::BedRegion& region) {
      regions[region.chromosome].emplace_back(region.start, region.end, mask);
    });
  }
  return regions;
}

// Build the comma-separated label list for an interval's label bit vector, in bit-index order.
std::string FormatLabels(const u16 label_bits, const RegionLabelsSet& labels) {
  std::string result;
  for (size_t bit = 0; bit < labels.Size(); ++bit) {
    if ((label_bits & (1U << bit)) != 0) {
      if (!result.empty()) {
        result += ',';
      }
      result += labels.NameAt(bit);
    }
  }
  return result;
}

void WriteLabeledInterval(std::ostream& out,
                          const std::string_view chromosome,
                          const Interval& interval,
                          const RegionLabelsSet& labels) {
  out << chromosome << '\t' << interval.start << '\t' << interval.end;
  const std::string label_list = FormatLabels(interval.labels, labels);
  if (!label_list.empty()) {
    out << '\t' << label_list;
  }
  out << '\n';
}

}  // namespace

void WriteLabeledBed(const RegionsByChromosome& regions,
                     const RegionLabelsSet& labels,
                     const fs::path& output,
                     const io::CommandLineInfo& command_line_info) {
  std::ofstream out(output);
  if (!out) {
    throw error::Error("Failed to open output BED file for writing: {}", output);
  }
  if (!command_line_info.name.empty()) {
    io::WriteTsvMetadata(out, command_line_info);
  }
  // `regions` is ordered by chromosome (std::map) and by position within each chromosome (merge
  // output), so writing in iteration order yields a sorted BED file.
  for (const auto& [chromosome, intervals] : regions) {
    for (const auto& interval : intervals) {
      WriteLabeledInterval(out, chromosome, interval, labels);
    }
  }
}

void MergeAndLabelBed(const MergeAndLabelBedOptions& options) {
  if (!options.overwrite && file::FileExists(options.output)) {
    throw error::Error("Output file already exists (use --overwrite to replace): {}", options.output);
  }
  const fs::path parent = options.output.parent_path();
  if (!parent.empty()) {
    file::CreateWritableDirectory(parent);
  }

  RegionLabelsSet labels;
  RegionsByChromosome regions = ReadLabeledRegions(options.inputs, labels);
  std::ofstream out(options.output);
  if (!out) {
    throw error::Error("Failed to open output BED file for writing: {}", options.output);
  }
  if (!options.command_line_info.name.empty()) {
    io::WriteTsvMetadata(out, options.command_line_info);
  }
  MergeRegionsByLabel(std::move(regions), [&out, &labels](const std::string_view chromosome, const Interval& interval) {
    WriteLabeledInterval(out, chromosome, interval, labels);
  });
  Logging::Info("Wrote merged labeled BED with {} label(s) to {}", labels.Size(), options.output);
}

}  // namespace xoos::alignment_metrics::merge_and_label_bed
