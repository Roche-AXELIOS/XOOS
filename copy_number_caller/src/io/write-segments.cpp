#include "io/write-segments.h"

#include <string>

#include <fmt/core.h>

#include <csv.hpp>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>

#include "segmentation/genomic-segment-column-names.h"

namespace xoos::cnc {

namespace {
vec<std::string> ChooseColumns(const SegmentType& mode) {
  using enum SegmentType;
  switch (mode) {
    case kLogROnly:
      return kGenomicSegColsRequiredLogRSegments;
    case kBaf:
      return kGenomicSegColsRequiredBafSegments;
    case kGermlineLikelihood:
      return kGenomicSegColsRequiredGermline;
    case kSomaticWithBafLikelihood:
      return kGenomicSegColsRequiredSomaticAlleleSpecific;
    case kSomaticNoBafLikelihood:
      return kGenomicSegColsRequiredSomaticNoAlleleSpecific;
    case kSomaticTumorOnlyWithBafLikelihood:
      return kGenomicSegColsRequiredSomaticTumorOnlyAlleleSpecific;
    default:
      Logging::Error("invalid SegmentType {}", enum_util::FormatEnumName(mode));
      throw std::runtime_error("invalid SegmentType");
  }
}
}  // namespace

void WriteSegments(const fs::path& out_fname,
                   const vec<GenomicSegment>& segments,
                   const SegmentType& mode,
                   const std::optional<io::CommandLineInfo>& command_line_info) {
  std::ofstream ofs(out_fname);
  if (!ofs.is_open()) {
    throw error::Error("Failed to open segments output file: {}", out_fname.string());
  }
  auto writer = csv::make_tsv_writer(ofs);
  if (command_line_info.has_value() && !command_line_info->version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info.value());
  }
  // check if segment has allele-specific information
  const vec<std::string> colnames = ChooseColumns(mode);
  writer << colnames;
  vec<std::string> row_to_write{};
  row_to_write.reserve(colnames.size());
  for (const auto& seg : segments) {
    row_to_write.clear();
    for (const std::string& col : colnames) {
      // convert genomic coordinates back to 1-based
      if (col == kGenomicSegColStart) {
        row_to_write.emplace_back(fmt::format("{}", seg.start + 1));
      } else {
        row_to_write.emplace_back(seg.GetFieldAsString(col));
      }
    }
    writer << row_to_write;
  }
}

void WriteLogRSegments(const vec<GenomicSegment>& segments,
                       const fs::path& logr_segments_out,
                       const Observations& logrs,
                       std::string_view sample_id,
                       const sex_predict::Sex sex,
                       const io::CommandLineInfo& command_line_info) {
  auto segments_to_write = segments;
  Logging::Info("Writing segments to {}", logr_segments_out.string());
  for (auto& seg : segments_to_write) {
    seg.id = sample_id;
    seg.sex = sex;
  }
  PopulateGenomicSegmentOptionalFields(segments_to_write, logrs);
  WriteSegments(logr_segments_out, segments_to_write, SegmentType::kLogROnly, command_line_info);
}

void WriteBafSegments(const vec<GenomicSegment>& segments,
                      const fs::path& baf_segments_out,
                      const io::CommandLineInfo& command_line_info) {
  Logging::Info("Writing segments to {}", baf_segments_out.string());
  WriteSegments(baf_segments_out, segments, SegmentType::kBaf, command_line_info);
}
}  // namespace xoos::cnc
