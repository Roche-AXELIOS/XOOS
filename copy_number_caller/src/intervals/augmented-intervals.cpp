#include "intervals/augmented-intervals.h"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <fstream>
#include <istream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include <fmt/core.h>

#include <csv.hpp>

#include <xoos/error/error.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/container-functions.h>
#include <xoos/util/parse-float.h>
#include <xoos/util/parse-int.h>
#include <xoos/util/string-functions.h>

#include "io/bigwig.h"
#include "io/column-names.h"
#include "io/fai.h"
#include "util/eigen-math.h"
#include "util/parse-util.h"

namespace xoos::cnc {

/**
 *  AugmentedIntervals file format described in README.md
 */
AugmentedIntervals::AugmentedIntervals(std::istream& ifs) {
  LoadIntervalsFile(ifs);
}

/**
 * @brief loads an intervals BED file.
 * Expected BED format: Contig\tStart\tEnd\tGCBias\tMappability\tOnTarget (0-based half-open).
 * `#`-prefixed metadata/comments and a `#`-prefixed column header are supported and skipped.
 * @param ifs input stream
 */
void AugmentedIntervals::LoadIntervalsFile(std::istream& ifs) {
  vec<f64> gc_biases;
  vec<f64> mappabilities;
  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::string contig;
    std::string start_str;
    std::string end_str;
    std::string gc_str;
    std::string mapp_str;
    std::string on_target_str;
    std::istringstream ss(line);
    if (!std::getline(ss, contig, '\t') || !std::getline(ss, start_str, '\t') || !std::getline(ss, end_str, '\t') ||
        !std::getline(ss, gc_str, '\t') || !std::getline(ss, mapp_str, '\t') ||
        !std::getline(ss, on_target_str, '\t')) {
      throw error::Error(
          "Intervals BED file must have 6 columns (Contig, Start, End, GCBias, Mappability, OnTarget): \"{}\"", line);
    }
    // reconstruct 1-based region string for internal use
    try {
      const u64 start = util::ParseU64(start_str);
      const u64 end = util::ParseU64(end_str);
      const std::string r = fmt::format("{}:{}-{}", contig, start + 1, end);
      _region.push_back(r);
      gc_biases.push_back(util::ParseF64(gc_str));
      mappabilities.push_back(util::ParseF64(mapp_str));
      _on_target.push_back(xoos::string::StringToBool(on_target_str));
    } catch (const std::exception& e) {
      throw error::Error("Intervals BED file parse error: {} (line: \"{}\")", e.what(), line);
    }
  }
  _gc_bias = eigen::FromStdVec(gc_biases);
  _mappability = eigen::FromStdVec(mappabilities);
}

/**
 * @brief get choice of on-target OR off-target intervals
 * @param o if true, return on-targets, if false return off-targets
 * @return an AugmentedIntervals object containing the desired target intervals
 */
vec<std::string> AugmentedIntervals::GetOnOrOffTargetRegions(bool o) {
  vec<std::string> ret;
  for (size_t i = 0; i < _region.size(); ++i) {
    if ((o && _on_target[i]) || (!o && !_on_target[i])) {
      ret.push_back(_region[i]);
    }
  }
  return ret;
}

/**
 * @brief set all targets to have on-target status
 */
void AugmentedIntervals::SetAllOnTarget() {
  _on_target = vec<bool>(_region.size(), true);
}

/**
 * @brief sort the object by contig and position, contig order defined by a provided FAI file
 * @param fai_file fai file containing desired order of contigs
 */
void AugmentedIntervals::SortByRegion(const fs::path& fai_file) {
  if (_region.empty()) {
    Logging::Error("cannot SortByRegion because there are no regions for this AugmentedIntervals object!");
    throw std::runtime_error("empty region std::vector");
  }
  StrUnorderedMap<size_t> contig_order = GetContigOrder(fai_file);
  // this will contain the "correct" order by which to reorder everything else
  vec<u64> argsort(_region.size());
  // Seed with sequential indices [0, size) so regions can be sorted indirectly.
  std::ranges::iota(argsort, 0);
  std::ranges::stable_sort(argsort, [this, &contig_order](size_t i, size_t j) {
    const auto& [lcontig, lstart, lend] = ParseRegionString(_region[i]);
    const auto& [rcontig, rstart, rend] = ParseRegionString(_region[j]);
    return std::tie(contig_order[lcontig], lstart) < std::tie(contig_order[rcontig], rstart);
  });
  // sort all the member std::vectors
  _region = util::container::Permute(_region, argsort);
  if (!_on_target.empty()) {
    _on_target = util::container::Permute(_on_target, argsort);
  }
  if (_gc_bias.size() != 0) {
    _gc_bias = eigen::Elem(_gc_bias, argsort);
  }
  if (_mappability.size() != 0) {
    _mappability = eigen::Elem(_mappability, argsort);
  }
}

/**
 * @brief annotate intervals with their respective GC contents
 * @param fa a FastaReader object
 */
void AugmentedIntervals::AnnotateGC(io::FastaReader& fa) {
  _gc_bias.resize(static_cast<Eigen::Index>(_region.size()));
  for (size_t i = 0; i < _region.size(); ++i) {
    const auto& [contig, start, end] = ParseRegionString(_region[i]);
    std::string seq(fa.GetSequence(contig, static_cast<s32>(start), static_cast<s32>(end), true));
    const size_t gc_count = std::ranges::count_if(seq.begin(), seq.end(), [](char c) { return c == 'G' || c == 'C'; });
    _gc_bias[static_cast<Eigen::Index>(i)] = static_cast<f64>(gc_count) / static_cast<f64>(seq.size());
  }
}

void AugmentedIntervals::FilterExtremeGC(f64 lower_threshold, f64 upper_threshold) {
  if (_gc_bias.size() == 0) {
    throw std::runtime_error("GC content not annotated. Call AnnotateGC() first.");
  }

  const size_t original_size = _region.size();

  // Filter by creating new std::vectors with only kept elements
  vec<std::string> filtered_region;
  vec<bool> filtered_on_target;
  vec<f64> filtered_gc;
  vec<f64> filtered_mapp;

  for (size_t i = 0; i < _region.size(); ++i) {
    const auto eigen_i = static_cast<Eigen::Index>(i);
    if (_gc_bias[eigen_i] >= lower_threshold && _gc_bias[eigen_i] <= upper_threshold) {
      filtered_region.push_back(_region[i]);
      filtered_on_target.push_back(_on_target[i]);
      filtered_gc.push_back(_gc_bias[eigen_i]);
      if (_mappability.size() != 0) {
        filtered_mapp.push_back(_mappability[eigen_i]);
      }
    }
  }

  // Replace with filtered std::vectors
  _region = std::move(filtered_region);
  _on_target = std::move(filtered_on_target);
  _gc_bias = eigen::FromStdVec(filtered_gc);
  if (_mappability.size() != 0) {
    _mappability = eigen::FromStdVec(filtered_mapp);
  }

  size_t removed = original_size - _region.size();
  Logging::Info("Removed {} intervals with extreme GC content (kept GC range {:.2f}-{:.2f})",
                removed,
                lower_threshold,
                upper_threshold);
}

/**
 * @brief annotate intervals with their respective mappabilities
 * @param mapp_bigwig
 */
void AugmentedIntervals::AnnotateMappability(BigWig& mapp_bigwig) {
  _mappability.resize(static_cast<Eigen::Index>(_region.size()));
  for (size_t i = 0; i < _region.size(); ++i) {
    const auto& [contig, start, end] = ParseRegionString(_region[i]);
    _mappability[static_cast<Eigen::Index>(i)] = mapp_bigwig.GetMean(contig, start, end);
  }
}

/**
 * @brief write AugmentedIntervals to BED file with metadata and column header.
 * @param ofs output file stream
 * @param command_line_info command line metadata including program version and rendered command line
 */
void AugmentedIntervals::Write(std::ofstream& ofs, const io::CommandLineInfo& command_line_info) const {
  auto writer = csv::make_tsv_writer(ofs);
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info);
  }
  const vec<std::string> column_names{
      kColumnContigAsComment, kColumnStart, kColumnEnd, kColumnGCBias, kColumnMappability, kColumnOnTarget};
  writer << column_names;
  for (size_t i = 0; i < _region.size(); ++i) {
    const auto eigen_i = static_cast<Eigen::Index>(i);
    const auto& [contig, start, end] = ParseRegionString(_region[i]);
    // Fixed six-decimal-place output matches the log-ratio and B-allele-frequency tracks.
    writer << std::make_tuple(contig,
                              start,
                              end,
                              FloatAsString(_gc_bias[eigen_i]),
                              FloatAsString(_mappability[eigen_i]),
                              xoos::string::BoolToString(_on_target[i]));
  }
}

}  // namespace xoos::cnc
