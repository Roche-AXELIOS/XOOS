#include "io/write-baf-values.h"

#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <tuple>

#include <csv.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/math.h>

#include "io/column-names.h"
#include "io/fai.h"
#include "io/write-bigwig.h"
#include "util/parse-util.h"

namespace xoos::cnc {
namespace {

/// Compute BAF from allele depths, returning NaN when total depth is zero.
f64 ComputeBaf(const f64 ref_ad, const f64 alt_ad) {
  if (math::IsCloseToZero(ref_ad + alt_ad)) {
    return std::numeric_limits<f64>::quiet_NaN();
  }
  return alt_ad / (ref_ad + alt_ad);
}

}  // namespace

void WriteBafValues(const Observations& ref_ads,
                    const Observations& alt_ads,
                    const fs::path& out_fname,
                    const io::CommandLineInfo& command_line_info) {
  std::ofstream ofs(out_fname);
  if (!ofs.is_open()) {
    Logging::Error("Failed to open output file: {}", out_fname.string());
    throw std::runtime_error("Failed to open output file");
  }
  auto writer = csv::make_tsv_writer(ofs);
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info);
  }
  const vec<std::string> column_names{
      kColumnContigAsComment, kColumnStart, kColumnEnd, kColumnBAF, kColumnRefAD, kColumnAltAD};
  writer << column_names;
  if (ref_ads.contigs.empty() || alt_ads.contigs.empty()) {
    Logging::Warn("WriteBafValues: No variants to write!");
    return;
  }
  for (size_t i = 0; i < ref_ads.contigs.size(); ++i) {
    const auto idx = static_cast<Eigen::Index>(i);
    const auto baf = ComputeBaf(ref_ads.obvs[idx], alt_ads.obvs[idx]);
    if (std::isnan(baf)) {
      Logging::Warn("WriteBafValues: zero total depth at {}:{}-{}, writing missing-value token",
                    ref_ads.contigs[i],
                    ref_ads.starts[idx],
                    ref_ads.ends[idx]);
    }
    writer << std::make_tuple(ref_ads.contigs[i],
                              ref_ads.starts[idx],
                              ref_ads.ends[idx],
                              FloatAsString(baf),
                              FloatAsIntString(ref_ads.obvs[idx]),
                              FloatAsIntString(alt_ads.obvs[idx]));
  }
  ofs.close();
}

void WriteBafValuesToBigWig(const Observations& ref_ads,
                            const Observations& alt_ads,
                            const fs::path& out_fname,
                            const fs::path& fai_filename) {
  if (ref_ads.contigs.empty() || alt_ads.contigs.empty()) {
    Logging::Warn("WriteBafValuesToBigWig: No variants to write!");
    std::ofstream ofs(out_fname);
    ofs.close();
    return;
  }
  vec<f64> baf_values;
  baf_values.reserve(static_cast<size_t>(ref_ads.obvs.size()));
  for (Eigen::Index i = 0; i < ref_ads.obvs.size(); ++i) {
    const auto baf = ComputeBaf(ref_ads.obvs[i], alt_ads.obvs[i]);
    if (std::isnan(baf)) {
      Logging::Warn("WriteBafValuesToBigWig: zero total depth at {}:{}-{}, writing NaN",
                    ref_ads.contigs[static_cast<size_t>(i)],
                    ref_ads.starts[i],
                    ref_ads.ends[i]);
    }
    baf_values.push_back(baf);
  }
  // figure out which contigs are where in the provided observations
  StrMap<std::tuple<size_t, size_t>> contig_intervals;
  std::string prev_contig{};
  size_t prev_i = 0;
  for (size_t i = 0; i < ref_ads.contigs.size(); ++i) {
    const auto& contig = ref_ads.contigs[i];
    if (!prev_contig.empty() && contig != prev_contig) {
      contig_intervals[prev_contig] = std::make_tuple(prev_i, i);
      prev_i = i;
    }
    prev_contig = contig;
  }
  // account for the last contig
  contig_intervals[ref_ads.contigs.back()] = std::make_tuple(prev_i, ref_ads.contigs.size());
  // finally, initialize and fill in the return object, in the order of the contigs provided from the FAI file
  Observations baf_observations;
  baf_observations.contigs.resize(ref_ads.contigs.size());
  baf_observations.starts.resize(ref_ads.starts.size());
  baf_observations.ends.resize(ref_ads.ends.size());
  // zero-initialize; entries left unset below default to a BAF of 0.
  baf_observations.obvs = Eigen::VectorXd::Zero(ref_ads.obvs.size());
  // order of contigs should follow from FAI file
  const auto contigs_in_order = GetContigsInOrder(fai_filename);
  size_t i = 0;
  for (const auto& contig : contigs_in_order) {
    if (!contig_intervals.contains(contig)) {
      Logging::Warn("Contig {} found in FAI file but not in BAF observations", contig);
      // this contig is not in the observations
      continue;
    }
    const auto [begin, end] = contig_intervals[contig];
    for (size_t j = begin; j < end; ++j) {
      const auto ei = static_cast<Eigen::Index>(i);
      const auto ej = static_cast<Eigen::Index>(j);
      baf_observations.contigs[i] = ref_ads.contigs[j];
      baf_observations.starts[ei] = ref_ads.starts[ej];
      baf_observations.ends[ei] = ref_ads.ends[ej];
      baf_observations.obvs[ei] = baf_values[j];
      ++i;
    }
  }
  if (i != baf_observations.contigs.size()) {
    Logging::Warn("Number of BAF observations written to BigWig ({}) does not match expected count ({})",
                  i,
                  baf_observations.contigs.size());
    throw std::runtime_error("BAF writing error");
  }
  // Write to BigWig using the existing function
  WriteObservationsToBigWig(baf_observations, out_fname, fai_filename);
}

void WriteBafFiles(const Observations& ref_obvs,
                   const Observations& alt_obvs,
                   const std::optional<fs::path>& baf_out,
                   const std::optional<fs::path>& baf_bw_out,
                   const std::optional<fs::path>& reference_genome_fai_fname,
                   const io::CommandLineInfo& command_line_info) {
  if (baf_out.has_value()) {
    Logging::Info("Writing BAFs BED to {}", baf_out->string());
    WriteBafValues(ref_obvs, alt_obvs, baf_out.value(), command_line_info);
  }
  if (baf_bw_out.has_value() && reference_genome_fai_fname.has_value()) {
    Logging::Info("Writing BAFs BigWig to {}", baf_bw_out->string());
    WriteBafValuesToBigWig(ref_obvs, alt_obvs, baf_bw_out.value(), reference_genome_fai_fname.value());
  } else if (baf_bw_out.has_value() && !reference_genome_fai_fname.has_value()) {
    Logging::Info("Reference FAI file was not provided and is required to write BAF values to BigWig");
  }
}

}  // namespace xoos::cnc
