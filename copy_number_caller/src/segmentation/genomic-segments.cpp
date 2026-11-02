#include "segmentation/genomic-segments.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <locale>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>

#include <fmt/core.h>
#include <fmt/format.h>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/missing-value.h>
#include <xoos/util/parse-float.h>
#include <xoos/util/parse-int.h>
#include <xoos/util/string-functions.h>

#include "likelihood/likelihood-flags.h"
#include "observations/observations.h"
#include "segmentation/genomic-segment-column-names.h"
#include "segmentation/interval-trees.h"
#include "util/parse-util.h"

namespace xoos::cnc {

/**
 * @brief constructor for GenomicSegment
 * @param segment a Segment object to be converted to GenomicSegmnt
 * @param obvs the corresponding observations for the Segmnet object
 */
GenomicSegment::GenomicSegment(const Segment& segment, const Observations& obvs)
    : contig(obvs.contigs[segment.start]),
      start(obvs.starts[static_cast<Eigen::Index>(segment.start)]),
      // segment.end is 0-based exclusive, so we have to substract one to access the last position
      end(obvs.ends[static_cast<Eigen::Index>(segment.end - 1)]) {
  arr_start = segment.start;
  arr_end = segment.end;
}

/**
 * @brief populate the .mean_logr and the .num_obs fields
 * @param obvs corresponding observations
 */
void GenomicSegment::PopulateOptionalFields(const Observations& obvs) {
  const auto seg_start = static_cast<Eigen::Index>(arr_start.value());
  const auto seg_count = static_cast<Eigen::Index>(arr_end.value() - arr_start.value());
  const Eigen::VectorXd obvs_vec = obvs.obvs.segment(seg_start, seg_count);
  mean_logr = eigen::Mean(obvs_vec);
  num_obs = static_cast<size_t>(obvs_vec.size());
}

namespace {
// Parse the Sex SEG column.
//
// The column is written as an enum name ("male"/"female"/"unknown"). The single-character codes "M", "F" and "N"
// written before that migration are still accepted so SEG files from older releases stay readable; drop that fallback
// once those files have aged out. An empty column parses as unknown -- not logged, because this runs once per row and
// a file with a blank Sex column would emit one line per record.
sex_predict::Sex ParseSexField(const std::string& val) {
  using enum sex_predict::Sex;

  if (val.empty()) {
    return kUnknown;
  }
  if (const auto parsed = enum_util::ParseEnumName<sex_predict::Sex>(val)) {
    return *parsed;
  }
  if (val.size() == 1) {
    switch (std::toupper(val[0], std::locale::classic())) {
      case 'M':
        return kMale;
      case 'F':
        return kFemale;
      case 'N':
        return kUnknown;
      default:
        break;
    }
  }
  throw error::Error("invalid sex value '{}': must be one of {}", val, enum_util::FormatEnumNames<sex_predict::Sex>());
}

// Set a core positional field (id, contig, coordinates). Returns true if the field was handled.
bool TrySetCoreField(GenomicSegment& seg, const std::string_view field, const std::string& val) {
  if (field == kGenomicSegColId) {
    seg.id = val;
  } else if (field == kGenomicSegColContig) {
    seg.contig = val;
  } else if (field == kGenomicSegColStart) {
    // the file format specifies 1-based, but we convert it to 0-based closed here for internal consistency
    const auto parsed_start = util::ParseU64(val);
    if (parsed_start == 0) {
      throw error::Error("start position in segment cannot be 0");
    }
    seg.start = parsed_start - 1;
  } else if (field == kGenomicSegColEnd) {
    // the file format specifies 1-based, but we convert it to 0-based open here for internal consistency
    seg.end = util::ParseU64(val);
  } else if (field == kGenomicSegColArrStart) {
    seg.arr_start = util::ParseU64(val);
  } else if (field == kGenomicSegColArrEnd) {
    seg.arr_end = util::ParseU64(val);
  } else {
    return false;
  }
  return true;
}

// Set a coverage/log-ratio metric field. Returns true if the field was handled.
bool TrySetMetricField(GenomicSegment& seg, const std::string_view field, const std::string& val) {
  if (field == kGenomicSegColNumLogr) {
    seg.num_obs = util::ParseU64(val);
  } else if (field == kGenomicSegColMeanLogr) {
    seg.mean_logr = util::ParseF64(val);
  } else if (field == kGenomicSegColNumSnp) {
    seg.num_snps = util::ParseU64(val);
  } else if (field == kGenomicSegColMBaf) {
    seg.mbaf = IsMissingValue(val) ? std::nullopt : std::optional<f64>(util::ParseF64(val));
  } else if (field == kGenomicSegColMeanDh) {
    seg.mean_dh = util::ParseF64(val);
  } else if (field == kGenomicSegColAvgMeanMapq) {
    seg.avg_mean_mapq = IsMissingValue(val) ? std::nullopt : std::optional<f64>(util::ParseF64(val));
  } else {
    return false;
  }
  return true;
}

// Set a copy-number field. Returns true if the field was handled.
bool TrySetCopyNumberField(GenomicSegment& seg, const std::string_view field, const std::string& val) {
  if (field == kGenomicSegColPloidy) {
    seg.ploidy = util::ParseF64(val);
  } else if (field == kGenomicSegColPurity) {
    seg.purity = util::ParseF64(val);
  } else if (field == kGenomicSegColTotalCopyNumber) {
    seg.total_copy_number = util::ParseU64(val);
  } else if (field == kGenomicSegColExpectedTotalCopyNumber) {
    seg.expected_total_copy_number = util::ParseF64(val);
  } else if (field == kGenomicSegColMinorCopyNumber) {
    seg.minor_copy_number = IsMissingValue(val) ? std::nullopt : std::optional<u64>(util::ParseU64(val));
  } else if (field == kGenomicSegColMajorCopyNumber) {
    seg.major_copy_number = IsMissingValue(val) ? std::nullopt : std::optional<u64>(util::ParseU64(val));
  } else {
    return false;
  }
  return true;
}

// Set a likelihood field. Returns true if the field was handled.
bool TrySetLikelihoodField(GenomicSegment& seg, const std::string_view field, const std::string& val) {
  if (field == kGenomicSegColLogrLikelihood) {
    seg.logr_likelihood = IsMissingValue(val) ? std::nullopt : std::optional<f64>(util::ParseF64(val));
  } else if (field == kGenomicSegColBafLikelihood) {
    seg.baf_likelihood = IsMissingValue(val) ? std::nullopt : std::optional<f64>(util::ParseF64(val));
  } else if (field == kGenomicSegColJointLikelihood) {
    seg.joint_likelihood = IsMissingValue(val) ? std::nullopt : std::optional<f64>(util::ParseF64(val));
  } else {
    return false;
  }
  return true;
}

// Set an annotation/flag field (sex, filter flags, allosome/germline booleans). Returns true if the field was handled.
bool TrySetAnnotationField(GenomicSegment& seg, const std::string_view field, const std::string& val) {
  if (field == kGenomicSegColSex) {
    seg.sex = ParseSexField(val);
  } else if (field == kGenomicSegColFilter) {
    seg.flags = StringToFlags(val);
  } else if (field == kGenomicSegColInAllosome) {
    seg.in_allosome = string::StringToBool(val);
  } else if (field == kGenomicSegColInPseudoAutosomalRegion) {
    seg.in_pseudo_autosomal_region = string::StringToBool(val);
  } else if (field == kGenomicSegColIsGermlineCnv) {
    seg.is_germline_cnv = string::StringToBool(val);
  } else {
    return false;
  }
  return true;
}
}  // namespace

void GenomicSegment::SetField(const std::string& field, const std::string& val) {
  if (TrySetCoreField(*this, field, val)) {
    return;
  }
  if (TrySetMetricField(*this, field, val)) {
    return;
  }
  if (TrySetCopyNumberField(*this, field, val)) {
    return;
  }
  if (TrySetLikelihoodField(*this, field, val)) {
    return;
  }
  if (TrySetAnnotationField(*this, field, val)) {
    return;
  }
  Logging::Error("GenomicSegCol {} not supported", field);
  throw std::runtime_error("unsupported column");
}

void GenomicSegment::SetFieldOrThrowError(const std::string& field, const std::string& val) {
  try {
    SetField(field, val);
  } catch (std::invalid_argument& e) {
    Logging::Error("Invalid value '{}' for field '{}': {}", val, field, e.what());
    throw std::runtime_error("invalid value for field");
  } catch (std::out_of_range& e) {
    Logging::Error("Out of range value '{}' for field '{}': {}", val, field, e.what());
    throw std::runtime_error("out of range value for field");
  } catch (std::runtime_error& e) {
    Logging::Error("Runtime error for field '{}': {}", field, e.what());
    throw;
  }
}

namespace {
// Format a core positional field as a string. Returns nullopt if the field is not in this group.
std::optional<std::string> TryGetCoreField(const GenomicSegment& seg, const std::string_view field) {
  if (field == kGenomicSegColId) {
    return fmt::format("{}", seg.id.value());
  }
  if (field == kGenomicSegColContig) {
    return seg.contig;
  }
  if (field == kGenomicSegColStart) {
    return fmt::format("{}", seg.start);
  }
  if (field == kGenomicSegColEnd) {
    return fmt::format("{}", seg.end);
  }
  if (field == kGenomicSegColArrStart) {
    return fmt::format("{}", seg.arr_start.value());
  }
  if (field == kGenomicSegColArrEnd) {
    return fmt::format("{}", seg.arr_end.value());
  }
  return std::nullopt;
}

// Format a coverage/log-ratio metric field as a string. Returns nullopt if the field is not in this group.
std::optional<std::string> TryGetMetricField(const GenomicSegment& seg, const std::string_view field) {
  if (field == kGenomicSegColNumLogr) {
    return fmt::format("{}", seg.num_obs.value());
  }
  if (field == kGenomicSegColMeanLogr) {
    return fmt::format("{:.2f}", seg.mean_logr.value());
  }
  if (field == kGenomicSegColNumSnp) {
    return fmt::format("{}", seg.num_snps.value());
  }
  if (field == kGenomicSegColMBaf) {
    return util::FormatOrNA(seg.mbaf, "{:.2f}");
  }
  if (field == kGenomicSegColMeanDh) {
    return fmt::format("{:.2f}", seg.mean_dh.value());
  }
  if (field == kGenomicSegColAvgMeanMapq) {
    return util::FormatOrNA(seg.avg_mean_mapq, "{:.2f}");
  }
  return std::nullopt;
}

// Format a copy-number field as a string. Returns nullopt if the field is not in this group.
std::optional<std::string> TryGetCopyNumberField(const GenomicSegment& seg, const std::string_view field) {
  if (field == kGenomicSegColPloidy) {
    return fmt::format("{:.2f}", seg.ploidy.value());
  }
  if (field == kGenomicSegColPurity) {
    return fmt::format("{:.2f}", seg.purity.value());
  }
  if (field == kGenomicSegColTotalCopyNumber) {
    return fmt::format("{}", seg.total_copy_number.value());
  }
  if (field == kGenomicSegColExpectedTotalCopyNumber) {
    return fmt::format("{:.2f}", seg.expected_total_copy_number.value());
  }
  const bool copy_number_missing = !seg.minor_copy_number.has_value() && !seg.major_copy_number.has_value();
  if (field == kGenomicSegColMinorCopyNumber) {
    return copy_number_missing ? std::string{util::kNA} : fmt::format("{}", seg.minor_copy_number.value());
  }
  if (field == kGenomicSegColMajorCopyNumber) {
    return copy_number_missing ? std::string{util::kNA} : fmt::format("{}", seg.major_copy_number.value());
  }
  return std::nullopt;
}

// Format a likelihood field as a string. Returns nullopt if the field is not in this group.
std::optional<std::string> TryGetLikelihoodField(const GenomicSegment& seg, const std::string_view field) {
  if (field == kGenomicSegColLogrLikelihood) {
    return util::FormatOrNA(seg.logr_likelihood, "{:.2f}");
  }
  if (field == kGenomicSegColBafLikelihood) {
    return util::FormatOrNA(seg.baf_likelihood, "{:.2f}");
  }
  if (field == kGenomicSegColJointLikelihood) {
    return util::FormatOrNA(seg.joint_likelihood, "{:.2f}");
  }
  return std::nullopt;
}

// Format an annotation/flag field as a string. Returns nullopt if the field is not in this group.
std::optional<std::string> TryGetAnnotationField(const GenomicSegment& seg, const std::string_view field) {
  if (field == kGenomicSegColSex) {
    return enum_util::FormatEnumName(seg.sex.value());
  }
  if (field == kGenomicSegColFilter) {
    return FlagsToString(seg.flags.value());
  }
  if (field == kGenomicSegColInAllosome) {
    return xoos::string::BoolToString(seg.in_allosome.value());
  }
  if (field == kGenomicSegColInPseudoAutosomalRegion) {
    return xoos::string::BoolToString(seg.in_pseudo_autosomal_region.value());
  }
  if (field == kGenomicSegColIsGermlineCnv) {
    return xoos::string::BoolToString(seg.is_germline_cnv.value());
  }
  return std::nullopt;
}
}  // namespace

std::string GenomicSegment::GetFieldAsString(const std::string& field) const {
  for (auto* const getter :
       {TryGetCoreField, TryGetMetricField, TryGetCopyNumberField, TryGetLikelihoodField, TryGetAnnotationField}) {
    if (const std::optional<std::string> result = getter(*this, field); result.has_value()) {
      return result.value();
    }
  }
  Logging::Error("{} is an unsupported SEG column", field);
  throw std::runtime_error("invalid column");
}

/**
 * @brief converts a Segment to a GenomicSegment (i.e. translates it into genomic coordinates). Requires original
 * observations to complete the conversion.
 * @param obvs original observations from which `segments` was produced
 * @param segments segments generated from obvs
 * @return vec<GenomicSegment>, which are `segments` translated to genomic coordinates
 */
vec<GenomicSegment> ConvertSegmentsToGenomicSegments(const Observations& obvs, const vec<Segment>& segments) {
  vec<GenomicSegment> genomic_segments;
  genomic_segments.reserve(segments.size());
  for (const auto& segment : segments) {
    genomic_segments.emplace_back(segment, obvs);
  }
  return genomic_segments;
}

/**
 * @brief same as above function, but appends results to a given std::vector
 * @param obvs original observations from which `segments` was produced
 * @param segments segments generated from obvs
 * @return vec<GenomicSegment>, which are `segments` translated to genomic coordinates
 */
vec<GenomicSegment>& ConvertSegmentsToGenomicSegments(const Observations& obvs,
                                                      const vec<Segment>& segments,
                                                      vec<GenomicSegment>& genomic_segments) {
  for (const auto& segment : segments) {
    genomic_segments.emplace_back(segment, obvs);
  }
  return genomic_segments;
}

/**
 * @brief adjusts the start and ends fields in each GenomicSegment in a list by setting the end to the midpoint of the
 * current GenomicSegment and the next, and setting the start to the midpoint of the current GenomicSegment and the
 * previous
 * @param genomic_segments std::vector of GenomicSegments
 * @return modifies the given GenomicSegments object
 */
vec<GenomicSegment>& AdjustGenomicSegmentBreakpoints(vec<GenomicSegment>& genomic_segments) {
  for (size_t i = 1; i < genomic_segments.size(); ++i) {
    // adjust [i-1].end [i].start
    const size_t prev_end = genomic_segments[i - 1].end;
    const size_t cur_start = genomic_segments[i].start;
    const size_t mid_point = prev_end + ((cur_start - prev_end) / 2);
    genomic_segments[i - 1].end = mid_point;
    genomic_segments[i].start = mid_point;
  }
  return genomic_segments;
}

/**
 * @brief populates the .num_obvs and .mean_logr fields of a list of GenomicSemgents. Also populate num_snps and
 * mean_dh if those fields are present
 * @param genomic_segments
 * @param logrs
 * @return same genomic_segments object, after populating specified fields
 */
vec<GenomicSegment>& PopulateGenomicSegmentOptionalFields(vec<GenomicSegment>& genomic_segments,
                                                          const Observations& obvs) {
  // we need to turn logrs into an interval tree
  // look up each genomic segment inside the interval tree
  // then update the mean and num logrs
  const IntervalTrees trees(obvs);
  for (auto& gseg : genomic_segments) {
    vec<size_t> idxs = trees.LookUp(gseg.contig, gseg.start, gseg.end);
    if (idxs.empty()) {
      Logging::Info("Segment {}:{}-{} not represented in observations! Will not update optional fields...",
                    gseg.contig,
                    gseg.start,
                    gseg.end);
      if (obvs.IsBAFSegObvs()) {
        gseg.num_snps = 0;
        gseg.mean_dh = NAN;
      } else {
        gseg.num_obs = 0;
        gseg.mean_logr = 0;
        gseg.arr_start = *(std::ranges::min_element(idxs));
        gseg.arr_end = *(std::ranges::max_element(idxs)) + 1;
      }
    } else {
      if (!obvs.IsBAFSegObvs()) {
        // Populate logr related fields
        f64 mean = 0;
        for (const size_t i : idxs) {
          mean += obvs.obvs[static_cast<Eigen::Index>(i)];
        }
        mean = mean / static_cast<f64>(idxs.size());
        gseg.num_obs = idxs.size();
        gseg.mean_logr = mean;
      } else {
        // Populate BAF related fields, DP-weighted DH mean
        f64 total_dp = 0;
        f64 af_dp_sum = 0;
        gseg.num_snps = idxs.size();
        for (const size_t i : idxs) {
          total_dp += sqrt(obvs.dps[i]);
          af_dp_sum += sqrt(obvs.dps[i]) * obvs.obvs[static_cast<Eigen::Index>(i)];
        }
        if (total_dp > 0) {
          gseg.mean_dh = af_dp_sum / total_dp;
        } else {
          gseg.mean_dh = NAN;
        }
      }
    }
  }
  return genomic_segments;
}

bool IsSortedGenomic(const vec<GenomicSegment>& v) {
  StrSet seen_contigs;
  for (size_t i = 1; i < v.size(); ++i) {
    if (v[i].contig != v[i - 1].contig) {
      if (seen_contigs.contains(v[i].contig)) {
        Logging::Warn("contigs not ordered");
        return false;
      }
      seen_contigs.insert(v[i - 1].contig);
    } else {
      if (v[i].start <= v[i - 1].start) {
        Logging::Warn("starts not ordered: [{}] {}:{}-{} vs [{}] {}:{}-{}",
                      i - 1,
                      v[i - 1].contig,
                      v[i - 1].start,
                      v[i - 1].end,
                      i,
                      v[i].contig,
                      v[i].start,
                      v[i].end);
        return false;
      }
      if (v[i].end <= v[i - 1].end) {
        Logging::Warn("ends not ordered: [{}] {}:{}-{} vs [{}] {}:{}-{}",
                      i - 1,
                      v[i - 1].contig,
                      v[i - 1].start,
                      v[i - 1].end,
                      i,
                      v[i].contig,
                      v[i].start,
                      v[i].end);
        return false;
      }
    }
  }
  return true;
}

vec<GenomicSegment>& SortSegments(vec<GenomicSegment>& segments, const vec<std::string>& contigs) {
  StrUnorderedMap<size_t> contig_to_order;
  for (size_t i = 0; i < contigs.size(); ++i) {
    contig_to_order[contigs[i]] = i;
  }
  std::ranges::sort(segments, [&contig_to_order](const GenomicSegment& lhs, const GenomicSegment& rhs) {
    return std::tie(contig_to_order[lhs.contig], lhs.start, lhs.end) <
           std::tie(contig_to_order[rhs.contig], rhs.start, rhs.end);
  });
  return segments;
}
}  // namespace xoos::cnc
