#pragma once

#include <optional>
#include <string>
#include <utility>

#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/vec.h>

#include "observations/observations.h"

namespace xoos::cnc {

// Copy number and ploidy of a normal diploid background (two homologous chromosomes).
constexpr size_t kDiploidCopyNumber = 2;
constexpr f64 kDiploidPloidy = 2.0;

struct Segment {
  // 0-based inclusive
  size_t start = 0;
  // 0-based exclusive
  size_t end = 0;
  // t-value for this segment compared to parent
  f64 t = 0;
  // p-value for this segment
  f64 p = 0;
  // number of values within this segment
  f64 n = 0;
  // mean value in this segment.
  f64 mean = 0;
};

struct GenomicSegment {
  GenomicSegment() = default;

  GenomicSegment(std::string c, const size_t s, const size_t e) : contig(std::move(c)), start(s), end(e) {
  }

  GenomicSegment(const Segment& segment, const Observations& obvs);

  void PopulateOptionalFields(const Observations& obvs);
  void SetField(const std::string& field, const std::string& val);
  void SetFieldOrThrowError(const std::string& field, const std::string& val);

  std::string GetFieldAsString(const std::string& field) const;
  std::string contig;
  // 0-based inclusive
  size_t start{};
  // 0-based exclusive
  size_t end{};
  std::optional<std::string> id;
  std::optional<size_t> arr_start;
  std::optional<size_t> arr_end;
  std::optional<size_t> num_obs;
  std::optional<f64> mean_logr;
  std::optional<size_t> num_snps;
  std::optional<f64> mean_dh;
  std::optional<f64> mbaf;
  std::optional<f64> avg_mean_mapq;
  std::optional<f64> ploidy;
  std::optional<f64> purity;
  std::optional<size_t> total_copy_number;
  std::optional<f64> expected_total_copy_number;
  std::optional<size_t> minor_copy_number;
  std::optional<size_t> major_copy_number;
  std::optional<f64> logr_likelihood;
  std::optional<f64> baf_likelihood;
  std::optional<f64> joint_likelihood;
  std::optional<sex_predict::Sex> sex;
  std::optional<vec<std::string>> flags;
  std::optional<bool> in_allosome;
  std::optional<bool> in_pseudo_autosomal_region;
  // Set only in tumor-only mode: true if the segment fits the germline log-ratio likelihood model better than the
  // somatic model. Written to the IsGermlineCnv SEG column and the VCF INFO/GERMLINE flag.
  std::optional<bool> is_germline_cnv;
};

vec<GenomicSegment> ConvertSegmentsToGenomicSegments(const Observations& obvs, const vec<Segment>& segments);
vec<GenomicSegment>& ConvertSegmentsToGenomicSegments(const Observations& obvs,
                                                      const vec<Segment>& segments,
                                                      vec<GenomicSegment>& genomic_segment);
vec<GenomicSegment>& PopulateGenomicSegmentOptionalFields(vec<GenomicSegment>& genomic_segments,
                                                          const Observations& obvs);
vec<GenomicSegment>& AdjustGenomicSegmentBreakpoints(vec<GenomicSegment>& genomic_segments);
bool IsSortedGenomic(const vec<GenomicSegment>& v);
vec<GenomicSegment>& SortSegments(vec<GenomicSegment>& segments, const vec<std::string>& contigs);
}  // namespace xoos::cnc
