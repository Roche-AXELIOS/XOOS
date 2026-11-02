#include "merge-segments.h"

#include <cstddef>
#include <optional>
#include <stdexcept>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/vec.h>
#include <xoos/util/math.h>

#include "merge-segments/merge-segments.h"
#include "segmentation/genomic-segment-column-names.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {

/**
 * @brief consequence: deletes all the internal fields for the segment. Does not support SNPs yet
 * @param segments
 * @param first_seg_idx
 * @param last_seg_idx
 * @return
 */
GenomicSegment MergeSegments(const vec<GenomicSegment>& segments_to_merge) {
  GenomicSegment ret;
  const GenomicSegment& first_seg = segments_to_merge[0];
  const GenomicSegment& last_seg = segments_to_merge[segments_to_merge.size() - 1];
  ret.start = first_seg.start;
  ret.end = last_seg.end;
  ret.contig = first_seg.contig;
  ret.arr_start = first_seg.arr_start;
  ret.arr_end = last_seg.arr_end;
  ret.id = first_seg.id;
  ret.sex = first_seg.sex;
  // calculate num_obs, num_segs, the weighted mean of the logrs, snps and mapqs fields
  f64 logr_total_weight = 0;
  f64 mean_logr = 0;
  f64 avg_mean_mapq = 0;
  f64 mapq_total_weight = 0;
  f64 new_logr_likelihood = 0;
  size_t total_num_obs = 0;
  bool in_allosome = segments_to_merge[0].in_allosome.value();
  bool in_par = segments_to_merge[0].in_pseudo_autosomal_region.value();
  for (const auto& seg : segments_to_merge) {
    if (seg.baf_likelihood.has_value() || seg.minor_copy_number.has_value() || seg.major_copy_number.has_value()) {
      throw error::Error("MAPQ merging not supported for segments with major/minor copy number estimations");
    }
    mean_logr += static_cast<f64>(seg.num_obs.value()) * seg.mean_logr.value();
    if (seg.avg_mean_mapq.has_value()) {
      avg_mean_mapq += static_cast<f64>(seg.num_obs.value()) * seg.avg_mean_mapq.value();
      mapq_total_weight += static_cast<f64>(seg.num_obs.value());
    }
    total_num_obs += seg.num_obs.value();
    if (seg.logr_likelihood.has_value()) {
      new_logr_likelihood += seg.logr_likelihood.value();
    }
    logr_total_weight += static_cast<f64>(seg.num_obs.value());
  }
  ret.mean_logr = mean_logr / logr_total_weight;
  ret.avg_mean_mapq = mapq_total_weight > 0 ? std::optional<f64>(avg_mean_mapq / mapq_total_weight) : std::nullopt;
  ret.num_obs = total_num_obs;
  ret.ploidy = first_seg.ploidy;
  ret.purity = first_seg.purity;
  ret.total_copy_number = first_seg.total_copy_number;
  ret.minor_copy_number = first_seg.minor_copy_number;
  // if there were no logr_likelihoods, then we want to keep the return value at nullopt
  if (!math::IsCloseToZero(new_logr_likelihood)) {
    ret.logr_likelihood = new_logr_likelihood / logr_total_weight;
    ret.joint_likelihood = ret.logr_likelihood;
  }
  ret.in_allosome = in_allosome;
  ret.in_pseudo_autosomal_region = in_par;
  return ret;
}

/**
 * @brief Merge segment into both of its flanking segments (A -> B -> C ----> ABC) if its mean avg MAPQ is less than a
 * specified cutoff AND if it is the same copy number as its neighbors AND the segment is not a homozygous deletion
 * (TCN==0)
 * @param mapqs
 * @param logrs
 * @param segments
 * @param cutoff
 */
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
vec<GenomicSegment> MergeLowMapqSegments(const vec<GenomicSegment>& segments, f64 min_mapq_cutoff) {
  if (min_mapq_cutoff <= 0) {
    Logging::Info("mapq threshold for merging either not specified or <= 0. No merging will be done");
    return segments;
  }
  Logging::Info("Finding low-mapq segments in a set of {} segments", segments.size());
  vec<GenomicSegment> new_segments;
  // check the first segment
  for (size_t i = 0; i < segments.size(); ++i) {
    const GenomicSegment& seg = segments[i];
    // skipped segments or segments with no MAPQ observations
    if (!seg.avg_mean_mapq.has_value()) {
      new_segments.emplace_back(seg);
      continue;
    }
    if (!seg.total_copy_number.has_value()) {
      Logging::Error("segment is missing the {} field!", kGenomicSegColTotalCopyNumber);
      throw std::runtime_error("missing column");
    }
    // cannot merge if the seg is a HOM del
    if (seg.total_copy_number.value() != 0 && seg.avg_mean_mapq.value() < min_mapq_cutoff) {
      vec<GenomicSegment> to_merge;
      bool merge_prev_seg = false;
      if (!new_segments.empty()) {
        const GenomicSegment& prev_seg = new_segments.back();
        // flanking segment must be from the same contig, same PAR status, and have MAPQ
        if (prev_seg.contig == seg.contig &&
            prev_seg.in_pseudo_autosomal_region.value() == seg.in_pseudo_autosomal_region.value() &&
            prev_seg.avg_mean_mapq.has_value()) {
          to_merge.emplace_back(prev_seg);
          merge_prev_seg = true;
        }
      }
      to_merge.emplace_back(seg);
      bool merge_next_seg = false;
      if (i + 1 < segments.size()) {
        const GenomicSegment& next_seg = segments[i + 1];
        // flanking segment must be from the same contig, same PAR status, and have MAPQ
        if (next_seg.contig == seg.contig &&
            next_seg.in_pseudo_autosomal_region.value() == seg.in_pseudo_autosomal_region.value() &&
            next_seg.avg_mean_mapq.has_value()) {
          to_merge.emplace_back(next_seg);
          merge_next_seg = true;
        }
      }
      if (!merge_prev_seg || !merge_next_seg) {
        Logging::Debug("Not merging {}:{}-{} because it does not have a flanking segments on one side",
                       seg.contig,
                       seg.start,
                       seg.end);
        new_segments.emplace_back(seg);
      } else if (to_merge.front().total_copy_number.value() != to_merge.back().total_copy_number.value()) {
        Logging::Debug("Not merging {}:{}-{} because its flanking segments have different copy numbers",
                       seg.contig,
                       seg.start,
                       seg.end);
        new_segments.emplace_back(seg);
      } else {
        const GenomicSegment new_seg = MergeSegments(to_merge);
        if (merge_prev_seg) {
          new_segments[new_segments.size() - 1] = new_seg;
        } else {
          new_segments.emplace_back(new_seg);
        }
        if (merge_next_seg) {
          // skip if the next segment has already been merged
          i += 1;
        }
      }
    } else {
      new_segments.emplace_back(seg);
    }
  }
  Logging::Info("after merging in low mapq segments, there are {} segments left", new_segments.size());
  return new_segments;
}

/**
 * @brief Merge adjacent segments if they have equal total copy number and are adjacent to each other
 * @param segments std::vector of GenomicSegment to objects to merge
 * @return std::vector of merged GenomicSegment objects
 */
vec<GenomicSegment> MergeAdjacentEqualCopyNumberSegments(const vec<GenomicSegment>& segments) {
  Logging::Info("Finding adjacent segments with equal copy numbers in a set of {} segments", segments.size());
  vec<GenomicSegment> new_segments;
  new_segments.reserve(segments.size());
  new_segments.push_back(segments[0]);
  for (size_t i = 1; i < segments.size(); ++i) {
    const auto& seg = segments[i];
    auto& last_seg = new_segments.back();
    // NOTE this does not take into account the seed segments from which these segments originated
    if (seg.contig == last_seg.contig && seg.total_copy_number == last_seg.total_copy_number &&
        seg.in_pseudo_autosomal_region == last_seg.in_pseudo_autosomal_region) {
      const vec<GenomicSegment> to_merge{last_seg, seg};
      const GenomicSegment new_seg = MergeSegments(to_merge);
      new_segments.pop_back();
      new_segments.push_back(new_seg);
    } else {
      new_segments.push_back(seg);
    }
  }
  Logging::Info("after merging in adjacent segments with equal copy number, there are {} segments left",
                new_segments.size());
  return new_segments;
}

}  // namespace xoos::cnc
