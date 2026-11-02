#include "util/mapq-util.h"

#include <cmath>
#include <cstddef>
#include <optional>

#include <xoos/error/error.h>
#include <xoos/types/float.h>
#include <xoos/types/vec.h>

#include "coverage/coverage.h"
#include "observations/observations.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/interval-trees.h"
#include "util/parse-util.h"

namespace xoos::cnc {

/**
 * @brief given coverage records with mean mapping qualities, return observations with those mapping qualities
 * @param coverage object with mean_mapping_quality field populated
 * @return Observations object where the observations consist of the mean_mapping_quality field from `coverage`
 */
Observations GetAvgMapqsFromCoverage(const CoverageRecords& coverage) {
  Observations mapqs;
  const auto n = static_cast<Eigen::Index>(coverage.region.size());
  mapqs.contigs.resize(coverage.region.size());
  mapqs.starts.resize(n);
  mapqs.ends.resize(n);
  mapqs.obvs.resize(n);
  for (size_t i = 0; i < coverage.region.size(); ++i) {
    const auto idx = static_cast<Eigen::Index>(i);
    const auto [contig, start, end] = ParseRegionString(coverage.region[i]);
    mapqs.contigs[i] = contig;
    mapqs.starts[idx] = start;
    mapqs.ends[idx] = end;
    mapqs.obvs[idx] = coverage.mean_mapping_quality[idx];
  }
  return mapqs;
}

/**
 * @brief given a set of segments (each segment represents a collection of observations), and a list of MAPQs per
 * observation, return a list of mean MAPQs per segment
 * @param mapqs
 * @param segments
 * @return
 */
Eigen::VectorXd GetAvgMeanMapqPerSegment(const vec<GenomicSegment>& segments, const Observations& mapqs) {
  Eigen::VectorXd avg_mean_mapqs(static_cast<Eigen::Index>(segments.size()));
  const IntervalTrees mapqs_trees(mapqs);
  Eigen::Index i = 0;
  for (const GenomicSegment& seg : segments) {
    if (!seg.num_obs.has_value() || seg.num_obs.value() == 0) {
      throw error::Error(
          "Segment {}:{}-{} must have non-zero num_obs to calculate mean MAPQ", seg.contig, seg.start, seg.end);
    }
    const vec<size_t> idxs = mapqs_trees.LookUp(seg.contig, seg.start, seg.end);
    f64 avg_mean_mapq = 0;
    for (const size_t idx : idxs) {
      avg_mean_mapq += mapqs.obvs[static_cast<Eigen::Index>(idx)];
    }
    avg_mean_mapq = avg_mean_mapq / static_cast<f64>(seg.num_obs.value());
    avg_mean_mapqs[i++] = avg_mean_mapq;
  }
  return avg_mean_mapqs;
}

/**
 * @brief given a set of segments (each segment represents a collection of observations), and a list of MAPQs per
 * observation, assign each segment's avg_mean_mapq field with the average mean mapq of all the targets in that segment
 * @param mapqs, non-const
 * @param segments, const
 * @return
 */
vec<GenomicSegment>& AssignAvgMeanMapqPerSegment(vec<GenomicSegment>& segments, const Observations& mapqs) {
  const IntervalTrees mapqs_trees(mapqs);
  for (GenomicSegment& seg : segments) {
    if (!seg.num_obs.has_value() || seg.num_obs.value() == 0) {
      throw error::Error(
          "Segment {}:{}-{} must have non-zero num_obs to calculate mean MAPQ", seg.contig, seg.start, seg.end);
    }
    const vec<size_t> idxs = mapqs_trees.LookUp(seg.contig, seg.start, seg.end);
    f64 avg_mean_mapq = 0;
    f64 n = 0;
    for (const size_t idx : idxs) {
      const auto eidx = static_cast<Eigen::Index>(idx);
      if (!std::isnan(mapqs.obvs[eidx])) {
        avg_mean_mapq += mapqs.obvs[eidx];
        n += 1;
      }
    }
    if (n > 0) {
      seg.avg_mean_mapq = avg_mean_mapq / n;
    } else {
      seg.avg_mean_mapq = std::nullopt;
    }
  }
  return segments;
}

}  // namespace xoos::cnc
