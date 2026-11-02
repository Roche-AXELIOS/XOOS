#include "segmentation/parent-specific-binary-segmentation-taskflow-graph.h"

#include <cstddef>
#include <stdexcept>

#include <taskflow/algorithm/for_each.hpp>

#include <xoos/log/logging.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/vec.h>

#include "segmentation/segment-pruning.h"
#include "segmentation/segmentation-options.h"
#include "segmentation/undo-segments.h"

namespace xoos::cnc {
/**
 * @brief build the Taskflow Graph for ParentSpecificBinarySegmentationTaskflowGraph. The logic of the PSCBS algorithm
 * is described here
 * @return
 */
tf::Graph ParentSpecificBinarySegmentationTaskflowGraph::BuildGraph() {
  tf::Graph graph;
  tf::FlowBuilder builder(graph);
  tf::Task init_task = builder.emplace([this]() {
    Logging::Info("Beginning Parent-Specific Binary Segmentation");
    if (_logrs.contigs.empty() || _logrs.starts.size() == 0 || _logrs.ends.size() == 0 || _logrs.obvs.size() == 0) {
      throw std::invalid_argument("logr observations are empty");
    }
    // get unique contigs (we assume logrs are sorted)
    for (size_t i = 0; i < _logrs.contigs.size(); ++i) {
      if (i == 0 || (i > 0 && _logrs.contigs[i] != _logrs.contigs[i - 1])) {
        _contigs.push_back(_logrs.contigs[i]);
      }
    }
  });
  // these are Taskflow Graphs, that must be converted to tasks of "taskflow"
  tf::Task seg_by_seg_task = builder.composed_of(_seg_by_seg_graph).name("segment_by_seg");
  tf::Task seg_by_seg_res_task = builder.emplace([this]() {
    Logging::Info("finished Segmenting by Seed Segments");
    _logr_segments = _seg_by_seg_graph.GetResult();
    Logging::Info("Log-Ratio segmentation yielded {} segments", _logr_segments.size());
    if (_options.undo_logr_mode != SegmentUndoMode::kDisabled) {
      const size_t prev_num_segments = _logr_segments.size();
      switch (_options.undo_logr_mode) {
        case SegmentUndoMode::kIterative:
          _logr_segments = IterativeUndoGenomicSegments(_seed_segments,
                                                        _logr_segments,
                                                        _logrs,
                                                        _options.undo_logr_segments_sd_factor,
                                                        _options.increment_undo_logr_segments_sd_factor,
                                                        _options.max_num_segments);
          break;
        case SegmentUndoMode::kNonIterative:
          _logr_segments = NonIterativeUndoGenomicSegments(
              _seed_segments, _logr_segments, _logrs, _options.undo_logr_segments_sd_factor);
          break;
        case SegmentUndoMode::kDisabled:
          break;
      }
      Logging::Info("Undid {} Log-Ratio segments. There are now {} Log-Ratio segments",
                    prev_num_segments - _logr_segments.size(),
                    _logr_segments.size());
    }
    return 0;
  });
  // determine if vcf is provided
  tf::Task baf_cond_task = builder.emplace([this]() { return _dh_vals.obvs.size() != 0 ? 1 : 0; });
  // convert the baf segmentation taskflow to a task of "taskflow"
  tf::Task baf_seg_task = builder.composed_of(_baf_seg_graph).name("baf_seg");
  // define the "stop" tasks
  tf::Task baf_seg_post_task = builder.emplace([this]() {
    _dh_segments = _baf_seg_graph.GetResult();
    Logging::Info("sorting BAF segments");
    SortSegments(_dh_segments, _contigs);
    _dh_segments = BreakLogrSegmentsByDHSegments(_logr_segments, _logrs, _dh_segments, _dh_vals);
    Logging::Info("BAF segmentation yielded {} segments", _dh_segments.size());
    if (!IsSortedGenomic(_dh_segments)) {
      throw std::runtime_error("segments are not sorted!!");
    }
    PopulateGenomicSegmentOptionalFields(_dh_segments, _logrs);
    PopulateGenomicSegmentOptionalFields(_dh_segments, _dh_vals);
    if (!_options.disable_hierarchical_pruning) {
      _res = PruneSegmentsByClusteringPerSeedSegment(
          _seed_segments, _dh_segments, true, _options.pruning_clustering_parameter, _options.baf_distance_weight);
      // Re-populate mean_dh from raw DH observations after pruning, since merging invalidates it
      PopulateGenomicSegmentOptionalFields(_res, _dh_vals);
    } else {
      Logging::Info("Hierarchical clustering disabled");
      _res = _dh_segments;
    }
  });
  // this will happen if there is no VCF provided
  tf::Task logr_seg_post_task = builder.emplace([this]() {
    Logging::Info("No VCF - returning segmentation from log-ratio values");
    PopulateGenomicSegmentOptionalFields(_logr_segments, _logrs);
    if (!_options.disable_hierarchical_pruning) {
      // logR-only path (no VCF): no BAF channel, so the BAF weight has no effect. Pass the default.
      _res = PruneSegmentsByClusteringPerSeedSegment(_seed_segments,
                                                     _logr_segments,
                                                     false,
                                                     _options.pruning_clustering_parameter,
                                                     kSegmentationDefaultBafDistanceWeight);
    } else {
      Logging::Info("Hierarchical clustering disabled");
      _res = _logr_segments;
    }
  });
  // define the graph
  init_task.precede(seg_by_seg_task);
  seg_by_seg_task.precede(seg_by_seg_res_task);
  seg_by_seg_res_task.precede(baf_cond_task);
  baf_cond_task.precede(logr_seg_post_task, baf_seg_task);
  baf_seg_task.precede(baf_seg_post_task);
  return graph;
}

/**
 * @brief modify seed segments depending on sex and mode (germline vs somatic)
 * Somatic (allele-specific):
 * For males:
 *   ChrX: Consider only PAR1/2 regions (only region with SNPs in it) and treat as diploid.
 *   ChrY: Ignore:
 * Remember that chrY PAR1/2 regions are masked. Thus all reads from this region map to chrX PAR1/2
 * For females:
 * ChrX: Treat like a normal autosomal (i.e. diploid)
 * ChrY: Ignore
 *
 * Germline (Total copy number)
 * For males:
 * ChrX:
 *   For PAR1/2 regions, treat them as diploid baseline when calculating total copy number.
 *   For non PAR1/2 regions, treat them as haploid
 * ChrY:
 *   For PAR1/2 regions, should be masked and thus can ignore.
 *   For non PAR1/2 regions, treat them as haploid
 * For females:
 * ChrX: Treat like a normal autosomal (i.e. diploid)
 * ChrY: Ignore
 *
 */
vec<GenomicSegment> ParentSpecificBinarySegmentationTaskflowGraph::LoadSeedSegments(
    const vec<GenomicSegment>& seed_segments, const sex_predict::Sex sex, const SegmentationMode mode) {
  using enum sex_predict::Sex;
  vec<GenomicSegment> ret = seed_segments;
  if (sex == kUnknown) {
    // remove allosomes from seed segments
    size_t num_removed = std::erase_if(ret, [](const GenomicSegment& seg) { return seg.in_allosome.value_or(false); });
    Logging::Info("Removed {} allosome segments from seed segments because sex is unknown", num_removed);
  } else if (sex == kFemale) {
    // remove chrom Y from seed segments
    size_t num_removed = std::erase_if(ret, [](const GenomicSegment& seg) {
      return seg.in_allosome.value_or(false) && sex_predict::IsInChromY(seg.contig);
    });
    Logging::Info("Removed {} chrY segments from seed segments because sex is female", num_removed);
  } else if (sex == kMale) {
    if (mode == SegmentationMode::kGermline) {
      // remove chrom Y PAR segments from seed segments
      size_t num_removed = std::erase_if(ret, [](const GenomicSegment& seg) {
        return (seg.in_allosome.value_or(false) && sex_predict::IsInChromY(seg.contig) &&
                seg.in_pseudo_autosomal_region.value_or(false));
      });
      Logging::Info("Removed {} chrY PAR segments from seed segments because assuming chrY PAR regions are masked",
                    num_removed);
    } else if (mode == SegmentationMode::kSomatic) {
      // remove any segments if they are in chrY, or if they are in chrX and non-PAR
      size_t num_removed = std::erase_if(ret, [](const GenomicSegment& seg) {
        return (seg.in_allosome.value_or(false) &&
                (sex_predict::IsInChromY(seg.contig) || !seg.in_pseudo_autosomal_region.value_or(false)));
      });
      Logging::Info("Removed {} chrY & non-PAR-chrX segments from seed segments because sex is male", num_removed);
      // Note: For germline mode and male, chrY PAR segments are removed; all other chrX and chrY segments are kept.
    }
  } else {
    Logging::Info("Keeping all allosome seed segments");
  }
  return ret;
}

}  // namespace xoos::cnc
