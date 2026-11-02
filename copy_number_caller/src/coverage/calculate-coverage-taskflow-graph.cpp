#include "coverage/calculate-coverage-taskflow-graph.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>

#include <htslib/sam.h>

#include <taskflow/algorithm/for_each.hpp>
#include <taskflow/core/flow_builder.hpp>
#include <taskflow/core/graph.hpp>
#include <taskflow/core/task.hpp>

#include <xoos/io/alignment-reader.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "coverage/calculate-coverage.h"
#include "util/parse-util.h"

namespace xoos::cnc {

// Target work per thread in base pairs, tuned for WGS (~1-2 million bp per thread).
constexpr s32 kTargetBasesPerThread = 1'500'000;

tf::Graph CalculateCoverageTaskflowGraph::BuildGraph() {
  tf::Graph graph;
  tf::FlowBuilder builder(graph);
  _res.region = _augmented_intervals.GetRegions();
  _res.on_target = _augmented_intervals.GetOnTargetStatus();
  const size_t total_regions = _res.region.size();
  _res.total_coverage = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(total_regions));
  _res.count = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(total_regions));
  _res.mean_mapping_quality = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(total_regions));
  _res.has_duplicates = false;
  _exclude_flags_int = 0;
  const s32 flags = bam_str2flag(_exclude_flags.c_str());
  constexpr u32 kMaxExcludeFlags = (static_cast<u32>(BAM_FSUPPLEMENTARY) << 1U) - 1U;
  if (flags < 0 || static_cast<u32>(flags) > kMaxExcludeFlags) {
    throw std::runtime_error("Error: Flag value not supported");
  }
  _exclude_flags_int |= static_cast<u32>(flags);
  auto region_tuple = ParseRegionString(_res.region[0]);
  const auto start = static_cast<s32>(std::get<1>(region_tuple));
  const auto end = static_cast<s32>(std::get<2>(region_tuple));
  s32 region_size = end - start;
  // step size is optimized for WGS, which is around 1 - 2 million bp per thread
  _step = kTargetBasesPerThread / region_size;
  _step = std::min(_step, static_cast<s32>(total_regions));
  _step = std::max(_step, 1);
  Logging::Debug("Step size {}, region size {}", _step, region_size);
  Logging::Debug("Exclude flags {} {}", _exclude_flags, _exclude_flags_int);
  tf::Task tf = builder.for_each_index(0, static_cast<s32>(total_regions), _step, [this](const s32 i) {
    const io::AlignmentReader reader = io::OpenAlignmentReader(_bam_file);
    CalculateCoverageRegion(reader, _res, _exclude_flags_int, _ignore_DN, i, _step);
  });
  const tf::Task empty_task = builder.emplace([]() {});
  tf.precede(empty_task);
  return graph;
}
}  // namespace xoos::cnc
