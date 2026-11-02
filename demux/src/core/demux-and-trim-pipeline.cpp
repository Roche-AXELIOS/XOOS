#include "demux-and-trim-pipeline.h"

#include <xoos/error/error.h>

#include "adapter-design/adapter-design.h"
#include "metrics/metric-filenames.h"
#include "task/flow-manager.h"

namespace xoos::demux {

namespace fs = std::filesystem;

/**
 * Check metrics validation and delete any existing metrics files if overwrite is enabled. If overwrite is not enabled
 * and any metrics files already exist, an error will be thrown.
 * Because metrics are written in the FlowManager destructor, we need to check for existing metrics files and handle
 * them before the FlowManager is created. This function is called in the DemuxAndTrimPipeline function before the
 * FlowManager is initialized.
 * @param param
 */
void HandleMetricsFiles(const DemuxAndTrimParam& param) {
  const auto metrics_dir = param.out_dir / kMetricsDirectory;
  if (!fs::exists(metrics_dir)) {
    fs::create_directories(metrics_dir);
    return;
  }
  if (!fs::is_directory(metrics_dir)) {
    throw error::Error(fmt::format(
        "Metrics path {} already exists but is not a directory. Please delete or rename this path, or choose a "
        "different output directory.",
        metrics_dir.string()));
  }
  const auto read_length_distr_dir = metrics_dir / kReadLengthDistributionDirectory;
  // Remove all generated distributions on overwrite so stale run- or sample-level files cannot survive.
  if (param.overwrite) {
    fs::remove_all(read_length_distr_dir);
    return;
  }
  if (fs::exists(read_length_distr_dir) &&
      (!fs::is_directory(read_length_distr_dir) || !fs::is_empty(read_length_distr_dir))) {
    throw error::Error(fmt::format(
        "Read-length metrics directory {} already contains output. Please delete it or specify --overwrite.",
        read_length_distr_dir.string()));
  }
  const auto metrics_files = std::vector{// duplex metrics file names used in demux-and-trim pipeline
                                         metrics_dir / kRunMetricsFile,
                                         metrics_dir / kSampleMetricsFile,
                                         metrics_dir / kSampleAssignmentMetrics,

                                         read_length_distr_dir / kTotalReadLengthDistr,
                                         read_length_distr_dir / kRawFailedReadLengthDistr,
                                         read_length_distr_dir / kTotalUntrimmedSimplexPassingReadLengthDistr,
                                         read_length_distr_dir / kTotalTrimmedSimplexPassingReadLengthDistr,
                                         read_length_distr_dir / kTotalUntrimmedDuplexPassingReadLengthDistr,
                                         read_length_distr_dir / kTotalTrimmedDuplexPassingReadLengthDistr,
                                         read_length_distr_dir / kTotalEndadapterPositionDistr,
                                         read_length_distr_dir / kNoHairpinReadLengthDistr};
  for (const auto& file : metrics_files) {
    if (fs::exists(file)) {
      throw error::Error(
          fmt::format("Metrics file {} already exists. Please delete it or specify the --overwrite flag to allow "
                      "overwriting existing metrics files.",
                      file.string()));
    }
  }
}

namespace {

/// @brief Resolve the requested adapter design from the DSL adapter config inside the bundle.
AdapterDesign ResolveAdapterDesign(const DemuxAndTrimParam& param) {
  return LoadAdapterDesign(param.adapter_design_bundle, param.adapter_design_name);
}

}  // namespace

void DemuxAndTrimPipeline(const DemuxAndTrimParam& param) {
  const AdapterDesign design = ResolveAdapterDesign(param);
  fs::create_directories(param.out_dir);
  HandleMetricsFiles(param);
  const FlowManager flow_manager(param, design);
}

}  // namespace xoos::demux
