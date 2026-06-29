#include "task/flow-manager.h"

#include <xoos/log/logging.h>

#include <algorithm>

#include "lut-bundle/lut-bundle.h"
#include "metrics/duplex-metrics.h"
#include "metrics/simplex-metrics.h"
#include "task/flow-context.h"
#include "utility/stop-watch.h"

namespace xoos::demux {

// The following code is taken from the Taskflow library. It is an incredible dirty hack to set the processor
// affinity as the library does not expose the native thread handles. So: I am creating a shadow class to
// get access to those handles by casting to a shadow class. Ugggggghhh, I feel so dirty.

struct ExecutorShadow {
  const size_t kSteals = 0;

  std::mutex wsq_mutex;
  std::mutex taskflows_mutex;

#ifdef __cpp_lib_atomic_wait
  std::atomic<size_t> num_topologies{0};
  std::atomic_flag all_spawned = ATOMIC_FLAG_INIT;
#else
  std::condition_variable topology_cv;
  std::mutex topology_mutex;
  size_t num_topologies{0};
#endif

  std::unordered_map<std::thread::id, size_t> wids;
  std::vector<std::thread> threads;
};

FlowManager::FlowManager(const DemuxAndTrimParam& param, const AdapterDesign& design)
    : adapter_type(design.type),
      is_duplex([this]() {
        switch (adapter_type) {
          using enum AdapterType;
          case kDuplex:
          case kDuplexUMI:
          case kDuplexStem:
            return true;
          default:
            return false;
        }
      }()),
      param(param),
      sid_pool(LoadSidPool(param, design, _lut_bundle_ysu_sl, _lut_bundle_ysu_te, _lut_bundle_ys, _lut_bundle_duplex,
                           _lut_bundle_duplex_umi, _lut_bundle_duplex_stem)) {
  // initialize Metrics configuration
  metrics_constraints::max_sid_id_index = sid_pool.size() - 1;
  metrics_constraints::max_logged_read_length = param.length_distribution_report_max;

  switch (design.type) {
    using enum AdapterType;
    case kYsuSl: {
      auto enable_partial = param.read_length_mode == ReadLengthMode::kAll;
      _demux_and_trim_ysu_sl = std::make_unique<DemuxAndTrimYsuSl>(enable_partial, *_lut_bundle_ysu_sl);
      break;
    }
    case kYsuTe: {
      auto enable_partial = param.read_length_mode == ReadLengthMode::kAll;
      _demux_and_trim_ysu_te = std::make_unique<DemuxAndTrimYsuTe>(enable_partial, *_lut_bundle_ysu_te);
      break;
    }
    case kYs: {
      auto enable_partial = param.read_length_mode == ReadLengthMode::kAll;
      _demux_and_trim_ys = std::make_unique<DemuxAndTrimYs>(enable_partial, *_lut_bundle_ys);
      break;
    }
    case kDuplex:
      _demux_and_trim_duplex = std::make_unique<DemuxAndTrimDuplex>(*_lut_bundle_duplex);
      break;
    case kDuplexUMI:
      _demux_and_trim_duplex_umi = std::make_unique<DemuxAndTrimDuplexUmi>(*_lut_bundle_duplex_umi);
      break;
    case kDuplexStem:
      _demux_and_trim_duplex_stem = std::make_unique<DemuxAndTrimDuplexStem>(*_lut_bundle_duplex_stem);
      break;
    default:
      throw error::Error("Unsupported adapter type: {}", static_cast<s32>(design.type));
  }
  _start_time = std::chrono::high_resolution_clock::now();

  // Some systems have plenty of cores but don't have enough memory to run all datasets at the same time. We therefore
  // will chunk up our input data assuming that every input dataset can keep 4 cores busy. This is a heuristic and
  // might not work for all systems.
  constexpr size_t kCoresPerDataset = 4;
  // Scale to the number of threads specified
  const size_t chunk_size = std::max(1UL, param.threads / kCoresPerDataset);
  Logging::Info("Demultiplexing...");
  _executor = std::make_unique<tf::Executor>(param.threads);

  for (size_t i = 0; i < param.inputs.size(); i += chunk_size) {
    for (size_t j = i; j < param.inputs.size() && j < (i + chunk_size); ++j) {
      // Creating a new flow context will add a graph to the executor.
      _flow_contexts.emplace_back(std::make_unique<FlowContext>(*this, param.inputs[j]));
      _total_file_size += std::filesystem::file_size(param.inputs[j]);
      ++_nr_files;
    }
    // While we're processing the input, TF will keep generating new tasks for every output file. Eventually,
    // those tasks will be done and the TF will shut down. We need to wait for that to happen.
    _executor->wait_for_all();
    for (const auto& flow_context : _flow_contexts) {
      _nr_sequences += flow_context->NumSequences();
      _nr_bases += flow_context->NumBases();
      _nr_bytes += flow_context->NumBytes();
    }
    // If we arrive here, all tasks are done. We can now safely remove the graphs from the executor.
    _flow_contexts.clear();
    if (Task::HasException()) {
      std::rethrow_exception(Task::GetException());
    }
    if (param.min_concord_dp_bases.has_value()) {
      auto min_concord_dp_bases = DuplexMetrics::MinConcordDupBases();

      // Check if we need to early stop
      if (param.min_concord_dp_bases.value() < min_concord_dp_bases) {
        Logging::Info("Early stopping due to minimum concordant duplex bases reached at {} bases",
                      min_concord_dp_bases);
        break;
      }
    }
  }
  _executor.reset();
  _end_time = std::chrono::high_resolution_clock::now();
}

// Allow SID pool to be specified in the constructor so it can exist a const object.
SidPool FlowManager::LoadSidPool(const DemuxAndTrimParam& param, const AdapterDesign& design,
                                 std::unique_ptr<LutBundleYsuSl>& lut_bundle_ysu_sl,
                                 std::unique_ptr<LutBundleYsuTe>& lut_bundle_ysu_te,
                                 std::unique_ptr<LutBundleYs>& lut_bundle_ys,
                                 std::unique_ptr<LutBundleDuplex>& lut_bundle_duplex,
                                 std::unique_ptr<LutBundleDuplexUmi>& lut_bundle_duplex_umi,
                                 std::unique_ptr<LutBundleDuplexStem>& lut_bundle_duplex_stem) {
  const auto load_bundle_sw = StopWatch{};
  const auto sid_pool_list =
      param.sample_sheet ? std::make_optional(detail::LoadSampleSheet(*param.sample_sheet)) : std::nullopt;
  SidPool sid_pool;
  switch (design.type) {
    using enum AdapterType;
    case kYsuSl:
      lut_bundle_ysu_sl = std::make_unique<LutBundleYsuSl>(
          LoadLutBundle<LutBundleYsuSl>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = detail::LoadSidPool(lut_bundle_ysu_sl->Sid5pPool());
      break;
    case kYsuTe:
      lut_bundle_ysu_te = std::make_unique<LutBundleYsuTe>(
          LoadLutBundle<LutBundleYsuTe>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = detail::LoadSidPool(lut_bundle_ysu_te->Sid5pPool());
      break;
    case kYs:
      lut_bundle_ys = std::make_unique<LutBundleYs>(
          LoadLutBundle<LutBundleYs>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = detail::LoadSidPool(lut_bundle_ys->Sid5pPool());
      break;
    case kDuplex:
      lut_bundle_duplex = std::make_unique<LutBundleDuplex>(
          LoadLutBundle<LutBundleDuplex>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = detail::LoadSidPool(lut_bundle_duplex->Sid5pPool());
      break;
    case kDuplexStem:
      lut_bundle_duplex_stem = std::make_unique<LutBundleDuplexStem>(
          LoadLutBundle<LutBundleDuplexStem>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = detail::LoadSidPool(lut_bundle_duplex_stem->Sid5pPool());
      break;
    case kDuplexUMI:
      lut_bundle_duplex_umi = std::make_unique<LutBundleDuplexUmi>(
          LoadLutBundle<LutBundleDuplexUmi>(param.adapter_design_bundle, design, sid_pool_list, param.threads));
      sid_pool = detail::LoadSidPool(lut_bundle_duplex_umi->Sid5pPool());
      break;
    default:
      throw error::Error("Unsupported adapter type: {}", static_cast<s32>(design.type));
  }
  const auto load_bundle_seconds = load_bundle_sw.ElapsedTime<std::chrono::milliseconds>();
  Logging::Info("Loaded adapter design bundle in {} ms", load_bundle_seconds);

  if (sid_pool.empty()) {
    throw error::Error("SID pool is empty after loading adapter design bundle; expected at least one SID.");
  }
  // copy elision should be happening here now
  return sid_pool;
}

FlowManager::~FlowManager() {
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(_end_time - _start_time).count();
  if (is_duplex) {
    // Aggregate metrics for duplex
    const auto global_results = DuplexMetrics::SumTotal();
    global_results.WriteMetrics(param, sid_pool);
    // if we have strand detection turn on log info
    if (param.strand_detector.has_value()) {
      ReportStrandMetrics(global_results);
    }
  } else {
    const auto simplex_metrics_total = SimplexMetrics::SumTotal();
    simplex_metrics_total.WriteMetrics(param, sid_pool);
  }
  // Calculate Gbases per minute.
  const auto gbases = static_cast<double>(_nr_bases) * 1e-9;
  const auto gbases_per_minute = 60000.0 * gbases / static_cast<double>(duration);
  const auto mega_sequences = 1e-6 * static_cast<double>(_nr_sequences);
  const auto mega_sequences_second = 1000.0 * mega_sequences / static_cast<double>(duration);

  Logging::Info(
      "\nThroughput: processed {} sequences, {} bases, {} bytes in {} ms (using {} workers)\nNet throughput: {:.5} "
      "Gbases/minute, {:.5} MSequences/s",
      _nr_sequences, _nr_bases, _nr_bytes, duration, param.threads, gbases_per_minute, mega_sequences_second);
}

const DemuxAndTrimDuplex& FlowManager::DemuxObjectDuplex() const {
  switch (adapter_type) {
    using enum AdapterType;
    case kDuplex:
      return *_demux_and_trim_duplex;
    case kDuplexStem:
      return *_demux_and_trim_duplex_stem;
    case kDuplexUMI:
      return *_demux_and_trim_duplex_umi;
    default:
      throw error::Error("FlowManager::DemuxObjectDuplex() called for non-duplex adapter type. Saw adapter type: {}",
                         static_cast<u32>(adapter_type));
  }
}

const DemuxAndTrimYsuSl& FlowManager::DemuxObjectYsuSl() const {
  if (!_demux_and_trim_ysu_sl) {
    throw error::Error("FlowManager::DemuxObjectYsuSl() called but YSU-SL demux object is not initialized");
  }
  return *_demux_and_trim_ysu_sl;
}

const DemuxAndTrimYsuTe& FlowManager::DemuxObjectYsuTe() const {
  if (!_demux_and_trim_ysu_te) {
    throw error::Error("FlowManager::DemuxObjectYsuTe() called but YSU-TE demux object is not initialized");
  }
  return *_demux_and_trim_ysu_te;
}

const DemuxAndTrimYs& FlowManager::DemuxObjectYs() const {
  if (!_demux_and_trim_ys) {
    throw error::Error("FlowManager::DemuxObjectYs() called but YS demux object is not initialized");
  }
  return *_demux_and_trim_ys;
}

}  // namespace xoos::demux
