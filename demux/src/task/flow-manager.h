#pragma once
#include <chrono>
#include <memory>
#include <unordered_map>
#include <vector>

#include "adapter-design/adapter-design.h"
#include "adapters/duplex/demux-and-trim-duplex.h"
#include "adapters/duplex/lut-bundle-duplex.h"
#include "adapters/duplex/stem/demux-and-trim-duplex-stem.h"
#include "adapters/duplex/stem/lut-bundle-duplex-stem.h"
#include "adapters/duplex/umi/demux-and-trim-duplex-umi.h"
#include "adapters/duplex/umi/lut-bundle-duplex-umi.h"
#include "adapters/ys/demux-and-trim-ys.h"
#include "adapters/ys/lut-bundle-ys.h"
#include "adapters/ysu-sl/demux-and-trim-ysu-sl.h"
#include "adapters/ysu-sl/lut-bundle-ysu-sl.h"
#include "adapters/ysu-te/demux-and-trim-ysu-te.h"
#include "adapters/ysu-te/lut-bundle-ysu-te.h"
#include "core/demux-and-trim-pipeline.h"
#include "task/sink.h"

namespace xoos::demux {
class FlowContext;

/// @brief FlowManager initializes and schedules tasks for the demuxing process. It is responsible for:
/// 1. Initializing the flow context for each input file
/// 2. Initialize the sequence writer used to write the demuxed data

using SidPool = std::unordered_map<uint, Barcode>;

class FlowManager {
 public:
  FlowManager(const DemuxAndTrimParam& param, const AdapterDesign& design);
  ~FlowManager();

  tf::Executor& Executor() const { return *_executor; }

  FlowManager(const FlowManager&) = delete;
  FlowManager& operator=(const FlowManager&) = delete;

  // For running the demuxing process
  const DemuxAndTrimDuplex& DemuxObjectDuplex() const;
  const DemuxAndTrimYsuSl& DemuxObjectYsuSl() const;
  const DemuxAndTrimYsuTe& DemuxObjectYsuTe() const;
  const DemuxAndTrimYs& DemuxObjectYs() const;

  // Init SID pool
  static SidPool LoadSidPool(const DemuxAndTrimParam& param, const AdapterDesign& design,
                             std::unique_ptr<LutBundleYsuSl>& lut_bundle_ysu_sl,
                             std::unique_ptr<LutBundleYsuTe>& lut_bundle_ysu_te,
                             std::unique_ptr<LutBundleYs>& lut_bundle_ys,
                             std::unique_ptr<LutBundleDuplex>& lut_bundle_duplex,
                             std::unique_ptr<LutBundleDuplexUmi>& lut_bundle_duplex_umi,
                             std::unique_ptr<LutBundleDuplexStem>& lut_bundle_duplex_stem);

  // useful control variable derived from AdapterDesign object
  const AdapterType adapter_type;
  const bool is_duplex;
  const DemuxAndTrimParam& param;

 private:
  std::unique_ptr<tf::Executor> _executor;
  std::vector<std::unique_ptr<FlowContext>> _flow_contexts;
  std::unique_ptr<DemuxAndTrimYsuSl> _demux_and_trim_ysu_sl;
  std::unique_ptr<DemuxAndTrimYsuTe> _demux_and_trim_ysu_te;
  std::unique_ptr<DemuxAndTrimYs> _demux_and_trim_ys;
  std::unique_ptr<DemuxAndTrimDuplex> _demux_and_trim_duplex;
  std::unique_ptr<DemuxAndTrimDuplexUmi> _demux_and_trim_duplex_umi;
  std::unique_ptr<DemuxAndTrimDuplexStem> _demux_and_trim_duplex_stem;
  // NOTE: LUT bundles must be declared before sid_pool since LoadSidPool() initializes them
  std::unique_ptr<LutBundleYsuSl> _lut_bundle_ysu_sl;
  std::unique_ptr<LutBundleYsuTe> _lut_bundle_ysu_te;
  std::unique_ptr<LutBundleYs> _lut_bundle_ys;
  std::unique_ptr<LutBundleDuplex> _lut_bundle_duplex;
  std::unique_ptr<LutBundleDuplexUmi> _lut_bundle_duplex_umi;
  std::unique_ptr<LutBundleDuplexStem> _lut_bundle_duplex_stem;

 public:
  // Initialized by LoadSidPool(), which uses the LUT bundles above thus needs to be defined they are initalized
  const SidPool sid_pool;

 private:
  std::chrono::time_point<std::chrono::high_resolution_clock> _start_time, _end_time;

  // For each SID, we'll create a sink buffer to store the demuxed records.
  std::vector<std::shared_ptr<Sink>> _sinks;
  size_t _nr_sequences{0};
  size_t _nr_bases{0};
  size_t _nr_bytes{0};
  size_t _total_file_size{0};
  size_t _nr_files{0};
};
}  // namespace xoos::demux
