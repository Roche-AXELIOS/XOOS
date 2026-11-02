#pragma once
#include <chrono>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>

#include <xoos/types/vec.h>

#include "adapter-design/adapter-design-bundle.h"
#include "adapter-design/adapter-design.h"
#include "adapters/duplex/demux-and-trim-duplex.h"
#include "adapters/duplex/hairpin-finder.h"
#include "adapters/duplex/lut-bundle-duplex.h"
#include "adapters/duplex/stem/demux-and-trim-duplex-stem.h"
#include "adapters/duplex/stem/lut-bundle-duplex-stem.h"
#include "adapters/duplex/umi/demux-and-trim-duplex-umi.h"
#include "adapters/duplex/umi/lut-bundle-duplex-umi.h"
#include "adapters/simplex-10x/demux-and-trim-simplex-10x.h"
#include "adapters/simplex-10x/lut-bundle-simplex-10x.h"
#include "adapters/simplex/demux-and-trim-simplex.h"
#include "adapters/ys/demux-and-trim-ys.h"
#include "adapters/ys/lut-bundle-ys.h"
#include "adapters/ysu/demux-and-trim-ysu.h"
#include "adapters/ysu/lut-bundle-ysu.h"
#include "core/demux-and-trim-pipeline.h"
#include "io/sample-sheet/sample-sheet.h"
#include "sequence/matcher/match-info.h"
#include "task/sink.h"

namespace xoos::demux {

/// (output_path, input_name_prefix)
using OutputFileEntry = std::pair<fs::path, std::string>;
using OutputFileMap = std::map<fs::path, vec<OutputFileEntry>>;

class FlowContext;

/**
 * @brief FlowManager initializes and schedules tasks for the demuxing process. It is responsible for:
 * 1. Initializing the flow context for each input file
 * 2. Initialize the sequence writer used to write the demuxed data
 */
class FlowManager {
 public:
  /**
   * @brief References to the per-adapter-type LUT bundle members that LoadSidPool populates.
   *
   * Bundling these into one struct keeps LoadSidPool's parameter list small; exactly one member is
   * assigned per call, selected by the adapter type.
   */
  struct LutBundleRefs {
    std::unique_ptr<LutBundleYsu>& ysu;
    std::unique_ptr<LutBundleYs>& ys;
    // Simplex consumes the adapter components directly (no named LUT bundle); the design bundle owns
    // the component vector for the lifetime of the demux object built from it.
    std::unique_ptr<AdapterDesignBundle>& simplex;
    std::unique_ptr<LutBundleSimplex10x>& simplex_10x;
    std::unique_ptr<LutBundleDuplex>& duplex;
    std::unique_ptr<LutBundleDuplexUmi>& duplex_umi;
    std::unique_ptr<LutBundleDuplexStem>& duplex_stem;
  };

  FlowManager(const DemuxAndTrimParam& param, const AdapterDesign& design);
  ~FlowManager();

  tf::Executor& Executor() const {
    return *_executor;
  }

  FlowManager(const FlowManager&) = delete;
  FlowManager& operator=(const FlowManager&) = delete;

  // For running the demuxing process
  const HairpinFinder& HairpinFinderObj() const;
  const DemuxAndTrimDuplex& DemuxObjectDuplex() const;
  const DemuxAndTrimYsu& DemuxObjectYsu() const;
  const DemuxAndTrimYs& DemuxObjectYs() const;
  const DemuxAndTrimSimplex& DemuxObjectSimplex() const;
  const DemuxAndTrimSimplex10x& DemuxObjectSimplex10x() const;

  AdapterType AdapterTypeValue() const {
    return _adapter_type;
  }

  bool IsDuplex() const {
    return _is_duplex;
  }

  const DemuxAndTrimParam& Param() const {
    return _param;
  }

  const SidPool& SidPoolValue() const {
    return _sid_pool;
  }

  // Init SID pool
  static SidPool LoadSidPool(const DemuxAndTrimParam& param,
                             const AdapterDesign& design,
                             LutBundleRefs lut_bundles,
                             DualSidMatchingPolicy& dual_sid_matching_policy);

  // Register a mapping from output file prefix to input file path for provenance tracking.
  void RegisterOutputPrefix(const std::string& output_prefix, std::string_view input_filepath) const;

  // Register an output file created by a sink worker.
  void RegisterOutputFile(const fs::path& output_path, const std::string& input_name_prefix) const;

 private:
  // Initialize the adapter-specific demux object based on the adapter type.
  void InitDemuxObject(const AdapterDesign& design);

  // Resolve the SID pool for a simplex design, wiring up dual-SID matching when the sample sheet
  // provides paired SIDs. Extracted from LoadSidPool's switch to keep each case small.
  static SidPool LoadSimplexSidPool(const DemuxAndTrimParam& param,
                                    const AdapterDesign& design,
                                    const std::optional<SampleSheetSidPools>& sample_sheet_pools,
                                    const std::optional<BarcodePool>& single_sid_pool,
                                    std::unique_ptr<AdapterDesignBundle>& design_bundle_simplex,
                                    DualSidMatchingPolicy& dual_sid_matching_policy);

  // useful control variable derived from AdapterDesign object
  const AdapterType _adapter_type;
  const bool _is_duplex;
  const DemuxAndTrimParam& _param;

  std::unique_ptr<tf::Executor> _executor;
  vec<std::unique_ptr<FlowContext>> _flow_contexts;
  std::unique_ptr<DemuxAndTrimYsu> _demux_and_trim_ysu;
  std::unique_ptr<DemuxAndTrimYs> _demux_and_trim_ys;
  std::unique_ptr<DemuxAndTrimSimplex> _demux_and_trim_simplex;
  std::unique_ptr<DemuxAndTrimSimplex10x> _demux_and_trim_simplex_10x;
  std::unique_ptr<DemuxAndTrimDuplex> _demux_and_trim_duplex;
  std::unique_ptr<DemuxAndTrimDuplexUmi> _demux_and_trim_duplex_umi;
  std::unique_ptr<DemuxAndTrimDuplexStem> _demux_and_trim_duplex_stem;
  // Hairpin finder for any duplex adapter subtype. Used by the demux step; the old DemuxAndTrimDuplex
  // hierarchy is retained for FindUMIPos/FindStartAdapterInConsensus used by the alignment step.
  std::unique_ptr<HairpinFinder> _hairpin_finder;
  // NOTE: LUT bundles must be declared before _sid_pool since LoadSidPool() initializes them
  std::unique_ptr<LutBundleYsu> _lut_bundle_ysu;
  std::unique_ptr<LutBundleYs> _lut_bundle_ys;
  // Simplex consumes the adapter components directly (no named LUT bundle); the design bundle owns
  // the component vector for the lifetime of the demux object built from it.
  std::unique_ptr<AdapterDesignBundle> _design_bundle_simplex;
  std::unique_ptr<LutBundleSimplex10x> _lut_bundle_simplex_10x;
  std::unique_ptr<LutBundleDuplex> _lut_bundle_duplex;
  std::unique_ptr<LutBundleDuplexUmi> _lut_bundle_duplex_umi;
  std::unique_ptr<LutBundleDuplexStem> _lut_bundle_duplex_stem;
  DualSidMatchingPolicy _dual_sid_matching_policy;

  // Initialized by LoadSidPool(), which uses the LUT bundles above thus needs to be defined after they
  // are initialized.
  const SidPool _sid_pool;

  std::chrono::time_point<std::chrono::high_resolution_clock> _start_time, _end_time;

  // For each SID, we'll create a sink buffer to store the demuxed records.
  vec<std::shared_ptr<Sink>> _sinks;
  size_t _nr_sequences{0};
  size_t _nr_bases{0};
  size_t _nr_bytes{0};
  size_t _total_file_size{0};
  size_t _nr_files{0};

  // Maps output file prefix to input file path, used for provenance tracking during rename.
  // Mutable because FlowContext holds a const reference to FlowManager but needs to register prefixes.
  mutable std::map<std::string, std::string> _output_prefix_to_input;
  // Tracks output files created by sink workers, grouped by sample directory.
  mutable OutputFileMap _output_files;
  mutable std::mutex _provenance_mutex;

  // Rename output FASTQ files to standard naming convention after processing.
  void RenameOutputFiles();

  // Guarantees a single valid-but-empty primary FASTQ for every sample that received zero reads across all inputs,
  // using the global _output_files map to find them. Gated on an explicit --sample-sheet. Must run after all inputs and
  // before RenameOutputFiles() so the file is renamed and recorded in the provenance TSV.
  void MaterializeEmptySamples() const;

  /**
   * @brief Writes a provenance TSV mapping renamed output files back to their source input files.
   *
   * @param[in] sample_dir  Directory in which to write the TSV.
   * @param[in] provenance  List of (new_filename, input_filepath) pairs.
   * @throws error::Error if the TSV already exists and overwrite is false.
   */
  void WriteProvenanceTsv(const fs::path& sample_dir, const vec<std::pair<std::string, std::string>>& provenance) const;
};
}  // namespace xoos::demux
