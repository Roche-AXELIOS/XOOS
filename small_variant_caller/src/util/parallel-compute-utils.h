#pragma once

#include <optional>
#include <string>
#include <vector>

#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

#include "compute-bam-features/alignment-reader.h"
#include "compute-bam-features/compute-bam-region-features.h"
#include "core/config.h"
#include "core/score-calculator.h"
#include "core/sex.h"
#include "util/region-util.h"

namespace xoos::svc {

/**
 * Synopsis:
 * This file contains utility functions and data structures for parallel feature computation for filtering variants.
 */

/**
 * @brief Manages state required for each thread to process variants in a region, this state is mutable
 * and cannot be accessed by multiple threads concurrently.
 */
struct WorkerContext {
  io::VcfReader filter_variants_reader;
  io::VcfReader vcf_feature_extraction_reader;
  std::optional<io::VcfReader> popaf_reader;
  vec<ScoreCalculator> calculators;
  vec<AlignmentReader> alignment_readers;

  WorkerContext(const std::string& vcf_file,
                const std::optional<std::string>& popaf_file,
                const std::vector<fs::path>& bam_inputs,
                const SVCConfig& model_config,
                AlignmentReaderCache& alignment_reader_cache);
};

/**
 * @brief The state required to filter variants in a region, this state is mutable and cannot be
 * accessed safely by multiple threads concurrently. A window of FlowContext are maintained to allow for
 * an unbalanced workload across threads at the cost of using additional memory.
 */
struct FlowContext {
  vec<io::VcfRecordPtr> out_records;
};

/**
 * @brief The global state required to filter variants, this state is immutable and shared across all threads.
 */
struct GlobalContext {
  fs::path genome{};
  ComputeBamFeaturesParams bam_feat_params{};
  StrMap<std::string> ref_seqs{};
  ChromIntervalsMap bed_regions{};
  ChromIntervalsMap interest_regions{};
  SVCConfig model_config{};
  io::VcfHeaderPtr hdr{};
  StrUnorderedMap<u32> normalize_targets{};
  Sex sex{};
  std::string chr_x_name{};
  std::string chr_y_name{};
  vec<Interval> chr_x_par{};
  vec<Interval> chr_y_par{};
  bool phased{};
  std::optional<vec<BedRegion>> force_calls{};
  std::optional<StrUnorderedSet> hotspots{};
  std::optional<StrUnorderedSet> block_list{};
  float min_allele_freq_threshold{};
  float weighted_counts_threshold{};
  float hotspot_weighted_counts_threshold{};
  float ml_threshold{};
  float somatic_tn_snv_ml_threshold{};
  float somatic_tn_indel_ml_threshold{};
  float hotspot_ml_threshold{};
  float germline_fail_ml_threshold{};
  float min_phased_allele_freq{};
  float max_phased_allele_freq{};
  u32 min_alt_counts{};
  std::optional<fs::path> skip_variants_vcf{};
  u32 tumor_support_threshold{};
};

std::tuple<ChromToVcfFeaturesMap, ChromToVariantInfoMap, RefInfoMap> ComputeBamAndVcfFeaturesForRegion(
    const GlobalContext& global_ctx,
    WorkerContext& worker_ctx,
    const TargetRegion& region,
    const ChromIntervalsMap& bed_regions,
    const ChromIntervalsMap& interest_regions);

vec<TargetRegion> PartitionVcfRegions(const fs::path& vcf_file, size_t threads);
vec<TargetRegion> PartitionRegionsForSomatic(u64 region_size,
                                             const io::HtsFilePtr& bam_file,
                                             const io::SamHdrPtr& header,
                                             const io::HtsIdxPtr& idx,
                                             const StrMap<vec<Interval>>& bed_regions);

}  // namespace xoos::svc
