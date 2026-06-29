#pragma once

#include <optional>
#include <string>
#include <vector>

#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

#include "compute-bam-features/alignment-reader.h"
#include "compute-bam-features/compute-bam-region-features.h"
#include "core/config.h"
#include "core/score-calculator.h"
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
  vec<io::InfoFieldMetadata> vcf_info_metadata{};
  vec<io::FormatFieldMetadata> vcf_fmt_metadata{};
  std::optional<s32> vcf_normal_index{};
  std::optional<s32> vcf_tumor_index{};
  ChromMedianDepth normalize_targets{};
  sex_predict::Sex sex{};
  std::string chr_x_name{};
  std::string chr_y_name{};
  vec<Interval> chr_x_par{};
  vec<Interval> chr_y_par{};
  bool phased{};
  std::optional<StrMap<vec<Interval>>> force_calls{};
  std::optional<StrUnorderedSet> hotspots{};
  std::optional<StrUnorderedSet> block_list{};
  f32 min_allele_freq_threshold{};
  f32 weighted_counts_threshold{};
  f32 hotspot_weighted_counts_threshold{};
  f32 ml_threshold{};
  f32 somatic_tn_snv_ml_threshold{};
  f32 somatic_tn_indel_ml_threshold{};
  f32 hotspot_ml_threshold{};
  f32 germline_fail_ml_threshold{};
  f32 min_phased_allele_freq{};
  f32 max_phased_allele_freq{};
  u32 min_alt_counts{};
  std::optional<fs::path> skip_variants_vcf{};
  u32 tumor_support_threshold{};
  bool is_germline_tagging{};
};

std::tuple<VarIdToVcfFeatures, BamRegionFeatureCollection> ComputeBamAndVcfFeaturesForRegion(
    const GlobalContext& global_ctx,
    WorkerContext& worker_ctx,
    const TargetRegion& region,
    const ChromIntervalsMap& bed_regions,
    const ChromIntervalsMap& interest_regions);

vec<TargetRegion> PartitionVcfRegions(const fs::path& vcf_file,
                                      size_t threads,
                                      std::optional<ChromIntervalsMap> bed_regions);

}  // namespace xoos::svc
