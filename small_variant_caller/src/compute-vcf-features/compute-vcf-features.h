#pragma once

#include <filesystem>
#include <map>
#include <optional>
#include <string>

#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

#include "core/config.h"
#include "core/variant-info.h"
#include "core/workflow.h"
#include "util/region-util.h"

namespace xoos::svc {

/**
 * @brief CLI parameters to compute features from a VCF file.
 */
struct ComputeVcfFeaturesParam {
  fs::path vcf_file{};                   // Path of input VCF file
  fs::path genome{};                     // Path of input genome FASTA file
  fs::path output_file{};                // Path of output TSV file for computed features
  std::optional<fs::path> pop_af_vcf{};  // Path of input population allele frequency VCF file
  SVCConfig config{};                    // Configuration from custom JSON file or default settings
  Workflow workflow{};                   // Workflow to execute (e.g., germline, etc.)
  std::optional<ChromIntervalsMap>
      target_regions{};  // Map of chromosome to target regions of interest for feature extraction
  std::optional<ChromIntervalsMap>
      interest_regions{};                // Map of chromosome to interest regions for `at_interest_region` feature
  size_t threads{};                      // Number of threads to use for parallel processing
  std::optional<fs::path> output_bed{};  // Path of output BED file for regions with computed features
  u64 left_pad{};       // Number of bases to include before the variant's start position for output BED
  u64 right_pad{};      // Number of bases to include after the variant's start position for output BED
  u32 collapse_dist{};  // Distance for collapsing nearby intervals for output BED
};

using ChromToIndexes = std::map<std::string, s32>;
using ChromToLengths = std::map<std::string, u64>;

/**
 * @brief Struct to hold information extracted from the VCF header.
 * It includes contig lengths, indexes, and flags indicating the presence of specific INFO/FORMAT fields.
 */
struct VcfHeaderInfo {
  ChromToLengths contig_lengths{};
  ChromToIndexes contig_indexes{};
  // flag whether INFO/FORMAT fields are present in the VCF header
  bool has_popaf{false};
  bool has_nalod{false};
  bool has_nlod{false};
  bool has_tlod{false};
  bool has_mpos{false};
  bool has_mmq{false};
  bool has_mbq{false};
  bool has_ad{false};
  bool has_af{false};
  bool has_dp{false};
  bool has_gq{false};
  bool has_gt{false};
  bool has_hapcomp{false};
  bool has_hapdom{false};
  bool has_ru{false};
  bool has_rpa{false};
  bool has_str{false};
  bool has_tumor_normal{false};  // flag whether the VCF file contains tumor-normal pair samples
  int normal_index{-1};          // index of the normal sample in the VCF records
  int tumor_index{-1};           // index of the tumor sample in the VCF records
};

/**
 * @brief Check feature column names against available resources and VCF fields; produce warnings for the absence
 * of resources or VCF fields.
 * @param feat_cols Vector of feature column enums
 * @param interest_regions Map of chromosome to vector of repeat region intervals
 * @param header VCF header shared pointer
 * @param has_popaf_vcf Flag whether POPAF VCF file is available
 * @return Number of warnings reported
 */
u16 CheckVcfFeatureResources(const vec<UnifiedFeatureCols>& feat_cols,
                             const ChromIntervalsMap& interest_regions,
                             const std::shared_ptr<io::VcfHeader>& header,
                             bool has_popaf_vcf);

/**
 * @brief Extract chromosome-level median FORMAT field DP from VCF.
 * @param vcf_path Path to VCF
 * @return Map of chromosome to median DP
 * @note Assumes that the VCF file is sorted by chromosome and position.
 */
StrUnorderedMap<u32> GetChromosomeMedianDP(const fs::path& vcf_path);

/**
 * @brief Structure to hold counts of different types of tandem repeats.
 */
struct RepeatCounts {
  u32 homopolymer;
  u32 direpeat;
  u32 trirepeat;
  u32 quadrepeat;
};

/**
 * @brief Extract repeat counts from a nucleotide sequence.
 * @param seq Nucleotide sequence
 * @return Repeat counts for homopolymer, di-repeat, tri-repeat, and quad-repeat
 * @detail Handle ambiguous repeat unit representations:
 *   Example: "AAATAAAT" has 3x "A" (homopolymer) and 2x "AAAT" (quad-repeat)
 *   - 2x "AAAT" covers 8 bases, while 3x "A" covers only 3 bases.
 *   - Choose the repeat type that covers the most bases (2x "AAAT" in this case).
 *   Example: "ATATATAT" can be represented as either 4x "AT" (di-repeat) or 2x "ATAT" (quad-repeat)
 *   - Both repeat types cover the same number of bases, i.e. 4 x 2bp == 2 x 4bp.
 *   - Choose the repeat type with the largest number of repetitions (4x "AT" in this case).
 */
RepeatCounts GetRepeatCounts(const std::string& seq);

/**
 * @brief Extract relevant information from VCF header.
 * @param header Shared pointer to VcfHeader
 * @return Extracted information
 */
VcfHeaderInfo GetVcfHeaderInfo(const std::shared_ptr<io::VcfHeader>& header);

/**
 * @brief Compute VCF features from a VCF file and reference genome FASTA file.
 * @param param Parameters for computing VCF features
 * @note This function extracts features from the input VCF file and reference genome FASTA file,
 * updates variant density for all extracted features, updates POPAF values if available,
 * writes feature positions to a BED file if specified, and writes features to an output TSV file if specified.
 */
void ComputeVcfFeatures(const ComputeVcfFeaturesParam& param);

/**
 * @brief Extract VCF features for target regions in the same chromosome, update POPAF and variant density.
 * @param vcf_reader VCF reader
 * @param genome_path Path to reference genome FASTA file
 * @param popaf_reader POPAF VCF reader
 * @param target_regions Map of chromosome to vector of target intervals
 * @param interest_regions Map of chromosome to vector of repeat intervals
 * @param chrom Chromsome name of target regions
 * @return PositionToVcfFeaturesMap mapping positions to VcfFeature structs
 * @throws std::runtime_error if the chromosome is not found in the VCF header
 * @note This function extracts features for a single chromosome, specified by `chrom`.
 */
PositionToVcfFeaturesMap ExtractFeaturesForRegion(io::VcfReader& vcf_reader,
                                                  const fs::path& genome_path,
                                                  std::optional<io::VcfReader>& popaf_reader,
                                                  const ChromIntervalsMap& target_regions,
                                                  const ChromIntervalsMap& interest_regions,
                                                  const std::string& chrom);

}  // namespace xoos::svc
