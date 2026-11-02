#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <string_view>

#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

#include "core/config.h"
#include "core/library-type.h"
#include "core/variant-info.h"
#include "core/workflow.h"
#include "util/region-util.h"

namespace xoos::svc {

/**
 * @brief CLI parameters to compute features from a VCF file.
 * @details Germline WGS requires an explicit library type because this command has no BAM from which to detect it.
 */
struct ComputeVcfFeaturesParam {
  // Path of input VCF file
  fs::path vcf_file;
  // Path of input genome FASTA file
  fs::path genome;
  // Path of output TSV file for computed features
  fs::path output_file;
  // Path of input population allele frequency VCF file
  std::optional<fs::path> pop_af_vcf;
  // Path of input Panel of Normals VCF file
  std::optional<fs::path> pon_vcf;
  // Path of input catalog VCF file for the `in_catalog` feature
  std::optional<fs::path> catalog_vcf;
  // Configuration from custom JSON file or default settings
  SVCConfig config;
  // Path of config JSON file
  std::optional<fs::path> config_file;
  // Workflow to execute (e.g., germline, etc.)
  Workflow workflow{};
  // Defaults consistently with BAM-based commands, but validation requires a concrete value because no BAM is read.
  LibraryType library_type{LibraryType::kAutoDetect};
  // Map of chromosome to target regions of interest for feature extraction
  std::optional<ChromIntervalsMap> target_regions;
  // Map of chromosome to interest regions for `at_interest_region` feature
  std::optional<ChromIntervalsMap> interest_regions;
  // Number of threads to use for parallel processing
  size_t threads{};
  // Path of output BED file for regions with computed features
  std::optional<fs::path> output_bed;
  // Number of bases to include before the variant's start position for output BED
  u64 left_pad{};
  // Number of bases to include after the variant's start position for output BED
  u64 right_pad{};
  // Distance for collapsing nearby intervals for output BED
  u32 collapse_dist{};
  std::optional<io::CommandLineInfo> command_line;
};

using ChromToIndexes = xoos::StrMap<s32>;
using ChromToLengths = xoos::StrMap<u64>;

/**
 * @brief Struct to hold information extracted from the VCF header.
 * It includes contig lengths, indexes, and flags indicating the presence of specific INFO/FORMAT fields.
 */
struct VcfHeaderInfo {
  ChromToLengths contig_lengths;
  ChromToIndexes contig_indexes;
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
  bool has_mq_rank_sum{false};
  bool has_sor{false};
  bool has_ru{false};
  bool has_rpa{false};
  bool has_str{false};
  bool has_tumor_normal{false};  // flag whether the VCF file contains tumor-normal pair samples
  s32 normal_index{-1};          // index of the normal sample in the VCF records
  s32 tumor_index{-1};           // index of the tumor sample in the VCF records
};

/**
 * @brief Check if a VCF record has the PASS filter.
 * @param record VCF record pointer
 * @return True if the record has the PASS filter, false otherwise
 */
bool HasPassFilter(const io::VcfRecordPtr& record);

/**
 * @brief Check feature column names against available resources and VCF fields; produce warnings for the absence
 * of resources or VCF fields.
 * @param feat_cols Vector of feature column
 * @param interest_regions Map of chromosome to vector of repeat region intervals
 * @param header VCF header shared pointer
 * @param has_popaf_vcf Flag whether POPAF VCF file is available
 * @param has_catalog_vcf Flag whether catalog VCF file is available
 * @return Number of warnings reported
 */
u16 CheckVcfFeatureResources(const vec<FeatureColumn>& feat_cols,
                             const ChromIntervalsMap& interest_regions,
                             const std::shared_ptr<io::VcfHeader>& header,
                             bool has_popaf_vcf,
                             bool has_catalog_vcf);

/**
 * @brief Compute the start and end coordinates of a GC content window centered on a position.
 * @param pos Variant position (0-based)
 * @param seq_len Length of the reference sequence (chromosome)
 * @param half_window Half-window size (number of bases upstream and downstream)
 * @return Pair of (start, end) where start is inclusive and end is exclusive, clamped to [0, seq_len]
 */
std::pair<u64, u64> GetGcContentWindow(u64 pos, u64 seq_len, u64 half_window);

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
RepeatCounts GetRepeatCounts(std::string_view seq);

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

// Forward declaration — full definition in parallel-compute-utils.h
struct VariantPreScanEntry;

/**
 * @brief Additional parameters for VCF feature extraction beyond the core reader and region inputs.
 * @details Bundles the germline tagging flag and optional pre-scan density data to keep the
 * ExtractFeaturesForRegion parameter count within limits.
 */
struct VcfFeatureExtractionParams {
  bool is_germline_tagging{false};
  const StrMap<vec<VariantPreScanEntry>>* pre_scan_variants{nullptr};
  const StrMap<vec<u32>>* pre_scan_density{nullptr};
};

/**
 * @brief Non-owning view of the optional external-VCF readers used to annotate features.
 * @details Bundles the POPAF, Panel of Normals (PON), and catalog readers into a single parameter to
 * avoid passing three adjacent same-type arguments. The referenced optionals remain owned by the caller
 * (e.g. WorkerContext) and stay mutable through the references. Each pointer must be non-null before
 * passing to ExtractFeaturesForRegion; all members are default-initialized to nullptr.
 */
struct AnnotationVcfReaders {
  std::optional<io::VcfReader>* popaf{nullptr};
  std::optional<io::VcfReader>* pon{nullptr};
  std::optional<io::VcfReader>* catalog{nullptr};
};

/**
 * @brief Extract VCF features for target regions in the same chromosome.
 * @details Computes per-variant features and variant density, and optionally annotates them from POPAF,
 * Panel of Normals (PON), and catalog VCFs when the corresponding readers are provided.
 * @param vcf_reader VCF reader
 * @param genome_path Path to reference genome FASTA file
 * @param annotation_readers Optional POPAF/PON/catalog readers used to annotate the extracted features
 * @param target_regions Map of chromosome to vector of target intervals
 * @param interest_regions Map of chromosome to vector of repeat intervals
 * @param chrom Chromosome name of target regions
 * @param params Extraction parameters (germline tagging flag and optional pre-scan density context)
 * @return VcfFeaturesMap mapping VariantId to VcfFeature structs
 * @throws std::runtime_error if the chromosome is not found in the VCF header
 * @note When params contains valid pre-scan pointers, variant density is looked up from the pre-scan data
 * instead of being computed per-region, avoiding boundary artifacts at region edges.
 */
VarIdToVcfFeatures ExtractFeaturesForRegion(io::VcfReader& vcf_reader,
                                            const fs::path& genome_path,
                                            const AnnotationVcfReaders& annotation_readers,
                                            const ChromIntervalsMap& target_regions,
                                            const ChromIntervalsMap& interest_regions,
                                            const std::string& chrom,
                                            const VcfFeatureExtractionParams& params);

}  // namespace xoos::svc
