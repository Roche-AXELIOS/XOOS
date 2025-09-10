#pragma once

#include <string>

#include <xoos/types/int.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "compute-bam-region-features.h"
#include "core/variant-info.h"
#include "region.h"
#include "xoos/types/float.h"

namespace xoos::svc {

// Maximum possible base quality and mapping quality values to calculate the `weighted_depth` feature
// TODO: confirm these values, how they were derived, and document them here
constexpr u8 kMaxPossibleBaseQual = 138;
constexpr u8 kMaxPossibleMapQual = 60;

// minimum base quality for concordant bases in a duplex read; works for concordant base in both 0,18,93 or 5,22,39
constexpr u8 kMinConcordantBaseQuality{30};

/**
 * @brief Sequence alignment operations describing the relationship between a read and a reference sequence.
 * @note Alignment clipping (e.g., soft and hard clips) are not included because they do not contribute to the
 *       feature extraction process.
 * @details Possible values:
 * - `kUnknown`: Unknown or unset alignment operation type.
 * - `kReferenceSkip`: A portion of the reference sequence is skipped in the alignment.
 * - `kMatch`: The read base matches the reference base.
 * - `kMismatch`: The read base does not match the reference base.
 * - `kInsertion`: An insertion in the read relative to the reference sequence.
 * - `kDeletion`: A deletion in the read relative to the reference sequence.
 */
enum class AlignOp {
  kUnknown,
  kReferenceSkip,
  kMatch,
  kMismatch,
  kInsertion,
  kDeletion,
};

/**
 * @brief Information about an alignment operation.
 * @see AlignOp
 */
struct AlignOpInfo {
  // alignment operation
  AlignOp op{AlignOp::kUnknown};
  // length of the operation
  u64 op_len{};
  // 0-based start position in the reference sequence for the operation
  u64 ref_pos{};
  // 0-based start position in the read for the operation
  u64 read_pos{};
  auto operator<=>(const AlignOpInfo&) const = default;
};

/**
 * @brief Read type to distinguish which feature values to be incremented.
 */
enum class ReadType {
  kDuplex,
  kSimplex,
  kUmiConsensus
};

/**
 * @brief Context for parameters and alignment operation information needed to determine whether an alignment
 * operation can be processed for feature extraction.
 * @details
 * The context is initialized with CLI parameters for feature extraction, alignment operation information, BAM record,
 * region, reference sequence, read ID, and VCF features.
 * Using these information, the following information are extracted and stored:
 * - family size for reads in the plus and minus strands
 * - read type (duplex, simplex, or UMI consensus)
 * - minimum base type required for processing this read, if applicable
 * - whether the read originated from a tumor sample, if applicable
 * - flags indicating whether there are no VCF features for the region
 * - vectors for alignment operation information, base types, and flags indicating whether each position is near
 *   a non-concordant base
 * It includes methods to check if an insertion, deletion, mismatch, or match can be processed
 * and whether they are near a non-concordant base.
 * @see ComputeBamFeaturesParams for CLI parameters used in feature extraction.
 * @see AlignOpInfo for alignment operation information.
 * @see bam1_t for the BAM record structure.
 * @see Region for the chromosomal region being processed.
 * @see PositionToVcfFeaturesMap for mapping positions to VCF features.
 * @see yc_decode::BaseType for base type information if yc tag decoding is used.
 * @see ReadId for the numeric ID of the read.
 * @see ReadType for the type of the read (duplex, simplex, or UMI consensus).
 * @see UnifiedVariantFeatures and UnifiedReferenceFeatures for variant and reference allele features.
 */
class AlignContext {
 public:
  AlignContext(const ComputeBamFeaturesParams& params,
               const vec<AlignOpInfo>& align_op_infos,
               const bam1_t* bam_record,
               const Region& region,
               const std::string& ref_seq,
               ReadId read_id,
               const PositionToVcfFeaturesMap& vcf_feats);

  // functions to check whether an alignment operation can be processed for feature extraction
  bool CanProcessInsertion(const AlignOpInfo& info) const;
  bool CanProcessDeletion(const AlignOpInfo& info) const;
  bool CanProcessMismatch(const AlignOpInfo& info) const;
  bool CanProcessMatch(u64 read_pos, u64 ref_pos) const;

  // functions to check whether an alignment operation is near a non-concordant base
  bool IsInsertionNearNonConcordantBase(const AlignOpInfo& info) const;
  bool IsDeletionNearNonConcordantBase(const AlignOpInfo& info) const;
  bool IsMismatchNearNonConcordantBase(const AlignOpInfo& info) const;
  bool IsMatchNearNonConcordantBase(u64 read_pos) const;

  // CLI parameters for feature extraction
  const ComputeBamFeaturesParams& params;
  // alignment operation information for the read
  const vec<AlignOpInfo>& align_op_infos;
  // BAM record for the read
  const bam1_t* bam_record;
  // chromosomal region being processed
  const Region& region;
  // chromosome sequence for the region
  const std::string& ref_seq;
  // numeric ID for the read
  const ReadId read_id;
  // Map of positions to VCF features for the region
  const PositionToVcfFeaturesMap& vcf_feats;
  // Flag indicating whether there are VCF features for the region
  const bool has_vcf_feats;
  // family size for read(s) in the plus strand
  u32 plus_counts = 0;
  // family size for read(s) in the minus strand
  u32 minus_counts = 0;
  // type of the read (duplex, simplex, or UMI consensus)
  ReadType read_type;
  // minimum base type required for processing this read; empty if not applicable
  std::optional<yc_decode::BaseType> min_base_type;
  // flag indicating whether this read originated from the tumor sample; empty if not applicable
  std::optional<bool> is_tumor_sample;
  // 0-based position of the first base after the last aligned reference base
  u64 ref_end = 0;
  // 0-based position of the first base of the homopolymer that contains the last aligned reference base
  std::optional<u64> homopolymer_start = std::nullopt;
  // vector of base types for each position in the read; empty if yc tag decoding is not used
  vec<yc_decode::BaseType> base_types;
  // vector of bool flagging whether each position is near a non-concordant base
  vec<bool> near_non_concordant;
  // flag indicating whether to skip processing this read
  bool skip_read = false;
};

/**
 * Parse the alignment CIGAR of a read to reference alignment to generate a vector of alignment operations that can be
 * easily referenced
 * @param record An htslib (bam1_t*) bam record struct
 * @param target_seq The reference sequence of the chromosome the alignment corresponds to
 * @return A vector of AlignInfo structs
 * @note The output vector will always start and end with a Match or Mismatch operation. All leading and trailing
 * insertion, deletion, soft-clip and hard clip are skipped.
 */
vec<AlignOpInfo> GetAlignOpInfos(const bam1_t* record, const std::string& target_seq);

/**
 * @brief Get the alignment length based on match, mismatch, and insertion only.
 * The is essentially the length of the read from the first aligned base to the last aligned base.
 * @param infos Vector of alignment operator information
 * @return Alignment length
 */
u64 GetAlignmentLength(const vec<AlignOpInfo>& infos);

/**
 * @brief Process alignment operations from a read and extract unified variant features.
 * @param align_ctx Context information for the read alignment
 * @param var_feats Map of variant features to be populated
 * @param ref_feats Map of reference allele features to be populated
 */
void ProcessAlignment(const AlignContext& align_ctx,
                      UnifiedVariantFeatures& var_feats,
                      UnifiedReferenceFeatures& ref_feats);

/**
 * @brief Update derived feature values based on the counts and sums in the UnifiedVariantFeature struct.
 * @param feature struct to be updated
 * @note Feature values depending on the REF allele or the other ALT allele at the same position are not updated.
 */
void UpdateDerivedFeatureValues(UnifiedVariantFeature& feature);

/**
 * @brief Update derived feature values based on the counts and sums in the UnifiedReferenceFeature struct.
 * @param feature struct to be updated
 */
void UpdateDerivedFeatureValues(UnifiedReferenceFeature& feature);

/**
 * @brief Increment unified variant feature.
 * @param feature A unified variant feature to be incremented
 * @param baseq Read base quality supporting the variant.
 * @param distance Minimum distance of variant to alignment ends.
 * @param near_non_concordant Flag whether this position is near a non-concordant base
 */
void IncrementFeature(
    UnifiedVariantFeature& feature, const AlignContext& align_ctx, f64 baseq, u64 distance, bool near_non_concordant);

/**
 * @brief Extract variant allele features for an alignment insertion operation.
 * @param features Map of position to variant features
 * @param align_ctx Context information for extracting features from the alignment
 * @param info Alignment info for the insertion operation
 */
void ProcessInsertion(UnifiedVariantFeatures& features, const AlignContext& align_ctx, const AlignOpInfo& info);

/**
 * @brief Extract variant allele features for an alignment deletion operation.
 * @param features Map of position to variant features
 * @param align_ctx Context information for extracting features from the alignment
 * @param info Alignment info for the deletion operation
 */
void ProcessDeletion(UnifiedVariantFeatures& features, const AlignContext& align_ctx, const AlignOpInfo& info);

/**
 * @brief Extract variant allele features for a mismatch base.
 * @param features Map of position to variant features
 * @param align_ctx Context information for extracting features from the alignment
 * @param info Alignment info for the mismatch base
 */
void ProcessMismatch(UnifiedVariantFeatures& features, const AlignContext& align_ctx, const AlignOpInfo& info);

/**
 * Extract reference allele features for a matching base.
 * @param features Map of chromosome position to reference allele feature
 * @param align_ctx Context information for extracting features from the alignment
 * @param read_pos 0-based position of the matching base on the read sequence
 * @param ref_pos 0-based position of the matching base on the reference sequence
 */
void ProcessMatch(UnifiedReferenceFeatures& features, const AlignContext& align_ctx, u64 read_pos, u64 ref_pos);

/**
 * @brief Count variants based on alignment operators.
 * @param infos Vector of alignment operator information
 * @return Number of variants
 */
u32 CountVariantsInRead(const vec<AlignOpInfo>& infos);

/**
 * @brief Count variants based on alignment operators, skipping any known variants.
 * @param infos Vector of alignment operator information
 * @param record BAM record
 * @param ref_seq Reference sequence
 * @param chrom Chromosome name
 * @param skip_variants Set of variant keys to skip
 * @return Number of variants
 */
u32 CountVariantsInRead(const vec<AlignOpInfo>& infos,
                        const bam1_t* record,
                        const std::string& ref_seq,
                        const std::string& chrom,
                        const StrUnorderedSet& skip_variants);

/**
 * Extract the 1-bp reference context before and after a given position.
 * @param ref_seq reference sequence
 * @param pos 0-based position on the reference sequence
 * @return 1-bp context concatenated without any delimiters
 * @note Context containing non-ACGT bases is not supported and will return an empty string.
 * @note Context bases must be within the bounds of the reference sequence (i.e., 0 < pos < ref_seq.size() - 1).
 */
std::string Get1bpContext(const std::string& ref_seq, u64 pos);

/**
 * Determine whether each base in the read is nearby a non-concordant base.
 * @param alignment BAM record
 * @return results as a boolean vector
 */
vec<bool> FindNearbyNonConcordantBase(const bam1_t* alignment);

/**
 * Determine whether each base in the read is nearby a non-concordant base.
 * @param alignment BAM record
 * @param base_types Vector of BaseType for each base in the read
 * @return results as a boolean vector
 */
vec<bool> FindNearbyNonConcordantBase(const bam1_t* alignment, const vec<yc_decode::BaseType>& base_types);

}  // namespace xoos::svc
