#pragma once
#include <string>

namespace xoos::svc {

/**
 * @brief Enum representing the type of genomic variant.
 * @note MNVs (multi-nucleotide variants) are treated as "unknown". They should be decomposed into multiple SNVs.
 */
enum class VariantType {
  kSNV,
  kDeletion,
  kInsertion,
  kUnknown
};

enum class VariantGroup {
  kAll,
  kSnvOnly,
  kIndelOnly
};

/**
 * @brief Determine the type of variant based on its reference and alternate alleles.
 * @param ref Reference allele
 * @param alt Alternate allele
 * @return VariantType indicating whether the variant is a substitution, insertion, deletion, or unknown
 * @note This function assumes that the input alleles are valid nucleotide sequences in the minimal representation.
 * @see `TrimVariant` for details on how to trim alleles.
 */
VariantType GetVariantType(const std::string& ref, const std::string& alt);

/**
 * @brief Check if the variant type matches the specified variant group.
  If var_group is kAll, all variant types are accepted.
  If var_group is kSnvOnly, only SNVs are accepted.
  If var_group is kIndelOnly, only indels are accepted.
 * @param var_group the type of variant to train
 * @param vid_type the type of the variant
 * @return True if the variant type and variant group align for training, false otherwise
 */
bool CheckVariantType(VariantGroup var_group, VariantType vid_type);

}  // namespace xoos::svc
