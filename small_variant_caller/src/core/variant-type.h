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

/**
 * @brief Determine the type of variant based on its reference and alternate alleles.
 * @param ref Reference allele
 * @param alt Alternate allele
 * @return VariantType indicating whether the variant is a substitution, insertion, deletion, or unknown
 * @note This function assumes that the input alleles are valid nucleotide sequences in the minimal representation.
 * @see `TrimVariant` for details on how to trim alleles.
 */
VariantType GetVariantType(const std::string& ref, const std::string& alt);

}  // namespace xoos::svc
