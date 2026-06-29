#pragma once
#include <functional>
#include <string>

#include "variant-type.h"
#include "xoos/types/int.h"

namespace xoos::svc {

/**
 * @brief Information required to uniquely identify a variant.
 * It contains the chromosome, position, reference allele, alternate allele, and variant type.
 */
struct VariantId {
 public:
  std::string chrom;
  u64 pos{0};
  std::string ref;
  std::string alt;
  VariantType type{VariantType::kUnknown};

  VariantId() = default;
  VariantId(std::string chrom, u64 pos, std::string ref, std::string alt);

  bool operator<(const VariantId& other) const;
  bool operator==(const VariantId& other) const = default;

  /**
   * @brief Converts the VariantId to a string representation in the format "chrom:pos:ref>alt".
   * @return string representation of the VariantId
   * @note Position is 1-based in the string representation.
   */
  std::string ToString() const;

  /**
   * @brief Returns the position of the reference feature for a given variant. The position is adjusted for deletions.
   * @return position of the reference feature
   */
  u64 GetRefFeaturePos() const;

  /**
   * @brief Maximum position of the reference allele in a variant.
   * @return Maximum position of the reference allele
   * @throws error::Error if the reference allele is empty
   */
  u64 GetMaxRefAllelePos() const;
};

}  // namespace xoos::svc

/**
 * @brief Custom hash function for VariantId to allow its use as a key in unordered containers.
 */
template <>
struct std::hash<xoos::svc::VariantId> {
  size_t operator()(const xoos::svc::VariantId& vid) const noexcept;
};
