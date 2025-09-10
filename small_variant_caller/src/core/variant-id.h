#pragma once
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

  u64 GetRefFeaturePos() const;
  u64 GetMaxRefAllelePos() const;
};

}  // namespace xoos::svc
