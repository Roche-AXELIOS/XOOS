#include "variant-id.h"

#include "xoos/error/error.h"

namespace xoos::svc {

VariantId::VariantId(std::string chrom, const u64 pos, std::string ref, std::string alt) {
  this->type = GetVariantType(ref, alt);
  this->chrom = std::move(chrom);
  this->pos = pos;
  this->ref = std::move(ref);
  this->alt = std::move(alt);
}

bool VariantId::operator<(const VariantId& other) const {
  return std::tie(chrom, pos, ref, alt, type) < std::tie(other.chrom, other.pos, other.ref, other.alt, other.type);
}

/**
 * @brief Returns the position of the reference feature for a given variant. The position is adjusted for deletions.
 * @return position of the reference feature
 */
u64 VariantId::GetRefFeaturePos() const {
  if (type == VariantType::kDeletion) {
    // To get reference features at the first variant base of a deletion, `id.pos` must be offset by 1.
    return pos + 1;
  }
  return pos;
}

/**
 * @brief Maximum position of the reference allele in a variant.
 * @return Maximum position of the reference allele
 * @throws error::Error if the reference allele is empty
 */
u64 VariantId::GetMaxRefAllelePos() const {
  if (ref.empty()) {
    throw error::Error("Reference allele is empty");
  }
  return pos + ref.length() - 1;
}

}  // namespace xoos::svc
