#include "variant-id.h"

#include "util/seq-util.h"
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

std::string VariantId::ToString() const {
  // Variant positions are 1-based in string representation
  return fmt::format("{}:{}:{}>{}", chrom, pos + 1, ref, alt);
}

u64 VariantId::GetRefFeaturePos() const {
  if (type == VariantType::kDeletion) {
    // To get reference features at the first variant base of a deletion, `id.pos` must be offset by 1.
    return pos + 1;
  }
  return pos;
}

u64 VariantId::GetMaxRefAllelePos() const {
  if (ref.empty()) {
    throw error::Error("Reference allele is empty");
  }
  return pos + ref.length() - 1;
}

}  // namespace xoos::svc

// Hash function implementation for VariantId using SeqToBits for ref and alt
size_t std::hash<xoos::svc::VariantId>::operator()(const xoos::svc::VariantId& vid) const noexcept {
  const size_t h1 = std::hash<std::string>{}(vid.chrom);
  const size_t h2 = std::hash<uint64_t>{}(vid.pos);
  const size_t h3 = xoos::svc::SeqToBits(vid.ref);
  const size_t h4 = xoos::svc::SeqToBits(vid.alt);
  const size_t h5 = std::hash<uint64_t>{}(static_cast<uint64_t>(vid.type));
  return (((h1 ^ (h2 << 1)) ^ (h3 << 2)) ^ (h4 << 3)) ^ (h5 << 4);
}
