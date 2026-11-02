#include "variant-id.h"

#include "util/seq-util.h"
#include "xoos/error/error.h"

namespace xoos::svc {

VariantId::VariantId(std::string chrom, const u64 pos, std::string ref, std::string alt)
    : chrom{std::move(chrom)},
      pos{pos},
      ref{std::move(ref)},
      alt{std::move(alt)},
      // Initialized last (after ref/alt members are set) so it derives from the moved-in values.
      type{GetVariantType(this->ref, this->alt)} {
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

u64 VariantId::GetAggFeaturePos() const {
  // Per-position aggregate features are keyed by the raw reported position for all variant types (no deletion offset),
  // so co-located SNVs, insertions, and deletions share the same aggregate values.
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
  constexpr size_t kHashShift1 = 1;
  constexpr size_t kHashShift2 = 2;
  constexpr size_t kHashShift3 = 3;
  constexpr size_t kHashShift4 = 4;
  return (((h1 ^ (h2 << kHashShift1)) ^ (h3 << kHashShift2)) ^ (h4 << kHashShift3)) ^ (h5 << kHashShift4);
}
