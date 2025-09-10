#include "variant-type.h"

namespace xoos::svc {

VariantType GetVariantType(const std::string& ref, const std::string& alt) {
  using enum VariantType;
  const auto ref_length = ref.length();
  const auto alt_length = alt.length();
  if (ref_length == 1 && alt_length == 1) {
    return kSNV;
  }
  if (ref_length < alt_length) {
    return kInsertion;
  }
  if (ref_length > alt_length) {
    return kDeletion;
  }
  return kUnknown;
}

}  // namespace xoos::svc
