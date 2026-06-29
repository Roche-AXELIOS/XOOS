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

bool CheckVariantType(const VariantGroup var_group, const VariantType vid_type) {
  if (var_group == VariantGroup::kAll) {
    return true;
  }
  if (vid_type == VariantType::kUnknown) {
    return false;
  }
  if (vid_type == VariantType::kSNV) {
    return var_group == VariantGroup::kSnvOnly;
  }
  if (vid_type == VariantType::kInsertion || vid_type == VariantType::kDeletion) {
    return var_group == VariantGroup::kIndelOnly;
  }
  return false;
}

}  // namespace xoos::svc
