#include "variant-type.h"

namespace xoos::svc {

VariantType GetVariantType(const std::string_view ref, const std::string_view alt) {
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
  using enum VariantGroup;
  using enum VariantType;
  if (var_group == kAll) {
    return true;
  }
  if (vid_type == kUnknown) {
    return false;
  }
  if (vid_type == kSNV) {  // NOSONAR: standard genomics abbreviation
    return var_group == kSnvOnly;
  }
  if (vid_type == kInsertion || vid_type == kDeletion) {
    return var_group == kIndelOnly;
  }
  return false;
}

}  // namespace xoos::svc
