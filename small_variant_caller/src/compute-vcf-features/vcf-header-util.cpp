#include "vcf-header-util.h"

namespace xoos::svc {

std::optional<TumorNormalSampleIndexes> GetTumorNormalSampleIndexes(const io::VcfHeaderPtr& hdr) {
  // extract sample indexes for `tumor_sample` and `normal_sample` if they are found in the header
  const std::map<std::string, s32>& sample_name_to_idx{hdr->GetSampleIndexes()};
  const char* const normal_name = hdr->GetValue("normal_sample");
  const char* const tumor_name = hdr->GetValue("tumor_sample");
  if (normal_name != nullptr && tumor_name != nullptr) {
    const s32 normal_idx = sample_name_to_idx.at(normal_name);
    const s32 tumor_idx = sample_name_to_idx.at(tumor_name);
    if (normal_idx >= 0 && tumor_idx >= 0 && tumor_idx != normal_idx) {
      return TumorNormalSampleIndexes{.tumor_sample_idx = tumor_idx, .normal_sample_idx = normal_idx};
    }
  }
  return std::nullopt;
}

}  // namespace xoos::svc
