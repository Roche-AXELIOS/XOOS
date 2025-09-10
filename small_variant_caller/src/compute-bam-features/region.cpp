#include "region.h"

#include <xoos/io/vcf/vcf-reader.h>

#include "core/filtering.h"
#include "util/seq-util.h"

namespace xoos::svc {

StrUnorderedSet ExtractVariantKeySet(const fs::path& vcf_file, const Region& region, u64 region_padding) {
  io::VcfReader reader(vcf_file);

  const auto chrom_len = reader.GetHeader()->GetContigLengths().at(region.chrom);
  const auto start = std::max(region.start, region_padding) - region_padding;
  const auto end = std::min(region.end + region_padding, chrom_len);
  reader.SetRegion(region.chrom, start, end);

  StrUnorderedSet variants;
  while (const auto& vcf_record = reader.GetNextRegionRecord(BCF_UN_STR)) {
    const auto& ref = vcf_record->Allele(0);
    if (IsAnyNotACTG(ref)) {
      continue;
    }

    for (int i = 1; i < vcf_record->NumAlleles(); ++i) {
      const auto& alt = vcf_record->Allele(i);
      if (IsAnyNotACTG(alt)) {
        continue;
      }
      variants.insert(GetVariantCorrelationKey(region.chrom, vcf_record->Position(), ref, alt, false));
    }
  }
  return variants;
}

s64 Region::EndInclusive() const {
  return static_cast<s64>(end + 1);
}

s64 Region::Start() const {
  return static_cast<s64>(start);
}
}  // namespace xoos::svc
