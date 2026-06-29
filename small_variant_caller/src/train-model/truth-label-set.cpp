#include "truth-label-set.h"

#include <ranges>

#include "core/genotype.h"
#include "util/seq-util.h"
#include "xoos/io/vcf/vcf-reader.h"

namespace xoos::svc {

GenotypeToVariantIds GetGenotypeToVariantIds(const fs::path& vcf) {
  GenotypeToVariantIds results;
  const io::VcfReader reader(vcf);
  while (const auto& vcf_record = reader.GetNextRecord()) {
    const auto& ref = vcf_record->Allele(0);
    if (!ContainsOnlyACTG(ref)) {
      continue;
    }
    // Unsupported genotypes in the ground truth are still kept because variants with these genotypes should be
    // excluded from the training data.
    const auto genotype = StringToGenotype(vcf_record->GetGTField());
    // Multi-allelics are rare compared to heterozygous REF-ALT and homozygous ALT.
    // To include as many true alleles as possible in the training data, store the ground truth GT label for each
    // individual ALT allele separately.
    // This is intended to handle two scenarios:
    // 1. Multi-allelics with 1 true allele and 1 false allele.
    // 2. Multi-allelics with 1 SNV and 1 indel. The germline workflow has separate training data for SNV and indel.
    const s32 num_alleles = vcf_record->NumAlleles();
    const auto chrom = vcf_record->Chromosome();
    const auto pos = static_cast<u64>(vcf_record->Position());
    for (s32 i = 1; i < num_alleles; ++i) {
      const auto& [ref_trimmed, alt_trimmed] = TrimVariant(ref, vcf_record->Allele(i));
      if (ContainsOnlyACTG(alt_trimmed)) {
        const VariantId vid(chrom, pos, ref_trimmed, alt_trimmed);
        results[genotype].emplace(vid);
      }
    }
  }
  return results;
}

std::optional<Genotype> GetGenotypeForVariant(const VariantId& vid, const GenotypeToVariantIds& genotypes) {
  for (const auto& [genotype, vids] : genotypes) {
    if (vids.contains(vid)) {
      return genotype;
    }
  }
  return std::nullopt;
}

}  // namespace xoos::svc
