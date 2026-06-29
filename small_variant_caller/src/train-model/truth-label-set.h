#pragma once

#include <filesystem>
#include <optional>
#include <unordered_set>

#include <xoos/types/fs.h>

#include "core/genotype.h"
#include "core/variant-id.h"

namespace xoos::svc {

/**
 * @brief Map of genotypes to their corresponding sets of VariantIds for a sample.
 * @note The key is the genotype label, and the value is a set of VariantIds that have that genotype label in the
 * sample.
 * @note Since a variant can only have one genotype label in a sample, the sets of VariantIds for different genotype
 * labels should be mutually exclusive.
 * @note This structure is more memory efficient than a map of VariantId to genotype label, which requires storing the
 * genotype label for each variant.
 */
using GenotypeToVariantIds = std::unordered_map<Genotype, std::unordered_set<VariantId>>;

/**
 * @brief Helper function to get the genotype label for a variant from the map of truth genotypes to variant IDs.
 * @param vid VariantId of the variant to get the genotype label for.
 * @param genotypes Map of ground truth genotypes to their corresponding sets of VariantIds.
 * @return Genotype label for the variant if found in the truth genotypes, std::nullopt otherwise.
 */
std::optional<Genotype> GetGenotypeForVariant(const VariantId& vid, const GenotypeToVariantIds& genotypes);

/**
 * @brief Helper function to extract the genotypes for ground truth variants and organize them in a map where the key
 * is the GT label and the value is a set of VariantIds with that GT label.
 * @param vcf VCF file containing ground truth variants and their genotypes.
 * @return Map of ground truth genotypes to their corresponding sets of VariantIds.
 */
GenotypeToVariantIds GetGenotypeToVariantIds(const fs::path& vcf);

}  // namespace xoos::svc
