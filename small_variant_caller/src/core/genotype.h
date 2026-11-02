#pragma once

#include <string>

#include <xoos/types/int.h>

namespace xoos::svc {

// Genotype representation for VCF FORMAT field "GT"
enum class Genotype : u8 {
  // diploid, FAIL or false positive ALT allele, GT=0/0
  kGT00,
  // diploid, heterozygous REF and ALT alleles, GT=0/1, 1/0, 0|1 or 1|0
  kGT01,
  // diploid, homozygous ALT allele, GT=1/1 or 1|1
  kGT11,
  // diploid, heterozygous ALT alleles, GT=1/2, 2/1, 1|2, or 2|1
  kGT12,
  // Unsupported diploid genotype
  kGTNA,
  // haploid, FAIL or false positive ALT allele, GT=0
  kGT0,
  // haploid, ALT allele, GT=1
  kGT1,
  // Unsupported haploid genotype
  kGTN
};

/**
 * @brief Convert genotype to VCF FORMAT field "GT" string.
 * @param gt Genotype
 * @return String representation
 */
std::string GenotypeToString(Genotype gt);

/**
 * @brief Convert genotype to integer class.
 * @param gt Genotype
 * @return Integer class
 */
u32 GenotypeToInt(Genotype gt);

/**
 * @brief Convert genotype to integer class label.
 * @param gt Genotype
 * @return Integer class label
 */
std::string GenotypeToIntString(Genotype gt);

/**
 * @brief Convert integer class to diploid genotype.
 * @param gt Integer class
 * @return Genotype
 */
Genotype IntToGenotype(u64 gt);

/**
 * @brief Convert VCF FORMAT field "GT" string to genotype.
 * Phasing symbols "/" and "|" are treated equivalently.
 * Both diploid genotypes (e.g. GT=0/1) and haploid genotypes (e.g. GT=1) are supported.
 * Unsupported or unrecognized genotype strings will be mapped to kGTNA.
 * @param gt FORMAT field "GT" string
 * @return Genotype
 */
Genotype StringToGenotype(const std::string& gt);

}  // namespace xoos::svc
