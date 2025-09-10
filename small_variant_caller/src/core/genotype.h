#pragma once

#include <string>

#include <xoos/types/int.h>

namespace xoos::svc {

// Genotype representation for VCF FORMAT field "GT"
enum class Genotype {
  kGT00,  // diploid, FAIL or false positive ALT allele, GT=0/0
  kGT01,  // diploid, heterozygous REF and ALT alleles, GT=0/1, 1/0, 0|1 or 1|0
  kGT11,  // diploid, homozygous ALT allele, GT=1/1 or 1|1
  kGT12,  // diploid, heterozygous ALT alleles, GT=1/2, 2/1, 1|2, or 2|1
  kGTNA,  // Unsupported genotype
  kGT0,   // haploid, FAIL or false positive ALT allele, GT=0
  kGT1    // haploid, ALT allele, GT=1
};

std::string GenotypeToString(Genotype gt);
u32 GenotypeToInt(Genotype gt);
std::string GenotypeToIntString(Genotype gt);
Genotype IntToGenotype(u64 gt);
Genotype StringToGenotype(const std::string& gt);

}  // namespace xoos::svc
