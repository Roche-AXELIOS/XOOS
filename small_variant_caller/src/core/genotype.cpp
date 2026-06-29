#include "genotype.h"

#include <xoos/types/str-container.h>

namespace xoos::svc {

using enum Genotype;

/**
 * @brief Convert genotype to VCF FORMAT field "GT" string.
 * @param gt Genotype
 * @return String representation
 */
std::string GenotypeToString(const Genotype gt) {
  switch (gt) {
    case kGT00:
      return "0/0";
    case kGT01:
      return "0/1";
    case kGT11:
      return "1/1";
    case kGT12:
      return "1/2";
    case kGT0:
      return "0";
    case kGT1:
      return "1";
    case kGTN:
      return ".";
    default:
      return "./.";
  }
}

/**
 * @brief Convert genotype to integer class.
 * @param gt Genotype
 * @return Integer class
 */
u32 GenotypeToInt(const Genotype gt) {
  // Only 4 classes are used for germline model training.
  switch (gt) {
    case kGT01:
      return 1;
    case kGT11:
    case kGT1:
      return 2;
    case kGT12:
      return 3;
    default:
      // unsupported genotypes: GT=0/0, GT=0, GT=./., GT=.
      return 0;
  }
}

/**
 * @brief Convert genotype to integer class label.
 * @param gt Genotype
 * @return Integer class label
 */
std::string GenotypeToIntString(const Genotype gt) {
  // Only 4 classes are used for germline model training.
  switch (gt) {
    case kGT01:
      return "1";
    case kGT11:
    case kGT1:
      return "2";
    case kGT12:
      return "3";
    default:
      // unsupported genotypes: GT=0/0, GT=0, GT=./., GT=.
      return "0";
  }
}

/**
 * @brief Convert integer class to diploid genotype.
 * @param gt Integer class
 * @return Genotype
 */
Genotype IntToGenotype(const u64 gt) {
  static constexpr u64 kGt00Int = 0;
  static constexpr u64 kGt01Int = 1;
  static constexpr u64 kGt11Int = 2;
  static constexpr u64 kGt12Int = 3;
  switch (gt) {
    case kGt00Int:
      return kGT00;
    case kGt01Int:
      return kGT01;
    case kGt11Int:
      return kGT11;
    case kGt12Int:
      return kGT12;
    default:
      return kGTNA;
  }
}

/**
 * @brief Map of VCF FORMAT field "GT" string to Genotype.
 */
static const StrMap<Genotype> kStringToGenotype{{"0/0", kGT00},
                                                {"0/1", kGT01},
                                                {"1/0", kGT01},
                                                {"1/1", kGT11},
                                                {"1/2", kGT12},
                                                {"2/1", kGT12},
                                                {"0|0", kGT00},
                                                {"0|1", kGT01},
                                                {"1|0", kGT01},
                                                {"1|1", kGT11},
                                                {"1|2", kGT12},
                                                {"2|1", kGT12},
                                                {"0", kGT0},
                                                {"1", kGT1},
                                                {"./.", kGTNA},
                                                {".", kGTN}};

/**
 * @brief Convert VCF FORMAT field "GT" string to genotype.
 * Phasing symbols "/" and "|" are treated equivalently.
 * Both diploid genotypes (e.g. GT=0/1) and haploid genotypes (e.g. GT=1) are supported.
 * Unsupported or unrecognized genotype strings will be mapped to kGTNA.
 * @param gt FORMAT field "GT" string
 * @return Genotype
 */
Genotype StringToGenotype(const std::string& gt) {
  const auto itr = kStringToGenotype.find(gt);
  if (itr != kStringToGenotype.end()) {
    return itr->second;
  }
  // Anything else is an unsupported genotype
  return kGTNA;
}

}  // namespace xoos::svc
