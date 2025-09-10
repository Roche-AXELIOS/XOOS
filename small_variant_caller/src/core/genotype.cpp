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
    default:
      return ".";
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
    case kGTNA:
    case kGT00:
    case kGT0:
      return 0;
    case kGT01:
      return 1;
    case kGT11:
    case kGT1:
      return 2;
    case kGT12:
      return 3;
    default:
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
    case kGTNA:
    case kGT00:
    case kGT0:
      return "0";
    case kGT01:
      return "1";
    case kGT11:
    case kGT1:
      return "2";
    case kGT12:
      return "3";
    default:
      return "0";
  }
}

/**
 * @brief Convert integer class to diploid genotype.
 * @param gt Integer class
 * @return Genotype
 */
Genotype IntToGenotype(const u64 gt) {
  switch (gt) {
    case 0:
      return kGT00;
    case 1:
      return kGT01;
    case 2:
      return kGT11;
    case 3:
      return kGT12;
    default:
      return kGTNA;
  }
}

/**
 * @brief Convert VCF FORMAT field "GT" string to genotype.
 * @param gt FORMAT field "GT" string
 * @return Genotype
 */
Genotype StringToGenotype(const std::string& gt) {
  // Phasing is ignored.
  // Ambiguous GT string (e.g. ".", "./.", "./1", "1/.") are not supported.
  static StrMap<Genotype> gt_str_to_int{{"0/0", kGT00},
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
                                        {"1", kGT1}};
  auto itr = gt_str_to_int.find(gt);
  if (itr != gt_str_to_int.end()) {
    return itr->second;
  }
  return kGTNA;  // unsupported genotype
}

}  // namespace xoos::svc
