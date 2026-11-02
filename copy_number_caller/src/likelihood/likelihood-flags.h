#pragma once
#include <string>

#include <xoos/types/vec.h>

namespace xoos::cnc {
static const std::string kLikelihoodFlagNone = ".";
static const std::string kLikelihoodFlagPass = "PASS";
static const std::string kLikelihoodFlagLowAvgMeanMapq = "lowMmMapq";
static const std::string kLikelihoodFlagNonIntegerTcn = "nonIntegerTcn";
static const std::string kLikelihoodFlagCnvLength = "cnvLength";
static const std::string kLikelihoodFlagChrYFemale = "chrYFemale";
static const vec<std::string> kLikelihoodFlagsAll = {kLikelihoodFlagNone,
                                                     kLikelihoodFlagLowAvgMeanMapq,
                                                     kLikelihoodFlagNonIntegerTcn,
                                                     kLikelihoodFlagCnvLength,
                                                     kLikelihoodFlagChrYFemale};
vec<std::string> StringToFlags(const std::string& str);
std::string FlagsToString(const vec<std::string>& flags);
}  // namespace xoos::cnc
