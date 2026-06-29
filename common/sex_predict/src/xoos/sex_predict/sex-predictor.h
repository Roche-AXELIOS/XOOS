#pragma once

#include <filesystem>
#include <optional>
#include <vector>

#include <xoos/io/htslib-util/htslib-ptr.h>

#include "sex.h"

namespace xoos::sex_predict {

namespace fs = std::filesystem;

constexpr float kMinRatioForSexDet = 25;
constexpr float kMinRatioForSexNA = 20;

struct BamInfo {
  io::HtsFilePtr bam_file;
  io::HtsIdxPtr idx;
  io::BamHdrPtr header;
};

struct RegionInfo {
  int tid;
  std::string name;
};

BamInfo LoadBamInfo(const fs::path& bam_input, const fs::path& bam_input_index);
std::optional<RegionInfo> GetXChromosome(const io::BamHdrPtr& bam_header);
std::optional<RegionInfo> GetYChromosome(const io::BamHdrPtr& bam_header);
std::optional<RegionInfo> GetChromosomeForStrings(const io::BamHdrPtr& bam_header,
                                                  const std::vector<std::string>& chromosomes);
Sex PredictSex(const fs::path& bam_file_path);

}  // namespace xoos::sex_predict
