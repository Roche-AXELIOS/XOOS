#pragma once

#include <string>

#include <xoos/types/fs.h>

namespace xoos::svc {

struct FeaturesToVcfParam {
  fs::path input_file;
  fs::path output_file;
  fs::path genome;
  std::string sample{"SAMPLE"};
};

void ConvertFeaturesToVcf(const FeaturesToVcfParam& param);

}  // namespace xoos::svc
