#pragma once

#include <xoos/io/bed-region.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/command-line-info.h"
#include "core/config.h"
#include "core/workflow.h"

namespace xoos::svc {

using BedRegion = io::BedRegion;

struct TrainModelParam {
  vec<fs::path> positive_features{};
  vec<fs::path> positive_vcf_features{};
  vec<fs::path> negative_features{};
  vec<fs::path> negative_vcf_features{};
  vec<fs::path> truth_vcfs{};
  vec<fs::path> output_file{};
  SVCConfig config{};
  vec<BedRegion> blocklist{};
  Workflow workflow{};
  u32 max_score{};
  u32 iterations{};
  u32 snv_iterations{};
  u32 indel_iterations{};
  size_t threads{};
  bool normalize_features{};
  bool write_training_data_tsv{};
  std::optional<CommandLineInfo> command_line;
};

void TrainModel(const TrainModelParam& param);

}  // namespace xoos::svc
