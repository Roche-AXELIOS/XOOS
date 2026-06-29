#include "sex-predictor.h"

#include <algorithm>
#include <filesystem>
#include <vector>

#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>

namespace xoos::sex_predict {

// get the most likely X chromosome string from the bam header
std::optional<RegionInfo> GetChromosomeForStrings(const io::BamHdrPtr& bam_header,
                                                  const std::vector<std::string>& chromosomes) {
  if (chromosomes.empty()) {
    return std::nullopt;
  }
  // iterate through all regions in the header and find the X chromosome
  for (int i = 0; i < bam_header->n_targets; i++) {
    // convert bam_header->target_name[i] to a lower-case string
    std::string target_name = bam_header->target_name[i];
    std::ranges::transform(target_name, target_name.begin(), ::tolower);
    if (std::ranges::find(chromosomes, target_name) != chromosomes.end()) {
      return RegionInfo{i, bam_header->target_name[i]};
    }
  }
  return std::nullopt;
}

std::optional<RegionInfo> GetXChromosome(const io::BamHdrPtr& bam_header) {
  return GetChromosomeForStrings(bam_header, std::vector<std::string>{"chrx", "x"});
}

std::optional<RegionInfo> GetYChromosome(const io::BamHdrPtr& bam_header) {
  return GetChromosomeForStrings(bam_header, std::vector<std::string>{"chry", "y"});
}

BamInfo LoadBamInfo(const fs::path& bam_input, const fs::path& bam_input_index) {
  auto bam_file = io::HtsFilePtr(hts_open(bam_input.c_str(), "r"));
  if (!bam_file) {
    throw error::Error("Could not open BAM file {}", bam_input.string());
  }
  auto header = io::BamHdrPtr(sam_hdr_read(bam_file.get()));
  if (!header) {
    throw error::Error("Could not read BAM header {}", bam_input.string());
  }

  auto idx = io::HtsIdxPtr(sam_index_load(bam_file.get(), bam_input_index.c_str()));
  if (!idx) {
    throw error::Error("could not load BAM index {}", bam_input_index);
  }
  return BamInfo{std::move(bam_file), std::move(idx), std::move(header)};
}

// Get all the reads from the bam file index and calculate the xy ration between X and Y.
// If the ratio X / Y > 25 then return female
// If the ratio X / Y > 20 then return unknown
// Else return male.
Sex PredictSex(const fs::path& bam_file_path) {
  const fs::path bam_file_index = bam_file_path.string() + ".bai";

  const auto [bam_file, idx, header] = LoadBamInfo(bam_file_path, bam_file_index);
  const auto alignment_reference_regions = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>(bam_init1(), bam_destroy1);

  // get the number of reads for X and Y

  // count X reads
  const auto region_x = GetXChromosome(header);
  u64 count_x = 0;
  if (region_x.has_value()) {
    const auto iter_x = io::HtsItrMultiPtr(sam_itr_querys(idx.get(), header.get(), region_x.value().name.c_str()));
    if (!iter_x) {
      throw error::Error("Could not parse region '{}'", region_x.value().name);
    }
    u64 unmapped_x = 0;
    // use hts_idx_get_stat to get the number of reads in the region
    hts_idx_get_stat(idx.get(), region_x.value().tid, &count_x, &unmapped_x);
  }

  // count Y reads
  const auto region_y = GetYChromosome(header);
  u64 count_y = 0;
  if (region_y.has_value()) {
    const auto iter_y = io::HtsItrMultiPtr(sam_itr_querys(idx.get(), header.get(), region_y.value().name.c_str()));
    if (!iter_y) {
      throw error::Error("Could not parse region '{}'", region_y.value().name.c_str());
    }

    u64 unmapped_y = 0;
    hts_idx_get_stat(idx.get(), region_y.value().tid, &count_y, &unmapped_y);
  }

  // and calculate the ratio between X and Y
  const double xy_ratio = static_cast<double>(count_x) / static_cast<double>(count_y);
  if (std::isnan(xy_ratio)) {
    return Sex::kUnknown;
  }
  if (xy_ratio > kMinRatioForSexDet) {
    return Sex::kFemale;
  }
  if (xy_ratio > kMinRatioForSexNA) {
    return Sex::kUnknown;
  }
  return Sex::kMale;
}

}  // namespace xoos::sex_predict
