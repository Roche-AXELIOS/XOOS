#include "file-util.h"

#include <filesystem>

#include <xoos/error/error.h>

namespace xoos::svc {

fs::path GetBamIndexPath(const fs::path& bam_path) {
  // NOTE: CRAM support for feature extraction not yet implemented
  const std::string index_extension = (bam_path.extension() == ".cram") ? ".crai" : ".bai";
  auto index_path = fs::path{bam_path.string() + index_extension};
  if (fs::exists(index_path)) {
    return index_path;
  }
  auto index_fallback_path = fs::path{bam_path}.replace_extension(index_extension);
  if (fs::exists(index_fallback_path)) {
    return index_fallback_path;
  }
  throw error::Error("No {} index exists for input {}", index_extension, bam_path);
}

void CreateParentDirectoryIfNotExists(const fs::path& file_path) {
  const auto parent_path = file_path.parent_path();
  if (!parent_path.empty()) {
    fs::create_directories(parent_path);
  }
}

}  // namespace xoos::svc
