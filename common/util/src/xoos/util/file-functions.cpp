#include "file-functions.h"

#include <unistd.h>

namespace xoos::file {

namespace fs = std::filesystem;

bool FileExists(const fs::path& p, fs::file_status s) {
  return fs::status_known(s) ? fs::exists(s) : fs::exists(p);
}

void CheckFileIsReadable(const fs::path& p) {
  if (access(p.c_str(), R_OK) != 0) {
    throw std::runtime_error("Path is not readable: " + p.string());
  }
}

void CheckFileIsWritable(const fs::path& p) {
  if (access(p.c_str(), W_OK) != 0) {
    throw std::runtime_error("Path is not writable: " + p.string());
  }
}

void CheckFilePermissions(const std::vector<fs::path>& input_paths, const std::vector<fs::path>& output_paths) {
  for (const auto& path : input_paths) {
    CheckFileIsReadable(path);
  }

  for (const auto& path : output_paths) {
    if (path.empty()) {
      throw std::runtime_error("Required output file not defined");
    }
    const fs::path directory = path.parent_path().empty() ? fs::current_path() : path.parent_path();
    CheckFileIsWritable(directory);
  }
}

void CheckFilePermissionsAndOutputPathExistence(const std::vector<fs::path>& input_paths,
                                                const std::vector<fs::path>& output_paths) {
  for (const auto& path : input_paths) {
    CheckFileIsReadable(path);
  }

  for (const auto& path : output_paths) {
    if (path.empty()) {
      throw std::runtime_error("Required output file not defined");
    }
    if (FileExists(path)) {
      throw std::runtime_error("Output file already exists: " + path.string());
    }
    const fs::path directory = path.parent_path().empty() ? fs::current_path() : path.parent_path();
    CheckFileIsWritable(directory);
  }
}

void CreateWritableDirectory(const fs::path& output_directory) {
  // if the  output directory already exists, we check if it is a directory and fail if not
  if (std::filesystem::exists(output_directory) && !std::filesystem::is_directory(output_directory)) {
    throw std::runtime_error("Output directory exists and is a file: " + output_directory.string());
  }

  // if it is a directory, we create it directly
  if (!FileExists(output_directory)) {
    create_directory(output_directory);
  }
  CheckFileIsWritable(output_directory);
}

}  // namespace xoos::file
