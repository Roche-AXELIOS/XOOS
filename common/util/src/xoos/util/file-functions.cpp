#include "file-functions.h"

#include <unistd.h>

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace xoos::file {

namespace fs = std::filesystem;

bool FileExists(const fs::path& p, fs::file_status s) {  // NOSONAR - exported helper used by downstream modules/tests
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

void CheckOutputFilePathsForWriting(const std::vector<fs::path>& output_paths, const bool overwrite) {
  // Output filenames may live in subfolders. Create each parent directory before the
  // existence/writability checks so callers do not fail only because a subdirectory is missing.
  for (const auto& path : output_paths) {
    const fs::path parent = path.parent_path();
    if (!parent.empty()) {
      CreateWritableDirectory(parent);
    }
  }
  CheckFilePermissions({}, output_paths);
  for (const auto& path : output_paths) {
    if (!FileExists(path)) {
      continue;
    }
    if (!overwrite) {
      throw std::runtime_error("Output file already exists: " + path.string());
    }
    if (!fs::is_regular_file(path)) {
      throw std::runtime_error("Output path is not a regular file: " + path.string());
    }
    CheckFileIsWritable(path);
  }
}

void CheckFilePermissionsAndOutputPathExistence(const std::vector<fs::path>& input_paths,
                                                const std::vector<fs::path>& output_paths) {
  for (const auto& path : input_paths) {
    CheckFileIsReadable(path);
  }
  CheckOutputFilePathsForWriting(output_paths, false);
}

void CreateWritableDirectory(const fs::path& output_directory) {
  // if the  output directory already exists, we check if it is a directory and fail if not
  if (std::filesystem::exists(output_directory) && !std::filesystem::is_directory(output_directory)) {
    throw std::runtime_error("Output directory exists and is a file: " + output_directory.string());
  }

  // if it is a directory, we create it directly
  if (!FileExists(output_directory)) {
    create_directories(output_directory);
  }
  CheckFileIsWritable(output_directory);
}

fs::path PrefixedOutputPath(const fs::path& output_dir,  // NOSONAR - exported helper used by downstream modules/tests
                            const std::string& output_prefix,
                            const fs::path& default_filename) {
  if (output_prefix.empty()) {
    return output_dir / default_filename;
  }
  // The default filename may embed a relative subfolder (e.g. "metrics/metrics.tsv"). Apply the
  // prefix to the basename so the subfolder is preserved: "metrics/metrics.tsv" with prefix
  // "SAMPLE" -> "metrics/SAMPLE_metrics.tsv". Filenames without a subfolder are unaffected.
  const fs::path prefixed_name = output_prefix + "_" + default_filename.filename().string();
  return output_dir / default_filename.parent_path() / prefixed_name;
}

}  // namespace xoos::file
