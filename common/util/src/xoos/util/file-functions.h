#pragma once

#include <filesystem>
#include <vector>

namespace xoos::file {

namespace fs = std::filesystem;

bool FileExists(const fs::path& p, fs::file_status s = fs::file_status{});

void CheckFileIsReadable(const fs::path& p);

void CheckFileIsWritable(const fs::path& p);

void CheckFilePermissions(const std::vector<fs::path>& input_paths, const std::vector<fs::path>& output_paths);

void CheckFilePermissionsAndOutputPathExistence(const std::vector<fs::path>& input_paths,
                                                const std::vector<fs::path>& output_paths);

// Creates an output directory if it doesn't exist and validates writing
// If it does exist it ensures that it is a directory and writable
void CreateWritableDirectory(const fs::path& output_directory);

}  // namespace xoos::file
