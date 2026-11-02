#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace xoos::file {

namespace fs = std::filesystem;

bool FileExists(const fs::path& p, fs::file_status s = fs::file_status{});

void CheckFileIsReadable(const fs::path& p);

void CheckFileIsWritable(const fs::path& p);

void CheckFilePermissions(const std::vector<fs::path>& input_paths, const std::vector<fs::path>& output_paths);

// Validate output file path permissions and reject existing files unless overwrite is enabled.
void CheckOutputFilePathsForWriting(const std::vector<fs::path>& output_paths, bool overwrite);

void CheckFilePermissionsAndOutputPathExistence(const std::vector<fs::path>& input_paths,
                                                const std::vector<fs::path>& output_paths);

// Creates an output directory if it doesn't exist and validates writing
// If it does exist it ensures that it is a directory and writable
void CreateWritableDirectory(const fs::path& output_directory);

// Build a default output path, optionally namespaced by a caller-supplied prefix.
// When the prefix is empty the default filename is used as-is, preserving backward-compatible
// output names. Otherwise the prefix and filename are joined with an underscore,
// e.g. prefix "SAMPLE" + "coverage.bed" -> "SAMPLE_coverage.bed".
// If the default filename embeds a relative subfolder, the prefix is applied to the basename and
// the subfolder is preserved, e.g. prefix "SAMPLE" + "metrics/metrics.tsv" -> "metrics/SAMPLE_metrics.tsv".
// The result is always returned as a path under output_dir.
fs::path PrefixedOutputPath(const fs::path& output_dir,
                            const std::string& output_prefix,
                            const fs::path& default_filename);

}  // namespace xoos::file
