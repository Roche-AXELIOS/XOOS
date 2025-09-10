#pragma once

#include <filesystem>
#include <optional>
#include <string>

#include "cli.h"

namespace xoos::cli {

namespace fs = std::filesystem;
/**
 * Add an output file option to the CLI which will be transformed to be relative to the output directory.
 * @param app The CLI app to add the option to
 * @param name The name of the option
 * @param value The value to store the option in
 * @param description The description of the option
 * @param default_value The default value of the option
 * @param out_dir The output directory to make the path relative to
 */
CLI::Option* AddOutputFileOption(AppPtr app,
                                 const std::string& name,
                                 fs::path& value,
                                 const std::string& description,
                                 const std::string& default_value,
                                 const fs::path& out_dir);

/**
 * Add an optional output file option to the CLI which will be transformed to be relative to the output directory.
 * @param app The CLI app to add the option to
 * @param name The name of the option
 * @param value The value to store the option in
 * @param description The description of the option
 * @param out_dir The output directory to make the path relative to
 */
CLI::Option* AddOutputFileOption(AppPtr app,
                                 const std::string& name,
                                 std::optional<fs::path>& value,
                                 const std::string& description,
                                 const fs::path& out_dir);

}  // namespace xoos::cli
