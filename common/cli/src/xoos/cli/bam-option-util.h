#pragma once

#include <filesystem>
#include <string>

#include "cli.h"

namespace xoos::cli {

namespace fs = std::filesystem;
/**
 * Adds an option to the CLI application for specifying a BAM file. This
 * will also validate the BAM index file.
 *
 * @param app The CLI application to which the option will be added.
 * @param name The name of the option.
 * @param bam_input Reference to a `fs::path` where the BAM file path will be stored.
 * @param bam_index Reference to a `fs::path` where the BAM index file path will be stored.
 * @param description A description of the option.
 * @return A pointer to the created `CLI::Option` object.
 */
CLI::Option* AddBamFileOption(
    AppPtr app, const std::string& name, fs::path& bam_input, fs::path& bam_index, const std::string& description);
}  // namespace xoos::cli
