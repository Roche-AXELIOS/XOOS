#pragma once

#include <filesystem>
#include <string>

#include <CLI/CLI.hpp>

namespace xoos::cli {

namespace fs = std::filesystem;

/**
 * Adds an option to the CLI application for specifying a FASTA file. This
 * will also validate the FASTA index file.
 *
 * @param app The CLI application to which the option will be added.
 * @param name The name of the option.
 * @param fasta_input Reference to a `fs::path` where the FASTA file path will be stored.
 * @param description A description of the option.
 * @return A pointer to the created `CLI::Option` object.
 */
CLI::Option* AddFastaFileOption(CLI::App* app,
                                const std::string& name,
                                fs::path& fasta_input,
                                const std::string& description);
}  // namespace xoos::cli
