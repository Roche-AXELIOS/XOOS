#include <filesystem>
#include <functional>
#include <string>

#include <fmt/format.h>

#include <xoos/cli/fasta-option-util.h>

#include "CLI/CLI.hpp"

namespace xoos::cli {

namespace {

std::function<void(const fs::path&)> FastaOptionFunction(fs::path& fasta_input) {
  return [&fasta_input](const fs::path& path) {
    fasta_input = path;
    if (!fs::exists(fasta_input)) {
      throw CLI::ValidationError(fmt::format("File '{}' does not exist", fasta_input.string()));
    }
    const auto fasta_index = fs::path(fasta_input.string() + ".fai");
    if (!fs::exists(fasta_index)) {
      throw CLI::ValidationError(fmt::format("FASTA index '{}' does not exist and is required", fasta_index.string()));
    }
  };
}

}  // namespace

CLI::Option* AddFastaFileOption(CLI::App* const app,
                                const std::string& name,
                                fs::path& fasta_input,
                                const std::string& description) {
  return app->add_option_function<fs::path>(name, FastaOptionFunction(fasta_input), description);
}

}  // namespace xoos::cli
