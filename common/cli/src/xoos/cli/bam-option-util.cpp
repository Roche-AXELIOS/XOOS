#include <xoos/cli/bam-option-util.h>

namespace xoos::cli {

static std::function<void(const fs::path&)> BamOptionFunction(fs::path& bam_input, fs::path& bam_index) {
  return [&bam_input, &bam_index](const fs::path& path) -> void {
    bam_input = path;
    if (!exists(bam_input)) {
      throw CLI::ValidationError(fmt::format("File '{}' does not exist", bam_input.string()));
    }
    bam_index = fs::path(bam_input.string() + ".bai");
    if (!exists(bam_index)) {
      throw CLI::ValidationError(fmt::format("BAM index '{}' does not exist and is required", bam_index.string()));
    }
  };
}

CLI::Option* AddBamFileOption(
    AppPtr app, const std::string& name, fs::path& bam_input, fs::path& bam_index, const std::string& description) {
  return app->add_option_function<fs::path>(name, BamOptionFunction(bam_input, bam_index), description);
}

}  // namespace xoos::cli
