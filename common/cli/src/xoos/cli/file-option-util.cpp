#include "file-option-util.h"

namespace xoos::cli {

CLI::Option* AddOutputFileOption(AppPtr app,
                                 const std::string& name,
                                 fs::path& value,
                                 const std::string& description,
                                 const std::string& default_value,
                                 const fs::path& out_dir) {
  return app->add_option(name, value, description)
      ->default_val(default_value)
      ->force_callback()
      ->transform([&out_dir](const std::string& path) -> std::string { return (out_dir / fs::path(path)).string(); })
      ->check(CLI::NonexistentPath);
}

CLI::Option* AddOutputFileOption(AppPtr app,
                                 const std::string& name,
                                 std::optional<fs::path>& value,
                                 const std::string& description,
                                 const fs::path& out_dir) {
  return app->add_option(name, value, description)
      ->transform([&out_dir](const std::string& path) -> std::string { return (out_dir / fs::path(path)).string(); })
      ->check(CLI::NonexistentPath);
}

}  // namespace xoos::cli
