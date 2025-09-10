#include <unistd.h>

#include <filesystem>

#include <fmt/format.h>

#include <xoos/cli/validators/file-permission-validator.h>
#include <xoos/types/fs.h>

namespace xoos::cli {

std::optional<std::string> CheckFileReadable(const fs::path& path) {
  if (!fs::exists(path)) {
    return fmt::format("Could not read file '{}'; does not exist", path.string());
  }
  if (!fs::is_regular_file(path)) {
    return fmt::format("Could not read file '{}'; not a regular file", path.string());
  }
  if (access(path.c_str(), R_OK) != 0) {
    return fmt::format("Could not read file '{}'; not readable", path.string());
  }
  return std::nullopt;
}

CLI::Validator FileReadableValidator() {
  const auto func = [](const std::string& input) -> std::string {
    const auto error_message = CheckFileReadable(input);
    return error_message ? *error_message : std::string{};
  };
  return CLI::Validator{func, "FILE_READABLE"};
}

std::optional<std::string> CheckFileWriteable(const fs::path& path) {
  fs::path directory = path.empty() ? fs::current_path() : fs::absolute(path).parent_path();
  if (!path.empty() && !path.parent_path().empty()) {
    while (!directory.empty() && !fs::exists(directory)) {
      directory = directory.parent_path();
    }
  }
  if (access(directory.c_str(), W_OK) != 0) {
    return fmt::format(
        "Cannot write to file '{}'; parent directory '{}' is not writeable", path.string(), directory.string());
  }
  return std::nullopt;
}

CLI::Validator FileWriteableValidator() {
  const auto func = [](const std::string& input) -> std::string {
    const auto error_message = CheckFileWriteable(input);
    return error_message ? *error_message : std::string{};
  };
  return CLI::Validator{func, "FILE_WRITABLE"};
}

std::variant<std::string, vec<fs::path>> CheckFileListReadable(const fs::path& path) {
  const auto error_message = CheckFileReadable(path);
  if (error_message) {
    return *error_message;
  }

  std::ifstream ifs(path);
  if (!ifs.is_open()) {
    return fmt::format("Could not open file list '{}'", path.string());
  }

  vec<fs::path> files;
  std::string line;
  while (std::getline(ifs, line)) {
    const auto error_message = CheckFileReadable(line);
    if (error_message) {
      return fmt::format("Cannot read file from list '{}'; {}", path.string(), *error_message);
    }
    files.emplace_back(line);
  }
  return files;
}

CLI::Validator FileListReadableValidator() {
  const auto func = [](std::string& input) -> std::string {
    const auto result = CheckFileListReadable(input);
    if (std::holds_alternative<std::string>(result)) {
      return std::get<std::string>(result);
    }
    return {};
  };
  return CLI::Validator({func, "FILE_LIST_READABLE"});
}

}  // namespace xoos::cli
