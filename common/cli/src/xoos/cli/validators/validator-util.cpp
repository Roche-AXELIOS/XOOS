#include "xoos/cli/validators/validator-util.h"

#include <exception>
#include <filesystem>
#include <functional>
#include <string>
#include <system_error>

namespace xoos::cli {

namespace {

/**
 * @brief Report a missing path or a directory before the content check runs.
 *
 * Every content validator opens the path and would fail on these anyway, but through the underlying
 * library, which reports an errno string ("Failed to open HTS file ...: No such file or directory") and
 * may also write its own message to stderr. Checking here means each option reports the same clear
 * reason without having to chain `CLI::ExistingFile` itself, and the messages match CLI11's own so the
 * two are indistinguishable to the user.
 *
 * @param path Filesystem path supplied for the option.
 * @return Error message, or an empty string when the path is an existing non-directory.
 */
std::string CheckPathIsExistingFile(const fs::path& path) {
  std::error_code error_code;
  const auto status = fs::status(path, error_code);
  if (!error_code) {
    if (!fs::exists(status)) {
      return "File does not exist: " + path.string();
    }
    if (fs::is_directory(status)) {
      return "File is actually a directory: " + path.string();
    }
    return {};
  }
  if (error_code == std::errc::no_such_file_or_directory) {
    return "File does not exist: " + path.string();
  }
  return {};
}

}  // namespace

std::function<std::string(const std::string&)> CreateValidatorFunction(
    const std::function<void(const fs::path&)>& validate) {
  return [validate](const std::string& str) {
    const fs::path path(str);
    if (auto message = CheckPathIsExistingFile(path); !message.empty()) {
      return message;
    }
    try {
      validate(path);
    } catch (const std::exception& e) {  // NOSONAR cpp:S1181 — broad catch; converts exception to a validator string
      // convert the exception message to a string and return it
      // a non-empty return value indicates validation failure
      return std::string(e.what());
    }
    return std::string();
  };
}

}  // namespace xoos::cli
