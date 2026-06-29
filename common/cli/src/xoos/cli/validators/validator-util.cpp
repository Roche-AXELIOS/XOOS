#include "xoos/cli/validators/validator-util.h"

namespace xoos::cli {

std::function<std::string(const std::string&)> CreateValidatorFunction(
    const std::function<void(const fs::path&)>& validate) {
  return [validate](const std::string& str) {
    try {
      validate(fs::path(str));
    } catch (const std::exception& e) {
      // convert the exception message to a string and return it
      // a non-empty return value indicates validation failure
      return std::string(e.what());
    }
    return std::string();
  };
}

}  // namespace xoos::cli
