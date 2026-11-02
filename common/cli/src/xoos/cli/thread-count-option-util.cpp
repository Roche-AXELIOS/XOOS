#include "xoos/cli/thread-count-option-util.h"

#include <cstddef>
#include <exception>
#include <string>
#include <thread>

#include <fmt/format.h>

#include "CLI/CLI.hpp"
#include "xoos/util/parse-int.h"

namespace xoos::cli {

CLI::Option* AddThreadCountOption(CLI::App* const app, const std::string& name, size_t& value) {
  return app->add_option(name, value, "Number of threads used (0=available hardware threads)")
      ->transform([name](const std::string& inner_value) {
        try {
          const size_t requested_threads = util::ParseU64(inner_value);
          return std::to_string(requested_threads == 0 ? std::thread::hardware_concurrency() : requested_threads);
        } catch (const std::exception&) {
          throw CLI::ConversionError(
              fmt::format("{}: '{}' cannot be converted to option, must be a non-negative integer "
                          "(0=available hardware threads).",
                          name,
                          inner_value));
        }
      })
      ->default_val(1);
}

}  // namespace xoos::cli
