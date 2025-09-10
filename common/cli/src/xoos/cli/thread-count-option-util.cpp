#include "xoos/cli/thread-count-option-util.h"

#include <thread>

#include "xoos/util/parser.h"

namespace xoos::cli {

CLI::Option* AddThreadCountOption(AppPtr app, const std::string& name, size_t& value) {
  return app->add_option(name, value, "Number of threads used (0=available hardware threads)")
      ->transform([name](const std::string& inner_value) {
        try {
          const size_t requested_threads = util::ParseU64(inner_value);
          return std::to_string(requested_threads == 0 ? std::thread::hardware_concurrency() : requested_threads);
        } catch (const std::exception& e) {
          const auto msg = fmt::format("{}: {}", name, e.what());
          throw CLI::ConversionError(msg);
        }
      })
      ->default_val(1);
}

}  // namespace xoos::cli
