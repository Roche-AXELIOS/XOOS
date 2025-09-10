#pragma once

#include <string>

#include "cli.h"

namespace xoos::cli {

/// Add a standard thread count option to the CLI which defaults to 1, allows 0 to mean all available threads
CLI::Option* AddThreadCountOption(AppPtr app, const std::string& name, size_t& value);

}  // namespace xoos::cli
