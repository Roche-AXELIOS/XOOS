#pragma once
#include <string>

#include "xoos/cli/cli.h"

namespace xoos::svc {

/**
 * @brief Struct to hold command line information: program name, version, and full command line.
 */
struct CommandLineInfo {
  std::string name;
  std::string version;
  std::string command_line;
};

CommandLineInfo GetCommandLineInfo(cli::ConstAppPtr app);
std::string GetVcfHeaderLine(const CommandLineInfo& info);

}  // namespace xoos::svc
