#include "command-line-info.h"

namespace xoos::svc {

/**
 * Extracts command line information from a CLI application instance.
 * @param app Pointer to the CLI application instance.
 * @return CommandLineInfo struct containing the program name, version, and full command line.
 */
CommandLineInfo GetCommandLineInfo(cli::ConstAppPtr app) {
  return CommandLineInfo{
      .name = app->get_name(), .version = app->version(), .command_line = cli::RenderCli(app, app->get_name())};
}

/**
 * Formats a CommandLineInfo struct into a VCF header line.
 * @param info CommandLineInfo struct containing the program name, version, and full command line.
 */
std::string GetVcfHeaderLine(const CommandLineInfo& info) {
  return fmt::format(
      "##RocheCommandLine=<ID={},Version={},CommandLine={}>", info.name, info.version, info.command_line);
}

}  // namespace xoos::svc
