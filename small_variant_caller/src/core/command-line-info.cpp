#include "command-line-info.h"

namespace xoos::svc {

/**
 * Extracts command line information from a CLI application instance.
 * @param app Pointer to the CLI application instance.
 * @return CommandLineInfo struct containing the program name, version, and full command line.
 */
CommandLineInfo GetCommandLineInfo(cli::ConstAppPtr app) {
  // If the application has subcommands, render the command-line arguments for each subcommand.
  // Common options (like log-level etc.) will be part of main application.
  // When a subcommand is triggered, then the subcommand-line arguments will be rendered along with common options.
  auto cli_args = cli::RenderCli(app, app->get_name());
  for (const auto* const sub : app->get_subcommands()) {
    cli_args += " " + cli::RenderCli(sub, sub->get_name());
  }
  return CommandLineInfo{.name = app->get_name(), .version = app->version(), .command_line = cli_args};
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
