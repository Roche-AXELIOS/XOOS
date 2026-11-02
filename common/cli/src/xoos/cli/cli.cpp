#include "xoos/cli/cli.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <exception>
#include <locale>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <fmt/format.h>

#include <xoos/enum/enum-util.h>
#include <xoos/log/logging.h>

#include "CLI/CLI.hpp"

namespace xoos::cli {

namespace {

/**
 * get_expected()  means the number of arguments expected for the option which is 0 for flag type.
 * @param opt - CLI11 option
 * @return true if the option is of boolean type, false otherwise
 */
bool IsKeyBoolType(const CLI::Option* const opt) {
  return opt->get_expected() == 0;
}

/**
 * if default value is "false", then the flag value will come as "true" if the flag is provided by the user.
 * Otherwise, the flag value will be 1 if the flag is provided by the user.
 * @param value - Value of the flag to be checked
 * @return true if the value is "true" or "1", false otherwise
 */
bool IsValueTrue(const std::string& value) {
  std::string lower_value = value;
  std::ranges::transform(
      lower_value, lower_value.begin(), [](const char c) { return std::tolower(c, std::locale::classic()); });
  return lower_value == "true" || lower_value == "1";
}

/**
 * Returns the CLI name of an option, handling hidden options (group == "").
 * CLI11's get_name() returns "" when the option group is empty (hidden), so we
 * fall back to lnames_[0] / snames_[0] directly via the public accessors.
 */
std::string GetOptionCliName(const CLI::Option* const opt) {
  auto name = opt->get_name();
  if (!name.empty()) {
    return name;
  }
  // Hidden option: get_name() returns "" because group is empty; use raw name lists
  if (!opt->get_lnames().empty()) {
    return "--" + opt->get_lnames()[0];
  }
  if (!opt->get_snames().empty()) {
    return "-" + opt->get_snames()[0];
  }
  return {};
}

std::string GetRenderedValue(const CLI::Option* const opt, const std::string& value) {
  const auto cli_name = GetOptionCliName(opt);
  if (cli_name.empty()) {
    return {};
  }
  if (!IsKeyBoolType(opt)) {
    return fmt::format(" {} {}", cli_name, value);
  }
  return IsValueTrue(value) ? fmt::format(" {}", cli_name) : "";
}

}  // namespace

/**
 * Renders the command-line string for a single application or subcommand level.
 * * This function iterates through all registered options of the provided app,
 * appending either the user-provided values or the default values (if available)
 * to the program name to reconstruct the command as it would appear in a shell.
 *
 * @param app A constant pointer to the CLI::App instance to render.
 * @param program_name The name of the command or subcommand to use as the base of the string.
 * @return A string representing the rendered command line for this specific level.
 */
std::string RenderCli(const CLI::App* const app, const std::string& program_name) {
  std::stringstream ss;
  ss << program_name;
  // Iterate over all options and print their names and values
  for (const auto* const opt : app->get_options()) {
    // Check if the option was actually provided on the command line
    if (opt->empty()) {
      const bool dependency_satisfied = std::ranges::all_of(
          opt->get_needs(), [](const CLI::Option* const dependency_opt) { return !dependency_opt->empty(); });
      // Option not provided by the user, print default value if available
      if (dependency_satisfied && !opt->get_default_str().empty()) {
        ss << GetRenderedValue(opt, opt->get_default_str());
      }
    } else {  // Print the user provided inputs
      for (const auto& value : opt->as<std::vector<std::string>>()) {
        ss << GetRenderedValue(opt, value);
      }
    }
  }
  return ss.str();
}

/**
 * Renders the complete command-line string, including parent commands if the app is a subcommand.
 * * This function checks for the existence of a parent application. If a parent exists,
 * it concatenates the rendered parent command with the current subcommand to provide
 * the full execution path. Otherwise, it renders the standalone application command.
 *
 * @param app A constant pointer to the CLI::App instance (either the main app or a subcommand).
 * @return A string containing the full reconstructed command line.
 */
std::string RenderFullCli(const CLI::App* const app) {
  std::string command_line;
  // If the app is a subcommand, we need to render the full command line including the parent command(s)
  if (const auto* const parent_app = app->get_parent(); parent_app != nullptr) {
    command_line = RenderCli(parent_app, parent_app->get_name()) + " " + RenderCli(app, app->get_name());
  } else {
    command_line = RenderCli(app, app->get_name());
  }
  return command_line;
}

std::string FullProgramName(const std::string& program_name, const std::string& version) {
  return program_name + " " + version;
}

std::shared_ptr<CLI::App> SetupDefaultCli(const std::string& program_name, const std::string& version) {
  Logging::Initialize();
  Logging::SetLevel(log::LogLevel::kInfo);

  auto app = std::make_shared<CLI::App>(FullProgramName(program_name, version), program_name);
  app->set_version_flag("-v,--version", version);

  const auto log_level_desc = fmt::format("Log level: {}", enum_util::FormatEnumNames<log::LogLevel>());
  const auto* const log_level_option =
      app->add_option<const std::string>("-l,--log-level", log_level_desc)->default_val("info");

  app->parse_complete_callback([app = app.get(), version, log_level_option]() {
    Logging::Info("Version: {}", version);
    Logging::Info("For Research Use Only. Not for use in diagnostic procedures.");

    // Render the full command line including any active subcommands.
    // Uses GetCommandLineInfo to share the rendering logic with VCF/TSV metadata output.
    const auto info = GetCommandLineInfo(app);
    Logging::Info("Running: {}", info.command_line);

    const auto log_level = log::ParseLogLevel(log_level_option->as<std::string>());
    if (!log_level) {
      throw std::runtime_error(fmt::format("Invalid log level: '{}'", log_level_option->as<std::string>()));
    }
    Logging::SetLevel(*log_level);
  });

  return app;
}

s32 RunCli(CLI::App* const app, const s32 argc, const char* const* const argv) {
  try {
    app->parse(argc, argv);
  } catch (const CLI::CallForHelp& e) {
    return app->exit(e);
  } catch (const CLI::CallForVersion& e) {
    return app->exit(e);
  } catch (const CLI::ParseError& e) {
    app->exit(CLI::CallForHelp());
    Logging::Error(e);
    return EXIT_FAILURE;
  } catch (const std::exception& e) {  // NOSONAR cpp:S1181 — top-level handler; consistent boundary for runtime errors
    Logging::Error(e);
    return EXIT_FAILURE;
  } catch (...) {  // NOSONAR cpp:S2738 — last-resort catch-all; clean exit instead of terminate
    Logging::Error("Unknown failure occurred");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/// Called from two contexts:
///  1. SetupDefaultCli's parse_complete_callback (app = root) — for the "Running: ..." log line.
///  2. Module pre-callbacks (app may be a subcommand) — for VCF/TSV metadata.
///
/// These run at different times: parse_complete_callback fires after parsing but before
/// subcommand callbacks, so the root app and all active subcommands are already resolved.
/// Module pre-callbacks fire later from within a subcommand callback, so app may point to
/// the subcommand rather than the root. The parent traversal handles both cases, always
/// producing the full command line with the root program name and version.
io::CommandLineInfo GetCommandLineInfo(const CLI::App* const app) {
  const auto* root = app;
  while (root->get_parent() != nullptr) {
    root = root->get_parent();
  }

  auto cli_args = RenderFullCli(app);
  // get_subcommands() with no arguments returns only parsed (active) subcommands, not all defined ones.
  for (const auto* const sub : app->get_subcommands()) {
    cli_args += " " + RenderCli(sub, sub->get_name());
  }

  return io::CommandLineInfo{.name = root->get_name(), .version = root->version(), .command_line = cli_args};
}

}  // namespace xoos::cli
