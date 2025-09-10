#include "xoos/cli/cli.h"

#include <xoos/enum/enum-util.h>
#include <xoos/util/uuid.h>

namespace xoos::cli {

/**
 * get_expected()  means the number of arguments expected for the option which is 0 for flag type.
 * @param opt - CLI11 option
 * @return true if the option is of boolean type, false otherwise
 */
static bool IsKeyBoolType(const CLI::Option* opt) {
  return opt->get_expected() == 0;
}

/**
 * if default value is "false", then the flag value will come as "true" if the flag is provided by the user.
 * Otherwise, the flag value will be 1 if the flag is provided by the user.
 * @param value - Value of the flag to be checked
 * @return true if the value is "true" or "1", false otherwise
 */
static bool IsValueTrue(const std::string& value) {
  std::string lower_value = value;
  std::transform(lower_value.begin(), lower_value.end(), lower_value.begin(), ::tolower);
  return lower_value == "true" || lower_value == "1";
}

static std::string GetRenderedValue(const CLI::Option* opt, const std::string& value) {
  if (!IsKeyBoolType(opt)) {
    return fmt::format(" {} {}", opt->get_name(), value);
  }
  return IsValueTrue(value) ? fmt::format(" {}", opt->get_name()) : "";
}

std::string RenderCli(ConstAppPtr app, const std::string& program_name) {
  std::stringstream ss;
  ss << program_name;
  // Iterate over all options and print their names and values
  for (const auto& opt : app->get_options()) {
    // Check if the option was actually provided on the command line
    if (opt->empty()) {
      // Option not provided by the user, print default value if available
      if (!opt->get_default_str().empty()) {
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

std::string FullProgramName(const std::string& program_name, const std::string& version) {
  return program_name + " " + version;
}

std::shared_ptr<CLI::App> SetupDefaultCli(const std::string& program_name, const std::string& version) {
  Logging::Initialize();
  Logging::SetLevel(log::LogLevel::kInfo);

  auto app = std::make_shared<CLI::App>(FullProgramName(program_name, version), program_name);
  app->set_version_flag("-v,--version", version);

  auto log_level_desc = fmt::format("Log level: {}", enum_util::FormatEnumNames<log::LogLevel>());
  auto* log_level_option = app->add_option<const std::string>("-l,--log-level", log_level_desc)->default_val("info");

  app->parse_complete_callback([app = app.get(), program_name, version, log_level_option]() {
    Logging::Info("Version: {}", version);
    Logging::Info("Execution ID: {}", uuid::ExecutionId().IsoDateTimeString());

    /// Render the command-line arguments for the main application. If the application has subcommands, render the
    /// command-line arguments for each subcommand. Few of the common options (like log-level etc.) will be part of main
    /// application. So when a subcommand will be triggered, then the program will render the subcommand-line arguments
    /// along with common options. For example, if the main application is `copy_number_caller` and it has a subcommand
    /// `GCCorrect`, the command-line arguments for `copy_number_caller` (mostly common options) and `GCCorrect` will be
    /// rendered.
    auto cli_args = RenderCli(app, program_name);
    for (const auto& subcommand : app->get_subcommands()) {
      auto sub_app = app->get_subcommand_ptr(subcommand);
      cli_args += " " + RenderCli(sub_app.get(), sub_app->get_name());
    }
    Logging::Info("Running: {}", cli_args);

    const auto log_level = log::ParseLogLevel(log_level_option->as<std::string>());
    if (!log_level) {
      throw std::runtime_error(fmt::format("Invalid log level: '{}'", log_level_option->as<std::string>()));
    }
    Logging::SetLevel(*log_level);
  });

  return app;
}

int RunCli(AppPtr app, int argc, const char* const* argv) {
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
  } catch (const std::exception& e) {
    Logging::Error(e);
    return EXIT_FAILURE;
  } catch (...) {
    Logging::Error("Unknown failure occurred");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

}  // namespace xoos::cli
