#pragma once

#include <memory>
#include <optional>
#include <string>

#include <fmt/format.h>

#include <CLI/CLI.hpp>

namespace xoos::cli {
using AppPtr = CLI::App*;
using ConstAppPtr = const CLI::App*;

std::string RenderCli(ConstAppPtr app, const std::string& program_name);

std::string FullProgramName(const std::string& program_name, const std::string& version);

// Setup a default CLI application with standard options and version information.
// The comments_callback, if provided, can be called with the version and CLI arguments.
std::shared_ptr<CLI::App> SetupDefaultCli(const std::string& program_name, const std::string& version);

int RunCli(AppPtr app, int argc, const char* const* argv);

template <class CliOpts>
using DefineOptions = std::function<void(AppPtr, std::shared_ptr<CliOpts>&)>;

template <class CliOpts>
using Main = std::function<void(const CliOpts&)>;

template <class CliOpts>
using PreCallback = std::function<void(ConstAppPtr, const std::shared_ptr<CliOpts>&)>;

/**
 * Adds a fallback command to the CLI application.
 * The callback command is executed if no subcommands are specified.
 *
 * @tparam CliOpts The type of the CLI options object.
 * @param app A shared pointer to the CLI::App instance.
 * @param define_options A function that configures the CLI options. It takes a shared pointer to the CLI::App and a
 * shared pointer to CliOpts.
 * @param cli_opts A shared pointer to the CliOpts instance, which holds the parsed command line options.
 * @param main A function that contains the main logic to be executed if no subcommands are specified.
 * @param pre_callback A function that is executed before the main function. Any extra validation or any preprocess
 * on the options which can not be done using the CLI library, can be done here.
 */
template <class CliOpts>
void AddFallbackCommand(AppPtr app,
                        const DefineOptions<CliOpts>& define_options,
                        std::shared_ptr<CliOpts>& cli_opts,
                        const Main<CliOpts>& main,
                        const PreCallback<CliOpts>& pre_callback = nullptr) {
  define_options(app, cli_opts);
  app->callback([app, cli_opts, main, pre_callback] {
    if (pre_callback) {
      pre_callback(app, cli_opts);
    }
    if (app->get_subcommands().empty()) {
      main(*cli_opts);
    }
  });
}

/**
 * Adds a subcommand to the CLI application with specified options and a callback function.
 *
 * @tparam CliOpts The type of the CLI options object.
 * @param app A shared pointer to the CLI::App instance to which the subcommand will be added.
 * @param subcommand_name A string representing the name of the subcommand.
 * @param define_sub_command_options A function that configures the options for the subcommand. It takes a shared
 * pointer to the subcommand's CLI::App and a shared pointer to CliOpts.
 * @param cli_opts A shared pointer to the CliOpts instance, which holds the parsed command line options for the
 * subcommand.
 * @param subcommand_main A function that contains the main logic to be executed when the subcommand is invoked.
 * @param subcommand_description Subcommand description, which will be displayed in the help message.
 * @param pre_callback A function that is executed before the main function. Any extra validation or any preprocess
 * on the options which can not be done using the CLI library, can be done here.
 */
template <class CliOpts>
void AddSubcommand(AppPtr app,
                   const std::string& subcommand_name,
                   const DefineOptions<CliOpts>& define_sub_command_options,
                   std::shared_ptr<CliOpts>& cli_opts,
                   const Main<CliOpts>& subcommand_main,
                   const std::optional<std::string>& subcommand_description = std::nullopt,
                   const PreCallback<CliOpts>& pre_callback = nullptr) {
  auto* const subcommand = app->add_subcommand(subcommand_name);
  if (subcommand_description.has_value()) {
    subcommand->description(subcommand_description.value());
  }
  define_sub_command_options(subcommand, cli_opts);
  subcommand->callback([subcommand, subcommand_main, cli_opts, pre_callback] {
    if (pre_callback) {
      pre_callback(subcommand, cli_opts);  // Execute pre_call_back if provided
    }
    subcommand_main(*cli_opts);
  });
}

template <class CliOpts>
struct StandardMainParam {
  std::string program_name{};
  std::string version{};
  std::shared_ptr<CliOpts> cli_opts{};
  DefineOptions<CliOpts> define_options{};
  Main<CliOpts> main{};
  PreCallback<CliOpts> pre_callback{};
};

template <class CliOpts>
int StandardMain(int argc, const char* const* argv, StandardMainParam<CliOpts>& param) {
  const auto app = SetupDefaultCli(param.program_name, param.version);
  AddFallbackCommand(app.get(), param.define_options, param.cli_opts, param.main, param.pre_callback);
  return RunCli(app.get(), argc, argv);
}

template <typename... Args>
CLI::ValidationError ValidationError(fmt::format_string<Args...> fmt, Args&&... args) {
  return CLI::ValidationError(fmt::format(fmt, std::forward<Args>(args)...));
}
}  // namespace xoos::cli
