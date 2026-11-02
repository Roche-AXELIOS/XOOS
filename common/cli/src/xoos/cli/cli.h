#pragma once

#include <memory>
#include <optional>
#include <string>

#include <fmt/format.h>

#include <CLI/CLI.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/types/int.h>

namespace xoos::cli {
std::string RenderCli(const CLI::App* app, const std::string& program_name);

std::string RenderFullCli(const CLI::App* app);

std::string FullProgramName(const std::string& program_name, const std::string& version);

/**
 * @brief Extract program name, version, and rendered command line from a CLI app.
 *
 * Traverses to the root app for name/version. Renders the full command line
 * including parent options (via RenderFullCli) and any active subcommands.
 */
io::CommandLineInfo GetCommandLineInfo(const CLI::App* app);  // NOSONAR — used by downstream modules

// Setup a default CLI application with standard options and version information.
// The comments_callback, if provided, can be called with the version and CLI arguments.
std::shared_ptr<CLI::App> SetupDefaultCli(const std::string& program_name, const std::string& version);

s32 RunCli(CLI::App* app, s32 argc, const char* const* argv);

template <class CliOpts>
using DefineOptions = std::function<void(CLI::App*, std::shared_ptr<CliOpts>&)>;

template <class CliOpts>
using Main = std::function<void(const CliOpts&)>;

// The app is passed as non-const so a pre-callback can finalize option state before the command line
// is rendered for provenance metadata (e.g. Option::default_str for a default resolved at parse time).
// Callbacks that only read the app may still take a `const CLI::App*`.
template <class CliOpts>
using PreCallback = std::function<void(CLI::App*, const std::shared_ptr<CliOpts>&)>;

/**
 * Adds a fallback command to the CLI application.
 *
 * @tparam CliOpts The type of the CLI options object.
 * @param app A non-owning raw pointer to the CLI::App instance.
 * @param define_options A function that configures the CLI options. It takes a raw pointer to the CLI::App and a
 * shared pointer to CliOpts.
 * @param cli_opts A shared pointer to the CliOpts instance, which holds the parsed command line options.
 * @param main A function that contains the main logic to be executed, even if subcommands are executed and fallthrough
 * is enabled.
 * @param pre_callback A function that is executed before the main function. Any extra validation or any preprocess
 * on the options which can not be done using the CLI library, can be done here.
 */
template <class CliOpts>
void AddFallbackCommand(CLI::App* const app,
                        const DefineOptions<CliOpts>& define_options,
                        std::shared_ptr<CliOpts>& cli_opts,
                        const Main<CliOpts>& main,
                        const PreCallback<CliOpts>& pre_callback) {
  define_options(app, cli_opts);
  app->callback([app, cli_opts, main, pre_callback] {
    if (pre_callback) {
      pre_callback(app, cli_opts);
    }
    if (main) {
      main(*cli_opts);
    }
  });
}

template <class CliOpts>
void AddFallbackCommand(CLI::App* const app,
                        const DefineOptions<CliOpts>& define_options,
                        std::shared_ptr<CliOpts>& cli_opts,
                        const Main<CliOpts>& main) {
  AddFallbackCommand(app, define_options, cli_opts, main, PreCallback<CliOpts>{});
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
 * @return A pointer to the newly added subcommand's CLI::App instance.
 */
template <class CliOpts>
CLI::App* AddSubcommand(CLI::App* const app,
                        const std::string& subcommand_name,
                        const DefineOptions<CliOpts>& define_sub_command_options,
                        std::shared_ptr<CliOpts>& cli_opts,
                        const Main<CliOpts>& subcommand_main,
                        const std::optional<std::string>& subcommand_description,
                        const PreCallback<CliOpts>& pre_callback) {
  auto* const subcommand = app->add_subcommand(subcommand_name);
  if (subcommand_description.has_value()) {
    subcommand->description(subcommand_description.value());
  }
  define_sub_command_options(subcommand, cli_opts);
  subcommand->callback([subcommand, subcommand_main, cli_opts, pre_callback] {
    if (pre_callback) {
      pre_callback(subcommand, cli_opts);
    }
    subcommand_main(*cli_opts);
  });
  return subcommand;
}

template <class CliOpts>
CLI::App* AddSubcommand(CLI::App* const app,
                        const std::string& subcommand_name,
                        const DefineOptions<CliOpts>& define_sub_command_options,
                        std::shared_ptr<CliOpts>& cli_opts,
                        const Main<CliOpts>& subcommand_main,
                        const std::optional<std::string>& subcommand_description) {
  return AddSubcommand(app,
                       subcommand_name,
                       define_sub_command_options,
                       cli_opts,
                       subcommand_main,
                       subcommand_description,
                       PreCallback<CliOpts>{});
}

template <class CliOpts>
CLI::App* AddSubcommand(CLI::App* const app,
                        const std::string& subcommand_name,
                        const DefineOptions<CliOpts>& define_sub_command_options,
                        std::shared_ptr<CliOpts>& cli_opts,
                        const Main<CliOpts>& subcommand_main) {
  return AddSubcommand(app, subcommand_name, define_sub_command_options, cli_opts, subcommand_main, std::nullopt);
}

template <class CliOpts>
struct StandardMainParam {
  std::string program_name;
  std::string version;
  std::shared_ptr<CliOpts> cli_opts{};
  DefineOptions<CliOpts> define_options{};
  Main<CliOpts> main{};
  PreCallback<CliOpts> pre_callback{};
};

template <class CliOpts>
s32 StandardMain(s32 argc, const char* const* argv, StandardMainParam<CliOpts>& param) {
  const auto app = SetupDefaultCli(param.program_name, param.version);
  AddFallbackCommand(app.get(), param.define_options, param.cli_opts, param.main, param.pre_callback);
  return RunCli(app.get(), argc, argv);
}

template <typename... Args>
CLI::ValidationError ValidationError(fmt::format_string<Args...> fmt, Args&&... args) {
  return CLI::ValidationError(fmt::format(fmt, std::forward<Args>(args)...));
}
}  // namespace xoos::cli
