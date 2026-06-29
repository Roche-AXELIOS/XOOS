#include <xoos/cli/cli.h>

#include "core/read-collapser-cli.h"
#include "core/read-collapser-options.h"

int main(int argc, char** argv) {
  xoos::cli::StandardMainParam<xoos::read_collapser::ReadCollapserOptions> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<xoos::read_collapser::ReadCollapserOptions>(),
      .define_options = xoos::read_collapser::DefineOptions,
      .pre_callback = xoos::read_collapser::SetCommandLineInfo,
  };
  return StandardMain(argc, argv, param);
}
