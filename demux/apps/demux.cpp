#include "cli/demux.h"

int main(int argc, char** argv) {
  const std::string version = VERSION;
  xoos::cli::StandardMainParam<xoos::demux::DemuxAndTrimParam> param = {
      .program_name = PROGRAM_NAME,
      .version = version,
      .cli_opts = std::make_shared<xoos::demux::DemuxAndTrimParam>(),
      .define_options = xoos::demux::DefineOptions,
      .main = xoos::demux::DemuxAndTrimPipeline,
      .pre_callback = xoos::demux::CreatePreCallback(version),
  };
  return xoos::cli::StandardMain(argc, argv, param);
}
