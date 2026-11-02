#include "cli/demux.h"

int main(int argc, char** argv) {
  const auto callbacks = xoos::demux::CreateCliCallbacks();
  xoos::cli::StandardMainParam<xoos::demux::DemuxAndTrimParam> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<xoos::demux::DemuxAndTrimParam>(),
      .define_options = callbacks.define_options,
      .main = xoos::demux::DemuxAndTrimPipeline,
      .pre_callback = callbacks.pre_callback,
  };
  return xoos::cli::StandardMain(argc, argv, param);
}
