#include "demux-strand-index/strand-build.h"

int main(int argc, char** argv) {
  xoos::cli::StandardMainParam<xoos::demux::strand::StrandBuildParam> param = {
      .cli_opts = std::make_shared<xoos::demux::strand::StrandBuildParam>(),
      .define_options = xoos::demux::strand::DefineOptions,
      .main = xoos::demux::strand::StrandBuildPipeline,
  };
  return StandardMain(argc, argv, param);
}
