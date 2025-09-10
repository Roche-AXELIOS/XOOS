#include "compute-bam-features/compute-bam-features-cli.h"

namespace xoos::svc {

static int Main(int argc, char** argv) {
  cli::StandardMainParam<ComputeBamFeaturesCliParams> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<ComputeBamFeaturesCliParams>(),
      .define_options = DefineOptions,
      .main = ParallelComputeBamFeatures,
  };
  return StandardMain(argc, argv, param);
}

}  // namespace xoos::svc

int main(int argc, char** argv) {
  return xoos::svc::Main(argc, argv);
}
