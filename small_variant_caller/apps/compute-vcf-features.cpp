#include "compute-vcf-features/compute-vcf-features-cli.h"

namespace xoos::svc {

static int Main(int argc, char** argv) {
  cli::StandardMainParam<ComputeVcfFeaturesParam> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<ComputeVcfFeaturesParam>(),
      .define_options = DefineOptionsComputeVcfFeatures,
      .main = ComputeVcfFeatures,
  };
  return StandardMain(argc, argv, param);
}

}  // namespace xoos::svc

int main(int argc, char** argv) {
  return xoos::svc::Main(argc, argv);
}
