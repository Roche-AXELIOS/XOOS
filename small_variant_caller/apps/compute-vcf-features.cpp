#include "compute-vcf-features/compute-vcf-features-cli.h"

namespace xoos::svc {

static int Main(int argc, char** argv) {
  cli::StandardMainParam<ComputeVcfFeaturesParam> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<ComputeVcfFeaturesParam>(),
      .define_options = compute_vcf_features::DefineOptions,
      .main = ComputeVcfFeatures,
      .pre_callback = compute_vcf_features::PreCallback,
  };
  return StandardMain(argc, argv, param);
}

}  // namespace xoos::svc

int main(int argc, char** argv) {
  return xoos::svc::Main(argc, argv);
}
