#include "features-to-vcf/features-to-vcf-cli.h"

namespace xoos::svc {

namespace {
int Main(int argc, char** argv) {
  cli::StandardMainParam<FeaturesToVcfParam> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<FeaturesToVcfParam>(),
      .define_options = DefineOptionsFeaturesToVcf,
      .main = ConvertFeaturesToVcf,
  };
  return StandardMain(argc, argv, param);
}
}  // namespace

}  // namespace xoos::svc

int main(int argc, char** argv) {
  return xoos::svc::Main(argc, argv);
}
