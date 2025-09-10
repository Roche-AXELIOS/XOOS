#include "vcf-to-bed/vcf-to-bed-cli.h"

namespace xoos::svc {

static int Main(int argc, char** argv) {
  cli::StandardMainParam<VcfToBedParam> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<VcfToBedParam>(),
      .define_options = DefineOptionsVcfToBed,
      .main = ConvertVcfToBed,
  };
  return StandardMain(argc, argv, param);
}

}  // namespace xoos::svc

int main(int argc, char** argv) {
  return xoos::svc::Main(argc, argv);
}
