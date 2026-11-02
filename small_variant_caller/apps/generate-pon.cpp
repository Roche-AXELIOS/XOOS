#include "generate-pon/generate-pon-cli.h"

namespace xoos::svc {

namespace {
int Main(int argc, char** argv) {
  cli::StandardMainParam<GeneratePonParam> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<GeneratePonParam>(),
      .define_options = DefineOptionsGeneratePon,
      .main = GeneratePon,
      .pre_callback = generate_pon::PreCallback,
  };
  return StandardMain(argc, argv, param);
}
}  // namespace

}  // namespace xoos::svc

int main(int argc, char** argv) {
  return xoos::svc::Main(argc, argv);
}
