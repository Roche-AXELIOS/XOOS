#include "phaser/phaser-cli.h"
#include "phaser/run.h"

namespace xoos::svc {

namespace {

int Main(int argc, char** argv) {
  cli::StandardMainParam<phaser::Params> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<phaser::Params>(),
      .define_options = phaser::DefineOptions,
      .main = phaser::Run,
      .pre_callback = phaser::PreCallback,
  };
  return StandardMain(argc, argv, param);
}

}  // namespace

}  // namespace xoos::svc

int main(int argc, char** argv) {
  return xoos::svc::Main(argc, argv);
}
