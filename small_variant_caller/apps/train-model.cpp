#include "train-model/train-model-cli.h"

namespace xoos::svc {

static int Main(int argc, char** argv) {
  cli::StandardMainParam<TrainModelParam> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<TrainModelParam>(),
      .define_options = train_model::DefineOptions,
      .main = TrainModel,
      .pre_callback = train_model::PreCallback,
  };
  return StandardMain(argc, argv, param);
}

}  // namespace xoos::svc

int main(int argc, char** argv) {
  return xoos::svc::Main(argc, argv);
}
