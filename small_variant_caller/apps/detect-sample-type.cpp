#include "file_type_inference/detect-sample-type-cli.h"

namespace xoos::svc {

namespace {
int Main(int argc, char** argv) {
  cli::StandardMainParam<DetectSampleTypeParam> param = {
      .program_name = "detect_sample_type",
      .version = VERSION,
      .cli_opts = std::make_shared<DetectSampleTypeParam>(),
      .define_options = DefineOptionsDetectSampleType,
      .main = DetectSampleType,
  };
  return StandardMain(argc, argv, param);
}
}  // namespace

}  // namespace xoos::svc

int main(int argc, char** argv) {
  return xoos::svc::Main(argc, argv);
}
