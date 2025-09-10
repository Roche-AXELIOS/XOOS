#include "filter-variants/filter-variants-cli.h"

namespace xoos::svc {

static int Main(int argc, char** argv) {
  cli::StandardMainParam<FilterVariantsParam> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<FilterVariantsParam>(),
      .define_options = DefineOptionsFilterVariants,
      .main = FilterVariants,
      .pre_callback = FilterVariantsPreCallback,
  };
  return StandardMain(argc, argv, param);
}

}  // namespace xoos::svc

int main(int argc, char** argv) {
  return xoos::svc::Main(argc, argv);
}
