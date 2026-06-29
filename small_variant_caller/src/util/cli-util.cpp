#include "cli-util.h"

#include "core/cli-option-names.h"
#include "log-util.h"
#include "xoos/cli/validators/bed-validator.h"
#include "xoos/cli/validators/indexed-bam-validator.h"
#include "xoos/cli/validators/indexed-fasta-validator.h"
#include "xoos/cli/validators/indexed-vcf-validator.h"
#include "xoos/cli/validators/nonempty-file-validator.h"

namespace xoos::svc {

CLI::Option* AddWarnAsErrorOption(CLI::App* const app) {
  return app->add_flag(cli_opt_name::kWarnAsError, warn_as_error, "Treat warn messages as errors");
}

CLI::Option* CheckNonEmptyFile(CLI::Option* const opt) {
  return opt->check(CLI::ExistingFile)->check(cli::NonEmptyFileValidator());
}

CLI::Option* CheckIndexedBamFile(CLI::Option* const opt) {
  return CheckNonEmptyFile(opt)->check(cli::IndexedBamFileValidator());
}

CLI::Option* CheckIndexedVcfFile(CLI::Option* const opt) {
  return CheckNonEmptyFile(opt)->check(cli::IndexedVcfFileValidator());
}

CLI::Option* CheckIndexedFastaFile(CLI::Option* const opt) {
  return CheckNonEmptyFile(opt)->check(cli::IndexedFastaFileValidator());
}

CLI::Option* CheckBedFile(CLI::Option* const opt) {
  return CheckNonEmptyFile(opt)->check(cli::BedFileValidator());
}

}  // namespace xoos::svc
