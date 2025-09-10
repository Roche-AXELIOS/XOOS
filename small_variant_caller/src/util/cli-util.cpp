#include "cli-util.h"

#include "file-util.h"
#include "log-util.h"

namespace xoos::svc {

/**
 * @brief a custom validator based on a provided validation function.
 * @param validate A function that takes a filesystem path and performs validation.
 * @return A function that takes a string and returns an error message if validation fails.
 */
static std::function<std::string(const std::string&)> MakeValidator(
    const std::function<void(const fs::path&)>& validate) {
  return [validate](const std::string& str) {
    try {
      validate(fs::path(str));
    } catch (const std::exception& e) {
      return std::string(e.what());
    }
    return std::string();
  };
}

/**
 * @brief Constructs a validator that checks if a file exists and is not empty.
 */
NonEmptyFileValidator::NonEmptyFileValidator() {
  name_ = "NONEMPTYFILE";
  func_ = MakeValidator(ValidateNonEmptyFile);
}

/**
 * @brief Constructs a validator that checks if a BAM file and its index exist and are valid.
 */
IndexedBamFileValidator::IndexedBamFileValidator() {
  name_ = "BAMFILE";
  func_ = MakeValidator(ValidateBamAndIndex);
}

/**
 * @brief Constructs a validator that checks if a VCF file and its index exist and are valid.
 */
IndexedVcfFileValidator::IndexedVcfFileValidator() {
  name_ = "VCFFILE";
  func_ = MakeValidator(ValidateVcfAndIndex);
}

/**
 * @brief Constructs a validator that checks if a FASTA file and its index exist and are valid.
 */
IndexedFastaFileValidator::IndexedFastaFileValidator() {
  name_ = "FASTAFILE";
  func_ = MakeValidator(ValidateFastaAndIndex);
}

/**
 * @brief Constructs a validator that checks if a BED file is valid.
 */
BedFileValidator::BedFileValidator() {
  name_ = "BEDFILE";
  func_ = MakeValidator(ValidateBed);
}

/**
 * @brief Adds a command line option to treat warnings as errors.
 * @param app Pointer to the CLI application instance.
 * @return Pointer to the added CLI option.
 */
CLI::Option* AddWarnAsErrorOption(cli::AppPtr app) {
  return app->add_flag("--warn-as-error", warn_as_error, "Treat warn messages as errors");
}

}  // namespace xoos::svc
