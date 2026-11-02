#include "demux-strand-index/strand-build.h"

#include <string>

#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/file-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/file-permission-validator.h>
#include <xoos/cli/validators/indexed-fasta-validator.h>

namespace xoos::demux::strand {
using cli::AddThreadCountOption;

namespace {

void AddRemovedOption(CLI::App* const app, const std::string& old_name, const std::string& replacement) {
  app->add_option_function<std::string>(
         old_name,
         [old_name, replacement](const std::string&) {
           throw CLI::ValidationError(
               fmt::format("The option '{}' is no longer supported; use '{}' instead.", old_name, replacement));
         },
         "")
      ->group("");
}

}  // namespace

void DefineOptions(CLI::App* const app, const StrandBuildParamPtr& param) {
  app->description(
      "This tool builds interleaved Bloom filters for strand detection of reads, to enable strand-aware quality base "
      "shifting for indels in duplex mode. All filters are powers of 2 in size, with the hash functions number "
      "determined by the max false positive rate (--max_false_positive_rate).");
  // Using powers of 2 speeds up modulo operations, so the filter is slightly faster.
  app->add_option("--reference", param->input, "Reference FASTA file. Must include a fai index.")
      ->required()
      ->check(cli::IndexedFastaFileValidator());
  // Negating ExistingPath rejects both files and directories; negating ExistingFile would allow an existing directory.
  app->add_option("--output-strand-index", param->output, "Output filename for the strand index.")
      ->default_val("strand.bloom")
      ->check(!CLI::IsMember({std::string{}}))
      ->check(!CLI::ExistingPath)
      ->check(cli::FileWriteableValidator());
  AddRemovedOption(app, "-i", "--reference");
  AddRemovedOption(app, "-o", "--output-strand-index");
  AddRemovedOption(app, "--input", "--reference");
  AddRemovedOption(app, "--output", "--output-strand-index");
  AddThreadCountOption(app, "--threads", param->threads);
  app->add_option("--kmer-size", param->kmer_size, "k-mer size for the output strand index.")
      ->default_val(kDefaultStrandBuildIndexKmerSize)
      ->check(CLI::Range(size_t{1}, kMaxSupportedKmerSize));
  // This default value is (1/2)^6 (i.e. the optimal FPR for 6 hash functions.)
  // This FPR is selected as the default because it happens to minimize the FPR of an exactly 2GB filter of hg38
  // (assuming k 19 and the power of 2 filter size).
  // To use 1GB on hg38 the best FPR would need to be 0.125 (3 hashes). This works, but has slightly worse performance.
  // Going any lower than 0.015625 for hg38 would require 4 gb of memory.
  app->add_option("--max-false-positive-rate",
                  param->max_false_positive_rate,
                  "Maximum false positive rate for the output interleaved bloom filter used as index.")
      ->default_val(kDefaultStrandBuildIndexMaxFalsePositiveRate)
      ->check(CLI::Range(0.0F, 1.0F) & !CLI::IsMember({0.0F, 1.0F}));
}

}  // namespace xoos::demux::strand
