#pragma once

#include <string>

#include <xoos/cli/cli.h>
#include <xoos/types/str-container.h>
#include <xoos/yc-decode/yc-decoder.h>

namespace xoos::svc {
// Default values for various parameters used in compute-bam-features
constexpr float kDefaultMinCTDNAAf = 0.0;
constexpr float kDefaultMinFFPEAf = 0.01;
const std::string kDefaultSampleType = "ctDNA";
const std::string kDefaultBamFeaturesFileName = "bam_features.txt";
constexpr u32 kDefaultMaxBamRegionSizePerThread = 16384;
constexpr u32 kDefaultMaxVcfRegionSizePerThread = 64000;
const std::string kDefaultMinBaseType = "simplex";

const StrMap<yc_decode::BaseType> kBaseTypeStrMap{
    {"discordant", yc_decode::BaseType::kDiscordant},
    {"simplex", yc_decode::BaseType::kSimplex},
    {"concordant", yc_decode::BaseType::kConcordant},
};

// Default values for parameters used in compute-vcf-features
const std::string kDefaultVcfFeaturesFileName = "vcf_features.txt";

// Default values for various parameters used in filter-variants
const std::string kDefaultChr1Name = "chr1";

// Default values for various parameters used in train-model
constexpr u32 kMaxScore = 4;

// Default values for various parameters used in vcf-to-bed
const std::string kDefaultVcfToBedFileName = "vcf-regions.bed";
constexpr u32 kDefaultLeftPad = 0;
constexpr u32 kDefaultRightPad = 0;
constexpr u32 kDefaultCollapsableDist = 0;

// inclusive range for MAPQ in SAM spec
static const CLI::Range kCliRangeMapq(0, 255);

// inclusive range for QUAL in SAM spec
static const CLI::Range kCliRangeBaseq(0, 93);

// inclusive range for fractional values
static const CLI::Range kCliRangeFraction(0.0, 1.0);

// CLI validator for non-empty input file
struct NonEmptyFileValidator : CLI::Validator {
  NonEmptyFileValidator();
};

// CLI validator for an indexed input BAM file
struct IndexedBamFileValidator : CLI::Validator {
  IndexedBamFileValidator();
};

// CLI validator for an indexed input VCF file
struct IndexedVcfFileValidator : CLI::Validator {
  IndexedVcfFileValidator();
};

// CLI validator for an indexed input FASTA file
struct IndexedFastaFileValidator : CLI::Validator {
  IndexedFastaFileValidator();
};

// CLI validator for an input BED file
struct BedFileValidator : CLI::Validator {
  BedFileValidator();
};

static const CLI::Validator kCliNonEmptyFile = CLI::ExistingFile & NonEmptyFileValidator();
static const CLI::Validator kCliIndexedBamFile = kCliNonEmptyFile & IndexedBamFileValidator();
static const CLI::Validator kCliIndexedVcfFile = kCliNonEmptyFile & IndexedVcfFileValidator();
static const CLI::Validator kCliIndexedFastaFile = kCliNonEmptyFile & IndexedFastaFileValidator();
static const CLI::Validator kCliBedFile = kCliNonEmptyFile & BedFileValidator();

CLI::Option* AddWarnAsErrorOption(cli::AppPtr app);

}  // namespace xoos::svc
