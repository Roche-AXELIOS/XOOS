#pragma once
#include <memory>
#include <optional>

#include <xoos/io/library-type.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "cli/copy-number-caller-modes.h"
#include "cli/sample-metadata-options.h"
#include "cli/vcf-parsing-options.h"
#include "coverage/calculate-coverage-options.h"
#include "gc-correct/gc-correct-options.h"
#include "intervals/augment-intervals-options.h"
#include "likelihood/likelihood-options.h"
#include "merge-segments/merge-segments-options.h"
#include "purity-ploidy-search/purity-ploidy-search-options.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc {

const size_t kCopyNumberCallerDefaultTumorMinCoverage = 0;
const size_t kCopyNumberCallerDefaultSomaticNormalMinCoverage = 30;
const size_t kCopyNumberCallerDefaultGermlineNormalMinCoverage = 0;
inline const fs::path kHg38MappabilityBigWig{"/resources/hg38_mappability.bw"};
inline const fs::path kHg38SbxSMappabilityBigWig{"/resources/hg38_500mer_mappability.bw"};
inline const fs::path kHg38SeedSegments{"/resources/hg38_seed.seg"};

struct CopyNumberCallerOptions {
  // Algorithm Submodule Options
  AugmentIntervalsOptions augment_intervals_options{};
  CalculateCoverageOptions calculate_coverage_options{};
  GCCorrectOptions gc_correct_options{};
  SegmentationOptions segmentation_options{};
  PurityPloidySearchOptions purity_ploidy_search_options{};
  LikelihoodOptions likelihood_options{};
  MergeSegmentsOptions merge_segments_options{};
  BafFilterOptions baf_filter_options{};
  SampleMetadataOptions sample_metadata_options{};

  size_t threads{1};
  CopyNumberCallerModes mode{CopyNumberCallerModes::kUnknown};

  // Somatic Coverage Thresholds
  size_t tumor_min_coverage{};
  size_t normal_min_coverage{};

  // Somatic Input File
  fs::path tumor_bam_fname;
  std::optional<fs::path> tumor_coverage_fname;

  // Sequencing library type of the input BAM. When unset, it is auto-detected
  // from the BAM; an explicit value takes precedence over auto-detection.
  std::optional<io::LibraryType> library_type;

  // Input Files
  fs::path normal_bam_fname;
  fs::path reference_genome_fname;
  fs::path reference_genome_fai_fname;
  fs::path mappability_bigwig_fname;
  fs::path seed_segments_fname;
  std::optional<fs::path> blocklist_bed_fname;
  std::optional<fs::path> augmented_intervals_fname;
  std::optional<fs::path> normal_coverage_fname;
  std::optional<fs::path> vcf_fname;
  std::optional<fs::path> ad_fname;
  std::optional<fs::path> segments_fname;
  std::optional<fs::path> logrs_fname;
  std::optional<fs::path> bafs_fname;
  std::optional<fs::path> mapping_qualities_fname;

  // Predict-sex input
  fs::path predict_sex_coverage_fname;
  // Output Directory
  fs::path output_dir;
  // Prefix prepended to all output filenames (e.g. a sample ID). Empty means filenames are unchanged.
  std::string output_prefix;
  // Overwrite output files if they already exist instead of erroring.
  bool overwrite{false};

  // Extract tumor-only BAF options
  std::optional<fs::path> population_vcf_fname;
  std::optional<fs::path> snp_blocklist_bed_fname;
  f64 min_population_af{};
  s32 min_tumor_depth{};

  // Output-control flags (hidden advanced options)
  bool save_taskflow_graph{false};
  bool save_log_ratio_segments{false};
  bool save_baf_segments{false};

  // Structured metadata for all output formats (TSV and VCF)
  io::CommandLineInfo command_line_info;
};

using CopyNumberCallerOptionsPtr = std::shared_ptr<CopyNumberCallerOptions>;

}  // namespace xoos::cnc
