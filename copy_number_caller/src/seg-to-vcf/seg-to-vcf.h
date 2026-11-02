#pragma once
#include <optional>
#include <string_view>

#include <xoos/io/metadata-util.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

#include "cli/copy-number-caller-modes.h"
#include "cli/copy-number-caller-options.h"
#include "io/fai.h"
#include "likelihood/total-copy-number-prior.h"
#include "seg-to-vcf/vcf-purity-source.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {

inline constexpr std::string_view kSampleWarningHighExtremeBafProportion = "##SampleWarning=highExtremeBafProportion";

struct WriteSegmentsToVcfParams {
  SeqLenMap seq_lengths;
  StrUnorderedMap<size_t> seq_order;
  std::string sample_id;
  std::string reference_file;
  sex_predict::Sex sex{sex_predict::Sex::kUnknown};
  CopyNumberCallerModes mode{CopyNumberCallerModes::kGermlineNormalWGS};
  std::optional<f64> purity;
  std::optional<f64> ploidy;
  std::optional<VcfPuritySource> vcf_purity_source;
  std::optional<TotalCopyNumberPrior> total_copy_number_prior;
  // Structured metadata for VCF ##RocheCommandLine header
  io::CommandLineInfo command_line_info;
  bool has_high_extreme_baf_proportion = false;
};

enum class CopyNumberStatus : u8 {
  kLoss,
  kRef,
  kCnLoh,
  kGain,
  kGainLoh,
  kUnknown
};

void WriteSegmentsToVcfMain(const CopyNumberCallerOptions& options);

void WriteSegmentsToVcf(const fs::path& fname,
                        const vec<GenomicSegment>& segments,
                        const WriteSegmentsToVcfParams& params);

CopyNumberStatus GetCopyNumberStatus(const GenomicSegment& seg);
CopyNumberStatus GetSomaticCopyNumberStatus(const GenomicSegment& seg);
std::string ConstructVcfFilterFieldFromGenomicSegment(const GenomicSegment& seg);
}  // namespace xoos::cnc
