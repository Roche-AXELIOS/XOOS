#pragma once
#include <fstream>
#include <optional>
#include <string>

#include <xoos/io/metadata-util.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "likelihood/total-copy-number-prior.h"

namespace xoos::cnc {

// Scalar configuration for a SegmentLikelihoodWriter. Identifies the pass (prior), the analysis mode (tumor_only), and
// the output destination. Purity/ploidy are NOT stored here: the diagnostic is emitted per grid cell, so each row
// carries its own purity/ploidy supplied at write time.
struct SegmentLikelihoodsParams {
  sex_predict::Sex sex = sex_predict::Sex::kUnknown;
  bool tumor_only = false;
  TotalCopyNumberPrior prior = TotalCopyNumberPrior::kDiploidBias;
  fs::path out_fname;
  io::CommandLineInfo command_line_info;
};

// Per-grid-cell, per-segment somatic-state MLE values for one diagnostic row.
struct SomaticLikelihoodRow {
  // somatic-state logR MLE at this grid cell (CN=2 included)
  f64 logr_likelihood = 0;
  // TCN at the somatic MLE
  size_t total_copy_number = 0;
  // somatic BAF MLE at the selected TCN; absent (N/A) on no-BAF cells
  std::optional<f64> baf_likelihood;
  // somatic joint (logR + BAF); equals logr_likelihood on no-BAF cells
  std::optional<f64> joint_likelihood;
  // minor CN at the selected TCN; absent on no-BAF cells
  std::optional<size_t> minor_copy_number;
};

// Per-segment germline-state MLE values (grid-invariant; tumor-only diploid-bias pass only).
struct GermlineLikelihoodRow {
  // germline-state logR MLE at fixed purity=0.99/ploidy=2 (CN=2 excluded)
  f64 logr_likelihood = 0;
  // TCN at the germline MLE (CN=2 excluded)
  size_t total_copy_number = 0;
};

// Positional identity of one diagnostic row: the segment locus and the purity/ploidy grid cell it was evaluated at.
// Bundled so WriteRow stays within the parameter limit.
struct SegmentLikelihoodRowKey {
  std::string contig;
  size_t start = 0;
  size_t end = 0;
  f64 purity = 0;
  f64 ploidy = 0;
};

/**
 * @brief Streaming, per-row writer for the per-segment purity/ploidy likelihood diagnostic.
 *
 * One writer corresponds to one output file (one purity/ploidy search pass). On construction it opens the file and
 * writes the metadata header, the prior tag, and the column header. Callers then emit one row per contributing segment
 * per grid cell via WriteRow(). Because rows are written inline from the search grid loop, segments excluded by the
 * search never reach the writer and are therefore absent from the file by construction.
 *
 * Column layout depends on the analysis mode (SegmentLikelihoodsParams::tumor_only):
 * - tumor-normal: somatic columns only.
 * - tumor-only: somatic columns plus germline likelihood/TCN, per-cell classification, and contribution likelihood.
 *
 * Germline columns, classification, and contribution are written as N/A on the uniform pass (prior == kUniform), where
 * germline classification is disabled.
 */
class SegmentLikelihoodWriter {
 public:
  explicit SegmentLikelihoodWriter(const SegmentLikelihoodsParams& params);

  // The writer owns an std::ofstream (move-only), so copying is disabled and moves are explicitly noexcept.
  // Constructors are grouped first, then assignment operators (see cpp:S4136).
  SegmentLikelihoodWriter(const SegmentLikelihoodWriter&) = delete;
  SegmentLikelihoodWriter(SegmentLikelihoodWriter&&) noexcept = default;
  SegmentLikelihoodWriter& operator=(const SegmentLikelihoodWriter&) = delete;
  SegmentLikelihoodWriter& operator=(SegmentLikelihoodWriter&&) noexcept = default;
  ~SegmentLikelihoodWriter() = default;

  /**
   * @brief Emits one diagnostic row for a segment at a single purity/ploidy grid cell.
   *
   * @param key Segment locus and grid-cell purity/ploidy that identify this row.
   * @param somatic Somatic-state MLE for this cell.
   * @param germline Per-segment germline-state MLE. Required (and only meaningful) for the tumor-only diploid-bias
   *        pass; pass std::nullopt for tumor-normal mode and the uniform pass, which emit N/A for germline columns.
   * @param classified_germline Whether the segment is classified germline at this cell (tumor-only diploid-bias pass);
   *        std::nullopt otherwise, emitting N/A for the classification and contribution columns.
   */
  void WriteRow(const SegmentLikelihoodRowKey& key,
                const SomaticLikelihoodRow& somatic,
                const std::optional<GermlineLikelihoodRow>& germline,
                std::optional<bool> classified_germline);

 private:
  std::ofstream _ofs;
  bool _tumor_only;
};

}  // namespace xoos::cnc
