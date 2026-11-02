#pragma once
#include <Eigen/Dense>
#include <string>
#include <utility>

#include <xoos/io/fasta-reader.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "io/bigwig.h"
#include "io/fai.h"

namespace xoos::cnc {
const f64 kCoverageMinOnTargetMapp = 0.6;
const f64 kCoverageMinOffTargetMapp = 0.1;
const f64 kCoverageGCQuantile = 0.001;
constexpr f64 kDefaultGCLowerThreshold = 0.15;
constexpr f64 kDefaultGCUpperThreshold = 0.75;

/**
 * @class AugmentedIntervals
 * @brief container for annotated genomic intervals
 *
 * Input file should have the following format:
 * line 1: header
 * line 2-n: <region>\t<gc bias>\t<mappability>\t<rep timing>\t<Gene>\t<on_target>\n
 */
class AugmentedIntervals {
 public:
  AugmentedIntervals() = default;
  AugmentedIntervals(const AugmentedIntervals&) = default;
  AugmentedIntervals(AugmentedIntervals&&) noexcept = default;
  AugmentedIntervals& operator=(const AugmentedIntervals&) = default;
  AugmentedIntervals& operator=(AugmentedIntervals&&) noexcept = default;
  ~AugmentedIntervals() = default;

  explicit AugmentedIntervals(std::istream& ifs);

  const vec<std::string>& GetRegions() const {
    return _region;
  }

  const Eigen::VectorXd& GetMappability() const {
    return _mappability;
  }

  const Eigen::VectorXd& GetGCBias() const {
    return _gc_bias;
  }

  const vec<bool>& GetOnTargetStatus() const {
    return _on_target;
  }

  void SetRegions(vec<std::string>&& regions) {
    _region = std::move(regions);
  }

  void SetMappability(Eigen::VectorXd&& mappability) {
    _mappability = std::move(mappability);
  }

  void SetGCBias(Eigen::VectorXd&& gc_bias) {
    _gc_bias = std::move(gc_bias);
  }

  void SetOnTargetStatus(const vec<bool>& on_target) {
    _on_target = on_target;
  }

  void SetOnTargetStatus(vec<bool>&& on_target) {
    _on_target = std::move(on_target);
  }

  void SetReferenceFile(const std::string& reference_name) {
    _reference_name = reference_name;
  }

  void SetAllOnTarget();

  void SetSeqLengths(SeqLenMap& seq_lengths) {
    _seq_lengths = seq_lengths;
  }

  const std::string& GetReferenceFile() const {
    return _reference_name;
  }

  vec<std::string> GetOnOrOffTargetRegions(bool /*o*/);

  const SeqLenMap& GetSeqLengths() const {
    return _seq_lengths;
  }

  void SortByRegion(const fs::path& fai_file);
  void AnnotateGC(io::FastaReader& fa);
  void AnnotateMappability(BigWig& bigwig);

  void Write(std::ofstream& ofs, const io::CommandLineInfo& command_line_info) const;

  // Remove intervals with GC content outside the specified thresholds
  void FilterExtremeGC(f64 lower_threshold = kDefaultGCLowerThreshold, f64 upper_threshold = kDefaultGCUpperThreshold);

 private:
  void LoadIntervalsFile(std::istream& ifs);
  vec<std::string> _region;
  Eigen::VectorXd _gc_bias;
  Eigen::VectorXd _mappability;
  vec<bool> _on_target;
  SeqLenMap _seq_lengths;
  std::string _reference_name;
};

}  // namespace xoos::cnc
