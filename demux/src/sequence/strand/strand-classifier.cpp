#include "strand-classifier.h"

#include <cmath>
#include <stdexcept>

#include <boost/math/distributions/binomial.hpp>

#include "sequence/kmer/kmerizer.h"

namespace xoos::demux::strand {

/**
 * @brief Constructor for StrandClassifier
 * @param ibf InterleavedBloomFilter object
 * @param noise Probability of false positive, between 0 and 1 (typically ibf average FPR)
 * @param critical Critical value for classification (Upper bound FPR of misclassification)
 * @param precompute_kmer_num Number of kmers to precompute thresholds for (if not enough, will calculate on-the-fly)
 **/
StrandClassifier::StrandClassifier(const InterleavedBloomFilter& ibf, double noise, double critical,
                                   int precompute_kmer_num)
    : _ibf(ibf),
      _noise(noise),
      _critical(critical),
      _min_kmer_count(ComputeMinCriticalThreshold(_noise, _critical)),
      _crit_thresholds(ComputeCriticalThresholds(precompute_kmer_num, noise, critical, _min_kmer_count)) {}

/**
 * @brief Classify the strand for a given read. Stop early if the difference in strand counts is statistically
 * significant. Otherwise, will return the strand with the highest count.
 * @param read a sequence as a string to classify (is k-merized on the fly)
 * @return The direction of the read (forward, reverse, or unknown) and if it is significant or not
 */
StrandType StrandClassifier::Classify(std::string_view read) const noexcept {
  // Basic idea, check if the difference between the forward and reverse counts could have arisen by noise
  // Noise is the probability of a false positive and is determined by filter
  kmer::Kmerizer kmerizer(read, _ibf.KmerSize());

  int fw_counts = 0;
  int rv_counts = 0;
  int total_kmers_seen = 0;

  // use kmerizer iterator to iterate over the read
  for (const auto& kmer_pair : kmerizer) {
    ++total_kmers_seen;
    // direction
    bool strand = kmer_pair.fw > kmer_pair.rv;
    // canonical k-mer
    u64 kmer = strand ? kmer_pair.fw : kmer_pair.rv;
    auto match = _ibf.Contains(kmer);

    if (match.first && match.second) {
      // Both direction match, do not increment any count
      // This is caused by at least one false positive in the filter
      continue;
    } else if (match.first) {
      // Only first filter matches (stores forward canonical k-mers in the relative to the reference)
      if (strand) {
        ++fw_counts;
      } else {
        ++rv_counts;
      }
    } else if (match.second) {
      // Only second filter matches (stores reverse canonical k-mers in the relative to the reference)
      if (strand) {
        ++rv_counts;
      } else {
        ++fw_counts;
      }
    }

    if (fw_counts > rv_counts) {
      if ((fw_counts - rv_counts) >= GetCriticalThreshold(total_kmers_seen)) {
        return StrandType::kForwardSig;
      }
    } else if (rv_counts > fw_counts) {
      if ((rv_counts - fw_counts) >= GetCriticalThreshold(total_kmers_seen)) {
        return StrandType::kReverseSig;
      }
    }
  }

  // we haven't reached significance, but we can still classify
  if (fw_counts > rv_counts) {
    return StrandType::kForward;
  } else if (rv_counts > fw_counts) {
    return StrandType::kReverse;
  }
  return StrandType::kUnknown;
}

/**
 * @brief Find the critical threshold for a given total number of kmers seen given the noise and critical value
 * @param total_kmers_seen Number of kmers seen until now
 * @param noise Probability of false positive
 * @param critical Critical value for classification (Upper bound FPR of misclassification)
 * @return Critical threshold of matches (required counts) for classification
 **/
int StrandClassifier::ComputeCriticalThreshold(const int total_kmers_seen, const double noise, const double critical) {
  boost::math::binomial_distribution<> dist(static_cast<double>(total_kmers_seen), noise);
  return std::floor(boost::math::quantile(complement(dist, critical)) + 1.0);
}

/**
 * @brief Compute the minimum number of kmers needed to significantly classify assuming only matches are seen
 * @param noise Probability of false positive (typically ibf average FPR)
 * @param critical Critical value for classification (Upper bound FPR of misclassification)
 * @return Minimum number of kmers needed to classify significantly
 * @throw std::invalid_argument if noise or critical is not in the range (0, 1)
 * @called by the constructor to compute the minimum number of kmers needed to classify and validate the input
 **/
int StrandClassifier::ComputeMinCriticalThreshold(const double noise, const double critical) {
  // Check if alpha and gamma are in the valid range
  if (noise <= 0.0 || noise >= 1.0 || critical <= 0.0 || critical >= 1.0) {
    throw std::invalid_argument("Invalid noise or critical_m2 value");
  }
  // Need smallest N such that noise^N <= critical
  // N >= log(critical) / log(noise)
  double min_n = std::log(critical) / std::log(noise);
  return std::ceil(min_n);
}

/**
 * @brief Precompute the critical thresholds for a range of different total kmers seen
 * @param precompute_kmer_num Number of kmers to precompute thresholds for
 * @param noise Probability of false positive (typically ibf average FPR)
 * @param critical Critical value for classification (Upper bound FPR of misclassification)
 * @param min_kmer_count Minimum number of kmers needed to classify
 * @return Vector of critical thresholds for each n value
 * @called by the constructor to precompute the thresholds for a range of n values
 **/
std::vector<int> StrandClassifier::ComputeCriticalThresholds(int precompute_kmer_num, double noise, double critical,
                                                             int min_kmer_count) {
  // initialize the vector with the minimum count just to be safe
  std::vector<int> m2_critical_count_thresholds(precompute_kmer_num, min_kmer_count);
  for (int kmer_count = min_kmer_count; kmer_count < precompute_kmer_num; ++kmer_count) {
    m2_critical_count_thresholds[kmer_count] = ComputeCriticalThreshold(kmer_count, noise, critical);
    boost::math::binomial_distribution<> binom_dist(kmer_count, noise);
  }
  return m2_critical_count_thresholds;
}

/**
 * @brief Get the critical threshold for a given total number of kmers seen
 * @param total_kmers_seen Number of kmers seen until now
 * @return Critical threshold (required counts) for classification
 * @details Computes the critical threshold for a given n if it is not precomputed on-the-fly
 **/
int StrandClassifier::GetCriticalThreshold(int total_kmers_seen) const {
  if (total_kmers_seen > static_cast<int>(_crit_thresholds.size())) {
    // Calculate on-the-fly
    return ComputeCriticalThreshold(total_kmers_seen, _noise, _critical);
  }
  return _crit_thresholds[total_kmers_seen];
}
}  // namespace xoos::demux::strand
