#pragma once

#include <vector>

#include "bloom-filter/interleaved-bloom-filter.h"
#include "strand.h"

namespace xoos::demux::strand {

/**
 * @Class StrandClassifier
 * @brief Classifies a read as either forward or reverse strand using canoncial k-mers from an Interleaved Bloom Filter.
 *
 * This class uses an Interleaved Bloom Filter to classify reads based on their strand. Assuming the noise specified
 * (typically Bloom filter FPR) is correct, the critical value will bound the upper bound FPR of assigning a read to
 * the wrong strand. It is an upper bound because it uses the difference in counts to each strand rather than just
 * count matches to a filter to more robustly handle cases where the read may actually contain k-mers from both strands.
 *
 * @details
 * We set up a Binomial distribution that can model if the number of matches seen to a strand could have arisen by the
 * noise of the IBF (i.e. FPR). We then use it to test if differences in counts (higher strand_count - lower
 * stand_count) is could have have arisen due to random chance during classification and end early if this happens to
 * be the case. This is a heuristic that is more stringent than just testing if enough k-mers match a strand (how one
 * might typical use of a binomial distribution modeled on the noise). The more typical approach would use all matches
 * (and not consider other strand), but taking the difference is a more robust to make sure to handle cases where the
 * read may actually contain k-mers from both strands while still leveraging low noise. Using this heuristic, the
 * critical value become an upper bound to the true FPR of misclassification as it it strictly always more stringent
 *than just considering the number of matches to a strand.
 *
 * @author Justin Chu
 * @date 2025-05-05
 **/
class StrandClassifier {
 public:
  StrandClassifier(const InterleavedBloomFilter& ibf, double noise, double critical, int precompute_kmer_num = 1024);

  // Disable copy/move semantics
  StrandClassifier(const StrandClassifier&) = delete;
  StrandClassifier& operator=(const StrandClassifier&) = delete;
  StrandClassifier(StrandClassifier&&) = delete;
  StrandClassifier& operator=(StrandClassifier&&) = delete;

  StrandType Classify(std::string_view read) const noexcept;

 private:
  const InterleavedBloomFilter& _ibf;  // This stores canonical the k-mers to each strand relative to the strand

  const double _noise;        // Noise/expected rate of randomly getting a match, typically the FPR of the Bloom filter
  const double _critical;     // Critical value for the Binomial distribution to see if the match counts the noise
  const int _min_kmer_count;  // Minimum number of k-mers seen to possibly classify a read

  const std::vector<int> _crit_thresholds;  // so we don't need to recompute critical thresholds (minimum matches)

  // Helper functions
  static int ComputeMinCriticalThreshold(double noise, double critical);
  static int ComputeCriticalThreshold(int total_kmers_seen, double noise, double critical);
  static std::vector<int> ComputeCriticalThresholds(int precompute_kmer_num, double noise, double critical,
                                                    int min_kmer_count);
  int GetCriticalThreshold(int total_kmers_seen) const;
};

}  // namespace xoos::demux::strand
