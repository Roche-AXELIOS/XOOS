#include "clustering/cluster-downsample.h"

#include <random>

namespace xoos::read_collapser {

void DownsampleReadsInCluster(vec<AlignmentPtr>& reads, const u32 max_reads) {
  // Partition reads into forward and reverse; and within those, into full and partial reads.
  if (reads.size() <= max_reads) {
    // If the number of reads is less than or equal to max_reads, return them as
    // no subsampling is needed.
    return;
  }
  // NOLINTNEXTLINE (cert-msc51-cpp)
  std::mt19937 random_generator(0);  // use constant seed to preserve determinacy for testing
  // Shuffle the reads to ensure randomness in subsampling.
  std::shuffle(reads.begin(), reads.end(), random_generator);
  // Partition reads into full and partial reads, and then into forward and reverse reads.
  const auto partial_reads_start =
      std::partition(reads.begin(), reads.end(), [](const AlignmentPtr& read) { return !read->IsPartial(); });
  const auto full_rev_start =
      std::partition(reads.begin(), partial_reads_start, [](const AlignmentPtr& read) { return read->IsForward(); });
  const auto partial_rev_start =
      std::partition(partial_reads_start, reads.end(), [](const AlignmentPtr& read) { return read->IsForward(); });
  // [reads.begin(), full_rev_start) are full forward reads
  // [full_rev_start, partial_reads_start) are full reverse reads
  // [partial_reads_start, partial_rev_start) are partial forward reads
  // [partial_rev_start, reads.end()) are partial reverse reads
  vec<AlignmentPtr> subsampled_reads;
  subsampled_reads.reserve(max_reads);
  // Try to copy a balanced number of full forward and full reverse reads first
  s32 total_selected_reads = 0;
  s32 full_reads_processed = 0;
  const auto full_forward_count = std::distance(reads.begin(), full_rev_start);
  const auto full_reverse_count = std::distance(full_rev_start, partial_reads_start);
  while (total_selected_reads < ToSigned(max_reads) && full_reads_processed < full_forward_count &&
         full_reads_processed < full_reverse_count) {
    subsampled_reads.push_back(*(reads.begin() + full_reads_processed));   // Full forward
    subsampled_reads.push_back(*(full_rev_start + full_reads_processed));  // Full reverse
    ++full_reads_processed;
    total_selected_reads += 2;
  }
  // If we still have space, fill with full reads
  while (total_selected_reads < ToSigned(max_reads) && full_reads_processed < full_forward_count) {
    subsampled_reads.push_back(*(reads.begin() + full_reads_processed));
    ++full_reads_processed;
    ++total_selected_reads;
  }
  while (total_selected_reads < ToSigned(max_reads) && full_reads_processed < full_reverse_count) {
    subsampled_reads.push_back(*(full_rev_start + full_reads_processed));
    ++full_reads_processed;
    ++total_selected_reads;
  }
  // If we still have space, fill with partial reads
  s32 partial_reads_processed = 0;
  const auto partial_fwd_count = std::distance(partial_reads_start, partial_rev_start);
  const auto partial_rev_count = std::distance(partial_rev_start, reads.end());
  while (total_selected_reads < ToSigned(max_reads) && partial_reads_processed < partial_fwd_count &&
         partial_reads_processed < partial_rev_count) {
    subsampled_reads.push_back(*(partial_reads_start + partial_reads_processed));  // Partial forward
    subsampled_reads.push_back(*(partial_rev_start + partial_reads_processed));    // Partial reverse
    ++partial_reads_processed;
    total_selected_reads += 2;
  }
  while (total_selected_reads < ToSigned(max_reads) && partial_reads_processed < partial_fwd_count) {
    subsampled_reads.push_back(*(partial_reads_start + partial_reads_processed));
    ++partial_reads_processed;
    ++total_selected_reads;
  }
  while (total_selected_reads < ToSigned(max_reads) && partial_reads_processed < partial_rev_count) {
    subsampled_reads.push_back(*(partial_rev_start + partial_reads_processed));
    ++partial_reads_processed;
    ++total_selected_reads;
  }
  // During balanced sampling, we copied reads in pairs so we may exceed max_reads by 1.
  if (subsampled_reads.size() > max_reads) {
    subsampled_reads.resize(max_reads);
  }
  // Swap the original reads with the subsampled reads.
  reads = std::move(subsampled_reads);
}

}  // namespace xoos::read_collapser
