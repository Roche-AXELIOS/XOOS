#pragma once

#include <xoos/io/fasta-reader.h>
#include <xoos/types/float.h>

#include <filesystem>

#include <gtl/phmap.hpp>
#include <taskflow/taskflow.hpp>

#include "bloom-filter/interleaved-bloom-filter.h"
#include "sequence/strand/strand-kmer-table.h"

namespace xoos::demux::strand {

// Basic Idea: Populate a hash table stores k-mers and their strand on the reference genome.
//             We use this table build a smaller data structure for strand detection.
//
// Store information need to keep track of each k-mers direction presence in the reference
// To simplify the logic, it may be better to use odd sized k-mers due to prevent palindromes
// Basic idea: Each k-mer in the reference genome can exist in the following states:
// - Occurs in canonical form in forward strand
// - Occurs in canonical form in reverse strand
// - Occurs in both forward and reverse strands (canonical form is present in both strands)
//
// Direction flags are stored in a single byte:
// - 0x01: canonical form is present in forward strand of reference
// - 0x02: canonical form is present in reverse strand of reference
// - The remaining bits are unused
// - A possible alternative scheme is to use 4 bits for each direction, allowing for counts up to 15
namespace fs = std::filesystem;

struct StrandBuildParam {
  /// @brief input fasta file, must have a fai index
  fs::path input;
  /// @brief output will be written to this filename
  fs::path output;
  /// @brief The number of threads to be used for parallel processing.
  size_t threads;
  /// @brief The maximum false positive rate for the bloom filter
  f32 max_false_positive_rate;
  /// @brief kmer size
  size_t kmer_size;
};

void BuildStrandMap(io::FastaReader& reader, tf::Executor& executor, StrandKmerTable& strand_map, u32 kmer_size);

size_t CountStrandKmers(const std::vector<std::atomic<u64>>& strand_table, tf::Executor& executor,
                        size_t& num_unique_kmers_fw_total, size_t& num_unique_kmers_rv_total);

void BuildStrandFilter(const std::vector<std::atomic<u64>>& strand_table, InterleavedBloomFilter& bf,
                       tf::Executor& executor);

void StrandBuildPipeline(const StrandBuildParam& param);

}  // namespace xoos::demux::strand
