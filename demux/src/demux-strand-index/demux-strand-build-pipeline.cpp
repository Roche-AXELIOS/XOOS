#include "demux-strand-index/demux-strand-build-pipeline.h"

#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "sequence/kmer/kmerizer.h"

namespace xoos::demux::strand {

/**
 * @brief Build the strand kmer map which determines unique strand k-mers using a hash table.
 * @param reader FastaReader object to read the reference genome.
 * @param executor Thread pool executor.
 * @param strand_map StrandKmerTable object used to store canonical k-mers (larger complement) and strand information.
 * @param kmer_size Size of the kmer.
 * @details This function reads the reference genome in chunks and builds the kmer strand map.
 */
void BuildStrandMap(io::FastaReader& reader, tf::Executor& executor, StrandKmerTable& strand_map, u32 kmer_size) {
  int fastq_chunk_size = 65536;
  tf::Taskflow taskflow;
  // iterate over the reference genome and build the kmer strand map
  for (const auto& record : reader) {
    const int end = record.second;
    for (int start = 0; start < end; start += fastq_chunk_size) {
      taskflow.emplace([&reader, record, kmer_size, start, end, fastq_chunk_size, &strand_map]() {
        // include overhang of kmer size bases
        const auto& seq = reader.GetSequence(record.first, start,
                                             std::min(end, start + fastq_chunk_size + static_cast<int>(kmer_size)));
        // kmerize the sequence and add to the kmer strand map
        // k-merizer iterator
        for (const auto& kmer_pair : kmer::Kmerizer(seq, kmer_size)) {
          // occur in both (forward and reverse) strands (palindromic), we don't need it at all
          if (kmer_pair.fw == kmer_pair.rv) {
            continue;
          }
          const bool strand = kmer_pair.fw > kmer_pair.rv;  // 1 for forward, 0 for reverse

          // determine canonical k-mer strand direction (larger of 2 k-mers is the canonical k-mer)
          // Note: Usually the smaller k-mer is the canonical k-mer but we use 0 as an empty key in custom table
          const u64 canon_kmer = strand ? kmer_pair.fw : kmer_pair.rv;

          // if k-mer is not in the map, then insert the strand. if it exists, then update the strand using or
          strand_map.InsertOrUpdate(canon_kmer, strand);
        }
      });
    }
  }

  executor.run(taskflow).get();
}

/**
 * @brief Count the number of unique k-mers in the raw strand table (k-mers and strand information, with empty entries).
 * @param strand_table StrandKmerTable raw array used to store canonical k-mers and strand information.
 * @param executor Thread pool executor.
 * @param num_unique_kmers_fw_total Total number of unique k-mers in the forward strand.
 * @param num_unique_kmers_rv_total Total number of unique k-mers in the reverse strand.
 * @return Total number of unique k-mers in both strands.
 */
size_t CountStrandKmers(const std::vector<std::atomic<u64>>& strand_table, tf::Executor& executor,
                        size_t& num_unique_kmers_fw_total, size_t& num_unique_kmers_rv_total) {
  const auto threads = executor.num_workers();
  const size_t map_chunk_size = strand_table.size() / threads + 1;

  std::vector<size_t> num_unique_kmers_fw(threads, 0);
  std::vector<size_t> num_unique_kmers_rv(threads, 0);
  std::vector<size_t> num_unique_kmers_both(threads, 0);

  tf::Taskflow taskflow;

  for (size_t i = 0; i < threads; ++i) {
    size_t start_index = i * map_chunk_size;
    size_t end_index = (i == threads - 1) ? strand_table.size() : (i + 1) * map_chunk_size;
    taskflow.emplace([start_index, end_index, &strand_table, &num_unique_kmers_fw, &num_unique_kmers_rv,
                      &num_unique_kmers_both, i]() {
      for (size_t j = start_index; j < end_index; ++j) {
        const u64 encoded_pair = strand_table.at(j).load(std::memory_order_relaxed);
        if (encoded_pair == StrandKmerTable::kEmptyEntry) {
          continue;
        }
        const auto& [kmer, strand] = StrandKmerTable::DecodeEntry(encoded_pair);
        if (strand == kForwardStrandBit) {
          num_unique_kmers_fw[i]++;
        } else if (strand == kReverseStrandBit) {
          num_unique_kmers_rv[i]++;
        } else {
          num_unique_kmers_both[i]++;
        }
      }
    });
  }

  executor.run(taskflow).get();

  // Accumulate thread-local counts into global totals
  size_t num_unique_kmers_both_total = 0;

  for (size_t i = 0; i < threads; ++i) {
    num_unique_kmers_fw_total += num_unique_kmers_fw[i];
    num_unique_kmers_rv_total += num_unique_kmers_rv[i];
    num_unique_kmers_both_total += num_unique_kmers_both[i];
  }

  return num_unique_kmers_both_total;
}

/**
 * @brief Build strand filters using the raw strand table (k-mers and strand information, with empty entries).
 * @param strand_table StrandKmerTable raw array used to store canonical k-mers and strand information.
 * @param bf Interleaved Bloom Filter with strand information
 * @param executor Thread pool executor.
 */
void BuildStrandFilter(const std::vector<std::atomic<u64>>& strand_table, InterleavedBloomFilter& bf,
                       tf::Executor& executor) {
  const auto threads = executor.num_workers();
  const size_t map_chunk_size = strand_table.size() / threads + 1;

  tf::Taskflow taskflow;
  for (size_t i = 0; i < threads; ++i) {
    size_t start_index = i * map_chunk_size;
    size_t end_index = (i == threads - 1) ? strand_table.size() : (i + 1) * map_chunk_size;
    taskflow.emplace([start_index, end_index, &strand_table, &bf]() {
      for (size_t j = start_index; j < end_index; ++j) {
        const u64 encoded_pair = strand_table.at(j).load(std::memory_order_relaxed);
        if (encoded_pair == StrandKmerTable::kEmptyEntry) {
          continue;
        }
        const auto& [kmer, strand] = StrandKmerTable::DecodeEntry(encoded_pair);
        if (strand == kForwardStrandBit) {
          bf.Insert(kmer, false);
        } else if (strand == kReverseStrandBit) {
          bf.Insert(kmer, true);
        }
      }
    });
  }
  executor.run(taskflow).get();
}

/**
 * @brief Execute the strand build pipeline.
 * @param param StrandBuildParam object containing input parameters for the pipeline.
 * @details This function validates input parameters, reads the reference genome, builds the k-mer strand map,
 * counts unique k-mers, and constructs Bloom filters for forward and reverse strands.
 */
void StrandBuildPipeline(const StrandBuildParam& param) {
  // check if the input file exists
  if (!fs::exists(param.input)) {
    throw std::runtime_error("Input file does not exist: " + param.input.string());
  }

  // check thread number
  if (param.threads == 0) {
    throw std::runtime_error("Number of threads must be greater than 0");
  }

  // check that out prefix has a folder in the path, and if so create directories if needed
  if (!param.output.parent_path().empty()) {
    fs::create_directories(param.output.parent_path());
  }

  // Check k-mer size
  if (param.kmer_size > 31) {
    throw std::runtime_error("Max kmer size is 31");
  }

  Logging::Info("Reading fasta and index file: {}", param.input.string());

  // start reading the reference genome using fasta reader
  io::FastaReader reader(param.input);

  // determine how much space to reserve for the kmer strand map
  size_t total_bases = 0;
  for (const auto& record : reader) {
    total_bases += record.second;
  }
  StrandKmerTable kmer_strand_map(total_bases);

  Logging::Info("Building kmer strand hash table. Total bases: {}. Allocated size: {} MB", total_bases,
                (kmer_strand_map.Size() * 8) / 1024 / 1024);

  tf::Executor executor(param.threads);

  BuildStrandMap(reader, executor, kmer_strand_map, param.kmer_size);

  size_t num_unique_kmers_fw_total = 0;
  size_t num_unique_kmers_rv_total = 0;

  const auto& table = kmer_strand_map.GetTable();

  const size_t num_unique_kmers_both_total =
      CountStrandKmers(table, executor, num_unique_kmers_fw_total, num_unique_kmers_rv_total);

  Logging::Info("Unique k-mers fw: {} rv: {} both: {} ", num_unique_kmers_fw_total, num_unique_kmers_rv_total,
                num_unique_kmers_both_total);

  // get optimal number of hash functions relative to max FPR
  const auto hash_count = InterleavedBloomFilter::CalcOptimalNumHashes(param.max_false_positive_rate);

  // Initialize the bloom filters
  InterleavedBloomFilter filter(num_unique_kmers_fw_total, num_unique_kmers_rv_total, hash_count, param.kmer_size,
                                param.max_false_positive_rate);

  Logging::Info("Building Interleaved Bloom filter with {} hash functions of size {} mb", hash_count,
                filter.Size() / 8388608);

  BuildStrandFilter(table, filter, executor);

  const auto fprs = filter.FalsePositiveRates();

  // Sanity check by reporting the FPR
  Logging::Info("Interleaved Bloom filter false positive rates. Forward: {} Reverse: {}", fprs.first, fprs.second);

  filter.Save(param.output);
  Logging::Info("Saved bloom filters to {}", param.output.string());
}

}  // namespace xoos::demux::strand
