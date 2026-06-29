#pragma once

#include <chrono>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "adapters/duplex/umi/trim-info-duplex-umi.h"
#include "adapters/ys/trim-info-ys.h"
#include "adapters/ysu-sl/trim-info-ysu-sl.h"
#include "adapters/ysu-te/trim-info-ysu-te.h"
#include "simd/simd-functions.h"
#include "xoos/types/int.h"

namespace xoos::demux {

// Match scores depict two values: a "raw" score (max 23) that only allows for substitutions and a Shifted Hamming
// Distance score that also allows for insertions and deletions.
struct MatchScore {
  uint8_t raw_score;
  uint8_t shd_score;
};

/**
 * KZ2024 - Original ReadRecord struct which uses four strings to store the id, comment, sequence, and
 * quality score of a read. Key issue: allocating and copying strings is expensive and shows up in the profiler.
 */
struct ReadRecord {
  std::string id;
  std::string comment;
  std::string seq;
  std::string qual;
};

// As we preallocate memory, we need to know the maximum length of a read. As we are doing short reads, I will
// use 2K as the maximum read length, which means that we allocate about 10KB per read record. This also means one
// needs to be careful with batch sizes; we'll need 1 MB of memory for a batch of 100 reads, I recommend a batch size
// of about 10,000 reads (100 MB per thread).

constexpr uint kMaxCommentLength = 128;
constexpr uint kMaxIdLength = 256;
constexpr uint kTwoBitPadding = 64;
constexpr uint kBlockSize = 64;    // Used by edlib
constexpr uint kMaxNrBlocks = 21;  // Used by edlib, 21 allows for maximum length of 1344 bases for 5p end of duplex
constexpr uint kMaxReadLength{kBlockSize * kMaxNrBlocks * 2};  // 3008 (use weird multiplier to avoid caching issues)
constexpr uint kBufferSize = 2 * kMaxReadLength + kMaxCommentLength + kMaxIdLength;
constexpr uint kMaxReadLength4 = kMaxReadLength >> 2;
constexpr uint kPeqTableSize{kMaxNrBlocks * 5};        // alphabet size is 4, add 1, multiply by nr blocks
constexpr uint kMaxAlignLength = kMaxReadLength >> 1;  // Maximum length for alignment

/**
 * ReadRecord struct which uses a fixed-size memory buffer to store the id, comment, sequence, and consensus sequence
 * This buffer support reads up to about 2700 bases long (assuming a small size of comment and id). If we encounter
 * longer reads, we will discard them and log it. If log shows far too many long reads, we can increase the buffer size.
 */
struct FixedReadRecord {                     // NOLINT - this is a struct, so we don't need to initialize the members
  constexpr static uint kUnassignedSID = 0;  // This is used to indicate that the file SID was not found or is invalid.

  char buf[kBufferSize];
  char consensus_buffer[kMaxAlignLength];  // Used to store the consensus sequence (duplex HD only)
  s32 consensus_seq_len;                   // Length of the consensus sequence
  u32 consensus_seq_offset{0};             // Offset of the consensus sequence in the consensus buffer
  u32 iupac_length;  // Length of the run-length encoded IUPAC string (which will overwrite the sequence)
  // last position in the consensus sequence
  u8 two_bit_seq[kTwoBitPadding + kTwoBitPadding + kMaxReadLength4];
  MatchScore match_values[kMaxReadLength + 64];  // 2 bytes per position, but add padding to allow out-of-bounds access

  uint8_t* TwoBitsSeq() { return two_bit_seq + kTwoBitPadding; }  // return pointer to the middle of the array

  const uint8_t* TwoBitsSeq() const {
    return two_bit_seq + kTwoBitPadding;
  }  // return pointer to the middle of the array

  // cannot use a union because of non-trivial constructors
  // TODO use c++17 std::variant instead of a union?
  //      Alternatively since this determined at the start of the program run and doesn't change could use a template?
  TrimInfoYsuSl trim_info_ysu_sl;
  TrimInfoYsuTe trim_info_ysu_te;
  TrimInfoYs trim_info_ys;
  // has umi information, even if we don't use it
  TrimInfoDuplexUMI trim_info_duplex;
  uint file_sid{kUnassignedSID};  // Store sid found separately

  // Offsets of the name, comment, sequence, and quality score in the buffers. FastQ is organized in that
  // way, so we assume that the buffer is contiguous and we can use offsets to access the data.
  uint name_offset;
  uint comment_offset;
  uint seq_offset;
  uint qual_offset;
  uint end_offset;

  uint SeqLen() const { return qual_offset - seq_offset; }

  uint QualLen() const { return end_offset - qual_offset; }

  uint NameLen() const { return comment_offset - name_offset; }

  uint CommentLen() const { return seq_offset - comment_offset; }

  char* Seq() { return buf + seq_offset; }

  char* Qual() { return buf + qual_offset; }

  char* Name() { return buf + name_offset; }

  char* Comment() { return buf + comment_offset; }

  char* ConsensusSeq() { return consensus_buffer + consensus_seq_offset; }

  const char* Seq() const { return buf + seq_offset; }

  const char* Qual() const { return buf + qual_offset; }

  const char* Name() const { return buf + name_offset; }

  const char* Comment() const { return buf + comment_offset; }

  const char* ConsensusSeq() const { return consensus_buffer + consensus_seq_offset; }

  // Used to get the maximum read length that can be stored in the buffer.
  uint MaximumReadLength() const { return (kBufferSize - seq_offset) >> 1; }

  // Status of the record, used to track the status of the record in the pipeline.
  enum class Status {
    // Valid states:
    // The record has not been read yet
    kNotRead,
    // The record has been read
    kRead,
    // The record has been demultiplexed (i.e. the SID has been found)
    kDemultiplexed,
    // The record has been written out
    kWritten,

    // Fail states (before assignment):
    // The raw record was too long and due to limited maximum read length
    kTooLongFail,
    // The raw read length is below the minimum read length
    kTooShortFail,
    // We couldn't find the midadapter for some reason
    kDuplexMidAdapterFail,

    // Fail states (after assignment):
    // Consensus errors exceeds edit distance threshold
    kDuplexEditDistanceFail,
    // The consensus sequence is too long for buffer during construction
    kDuplexTooLongFail,
    // The trimmed read length is below the minimum read length
    kTrimmedTooShortFail,
    // There was extra expected loop sequence but we failed to trim it
    kFailedMidadapterTrimFail
  };

  /**
   * @brief Sets the status of the FixedReadRecord and updates the file_sid if necessary.
   *
   * This function modifies the status of the FixedReadRecord and performs additional actions
   * based on the new status. If the new status indicates that the record is invalid (e.g., too long,
   * too short, or failed due to duplex-specific issues), the file_sid is reset to 0. For valid statuses,
   * the file_sid remains unchanged.
   *
   * Try to use this function to set the status of the FixedReadRecord instead of directly modifying the status member.
   *
   * @param new_status The new status to set for the FixedReadRecord.
   * @throws std::invalid_argument If the provided status is invalid.
   */
  void SetStatus(const Status new_status) {
    switch (new_status) {
      using enum Status;
      case kNotRead:
      case kTooLongFail:
      case kTooShortFail:
      case kTrimmedTooShortFail:
      case kDuplexMidAdapterFail:
      case kDuplexEditDistanceFail:
      case kDuplexTooLongFail:
      case kFailedMidadapterTrimFail:
        // no longer a valid read so set the file_sid failed state
        file_sid = kUnassignedSID;
        break;
      case kRead:
      case kWritten:
      case kDemultiplexed:
        // A valid read, so we can keep the file_sid
        break;
      default:
        // Invalid status
        throw std::invalid_argument("Invalid status for FixedReadRecord");
        return;
    }
    _status = new_status;
  }

  Status GetStatus() const { return _status; }

  int hairpin_pos;
  int error_metric;

  // Called to clear the record before reusing it.
  inline void Clear() {
    name_offset = comment_offset = seq_offset = qual_offset = end_offset = 0;
    trim_info_ysu_sl.Clear();
    trim_info_ysu_te.Clear();
    trim_info_duplex.Clear();
    SetStatus(Status::kNotRead);
    hairpin_pos = -1;
    error_metric = 51360;  // why not (no sid found)
  }

  // To make conversion easier from old-style ReadRecord
  void Initialize(const std::string& name, const std::string& comment, const std::string& seq,
                  const std::string& qual) {
    Clear();
    // Copy the name, comment, sequence, and quality score to the buffer.
    name_offset = 0;
    comment_offset = name_offset + name.copy(buf + name_offset, kMaxIdLength);
    seq_offset = comment_offset + comment.copy(buf + comment_offset, kMaxCommentLength);
    qual_offset = seq_offset + seq.copy(buf + seq_offset, kMaxReadLength4);
    end_offset = qual_offset + qual.copy(buf + qual_offset, kMaxReadLength4);
    InitTwoBitSeq();
    SetStatus(Status::kRead);
  }

  void InitTwoBit(const u8* seq, u32 seq_len) {
    simd::ConvertTo2Bit(seq, 0, seq_len, TwoBitsSeq());
    // fill padding regions with 0s
    std::memset(TwoBitsSeq() - kTwoBitPadding, 0, kTwoBitPadding);
    uint num_bytes = (seq_len + 3) / 4;
    std::memset(TwoBitsSeq() + num_bytes, 0, kTwoBitPadding);
  }

  void InitTwoBitSeq() { InitTwoBit(reinterpret_cast<u8*>(Seq()), SeqLen()); }

 private:
  Status _status{Status::kNotRead};  // Status of the record, default is not read
};

// This structure is used to store the state information for a sink task. I opted to make this data part of the
// formatted output data rather than a task because the corresponding asynchronous task requires this data, but
// will not be visible as part of a task.
class FlowContext;

struct SinkData {
  // need to correctly write out the data.
  FlowContext* flow_context{nullptr};
  size_t batch_id{0};
  size_t sid_id{0};
  char* p_data{nullptr};
  size_t length{0};

  // Function implemented in flow-context.cpp.
  void operator()() const;
};

// Container for the formatted representation of a batch of reads after demuxing and trimming.
// Holds a single contiguous buffer of formatted records plus per-sink state used by asynchronous writers.
struct FormattedOutput {
  // Worst-case, formatted records are a little bit longer due to things like tags and UMI information
  // In practice this shouldn't happen due to trimming
  // TODO: As requirements change (more tags) we should perhaps add some checks to make sure the above is true

  // The formatted records (all data in one contiguous block)
  std::shared_ptr<std::vector<char>> formatted_records{};

  explicit FormattedOutput(uint64_t batch_size) {
    formatted_records = std::make_shared<std::vector<char>>(batch_size * kBufferSize);
    p_data = formatted_records->data();
    capacity = formatted_records->size();
  }

  // Number of bytes of valid formatted data currently stored in the buffer.
  size_t nr_bytes{0};
  // number of bytes allocated in the buffer
  size_t capacity{0};
  char* p_data{nullptr};
  // Maximum number of sinks that we can have in a batch
  static constexpr int kMaxNumberSinks = 200;
  // This state information is used by the async writer tasks.
  SinkData file_sinks[kMaxNumberSinks];
};

struct FixedReadRecordBatch {
  // We will reuse these batches and once allocated, we do not want to resize the batch. Hence, the allocated
  // size is fixed and we need to keep track of the number of records in the batch.
  std::shared_ptr<std::vector<FixedReadRecord>> records{};
  FormattedOutput formatted_output;
  uint64_t num_records{0UL};
  uint64_t max_num_records{0UL};
  uint64_t num_bytes{0UL};
  uint64_t num_bases{0UL};
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
  std::chrono::time_point<std::chrono::high_resolution_clock> end_time;

  // These objects will never be copied.
  FixedReadRecordBatch(const FixedReadRecordBatch&) = delete;
  FixedReadRecordBatch& operator=(const FixedReadRecordBatch&) = delete;

  uint64_t Capacity() const { return max_num_records; }

  uint64_t Size() const { return num_records; }

  explicit FixedReadRecordBatch(uint64_t size) : formatted_output(size), max_num_records(size) {
    records = std::make_shared<std::vector<FixedReadRecord>>(size);
  }
};

}  // namespace xoos::demux
