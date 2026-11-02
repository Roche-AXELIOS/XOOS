#pragma once

#include <array>
#include <chrono>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

#include <xoos/types/vec.h>

#include "adapters/duplex/umi/trim-info-duplex-umi.h"
#include "adapters/simplex/trim-info-simplex.h"
#include "adapters/ysu/trim-info-ysu.h"
#include "sequence/matcher/match-position-states.h"
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
 * Original ReadRecord struct which uses four strings to store the id, comment, sequence, and
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
// Used by edlib
constexpr uint kBlockSize = 64;
// Used by edlib, 21 allows for maximum length of 1344 bases for 5p end of duplex
constexpr uint kMaxNrBlocks = 21;
// 3008 (use weird multiplier to avoid caching issues)
constexpr uint kMaxReadLength{kBlockSize * kMaxNrBlocks * 2};
constexpr uint kBufferSize = (2 * kMaxReadLength) + kMaxCommentLength + kMaxIdLength;
constexpr uint kMaxReadLength4 = kMaxReadLength >> 2;
// alphabet size is 4, add 1, multiply by nr blocks
constexpr uint kPeqTableSize{kMaxNrBlocks * 5};
// Maximum length for alignment
constexpr uint kMaxAlignLength = kMaxReadLength >> 1;

/**
 * ReadRecord struct which uses a fixed-size memory buffer to store the id, comment, sequence, and consensus sequence
 * This buffer support reads up to about 2700 bases long (assuming a small size of comment and id). If we encounter
 * longer reads, we will discard them and log it. If log shows far too many long reads, we can increase the buffer size.
 */
// NOLINTNEXTLINE - this is a struct, so we don't need to initialize the members
// This is a deliberately fixed-buffer, zero-allocation record reused across a hot processing loop. Its inline
// C-arrays, large field count, and the single trailing private _status member are intentional design choices for
// cache locality; members are not initialized on construction because Clear()/Initialize() populate them before
// use. The related structural rules are scoped out for this file in spec.yaml.
struct FixedReadRecord {
  // This is used to indicate that the file SID was not found or is invalid.
  static constexpr u32 kUnassignedSID = 0;

  // Sentinel error metric used to mark a record with no assigned SID.
  static constexpr s32 kUnassignedErrorMetric = 51'360;

  char buf[kBufferSize];
  // Used to store the consensus sequence (duplex HD only)
  char consensus_buffer[kMaxAlignLength];
  // Length of the consensus sequence
  s32 consensus_seq_len;
  // Offset of the consensus sequence in the consensus buffer
  u32 consensus_seq_offset{0};
  // Length of the run-length encoded IUPAC string (which will overwrite the sequence)
  u32 iupac_length;
  // last position in the consensus sequence
  u8 two_bit_seq[kTwoBitPadding + kTwoBitPadding + kMaxReadLength4];
  // 2 bytes per position, but add padding to allow out-of-bounds access
  MatchScore match_values[kMaxReadLength + 64];

  // return pointer to the middle of the array
  uint8_t* TwoBitsSeq() {
    return two_bit_seq + kTwoBitPadding;
  }

  // return pointer to the middle of the array
  const uint8_t* TwoBitsSeq() const {
    return two_bit_seq + kTwoBitPadding;
  }

  // cannot use a union because of non-trivial constructors
  // NOTE: use c++17 std::variant instead of a union?
  //       Alternatively since this determined at the start of the program run and doesn't change could use a template?
  TrimInfoYsu trim_info_ysu;
  TrimInfoSimplex trim_info_simplex;
  // has umi information, even if we don't use it
  TrimInfoDuplexUMI trim_info_duplex;
  // Stable one-based SID identifier written in the FASTQ read name.
  u32 sid_index{kUnassignedSID};
  // Sink identifier used for output routing. In all-split mode, partial reads use a separate sink
  // from full reads while retaining the same sid_index.
  u32 sink_id{kUnassignedSID};

  // Offsets of the name, comment, sequence, and quality score in the buffers. FastQ is organized in that
  // way, so we assume that the buffer is contiguous and we can use offsets to access the data.
  uint name_offset;
  uint comment_offset;
  uint seq_offset;
  uint qual_offset;
  uint end_offset;

  uint SeqLen() const {
    return qual_offset - seq_offset;
  }

  uint QualLen() const {
    return end_offset - qual_offset;
  }

  uint NameLen() const {
    return comment_offset - name_offset;
  }

  uint CommentLen() const {
    return seq_offset - comment_offset;
  }

  char* Seq() {
    return buf + seq_offset;
  }

  const char* Seq() const {
    return buf + seq_offset;
  }

  char* Qual() {
    return buf + qual_offset;
  }

  const char* Qual() const {
    return buf + qual_offset;
  }

  char* Name() {
    return buf + name_offset;
  }

  const char* Name() const {
    return buf + name_offset;
  }

  char* Comment() {
    return buf + comment_offset;
  }

  const char* Comment() const {
    return buf + comment_offset;
  }

  char* ConsensusSeq() {
    return consensus_buffer + consensus_seq_offset;
  }

  const char* ConsensusSeq() const {
    return consensus_buffer + consensus_seq_offset;
  }

  // Used to get the maximum read length that can be stored in the buffer.
  uint MaximumReadLength() const {
    return (kBufferSize - seq_offset) >> 1;
  }

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
    kFailedMidadapterTrimFail,
    // The 5' and 3' SIDs disagree and the discordant-sid-mode requires discard
    kDiscordantSidFail,
    // The read is partial (single-ended SID) and full-only mode requires discard
    kFullOnlyFail,
    // A UMI-simplex (YSU) read was assigned by SID but lacks a SID+UMI combination on any end
    kUmiFail
  };

  /**
   * @brief Sets the status of the FixedReadRecord and updates the SID fields if necessary.
   *
   * This function modifies the status of the FixedReadRecord and performs additional actions
   * based on the new status. If the new status indicates that the record is invalid (e.g., too long,
   * too short, or failed due to duplex-specific issues), both SID fields are reset to 0. For valid statuses,
   * the SID fields remain unchanged.
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
      case kDiscordantSidFail:
      case kFullOnlyFail:
      case kUmiFail:
        // No longer a valid read, so clear the output header and routing SID fields.
        sid_index = kUnassignedSID;
        sink_id = kUnassignedSID;
        break;
      case kRead:
      case kWritten:
      case kDemultiplexed:
        // A valid read, so we can keep the SID fields.
        break;
      default:
        // Invalid status
        throw std::invalid_argument("Invalid status for FixedReadRecord");
        return;
    }
    _status = new_status;
  }

  Status GetStatus() const {
    return _status;
  }

  // Start position (inclusive) of the duplex loop sequence in the raw read, or -1 if not found
  s32 loop_start_pos;
  // End position (inclusive) of the duplex loop sequence in the raw read, or -1 if not found
  s32 loop_end_pos;
  s32 error_metric;

  // Raw sequence snapshot for duplex reads, captured before alignment overwrites Seq().
  // Only populated when the kRawSeq trimming debug mode is active.
  std::string duplex_raw_seq;

  // Trimmed R1 (5' side) and R2 (3' side) insert lengths from duplex alignment.
  // Set after GenerateConsensusQualAndIUPAC; used by FormatDuplexRawSeqTag to
  // derive adapter boundaries relative to the midadapter position.
  u32 duplex_r1_len{0};
  u32 duplex_r2_len{0};

  // Per-read state bitfield carried through from the RDB `molecular_traces_read_state` dataset.
  // Holds the upstream threaded-state bits; propagated into the output read bitflag by the
  // formatter. Defaults to 0 when the input RDB predates the dataset.
  u8 read_state{0};

  // Called to clear the record before reusing it.
  void Clear() {
    name_offset = 0;
    comment_offset = 0;
    seq_offset = 0;
    qual_offset = 0;
    end_offset = 0;
    trim_info_ysu.Clear();
    trim_info_simplex.Clear();
    trim_info_duplex.Clear();
    SetStatus(Status::kNotRead);
    loop_start_pos = kNoMatchPosition;
    loop_end_pos = kNoMatchPosition;
    // sentinel value used when no SID has been assigned
    error_metric = kUnassignedErrorMetric;
    duplex_raw_seq.clear();
    duplex_r1_len = 0;
    duplex_r2_len = 0;
    read_state = 0;
  }

  // Test-support helper for constructing a record from string fields; used by the unit tests, which
  // SonarQube does not scan, so it appears unused from src/ alone.
  void Initialize(const std::string_view name,
                  const std::string_view comment,
                  const std::string_view seq,
                  const std::string_view qual) {
    Clear();
    // Copy the name, comment, sequence, and quality score to the buffer.
    name_offset = 0;
    comment_offset = name_offset + static_cast<uint>(name.copy(buf + name_offset, kMaxIdLength));
    seq_offset = comment_offset + static_cast<uint>(comment.copy(buf + comment_offset, kMaxCommentLength));
    // seq and qual hold base/quality characters (one byte each), so they are bounded by kMaxReadLength, not by the
    // 2-bit-packed buffer size kMaxReadLength4. kBufferSize reserves kMaxReadLength bytes for each; using the smaller
    // limit here would truncate long reads and diverge from the production fastq/rdb readers.
    qual_offset = seq_offset + static_cast<uint>(seq.copy(buf + seq_offset, kMaxReadLength));
    end_offset = qual_offset + static_cast<uint>(qual.copy(buf + qual_offset, kMaxReadLength));
    InitTwoBitSeq();
    SetStatus(Status::kRead);
  }

  // Low-level helper that bridges to the SIMD 2-bit encoder, which operates on a raw (pointer, length) pair;
  // a std::string_view cannot be used because ConvertTo2Bit needs a raw byte pointer and explicit length
  // semantics for out-of-bounds padding. The related rules are scoped out for this file in spec.yaml.
  void InitTwoBit(const u8* seq, u32 seq_len) {
    simd::ConvertTo2Bit(seq, 0, seq_len, TwoBitsSeq());
    // fill padding regions with 0s
    std::memset(TwoBitsSeq() - kTwoBitPadding, 0, kTwoBitPadding);
    const uint num_bytes = (seq_len + 3) / 4;
    std::memset(TwoBitsSeq() + num_bytes, 0, kTwoBitPadding);
  }

  void InitTwoBitSeq() {
    // NOSONAR(S5356) -- reinterpret_cast from char* to u8* is required to feed the raw sequence bytes to the
    // SIMD 2-bit encoder; both are byte-sized types and the reinterpretation is well defined.
    InitTwoBit(reinterpret_cast<u8*>(Seq()), SeqLen());  // NOSONAR(S5356)
  }

 private:
  // Status of the record, default is not read
  Status _status{Status::kNotRead};
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
  // NOTE: As requirements change (more tags) we should perhaps add some checks to make sure the above is true
  // The formatted records (all data in one contiguous block)
  std::shared_ptr<vec<char>> formatted_records;

  // Convenience constructor used by unit tests; non-test code uses the two-arg form.
  explicit FormattedOutput(const uint64_t batch_size)  // NOSONAR(cpp:S5536)
      : FormattedOutput(batch_size, kBufferSize) {
  }

  FormattedOutput(const uint64_t batch_size, const uint64_t record_budget) {
    formatted_records = std::make_shared<vec<char>>(batch_size * record_budget);
    p_data = formatted_records->data();
    capacity = formatted_records->size();
  }

  // Number of bytes of valid formatted data currently stored in the buffer.
  size_t nr_bytes{0};
  // number of bytes allocated in the buffer
  size_t capacity{0};
  char* p_data{nullptr};
  static constexpr u32 kMaxNumberSids = 256;
  // Two sinks per SID for full/partial split mode plus the reserved failed-read sink.
  static constexpr u32 kMaxNumberSinks = (2 * kMaxNumberSids) + 1;
  // This state information is used by the async writer tasks.
  std::array<SinkData, kMaxNumberSinks> file_sinks;
};

struct FixedReadRecordBatch {
  // We will reuse these batches and once allocated, we do not want to resize the batch. Hence, the allocated
  // size is fixed and we need to keep track of the number of records in the batch.
  std::shared_ptr<vec<FixedReadRecord>> records;
  FormattedOutput formatted_output;
  uint64_t num_records{0UL};
  uint64_t max_num_records{0UL};
  uint64_t num_bytes{0UL};
  uint64_t num_bases{0UL};
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
  std::chrono::time_point<std::chrono::high_resolution_clock> end_time;

  // Convenience constructor used by unit tests; non-test code uses the two-arg form.
  explicit FixedReadRecordBatch(const uint64_t size)  // NOSONAR(cpp:S5536)
      : FixedReadRecordBatch(size, kBufferSize) {
  }

  FixedReadRecordBatch(const uint64_t size, const uint64_t record_budget)
      : formatted_output(size, record_budget), max_num_records(size) {
    records = std::make_shared<vec<FixedReadRecord>>(size);
  }

  // These objects will never be copied.
  FixedReadRecordBatch(const FixedReadRecordBatch&) = delete;
  FixedReadRecordBatch& operator=(const FixedReadRecordBatch&) = delete;

  uint64_t Capacity() const {
    return max_num_records;
  }

  uint64_t Size() const {
    return num_records;
  }
};

}  // namespace xoos::demux
