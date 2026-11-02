#pragma once
#include <algorithm>
#include <charconv>
#include <cstddef>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include <xoos/error/error.h>

#include "core/trimming-debug-mode.h"
#include "io/read-record.h"
#include "task/task.h"
#include "utility/alignment-util.h"

namespace xoos::demux {

constexpr char kOriginalIdPrefix = '@';
constexpr char kReadNameSeparator = ':';
constexpr char kMissingBarcodeChar = '*';
constexpr char kQualitySeparator = '+';
constexpr char kBitFlagSidSeparator = '|';
constexpr char kCommentSeparator = ' ';

// Constants for the --trimming-debug tags.
namespace debug {
constexpr std::string_view kSqTag = "sq:Z:";
constexpr std::string_view kSqSimplexInsertDelimiter = "<R>";
constexpr std::string_view kSqDuplexR1Delimiter = "<R1>";
constexpr std::string_view kSqDuplexR2Delimiter = "<R2>";
}  // namespace debug

// FASTQ header bitflag layout: a single byte encoded as a decimal integer between the read name and
// the SID, e.g. `@<name>:<bitflag>|<sid>`. A set bit means the corresponding feature is present.
// See the Everest read specification for the full 8-bit layout.
// Bit 7: read started in a threaded state.
inline constexpr auto kBitFlagStartsThreaded = static_cast<std::byte>(0b1000'0000);
// Bit 6: read ended in a threaded state.
inline constexpr auto kBitFlagEndsThreaded = static_cast<std::byte>(0b0100'0000);
// 5' SID-present bit (0b1000) in the FASTQ header bitflag.
inline constexpr auto kBitFlagSid5p = static_cast<std::byte>(0b1000);
// 3' SID-present bit (0b0100) in the FASTQ header bitflag.
inline constexpr auto kBitFlagSid3p = static_cast<std::byte>(0b0100);
// 5' UMI-present bit (0b0010) in the FASTQ header bitflag.
inline constexpr auto kBitFlagUmi5p = static_cast<std::byte>(0b0010);
// 3' UMI-present bit (0b0001) in the FASTQ header bitflag.
inline constexpr auto kBitFlagUmi3p = static_cast<std::byte>(0b0001);
// Convenience mask: duplex reads always carry both SIDs.
inline constexpr auto kBitFlagDuplexSids = kBitFlagSid5p | kBitFlagSid3p;
// Threaded-state bits (7 and 6), carried through unchanged from the RDB per-read state byte.
inline constexpr auto kBitFlagThreadedMask = kBitFlagStartsThreaded | kBitFlagEndsThreaded;

// Extract the threaded-state bits from the RDB per-read state byte for OR-ing into the output bitflag.
// The upstream read-state byte uses the same bit positions (7/6) as the output bitflag, so the
// threaded bits are copied through unchanged. Bits default to 0 when the source RDB predates the
// molecular_traces_read_state dataset.
constexpr std::byte ThreadedBitsFromReadState(const u8 read_state) {
  return static_cast<std::byte>(read_state) & kBitFlagThreadedMask;
}

// Maximum number of decimal digits required to render a u8 (the FASTQ header bitflag).
inline constexpr size_t kMaxU8DecimalDigits = std::numeric_limits<u8>::digits10 + 1;
// Maximum number of decimal digits required to render a u32 (SID, raw length, trim positions).
inline constexpr size_t kMaxU32DecimalDigits = std::numeric_limits<u32>::digits10 + 1;

/**
 * @brief Options controlling how simplex reads are formatted into FASTQ output.
 *
 * Bundled so callers pass a named-field struct instead of positional arguments — self-documenting at the
 * call site and extensible without touching every signature in the formatting call chain.
 */
struct SimplexFormatOptions {
  /// If true, override quality scores with `simplex_qual_score`. By default the original scores are preserved.
  bool override_qual_score = false;
  /// Quality score (phred+33) written for every base when `override_qual_score` is set.
  char simplex_qual_score = BaseQualityScores::kSimplexBaseQual;
  /// Active trimming debug modes. Empty means no debug tags.
  std::vector<TrimmingDebugMode> trimming_debug;
};

/**
 * @brief Writes demultiplexed data into appropriate output files.
 *
 * The Formatter class converts read records into various output formats for each adapter and writes them to disk.
 * It operates as a task within the taskflow scheduler framework.
 */
class Formatter : public Task {
 public:
  /**
   * @brief Constructs a Formatter task.
   * @param[in] exec The flow context providing execution environment.
   * @param[in] batch_nr The batch number to process.
   */
  Formatter(FlowContext& exec, size_t batch_nr);
  ~Formatter() override = default;

  /**
   * @brief Executes the formatter task.
   *
   * Called by the task scheduler to write a batch of records from the demux data into the appropriate output files.
   */
  void operator()();

  /**
   * @brief Converts a FixedReadRecord to processed FASTQ format for duplex adapter data. Outputs processed format
   * including YC tag, consensus sequence, etc.
   *
   * Exposed as static for projects that use demux like an API in another context.
   *
   * @param read The read record to convert.
   * @param p_data Output buffer for the formatted data.
   * @param offset Current position in buffer; updated after write.
   * @details Assumes the FixedReadRecord provides any additional metadata (for example a YC tag) via its
   * comment/header fields, and that ConsensusSeq() contains the duplex consensus sequence.
   */
  // Convenience overload used by unit tests; unused in non-test code.
  static void ToFastqDuplex(const FixedReadRecord& read, char* const p_data, size_t& offset) {  // NOSONAR(cpp:S5536)
    ToFastqDuplex(read, p_data, offset, {});
  }

  static void ToFastqDuplex(const FixedReadRecord& read,
                            char* p_data,
                            size_t& offset,
                            const std::vector<TrimmingDebugMode>& trimming_debug);

  /**
   * @brief Converts a FixedReadRecord to original raw FASTQ format. Outputs the original unprocessed read data. Usually
   * used for failed reads.
   *
   * Exposed as static for projects that use demux like an API in another context.
   *
   * @param[in] read The read record to convert.
   * @param[out] p_data Output buffer for the formatted data.
   * @param[in,out] offset Current position in buffer; updated after write.
   * @warning Assumes Seq(), Qual(), Name(), Comment() are unaltered, which may not hold true due to partial processing
   */
  static void ToFastqRaw(const FixedReadRecord& read, char* p_data, size_t& offset);

 private:
  /**
   * @brief Determine whether a record should be skipped during the initial scan for unwritten records.
   *
   * @param record The record to evaluate.
   * @return true if the record should be skipped (already processed or not eligible).
   */
  bool ShouldSkipRecord(const FixedReadRecord& record) const;

  /**
   * @brief Format a single record into the output buffer, dispatching on its status.
   *
   * Encapsulates the status-based switch logic previously inlined in operator().
   *
   * @param[in,out] record    The record to format; status is updated to kWritten on success.
   * @param         sink_sid  The output sink SID currently being collected.
   * @param[in,out] mem       The output buffer.
   * @return true if bytes were written to mem (caller should bounds-check), false if skipped.
   */
  bool FormatRecord(FixedReadRecord& record, u32 sink_sid, FormattedOutput& mem) const;

  /**
   * @brief Format a failed read for raw output, honoring --write-failed-reads and simplex quality overrides.
   *
   * Extracted from FormatRecord so the shared body of all failure states lives in one place.
   *
   * @param[in,out] record The failed record to format.
   * @param         sid    The SID currently being collected.
   * @param[in,out] mem    The output buffer.
   * @return true if bytes were written to mem, false if the read was skipped.
   */
  bool FormatFailedRecord(FixedReadRecord& record, u32 sid, FormattedOutput& mem) const;

  /**
   * @brief Schedule async write tasks for the collected sinks and mark the batch complete.
   *
   * @param sinks Span of sink descriptors populated during formatting.
   */
  void ScheduleWriteTasks(std::span<SinkData> sinks) const;

  /**
   * @brief Write a successfully-demultiplexed record to the buffer, dispatching on the configured adapter type.
   *
   * Reads the simplex formatting options and adapter type from the enclosing flow context, so they do not have to
   * be threaded through from `operator()`.
   */
  void PassingReadToBuffer(FixedReadRecord& read, FormattedOutput& mem) const;

  /**
   * @brief Write a failed record to the buffer as raw FASTQ.
   *
   * @warning Assumes data is unaltered, which may not hold true due to partial processing.
   */
  static void FailedReadToBuffer(FixedReadRecord& read, FormattedOutput& mem);
};

void WriteUintToBuffer(char* p_data, size_t& offset, size_t buffer_end, u32 value, std::string_view error_msg);
void WriteUmiValue(char* p_data, size_t& offset, size_t buffer_end, const std::optional<u32>& umi_value);
void WriteFastqHeaderWithBitflagAndSid(const FixedReadRecord& read, char* p_data, size_t& offset, std::byte bitflag);

/**
 * @brief Append a contiguous byte range to the output buffer and advance the offset.
 *
 * For trivially-copyable contiguous ranges the standard library lowers `std::copy_n` to
 * `__builtin_memmove`, which on modern libc shares fast paths with `memcpy`; functionally
 * identical here because the destination never overlaps the source. Caller is responsible
 * for ensuring the destination has sufficient capacity.
 */
inline void AppendBytes(char* const p_data, size_t& offset, const std::string_view src) {
  std::copy_n(src.data(), src.size(), p_data + offset);
  offset += src.size();
}

/// @brief Format raw sequence tag: "sq:Z:<5p_trimmed><R><3p_trimmed>"
std::string FormatSimplexRawSeqTag(const FixedReadRecord& read, const LociRange& insert);

// Format the duplex sq:Z tag from the saved raw sequence and hairpin range.
// Format: sq:Z:<5p_trimmed><R1><hairpin_trimmed><R2><3p_trimmed>
std::string FormatDuplexRawSeqTag(const FixedReadRecord& read);

void ToFastqDuplexUMI(const FixedReadRecord& read,
                      FormattedOutput& mem,
                      const std::vector<TrimmingDebugMode>& trimming_debug);

// Convenience overload used by unit tests; unused in non-test code.
inline void ToFastqDuplexUMI(const FixedReadRecord& read, FormattedOutput& mem) {  // NOSONAR(cpp:S5536)
  ToFastqDuplexUMI(read, mem, {});
}

/**
 * @brief Write a simplex read to the buffer without UMI information in the read name.
 *
 * @tparam TrimInfo The type containing trimming information (e.g., TrimInfoYsu, TrimInfoSimplex).
 * @param read The read to write.
 * @param mem The buffer to write to.
 * @param trim_info The trimming information for the read, used to extract the trimmed sequence if applicable.
 * @param opts Formatting options (quality override, trim debug tags).
 *
 * @details Output format: @<original id> <original comment>
 *          This function writes the read name without any UMI sequences, followed by the comment,
 *          sequence (trimmed if applicable), quality separator, and quality scores.
 *          Used for adapter types that don't capture UMI information (e.g., kYs).
 */
template <typename TrimInfo>
void ToFastqSimplex(const FixedReadRecord& read,
                    FormattedOutput& mem,
                    const TrimInfo& trim_info,
                    SimplexFormatOptions opts) {
  auto& offset = mem.nr_bytes;
  auto* const p_data = mem.p_data;

  // Write header: @<name>:<bitflag>|<sid>
  auto bitflag = GenerateSidBitFlagFromTrimInfo(trim_info);
  WriteFastqHeaderWithBitflagAndSid(read, p_data, offset, bitflag);

  // Write the body (comment, sequence, quality scores, etc.)
  AppendFastqSimplexBody(read, p_data, offset, trim_info, opts);
}

/// @brief Non-template implementation of the simplex body writer.
/// Called by the thin template wrapper with extracted TrimInfo fields.
void AppendFastqSimplexBodyImpl(const FixedReadRecord& read,
                                char* p_data,
                                size_t& offset,
                                const LociRange& insert,
                                bool is_trimmed,
                                const SimplexFormatOptions& opts);

/**
 * @brief Writes the read body (comment, sequence, and quality scores) to the buffer,
 *        assuming the read name has already been written.
 * Useful to reduce code duplication between simplex UMI and non-UMI output.
 */
template <typename TrimInfo>
void AppendFastqSimplexBody(const FixedReadRecord& read,
                            char* const p_data,
                            size_t& offset,
                            const TrimInfo& trim_info,
                            const SimplexFormatOptions& opts) {
  AppendFastqSimplexBodyImpl(read, p_data, offset, trim_info.insert, trim_info.sid.has_value(), opts);
}

/**
 * @brief Write the read to the buffer, using the TrimInfo to get UMI sequences and trimming information.
 * @param read The read to write.
 * @param mem The buffer to write to.
 * @param trim_info The trimming and UMI information for the read.
 * @param opts Formatting options (quality override, trim debug tags).
 */
template <typename TrimInfo>
void ToFastqSimplexUMI(const FixedReadRecord& read,
                       FormattedOutput& mem,
                       const TrimInfo& trim_info,
                       SimplexFormatOptions opts) {
  auto& offset = mem.nr_bytes;
  auto* const p_data = mem.p_data;

  // Write header: @<name>:<bitflag>|<sid>
  auto bitflag = GenerateUmiBitFlagFromTrimInfo(trim_info);
  WriteFastqHeaderWithBitflagAndSid(read, p_data, offset, bitflag);

  // UMI name portion, parsers look for a pipe '|' character currently
  p_data[offset++] = kReadNameSeparator;
  if (trim_info.umi_5p.has_value()) {
    auto [ptr, ec] = std::to_chars(p_data + offset, p_data + mem.capacity, *trim_info.umi_5p);
    offset = static_cast<size_t>(ptr - p_data);
    if (ec != std::errc()) {
      throw error::Error("UMI conversion error. Either buffer too small or value invalid.\n{}",
                         std::make_error_code(ec).message());
    }
  } else {
    p_data[offset++] = kMissingBarcodeChar;
  }
  p_data[offset++] = kReadNameSeparator;
  if (trim_info.umi_3p.has_value()) {
    auto [ptr, ec] = std::to_chars(p_data + offset, p_data + mem.capacity, *trim_info.umi_3p);
    offset = static_cast<size_t>(ptr - p_data);
    if (ec != std::errc()) {
      throw error::Error("UMI conversion error. Either buffer too small or value invalid.\n{}",
                         std::make_error_code(ec).message());
    }
  } else {
    p_data[offset++] = kMissingBarcodeChar;
  }

  // Write the body (comment, sequence, quality scores, etc.)
  AppendFastqSimplexBody(read, p_data, offset, trim_info, opts);
}

/**
 * @brief Generates a SID-only bitflag indicating which SIDs are present.
 * UMI bits are always zero; use GenerateUmiBitFlag when UMI presence matters.
 * @param has_sid_5p Whether the 5' SID is present.
 * @param has_sid_3p Whether the 3' SID is present.
 * @return Bitflag with 5' SID = 0b1000 and 3' SID = 0b0100.
 */
std::byte GenerateSidBitFlag(bool has_sid_5p, bool has_sid_3p);

/**
 * @brief Generates a bitflag indicating which SIDs and UMIs are present.
 * @param has_sid_5p Whether the 5' SID is present.
 * @param has_sid_3p Whether the 3' SID is present.
 * @param has_umi_5p Whether the 5' UMI is present.
 * @param has_umi_3p Whether the 3' UMI is present.
 * @return Bitflag with 5' SID = 0b1000, 3' SID = 0b0100, 5' UMI = 0b0010, and 3' UMI = 0b0001.
 */
std::byte GenerateUmiBitFlag(bool has_sid_5p, bool has_sid_3p, bool has_umi_5p, bool has_umi_3p);

/**
 * @brief Generates a SID-only bitflag indicating which SIDs are present.
 * UMI bits are always zero; use GenerateUmiBitFlagFromTrimInfo when UMI presence matters.
 * @tparam TrimInfo Any trim-info type with optional sid_5p / sid_3p fields
 *         (e.g., TrimInfoYsu, TrimInfoSimplex).
 * @param trim_info The trimming information to check.
 * @return Bitflag: 5' sid (0b1000), 3' sid (0b0100), 5' umi (0b0010), 3' umi (0b0001)
 */
template <typename TrimInfo>
std::byte GenerateSidBitFlagFromTrimInfo(const TrimInfo& trim_info) {
  return GenerateSidBitFlag(trim_info.sid_5p.has_value(), trim_info.sid_3p.has_value());
}

/**
 * @brief Generates a bitflag based on which SIDs and UMIs are present in the given TrimInfo.
 * @tparam TrimInfo Any trim-info type with optional sid_5p / sid_3p / umi_5p / umi_3p fields
 *         (e.g., TrimInfoYsu, TrimInfoDuplexUMI).
 * @param trim_info The trimming and UMI information to check.
 * @return Bitflag: 5' sid (0b1000), 3' sid (0b0100), 5' umi (0b0010), 3' umi (0b0001)
 */
template <typename TrimInfo>
std::byte GenerateUmiBitFlagFromTrimInfo(const TrimInfo& trim_info) {
  return GenerateUmiBitFlag(trim_info.sid_5p.has_value(),
                            trim_info.sid_3p.has_value(),
                            trim_info.umi_5p.has_value(),
                            trim_info.umi_3p.has_value());
}

std::byte GenerateUmiBitFlagFromFixedDuplexRead(const FixedReadRecord& read);

}  // namespace xoos::demux
