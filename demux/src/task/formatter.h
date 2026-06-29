#pragma once
#include <xoos/error/error.h>

#include <charconv>
#include <cstddef>
#include <system_error>

#include "io/read-record.h"
#include "task/task.h"

namespace xoos::demux {

constexpr char kOriginalIdPrefix = '@';
constexpr char kReadNameSeparator = ':';
constexpr char kMissingBarcodeChar = '*';
constexpr char kQualitySeparator = '+';
constexpr char kBitFlagSidSeparator = '|';
constexpr char kCommentSeparator = ' ';

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
  static void ToFastqDuplex(const FixedReadRecord& read, char* p_data, size_t& offset);

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
};

void WriteUintToBuffer(char* p_data, size_t& offset, size_t buffer_end, u32 value, std::string_view error_msg);
void WriteUmiValue(char* p_data, size_t& offset, size_t buffer_end, const std::optional<u32>& umi_value);
void WriteFastqHeaderWithBitflagAndSid(const FixedReadRecord& read, char* p_data, size_t& offset, u8 bitflag);

void ToFastqDuplexUMI(const FixedReadRecord& read, FormattedOutput& mem);

/**
 * @brief Write a simplex read to the buffer without UMI information in the read name.
 *
 * @tparam TrimInfo The type containing trimming information (e.g., TrimInfoYs, TrimInfoYsuSl, TrimInfoYsuTe).
 * @param read The read to write.
 * @param mem The buffer to write to.
 * @param trim_info The trimming information for the read, used to extract the trimmed sequence if applicable.
 *
 * @details Output format: @<original id> <original comment>
 *          This function writes the read name without any UMI sequences, followed by the comment,
 *          sequence (trimmed if applicable), quality separator, and quality scores.
 *          Used for adapter types that don't capture UMI information (e.g., kYs).
 */
template <typename TrimInfo>
void ToFastqSimplex(const FixedReadRecord& read, FormattedOutput& mem, const TrimInfo& trim_info) {
  // Copy over comment prefixed with '@' character.
  auto& offset = mem.nr_bytes;
  auto* const p_data = mem.p_data;

  // Write header: @<name>:<bitflag>|<sid>
  auto bitflag = GenerateYsuSlBitFlagFromTrimInfo(trim_info);
  WriteFastqHeaderWithBitflagAndSid(read, p_data, offset, bitflag);

  // Write the body (comment, sequence, quality scores, etc.)
  AppendFastqSimplexBody(read, p_data, offset, trim_info);
}

/**
 * @brief Writes the read body (comment, sequence, and quality scores) to the buffer,
 *        assuming the read name has already been written.
 * Useful to reduce code duplication between simplex UMI and non-UMI output.
 *
 * @param read The read to write.
 * @param p_data Output buffer for the formatted data.
 * @param offset Current position in buffer; updated after write.
 * @param trim_info The trimming and UMI information for the read.
 */
template <typename TrimInfo>
void AppendFastqSimplexBody(const FixedReadRecord& read, char* const p_data, size_t& offset,
                            const TrimInfo& trim_info) {
  // Add space between comment and name
  p_data[offset++] = ' ';

  // Copy the comment as-is - and don't forget to add the '\n' at the end.
  std::memcpy(p_data + offset, read.Comment(), read.CommentLen());
  offset += read.CommentLen();
  // Add a newline after the comment.
  p_data[offset++] = '\n';

  // Copy the sequence.
  if (trim_info.sid) {
    // If the sid was trimmed, we need to copy the trimmed sequence
    const auto length = trim_info.insert.Length();
    std::memcpy(p_data + offset, read.Seq() + trim_info.insert.spos, length);
    offset += length;
  } else {
    // copy untrimmed sequence.
    std::memcpy(p_data + offset, read.Seq(), read.SeqLen());
    offset += read.SeqLen();
  }
  p_data[offset++] = '\n';
  p_data[offset++] = kQualitySeparator;
  p_data[offset++] = '\n';

  // Copy the quality scores
  if (trim_info.sid) {
    // If the sid was trimmed, we need to copy the trimmed quality scores
    const auto length = trim_info.insert.Length();
    std::memcpy(p_data + offset, read.Qual() + trim_info.insert.spos, length);
    offset += length;
  } else {
    // copy untrimmed quality scores.
    std::memcpy(p_data + offset, read.Qual(), read.QualLen());
    offset += read.QualLen();
  }
  // Add a newline at the end of the quality scores.
  p_data[offset++] = '\n';
}

/**
 * @brief Write the read to the buffer, using the TrimInfo to get UMI sequences and trimming information.
 * @param read The read to write.
 * @param mem The buffer to write to.
 * @param trim_info The trimming and UMI information for the read.
 * @details Output format: @<original id>|7|3 <original suffix>
 * '*' is used to indicate missing barcodes.
 */
template <typename TrimInfo>
void ToFastqSimplexUMI(const FixedReadRecord& read, FormattedOutput& mem, const TrimInfo& trim_info) {
  // Copy over comment prefixed with '@' character.
  auto& offset = mem.nr_bytes;
  auto* const p_data = mem.p_data;

  // Write header: @<name>:<bitflag>|<sid>
  auto bitflag = GenerateUmiBitFlagFromTrimInfo(trim_info);
  WriteFastqHeaderWithBitflagAndSid(read, p_data, offset, bitflag);

  // UMI name portion, parsers look for a pipe '|' character currently
  p_data[offset++] = kReadNameSeparator;
  if (trim_info.umi_5p.has_value()) {
    auto [ptr, ec] = std::to_chars(p_data + offset, p_data + mem.capacity, trim_info.umi_5p.value());
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
    auto [ptr, ec] = std::to_chars(p_data + offset, p_data + mem.capacity, trim_info.umi_3p.value());
    offset = static_cast<size_t>(ptr - p_data);
    if (ec != std::errc()) {
      throw error::Error("UMI conversion error. Either buffer too small or value invalid.\n{}",
                         std::make_error_code(ec).message());
    }
  } else {
    p_data[offset++] = kMissingBarcodeChar;
  }

  // Write the body (comment, sequence, quality scores, etc.)
  AppendFastqSimplexBody(read, p_data, offset, trim_info);
}

/**
 * @brief Generates a bitflag based on which SIDs are present
 * The bitflag is used to indicate the presence of 5' and 3' SIDs and UMIs in the output FASTQ header.
 * In this case no UMI information is captured in the read name, so the bitflag only reflects which SIDs are present.
 * For duplex without a UMI this is a constant as both SIDs are present, but here we have to look for both SID values.
 * @tparam TrimInfo The type of trim info (TrimInfoYsuSl, TrimInfoYsuTe, or TrimInfoDuplexUMI)
 * @param trim_info The trimming and UMI information to check.
 * @return Bitflag: 5' sid (0b1000), 3' sid (0b0100), 5' umi (0b0010), 3' umi (0b0001)
 */
template <typename TrimInfo>
u8 GenerateYsuSlBitFlagFromTrimInfo(const TrimInfo& trim_info) {
  // 5' sid: 0b1000
  // 3' sid: 0b0100
  // 5' umi: 0b0010
  // 3' umi: 0b0001
  std::byte bitflag{0};

  if (trim_info.sid_5p.has_value()) {
    bitflag |= std::byte{0b1000};
  }
  if (trim_info.sid_3p.has_value()) {
    bitflag |= std::byte{0b0100};
  }
  return std::to_integer<u8>(bitflag);
}

/**
 * @brief Generates a bitflag based on which SIDs and UMIs are present in the given TrimInfo.
 * The bitflag is used to indicate the presence of 5' and 3' SIDs and UMIs in the output FASTQ header.
 * @tparam TrimInfo The type of trim info (TrimInfoYsuSl, TrimInfoYsuTe, or TrimInfoDuplexUMI)
 * @param trim_info The trimming and UMI information to check.
 * @return Bitflag: 5' sid (0b1000), 3' sid (0b0100), 5' umi (0b0010), 3' umi (0b0001)
 */
template <typename TrimInfo>
u8 GenerateUmiBitFlagFromTrimInfo(const TrimInfo& trim_info) {
  // 5' sid: 0b1000
  // 3' sid: 0b0100
  // 5' umi: 0b0010
  // 3' umi: 0b0001
  std::byte bitflag{0};

  if (trim_info.sid_5p.has_value()) {
    bitflag |= std::byte{0b1000};
  }
  if (trim_info.sid_3p.has_value()) {
    bitflag |= std::byte{0b0100};
  }
  if (trim_info.umi_5p.has_value()) {
    bitflag |= std::byte{0b0010};
  }
  if (trim_info.umi_3p.has_value()) {
    bitflag |= std::byte{0b0001};
  }
  return std::to_integer<u8>(bitflag);
}

u8 GenerateUmiBitFlagFromFixedDuplexRead(const FixedReadRecord& read);

}  // namespace xoos::demux
