#include "task/formatter.h"

#include <utility>

#include <fmt/format.h>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include "task/flow-context.h"
#include "task/flow-manager.h"
#include "utility/alignment-util.h"
#include "utility/stop-watch.h"

namespace xoos::demux {
Formatter::Formatter(FlowContext& exec, size_t batch_nr) : Task(exec, batch_nr, fmt::format("Writer{}", batch_nr)) {
}

/**
 * @brief Generates a bitflag for the read based on which UMIs are present in the TrimInfoDuplexUMI. The bitflag is used
 * to indicate the presence of 5' and 3' SIDs and UMIs in the output FASTQ header for duplex reads.
 * @param[in] read The read to generate the bitflag for. We will look at the TrimInfoDuplexUMI to determine which UMIs
 * were present and set the appropriate bits in the bitflag.
 * @return
 */
std::byte GenerateUmiBitFlagFromFixedDuplexRead(const FixedReadRecord& read) {
  auto bitflag = kBitFlagDuplexSids;
  if (read.trim_info_duplex.umi_5p.has_value()) {
    bitflag |= kBitFlagUmi5p;
  }
  if (read.trim_info_duplex.umi_3p.has_value()) {
    bitflag |= kBitFlagUmi3p;
  }
  return bitflag;
}

/**
 * @brief Writes the read body (YC tag, comment, consensus sequence, and quality scores) to the buffer,
 *        assuming the read name has already been written.
 * Useful to reduce code duplication between UMI and non-UMI output.
 *
 * @param read The read to write.
 * @param p_data Output buffer for the formatted data.
 * @param offset Current position in buffer; updated after write.
 */
static void AppendFastqDuplexBody(const FixedReadRecord& read,
                                  char* const p_data,
                                  size_t& offset,
                                  const vec<TrimmingDebugMode>& trimming_debug) {
  // Add the YC tag (starts with space, e.g. " YC:Z:...")
  AppendBytes(p_data, offset, {read.Seq(), read.iupac_length});

  // Append sq:Z tag if raw-seq debug mode is active.
  if (HasTrimmingDebugMode(trimming_debug, TrimmingDebugMode::kRawSeq)) [[unlikely]] {  // NOLINT(readability/braces)
    const auto tag = FormatDuplexRawSeqTag(read);
    if (tag.size() > debug::kSqTag.size()) {
      p_data[offset++] = '\t';
      AppendBytes(p_data, offset, {tag.data(), tag.size()});
    }
  }

  // Copy the comment as-is.
  AppendBytes(p_data, offset, {read.Comment(), read.CommentLen()});
  // Add a newline after the comment.
  p_data[offset++] = '\n';

  const auto consensus_len = static_cast<size_t>(read.consensus_seq_len);

  // The sequence is now a consensus sequence.
  AppendBytes(p_data, offset, {read.ConsensusSeq(), consensus_len});

  // Add a newline at the end of the sequence.
  p_data[offset++] = '\n';
  // Add the '+' character after the sequence.
  p_data[offset++] = kQualitySeparator;
  // Add a newline after the '+' character.
  p_data[offset++] = '\n';

  // Copy the quality scores; these artificial scores have replaced the sequence.
  AppendBytes(p_data, offset, {read.Qual(), consensus_len});

  // Add a newline at the end of the quality scores.
  p_data[offset++] = '\n';
}

/**
 * @brief Write the read to the buffer, using the DuplexUmi bundle to get the UMI sequences.
 * @param read The read to write.
 * @param mem The buffer to write to.
 * @details Implementation is similar to the Duplex version, but we add the UMI sequences to the read name.
 * '*' is used to indicate missing barcodes.
 */
void ToFastqDuplexUMI(const FixedReadRecord& read, FormattedOutput& mem, const vec<TrimmingDebugMode>& trimming_debug) {
  auto& offset = mem.nr_bytes;
  auto* const p_data = mem.p_data;

  // Write header: @<name>:<bitflag>|<sid>
  const auto bitflag = GenerateUmiBitFlagFromFixedDuplexRead(read);
  WriteFastqHeaderWithBitflagAndSid(read, p_data, offset, bitflag);

  // Write UMI values: :<umi_5p>:<umi_3p>
  p_data[offset++] = kReadNameSeparator;
  WriteUmiValue(p_data, offset, mem.capacity, read.trim_info_duplex.umi_5p);
  p_data[offset++] = kReadNameSeparator;
  WriteUmiValue(p_data, offset, mem.capacity, read.trim_info_duplex.umi_3p);

  AppendFastqDuplexBody(read, p_data, offset, trimming_debug);
}

/**
 * @brief Write original read to buffer. Intended for failed reads.
 * @param read The read to write.
 * @param mem The buffer to write to.
 * @warning This function assumes data is unaltered, which may not hold true due to partial processing.
 */
void Formatter::FailedReadToBuffer(FixedReadRecord& read, FormattedOutput& mem) {
  auto& offset = mem.nr_bytes;
  auto* const p_data = mem.p_data;

  Formatter::ToFastqRaw(read, p_data, offset);
  read.SetStatus(FixedReadRecord::Status::kWritten);
}

/**
 * @brief Helper function to write a uint value to buffer using std::to_chars with error handling.
 * @param p_data Output buffer for the formatted data.
 * @param offset Current position in buffer; updated after write.
 * @param buffer_end End of the buffer.
 * @param value The value to write.
 * @param error_msg Error message prefix if conversion fails.
 */
void WriteUintToBuffer(
    char* const p_data, size_t& offset, const size_t buffer_end, const u32 value, const std::string_view error_msg) {
  const auto [ptr, ec] = std::to_chars(p_data + offset, p_data + buffer_end, value);
  if (ec != std::errc()) {
    throw error::Error("{} Either buffer too small or value invalid.\n", error_msg, std::make_error_code(ec).message());
  }
  offset = static_cast<size_t>(ptr - p_data);
}

/**
 * @brief Helper function to write an optional UMI value to buffer.
 * Writes the UMI value if present, otherwise writes '*'.
 * @param p_data Output buffer for the formatted data.
 * @param offset Current position in buffer; updated after write.
 * @param buffer_end End of the buffer.
 * @param umi_value Optional UMI value to write.
 */
void WriteUmiValue(char* const p_data, size_t& offset, const size_t buffer_end, const std::optional<u32>& umi_value) {
  if (umi_value.has_value()) {
    WriteUintToBuffer(p_data, offset, buffer_end, *umi_value, "UMI conversion error.");
  } else {
    p_data[offset++] = kMissingBarcodeChar;
  }
}

/**
 * @brief Helper function to write the FASTQ header with bitflag and SID.
 * Writes: @<name>:<bitflag>|<sid>
 * @param read The read record.
 * @param p_data Output buffer for the formatted data.
 * @param offset Current position in buffer; updated after write.
 * @param bitflag The bitflag value to write.
 */
void WriteFastqHeaderWithBitflagAndSid(const FixedReadRecord& read,
                                       char* const p_data,
                                       size_t& offset,
                                       const std::byte bitflag) {
  // Write '@' and name
  p_data[offset++] = kOriginalIdPrefix;
  AppendBytes(p_data, offset, {read.Name(), read.NameLen()});

  // Merge the threaded-state bits from the RDB read-state byte into the SID/UMI bitflag. This is the
  // single choke point for all read types, so every output path carries the threaded bits.
  const auto full_bitflag = bitflag | ThreadedBitsFromReadState(read.read_state);

  // Write ':' and bitflag
  p_data[offset++] = kReadNameSeparator;
  WriteUintToBuffer(p_data,
                    offset,
                    offset + kMaxU8DecimalDigits,
                    static_cast<u32>(std::to_underlying(full_bitflag)),
                    "Bitflag conversion error.");

  // Write '|' and SID.
  p_data[offset++] = kBitFlagSidSeparator;
  WriteUintToBuffer(p_data, offset, offset + kMaxU32DecimalDigits, read.sid_index, "SID conversion error.");
}

/**
 * @brief Format raw sequence tag for simplex reads.
 * Format: "sq:Z:<5p_trimmed><R><3p_trimmed>"
 *
 * The 5' trimmed region is [0, insert.spos), the 3' trimmed region is [insert.epos, SeqLen()).
 * The insert itself is replaced by the literal "<R>".
 * @param read The read record containing the raw sequence.
 * @param insert The insert range from trimming.
 * @return Formatted tag string.
 */
std::string FormatSimplexRawSeqTag(const FixedReadRecord& read, const LociRange& insert) {
  std::string result;
  result += debug::kSqTag;
  result.append(read.Seq(), insert.spos);
  result += debug::kSqSimplexInsertDelimiter;
  result.append(read.Seq() + insert.epos, read.SeqLen() - insert.epos);
  return result;
}

std::string FormatDuplexRawSeqTag(const FixedReadRecord& read) {
  const auto& raw = read.duplex_raw_seq;
  const auto& mid = read.trim_info_duplex.midadapter_range;
  std::string result;
  result += debug::kSqTag;
  if (raw.empty() || !mid.has_value()) {
    return result;
  }
  const auto hairpin_spos = mid->spos;
  const auto hairpin_epos_inclusive = mid->epos + 1;

  // Derive the 5' and 3' trimmed boundaries in the raw read.
  //
  // Raw read layout: [5p_trimmed][S1 insert][hairpin_trimmed][S2_rc insert][3p_trimmed]
  //
  // The hairpin position (hairpin_spos, hairpin_epos) is the ground truth found
  // by exact matching.  duplex_r1_len / duplex_r2_len (set during alignment) give
  // the trimmed insert lengths.  Trimmed boundaries are derived relative to the
  // hairpin:
  //   trim_5p_end   = hairpin_spos - R1_len
  //   trim_3p_start = hairpin_epos + 1 + R2_len
  const auto r1_len = read.duplex_r1_len;
  const auto r2_len = read.duplex_r2_len;
  const auto trim_5p_end = (r1_len <= hairpin_spos) ? hairpin_spos - r1_len : u32{0};
  const auto trim_3p_start = std::min(hairpin_epos_inclusive + r2_len, static_cast<u32>(raw.size()));

  result.append(raw, 0, trim_5p_end);
  result += debug::kSqDuplexR1Delimiter;
  result.append(raw, hairpin_spos, hairpin_epos_inclusive - hairpin_spos);
  result += debug::kSqDuplexR2Delimiter;
  result.append(raw, trim_3p_start, raw.size() - trim_3p_start);
  return result;
}

// bool args map 1:1 to TrimInfo SID flags
std::byte GenerateSidBitFlag(const bool has_sid_5p, const bool has_sid_3p) {  // NOSONAR(cpp:S5422)
  auto bitflag = std::byte{0};
  if (has_sid_5p) {
    bitflag |= kBitFlagSid5p;
  }
  if (has_sid_3p) {
    bitflag |= kBitFlagSid3p;
  }
  return bitflag;
}

// bool args map 1:1 to TrimInfo SID/UMI flags
std::byte GenerateUmiBitFlag(const bool has_sid_5p,  // NOSONAR(cpp:S5422)
                             const bool has_sid_3p,
                             const bool has_umi_5p,
                             const bool has_umi_3p) {
  auto bitflag = GenerateSidBitFlag(has_sid_5p, has_sid_3p);
  if (has_umi_5p) {
    bitflag |= kBitFlagUmi5p;
  }
  if (has_umi_3p) {
    bitflag |= kBitFlagUmi3p;
  }
  return bitflag;
}

/**
 * @brief Writes the read body (comment, sequence, and quality scores) to the buffer,
 *        assuming the read name has already been written.
 */
void AppendFastqSimplexBodyImpl(const FixedReadRecord& read,
                                char* const p_data,
                                size_t& offset,
                                const LociRange& insert,
                                const bool is_trimmed,
                                const SimplexFormatOptions& opts) {
  // Add space between name and comment/tags
  p_data[offset++] = ' ';

  // Copy the comment as-is.
  AppendBytes(p_data, offset, {read.Comment(), read.CommentLen()});

  // Append the raw-seq debug tag when that mode is active.
  if (HasTrimmingDebugMode(opts.trimming_debug, TrimmingDebugMode::kRawSeq)) [[unlikely]] {
    if (read.CommentLen() > 0) {
      p_data[offset++] = '\t';
    }
    AppendBytes(p_data, offset, FormatSimplexRawSeqTag(read, insert));
  }

  // Add a newline after the comment.
  p_data[offset++] = '\n';

  // Determine the source pointers and length based on whether the read was trimmed.
  const auto insert_offset = is_trimmed ? insert.spos : 0;
  const auto length = is_trimmed ? insert.Length() : read.SeqLen();

  // Copy the sequence.
  AppendBytes(p_data, offset, {read.Seq() + insert_offset, length});

  p_data[offset++] = '\n';
  p_data[offset++] = kQualitySeparator;
  p_data[offset++] = '\n';

  // When override_qual_score is true, replace them with the configured simplex quality score.
  if (opts.override_qual_score) {
    std::ranges::fill_n(p_data + offset, static_cast<std::ptrdiff_t>(length), opts.simplex_qual_score);
    offset += length;
  } else {
    AppendBytes(p_data, offset, {read.Qual() + insert_offset, length});
  }

  // Add a newline at the end of the quality scores.
  p_data[offset++] = '\n';
}

/**
 * @brief Write the record to the buffer based on the configured adapter type.
 * Reads adapter type and simplex formatting options from the enclosing flow context.
 * @param read The read record to write.
 * @param mem The buffer to write to.
 */
void Formatter::PassingReadToBuffer(FixedReadRecord& read, FormattedOutput& mem) const {
  using enum AdapterType;
  const auto& params = context.DemuxParam();
  const SimplexFormatOptions simplex_opts{.override_qual_score = params.override_simplex_qual_score,
                                          .simplex_qual_score = params.base_quality_scores.simplex,
                                          .trimming_debug = params.trimming_debug};
  switch (context.GetManager().AdapterTypeValue()) {
    case kYsu:
      ToFastqSimplexUMI(read, mem, read.trim_info_ysu, simplex_opts);
      break;
    case kYs:
    case kSimplex:
    case kSimplex10x:
      // These use the same TrimInfoSimplex layout and FASTQ formatting logic.
      ToFastqSimplex(read, mem, read.trim_info_simplex, simplex_opts);
      break;
    case kDuplexStem:
    case kDuplex:
      Formatter::ToFastqDuplex(read, mem.p_data, mem.nr_bytes, params.trimming_debug);
      if (mem.nr_bytes > mem.capacity) {
        throw error::Error("Sequence too long to store in internal buffer!");
      }
      break;
    case kDuplexUMI:
      ToFastqDuplexUMI(read, mem, params.trimming_debug);
      break;
    default:
      throw error::Error("Unknown adapter type {} encountered in Formatter",
                         std::to_underlying(context.GetManager().AdapterTypeValue()));
  }
  read.SetStatus(FixedReadRecord::Status::kWritten);
}

void Formatter::ToFastqDuplex(const FixedReadRecord& read,
                              char* const p_data,
                              size_t& offset,
                              const vec<TrimmingDebugMode>& trimming_debug) {
  // Write header: @<name>:<bitflag>|<sid>
  const auto bitflag = GenerateUmiBitFlagFromFixedDuplexRead(read);
  WriteFastqHeaderWithBitflagAndSid(read, p_data, offset, bitflag);

  // Write the body (YC tag, comment, consensus sequence, quality scores, etc.)
  AppendFastqDuplexBody(read, p_data, offset, trimming_debug);
}

void Formatter::ToFastqRaw(const FixedReadRecord& read, char* const p_data, size_t& offset) {
  p_data[offset++] = kOriginalIdPrefix;
  AppendBytes(p_data, offset, {read.Name(), read.NameLen()});
  // copy the comment
  p_data[offset++] = kCommentSeparator;
  AppendBytes(p_data, offset, {read.Comment(), read.CommentLen()});
  // Add a newline after the name.
  p_data[offset++] = '\n';

  // Copy the sequence.
  AppendBytes(p_data, offset, {read.Seq(), read.SeqLen()});

  // Separator line between the sequence and the quality scores.
  p_data[offset++] = '\n';
  p_data[offset++] = kQualitySeparator;
  p_data[offset++] = '\n';

  // Copy the quality scores.
  AppendBytes(p_data, offset, {read.Qual(), read.QualLen()});
  // Add a newline after the quality scores
  p_data[offset++] = '\n';
}

bool Formatter::ShouldSkipRecord(const FixedReadRecord& record) const {
  if (context.DemuxParam().output_failed_reads) {
    return record.GetStatus() == FixedReadRecord::Status::kWritten;
  }
  return record.GetStatus() != FixedReadRecord::Status::kDemultiplexed;
}

bool Formatter::FormatRecord(FixedReadRecord& record, const u32 sink_sid, FormattedOutput& mem) const {
  using enum FixedReadRecord::Status;

  switch (record.GetStatus()) {
    case kNotRead:
      throw error::Error("Internal error: record in Formatter with status {} encountered.",
                         static_cast<s32>(record.GetStatus()));
    case kWritten:
      // Nothing to do here, but skip buffer checking
    case kTooLongFail:
      // The read being too long means we don't have it in the buffer, so we can't write it out.
      // TODO: Find some way of preserving this read
      return false;
    case kRead:
      // shouldn't be possible but we don't always track status for simplex reads and will be considered failed
      // TODO: Make simplex adapters use failure states
    case kTooShortFail:
    case kTrimmedTooShortFail:
    case kDuplexMidAdapterFail:
    case kDuplexEditDistanceFail:
    case kDuplexTooLongFail:
    case kFailedMidadapterTrimFail:
    case kDiscordantSidFail:
    case kFullOnlyFail:
    case kUmiFail:
      return FormatFailedRecord(record, sink_sid, mem);
    case kDemultiplexed:
      if (sink_sid != record.sink_id) {
        return false;
      }
      PassingReadToBuffer(record, mem);
      return true;
    default:
      throw error::Error("Unknown FixedReadRecord status {} encountered in formatter.",
                         static_cast<s32>(record.GetStatus()));
  }
}

bool Formatter::FormatFailedRecord(FixedReadRecord& record, const u32 sid, FormattedOutput& mem) const {
  // Failed reads are only written when --write-failed-reads is set and the read is unassigned.
  if (!context.DemuxParam().output_failed_reads || sid != FixedReadRecord::kUnassignedSID) {
    return false;
  }
  // For simplex adapters, override quality scores in-place before writing raw output when requested.
  if (!context.GetManager().IsDuplex() && context.DemuxParam().override_simplex_qual_score) {
    std::ranges::fill_n(record.Qual(), record.QualLen(), context.DemuxParam().base_quality_scores.simplex);
  }
  FailedReadToBuffer(record, mem);
  return true;
}

void Formatter::ScheduleWriteTasks(const std::span<SinkData> sinks) const {
  auto& executor{context.GetManager().Executor()};

  if (!sinks.empty()) {
    if (sinks.size() > FormattedOutput::kMaxNumberSinks) {
      throw error::Error("Too many sinks. Expected at most {} sinks, but got {} sinks.",
                         FormattedOutput::kMaxNumberSinks,
                         sinks.size());
    }

    // Now create the write tasks. These tasks will not need to be stored, so creating them on the stack.
    std::array<tf::AsyncTask, FormattedOutput::kMaxNumberSinks> write_tasks{};
    for (size_t i = 0; i < sinks.size(); ++i) {
      write_tasks[i] = context.CreateWriteTask(sinks[i]);
    }
    context.nr_writes_scheduled += sinks.size();

    // Join task: completes once all per-sink write_tasks finish; no work of its own.
    context.SetWriteTask(
        batch_nr,
        executor.silent_dependent_async(
            [] { /* barrier task, no work */ }, write_tasks.begin(), write_tasks.begin() + sinks.size()));
  } else {
    // No sinks created, so we need to create a dummy write task.
    context.SetWriteTask(batch_nr, executor.silent_dependent_async([] { /* no-op placeholder write task */ }));
  }
}

void Formatter::operator()() {
  try {
    const StopWatch sw;
    auto& batch{context.GetBatchData(batch_nr)};

    size_t nr_sinks = 0;
    size_t current_index = 0;
    auto& mem{batch.formatted_output};
    mem.nr_bytes = 0;
    auto* const sink_data = batch.formatted_output.file_sinks.data();

    while (current_index != batch.num_records) {
      // Advance past records that are already handled or not eligible.
      while (current_index < batch.num_records && ShouldSkipRecord((*batch.records)[current_index])) {
        ++current_index;
      }
      if (current_index >= batch.num_records) {
        break;
      }

      if (nr_sinks >= FormattedOutput::kMaxNumberSinks) {
        throw error::Error(
            "Error condition detected: max number output sinks in a batch ({}) exceeded, generated output would be "
            "invalid.\n"
            "Saw {} sinks so far.",
            FormattedOutput::kMaxNumberSinks,
            nr_sinks);
      }

      const auto sink_sid = (*batch.records)[current_index].sink_id;
      const auto start_offset = mem.nr_bytes;
      auto& this_sink = sink_data[nr_sinks];
      this_sink.p_data = mem.p_data + start_offset;
      this_sink.sid_id = sink_sid;
      this_sink.batch_id = batch_nr;
      this_sink.flow_context = &context;

      // Format all records matching this SID into the buffer.
      for (auto index = current_index; index < batch.num_records; ++index) {
        if (!FormatRecord((*batch.records)[index], sink_sid, mem)) {
          continue;
        }
        // Written something — bounds-check the buffer.
        if (start_offset + mem.nr_bytes > mem.capacity) {
          throw error::Error("Write buffer not large enough. Bytes available: {}, bytes needed: {}",
                             mem.capacity,
                             start_offset + mem.nr_bytes);
        }
      }

      this_sink.length = mem.nr_bytes - start_offset;
      ++nr_sinks;
    }

    ScheduleWriteTasks({sink_data, nr_sinks});

    // Mark batch as completed, but this also schedules a new batch if input is available.
    context.MarkBatchCompleted();
    context.AddToFormatTime(static_cast<u64>(sw.ElapsedTime()));
  } catch (const std::exception& e) {
    Logging::Error("Formatter::operator() failed: {}", e.what());
    SetTaskException(std::current_exception());
  }
}
}  // namespace xoos::demux
