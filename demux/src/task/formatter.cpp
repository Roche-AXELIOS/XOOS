#include "task/formatter.h"

#include <fmt/format.h>
#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include <cstddef>

#include "task/flow-context.h"
#include "task/flow-manager.h"
#include "utility/stop-watch.h"

namespace xoos::demux {
Formatter::Formatter(FlowContext& exec, size_t batch_nr) : Task(exec, batch_nr, fmt::format("Writer{}", batch_nr)) {}

/**
 * @brief Generates a bitflag for the read based on which UMIs are present in the TrimInfoDuplexUMI. The bitflag is used
 * to indicate the presence of 5' and 3' SIDs and UMIs in the output FASTQ header for duplex reads.
 * @param read  The read to generate the bitflag for. We will look at the TrimInfoDuplexUMI to determine which UMIs were
 * present and set the appropriate bits in the bitflag.
 * @return
 */
u8 GenerateUmiBitFlagFromFixedDuplexRead(const FixedReadRecord& read) {
  // 5' sid: 0b1000 - always true for duplex
  // 3' sid: 0b0100 - always true for duplex
  // 5' umi: 0b0010
  // 3' umi: 0b0001
  // we assume both SIDS are present for duplex and create the appropriate flag
  std::byte bitflag{0b1000 | 0b0100};

  if (read.trim_info_duplex.umi_5p.has_value()) {
    bitflag |= std::byte{0b0010};
  }
  if (read.trim_info_duplex.umi_3p.has_value()) {
    bitflag |= std::byte{0b0001};
  }
  return std::to_integer<u8>(bitflag);
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
void AppendFastqDuplexBody(const FixedReadRecord& read, char* const p_data, size_t& offset) {
  // Add the YC tag
  std::memcpy(p_data + offset, read.Seq(), read.iupac_length);
  offset += read.iupac_length;

  // Copy the comment as-is - and don't forget to add the '\n' at the end.
  std::memcpy(p_data + offset, read.Comment(), read.CommentLen());
  offset += read.CommentLen();
  // Add a newline after the comment.
  p_data[offset++] = '\n';

  const auto consensus_len = static_cast<size_t>(read.consensus_seq_len);

  // The sequence is now a consensus sequence.
  std::memcpy(p_data + offset, read.ConsensusSeq(), consensus_len);
  offset += consensus_len;

  // Add a newline at the end of the sequence.
  p_data[offset++] = '\n';
  // Add the '+' character after the sequence.
  p_data[offset++] = kQualitySeparator;
  // Add a newline after the '+' character.
  p_data[offset++] = '\n';

  // Copy the quality scores; these artificial scores have replaced the sequence.
  std::memcpy(p_data + offset, read.Qual(), consensus_len);
  offset += consensus_len;

  // Add a newline at the end of the quality scores.
  p_data[offset++] = '\n';
}

/**
 * @brief Write the read to the buffer, using the DuplexUmi bundle to get the UMI sequences.
 * @param read The read to write.
 * @param mem The buffer to write to.
 * @details Implementation is similar to the Duplex version, but we add the UMI sequences to the read name.
 * Output format: @<original id>|7|3 <original suffix>
 * '*' is used to indicate missing barcodes.
 */
void ToFastqDuplexUMI(const FixedReadRecord& read, FormattedOutput& mem) {
  // Copy over comment prefixed with '@' character.
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

  AppendFastqDuplexBody(read, p_data, offset);
}

/**
 * @brief Write original read to buffer. Intended for failed reads.
 * @param read The read to write.
 * @param mem The buffer to write to.
 * @warning This function assumes data is unaltered, which may not hold true due to partial processing.
 */
void FailedReadToBuffer(FixedReadRecord& read, FormattedOutput& mem) {
  // Copy over comment prefixed with '@' character.
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
void WriteUintToBuffer(char* const p_data, size_t& offset, const size_t buffer_end, const u32 value,
                       const std::string_view error_msg) {
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
    WriteUintToBuffer(p_data, offset, buffer_end, umi_value.value(), "UMI conversion error.");
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
void WriteFastqHeaderWithBitflagAndSid(const FixedReadRecord& read, char* const p_data, size_t& offset,
                                       const u8 bitflag) {
  // Write '@' and name
  p_data[offset++] = kOriginalIdPrefix;
  std::memcpy(p_data + offset, read.Name(), read.NameLen());
  offset += read.NameLen();

  // Write ':' and bitflag
  p_data[offset++] = kReadNameSeparator;
  WriteUintToBuffer(p_data, offset, offset + 3, bitflag, "Bitflag conversion error.");

  // Write '|' and SID
  p_data[offset++] = kBitFlagSidSeparator;
  const auto read_sid = std::to_string(read.file_sid);
  std::memcpy(p_data + offset, read_sid.c_str(), read_sid.size());
  offset += read_sid.size();
}

/**
 * @brief Write the record to the buffer based on adapter type.
 * @param read The read record to write.
 * @param mem The buffer to write to.
 * @param adapter_type The type of adapter used for demultiplexing.
 */
void PassingReadToBuffer(FixedReadRecord& read, FormattedOutput& mem, const AdapterType adapter_type) {
  using enum AdapterType;
  switch (adapter_type) {
    case kYsuSl:
      ToFastqSimplexUMI(read, mem, read.trim_info_ysu_sl);
      break;
    case kYsuTe:
      ToFastqSimplexUMI(read, mem, read.trim_info_ysu_te);
      break;
    case kYs:
      ToFastqSimplex(read, mem, read.trim_info_ys);
      break;
    case kDuplexStem:
    case kDuplex:
      Formatter::ToFastqDuplex(read, mem.p_data, mem.nr_bytes);
      if (mem.nr_bytes > mem.capacity) {
        throw error::Error("Sequence too long to store in internal buffer!");
      }
      break;
    case kDuplexUMI:
      ToFastqDuplexUMI(read, mem);
      break;
    default:
      throw error::Error("Unknown adapter type {} encountered in Formatter", static_cast<s32>(adapter_type));
  }
  read.SetStatus(FixedReadRecord::Status::kWritten);
}

void Formatter::ToFastqDuplex(const FixedReadRecord& read, char* const p_data, size_t& offset) {
  // Write header: @<name>:<bitflag>|<sid>
  const auto bitflag = GenerateUmiBitFlagFromFixedDuplexRead(read);
  WriteFastqHeaderWithBitflagAndSid(read, p_data, offset, bitflag);

  // Write the body (YC tag, comment, consensus sequence, quality scores, etc.)
  AppendFastqDuplexBody(read, p_data, offset);
}

void Formatter::ToFastqRaw(const FixedReadRecord& read, char* const p_data, size_t& offset) {
  p_data[offset++] = kOriginalIdPrefix;
  std::memcpy(p_data + offset, read.Name(), read.NameLen());
  offset += read.NameLen();
  // copy the comment
  p_data[offset++] = kCommentSeparator;
  std::memcpy(p_data + offset, read.Comment(), read.CommentLen());
  offset += read.CommentLen();
  // Add a newline after the name.
  p_data[offset++] = '\n';

  // Copy the sequence.
  std::memcpy(p_data + offset, read.Seq(), read.SeqLen());
  offset += read.SeqLen();

  // Separator line between the sequence and the quality scores.
  p_data[offset++] = '\n';
  p_data[offset++] = kQualitySeparator;
  p_data[offset++] = '\n';

  // Copy the quality scores.
  std::memcpy(p_data + offset, read.Qual(), read.QualLen());
  offset += read.QualLen();
  // Add a newline after the quality scores
  p_data[offset++] = '\n';
}

void Formatter::operator()() {
  try {
    const StopWatch sw;
    auto& batch{context.GetBatchData(batch_nr)};

    // Loop over all the records in the batch and convert them to output format but cluster them by SID so we can
    // hand them over to the sink as a single block of memory.
    auto& mgr{context.GetManager()};

    size_t nr_sinks = 0;

    size_t current_index = 0;
    auto& mem{batch.formatted_output};
    mem.nr_bytes = 0;
    auto& sink_data{batch.formatted_output.file_sinks};

    while (current_index != batch.num_records) {
      // Find the first record that has not been written out yet.

      while (current_index < batch.num_records &&
             (context.DemuxParam().output_failed_reads
                  ? (*batch.records)[current_index].GetStatus() == FixedReadRecord::Status::kWritten
                  : (*batch.records)[current_index].GetStatus() != FixedReadRecord::Status::kDemultiplexed)) {
        ++current_index;
      }
      if (current_index < batch.num_records) {
        if (nr_sinks >= FormattedOutput::kMaxNumberSinks) {
          throw error::Error(
              "Error condition detected: max number SIDs in a batch ({}) exceeded, generated output would be invalid.\n"
              "Saw {} sinks so far.",
              FormattedOutput::kMaxNumberSinks, nr_sinks);
        }
        // Get the SID of the current record
        const auto sid = (*batch.records)[current_index].file_sid;

        // Collect all records with the same SID. We will write them out in a single go.
        const auto start_offset = mem.nr_bytes;
        auto& this_sink = sink_data[nr_sinks];
        this_sink.p_data = mem.p_data + start_offset;
        this_sink.sid_id = sid;
        this_sink.batch_id = batch_nr;
        this_sink.flow_context = &context;

        for (auto index = current_index; index < batch.num_records; ++index) {
          auto& record = (*batch.records)[index];
          using enum FixedReadRecord::Status;
          switch (record.GetStatus()) {
            case kNotRead:
              // We should never see these states here
              throw error::Error("Internal error: record in Formatter with status {} encountered.",
                                 static_cast<s32>(record.GetStatus()));
            case kWritten:
              // Nothing to do here, but skip buffer checking
            case kTooLongFail:
              // The read being too long means we don't have it in the buffer, so we can't write it out.
              // TODO: Find some way of preserving this read
              continue;
            case kRead:
              // shouldn't be possible but we don't always track status for simplex reads and will be considered failed
              // TODO: Make simplex adapters use failure states
            case kTooShortFail:
            case kTrimmedTooShortFail:
            case kDuplexMidAdapterFail:
            case kDuplexEditDistanceFail:
            case kDuplexTooLongFail:
            case kFailedMidadapterTrimFail:
              // Failed reads should only be written if output_failed_reads is set
              if (context.DemuxParam().output_failed_reads && sid == FixedReadRecord::kUnassignedSID) {
                FailedReadToBuffer(record, mem);
                // We have written, do a buffer check
                break;
              }
              // Skip buffer checking
              continue;
            case kDemultiplexed:
              if (sid == record.file_sid) {
                // This record belongs to the current SID, so write it to the buffer.
                PassingReadToBuffer(record, mem, mgr.adapter_type);
                // We have written, do a buffer check
                break;
              }
              // Different SID, so nothing written; skip check
              continue;
            default:
              throw error::Error("Unknown FixedReadRecord status {} encountered in formatter.",
                                 static_cast<s32>(record.GetStatus()));
          }
          // If we reach here, we have written something to the buffer, so we need to check that we didn't overflow.
          if (start_offset + mem.nr_bytes > mem.capacity) {
            throw error::Error("Write buffer not large enough. Bytes available: {}, bytes needed: {}", mem.capacity,
                               start_offset + mem.nr_bytes);
          }
        }
        this_sink.length = mem.nr_bytes - start_offset;
        ++nr_sinks;
      }
    }

    auto& executor{context.GetManager().Executor()};

    if (nr_sinks > 0) {
      if (nr_sinks > FormattedOutput::kMaxNumberSinks) {
        throw error::Error("Too many sinks. Expected at most {} sinks, but got {} sinks.",
                           FormattedOutput::kMaxNumberSinks, nr_sinks);
      }

      // Now create the write tasks. These tasks will not need to be stored, so creating them on the stack.
      // TODO: need to look into dynamic allocation
      tf::AsyncTask write_tasks[FormattedOutput::kMaxNumberSinks] = {};
      for (size_t i = 0; i < nr_sinks; ++i) {
        // Every sink has a write task that is dependent on the preceding write task, but the sink figured that one out.
        write_tasks[i] = context.CreateWriteTask(sink_data[i]);
      }
      context.nr_writes_scheduled += nr_sinks;

      // Create a write task that is dependent on all the sinks that were created.
      context.SetWriteTask(batch_nr, executor.silent_dependent_async([]() {}, write_tasks, write_tasks + nr_sinks));
    } else {
      // No sinks created, so we need to create a dummy write task.
      context.SetWriteTask(batch_nr, executor.silent_dependent_async([]() {}));
    }

    // Mark batch as completed, but this also schedules a new batch if input is available.
    context.MarkBatchCompleted();
    context.AddToFormatTime(static_cast<u64>(sw.ElapsedTime()));
  } catch (const std::exception& e) {
    Logging::Error("Formatter::operator() failed: {}", e.what());
    SetTaskException(std::current_exception());
  }
}
}  // namespace xoos::demux
