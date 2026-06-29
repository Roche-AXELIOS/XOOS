#include "rdb-sequence-reader.h"

#include <fmt/format.h>
#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/util/string-functions.h>

#include <cstring>
#include <limits>
#include <tuple>

#include <nlohmann/json.hpp>

#include "rdb-2bit-utils.h"
#include "utility/math-util.h"
#include "xoos/enum/enum-util.h"

namespace xoos::demux {

// Constants for RDB file format
constexpr u32 kRunIdOffset = 40;
constexpr u32 kRunIdLength = 64;
constexpr u32 kBatchIndexOffset = 104;
constexpr u32 kBlockTypeOffset = 110;
constexpr u32 kNumDatasetsOffset = 112;
constexpr u32 kTotalDatasetBytesOffset = 128;
constexpr u32 kDatasetSizeOffset = 66;

constexpr u32 kNormalBlock = 0;
constexpr u32 kParameterBlock = 1;
constexpr u32 kChunkSize = 4096;
constexpr u32 kBlockInfoSize = 160;
constexpr u32 kDatasetInfoSize = 82;
// Maximum number of datasets that fit in a single RDB block header.
// The header is kChunkSize (4096) bytes; the first kBlockInfoSize (160) bytes are block info,
// and each dataset descriptor is kDatasetInfoSize (82) bytes, giving floor((4096-160)/82) = 48.
// Reading beyond this count would walk past the end of the header buffer.
constexpr u32 kMaxNumDatasets = (kChunkSize - kBlockInfoSize) / kDatasetInfoSize;  // 48 per RDB spec

constexpr u32 kMaxBase64UIDLength = 6;

// Datapath 2.0 format has 5 sections: date_group_num_instr_cd
// NS3_Gamma1 format has 5 sections.
constexpr u32 kNumPartsDp2OrGamma1 = 5;
// HTP 1.0 format has 6 sections: date_group_num_instr_cd_cycle
// NS3_Gamma2 format has 6 sections. NOTE: order name can have underscores.
constexpr u32 kNumPartsHtp1OrGamma2 = 6;

// if values are not specified by in the prefix queue num
constexpr auto kUnspecifiedQueueNumber = "Q*";

// constant to help logging format
constexpr std::string_view kRunIdTypeSuffix = "-run-id-type";

template <class T>
T GetField(const char* base, u16 offset) {
  return *(reinterpret_cast<const T*>(base + offset));
}

/**
 * Here we build the run prefix portion of the read name based on the run_id and cycle_id.
 *
 * The returned prefix will be:
 * {date}:{sequencer}:Q{queue_num}:R{run_num}:{cycle_id}
 *
 *
 *    HTP 1.0 format.
 *    220421_ENG-SYS-HTP_01_cloverleaf_WWX06R05C05_cycle01
 *    0      1           2  3          4            5
 *    date   group       │  instrument chip         cycle
 *                       run
 *   Datapath 2.0 format.
 *    220421_ENG-SYS-HTP_01_cloverleaf_WWX06R05C05
 *    0      1           2  3          4
 *    date   group       │  instrument chip
 *                       run
 *  NS3_Gamma1 format.
 *    20230526_2-00063_13_123456A_1
 *    0        1       2  3       4
 *    date     |       |  tube    cycle
 *             |       run
 *             instrument
 *  NS3_Gamma2 format.
 *    20230526_2-00063_Q13_R12_123456A_1
 *    0        1       2   3   4       5
 *    date     |       |   run order   tube
 *             |       queue
 *             instrument
 *
 * @brief Parse the run_id and determine the RDB run id type.
 *
 * If the run_id doesn't match expected formats, returns the placeholder
 * prefix "instrument:complex:chip:cycle" and kUnknownRdbRunIdType.
 *
 * Note that some formats include the cycle_id and some do not.
 *
 * @param run_id Prefix input.
 * @return The formatted prefix as
 *         {date}:{instrument}:Q{queue_num}:R{run_num}:{cycle_id, optional per format}
 *         and the corresponding RdbRunIdType.
 */
std::tuple<std::string, RdbRunIdType> BuildNamePrefix(const std::string& run_id) {
  // HTP 1.1 format has 5 sections: date_instr_cd_group_num
  const auto parts = string::Split(run_id, "_");
  const auto num_parts = parts.size();
  if (kNumPartsHtp1OrGamma2 <= num_parts) {
    // this is the Gamma2 format, which is in the same order and format as the desired read output, but wihtout the
    // cycle
    if (parts[2].starts_with("Q") && parts[3].starts_with("R")) {
      return {fmt::format("{}:{}:{}:{}", parts[0], parts[1], parts[2], parts[3]), RdbRunIdType::kGamma2RdbRunIdType};
    }
    // this is the HTP format, date, group, run, instrument, chip, cycle
    //    220421_ENG-SYS-HTP_01_cloverleaf_WWX06R05C05_cycle01
    //    0      1           2  3          4            5
    //    date   group       │  instrument chip         cycle
    //                       run
    const auto cycle_number = std::stoi(parts[5].substr(5, 2));
    return {fmt::format("{}:{}:{}:R{}:{}", parts[0], parts[3], kUnspecifiedQueueNumber, parts[2], cycle_number),
            RdbRunIdType::kHtp1RdbRunIdType};
  }
  if (parts.size() == kNumPartsDp2OrGamma1) {
    if (std::ranges::all_of(parts[4], [](const auto x) { return std::isdigit(x); })) {
      // If the 5th field is an integer, we have NS3_Gamma1 format.
      //    20230526_2-00063_13_123456A_1
      //    0        1       2  3       4
      //    date     |       |  tube    cycle
      //             |       run
      //             instrument
      return {fmt::format("{}:{}:{}:R{}:{}", parts[0], parts[1], kUnspecifiedQueueNumber, parts[2], parts[4]),
              RdbRunIdType::kGamma1RdbRunIdType};
    }
    // If the 5th field is not an integer, we have DP2 format.
    //    220421_ENG-SYS-HTP_01_cloverleaf_WWX06R05C05
    //    0      1           2  3          4
    //    date   group       │  instrument chip
    //                       run
    return {fmt::format("{}:{}:{}:R{}", parts[0], parts[3], kUnspecifiedQueueNumber, parts[2]),
            RdbRunIdType::kDataPath2RdbRunIdType};
  }
  // using placeholder that adheres to the documented output format {date}:{instrument}:Q*:R*
  // cycle_id will be appended by LoadBlock() for kUnknownRdbRunIdType
  return {"date:instrument:Q*:R*", RdbRunIdType::kUnknownRdbRunIdType};
}

constexpr auto kMainRdbFilename = "main_000001.rdb";
constexpr auto kCycleIdKey = "expstate::cycle_id";

/**
 *
 * cell-id is to be extracted from the parameters block in the
 * main_000001.rdb in the same folder as the sequence_file_path
 * e.g.
 *   "expstate::cycle_id": {
 *       "type": "uint32",
 *       "value": "1"
 *   },
 * @param sequence_file_path Path to the current RDB filename that should be parallel to the main_000001.rdb file that
 * contains the metadata
 * @return The cell_id from the metadata
 */
static u32 GetDataPath2CycleId(const fs::path& sequence_file_path) {
  auto main_file_path = sequence_file_path;
  main_file_path.replace_filename(kMainRdbFilename);
  // only the DataPAth2 objects have this.  If it is DataPath2 and this is missing then there is a formatting error
  if (!fs::exists(main_file_path)) {
    throw error::Error("Missing file '{}' for reading DataPath2 RDB format in path '{}'", kMainRdbFilename,
                       main_file_path.parent_path().string());
  }
  std::ifstream ifs(main_file_path, std::ios::binary);
  char header[kChunkSize];
  std::streamoff cur_block_start = 0;
  // Cache file size up front so we can validate untrusted fields (num_datasets, dataset_size)
  // against the actual file bounds before using them.
  const auto main_file_size_u64 = fs::file_size(main_file_path);
  const auto max_streamsize_u64 = static_cast<u64>(std::numeric_limits<std::streamsize>::max());
  const auto max_streamoff_u64 = static_cast<u64>(std::numeric_limits<std::streamoff>::max());
  if (main_file_size_u64 > max_streamoff_u64) {
    throw error::Error("File '{}' is too large to read with std::streamoff ({})", main_file_path.string(),
                       main_file_size_u64);
  }
  ifs.read(header, kChunkSize);
  while (ifs) {
    if (GetField<u16>(header, kBlockTypeOffset) == kParameterBlock) {
      auto dataset_info = header + kBlockInfoSize;
      auto dataset_offset = cur_block_start + kChunkSize;
      if (dataset_offset < 0) {
        ifs.close();
        throw error::Error("Corrupted RDB block in '{}': dataset offset ({}) is negative", main_file_path.string(),
                           dataset_offset);
      }
      auto dataset_offset_u64 = static_cast<u64>(dataset_offset);
      auto num_datasets = GetField<u16>(header, kNumDatasetsOffset);
      // FIX: Validate num_datasets against the maximum that fits in the header (48 per RDB spec).
      // Without this check, a corrupted num_datasets value causes the dataset_info pointer to walk
      // past the end of the stack-allocated header[] buffer, resulting in out-of-bounds reads.
      if (num_datasets > kMaxNumDatasets) {
        ifs.close();
        throw error::Error("Corrupted RDB block in '{}': num_datasets ({}) exceeds maximum ({})",
                           main_file_path.string(), num_datasets, kMaxNumDatasets);
      }
      for (u16 i = 0; i < num_datasets; ++i, dataset_info += kDatasetInfoSize) {
        auto dataset_name = dataset_info;
        const auto dataset_size_u64 = GetField<u64>(dataset_info, kDatasetSizeOffset);
        // Validate in unsigned space first, then cast for allocation/read().
        if (dataset_offset_u64 > main_file_size_u64) {
          ifs.close();
          throw error::Error("Corrupted dataset offset in '{}': dataset '{}' starts at {} beyond file size {}",
                             main_file_path.string(), dataset_name, dataset_offset_u64, main_file_size_u64);
        }
        const auto remaining_bytes = main_file_size_u64 - dataset_offset_u64;
        if (dataset_size_u64 > remaining_bytes || dataset_size_u64 > max_streamsize_u64) {
          ifs.close();
          throw error::Error("Corrupted dataset size in '{}': dataset '{}' claims {} bytes but only {} remain",
                             main_file_path.string(), dataset_name, dataset_size_u64, remaining_bytes);
        }
        const auto dataset_size = static_cast<std::streamsize>(dataset_size_u64);
        if (strcmp(dataset_name, "parameters") == 0) {
          ifs.seekg(static_cast<std::streamoff>(dataset_offset_u64), std::ios::beg);
          std::string parameters_json_string(dataset_size, '\0');
          ifs.read(parameters_json_string.data(), dataset_size);
          auto json = nlohmann::json::parse(parameters_json_string);
          if (json.contains(kCycleIdKey)) {
            auto cycle_id = std::stoi(static_cast<std::string>(json[kCycleIdKey]["value"]));
            ifs.close();
            return cycle_id;
          }
        }
        const auto rounded_dataset_size = RoundUp(dataset_size_u64, 8u);
        if (rounded_dataset_size > main_file_size_u64 - dataset_offset_u64) {
          ifs.close();
          throw error::Error("Corrupted dataset size in '{}': dataset '{}' rounded size {} exceeds {} remaining bytes",
                             main_file_path.string(), dataset_name, rounded_dataset_size,
                             main_file_size_u64 - dataset_offset_u64);
        }
        dataset_offset_u64 += rounded_dataset_size;
      }
    }
    cur_block_start +=
        kChunkSize + static_cast<std::streamoff>(RoundUp(GetField<u64>(header, kTotalDatasetBytesOffset), kChunkSize));
    ifs.seekg(cur_block_start, std::ios::beg);
    ifs.read(header, kChunkSize);
  }
  ifs.close();
  throw error::Error("Missing cycle-id value in '{}' file", main_file_path.string());
}

// NSA numeric base-64 alphabet – matches formatInteger<6, char*, uint32, 64, UNIFORM_UID_SIZE>.
// Order: 0-9, A-Z, a-z, +, _
static constexpr char kBase64Chars[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+_";

/**
 * Encode a u32 value to a base-64 string.
 *
 *   - Alphabet: 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+_
 *   - Pure numeric base-64 conversion (most-significant digit first)
 *   - Always produces exactly kMaxBase64UIDLength (6) characters, zero-padded
 *
 * 64^6 = 68,719,476,736 > 2^32, so every u32 value maps to a unique 6-character string.
 *
 * @param value The u32 cell index to encode
 * @return 6-character base-64 UID string
 */
std::string EncodeU32ToBase64(const u32 value) {
  std::string result(kMaxBase64UIDLength, '\0');
  u32 v = value;
  u32 divisor = 1073741824u;  // 64^5 – start with the most-significant digit
  for (u32 i = 0; i < kMaxBase64UIDLength; ++i) {
    result[i] = kBase64Chars[v / divisor];
    v %= divisor;
    divisor /= 64u;
  }
  return result;
}

RdbSequenceReader::RdbSequenceReader(const fs::path& sequence_file_path)
    : _file_path(sequence_file_path),
      _read_seq_stream(_file),
      _read_qual_stream(_file),
      _read_len_stream(_file),
      _read_cell_index_stream(_file) {
  _file.exceptions(std::ifstream::failbit | std::ifstream::badbit | std::ifstream::eofbit);
  _file.open(sequence_file_path, std::ios::binary);

  _file.seekg(0, std::ios::end);
  _file_size_bytes = _file.tellg();
  if (_file_size_bytes <= 0) {
    // We are logging here in case there is a crash prior to the error from another thread
    const auto error_message = fmt::format("Error: '{}' is empty or malformed.", sequence_file_path.string());
    Logging::Error(error_message);
    throw error::Error(error_message);
  }

  _file.seekg(0, std::ios::beg);

  _read_seq_stream.SetFileSize(_file_size_bytes);
  _read_qual_stream.SetFileSize(_file_size_bytes);
  _read_len_stream.SetFileSize(_file_size_bytes);
  _read_cell_index_stream.SetFileSize(_file_size_bytes);
  _cur_block_start = 0;

  // load first block
  _cur_block_size = 0;
  _num_reads_in_block = 0;
  _cur_read_index = 0;
  _cur_batch_index = 0;
  LoadBlock();
}

void RdbSequenceReader::LoadBlock() {
  char header[kChunkSize];

  try {
    _file.seekg(_cur_block_start, std::ios::beg);
    _file.read(header, kChunkSize);
  } catch (const std::ios_base::failure& e) {
    const auto error_message =
        fmt::format("Failed to read header from file '{}'.  Error: '{}'", _file_path.string(), e.what());
    // logging the error here in case there is a crash before the error message is displayed due to another thread's
    // failure
    Logging::Error(error_message);
    throw error::Error(error_message);
  }

  if (GetField<u16>(header, kBlockTypeOffset) != kNormalBlock) {
    throw error::Error("Not a normal block at {}", _cur_block_start);
  }

  const auto max_streamoff_u64 = static_cast<u64>(std::numeric_limits<std::streamoff>::max());
  if (_file_size_bytes < 0 || _cur_block_start < 0) {
    throw error::Error("Corrupted RDB reader state for '{}': negative file size ({}) or block offset ({})",
                       _file_path.string(), _file_size_bytes, _cur_block_start);
  }
  const auto file_size_u64 = static_cast<u64>(_file_size_bytes);
  const auto cur_block_start_u64 = static_cast<u64>(_cur_block_start);
  if (cur_block_start_u64 > file_size_u64 || file_size_u64 - cur_block_start_u64 < kChunkSize) {
    throw error::Error("Corrupted RDB block at offset {}: header exceeds remaining file size", _cur_block_start);
  }

  _num_reads_in_block = 0;
  _cur_read_index = 0;
  _cur_batch_index = GetField<u32>(header, kBatchIndexOffset);

  // We load the run_id and use that to identify the run type for collecting the cycle_id and read_name_prefix
  const std::string run_id(header + kRunIdOffset,
                           std::find(header + kRunIdOffset, header + kRunIdOffset + kRunIdLength, '\0'));
  const auto [prefix, run_id_type] = BuildNamePrefix(run_id);
  // Log the type of the run id once per session
  // remove the suffix -run-id-type from the formatted name for cleaner logging
  auto run_id_formatted = enum_util::FormatEnumName(run_id_type);
  if (run_id_formatted.ends_with(kRunIdTypeSuffix)) {
    run_id_formatted.erase(run_id_formatted.size() - kRunIdTypeSuffix.size());
  }
  XOOS_LOG_INFO_ONCE("Detected RDB input of type '{}'", run_id_formatted);
  _read_name_prefix = prefix;
  // For some run_id types, the cycle_id is not encoded in the run_id itself, so we must derive it from
  // other sources (e.g., file path for DataPath2 or batch index for Gamma2/unknown) and append it to the prefix.
  if (run_id_type == RdbRunIdType::kDataPath2RdbRunIdType) {
    // NOTE: this is a slower way to do it, but this should be a legacy RDB type
    _cycle_id = GetDataPath2CycleId(_file_path);
    // Append the cycle_id to the read name prefix.
    _read_name_prefix += fmt::format(":{}", _cycle_id.value());
  }
  // For Gamma2 or unknown run_id formats, derive cycle_id from the batch index stored in the block header.
  if (run_id_type == RdbRunIdType::kGamma2RdbRunIdType || run_id_type == RdbRunIdType::kUnknownRdbRunIdType) {
    _cycle_id = _cur_batch_index;
    // Append the cycle_id to the read name prefix.
    _read_name_prefix += fmt::format(":{}", _cycle_id.value());
  }
  // NOTE: Gamma1 and HTP1 run_id formats already have the cycle_id included in the read name prefix produced above.

  const auto total_dataset_bytes = GetField<u64>(header, kTotalDatasetBytesOffset);
  if (total_dataset_bytes == 0) {
    _cur_block_size = kChunkSize;
    // no data in this block
    return;
  }

  const auto dataset_section_start_u64 = cur_block_start_u64 + kChunkSize;
  const auto remaining_dataset_bytes = file_size_u64 - dataset_section_start_u64;
  // FIX: Sanity-check total_dataset_bytes against the remaining file size.  A corrupted value
  // could cause _cur_block_size to overshoot, making the next LoadBlock() seek past EOF and
  // producing undefined results or file read errors.
  if (total_dataset_bytes > remaining_dataset_bytes) {
    throw error::Error("Corrupted RDB block at offset {}: total_dataset_bytes ({}) exceeds remaining file size",
                       _cur_block_start, total_dataset_bytes);
  }
  const auto rounded_total_dataset_bytes = RoundUp(total_dataset_bytes, static_cast<u64>(kChunkSize));
  if (rounded_total_dataset_bytes > remaining_dataset_bytes ||
      rounded_total_dataset_bytes > max_streamoff_u64 - kChunkSize) {
    throw error::Error("Corrupted RDB block at offset {}: rounded total_dataset_bytes ({}) exceeds remaining file size",
                       _cur_block_start, rounded_total_dataset_bytes);
  }
  _cur_block_size = kChunkSize + static_cast<std::streamoff>(rounded_total_dataset_bytes);

  const auto num_datasets = GetField<u16>(header, kNumDatasetsOffset);
  // FIX: Validate num_datasets against the maximum that can fit in the header (48 per RDB spec).
  // Without this check, the for-loop below advances `dataset_info` past the end of the
  // stack-allocated header[] buffer, causing out-of-bounds reads.
  if (num_datasets > kMaxNumDatasets) {
    throw error::Error("Corrupted RDB block at offset {}: num_datasets ({}) exceeds maximum ({})", _cur_block_start,
                       num_datasets, kMaxNumDatasets);
  }
  auto dataset_offset_u64 = dataset_section_start_u64;
  auto dataset_info = header + kBlockInfoSize;
  constexpr auto kNumDatasetsRequired = 5;
  int num_datasets_found = 0;
  for (u16 i = 0; i < num_datasets; ++i, dataset_info += kDatasetInfoSize) {
    const auto dataset_name = dataset_info;
    if (dataset_offset_u64 > file_size_u64 || dataset_offset_u64 > max_streamoff_u64) {
      throw error::Error("Corrupted dataset offset in '{}': dataset '{}' starts at {} beyond file size {}",
                         _file_path.string(), dataset_name, dataset_offset_u64, file_size_u64);
    }
    const auto dataset_size_u64 = GetField<u64>(dataset_info, kDatasetSizeOffset);
    const auto remaining_bytes = file_size_u64 - dataset_offset_u64;
    if (dataset_size_u64 > remaining_bytes || dataset_size_u64 > max_streamoff_u64) {
      throw error::Error("Corrupted dataset size in '{}': dataset '{}' claims {} bytes but only {} remain",
                         _file_path.string(), dataset_name, dataset_size_u64, remaining_bytes);
    }
    const auto dataset_offset = static_cast<std::streamoff>(dataset_offset_u64);
    if (strcmp(dataset_name, "molecular_traces_count") == 0) {
      _file.seekg(dataset_offset, std::ios::beg);
      _file.read(reinterpret_cast<char*>(&_num_reads_in_block), sizeof(u32));
      ++num_datasets_found;
    } else if (strcmp(dataset_name, "molecular_traces") == 0) {
      _read_seq_stream.Reset(dataset_offset);
      ++num_datasets_found;
    } else if (strcmp(dataset_name, "molecular_traces_qscores") == 0) {
      _read_qual_stream.Reset(dataset_offset);
      ++num_datasets_found;
    } else if (strcmp(dataset_name, "molecular_traces_read_length") == 0) {
      _read_len_stream.Reset(dataset_offset);
      ++num_datasets_found;
    } else if (strcmp(dataset_name, "mt_cell_index") == 0) {
      _read_cell_index_stream.Reset(dataset_offset);
      ++num_datasets_found;
    }
    const auto rounded_dataset_size = RoundUp(dataset_size_u64, 8u);
    if (rounded_dataset_size > file_size_u64 - dataset_offset_u64) {
      throw error::Error("Corrupted dataset size in '{}': dataset '{}' rounded size {} exceeds {} remaining bytes",
                         _file_path.string(), dataset_name, rounded_dataset_size, file_size_u64 - dataset_offset_u64);
    }
    dataset_offset_u64 += rounded_dataset_size;
  }
  if (num_datasets_found != kNumDatasetsRequired) {
    throw error::Error("Expected datasets are not available for block at {}", _cur_block_start);
  }
}

BatchStatistics RdbSequenceReader::ReadBatchIntoArena(FixedReadRecordBatch& batch) {
  batch.start_time = std::chrono::high_resolution_clock::now();
  batch.num_bytes = 0ul;
  batch.num_bases = 0ul;
  for (batch.num_records = 0ul; batch.num_records < batch.Capacity() && HasMoreReads(); ++batch.num_records) {
    auto& record = (*batch.records)[batch.num_records];
    batch.num_bytes += ReadSingleFixed(record);
    batch.num_bases += record.SeqLen();
  }
  batch.end_time = std::chrono::high_resolution_clock::now();
  return {batch.num_records, batch.num_bases};
}

u32 RdbSequenceReader::ReadSingleFixed(FixedReadRecord& rec) {
  if (_cur_read_index >= _num_reads_in_block) {
    _cur_block_start += _cur_block_size;
    LoadBlock();
  }

  rec.Clear();

  auto read_length = _read_len_stream.Read<u32>();
  // Guard: a zero-length read is never valid in practice and is a strong indicator of stream
  // desynchronization or file corruption.  Fail fast instead of silently producing empty records
  // that waste I/O cycles (which on networked storage can manifest as long CPU-near-0 stalls).
  if (read_length == 0) {
    auto error_string = fmt::format(
        "Read length of 0 encountered at read index {} in block at offset {} in file '{}'. "
        "This likely indicates file corruption or an incomplete RDB file.",
        _cur_read_index, _cur_block_start, _file_path.string());
    Logging::Error(error_string);
    throw error::Error(error_string);
  }
  auto read_length_packed = CeilDiv(read_length, 4u);
  u32 bytes_read = 0;

  if (read_length <= kMaxReadLength) {
    const auto read_cell_index = _read_cell_index_stream.Read<u32>();
    // We take the u32 and encode it to Base64. The buffer will always be there on the last 2 so we always exclude it by
    // truncating to 6 characters, providing a fixed length string for the read name.
    std::string b64_uid = EncodeU32ToBase64(read_cell_index);

    const auto [name_out, name_size] =
        fmt::format_to_n(rec.Name(), kBufferSize - 1, "{}:{}", _read_name_prefix, b64_uid);

    rec.comment_offset = name_out - rec.Name();

    rec.seq_offset = rec.comment_offset;
    _read_seq_stream.Read(rec.TwoBitsSeq(), read_length_packed);
    Reverse2BitOrder(rec.TwoBitsSeq(), read_length_packed);

    // fill padding regions with 0s
    std::memset(rec.TwoBitsSeq() - kTwoBitPadding, 0, kTwoBitPadding);
    // FIX: Was `+ read_length_packed + 1` which is off-by-one.  The +1 skips a byte and at
    // kMaxReadLength writes 1 byte past the end of the two_bit_seq[] array.
    // InitTwoBit() in read-record.h correctly uses `+ num_bytes` without +1.
    std::memset(rec.TwoBitsSeq() + read_length_packed, 0, kTwoBitPadding);

    DecodeDnaBases(rec.TwoBitsSeq(), rec.Seq(), read_length);

    rec.qual_offset = rec.seq_offset + read_length;
    _read_qual_stream.Read(_qual_buf, read_length_packed);
    DecodeQualScores(_qual_buf, rec.Qual(), read_length);

    rec.end_offset = rec.qual_offset + read_length;

    bytes_read = read_length_packed * 2;
    rec.SetStatus(FixedReadRecord::Status::kRead);
  } else {
    // FIX: Even though we discard this read, we must advance all per-read streams past its
    // data so that the next call to ReadSingleFixed() reads from the correct offsets.
    // Without these skips, every subsequent read in the block parses garbled data because
    // the cell-index, sequence, and quality streams are still pointing at *this* read's bytes.
    _read_cell_index_stream.Skip(sizeof(u32));
    _read_seq_stream.Skip(read_length_packed);
    _read_qual_stream.Skip(read_length_packed);
    rec.SetStatus(FixedReadRecord::Status::kTooLongFail);
  }
  ++_cur_read_index;
  return bytes_read;
}

bool RdbSequenceReader::HasMoreReads() const {
  return (_cur_read_index < _num_reads_in_block) || (_cur_block_start + _cur_block_size < _file_size_bytes);
}

}  // namespace xoos::demux
