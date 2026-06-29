#include "sequence-reader.h"

#include <xoos/error/error.h>

#include <algorithm>
#include <string>
#include <vector>

#include "fastq/kseqpp-sequence-reader.h"
#include "rdb/rdb-sequence-reader.h"

namespace xoos::demux {
bool EndsWithOneOf(const fs::path& p, const std::vector<std::string>& suffixes) {
  std::string ps{p.string()};
  std::transform(std::cbegin(ps), std::cend(ps), std::begin(ps), [](auto c) { return std::tolower(c); });
  return std::any_of(suffixes.cbegin(), suffixes.cend(), [&ps](const auto& suffix) { return ps.ends_with(suffix); });
}

bool FileNameBeginsWith(const fs::path& p, const std::string& prefix) {
  std::string ps{p.filename().string()};
  std::transform(std::cbegin(ps), std::cend(ps), std::begin(ps), [](auto c) { return std::tolower(c); });
  return ps.starts_with(prefix);
}

bool ParentPathEndsWith(const fs::path& p, const std::string& suffix) {
  std::string ps{p.parent_path().string()};
  std::transform(std::cbegin(ps), std::cend(ps), std::begin(ps), [](auto c) { return std::tolower(c); });
  return ps.ends_with(suffix);
}

bool IsFastqSequenceFileFormat(const fs::path& sequence_file_path) {
  return EndsWithOneOf(sequence_file_path, {".fastq", ".fastq.gz", ".fq", ".fq.gz"});
}

bool IsRdbSequenceFileFormat(const fs::path& sequence_file_path) {
  return (sequence_file_path.filename().string().starts_with("PrimaryAnalysis_") ||
          sequence_file_path.filename().string().starts_with("ExtractMt_")) &&
         sequence_file_path.string().ends_with(".rdb");
}

bool IsSequenceFileFormat(const fs::path& sequence_file_path) {
  return IsFastqSequenceFileFormat(sequence_file_path) || IsRdbSequenceFileFormat(sequence_file_path);
}

std::shared_ptr<SequenceReader> OpenSequenceFile(const fs::path& sequence_file_path) {
  if (IsFastqSequenceFileFormat(sequence_file_path)) {
    return std::make_shared<KSeqppSequenceReader>(sequence_file_path);
  }
  if (IsRdbSequenceFileFormat(sequence_file_path)) {
    return std::make_shared<RdbSequenceReader>(sequence_file_path);
  }
  throw error::Error("Unrecognized sequence file type: '{}'", sequence_file_path.string());
}
}  // namespace xoos::demux
