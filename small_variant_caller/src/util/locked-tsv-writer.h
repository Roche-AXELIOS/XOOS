#pragma once

#include <fstream>
#include <mutex>

#include <csv.hpp>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "xoos/io/metadata-util.h"

namespace xoos::svc {

/**
 * @brief Thread-safe TSV writer class that manages file output with locking.
 */
class LockedTsvWriter {
 public:
  explicit LockedTsvWriter(const fs::path& file_path);

  /**
   * @brief Write a row to the TSV file.
   * @param row Vector of strings representing a single row
   */
  void AppendRow(const vec<std::string>& row);

  /**
   * @brief Append rows to the TSV file in a thread-safe manner.
   * @param rows Vector of rows, where each row is a vector of strings
   */
  void AppendRows(const vec<vec<std::string>>& rows);

  /**
   * @brief Write version and command line metadata to the TSV file.
   * @param info Command line info containing version and command line
   */
  void AppendMetadata(const io::CommandLineInfo& info);

  /**
   * @brief Write only a ##Sample=<...> metadata line to the TSV file.
   * @param read_group Read group whose fields populate the sample line
   *
   * Written directly to the underlying stream to avoid CSV quoting. Use when the command line
   * metadata is written separately or omitted.
   */
  void AppendSampleMetadata(const io::ReadGroup& read_group);

  void Flush();

 private:
  std::mutex _mutex;

  std::ofstream _ofs;
  csv::TSVWriter<std::ofstream> _writer;
};

}  // namespace xoos::svc
