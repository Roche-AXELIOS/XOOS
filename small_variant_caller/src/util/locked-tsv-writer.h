#pragma once

#include <fstream>
#include <mutex>

#include <csv.hpp>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

namespace xoos::svc {

/**
 * @brief Thread-safe TSV writer class that manages file output with locking.
 */
class LockedTsvWriter {
 public:
  explicit LockedTsvWriter(const fs::path& file_path);

  /**
   * @brief Write the header row to the TSV file.
   * @param header Vector of strings representing the header row
   */
  void AppendRow(const vec<std::string>& row);

  /**
   * @brief Append rows to the TSV file in a thread-safe manner.
   * @param rows Vector of rows, where each row is a vector of strings
   */
  void AppendRows(const vec<vec<std::string>>& rows);

  void Flush();

 private:
  std::mutex _mutex;

  std::ofstream _ofs;
  csv::TSVWriter<std::ofstream> _writer;
};

}  // namespace xoos::svc
