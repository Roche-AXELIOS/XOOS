#include "locked-tsv-writer.h"

#include <xoos/error/error.h>

namespace xoos::svc {

static std::ofstream OpenStream(std::mutex& mutex, const fs::path& file_path) {
  std::scoped_lock lock{mutex};
  if (!file_path.parent_path().empty()) {
    fs::create_directories(file_path.parent_path());
  }
  std::ofstream ofs(file_path, std::ios::out);
  if (!ofs.is_open()) {
    throw error::Error("Failed to open file: {}", file_path);
  }
  return ofs;
}

LockedTsvWriter::LockedTsvWriter(const fs::path& file_path) : _ofs(OpenStream(_mutex, file_path)), _writer(_ofs) {
}

void LockedTsvWriter::AppendRow(const vec<std::string>& row) {
  std::scoped_lock lock{_mutex};
  _writer << row;
}

void LockedTsvWriter::AppendRows(const vec<vec<std::string>>& rows) {
  std::scoped_lock lock{_mutex};
  for (const auto& row : rows) {
    _writer << row;
  }
}

void LockedTsvWriter::Flush() {
  std::scoped_lock lock{_mutex};
  _ofs.flush();
}
}  // namespace xoos::svc
