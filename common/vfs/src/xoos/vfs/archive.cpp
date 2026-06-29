#include "xoos/vfs/archive.h"

#include <archive.h>
#include <archive_entry.h>

#include <cstring>

namespace xoos::vfs {

ArchiveVirtualFileHandle::ArchiveVirtualFileHandle(const std::shared_ptr<struct archive>& archive) : _archive{archive} {
}

std::string ArchiveVirtualFileHandle::GetError() const {
  return archive_error_string(_archive.get());
}

ptrdiff_t ArchiveVirtualFileHandle::Read(char* buffer, ptrdiff_t buffer_size) {
  return archive_read_data(_archive.get(), buffer, buffer_size);
}

ArchiveVirtualFilesystem::ArchiveVirtualFilesystem(fs::path archive, size_t block_size)
    : _archive{std::move(archive)}, _block_size{block_size} {
}

std::vector<fs::path> ArchiveVirtualFilesystem::ListRecursiveImpl(size_t max_results) const {
  const ArchivePtr archive = OpenArchive();
  std::vector<fs::path> entries;
  while (entries.size() != max_results) {
    struct archive_entry* entry;
    const int r = archive_read_next_header(archive.get(), &entry);
    if (r == ARCHIVE_EOF) {
      break;
    }
    if (r != ARCHIVE_OK) {
      throw std::runtime_error("Unable to read archive entry: " + _archive.generic_string() + ": " + GetError(archive));
    }
    entries.emplace_back(archive_entry_pathname(entry));
    archive_read_data_skip(archive.get());
  }
  return entries;
}

VirtualFileHandlePtr ArchiveVirtualFilesystem::Open(const fs::path& path) const {
  const ArchivePtr archive = OpenArchive();

  struct archive_entry* entry;
  while (archive_read_next_header(archive.get(), &entry) == ARCHIVE_OK) {
    if (std::strcmp(path.c_str(), archive_entry_pathname(entry)) == 0) {
      return std::make_shared<ArchiveVirtualFileHandle>(archive);
    }
    archive_read_data_skip(archive.get());
  }

  return nullptr;
}

std::string ArchiveVirtualFilesystem::GetName() const {
  return _archive.string();
}

std::string ArchiveVirtualFilesystem::GetError(const ArchiveVirtualFilesystem::ArchivePtr& archive) {
  return archive_error_string(archive.get());
}

ArchiveVirtualFilesystem::ArchivePtr ArchiveVirtualFilesystem::OpenArchive() const {
  ArchivePtr archive = ArchivePtr(archive_read_new(), &archive_read_free);
  archive_read_support_filter_all(archive.get());
  archive_read_support_format_all(archive.get());

  const int r = archive_read_open_filename(archive.get(), _archive.c_str(), _block_size);
  if (r != ARCHIVE_OK) {
    throw std::runtime_error("Unable to open archive file: " + _archive.generic_string() + ": " + GetError(archive));
  }
  return archive;
}

}  // namespace xoos::vfs
