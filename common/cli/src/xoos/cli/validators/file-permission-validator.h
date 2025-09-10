#pragma once

#include <optional>
#include <variant>

#include <fmt/format.h>

#include <CLI/CLI.hpp>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

namespace xoos::cli {

std::optional<std::string> CheckFileReadable(const fs::path& path);

CLI::Validator FileReadableValidator();

std::optional<std::string> CheckFileWriteable(const fs::path& path);

CLI::Validator FileWriteableValidator();

std::variant<std::string, vec<fs::path>> CheckFileListReadable(const fs::path& path);

CLI::Validator FileListReadableValidator();

}  // namespace xoos::cli
