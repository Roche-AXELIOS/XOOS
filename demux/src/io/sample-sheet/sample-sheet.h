#pragma once

#include <filesystem>
#include <string>
#include <unordered_map>

namespace xoos::demux {

namespace fs = std::filesystem;

std::unordered_map<std::string, std::string> Read(const fs::path& path);

}  // namespace xoos::demux
