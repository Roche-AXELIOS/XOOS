#include "ios-util.h"

#include <fstream>
#include <ios>

#include <xoos/error/error.h>

namespace xoos::svc {

std::ofstream OpenOutput(const fs::path& output, std::ios::openmode mode) {
  std::ofstream out(output, mode);
  if (!out) {
    throw error::Error("Failed to open output file: {}", output);
  }
  return out;
}

std::ofstream OpenOutput(const fs::path& output) {
  return OpenOutput(output, std::ios::out | std::ios::trunc);
}

std::ifstream OpenInput(const fs::path& input, std::ios::openmode mode) {
  std::ifstream in(input, mode);
  if (!in) {
    throw error::Error("Failed to open input file: {}", input);
  }
  return in;
}

std::ifstream OpenInput(const fs::path& input) {
  return OpenInput(input, std::ios::in);
}

void RequireIoNotBad(const fs::path& path, std::ios& ios) {
  if (ios.bad()) {
    throw error::Error("Failed to read from file: {}", path);
  }
}

void RequireIoNotBad(std::ios& ios) {
  if (ios.bad()) {
    throw error::Error("Failed to read from stream");
  }
}

size_t Read(std::istream& in, vec<char>& buffer) {
  in.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
  return in.gcount();
}

std::string Read(std::istream& in) {
  const std::string content{(std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>()};
  RequireIoNotBad(in);
  return content;
}

std::string Read(const fs::path& path) {
  auto in = OpenInput(path, std::ios::binary);
  const auto content = Read(in);
  RequireIoNotBad(path, in);
  return content;
}

vec<std::string> ReadLines(std::istream& in) {
  vec<std::string> lines;
  std::string line;
  while (std::getline(in, line)) {
    lines.emplace_back(line);
  }
  RequireIoNotBad(in);
  return lines;
}

vec<std::string> ReadLines(const fs::path& path) {
  auto in = OpenInput(path);
  const auto lines = ReadLines(in);
  RequireIoNotBad(path, in);
  return lines;
}

void Write(std::ostream& out, const std::string& content) {
  out << content;
  RequireIoNotBad(out);
}

void Write(const fs::path& path, std::ios::openmode mode, const std::string& content) {
  auto out = OpenOutput(path, mode);
  Write(out, content);
}

void Write(const fs::path& path, const std::string& content) {
  Write(path, std::ios::out | std::ios::trunc, content);
}

void WriteLines(std::ostream& out, const vec<std::string>& lines) {
  for (const auto& line : lines) {
    out << line << '\n';
  }
  RequireIoNotBad(out);
}

void WriteLines(const fs::path& path, const vec<std::string>& lines) {
  auto out = OpenOutput(path, std::ios::out | std::ios::trunc);
  WriteLines(out, lines);
  RequireIoNotBad(path, out);
}
}  // namespace xoos::svc
