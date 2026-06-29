#include "xoos/io/metadata-util.h"

#include <optional>

#include <fmt/format.h>

namespace xoos::io {

const std::string kTsvCommentKeyValueSeparator = "=";
const std::string kTsvVersionKey = "Version";
const std::string kTsvCommandLineKey = "CommandLine";

// insert a comment into the comments vector
void AddComment(Comments& comments, const std::string& comment) {
  comments.emplace_back(comment);
}

// insert the key and value into the comments vector as a key=value pair
void AddCommentKeyPair(Comments& comments, const std::string& key, const std::string& value) {
  comments.emplace_back(key + kTsvCommentKeyValueSeparator + value);
}

void AddVersionAndCommandLineComment(Comments& comments, const std::string& version, const std::string& command_line) {
  AddCommentKeyPair(comments, kTsvVersionKey, version);
  AddCommentKeyPair(comments, kTsvCommandLineKey, command_line);
}

/**
 * @brief Extract the version from the comments.
 * @param comments The vector of comments.
 * @return An optional string containing the version if found, otherwise std::nullopt.
 */
std::optional<std::string> GetVersion(const Comments& comments) {
  // look for the version key in the comments
  for (const auto& comment : comments) {
    if (comment.rfind(kTsvVersionKey + kTsvCommentKeyValueSeparator, 0) == 0) {
      return comment.substr(kTsvVersionKey.size() + kTsvCommentKeyValueSeparator.size());
    }
  }
  return std::nullopt;
}

/**
 * @brief Extract the command line from the comments.
 * @param comments The vector of comments.
 * @return An optional string containing the command line if found, otherwise std::nullopt.
 */
std::optional<std::string> GetCommandLine(const Comments& comments) {
  // look for the command line key in the comments
  for (const auto& comment : comments) {
    if (comment.rfind(kTsvCommandLineKey + kTsvCommentKeyValueSeparator, 0) == 0) {
      return comment.substr(kTsvCommandLineKey.size() + kTsvCommentKeyValueSeparator.size());
    }
  }
  return std::nullopt;
}

/**
 * @brief Extract the version from the comments that is formatted for VCF with an upper-case key value to satisfy more
 * strict parsers.
 * @param comments The vector of comments.
 * @return An optional string containing the version if found, otherwise std::nullopt.
 */
std::optional<std::string> GetVcfVersion(const Comments& comments) {
  // look for the version key in the comments
  auto version = GetVersion(comments);
  if (version.has_value()) {
    return fmt::format("##VERSION={}", version.value());
  }
  return std::nullopt;
}

/**
 * @brief Extract the command line from the comments that is formatted for VCF with an upper-case key value to satisfy
 * more strict parsers.
 * @param comments The vector of comments.
 * @return An optional string containing the command line if found, otherwise std::nullopt.
 */
std::optional<std::string> GetVcfCommandLine(const Comments& comments) {
  auto command_line = GetCommandLine(comments);
  if (command_line.has_value()) {
    return fmt::format("##COMMANDLINE={}", command_line.value());
  }
  return std::nullopt;
}

}  // namespace xoos::io
