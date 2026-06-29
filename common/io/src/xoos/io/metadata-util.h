#pragma once
#include <optional>
#include <string>
#include <vector>

#include <csv.hpp>

namespace xoos::io {
using Comments = std::vector<std::string>;

void AddComment(Comments& comments,
                const std::string& comment);  // now comments has {"a", "b", "c", "x=y"}

void AddCommentKeyPair(Comments& comments,
                       const std::string& key,
                       const std::string& value);  // now comments has {"a", "b", "c", "x=y"}

void AddVersionAndCommandLineComment(Comments& comments, const std::string& version, const std::string& command_line);

const std::string kTsvCommentLinePrefix = "#";
const std::string kVcfCommentLinePrefix = "##";

/**
 * @brief Write comments to a TSV writer, prefixing each comment with a comment character.
 * @tparam TSVWriter The type of the TSV writer, which must support the '<<' operator for writing rows.
 * @param writer The TSV writer to which comments will be written.
 * @param comments A vector of comment strings to be written.
 */
template <typename TSVWriter>
void WriteTsvComments(TSVWriter& writer, const Comments& comments) {
  for (const auto& comment : comments) {
    writer << std::vector{kTsvCommentLinePrefix + comment};
  }
}

/**
 * @brief Generate a string representation for the version key from comments.
 */
std::optional<std::string> GetVersion(const Comments& comments);

/**
 * @brief Generate a string representation for the command line key from comments.
 */
std::optional<std::string> GetCommandLine(const Comments& comments);

// helper functions to extract version and command line from comments
std::optional<std::string> GetVcfVersion(const Comments& comments);
std::optional<std::string> GetVcfCommandLine(const Comments& comments);

}  // namespace xoos::io
