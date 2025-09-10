#include "xoos/io/metadata-util.h"

namespace xoos::io {

const std::string kTsvCommentLinePrefix = "#";
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

// we write TSV comments to the writer, prepending with a comment
void WriteTsvComments(csv::TSVWriter<std::ofstream, false> writer, const Comments& comments) {
  for (const auto& comment : comments) {
    writer << std::vector{kTsvCommentLinePrefix + comment};
  }
}

}  // namespace xoos::io
