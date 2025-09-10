#pragma once
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

void WriteTsvComments(csv::TSVWriter<std::ofstream, false> writer, const Comments& comments);
}  // namespace xoos::io
