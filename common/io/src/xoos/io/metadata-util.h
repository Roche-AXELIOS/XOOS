#pragma once
#include <optional>
#include <string>
#include <string_view>

#include <xoos/io/bam/read-group.h>

namespace xoos::io {

const std::string kTsvCommentLinePrefix = "#";
const std::string kVcfCommentLinePrefix = "##";
const std::string kRocheCommandLinePrefix = "##RocheCommandLine=<";
const std::string kKeyId = "ID";
const std::string kKeyVersion = "Version";
const std::string kKeyCommandLine = "CommandLine";

const std::string kSampleLinePrefix = "##Sample=<";
const std::string kKeyName = "Name";
const std::string kKeyLibraryId = "LibraryId";
const std::string kKeyPlatform = "Platform";
const std::string kKeyPlatformModel = "PlatformModel";

/**
 * @brief Holds program name, version, and full command line for metadata output.
 */
struct CommandLineInfo {
  std::string name;
  std::string version;
  std::string command_line;
};

/**
 * @brief Format a CommandLineInfo as a RocheCommandLine metadata line.
 *
 * Produces: ##RocheCommandLine=<ID=name,Version="version",CommandLine="command_line">
 *
 * Used for both VCF and TSV output.
 */
std::string GetCommandInfo(const CommandLineInfo& info);

/**
 * @brief Parse a RocheCommandLine metadata line back into a CommandLineInfo struct.
 *
 * Expected format:
 *   ##RocheCommandLine=<ID=name,Version="version",CommandLine="command_line">
 *
 * Key-value pairs inside the angle brackets are parsed in any order.
 * Handles escaped quotes and backslashes in quoted values (inverse of EscapeMetadataValue).
 * Returns std::nullopt silently if the line is not a RocheCommandLine line, or with a
 * warning log if the line starts with ##RocheCommandLine= but is malformed.
 */
std::optional<CommandLineInfo> ParseCommandInfo(std::string_view line);

/**
 * @brief Format the sample information from a ReadGroup as a Sample metadata line.
 *
 * Produces: ##Sample=<Name=SAMPLE_NAME,LibraryId=LIBRARY_ID,Platform=PLATFORM,PlatformModel=PLATFORM_MODEL>
 *
 * Fields are populated from the read group's sample name (SM), library (LB),
 * platform (PL), and platform model (PM) tags respectively. Keys whose value is
 * empty are omitted, so a read group with only a sample name renders as
 * ##Sample=<Name=SAMPLE_NAME>.
 * Values containing a comma (the field delimiter) or a double-quote are wrapped
 * in double quotes so they remain unambiguous; embedded double-quotes and
 * backslashes are escaped.
 *
 * Used to annotate TSV metrics files with the sample they were computed from.
 */
std::string GetSampleInfo(const ReadGroup& read_group);

/**
 * @brief Bundles the metadata written to the comment header of a TSV file.
 *
 * The read group is optional because not all modules have sample-level information.
 */
struct TsvMetadata {
  CommandLineInfo command_line_info;
  std::optional<ReadGroup> read_group;
};

/**
 * @brief Write a RocheCommandLine metadata comment to an output stream.
 *
 * Writes a single ##RocheCommandLine=<...> comment line directly to the stream
 * to avoid CSV quoting of special characters (commas, double-quotes).
 */
template <typename OutputStream>
void WriteTsvMetadata(OutputStream& stream, const CommandLineInfo& info) {
  stream << GetCommandInfo(info) << "\n";
}

/**
 * @brief Write the TSV comment header for the given metadata to an output stream.
 *
 * Writes the ##RocheCommandLine=<...> line, followed by a ##Sample=<...> line
 * when sample information is present. Both lines are written directly to the
 * stream to avoid CSV quoting of special characters (commas, double-quotes).
 */
template <typename OutputStream>
void WriteTsvMetadata(OutputStream& stream, const TsvMetadata& metadata) {
  stream << GetCommandInfo(metadata.command_line_info) << "\n";
  if (metadata.read_group.has_value()) {
    stream << GetSampleInfo(*metadata.read_group) << "\n";
  }
}

}  // namespace xoos::io
