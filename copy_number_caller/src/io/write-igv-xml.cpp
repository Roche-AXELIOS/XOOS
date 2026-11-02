#include "io/write-igv-xml.h"

#include <fstream>
#include <optional>
#include <stdexcept>

#include <xoos/log/logging.h>
#include <xoos/types/fs.h>

namespace xoos::cnc {
void WriteIGVXML(const fs::path& output_xml_path,
                 const std::optional<fs::path>& logr_file,
                 const std::optional<fs::path>& bigwig_file) {
  if (!logr_file.has_value() && !bigwig_file.has_value()) {
    Logging::Error("No LogR or BAF bigWig files generated, skipping IGV XML writing");
    throw std::runtime_error("No LogR or BAF bigWig files generated. IGV XML cannot be written.");
  }
  if (!logr_file.has_value()) {
    Logging::Warn("No LogR bigWig file generated, IGV XML will be written without LogR track");
  }
  if (!bigwig_file.has_value()) {
    Logging::Warn("No BAF bigWig file generated, IGV XML will be written without BAF track");
  }
  std::ofstream xml_file(output_xml_path);
  if (!xml_file.is_open()) {
    throw std::runtime_error("Could not open XML output file for writing: " + output_xml_path.string());
  }
  // note: filesystem::path objects are automatically quoted when streamed:
  // https://en.cppreference.com/w/cpp/filesystem/path.html
  // header
  xml_file << R"(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)" << "\n"
           << R"(<Session version="8">)"
           << "\n"
           // start Resources section
           << R"(    <Resources>)" << "\n";
  if (logr_file.has_value()) {
    xml_file << R"(        <Resource path=)" << logr_file.value().filename() << R"( type="bw"/>)" << "\n";
  }
  if (bigwig_file.has_value()) {
    xml_file << R"(        <Resource path=)" << bigwig_file.value().filename() << R"( type="bw"/>)" << "\n";
  }
  xml_file << R"(    </Resources>)" << "\n" << R"(    <Panel name="DataPanel">)" << "\n";
  // start Tracks sections
  if (logr_file.has_value()) {
    xml_file << R"(        <Track attributeKey=)" << logr_file.value().filename()
             << R"( autoScale="true" clazz="org.broad.igv.track.DataSourceTrack" id=)" << logr_file.value().filename()
             << R"( name=)" << logr_file.value().filename()
             << R"( renderer="SCATTER_PLOT" visible="true" windowFunction="none">)" << "\n"
             << R"(            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" type="LINEAR"/>)" << "\n"
             << R"(        </Track>)" << "\n";
  }
  if (bigwig_file.has_value()) {
    xml_file
        << R"(        <Track attributeKey=)" << bigwig_file.value().filename()
        << R"( autoScale="true" clazz="org.broad.igv.track.DataSourceTrack" id=)" << bigwig_file.value().filename()
        << R"( name=)" << bigwig_file.value().filename()
        << R"( renderer="SCATTER_PLOT" visible="true" windowFunction="none">)" << "\n"
        << R"(            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="1.0" minimum="0.0" type="LINEAR"/>)"
        << "\n"
        << R"(        </Track>)" << "\n";
  }
  xml_file << R"(    </Panel>)" << "\n"
           << R"(</Session>)" << "\n";
  xml_file.close();
}

}  // namespace xoos::cnc
