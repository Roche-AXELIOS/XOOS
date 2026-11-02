#include "io/fai.h"

#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <xoos/log/logging.h>
#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

namespace xoos::cnc {

SeqLenMap ParseFai(const fs::path& fname) {
  SeqLenMap ret;

  if (!fs::exists(fname)) {
    Logging::Error("ParseFai: fai file does not exist");
    throw std::runtime_error("Fai file does not exist: " + fname.string());
  }

  if (fs::file_size(fname) == 0) {
    Logging::Error("ParseFai: fai file is empty");
    throw std::runtime_error("Empty fai file: " + fname.string());
  }

  std::ifstream ifs(fname);
  std::string line;
  while (std::getline(ifs, line)) {
    std::istringstream is(line);
    std::string key;
    size_t val = 0;
    is >> key;
    is >> val;
    // this will happen if >>key or >>val fails
    if (!is) {
      Logging::Error("ParseFai: Either found a 0 length contig or is invalid .fai file");
      throw std::runtime_error("Invalid fai file: " + fname.string());
    }
    ret[key] = val;
  }
  return ret;
}

StrUnorderedMap<size_t> GetContigOrder(const fs::path& fname) {
  std::ifstream ifs(fname);
  std::string line;
  // this will hold the correct ordering of the contigs, determined by fai_file
  StrUnorderedMap<size_t> contig_order;
  size_t i = 0;
  while (std::getline(ifs, line)) {
    std::istringstream is(line);
    std::string key;
    is >> key;
    contig_order[key] = i++;
  }
  return contig_order;
}

/**
 * @brief return a std::vector of contig names in the order they appear in the fai file
 * @param fname path to fai file
 */
vec<std::string> GetContigsInOrder(const fs::path& fname) {
  std::ifstream ifs(fname);
  std::string line;
  // this will hold the correct ordering of the contigs, determined by fai_file
  vec<std::string> contigs_in_order;
  while (std::getline(ifs, line)) {
    std::istringstream is(line);
    std::string contig;
    is >> contig;
    contigs_in_order.emplace_back(contig);
  }
  return contigs_in_order;
}
}  // namespace xoos::cnc
