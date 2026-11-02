#include "io/write-bigwig.h"

#include <cstddef>
#include <exception>
#include <fstream>
#include <ios>
#include <limits>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/adopt.h>

#include "io/bigwig.h"
#include "observations/observations.h"

extern "C" {
#include <libbigwig/bigWig.h>
}

namespace xoos::cnc {

namespace {
// Seven zoom levels balance file size and reduced-resolution detail for observation tracks.
constexpr s32 kObservationZoomLevels = 7;
// An empty bigWig carries no intervals, so no zoom summaries are needed.
constexpr s32 kEmptyZoomLevels = 0;
}  // namespace

BigWigWriter::BigWigWriter(const fs::path& path) {
  // Probe with append mode because libbigwig cannot safely handle a write-open failure. Append mode
  // preserves an existing file while verifying that the destination is writable.
  std::error_code exists_ec;
  const bool existed_before = fs::exists(path, exists_ec);
  if (exists_ec) {
    throw error::Error("Failed to inspect BigWig output path '{}': {}", path.string(), exists_ec.message());
  }
  if (const std::ofstream probe(path, std::ios::app); !probe) {
    throw error::Error("Failed to open BigWig file for writing: {}", path.string());
  }

  // Remove an empty probe artifact without deleting a pre-existing destination.
  const auto remove_probe_artifact = [&existed_before, &path]() {
    if (!existed_before) {
      std::error_code remove_ec;
      fs::remove(path, remove_ec);
    }
  };

  _fp = AdoptUnique<bigWigFile_t, BigWigDeleter>(bwOpen(path.c_str(), nullptr, "w"));  // NOSONAR(M23_286): C adopt
  if (!_fp) {
    remove_probe_artifact();
    throw error::Error("Failed to open BigWig file for writing: {}", path.string());
  }
}

void BigWigWriter::SetChromList(const vec<const char*>& chrom_names, const vec<u32>& chrom_lengths) {
  if (chrom_names.size() != chrom_lengths.size()) {
    throw error::Error("Chromosome name and length vectors must have the same size ({} vs {})",
                       chrom_names.size(),
                       chrom_lengths.size());
  }
  _fp->cl = bwCreateChromList(chrom_names.data(), chrom_lengths.data(), static_cast<s64>(chrom_names.size()));
  if (_fp->cl == nullptr) {
    throw error::Error("Failed to create chromosome list");
  }
}

void BigWigWriter::WriteHeader(const s32 zoom_levels) {  // NOSONAR(S5817)
  if (bwCreateHdr(_fp.get(), zoom_levels) != 0) {
    throw error::Error("Failed to create BigWig header");
  }
  if (bwWriteHdr(_fp.get()) != 0) {
    throw error::Error("Failed to write BigWig header");
  }
}

void BigWigWriter::AddIntervals(const vec<const char*>& chroms,  // NOSONAR(S5817)
                                const vec<u32>& starts,
                                const vec<u32>& ends,
                                const vec<f32>& values,
                                const size_t count) {
  if (count == 0) {
    return;
  }
  if (chroms.size() < count || starts.size() < count || ends.size() < count || values.size() < count) {
    throw error::Error(
        "Interval vectors are smaller than the requested count ({}): chroms={}, starts={}, ends={}, values={}",
        count,
        chroms.size(),
        starts.size(),
        ends.size(),
        values.size());
  }
  // libbigwig's bwAddIntervals takes the interval count as a uint32_t. Guard the narrowing conversion
  // so an oversized count fails loudly instead of silently truncating and writing fewer intervals.
  if (count > std::numeric_limits<u32>::max()) {
    throw error::Error(
        "Too many intervals for a single BigWig write ({}); maximum is {}", count, std::numeric_limits<u32>::max());
  }
  if (bwAddIntervals(_fp.get(), chroms.data(), starts.data(), ends.data(), values.data(), static_cast<u32>(count)) !=
      0) {
    throw error::Error("Failed to add intervals to BigWig");
  }
}

std::pair<StrMap<u32>, vec<std::string>> GetChromLengthsFromFaiFile(const std::string& fai_filename) {
  StrMap<u32> chrom_lengths;
  // Preserve .fai order
  vec<std::string> chrom_order;
  std::ifstream fai_file(fai_filename);

  if (!fai_file.is_open()) {
    throw error::Error("Could not open .fai file: {}", fai_filename);
  }

  std::string line;
  while (std::getline(fai_file, line)) {
    std::istringstream iss(line);
    std::string chrom_name;
    u32 length = 0;

    if (iss >> chrom_name >> length) {
      chrom_lengths[chrom_name] = length;
      // Preserve original order
      chrom_order.push_back(chrom_name);
    }
  }

  return {chrom_lengths, chrom_order};
}

bool ValidateIntervalsSorted(const vec<const char*>& chroms,
                             const vec<u32>& starts,
                             const vec<std::string>& chrom_order) {
  if (chroms.empty()) {
    return true;
  }

  // Map chromosomes to their positions in the .fai file.
  StrMap<s32> chrom_order_map;
  for (size_t i = 0; i < chrom_order.size(); ++i) {
    chrom_order_map[chrom_order[i]] = static_cast<s32>(i);
  }

  for (size_t i = 1; i < chroms.size(); ++i) {
    std::string prev_chrom(chroms[i - 1]);
    std::string curr_chrom(chroms[i]);

    const auto prev_it = chrom_order_map.find(prev_chrom);
    const auto curr_it = chrom_order_map.find(curr_chrom);
    if (prev_it == chrom_order_map.end() || curr_it == chrom_order_map.end()) {
      Logging::Error("Unknown chromosome found: {} or {}", prev_chrom, curr_chrom);
      return false;
    }

    s32 prev_order = prev_it->second;
    s32 curr_order = curr_it->second;

    if (prev_order > curr_order) {
      Logging::Error("Chromosomes out of order at interval {}: {} (order {}) comes after {} (order {})",
                     i,
                     curr_chrom,
                     curr_order,
                     prev_chrom,
                     prev_order);
      return false;
    }
    if (prev_order == curr_order) {
      if (starts[i] < starts[i - 1]) {
        Logging::Error(
            "Positions out of order at interval {} on {}: {} comes after {}", i, curr_chrom, starts[i], starts[i - 1]);
        return false;
      }
    }
  }

  return true;
}

void WriteObservationsToBigWig(const Observations& obs, const fs::path& filename, const fs::path& fai_filename) {
  // Read chromosome information from .fai file
  StrMap<u32> chrom_lengths;
  vec<std::string> chrom_order;
  try {
    auto [cl, co] = GetChromLengthsFromFaiFile(fai_filename.string());
    chrom_lengths = std::move(cl);
    chrom_order = std::move(co);
    Logging::Info("Read {} chromosomes from {}", chrom_lengths.size(), fai_filename.string());
  } catch (const std::exception& e) {
    throw error::Error("Error reading .fai file '{}': {}", fai_filename.string(), e.what());
  }

  {
    // Closing the writer finalizes the file before success is reported.
    BigWigWriter writer(filename);

    // Store chromosome names in .fai order before taking C-string pointers.
    vec<std::string> chrom_name_strings;
    chrom_name_strings.reserve(chrom_order.size());

    for (const std::string& chrom_name : chrom_order) {
      chrom_name_strings.push_back(chrom_name);
    }

    // Create the C arrays after string storage is stable.
    vec<const char*> chrom_names_c;
    vec<u32> chrom_lens_c;
    for (const auto& chrom_name : chrom_name_strings) {
      chrom_names_c.push_back(chrom_name.c_str());
      chrom_lens_c.push_back(chrom_lengths[chrom_name]);
    }
    writer.SetChromList(chrom_names_c, chrom_lens_c);

    // Create and write header
    writer.WriteHeader(kObservationZoomLevels);

    // Store all chromosome names before taking C-string pointers.
    vec<std::string> chrom_strings;
    chrom_strings.reserve(obs.contigs.size());
    for (const auto& contig : obs.contigs) {
      chrom_strings.push_back(contig);
    }

    // Create interval arrays after string storage is stable.
    vec<const char*> chroms;
    vec<u32> starts_vec;
    vec<u32> ends_vec;
    vec<f32> values_vec;

    for (size_t i = 0; i < obs.contigs.size(); ++i) {
      if (!chrom_lengths.contains(chrom_strings[i])) {
        throw error::Error("Chromosome '{}' not found in .fai file '{}'", chrom_strings[i], fai_filename.string());
      }
      const auto idx = static_cast<Eigen::Index>(i);
      chroms.push_back(chrom_strings[i].c_str());
      starts_vec.push_back(static_cast<u32>(obs.starts[idx]));
      ends_vec.push_back(static_cast<u32>(obs.ends[idx]));
      values_vec.push_back(static_cast<f32>(obs.obvs[idx]));
    }

    // Validate sorting before writing
    Logging::Info("Validating interval sorting for bigWig writing...");
    if (!ValidateIntervalsSorted(chroms, starts_vec, chrom_order)) {
      throw error::Error("Intervals are not properly sorted for BigWig format");
    }

    // Write intervals to BigWig
    writer.AddIntervals(chroms, starts_vec, ends_vec, values_vec, obs.contigs.size());
  }

  Logging::Info("BigWig file '{}' created successfully with {} intervals", filename.string(), obs.contigs.size());
}

void WriteEmptyBigWig(const fs::path& filename, const fs::path& fai_filename) {
  auto [chrom_lengths, chrom_order] = GetChromLengthsFromFaiFile(fai_filename.string());

  // The writer owns the file handle and releases it on scope exit, including if any step below throws.
  BigWigWriter writer(filename);

  vec<const char*> chrom_names_c;
  vec<u32> chrom_lens_c;
  for (const auto& name : chrom_order) {
    chrom_names_c.push_back(name.c_str());
    chrom_lens_c.push_back(chrom_lengths[name]);
  }

  writer.SetChromList(chrom_names_c, chrom_lens_c);
  writer.WriteHeader(kEmptyZoomLevels);
}

}  // namespace xoos::cnc
