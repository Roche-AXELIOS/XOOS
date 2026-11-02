#pragma once
#include <cstddef>
#include <string>
#include <utility>

#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

extern "C" {
#include <libbigwig/bigWig.h>
}

#include "io/bigwig.h"

namespace xoos::cnc {

// Forward declaration
struct Observations;

/**
 * @brief RAII wrapper owning a libbigwig write-mode file handle.
 *
 * Owns the `bigWigFile_t*` opened for writing and releases it via BigWigPtr on
 * scope exit. Each member wraps a libbigwig write step and throws
 * xoos::error::Error on a libbigwig error return.
 */
class BigWigWriter {
 public:
  BigWigWriter() = delete;

  /**
   * @brief Open @p path for writing.
   * @throws xoos::error::Error if the path is not writable or the handle cannot be opened.
   */
  explicit BigWigWriter(const fs::path& path);

  BigWigWriter(const BigWigWriter&) = delete;
  BigWigWriter(BigWigWriter&&) noexcept = default;
  ~BigWigWriter() = default;

  BigWigWriter& operator=(const BigWigWriter&) = delete;
  BigWigWriter& operator=(BigWigWriter&&) noexcept = default;

  /**
   * @brief Set the chromosome list on the write handle.
   * @param chrom_names chromosome names (C strings), same size as @p chrom_lengths
   * @param chrom_lengths chromosome lengths, same size as @p chrom_names
   * @throws xoos::error::Error on size mismatch or if the list cannot be created.
   */
  void SetChromList(const vec<const char*>& chrom_names, const vec<u32>& chrom_lengths);

  /**
   * @brief Create and write the bigWig header.
   * @param zoom_levels number of zoom (reduced-resolution) levels to store
   * @throws xoos::error::Error if the header cannot be created or written.
   */
  void WriteHeader(s32 zoom_levels);

  /**
   * @brief Append @p count intervals to the bigWig.
   * @param chroms chromosome names (C strings) per interval
   * @param starts 0-based start positions per interval
   * @param ends end positions per interval
   * @param values values per interval
   * @param count number of intervals to write (no-op when 0)
   * @throws xoos::error::Error if any input vector is smaller than @p count or the write fails.
   */
  void AddIntervals(const vec<const char*>& chroms,
                    const vec<u32>& starts,
                    const vec<u32>& ends,
                    const vec<f32>& values,
                    size_t count);

 private:
  BigWigPtr _fp;
};

/**
 * @brief Reads chromosome lengths from a FASTA index (.fai) file
 *
 * This function parses a FASTA index file to extract chromosome names and their
 * corresponding lengths. The .fai file format contains tab-separated values where
 * the first column is the chromosome name and the second column is its length.
 *
 * @param fai_filename fs::path to the FASTA index (.fai) file to read
 * @return std::pair containing:
 *         - first: A map with chromosome names as keys and their lengths as values
 *         - second: A std::vector of chromosome names in the order they appear in the file
 * @throws std::runtime_error if the file cannot be opened or parsed
 */
std::pair<StrMap<u32>, vec<std::string>> GetChromLengthsFromFaiFile(const std::string& fai_filename);

/**
 * @brief Validates that genomic intervals are sorted according to chromosome order and genomic positions.
 *
 * This function checks whether a collection of genomic intervals is properly sorted first by
 * chromosome order (as specified in the reference genome index) and then by genomic position
 * within each chromosome. The validation ensures that intervals can be processed sequentially
 * for formats that require sorted input (e.g., BigWig files).
 *
 * @param chroms Vector of chromosome names as C-strings for each interval
 * @param starts Vector of start positions (0-based) corresponding to each interval
 * @param chrom_order Vector defining the canonical chromosome ordering (typically from .fai file)
 *
 * @return true if intervals are properly sorted, false otherwise
 *
 * @details The function performs the following validation steps:
 * - Returns true immediately for empty input
 * - Creates a mapping from chromosome names to their canonical order indices
 * - Iterates through consecutive interval pairs checking:
 *   - Both chromosomes exist in the reference order
 *   - Current chromosome order >= previous chromosome order
 *   - For same chromosome: current position >= previous position
 * - Logs detailed error messages for any ordering violations
 *
 * @note This function assumes chroms and starts std::vectors have the same size
 * @warning Will return false and log errors if any chromosome is not found in chrom_order
 */
bool ValidateIntervalsSorted(const vec<const char*>& chroms,
                             const vec<u32>& starts,
                             const vec<std::string>& chrom_order);

/**
 * @brief Write the observations to a BigWig file
 * @param obs The observations object to write
 * @param filename path to output BigWig file
 * @param fai_filename path to FASTA index file containing chromosome lengths
 */
void WriteObservationsToBigWig(const Observations& obs, const fs::path& filename, const fs::path& fai_filename);

/**
 * @brief Write a valid but data-empty BigWig file (header + chromosome list, no intervals).
 *
 * Used by the low-coverage early-exit path to produce structurally valid BigWig
 * files that downstream consumers (IGV, pyBigWig) can open without error.
 *
 * @param filename path to output BigWig file
 * @param fai_filename path to FASTA index file containing chromosome lengths
 */
void WriteEmptyBigWig(const fs::path& filename, const fs::path& fai_filename);

}  // namespace xoos::cnc
