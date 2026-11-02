#include "detect-sample-type.h"

#include <limits>
#include <mutex>
#include <string>
#include <utility>

#include <fmt/format.h>
#include <htslib/vcf.h>

#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/for_each.hpp>

#include <xoos/error/error.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/log/logging.h>
#include <xoos/types/vec.h>

#include "compute-vcf-features/vcf-header-util.h"
#include "file_type_inference/sample-type-metrics.h"
#include "util/region-util.h"

namespace xoos::svc {

namespace {
constexpr s32 kDetectSampleTypeUnpackLevel =
    static_cast<s32>(static_cast<u32>(BCF_UN_STR) | static_cast<u32>(BCF_UN_FMT));

/**
 * @brief Indexed VCF region assigned to one sample-type detection worker.
 * @details The region is represented with htslib contig index and half-open coordinates so each worker can
 * call VcfReader::SetRegion() directly without shared reader state.
 */
struct DetectSampleTypeTask {
  /// Chromosome name, retained for logging and target interval lookup.
  std::string chrom;
  /// htslib/VCF header contig index used for region queries.
  s32 chrom_index{};
  /// Inclusive zero-based region start.
  u64 start{};
  /// Exclusive zero-based region end.
  u64 end{};
};

/**
 * @brief Prepared sample-type detection inputs shared by sequential and parallel scans.
 * @details Separates setup from scanning so sample index selection, denominator length calculation, and
 * region task construction happen exactly once before choosing the execution mode.
 */
struct DetectSampleTypeSetup {
  /// Zero-based VCF sample index evaluated for sample-type metrics.
  s32 sample_idx{};
  /// Denominator length used by FinalizeSampleTypeMetrics().
  u64 total_length{};
  /// True when sample-type detection should avoid indexed parallel queries.
  bool use_single_thread{};
  /// Per-region work items used only when use_single_thread is false.
  vec<DetectSampleTypeTask> tasks;
};

class TargetRegionCursor {
 public:
  /**
   * @brief Create a cursor over BED intervals grouped by chromosome.
   * @details The cursor stores a per-chromosome next-interval position so repeated calls can advance
   * through sorted VCF records without re-scanning intervals from the beginning.
   * @param targets BED-style intervals grouped by chromosome. The referenced map must outlive the cursor.
   */
  explicit TargetRegionCursor(const ChromIntervalsMap& targets) : _targets(&targets) {
  }

  /**
   * @brief Test whether a VCF record overlaps the target intervals.
   * @details Advances the cursor for the record chromosome until the next interval could overlap the
   * record position. This is efficient for records visited in coordinate order.
   * @param record VCF record to test.
   * @return True when the record position is inside a target interval; otherwise false.
   */
  bool Overlaps(const io::VcfRecordPtr& record) {
    const auto chrom = record->Chromosome();
    const auto intervals_it = _targets->find(chrom);
    if (intervals_it == _targets->end()) {
      return false;
    }

    const auto pos = static_cast<u64>(record->Position());
    const auto& intervals = intervals_it->second;
    size_t& next_region_idx = _next_region_by_chromosome[chrom];
    while (next_region_idx < intervals.size() && intervals.at(next_region_idx).end <= pos) {
      ++next_region_idx;
    }
    if (next_region_idx >= intervals.size()) {
      return false;
    }
    const auto& interval = intervals.at(next_region_idx);
    return interval.start <= pos && pos < interval.end;
  }

 private:
  const ChromIntervalsMap* _targets;
  StrMap<size_t> _next_region_by_chromosome;
};

/**
 * @brief Build per-contig work items for indexed sample-type detection.
 * @details When target regions are present, each task spans the first through last target interval for a
 * chromosome. Without targets, each task spans the full contig length from the VCF header, falling back to
 * the maximum signed htslib coordinate when the contig length is unavailable.
 * @param reader Open VCF reader used to access contig indexes and header lengths.
 * @param target_regions Optional target intervals limiting task ranges.
 * @return Task list suitable for parallel region queries.
 */
vec<DetectSampleTypeTask> BuildDetectSampleTypeTasks(const io::VcfReader& reader,
                                                     const std::optional<ChromIntervalsMap>& target_regions) {
  vec<DetectSampleTypeTask> tasks;
  const auto contig_indexes = reader.GetContigIndexes();
  const auto contig_lengths = reader.GetHeader()->GetContigLengths();
  if (target_regions.has_value()) {
    for (const auto& [chrom, intervals] : target_regions.value()) {
      if (!contig_indexes.contains(chrom) || intervals.empty()) {
        continue;
      }
      // Query the full target span for the chromosome; TargetRegionCursor still filters individual
      // records to the actual intervals while avoiding one htslib query per BED row.
      tasks.emplace_back(chrom, contig_indexes.at(chrom), intervals.front().start, intervals.back().end);
    }
    return tasks;
  }

  for (const auto& [chrom, chrom_index] : contig_indexes) {
    const u64 end = contig_lengths.contains(chrom) && contig_lengths.at(chrom) > 0
                        ? contig_lengths.at(chrom)
                        : static_cast<u64>(std::numeric_limits<s64>::max());
    tasks.emplace_back(chrom, chrom_index, 0, end);
  }
  return tasks;
}

/**
 * @brief Tally sample-type observations from an already-open sequential VCF reader.
 * @details Iterates records in reader order, optionally applies a target BED cursor, and adds passing
 * records to @p accumulator. This helper backs the single-threaded path.
 * @param reader Open VCF reader positioned at the start of the scan.
 * @param sample_idx Zero-based VCF sample index used for FORMAT fields.
 * @param max_af Maximum allele frequency for low-frequency observations.
 * @param has_bed Whether @p target_regions_by_chromosome should be applied.
 * @param target_regions_by_chromosome Target intervals grouped by chromosome; ignored when @p has_bed is false.
 * @param accumulator Accumulator updated in place.
 */
void TallyRecords(const io::VcfReader& reader,
                  const s32 sample_idx,
                  const f64 max_af,
                  const bool has_bed,
                  const ChromIntervalsMap& target_regions_by_chromosome,
                  SampleTypeMetricAccumulator& accumulator) {
  std::optional<TargetRegionCursor> target_region_cursor;
  if (has_bed) {
    target_region_cursor.emplace(target_regions_by_chromosome);
  }
  while (const auto& record = reader.GetNextRecord(kDetectSampleTypeUnpackLevel)) {
    if (target_region_cursor.has_value() && !target_region_cursor->Overlaps(record)) {
      continue;
    }
    TallySampleTypeRecord(record, sample_idx, max_af, accumulator);
  }
}

/**
 * @brief Run sample-type metric accumulation in one sequential VCF pass.
 * @details Opens a VCF reader and delegates record filtering and metric updates to TallyRecords().
 * @param vcf_file Path to the VCF to scan.
 * @param sample_idx Zero-based VCF sample index used for FORMAT fields.
 * @param max_af Maximum allele frequency for low-frequency observations.
 * @param has_bed Whether target interval filtering is active.
 * @param target_regions_by_chromosome Target intervals grouped by chromosome; ignored when @p has_bed is false.
 * @return Accumulated sample-type observations from the full scan.
 */
SampleTypeMetricAccumulator TallyRecordsSequential(const fs::path& vcf_file,
                                                   const s32 sample_idx,
                                                   const f64 max_af,
                                                   const bool has_bed,
                                                   const ChromIntervalsMap& target_regions_by_chromosome) {
  const io::VcfReader reader(vcf_file);
  SampleTypeMetricAccumulator accumulator;
  TallyRecords(reader, sample_idx, max_af, has_bed, target_regions_by_chromosome, accumulator);
  return accumulator;
}

/**
 * @brief Run sample-type metric accumulation for one indexed region task.
 * @details Opens a thread-local VCF reader, seeks to @p task, applies optional target filtering, and
 * accumulates only records from that region. It returns an empty accumulator when htslib cannot set the
 * requested region.
 * @param vcf_file Path to the indexed VCF to scan.
 * @param task Chromosome region assigned to this worker.
 * @param sample_idx Zero-based VCF sample index used for FORMAT fields.
 * @param max_af Maximum allele frequency for low-frequency observations.
 * @param has_bed Whether target interval filtering is active.
 * @param target_regions_by_chromosome Target intervals grouped by chromosome; ignored when @p has_bed is false.
 * @return Accumulated observations from the assigned region.
 */
SampleTypeMetricAccumulator TallyRecordsForTask(const fs::path& vcf_file,
                                                const DetectSampleTypeTask& task,
                                                const s32 sample_idx,
                                                const f64 max_af,
                                                const bool has_bed,
                                                const ChromIntervalsMap& target_regions_by_chromosome) {
  SampleTypeMetricAccumulator accumulator;
  io::VcfReader reader(vcf_file);
  if (!reader.SetRegion(task.chrom_index, task.start, task.end)) {
    return accumulator;
  }
  std::optional<TargetRegionCursor> target_region_cursor;
  if (has_bed) {
    target_region_cursor.emplace(target_regions_by_chromosome);
  }

  while (const auto& record = reader.GetNextRegionRecord(kDetectSampleTypeUnpackLevel)) {
    if (target_region_cursor.has_value() && !target_region_cursor->Overlaps(record)) {
      continue;
    }
    TallySampleTypeRecord(record, sample_idx, max_af, accumulator);
  }
  return accumulator;
}

/**
 * @brief Run sample-type metric accumulation across indexed VCF regions in parallel.
 * @details Dispatches one taskflow job per region task. Each worker uses its own VCF reader and accumulator;
 * results are merged under a mutex after the worker completes.
 * @param vcf_file Path to the indexed VCF to scan.
 * @param tasks Region tasks built from contigs or target extents.
 * @param sample_idx Zero-based VCF sample index used for FORMAT fields.
 * @param options Detector options supplying the thread count and AF threshold.
 * @param has_bed Whether target interval filtering is active.
 * @param target_regions_by_chromosome Target intervals grouped by chromosome; ignored when @p has_bed is false.
 * @return Combined observations from all region tasks.
 */
SampleTypeMetricAccumulator TallyRecordsParallel(const fs::path& vcf_file,
                                                 const vec<DetectSampleTypeTask>& tasks,
                                                 const s32 sample_idx,
                                                 const DetectSampleTypeOptions& options,
                                                 const bool has_bed,
                                                 const ChromIntervalsMap& target_regions_by_chromosome) {
  SampleTypeMetricAccumulator combined;
  std::mutex combined_mutex;
  tf::Executor executor(options.threads);
  tf::Taskflow flow;
  flow.for_each(tasks.begin(),
                tasks.end(),
                [&vcf_file, sample_idx, &options, has_bed, &target_regions_by_chromosome, &combined, &combined_mutex](
                    const DetectSampleTypeTask& task) {
                  auto task_accumulator = TallyRecordsForTask(
                      vcf_file, task, sample_idx, options.max_af, has_bed, target_regions_by_chromosome);
                  if (task_accumulator.error_count == 0 && task_accumulator.dps.empty()) {
                    return;
                  }
                  const std::scoped_lock lock{combined_mutex};
                  MergeSampleTypeAccumulators(combined, std::move(task_accumulator));
                })
      .name("Parallel sample-type detection VCF scan");
  executor.run(flow).get();
  return combined;
}

/**
 * @brief Prepare common inputs for sample-type detection.
 * @details Determines the sample index from tumor-normal VCF header metadata, falling back to sample 0
 * when the header lacks `tumor_sample`; computes the denominator length from BED intervals or VCF contigs;
 * and builds indexed region tasks when multithreading can be used.
 * @param vcf_file Path to the VCF to inspect.
 * @param target_regions Optional BED intervals limiting detection.
 * @param options Detection options used to decide whether to build parallel tasks.
 * @return Detection setup shared by sequential and parallel scans.
 */
DetectSampleTypeSetup BuildDetectSampleTypeSetup(const fs::path& vcf_file,
                                                 const std::optional<ChromIntervalsMap>& target_regions,
                                                 const DetectSampleTypeOptions& options) {
  io::VcfReader reader(vcf_file);
  const auto header = reader.GetHeader();
  const auto sample_indexes = GetTumorNormalSampleIndexes(header);
  s32 sample_idx = 0;
  if (sample_indexes.has_value()) {
    sample_idx = sample_indexes->tumor_sample_idx;
  } else {
    Logging::Warn("VCF header '{}' does not define tumor_sample; using sample index 0 for sample-type detection",
                  vcf_file);
  }

  u64 total_length =
      target_regions.has_value() ? SumChromIntervalLengths(*target_regions) : SumVcfHeaderContigLengths(header);
  if (total_length == 0 && !target_regions.has_value()) {
    total_length = kDefaultWholeGenomeLengthForSampleTypeDetection;
  }

  const bool use_single_thread = options.threads <= 1 || !reader.HasIndex();
  vec<DetectSampleTypeTask> tasks;
  if (!use_single_thread) {
    // Only indexed VCFs can be split into independent region queries. Unindexed input falls back to one
    // streaming pass so standalone detection still works for plain VCF files.
    tasks = BuildDetectSampleTypeTasks(reader, target_regions);
  }

  return {
      .sample_idx = sample_idx,
      .total_length = total_length,
      .use_single_thread = use_single_thread,
      .tasks = std::move(tasks),
  };
}
}  // namespace

SampleTypeMetrics ComputeSampleTypeMetrics(const fs::path& vcf_file,
                                           const std::optional<fs::path>& bed_file,
                                           const DetectSampleTypeOptions& options) {
  const auto target_regions = GetChromIntervalMap(bed_file);
  const auto setup = BuildDetectSampleTypeSetup(vcf_file, target_regions, options);

  const bool has_bed = target_regions.has_value() && !target_regions->empty();
  const ChromIntervalsMap empty_bed;
  const auto& target_regions_by_chromosome = has_bed ? target_regions.value() : empty_bed;
  SampleTypeMetricAccumulator accumulator;
  if (setup.use_single_thread) {
    accumulator =
        TallyRecordsSequential(vcf_file, setup.sample_idx, options.max_af, has_bed, target_regions_by_chromosome);
  } else {
    Logging::Info("Scanning VCF for sample-type metrics with {} threads across {} chromosomes",
                  options.threads,
                  setup.tasks.size());
    accumulator =
        TallyRecordsParallel(vcf_file, setup.tasks, setup.sample_idx, options, has_bed, target_regions_by_chromosome);
  }

  return FinalizeSampleTypeMetrics(accumulator, setup.total_length, vcf_file, options.max_af);
}

SampleType ClassifySampleType(const f64 phred, const f64 max_ffpe_phred) {
  using enum SampleType;
  return phred <= max_ffpe_phred ? kFfpe : kCellLine;
}

}  // namespace xoos::svc
