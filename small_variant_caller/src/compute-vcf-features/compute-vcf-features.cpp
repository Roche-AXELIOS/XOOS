#include "compute-vcf-features.h"

#include <algorithm>
#include <span>
#include <unordered_set>

#include <csv.hpp>

#include <taskflow/algorithm/for_each.hpp>

#include <xoos/error/error.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>

#include "core/config.h"
#include "core/filtering.h"
#include "core/variant-info-serializer.h"
#include "core/vcf-fields.h"
#include "util/file-util.h"
#include "util/locked-tsv-writer.h"
#include "util/log-util.h"
#include "util/seq-util.h"
#include "vcf-field-util.h"
#include "vcf-header-util.h"
#include "xoos/types/float.h"

namespace xoos::svc {

// Length of the sequence context for manual counting of tandem repeats after the variant's start position
static constexpr u64 kTandemRepeatSearchDistance{100};
// Upstream window size before variant's start position for calculating variant density
static constexpr u64 kVariantDensityUpstreamWindow{100};

/**
 * @brief Store the population allele frequency and quality for a variant.
 */
struct RefAltPopaf {
  std::string ref;
  std::string alt;
  f32 popaf;
  f32 qual;
};

/**
 * @brief Check if the VCF record has a negative position and warn if it does.
 * @param record VCF record to check
 * @return true if the position is negative, false otherwise
 * @note This function is used to ensure that all VCF records have valid positions.
 * @thows std::runtime_error if the position is negative and the warning is set as an error.
 */
static bool CheckAndWarnNegativePosition(const io::VcfRecordPtr& record) {
  if (record->Position() < 0) {
    WarnAsErrorIfSet("Found VCF record with position less than 0");
    return true;
  }
  return false;
}

/**
 * @brief Helper function to update popaf values of features at a given chromosome position.
 * @param feature_itr Iterator of vector of features in the same chromosome
 * @param feature_itr_end End iterator of vector of features
 * @param pos Chromosome position
 * @param ref_alt_af Map of reference and alt allele to allele frequency value for the position
 */
static void UpdatePopafAtPosition(std::span<VcfFeature>::iterator& feature_itr,
                                  const std::span<VcfFeature>::iterator& feature_itr_end,
                                  const u64 pos,
                                  const vec<RefAltPopaf>& ref_alt_af) {
  // Advance the iterator to the first feature at or after the given position.
  while (feature_itr != feature_itr_end && feature_itr->pos < pos) {
    ++feature_itr;
  }
  // Find all features at the given position and update their popaf values if a matching ref/alt allele is found.
  while (feature_itr != feature_itr_end && feature_itr->pos == pos) {
    f32 qual_found = -1;
    for (const auto& [ref, alt, popaf, qual] : ref_alt_af) {
      if (feature_itr->ref == ref && feature_itr->alt == alt && qual > qual_found) {
        qual_found = qual;
        feature_itr->popaf = popaf;
        // Do not break out of loop here because gnomAD VCF can potentially have >1 entries for the same variant, and
        // we only use the AF with higher QUAL.
      }
    }
    ++feature_itr;
  }
}

/**
 * @brief Block of VCF features for a specific chromosome.
 * @note This struct is used to process features in parallel for a specific chromosome.
 * Each block contains a span of VcfFeature objects, the chromosome index, and the chromosome length.
 * The span allows for efficient access to the features without copying them.
 * @note The `chrom_index` is the index of the chromosome in the VCF index,
 * and `chrom_length` is the length of the chromosome in base pairs.
 * @note The `features` span contains the VcfFeature objects for the chromosome.
 * It is expected that the features are sorted by position within the span.
 */
struct ChromFeatureBlock {
  std::span<VcfFeature> features;
  s32 chrom_index;
  u64 chrom_length;
};

/**
 * @brief Update the feature block with allele frequency values from a POPAF VCF file.
 * @param block Feature block to be updated
 * @param reader VCF reader for the POPAF VCF file
 */
static void UpdatePopafAtBlock(ChromFeatureBlock& block, io::VcfReader& reader) {
  auto feature_itr = block.features.begin();
  const auto feature_itr_end = block.features.end();

  const u64 start = std::max<u64>(0, feature_itr->pos);
  const u64 end = std::min<u64>((feature_itr_end - 1)->pos + 1, block.chrom_length);
  // +1 to the position of last feature because region is end-exclusive
  if (!reader.SetRegion(block.chrom_index, start, end)) {
    return;
  }
  u64 prev_pos = 0;
  vec<RefAltPopaf> ref_alt_af;
  while (const auto& record = reader.GetNextRegionRecord(BCF_UN_INFO)) {
    if (CheckAndWarnNegativePosition(record)) {
      continue;
    }
    const auto pos = static_cast<u64>(record->Position());
    if (pos < start || pos >= end) {
      // skip records outside the region of interest
      continue;
    }
    if (pos > prev_pos) {
      UpdatePopafAtPosition(feature_itr, feature_itr_end, prev_pos, ref_alt_af);
      if (feature_itr >= feature_itr_end) {
        break;
      }
      prev_pos = pos;
      ref_alt_af = {};
    }
    if (feature_itr->pos > pos) {
      continue;
    }
    const auto ref = record->Allele(0);
    const f32 qual = record->GetQuality(0);
    // start from 1 to skip the REF allele
    s32 allele_idx = 1;
    for (const auto af : record->GetInfoFieldNoCheck<f32>(kFieldAf)) {
      ref_alt_af.emplace_back(ref, record->Allele(allele_idx), af, qual);
      ++allele_idx;
    }
  }
  // update the last position
  UpdatePopafAtPosition(feature_itr, feature_itr_end, prev_pos, ref_alt_af);
}

// minimum number of features per block
constexpr size_t kMinBlockSize{100};

/**
 * @brief Get contig indexes and lengths from the VCF reader.
 * @param reader VCF reader
 * @return Tuple of contig indexes and lengths
 */
static std::tuple<ChromToIndexes, ChromToLengths> GetContigInfo(io::VcfReader& reader) {
  const auto& indexes = reader.GetContigIndexes();
  const auto& lengths = reader.GetHeader()->GetContigLengths();
  return std::make_tuple(indexes, lengths);
}

/**
 * @brief Get contig indexes and lengths from the VCF file.
 * @param vcf_path Path of the VCF file
 * @return Tuple of contig indexes and lengths
 */
static std::tuple<ChromToIndexes, ChromToLengths> GetContigInfo(const fs::path& vcf_path) {
  auto reader = io::VcfReader(vcf_path);
  return GetContigInfo(reader);
}

/**
 * @brief Update all features with allele frequency values from an indexed POPAF VCF.
 * @param all_features Map of chromosome to vector of features
 * @param vcf_path Path of POPAF VCF
 * @param threads Number of threads
 * @note VCF file must be indexed and contain the INFO/FORMAT field `AF` for population allele frequency.
 */
static void UpdatePopaf(StrMap<vec<VcfFeature>>& all_features, const fs::path& vcf_path, const u32 threads) {
  // This method annotates VCF feature entries with a population based allele
  // frequency value if a matching variant is found in a provided population af file.

  const auto& [chrom_indexes, chrom_lengths] = GetContigInfo(vcf_path);

  // Break up the set of VCF features from the sample to assign one feature block per task
  vec<ChromFeatureBlock> feat_blocks;
  for (auto& [chrom, features] : all_features) {
    if (features.empty() || !chrom_indexes.contains(chrom) || !chrom_lengths.contains(chrom)) {
      // skip empty features or chromosomes not in the VCF header
      continue;
    }
    // this chromosome is also in the POPAF VCF header
    const s32 chrom_index = chrom_indexes.at(chrom);
    const u64 chrom_length = chrom_lengths.at(chrom);
    const u64 num_feat = features.size();
    const u64 num_feat_per_thread = num_feat / threads + (num_feat % threads != 0 ? 1 : 0);
    const u64 block_size = std::max<u64>(kMinBlockSize, num_feat_per_thread);
    for (u64 i = 0; i < num_feat; i += block_size) {
      // each task will be responsible for a specific range of the vector of `VcfFeature`
      auto feat_span = std::span(features.data() + i, std::min(num_feat - i, block_size));
      if (!feat_span.empty()) {
        feat_blocks.emplace_back(feat_span, chrom_index, chrom_length);
      }
    }
  }

  std::atomic_uint32_t num_chrom_processed{0};

  // Setups up a task that iterates through a given block of the popaf vcf and updates the features for the matching
  // chromosomal positions
  auto task = [&vcf_path, &num_chrom_processed](ChromFeatureBlock& block) {
    try {
      auto reader = io::VcfReader(vcf_path);
      UpdatePopafAtBlock(block, reader);
      num_chrom_processed.fetch_add(1);
    } catch (const std::exception& e) {
      const auto& first = block.features.begin();
      throw error::Error("Error extracting POPAF values in feature block '{}:{}-{}': {}",
                         first->chrom,
                         first->pos,
                         block.features.back().pos,
                         e.what());
    }
  };

  tf::Executor executor(threads);
  tf::Taskflow flow;
  for (auto& block : feat_blocks) {
    flow.emplace([&block, &task] { task(block); });
  }
  executor.run(flow).get();

  if (feat_blocks.size() != num_chrom_processed) {
    throw error::Error("Number of feature blocks {} and tasks executed {} are not the same!",
                       feat_blocks.size(),
                       std::to_string(num_chrom_processed));
  }
}

VcfHeaderInfo GetVcfHeaderInfo(const io::VcfHeaderPtr& header) {
  // This method reads the VCF header and extracts relevant information for feature extraction:
  // - reference contig indexes and lengths
  // - tumor/normal sample indexes
  // - presence of INFO/FORMAT fields of interest
  // Since these are invariant to individual VCF records, these are stored to avoid repeated lookups in the VCF header.
  VcfHeaderInfo info;
  info.contig_lengths = header->GetContigLengths();
  info.contig_indexes = header->GetContigIndexes();

  // extract sample indexes for `tumor_sample` and `normal_sample` if they are found in the header
  const auto tn_sample_idx = GetTumorNormalSampleIndexes(header);
  if (tn_sample_idx.has_value()) {
    info.tumor_index = tn_sample_idx->tumor_sample_idx;
    info.normal_index = tn_sample_idx->normal_sample_idx;
    info.has_tumor_normal = true;
  } else {
    info.tumor_index = -1;
    info.normal_index = -1;
    info.has_tumor_normal = false;
  }

  // check header for field names
  info.has_popaf = header->HasInfoField(kFieldPopaf);
  info.has_nalod = header->HasInfoField(kFieldNalod);
  info.has_nlod = header->HasInfoField(kFieldNlod);
  info.has_tlod = header->HasInfoField(kFieldTlod);
  info.has_mpos = header->HasInfoField(kFieldMpos);
  info.has_mmq = header->HasInfoField(kFieldMmq);
  info.has_mbq = header->HasInfoField(kFieldMbq);
  info.has_ad = header->HasInfoField(kFieldAd);
  info.has_af = header->HasInfoField(kFieldAf);
  info.has_dp = header->HasInfoField(kFieldDp);
  info.has_gq = header->HasInfoField(kFieldGq);
  info.has_gt = header->HasInfoField(kFieldGt);
  info.has_hapcomp = header->HasInfoField(kFieldHapcomp);
  info.has_hapdom = header->HasInfoField(kFieldHapdom);
  info.has_str = header->HasInfoField(kFieldStr);
  info.has_ru = header->HasInfoField(kFieldRu);
  info.has_rpa = header->HasInfoField(kFieldRpa);
  return info;
}

/**
 * @brief Update the variant density values for all features.
 * @param chrom_to_features Map of chromosome to vector of VcfFeature
 */
static void UpdateVariantDensity(StrMap<vec<VcfFeature>>& chrom_to_features) {
  // Variant density is defined as the number of variants within a window upstream and downstream of the variant's
  // start position. If there is only one variant in the window, the density is 1.
  for (auto& [chrom, features] : chrom_to_features) {
    const u64 num_feat = features.size();
    for (u64 i = 0; i < num_feat; ++i) {
      auto& feat = features[i];
      const auto max_pos = feat.pos + kVariantDensityUpstreamWindow;
      for (u64 j = i + 1; j < num_feat && features[j].pos <= max_pos; ++j) {
        ++feat.variant_density;
        ++features[j].variant_density;
      }
    }
  }
}

/**
 * @brief Append one variant's VCF features as one line to a TSV file.
 * @param feat VCF features for one variant
 * @param feature_cols Vector of feature column enums
 */
static vec<std::string> FeatureToString(const VcfFeature& feat, const vec<FeatureColumn>& feature_cols) {
  // convert coordinates back to 1-based
  const auto vcf_feature_serialized = VariantInfoSerializer::SerializeVcfFeature(feat);
  vec<std::string> feature_strings;
  feature_strings.reserve(feature_cols.size());
  for (const auto& [enum_val, prefix] : feature_cols) {
    feature_strings.push_back(vcf_feature_serialized.at(enum_val));
  }
  return feature_strings;
}

/**
 * @brief Write all variants' VCF features to a TSV file.
 * @param writer LockedTsvWriter instance for thread-safe writing
 * @param feature_names Vector of feature names
 * @param feature_cols Vector of feature columns
 * @param chrom_to_features Map of chromosome to vector of VcfFeature
 * @param threads Number of threads
 * @param cmd_info Optional command line information for adding comments to the output file
 */
static void WriteAllFeatures(LockedTsvWriter& writer,
                             const vec<std::string>& feature_names,
                             const vec<FeatureColumn>& feature_cols,
                             StrMap<vec<VcfFeature>>& chrom_to_features,
                             const u32 threads,
                             const std::optional<CommandLineInfo>& cmd_info) {
  // write version and command line as comments
  if (cmd_info.has_value()) {
    io::Comments cmt;
    io::AddVersionAndCommandLineComment(cmt, cmd_info.value().version, cmd_info.value().command_line);
    writer.AppendComments(cmt);
  }
  // write header row
  writer.AppendRow(feature_names);

  if (threads <= 1) {
    // serialize all features to an output TSV file in a single thread
    for (const auto& [chrom, features] : chrom_to_features) {
      for (const auto& f : features) {
        writer.AppendRow(FeatureToString(f, feature_cols));
      }
    }
    return;
  }
  // use multiple threads to write features in parallel
  // determine the block size for each thread
  u64 total_num_features = 0;
  for (const auto& [chrom, features] : chrom_to_features) {
    total_num_features += features.size();
  }
  const u64 num_features_per_thread = total_num_features / threads + (total_num_features % threads != 0 ? 1 : 0);
  const u64 block_size = std::max<u64>(kMinBlockSize, num_features_per_thread);

  // assign feature blocks
  vec<std::span<VcfFeature>> blocks;
  for (auto& [chrom, features] : chrom_to_features) {
    const u64 num_feat = features.size();
    for (u64 i = 0; i < num_feat; i += block_size) {
      auto feat_span = std::span(features.data() + i, std::min(num_feat - i, block_size));
      if (!feat_span.empty()) {
        blocks.emplace_back(feat_span);
      }
    }
  }

  // set up one task per feature block
  auto task = [&feature_cols, &writer](const std::span<VcfFeature>& block) {
    try {
      vec<vec<std::string>> feature_strings;
      for (const auto& feat : block) {
        feature_strings.emplace_back(FeatureToString(feat, feature_cols));
      }
      writer.AppendRows(feature_strings);
    } catch (const std::exception& e) {
      throw error::Error("Error in serializing feature block: {}", e.what());
    }
  };

  tf::Executor executor(threads);
  tf::Taskflow flow;
  for (const auto& block : blocks) {
    flow.emplace([&task, &block] { task(block); });
  }
  executor.run(flow).get();
}

/**
 * @brief Update the repeat counts in a VCF feature struct.
 * @param feat VCF feature struct
 */
static void UpdateRepeatCounts(VcfFeature& feat) {
  if (!feat.str || feat.ru.empty()) {
    // no repeat unit or STR information available
    return;
  }

  // repeat unit lengths for tandem repeats
  static constexpr size_t kHomoPolymerRepeatUnitLength{1};
  static constexpr size_t kDiRepeatUnitLength{2};
  static constexpr size_t kTriRepeatUnitLength{3};
  static constexpr size_t kQuadRepeatUnitLength{4};

  // RPA is the number of times tandem repeat unit is repeated.
  // There are values for both REF and ALT alleles, but we are only interested in REF.
  switch (feat.ru.length()) {
    case kHomoPolymerRepeatUnitLength:
      feat.homopolymer = feat.rpa_ref;
      break;
    case kDiRepeatUnitLength:
      feat.direpeat = feat.rpa_ref;
      break;
    case kTriRepeatUnitLength:
      feat.trirepeat = feat.rpa_ref;
      break;
    case kQuadRepeatUnitLength:
      feat.quadrepeat = feat.rpa_ref;
      break;
    default:
      break;
  }
}

/**
 * @brief Update a VcfFeature struct based on Mutect2 fields in a VCF record.
 * @pre The VCF record must have been called by Mutect2 and contain the relevant INFO/FORMAT fields.
 * @pre Allele index must be provided to indicate which ALT allele's features to extract for multi-allelic records.
 * @param feat VcfFeature struct to be updated
 * @param record VcfRecordPtr containing the variant information
 * @param header_info Extracted information from VCF header
 * @param use_popaf_field Flag whether to use POPAF VCF field
 * @param allele_index Index of ALT allele in VcfRecord
 */
static void UpdateMutect2Features(VcfFeature& feat,
                                  const io::VcfRecordPtr& record,
                                  const VcfHeaderInfo& header_info,
                                  const bool use_popaf_field,
                                  const s32 allele_index) {
  if (allele_index < 1) {
    throw error::Error("Error extracting VCF features at {}:{}. ALT allele index must be greater than 0 ({})",
                       record->Chromosome(),
                       record->Position(),
                       allele_index);
  }
  const auto allele_idx_u = static_cast<u32>(allele_index);
  const auto value_idx_u = allele_idx_u - 1;
  feat.popaf = GetInfo<f32, f32>(use_popaf_field && header_info.has_popaf, record, value_idx_u, kFieldPopaf);
  feat.nalod = GetInfo<f32, f32>(header_info.has_nalod, record, value_idx_u, kFieldNalod);
  feat.nlod = GetInfo<f32, f32>(header_info.has_nlod, record, value_idx_u, kFieldNlod);
  feat.tlod = GetInfo<f32, f32>(header_info.has_tlod, record, value_idx_u, kFieldTlod);
  feat.mpos = GetInfo<s32, u64>(header_info.has_mpos, record, value_idx_u, kFieldMpos);
  std::tie(feat.mmq_ref, feat.mmq_alt) = GetInfoRefAlt<s32, u8>(header_info.has_mmq, record, allele_idx_u, kFieldMmq);
  std::tie(feat.mbq_ref, feat.mbq_alt) = GetInfoRefAlt<s32, u8>(header_info.has_mbq, record, allele_idx_u, kFieldMbq);
  if (header_info.has_tumor_normal) {
    if (header_info.has_ad) {
      const auto& values = record->GetFormatFieldNoCheck<s32>(kFieldAd);
      if (values.size() == 4) {
        // There are AD values for both REF and ALT
        // REF values are at indexes 0 and 2; ALT values are at indexes 1 and 3
        feat.tumor_alt_ad = values.at(2 * header_info.tumor_index + 1);    // +1 for ALT
        feat.normal_alt_ad = values.at(2 * header_info.normal_index + 1);  // +1 for ALT
        feat.alt_ad = feat.tumor_alt_ad + feat.normal_alt_ad;
        feat.ref_ad = values.at(0) + values.at(2);
      }
    }

    if (header_info.has_af) {
      // Note that this is FORMAT field AF
      const auto& values = record->GetFormatFieldNoCheck<f32>(kFieldAf);
      if (values.size() == 2) {
        feat.tumor_af = values.at(header_info.tumor_index);
        feat.normal_af = values.at(header_info.normal_index);
        if (feat.normal_af > 0) {
          feat.tn_af_ratio = feat.tumor_af / feat.normal_af;
        }
      }
    }

    if (header_info.has_dp) {
      const auto& values = record->GetFormatFieldNoCheck<s32>(kFieldDp);
      if (values.size() == 2) {
        feat.tumor_dp = values.at(header_info.tumor_index);
        feat.normal_dp = values.at(header_info.normal_index);
      }
    }
  }
}

/**
 * @brief Extract features from VcfRecord.
 * @param record VcfRecord
 * @param header_info Extracted information from VCF header
 * @param use_popaf_field Flag whether to use POPAF field from the VCF record
 * @param allele_index Index of allele in VcfRecord, must be >= 1
 * @return VcfFeature struct
 */
static std::optional<VcfFeature> RecordToFeature(const io::VcfRecordPtr& record,
                                                 const VcfHeaderInfo& header_info,
                                                 const bool use_popaf_field,
                                                 const s32 allele_index) {
  // Overview:
  // 1. Trim the variant representation into its shortest form as needed.
  // 2. Create a VcfFeature struct and populate it with the variant ID information.
  // 3. Extract INFO/FORMAT fields from the record and populate the VcfFeature struct.

  if (CheckAndWarnNegativePosition(record)) {
    return std::nullopt;
  }
  if (allele_index < 1) {
    throw error::Error("Error extracting VCF features at {}:{}. ALT allele index must be greater than 0 ({})",
                       record->Chromosome(),
                       record->Position(),
                       allele_index);
  }
  const std::string alt = record->Allele(allele_index);
  if (alt == "*") {
    // do not extract features for wildcard allele
    return std::nullopt;
  }

  // typically, variant representation is only trimmed for multi-allelic indel records
  const auto& [ref_trimmed, alt_trimmed] = TrimVariant(record->Allele(0), alt);
  VcfFeature feat{.chrom = record->Chromosome(),
                  .pos = static_cast<u64>(record->Position()),
                  .ref = ref_trimmed,
                  .alt = alt_trimmed,
                  .qual = record->GetQuality(0)};

  // `value_index` is intended for fields that do not account for the reference allele
  // therefore, it is equal to the allele index minus one
  const auto allele_idx_u = static_cast<u32>(allele_index);
  const auto value_index = allele_idx_u - 1;

  // extract INFO and FORMAT field values from the record
  feat.hapcomp = GetInfo<s32, u32>(header_info.has_hapcomp, record, value_index, kFieldHapcomp);
  feat.hapdom = GetInfo<f32, f32>(header_info.has_hapdom, record, value_index, kFieldHapdom);

  if (header_info.has_str) {
    feat.str = record->HasInfoFieldNoCheck(kFieldStr);
    if (feat.str) {
      // Although STR, RU, and RPA may all appear in the header, RU and RPA are only present if STR is in the record
      if (header_info.has_ru) {
        feat.ru = record->GetInfoFieldStringNoCheck(kFieldRu);
      }
      std::tie(feat.rpa_ref, feat.rpa_alt) =
          GetInfoRefAlt<s32, u32>(header_info.has_rpa, record, allele_index, kFieldRpa);
      UpdateRepeatCounts(feat);
    }
  }

  if (!header_info.has_tumor_normal) {
    feat.gq = GetFormat<s32, u32>(header_info.has_gq, record, 0, kFieldGq);
    // Note that this is INFO field AF
    feat.normal_af = GetInfo<f32, f32>(header_info.has_af, record, 0, kFieldAf);
    feat.normal_dp = GetFormat<s32, u32>(header_info.has_dp, record, 0, kFieldDp);

    if (header_info.has_gt) {
      const std::string& value = record->GetGTField();
      if (!value.empty()) {
        feat.genotype = StringToGenotype(value);
      }
    }

    // Maximum number of alleles in a germline VCF record, e.g. 1 REF + 2 ALT
    static constexpr u32 kGermlineRecordMaxAlleles{3};
    if (std::cmp_greater_equal(record->NumAlleles(), kGermlineRecordMaxAlleles)) {
      // multi-allelic variant
      std::tie(feat.ref_ad, feat.alt_ad, feat.alt_ad2) =
          GetFormatThreeAlleles<s32, u32>(header_info.has_ad, record, allele_idx_u, kFieldAd);
      if (feat.normal_dp > 0) {
        const auto dp = static_cast<f64>(feat.normal_dp);
        feat.ref_ad_af = feat.ref_ad / dp;
        feat.alt_ad_af = feat.alt_ad / dp;
        feat.alt_ad2_af = feat.alt_ad2 / dp;
      }
    } else {
      std::tie(feat.ref_ad, feat.alt_ad) =
          GetFormatRefAlt<s32, u32>(header_info.has_ad, record, allele_idx_u, kFieldAd);
      if (feat.normal_dp > 0) {
        const auto dp = static_cast<f64>(feat.normal_dp);
        feat.ref_ad_af = feat.ref_ad / dp;
        feat.alt_ad_af = feat.alt_ad / dp;
      }
    }
  }

  // Update features based on Mutect2 specific fields
  UpdateMutect2Features(feat, record, header_info, use_popaf_field, allele_index);

  return feat;
}

RepeatCounts GetRepeatCounts(const std::string& seq) {
  // Count the number of repetitions for each repeat type
  u32 homopolymer = GetHomopolymerLength(seq);
  u32 direpeat = CountRepeats(seq, 2);
  u32 trirepeat = CountRepeats(seq, 3);
  u32 quadrepeat = CountRepeats(seq, 4);
  // Calculate the length covered by each repeat type
  const u32 len2 = 2 * direpeat;
  const u32 len3 = 3 * trirepeat;
  const u32 len4 = 4 * quadrepeat;
  // Handle ambiguous repeat unit representations
  // Choose the repeat unit that covers the most sequence length
  // If there is a tie, chose the shorter repeat unit
  if (homopolymer > 0 && homopolymer >= len2 && homopolymer >= len3 && homopolymer >= len4) {
    direpeat = 0;
    trirepeat = 0;
    quadrepeat = 0;
  } else if (direpeat > 0 && len2 > homopolymer && len2 >= len3 && len2 >= len4) {
    homopolymer = 0;
    trirepeat = 0;
    quadrepeat = 0;
  } else if (trirepeat > 0 && len3 > homopolymer && len3 > len2 && len3 >= len4) {
    homopolymer = 0;
    direpeat = 0;
    quadrepeat = 0;
  } else if (quadrepeat > 0 && len4 > homopolymer && len4 > len2 && len4 > len3) {
    homopolymer = 0;
    direpeat = 0;
    trirepeat = 0;
  }
  return RepeatCounts{homopolymer, direpeat, trirepeat, quadrepeat};
}

/**
 * @brief Information about a genomic region for parallel processing.
 */
struct RegionInfo {
  // chromosome name
  const std::string chrom;
  // chromosome index in the VCF index
  const s32 chrom_index;
  // 0-based, inclusive start position
  const u64 start;
  // 0-based, exclusive end position
  const u64 end;
  // flag whether the POPAF VCF is provided
  const bool has_popaf_vcf;
  // reference sequence for the chromosome
  const std::string& chrom_seq;
  // flag whether to extract features for germline-tagging workflow
  const bool is_germline_tagging;
};

/**
 * @brief Features extracted from the reference sequence.
 */
struct RefSeqFeatures {
  std::string pre_2bp_context{};
  std::string post_2bp_context{};
  std::string post_30bp_context{};
  u32 uniq_3mers{0};
  u32 uniq_4mers{0};
  u32 uniq_5mers{0};
  u32 uniq_6mers{0};
  u32 homopolymer{0};
  u32 direpeat{0};
  u32 trirepeat{0};
  u32 quadrepeat{0};
};

/**
 * @brief Extract sequence context features from the reference sequence.
 * @param ref_seq Reference sequence
 * @param pos Position of the variant in the reference sequence (0-based)
 * @param count_repeats Flag whether to count tandem repeats in the sequence context
 * @return RefSeqFeatures struct containing extracted sequence context features
 */
static RefSeqFeatures GetRefSeqFeatures(const std::string& ref_seq, const u64 pos, const bool count_repeats) {
  static constexpr u64 kShortContext{2};
  static constexpr u64 kLongContext{30};
  // start position of the sequence context after the variant position
  const auto post_context_start = pos + 1;
  const auto seq_len = ref_seq.size();

  RefSeqFeatures feat;
  feat.pre_2bp_context = pos >= kShortContext ? ref_seq.substr(pos - kShortContext, kShortContext) : "";
  feat.post_2bp_context =
      post_context_start + kShortContext <= seq_len ? ref_seq.substr(post_context_start, kShortContext) : "";
  if (post_context_start + kLongContext <= seq_len) {
    // has enough sequence context after the variant position
    feat.post_30bp_context = ref_seq.substr(post_context_start, kLongContext);
    feat.uniq_3mers = CountUniqueKmers(feat.post_30bp_context, 3);
    feat.uniq_4mers = CountUniqueKmers(feat.post_30bp_context, 4);
    feat.uniq_5mers = CountUniqueKmers(feat.post_30bp_context, 5);
    feat.uniq_6mers = CountUniqueKmers(feat.post_30bp_context, 6);
  }
  if (count_repeats) {
    // Counting starts at position + 1 to match GATK's behavior for `RU` and `RPA`.
    const auto seq = ref_seq.substr(pos + 1, kTandemRepeatSearchDistance);
    const auto& [homopolymer, direpeat, trirepeat, quadrepeat] = GetRepeatCounts(seq);
    feat.homopolymer = homopolymer;
    feat.direpeat = direpeat;
    feat.trirepeat = trirepeat;
    feat.quadrepeat = quadrepeat;
  }
  return feat;
}

/**
 * @brief Update a VcfFeature struct with extracted reference sequence features.
 * @param feat VcfFeature struct to be updated
 * @param ref_seq_feat RefSeqFeatures struct containing the extracted features
 */
static void UpdateFeatureWithRefSeqFeatures(VcfFeature& feat,
                                            const RefSeqFeatures& ref_seq_feat,
                                            const bool count_repeats) {
  // Update the VcfFeature struct with the extracted reference sequence features
  feat.pre_2bp_context = ref_seq_feat.pre_2bp_context;
  feat.post_2bp_context = ref_seq_feat.post_2bp_context;
  feat.post_30bp_context = ref_seq_feat.post_30bp_context;
  feat.uniq_3mers = ref_seq_feat.uniq_3mers;
  feat.uniq_4mers = ref_seq_feat.uniq_4mers;
  feat.uniq_5mers = ref_seq_feat.uniq_5mers;
  feat.uniq_6mers = ref_seq_feat.uniq_6mers;
  if (count_repeats) {
    feat.homopolymer = ref_seq_feat.homopolymer;
    feat.direpeat = ref_seq_feat.direpeat;
    feat.trirepeat = ref_seq_feat.trirepeat;
    feat.quadrepeat = ref_seq_feat.quadrepeat;
  }
}

/**
 * @brief Check if a variant position overlaps with an interval.
 * @param intervals Vector of intervals
 * @param start Start position of the variant (0-based, inclusive)
 * @param end End position of the variant (0-based, exclusive)
 * @param itr Iterator to `intervals`
 * @param padding Padding to consider for the overlap check
 * @return True if the variant overlaps with a repeat region, false otherwise
 * @note `itr` will be advanced to the variant position.
 */
static bool FindInterval(
    const vec<Interval>& intervals, const u64 start, const u64 end, vec<Interval>::iterator& itr, const u64 padding) {
  if (intervals.empty() || itr == intervals.end()) {
    return false;
  }
  // keep advancing the iterator until the current interval's end is either at or
  // greater than this variant's start position
  while (itr != intervals.end() && itr->end < start) {
    ++itr;
  }
  if (itr != intervals.end()) {
    // determine whether a variant is within an interval
    return IntervalOverlap(*itr, start >= padding ? start - padding : start, end + padding);
  }
  return false;
}

bool HasPassFilter(const io::VcfRecordPtr& record) {
  const auto& filters = record->GetFilters();
  return std::ranges::find(filters, kFilteringPassId) != filters.end();
}

/**
 * @brief Helper function to compute features for target regions in the same chromosome.
 * @param region VCF region information
 * @param vcf_reader VCF reader
 * @param target_regions Map of chromosome to target intervals
 * @param interest_regions Map of chromosome to repeat intervals
 * @param header_info Extracted information from VCF header
 * @return Vector of VcfFeature structs
 */
static vec<VcfFeature> ComputeVcfFeaturesForRegion(const RegionInfo& region,
                                                   io::VcfReader& vcf_reader,
                                                   const ChromIntervalsMap& target_regions,
                                                   const ChromIntervalsMap& interest_regions,
                                                   const VcfHeaderInfo& header_info) {
  // This function does not update variant density for each feature struct; the update is done outside of this function
  // once all features in the same chromosome has been extracted.
  // This function does not read the POPAF VCF to update POPAF values.

  // Overview:
  // 1. Set the region in the VCF reader
  // 2. Iterate over the VCF records in the region
  //   2.1. Check if the record overlaps with target regions; advance the iterator for target regions as necessary
  //   2.2. Check if the record overlaps with repeat regions; advance the iterator for repeat regions as necessary
  //   2.3. Extract sequence context features from the reference genome
  //   2.4. Extract repeat counts from the VCF record; otherwise, count repeats using the reference genome sequence
  //   2.5. Extract INFO/FORMAT field features from the VCF record
  // 3. Return the features extracted for the region

  vec<VcfFeature> features;
  if (!vcf_reader.SetRegion(region.chrom_index, region.start, region.end)) {
    return features;
  }

  const bool has_target_regions = !target_regions.empty() && target_regions.contains(region.chrom);
  vec<Interval> target_intervals;
  vec<Interval>::iterator target_intervals_itr;
  if (has_target_regions) {
    target_intervals = target_regions.at(region.chrom);
    target_intervals_itr = target_intervals.begin();
  }

  vec<Interval> interest_intervals;
  vec<Interval>::iterator interest_intervals_itr;
  if (!interest_regions.empty() && interest_regions.contains(region.chrom)) {
    interest_intervals = interest_regions.at(region.chrom);
    interest_intervals_itr = interest_intervals.begin();
  }

  // padding of 0 to ensure that the variant is within the target region
  static constexpr u64 kTargetPadding{0};
  // padding of 1 for variants adjacent to the interest region
  static constexpr u64 kAtInterestPadding{1};

  const bool has_repeat_fields = header_info.has_str && header_info.has_ru && header_info.has_rpa;
  while (const auto& record = vcf_reader.GetNextRegionRecord(BCF_UN_ALL)) {
    if (CheckAndWarnNegativePosition(record)) {
      continue;
    }
    if (region.is_germline_tagging && HasPassFilter(record)) {
      continue;
    }
    const auto pos = static_cast<u64>(record->Position());
    const u64 pos_end = pos + record->Allele(0).length();
    if (has_target_regions && !FindInterval(target_intervals, pos, pos_end, target_intervals_itr, kTargetPadding)) {
      // target region not found, skip this record
      continue;
    }
    const bool at_interest = FindInterval(interest_intervals, pos, pos_end, interest_intervals_itr, kAtInterestPadding);
    const bool count_repeats = !has_repeat_fields && pos + 1 + kTandemRepeatSearchDistance <= region.chrom_seq.size();
    const RefSeqFeatures ref_seq_feat = GetRefSeqFeatures(region.chrom_seq, pos, count_repeats);

    // each record can have more than one ALT
    const s32 num_alleles = record->NumAlleles();
    if (num_alleles < 1) {
      // must have at least 2 alleles: REF and ALT
      continue;
    }
    for (s32 i = 1; i < num_alleles; ++i) {
      auto ret = RecordToFeature(record, header_info, !region.has_popaf_vcf, i);
      if (ret.has_value()) {
        VcfFeature feat = ret.value();
        // variant density is initialized to 1 for itself
        feat.variant_density = 1;
        feat.at_interest = at_interest;
        UpdateFeatureWithRefSeqFeatures(feat, ref_seq_feat, count_repeats);
        features.emplace_back(feat);
      }
    }
  }

  return features;
}

/**
 * @brief Sort a vector of VcfFeature structs by position in ascending order.
 * @param features Vector of VcfFeature structs for the same chromosome
 */
static void SortByPosition(vec<VcfFeature>& features) {
  std::ranges::sort(features, [](const VcfFeature& a, const VcfFeature& b) { return a.pos < b.pos; });
}

VarIdToVcfFeatures ExtractFeaturesForRegion(io::VcfReader& vcf_reader,
                                            const fs::path& genome_path,
                                            std::optional<io::VcfReader>& popaf_reader,
                                            const ChromIntervalsMap& target_regions,
                                            const ChromIntervalsMap& interest_regions,
                                            const std::string& chrom,
                                            const bool is_germline_tagging) {
  // Overview:
  // 1. Extract relevant information from VCF header
  // 2. Extract features for target regions in the same chromosome
  // 3. Update POPAF values for extracted features if a POPAF VCF is provided
  // 4. Update variant density for each extracted feature struct
  // 5. Convert features vector to PositionToVcfFeaturesMap

  const ChromToIndexes& chrom_indexes = vcf_reader.GetContigIndexes();
  if (!chrom_indexes.contains(chrom)) {
    throw error::Error("Error finding VCF records for region '{}:{}-{}'",
                       chrom,
                       target_regions.at(chrom).begin()->start,
                       (target_regions.at(chrom).end() - 1)->end);
  }

  io::FastaReader fasta_reader(genome_path);
  const s32 cid = chrom_indexes.at(chrom);
  const auto& chrom_seq = fasta_reader.GetSequence(chrom);
  const auto vcf_region = RegionInfo{.chrom = chrom,
                                     .chrom_index = cid,
                                     .start = target_regions.at(chrom).begin()->start,
                                     .end = (target_regions.at(chrom).end() - 1)->end,
                                     .has_popaf_vcf = popaf_reader.has_value(),
                                     .chrom_seq = chrom_seq,
                                     .is_germline_tagging = is_germline_tagging};

  const VcfHeaderInfo& header_info = GetVcfHeaderInfo(vcf_reader.GetHeader());
  auto features = ComputeVcfFeaturesForRegion(vcf_region, vcf_reader, target_regions, interest_regions, header_info);
  if (features.empty()) {
    return {};
  }

  // sort the features by position so POPAF and variant density can be updated properly
  SortByPosition(features);

  if (popaf_reader.has_value()) {
    const auto& [popaf_chrom_indexes, popaf_chrom_lengths] = GetContigInfo(popaf_reader.value());
    if (popaf_chrom_indexes.contains(chrom) && popaf_chrom_lengths.contains(chrom)) {
      // only update POPAF if this chromosome is also in the POPAF VCF header
      ChromFeatureBlock feat_block{features, popaf_chrom_indexes.at(chrom), popaf_chrom_lengths.at(chrom)};
      UpdatePopafAtBlock(feat_block, popaf_reader.value());
    }
  }

  StrMap<vec<VcfFeature>> chrom_features;
  chrom_features[chrom] = features;
  UpdateVariantDensity(chrom_features);

  // Return a map of VariantId to VcfFeature for easy lookup
  VarIdToVcfFeatures features_map;
  for (const auto& feat : chrom_features.at(chrom)) {
    const VariantId id(feat.chrom, feat.pos, feat.ref, feat.alt);
    features_map[id] = feat;
  }

  return features_map;
}

/**
 * @brief Warn if a feature column is specified but the corresponding field is not found in the VCF header.
 * @param feat_col_set Set of all feature columns
 * @param feat_col Feature column
 * @param feat_name Feature name
 * @param field_type VCF field type, e.g. INFO or FORMAT
 * @param field_name VCF field name
 * @param found_in_header Flag whether VCF field is found in VCF header
 * @return Flag whether a warning was issued
 */
static u16 WarnIfMissing(const std::unordered_set<UnifiedFeatureCols>& feat_col_set,
                         const UnifiedFeatureCols feat_col,
                         const std::string& feat_name,
                         const std::string& field_type,
                         const std::string& field_name,
                         const bool found_in_header) {
  if (!found_in_header && feat_col_set.contains(feat_col)) {
    WarnAsErrorIfSet("VCF feature '{}' is specified, but {} field '{}' is not found in the VCF header!",
                     feat_name,
                     field_type,
                     field_name);
    return 1;
  }
  return 0;
}

/**
 * @brief Check feature column names against VCF INFO fields.
 * @param cols Set of feature columns
 * @param info Extracted information from VCF header
 * @return Number of warnings reported
 */
static u16 CheckVcfFeatureInfoFields(const std::unordered_set<UnifiedFeatureCols>& cols, const VcfHeaderInfo& info) {
  using enum UnifiedFeatureCols;
  static const std::string kType{"INFO"};
  return WarnIfMissing(cols, kVcfNormalAf, kNameNormalAF, kType, kFieldAf, info.has_af) +
         WarnIfMissing(cols, kVcfHapcomp, kNameHapcomp, kType, kFieldHapcomp, info.has_hapcomp) +
         WarnIfMissing(cols, kVcfHapdom, kNameHapdom, kType, kFieldHapdom, info.has_hapdom) +
         WarnIfMissing(cols, kVcfMbqRef, kNameMbqRef, kType, kFieldMbq, info.has_mbq) +
         WarnIfMissing(cols, kVcfMbqAlt, kNameMbqAlt, kType, kFieldMbq, info.has_mbq) +
         WarnIfMissing(cols, kVcfMmqRef, kNameMmqRef, kType, kFieldMmq, info.has_mmq) +
         WarnIfMissing(cols, kVcfMmqAlt, kNameMmqAlt, kType, kFieldMmq, info.has_mmq) +
         WarnIfMissing(cols, kVcfMpos, kNameMpos, kType, kFieldMpos, info.has_mpos) +
         WarnIfMissing(cols, kVcfNalod, kNameNalod, kType, kFieldNalod, info.has_nalod) +
         WarnIfMissing(cols, kVcfNlod, kNameNlod, kType, kFieldNlod, info.has_nlod) +
         WarnIfMissing(cols, kVcfRu, kNameRu, kType, kFieldRu, info.has_ru) +
         WarnIfMissing(cols, kVcfRpaRef, kNameRpaRef, kType, kFieldRpa, info.has_rpa) +
         WarnIfMissing(cols, kVcfRpaAlt, kNameRpaAlt, kType, kFieldRpa, info.has_rpa) +
         WarnIfMissing(cols, kVcfStr, kNameStr, kType, kFieldStr, info.has_str) +
         WarnIfMissing(cols, kVcfTlod, kNameTlod, kType, kFieldTlod, info.has_tlod);
}

/**
 * @brief Check feature column names against VCF FORMAT fields.
 * @param cols Set of feature columns
 * @param info Extracted information from VCF header
 * @return Number of warnings reported
 */
static u16 CheckVcfFeatureFormatFields(const std::unordered_set<UnifiedFeatureCols>& cols, const VcfHeaderInfo& info) {
  using enum UnifiedFeatureCols;
  static const std::string kType{"FORMAT"};
  return WarnIfMissing(cols, kVcfAltAd, kNameAltAD, kType, kFieldAd, info.has_ad) +
         WarnIfMissing(cols, kVcfRefAd, kNameRefAD, kType, kFieldAd, info.has_ad) +
         WarnIfMissing(cols, kVcfNormalDp, kNameNormalDP, kType, kFieldDp, info.has_dp) +
         WarnIfMissing(cols, kVcfGq, kNameGq, kType, kFieldGq, info.has_gq) +
         WarnIfMissing(cols, kVcfGenotype, kNameGenotype, kType, kFieldGt, info.has_gt);
}

/**
 * @brief Check feature column names against available resources and VCF fields; produce warnings for the absence
 * of resources or VCF fields.
 * @param feat_cols Vector of feature column
 * @param interest_regions Map of chromosome to vector of repeat region intervals
 * @param header_info VCF header information
 * @param has_popaf_path Flag whether POPAF VCF file is available
 * @return Number of warnings reported
 */
static u16 CheckVcfFeatureResources(const vec<FeatureColumn>& feat_cols,
                                    const ChromIntervalsMap& interest_regions,
                                    const VcfHeaderInfo& header_info,
                                    const bool has_popaf_path) {
  using enum UnifiedFeatureCols;
  std::unordered_set<UnifiedFeatureCols> feat_cols_set;
  for (const auto& [enum_val, prefix] : feat_cols) {
    // insert into set for easy lookup
    // ignore prefix because it is not relevant for resource checking
    feat_cols_set.insert(enum_val);
  }

  u16 num_warnings{0};
  if (interest_regions.empty() && feat_cols_set.contains(kVcfAtInterest)) {
    WarnAsErrorIfSet("VCF feature '{}' is specified, but there are no Interest regions!", kNameAtInterest);
    ++num_warnings;
  }
  if (!interest_regions.empty() && !feat_cols_set.contains(kVcfAtInterest)) {
    WarnAsErrorIfSet("Interest regions are specified, but VCF feature '{}' is not specified!", kNameAtInterest);
    ++num_warnings;
  }
  if (!has_popaf_path && !header_info.has_popaf && feat_cols_set.contains(kVcfPopAf)) {
    WarnAsErrorIfSet(
        "VCF feature '{}' is specified, but POPAF VCF file is not provided and the INFO field '{}' is not "
        "found in the VCF header!",
        kNamePopAF,
        kFieldPopaf);
    ++num_warnings;
  }
  num_warnings += CheckVcfFeatureInfoFields(feat_cols_set, header_info);
  num_warnings += CheckVcfFeatureFormatFields(feat_cols_set, header_info);
  return num_warnings;
}

u16 CheckVcfFeatureResources(const vec<FeatureColumn>& feat_cols,
                             const ChromIntervalsMap& interest_regions,
                             const io::VcfHeaderPtr& header,
                             const bool has_popaf_vcf) {
  const auto header_info = GetVcfHeaderInfo(header);
  return CheckVcfFeatureResources(feat_cols, interest_regions, header_info, has_popaf_vcf);
}

/**
 * @brief Extract features from VCF file using one or more threads.
 * @param param Parameters for VCF feature extraction
 * @param feat_cols Vector of feature column
 * @return Map of chromosome to vector of VcfFeature structs sorted by position
 * @note Input VCF file must be indexed.
 */
static StrMap<vec<VcfFeature>> ExtractFeatures(const ComputeVcfFeaturesParam& param,
                                               const vec<FeatureColumn>& feat_cols) {
  // Overview:
  // 1. Extract relevant information from VCF header
  // 2. check if the VCF file has the required INFO and FORMAT fields and whether resource files are available
  // 3. load reference sequences for each chromosome in the target regions
  // 4. create a task for each chromosome in the target regions
  // 5. run the tasks in parallel using a Taskflow thread pool
  // 6. sort the extracted features for each chromosome
  // 7. return the map of chromosome to vector of VcfFeature structs

  Logging::Info("Begin VCF feature extraction");
  const bool has_popaf_vcf{param.pop_af_vcf.has_value()};
  const SVCConfig model_config = param.config;

  const auto& target_regions = param.target_regions.has_value() ? param.target_regions.value() : ChromIntervalsMap{};
  const auto& interest_regions =
      param.interest_regions.has_value() ? param.interest_regions.value() : ChromIntervalsMap{};

  const auto& vcf_path = param.vcf_file;
  VcfHeaderInfo header_info;
  ChromToIndexes chrom_indexes;
  {
    auto reader = io::VcfReader(vcf_path);
    // if vcf.gz, then indexes may be different from the header
    chrom_indexes = reader.GetContigIndexes();
    header_info = GetVcfHeaderInfo(reader.GetHeader());
  }
  const ChromToLengths chrom_lengths = header_info.contig_lengths;
  if (auto num_warnings = CheckVcfFeatureResources(feat_cols, interest_regions, header_info, has_popaf_vcf);
      num_warnings > 0) {
    WarnAsErrorIfSet(
        "There were {} warnings while checking VCF feature resources; VCF feature extraction may not be accurate!",
        num_warnings);
  }

  std::map<s32, std::string> chrom_seqs;
  vec<RegionInfo> region_infos;
  const bool is_germline_tagging = model_config.workflow == Workflow::kGermlineTagging;
  if (target_regions.empty()) {
    // no target regions were specified
    // will assign one parallel task for each chromosome
    // load all chromosome sequences in the reference genome
    io::FastaReader fasta_reader(param.genome);
    for (const auto& [chrom, chrom_length] : chrom_lengths) {
      if (chrom_indexes.contains(chrom)) {
        const s32 cid = chrom_indexes.at(chrom);
        chrom_seqs[cid] = fasta_reader.GetSequence(chrom);
        region_infos.emplace_back(chrom, cid, 0, chrom_length, has_popaf_vcf, chrom_seqs[cid], is_germline_tagging);
      }
    }
  } else {
    // one or more target intervals were specified
    // will assign one parallel task for each target chromosome
    // load only target chromosome sequences in the reference genome
    io::FastaReader fasta_reader(param.genome);
    for (const auto& [chrom, intervals] : target_regions) {
      if (chrom_indexes.contains(chrom)) {
        const s32 cid = chrom_indexes.at(chrom);
        chrom_seqs[cid] = fasta_reader.GetSequence(chrom);
        region_infos.emplace_back(chrom,
                                  cid,
                                  intervals.begin()->start,
                                  (intervals.end() - 1)->end,
                                  has_popaf_vcf,
                                  chrom_seqs[cid],
                                  is_germline_tagging);
      }
    }
  }

  StrMap<vec<VcfFeature>> chrom_to_features;
  std::mutex chrom_features_mutex;

  auto task = [&vcf_path, &target_regions, &interest_regions, &header_info, &chrom_features_mutex, &chrom_to_features](
                  const RegionInfo& info) {
    try {
      auto reader = io::VcfReader(vcf_path);
      const auto features = ComputeVcfFeaturesForRegion(info, reader, target_regions, interest_regions, header_info);
      if (!features.empty()) {
        // do not store the vector if it is empty
        const std::scoped_lock lock{chrom_features_mutex};
        chrom_to_features[info.chrom] = features;
      }
    } catch (const std::exception& e) {
      throw error::Error(
          "Error extracting feature in region '{}:{}-{}': {}", info.chrom, info.start, info.end, e.what());
    }
  };

  const size_t num_workers = param.threads <= 1 ? 1 : param.threads;
  Logging::Info("use {} threads to extract features", num_workers);
  tf::Executor executor(num_workers);
  tf::Taskflow flow;
  tf::Task map_task = flow.for_each(region_infos.begin(), region_infos.end(), task);
  map_task.name("Map VCF feature extraction task to region");
  executor.run(flow).get();

  // Sort features for each chromosome. Calculating variant density would require the vector of features to be sorted.
  for (auto& [chrom, features] : chrom_to_features) {
    SortByPosition(features);
  }

  return chrom_to_features;
}

/**
 * @brief Extract a map of chromosome to vector of intervals based on VCF feature positions.
 * @param chrom_to_features Map of chromosome to vector of VCF features
 * @param left_pad Padding to the left of variant start position
 * @param right_pad Padding to the right of variant end position
 * @param collapse_dist Distance to collapse neighboring intervals
 * @return Map of chromosome to vector of intervals
 */
static ChromIntervalsMap ExtractFeatureIntervalMap(const StrMap<vec<VcfFeature>>& chrom_to_features,
                                                   const u64 left_pad,
                                                   const u64 right_pad,
                                                   const u64 collapse_dist) {
  // Overview:
  // 1. extract feature positions from VCF features and create intervals
  // 2. merge overlapping intervals and collapse neighboring intervals
  // 3. store all intervals in a map of chromosome to vector of intervals
  ChromIntervalsMap output_chrom_intervals;
  for (const auto& [chrom, features] : chrom_to_features) {
    vec<Interval> intervals;
    intervals.reserve(features.size());
    for (const auto& f : features) {
      const u64 start = f.pos >= left_pad ? f.pos - left_pad : 0;
      const u64 end = f.pos + f.ref.length() + right_pad;
      intervals.emplace_back(start, end);
    }
    std::ranges::sort(intervals);
    output_chrom_intervals[chrom] = MergeIntervals(intervals, collapse_dist);
  }
  return output_chrom_intervals;
}

void ComputeVcfFeatures(const ComputeVcfFeaturesParam& param) {
  // Overview:
  // 1. extract features from the input VCF file and reference genome FASTA file
  // 2. update variant density for all extracted features
  // 3. update POPAF values for all extracted features, if POPAF VCF file is available
  // 4. write feature positions to a BED file, if BED file path is specified
  // 5. write features to output TSV file, if output file path is specified

  const SVCConfig model_config = param.config;
  const auto threads = static_cast<u32>(param.threads);

  auto chrom_to_features = ExtractFeatures(param, model_config.vcf_feature_cols);
  if (chrom_to_features.empty()) {
    WarnAsErrorIfSet("There are no variants with extracted features!");
  } else {
    UpdateVariantDensity(chrom_to_features);

    if (param.pop_af_vcf.has_value()) {
      Logging::Info("updating POPAF values...");
      UpdatePopaf(chrom_to_features, param.pop_af_vcf.value(), threads);
    }
  }

  if (param.output_bed.has_value()) {
    CreateParentDirectoryIfNotExists(param.output_bed.value());
    Logging::Info("writing features BED file...");
    const auto chrom_to_intervals =
        ExtractFeatureIntervalMap(chrom_to_features, param.left_pad, param.right_pad, param.collapse_dist);
    WriteBed(chrom_to_intervals, param.output_bed.value());
  }

  if (!param.output_file.empty()) {
    CreateParentDirectoryIfNotExists(param.output_file);
    Logging::Info("writing features file...");
    LockedTsvWriter writer(param.output_file);
    WriteAllFeatures(writer,
                     model_config.vcf_feature_names,
                     model_config.vcf_feature_cols,
                     chrom_to_features,
                     threads,
                     param.command_line);
  }

  Logging::Info("Completed VCF feature extraction");
}

}  // namespace xoos::svc
