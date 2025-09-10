#include "compute-bam-region-features.h"

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>

#include "compute-bam-features/alignment-reader.h"
#include "compute-bam-features/read-processing.h"
#include "core/config.h"
#include "core/variant-info-serializer.h"
#include "util/log-util.h"

namespace xoos::svc {

ComputeBamRegionFeatures::ComputeBamRegionFeatures(const vec<AlignmentReader>& alignment_readers,
                                                   ComputeBamFeaturesParams params,
                                                   LockedTsvWriter* writer)
    : _alignment_readers(alignment_readers), _params(std::move(params)), _writer(writer) {
}

void ComputeBamRegionFeatures::ResetForRegion() {
  _read_names.clear();
  _current_read_id = 0;
  _var_features = {};
  _ref_features.clear();
  _vcf_features.clear();
  _vcf_positions.clear();
}

BamFeatures ComputeBamRegionFeatures::ComputeBamFeatures(
    const Region& region,
    const std::string& ref_sequence,
    const std::optional<ChromToVcfFeaturesMap>& chrom_to_vcf_features,
    const std::optional<StrUnorderedSet>& skip_variants) {
  // This function is designed to be used for feature extraction for a given genomic region. It can be run in parallel
  // or in sequence to compute features for different genomic regions with a region being as large as a full chromosome.
  // This function is also designed to take in multiple BAM files, with each BAM file's alignments to the given region
  // checked to identify all alignments that overlap the region over all input BAM files.
  ResetForRegion();

  // Check if there are VCF Features provided. If VCF features are passed we want to limit the feature computation
  // only to positions where there are observed variants, as only features at these positions will be used in
  // subsequent training and filtering applications.
  if (chrom_to_vcf_features.has_value()) {
    _vcf_features = chrom_to_vcf_features->at(region.chrom);
    if (_vcf_features.empty()) {
      // VCF Features were specified, but no features for this chromosome. No need to parse region
      return std::make_pair(_var_features, _ref_features);
    } else {
      // Have VCF Features for the chromosome, check if they intersect the region at all
      // Pad by one to ensure that edge cases are covered
      for (const auto& [pos, feats] : _vcf_features) {
        if (region.start <= pos - 1 && region.end >= pos + 1) {
          _vcf_positions.insert(pos);
        }
      }
      if (_vcf_positions.empty()) {
        // VCF Features were specified, but no features for this region. No need to parse region
        return std::make_pair(_var_features, _ref_features);
      }
    }
  }

  // loop through all BAM files and extract features from the alignments overlapping the target region. BAM files are
  // processed individually but feature values are cumulatively incremented with alignments across all BAM files.
  for (const auto& alignment_reader : _alignment_readers) {
    ComputeFeaturesForBam(region, ref_sequence, skip_variants, alignment_reader);
  }
  ProcessVariantFeatures();
  return std::make_pair(_var_features, _ref_features);
}

void ComputeBamRegionFeatures::ComputeFeaturesForBam(const Region& region,
                                                     const std::string& ref_sequence,
                                                     const std::optional<StrUnorderedSet>& skip_variants,
                                                     const AlignmentReader& alignment_reader) {
  // Add one to the end position to make it end inclusive. This is so that we correctly include reads that start
  // at the end position computing depth for insertion at the end
  const auto tid = io::SamHdrName2Tid(alignment_reader.hdr.get(), region.chrom);
  const auto alignments = io::SamItrQueryI(alignment_reader.idx.get(), tid, region.Start(), region.EndInclusive());
  if (alignments == nullptr) {
    throw error::Error("Failed to query region {}:{}-{}", region.chrom, region.start, region.end);
  }

  // Iterate through the specified region of bam file
  io::Bam1Ptr b(bam_init1());
  if (b == nullptr) {
    throw error::Error("Failed to allocate memory for bam record");
  }
  while (true) {
    const bool has_next = io::SamItrNext(alignment_reader.fp.get(), alignments.get(), b.get());
    if (!has_next) {
      break;
    }
    // skip reads that are unmapped, secondary, supplementary, duplicates, or have low mapping quality
    if ((b->core.flag & BAM_FUNMAP) != 0 || (b->core.flag & BAM_FSECONDARY) != 0 ||
        (b->core.flag & BAM_FSUPPLEMENTARY) != 0 || (b->core.flag & BAM_FDUP) != 0 || b->core.qual < _params.min_mapq ||
        b->core.tid < 0 || b->core.pos < 0) {
      continue;
    }
    // Sanity check to skip any alignments that start and end outside the region of interest.
    if (std::cmp_greater(b->core.pos, region.end)) {
      continue;
    }

    const std::string read_name = bam_get_qname(b.get());
    // exclude artificial haplotypes generated by GATK
    if (read_name.starts_with("HC")) {
      continue;
    }

    ComputeFeaturesForRead(region, ref_sequence, skip_variants, b.get());
  }
}

void ComputeBamRegionFeatures::ComputeFeaturesForRead(const Region& region,
                                                      const std::string& ref_sequence,
                                                      const std::optional<StrUnorderedSet>& skip_variants,
                                                      const bam1_t* b) {
  if (!_params.decode_yc) {
    // do not decode the YC tag; determine base types using base quality
    ComputeBamFeaturesHelper(b, ref_sequence, region, _read_names, _current_read_id, skip_variants);
  } else {
    // decode the YC tag
    auto decoded_records = yc_decode::DecodeToBamRecords(b);
    // TODO: update `plus_counts` and `minus_counts`?
    // This can get complicated in function `IncrementUnifiedFeature` within read-processing.cpp
    if (std::holds_alternative<std::monostate>(decoded_records)) {
      // no decoded records; the YC tag is either absent or empty
      // if input data is a hybrid of duplex reads and SBX-S reads, then this is an SBX-S read
      // we still want to decode YC tag (even if absent) to avoid using base quality to determine base types
      if (kSimplexReadFamilySize >= _params.min_family_size) {
        ComputeBamFeaturesHelper(b, ref_sequence, region, _read_names, _current_read_id, skip_variants);
      }
    } else if (std::holds_alternative<io::Bam1Ptr>(decoded_records)) {
      // one decoded record
      if (kSimplexReadFamilySize >= _params.min_family_size) {
        const auto r1 = std::move(std::get<io::Bam1Ptr>(decoded_records));
        ComputeBamFeaturesHelper(r1.get(), ref_sequence, region, _read_names, _current_read_id, skip_variants);
      }
    } else {
      // two decoded records
      if (kDuplexReadFamilySize >= _params.min_family_size) {
        const auto [r1, r2] = std::move(std::get<std::pair<io::Bam1Ptr, io::Bam1Ptr>>(decoded_records));
        ComputeBamFeaturesHelper(r1.get(), ref_sequence, region, _read_names, _current_read_id, skip_variants);
        ComputeBamFeaturesHelper(r2.get(), ref_sequence, region, _read_names, _current_read_id, skip_variants);
      }
    }
  }
}

void ComputeBamRegionFeatures::ProcessVariantFeatures() {
  // We only want to compute certain somatic tumor normal related features if we are computing features in the somatic
  // tumor normal workflow and have a set tumor read group.
  const bool tally_tumor_normal_features = _params.tumor_read_group.has_value() && !_params.tumor_read_group->empty();

  // At this stage, the features for all alignments in the target region have been computed. Although all features for
  // the variant/reference alleles are computed at each alignment position, only the user-selected feature columns
  // would be converted to their string representations.
  vec<vec<std::string>> region_feature_strings;
  if (_writer != nullptr) {
    region_feature_strings.reserve(_var_features.size());
  }
  // Gather variants at their reference feature positions. Since ref feature positions are offset by 1 for deletions,
  // gathering variants at their positions (instead of ref feature position) would lead to an incorrect tally.
  std::map<u64, vec<VariantId>> ref_pos_to_vids{};
  for (const auto& [id, f] : _var_features) {
    ref_pos_to_vids[id.GetRefFeaturePos()].emplace_back(id);
  }
  // Process variant features at each reference position and use helper function to accumulate and increment feature
  // values
  for (auto& [pos, vids] : ref_pos_to_vids) {
    if (vids.size() == 1) {
      const auto& vid = vids.front();
      auto& var_feat = _var_features.at(vid);
      UpdateDerivedFeatureValues(var_feat);
      auto& ref_feat = _ref_features[vid.GetRefFeaturePos()];
      ref_feat.num_alt = 1;
      UpdateDerivedFeatureValues(ref_feat);
      // Compute the total duplex dp here over all variant and ref positions with regular and lowbq
      ref_feat.duplex_dp = var_feat.duplex_lowbq + var_feat.duplex + ref_feat.duplex_lowbq + ref_feat.support;
      UpdateFeaturesHelper(var_feat,
                           ref_feat,
                           var_feat.mapq_sum + ref_feat.nonhomopolymer_mapq_sum,
                           var_feat.baseq_sum + ref_feat.nonhomopolymer_baseq_sum,
                           tally_tumor_normal_features);
      if (_writer != nullptr) {
        auto feature_strings = SerializeFeatureHelper(vid, var_feat, ref_feat, _params.feature_cols);
        region_feature_strings.emplace_back(feature_strings);
      }
    } else {
      // Multiple variants at ths position
      // Tally total MAPQ sum and BASEQ sum for all variant and reference alleles.
      auto& ref_feat = _ref_features[vids[0].GetRefFeaturePos()];
      UpdateDerivedFeatureValues(ref_feat);
      u32 total_mapq_sum = ref_feat.nonhomopolymer_mapq_sum;
      double total_baseq_sum = ref_feat.nonhomopolymer_baseq_sum;
      ref_feat.duplex_dp = ref_feat.duplex_lowbq + ref_feat.support;
      ref_feat.num_alt = static_cast<u32>(vids.size());
      for (const auto& vid : vids) {
        const auto& var_feat = _var_features.at(vid);
        total_mapq_sum += var_feat.mapq_sum;
        total_baseq_sum += var_feat.baseq_sum;
        // Compute the total duplex dp here over all variant and ref positions with regular and lowbq
        ref_feat.duplex_dp += var_feat.duplex_lowbq + var_feat.duplex;
      }
      // Update features for each variant and convert feature values to their string representations.
      for (const auto& vid : vids) {
        auto& var_feat = _var_features.at(vid);
        UpdateDerivedFeatureValues(var_feat);
        UpdateFeaturesHelper(var_feat, ref_feat, total_mapq_sum, total_baseq_sum, tally_tumor_normal_features);
        if (_writer != nullptr) {
          auto feature_strings = SerializeFeatureHelper(vid, var_feat, ref_feat, _params.feature_cols);
          region_feature_strings.emplace_back(feature_strings);
        }
      }
    }
  }
  if (_writer != nullptr) {
    _writer->AppendRows(region_feature_strings);
  }
}

void ComputeBamRegionFeatures::PopulateRegionsQueue(const StrMap<RefRegion>& contig_map,
                                                    const StrMap<vec<Interval>>& bed_regions,
                                                    vec<Region>& regions_queue,
                                                    const u64 max_region_size_per_thread) {
  if (!bed_regions.empty()) {
    GetGenomicRegionsForBedFile(contig_map, bed_regions, regions_queue, max_region_size_per_thread);
  } else {
    GetGenomicRegionsForChromosome(contig_map, regions_queue, max_region_size_per_thread);
  }
}

void ComputeBamRegionFeatures::GetGenomicRegionsForBedFile(const StrMap<RefRegion>& contig_map,
                                                           const StrMap<vec<Interval>>& bed_regions,
                                                           vec<Region>& regions_queue,
                                                           const u64 max_region_size_per_thread) {
  // This function breaks up BED entries per contig into smaller regions of at most a set size that can be used for
  // parallel processing.
  std::optional<std::string> prev_chrom = std::nullopt;
  std::optional<Interval> prev_interval = std::nullopt;
  for (const auto& [chrom, intervals] : bed_regions) {
    // skip if the contig has no alignments
    if (!contig_map.contains(chrom)) {
      continue;
    }

    if (!prev_chrom.has_value() || prev_chrom.value() != chrom) {
      // this is the first interval in the chromosome; no previous interval
      prev_interval = std::nullopt;
    }
    prev_chrom = chrom;

    // Add intervals to the list of regions. If the interval is too big break up the interval into smaller chunks
    for (const auto& interval : intervals) {
      const u64 region_start = interval.start;
      const u64 region_end = interval.end;
      for (u64 start = region_start; start <= region_end + 1; start += max_region_size_per_thread) {
        const u64 end = std::min(start + max_region_size_per_thread, region_end + 1);
        regions_queue.emplace_back(chrom, start, end, prev_interval);
        prev_interval = Interval{start, end};
      }
    }
  }
}

void ComputeBamRegionFeatures::GetGenomicRegionsForChromosome(const StrMap<RefRegion>& contig_map,
                                                              vec<Region>& regions_queue,
                                                              const u64 max_region_size_per_thread) {
  // This function breaks up contigs into smaller regions up to a certain size, which are used to compute features in
  // parallel
  std::optional<std::string> prev_chrom = std::nullopt;
  std::optional<Interval> prev_interval = std::nullopt;
  for (const auto& [chrom, region] : contig_map) {
    // leftmost aligned position of the contig
    const u64 region_start = region.start_position;
    // rightmost aligned position of the contig
    const u64 region_end = region.end_position;

    if (!prev_chrom.has_value() || prev_chrom.value() != chrom) {
      // this is the first interval in the chromosome; no previous interval
      prev_interval = std::nullopt;
    }
    prev_chrom = chrom;

    // Add regions to the list up to the max region size. Individual threads will run BAM feature computation per each
    // of these regions.
    for (u64 start = region_start; start <= region_end + 1; start += max_region_size_per_thread) {
      const u64 end = std::min(start + max_region_size_per_thread, region_end + 1);
      regions_queue.emplace_back(chrom, start, end, prev_interval);
      prev_interval = Interval{start, end};
    }
    // TODO: add dynamic partitioning that breaks up dense/high coverage regions
  }
}

void ComputeBamRegionFeatures::UpdateFeaturesHelper(UnifiedVariantFeature& var_feat,
                                                    UnifiedReferenceFeature& ref_feat,
                                                    const u32 total_mapq_sum,
                                                    const double total_baseq_sum,
                                                    const bool tally_tumor_normal_features) {
  // Helper function to update feature values. This helper function is reliant on both the variant specific and
  // reference
  // position specific values and totals. Thus, it should be run only once all alignments have been processed and read
  // alignments supporting the reference or given variant at the position noted and used to update the respective
  // feature values. Begin by updating outstanding features
  var_feat.mq_af = total_mapq_sum > 0 ? static_cast<double>(var_feat.mapq_sum) / total_mapq_sum : 0;
  var_feat.bq_af = total_baseq_sum > 0 ? var_feat.baseq_sum / total_baseq_sum : 0;
  ref_feat.mq_af = total_mapq_sum > 0 ? static_cast<double>(ref_feat.nonhomopolymer_mapq_sum) / total_mapq_sum : 0;
  ref_feat.bq_af = total_baseq_sum > 0 ? static_cast<double>(ref_feat.nonhomopolymer_baseq_sum) / total_baseq_sum : 0;
  var_feat.duplex_af = ref_feat.duplex_dp > 0 ? (var_feat.duplex_lowbq + var_feat.duplex) / ref_feat.duplex_dp : 0;
  ref_feat.duplex_af = ref_feat.duplex_dp > 0 ? (ref_feat.duplex_lowbq + ref_feat.support) / ref_feat.duplex_dp : 0;
  if (var_feat.support > 0) {
    const double f = static_cast<double>(var_feat.support_reverse) / var_feat.support;
    var_feat.alignmentbias = 0.25 - f * (1 - f);
  }
  const u32 duplex_sum = var_feat.duplex + ref_feat.nonhomopolymer_support;
  if (duplex_sum > 0) {
    auto duplex = static_cast<double>(var_feat.duplex);
    var_feat.adt = std::abs(duplex - ref_feat.nonhomopolymer_support) / duplex_sum;
    var_feat.adtl = std::log10(duplex_sum) * var_feat.adt;
    var_feat.indel_af = duplex / duplex_sum;
  }
#ifdef SOMATIC_ENABLE
  // Only compute these tumor specific features if running feature computation for somatic tumor-normal
  if (tally_tumor_normal_features) {
    auto total_tumor_support = var_feat.tumor_support + ref_feat.tumor_support;
    var_feat.tumor_af = total_tumor_support > 0 ? static_cast<double>(var_feat.tumor_support) / total_tumor_support : 0;
    auto total_normal_support = var_feat.normal_support + ref_feat.normal_support;
    var_feat.normal_af =
        total_normal_support > 0 ? static_cast<double>(var_feat.normal_support) / total_normal_support : 0;
    var_feat.rat =
        var_feat.tumor_af + var_feat.normal_af > 0 ? var_feat.tumor_af / (var_feat.tumor_af + var_feat.normal_af) : 0;
    if (var_feat.tumor_support > 0) {
      double f = static_cast<double>(var_feat.tumor_support_reverse) / var_feat.tumor_support;
      var_feat.tumor_alignmentbias = 0.25 - f * (1 - f);
    }
    if (var_feat.normal_support > 0) {
      double f = static_cast<double>(var_feat.normal_support_reverse) / var_feat.normal_support;
      var_feat.normal_alignmentbias = 0.25 - f * (1 - f);
    }
  }
#endif  // SOMATIC_ENABLE
}

vec<std::string> ComputeBamRegionFeatures::SerializeFeatureHelper(const VariantId& vid,
                                                                  const UnifiedVariantFeature& var_feat,
                                                                  const UnifiedReferenceFeature& ref_feat,
                                                                  const vec<UnifiedFeatureCols>& feature_cols) {
  // Helper function to convert selected feature values to strings. This function is used primarily for generating
  // strings to be output to file. First convert all feature structs to their string representations, and combine into a
  // single map.
  auto ref_feature_serialized = VariantInfoSerializer::SerializeReferenceFeature(ref_feat);
  auto variant_id_serialized = VariantInfoSerializer::SerializeVariantId(vid);
  auto variant_feature_serialized = VariantInfoSerializer::SerializeVariantFeature(var_feat);
  variant_feature_serialized.insert(variant_id_serialized.begin(), variant_id_serialized.end());
  variant_feature_serialized.insert(ref_feature_serialized.begin(), ref_feature_serialized.end());
  // Only the feature columns selected by the user are stored. Ensure output vector matches specified feature order
  vec<std::string> feature_strings;
  feature_strings.reserve(feature_cols.size());
  for (const auto& col : feature_cols) {
    feature_strings.push_back(variant_feature_serialized[col]);
  }
  return feature_strings;
}

void ComputeBamRegionFeatures::ComputeBamFeaturesHelper(const bam1_t* record,
                                                        const std::string& ref_seq,
                                                        const Region& region,
                                                        ReadNameToId& read_names,
                                                        ReadId& current_read_id,
                                                        const std::optional<StrUnorderedSet>& skip_variants) {
  // This function extracts features for an individual BAM record. It assumes that the record has already been checked
  // for mapping quality and flag stats. This function should only be used for records that have passed checks for
  // minimum mapping quality, flag stats, and are aligned within the region of interest. First extract alignment
  // positions for all CIGAR operators in the BAM record determine the alignment's end position on the reference

  if (record->core.pos < 0) {
    WarnAsErrorIfSet("Found alignment with start position < 0: '{}'", bam_get_qname(record));
    return;
  }

  const auto& align_ops = GetAlignOpInfos(record, ref_seq);
  if (align_ops.empty()) {
    // no alignment information; skip this record
    return;
  }

  const auto bam_pos = static_cast<u64>(record->core.pos);
  const auto reference_end = align_ops.back().ref_pos + align_ops.back().op_len;
  // Sanity check to make sure the read doesn't start and end before the region of interest
  if (reference_end < region.start || bam_pos >= region.end) {
    return;
  }
  if (!_vcf_positions.empty()) {
    // check whether the alignment intersects the VCF positions for this region. We only want to extract features from
    // reads that overlap an observed variant in the VCF positions set
    auto lower = _vcf_positions.lower_bound(bam_pos);
    auto upper = _vcf_positions.upper_bound(reference_end);
    if (lower == upper) {
      return;
    }
  }

  // check whether thresholds are set for max number of variants allowed per alignment
  const bool check_num_var = _params.max_read_variant_count.has_value() && _params.max_read_variant_count.value() > 0;
  const bool check_num_var_norm =
      _params.max_read_variant_count_normalized.has_value() && _params.max_read_variant_count_normalized.value() > 0;
  if (check_num_var || check_num_var_norm) {
    u32 num_var{0};
    if (skip_variants.has_value() && !skip_variants->empty()) {
      // count variants in the alignment, ignoring any known variants
      num_var = CountVariantsInRead(align_ops, record, ref_seq, region.chrom, skip_variants.value());
    } else {
      // naively count variants in the alignment
      num_var = CountVariantsInRead(align_ops);
    }
    if ((check_num_var && num_var > _params.max_read_variant_count.value()) ||
        (check_num_var_norm && (static_cast<float>(num_var) / static_cast<float>(GetAlignmentLength(align_ops))) >
                                   _params.max_read_variant_count_normalized.value())) {
      // the alignment has too many variants; do not extract any features from it
      return;
    }
  }
  // extract a numeric ID from the read name
  const auto read_id = GetReadId(read_names, bam_get_qname(record), current_read_id);

  // extract features for variants and reference alleles using the extracted CIGAR operator alignment info. Computed
  // features are added to the maps that has been passed in.
  const AlignContext align_ctx(_params, align_ops, record, region, ref_seq, read_id, _vcf_features);
  ProcessAlignment(align_ctx, _var_features, _ref_features);
}

ReadId ComputeBamRegionFeatures::GetReadId(ReadNameToId& read_names,
                                           const std::string& read_name,
                                           ReadId& current_read_id) {
  auto it = read_names.find(read_name);
  if (it != read_names.end()) {
    return it->second;
  }
  auto read_id = current_read_id;
  current_read_id++;
  read_names.try_emplace(read_name, read_id);
  return read_id;
}

}  // namespace xoos::svc
