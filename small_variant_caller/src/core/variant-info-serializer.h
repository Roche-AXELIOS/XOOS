#pragma once

#include <fstream>
#include <map>
#include <string>

#include <xoos/types/fs.h>

#include "bam-feature-collection.h"
#include "xoos/types/float.h"
#include "xoos/types/vec.h"

namespace xoos::svc {

/**
 * Synopsis:
 * This file provides serialization and deserialization methods for variant features,
 * reference features, and variant IDs. It also includes methods to load features from files
 * and numericalize feature values.
 */

class VariantInfoSerializer {
 public:
  /**
   * @brief Deserialize a row of feature value strings into a VariantId struct
   * @param header a vector of FeatureColumn structs
   * @param fields a vector of string values corresponding to the FeatureColumn structs
   * @return a VariantId Struct
   */
  static VariantId DeserializeVariantId(const vec<FeatureColumn>& header, const vec<std::string>& fields);

  /**
   * @brief Deserialize a row of feature value strings into a BamFeatureTuple struct
   * @param header a vector of FeatureColumn structs
   * @param fields a vector of string values corresponding to the FeatureColumn structs
   * @return a BamFeatureTuple Struct
   */
  static BamFeatureTuple DeserializeBamFeatureRow(const vec<FeatureColumn>& header, const vec<std::string>& fields);

  /**
   * @brief Deserialize a row of feature value strings into a TumorNormalBamFeatureTuple struct
   * @param header a vector of FeatureColumn structs
   * @param fields a vector of string values corresponding to the FeatureColumn structs
   * @return a TumorNormalBamFeatureTuple Struct
   */
  static TumorNormalBamFeatureTuple DeserializeTumorNormalBamFeatureRow(const vec<FeatureColumn>& header,
                                                                        const vec<std::string>& fields);

  /**
   * @brief Parse the header of a BAM features file and return a vector of FeatureColumn
   * @param input Input file stream
   * @param features_file Path to BAM features file
   * @param infer_sample_context Whether to infer sample context from feature name prefix (e.g., "tumor_", "normal_")
   * @return Vector of FeatureColumn
   */
  static vec<FeatureColumn> ParseBamFeaturesFileHeader(std::ifstream& input,
                                                       const fs::path& features_file,
                                                       bool infer_sample_context);

  /**
   * @brief Parse the header of a VCF features file and return a vector of FeatureColumn
   * @param input Input file stream
   * @param features_file Path to VCF features file
   * @param infer_sample_context Whether to infer sample context from feature name prefix (e.g., "tumor_", "normal_")
   * @return Vector of FeatureColumn
   */
  static vec<FeatureColumn> ParseVcfFeaturesFileHeader(std::ifstream& input,
                                                       const fs::path& features_file,
                                                       bool infer_sample_context);

  /**
   * @brief Validate the header of a features file for feature normalization.
   * @param features_file Path to features file
   * @param has_sample_context Whether the features file has sample context in its feature columns (e.g., "tumor_",
   * "normal_")
   * @throws error::Error if required feature(s) for normalization are missing from the header
   */
  static void ValidateFeatureFileHeaderForNormalization(const fs::path& features_file, bool has_sample_context);

  /**
   * @brief Read a line from a TSV file, skip comments and empty lines, and validate the number of columns
   * @param input Input file stream
   * @param num_cols Expected number of columns
   * @return Vector of field values, or empty vector if end of file is reached
   */
  static vec<std::string> ParseNextTsvRow(std::ifstream& input, size_t num_cols);

  /**
   * @brief Load VCF features from a TSV file, filtering by a container for VariantId.
   * @tparam VariantIdContainer A container type that supports the `contains` method for VariantId.
   * @param features_file a text file containing BAM features
   * @param var_ids a container of VariantId to filter the loaded features
   * @param infer_sample_context Whether to infer sample context from feature name prefix (e.g., "tumor_", "normal_")
   * @return VcfFeaturesMap containing the loaded VCF features
   */
  template <typename VariantIdContainer>
    requires requires(const VariantIdContainer& m, const VariantId& id) { m.contains(id); }
  static VarIdToVcfFeatures LoadVcfFeatures(const fs::path& features_file,
                                            VariantIdContainer& var_ids,
                                            const bool infer_sample_context) {
    VarIdToVcfFeatures features;
    std::ifstream input(features_file);
    vec<std::string> fields;
    // extract header columns
    const vec<FeatureColumn>& header = ParseVcfFeaturesFileHeader(input, features_file, infer_sample_context);
    const auto num_cols = header.size();
    // parse rest of the file to extract VCF features
    while (!(fields = ParseNextTsvRow(input, num_cols)).empty()) {
      const VariantId vid = DeserializeVariantId(header, fields);
      // Skip this line if VariantId is not in the filter container
      if (!var_ids.empty() && !var_ids.contains(vid)) {
        continue;
      }
      features[vid] = DeserializeVcfFeatureRow(header, fields);
    }
    return features;
  }

  /**
   * @brief Load VCF features from a features file.
   * @param features_file a text file containing VCF features
   * @param infer_sample_context Whether to infer sample context from feature name prefix (e.g., "tumor_", "normal_")
   * @return a map of VariantIds to VcfFeatures
   */
  static VarIdToVcfFeatures LoadVcfFeatures(const fs::path& features_file, bool infer_sample_context);

  /**
   * @brief Serialize a UnifiedVariantFeature struct into a map of feature enum to string value
   * @param variant_info a UnifiedVariantFeature struct
   * @return a map of UnifiedFeatureCols to std::strings
   */
  static std::map<UnifiedFeatureCols, std::string> SerializeVariantFeature(const UnifiedVariantFeature& variant_info);

  /**
   * @brief Serialize a UnifiedReferenceFeature struct into a map of feature enum to string value
   * @param ref_info a UnifiedReferenceFeature struct
   * @return a map of UnifiedFeatureCols to std::strings
   */
  static std::map<UnifiedFeatureCols, std::string> SerializeReferenceFeature(const UnifiedReferenceFeature& ref_info);

  /**
   * @brief Serialize a VariantId struct into a map of feature enum to string value
   * @param vid a VariantId struct
   * @return a map of UnifiedFeatureCols to std::strings
   */
  static std::map<UnifiedFeatureCols, std::string> SerializeVariantId(const VariantId& vid);

  /**
   * @brief Serialize a BamFeatureTuple struct into a vector of string values based on the provided header
   * @param header a vector of FeatureColumn structs
   * @param vid a VariantId struct
   * @param feat_row a BamFeatureTuple struct
   * @return a vector of string values corresponding to the FeatureColumn structs
   */
  static vec<std::string> SerializeBamFeatureRow(const vec<FeatureColumn>& header,
                                                 const VariantId& vid,
                                                 const BamFeatureTuple& feat_row);

  /**
   * @brief Serialize a TumorNormalBamFeatureTuple struct into a vector of string values based on the provided header
   * @param header a vector of FeatureColumn structs
   * @param vid a VariantId struct
   * @param feat_row a TumorNormalBamFeatureTuple struct
   * @param has_tumor_feat whether tumor-specific features are present
   * @param has_normal_feat whether normal-specific features are present
   * @return a vector of string values corresponding to the FeatureColumn structs
   */
  static vec<std::string> SerializeTumorNormalBamFeatureRow(const vec<FeatureColumn>& header,
                                                            const VariantId& vid,
                                                            const TumorNormalBamFeatureTuple& feat_row,
                                                            bool has_tumor_feat,
                                                            bool has_normal_feat);

  /**
   * @brief Numericalize an attribute in the feature structs.
   * @param col Enum of the attribute
   * @param vid VariantId struct
   * @param bam_feat UnifiedVariantFeature struct
   * @param vcf_feat VcfFeature struct
   * @param ref_feat UnifiedReferenceFeature struct
   * @return double value of the attribute
   */
  static f64 NumericalizeFeature(UnifiedFeatureCols col,
                                 const VariantId& vid,
                                 const UnifiedVariantFeature& bam_feat,
                                 const VcfFeature& vcf_feat,
                                 const UnifiedReferenceFeature& ref_feat);

  /**
   * @brief Numericalize an attribute in a VariantId struct.
   * @param col Enum of the attribute
   * @param vid VariantId struct
   * @return double value of the attribute
   */
  static f64 NumericalizeFeature(UnifiedFeatureCols col, const VariantId& vid);

  /**
   * @brief Numericalize an attribute in a UnifiedVariantFeature struct.
   * @param col Enum of the attribute
   * @param feat UnifiedVariantFeature struct
   * @return double value of the attribute
   */
  static f64 NumericalizeFeature(UnifiedFeatureCols col, const UnifiedVariantFeature& feat);

  /**
   * @brief Numericalize an attribute in a VcfFeature struct.
   * @param col Enum of the attribute
   * @param feat VcfFeature struct
   * @return double value of the attribute
   */
  static f64 NumericalizeFeature(UnifiedFeatureCols col, const VcfFeature& feat);

  /**
   * @brief Numericalize an attribute in a UnifiedReferenceFeature struct.
   * @param col Enum of the attribute
   * @param feat UnifiedReferenceFeature struct
   * @return double value of the attribute
   */
  static f64 NumericalizeFeature(UnifiedFeatureCols col, const UnifiedReferenceFeature& feat);

  /**
   * @brief Deserialize a row of VCF feature file into a VcfFeature struct
   * @param header Vector of FeatureColumn structs extracted from the features file header line
   * @param fields Vector of string values corresponding to the FeatureColumn structs
   * @return VcfFeature struct
   */
  static VcfFeature DeserializeVcfFeatureRow(const vec<FeatureColumn>& header, const vec<std::string>& fields);

  /**
   * @brief Serialize a VcfFeatures struct into a map of VcfFeatureCols to string
   * @param feature a VcfFeatures struct
   * @return a map of VcfFeatureCols to string
   */
  static std::map<UnifiedFeatureCols, std::string> SerializeVcfFeature(const VcfFeature& feature);

  /**
   * @brief Generator function that yields (VariantId, BamFeatureTuple) pairs from a BAM feature file.
   *        The generator parses the file using the same logic as LoadBamFeatures, but yields one tuple at a time.
   * @param features_file Path to the BAM features file (TSV).
   * @return Generator yielding std::pair<VariantId, BamFeatureTuple> for each record.
   * @throws std::runtime_error if the file cannot be opened or parsing fails.
   */
  static auto BamFeatureTupleGenerator(const fs::path& features_file)
      -> std::function<std::optional<std::pair<VariantId, BamFeatureTuple>>()>;

  /**
   * @brief Generator function that yields (VariantId, TumorNormalBamFeatureTuple) pairs from a BAM feature file.
   *        The generator parses the file using the same logic as LoadBamFeatures, but yields one tuple at a time.
   * @param features_file Path to the BAM features file (TSV).
   * @return Generator yielding std::pair<VariantId, TumorNormalBamFeatureTuple> for each record.
   * @throws error::Error if the file cannot be opened or parsing fails.
   */
  static auto TumorNormalBamFeatureTupleGenerator(const fs::path& features_file)
      -> std::function<std::optional<std::pair<VariantId, TumorNormalBamFeatureTuple>>()>;
};

}  // namespace xoos::svc
