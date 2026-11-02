#include "features-to-vcf.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include <xoos/error/error.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-writer.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/string-functions.h>

#include "core/column-names.h"
#include "core/vcf-fields.h"
#include "util/file-util.h"

namespace xoos::svc {
namespace {

constexpr s32 kDummyAlleleCount = 1;
constexpr f32 kDummyAlleleFrequency = 0.5F;
constexpr s32 kDummyTotalAlleles = 2;
constexpr s32 kDummyReadDepth = 2;
constexpr s32 kDummyRefAlleleDepth = 1;
constexpr s32 kDummyAltAlleleDepth = 1;
constexpr s32 kDummyGenotypeQuality = 50;
const std::string kDummyGenotype = "0/1";
constexpr size_t kRequiredColumnCount = 4;
constexpr auto* kColumnDelimiter = "\t";

struct Allele {
  std::string ref;
  std::string alt;
};

using VariantMap = std::map<size_t, std::map<s32, vec<Allele>>>;

size_t FindColumnIndex(const vec<std::string>& header_fields, const std::string_view column_name) {
  const auto iter = std::ranges::find(header_fields, column_name);
  if (iter == header_fields.end()) {
    throw error::Error("Required column '{}' not found in input file header", column_name);
  }
  return static_cast<size_t>(std::distance(header_fields.begin(), iter));
}

std::array<size_t, kRequiredColumnCount> FindColumnIndices(const vec<std::string>& header_fields) {
  return {FindColumnIndex(header_fields, kNameChrom),
          FindColumnIndex(header_fields, kNamePos),
          FindColumnIndex(header_fields, kNameRef),
          FindColumnIndex(header_fields, kNameAlt)};
}

s32 ParseOneBasedPosition(const std::string& value, const size_t line_num) {
  try {
    size_t parsed_chars = 0;
    const auto position = std::stoll(value, &parsed_chars);
    if (parsed_chars != value.size() || position < 1 || position > std::numeric_limits<s32>::max()) {
      throw error::Error("Invalid 1-based position '{}' at line {}", value, line_num);
    }
    return static_cast<s32>(position - 1);
  } catch (const std::invalid_argument&) {
    throw error::Error("Invalid 1-based position '{}' at line {}", value, line_num);
  } catch (const std::out_of_range&) {
    throw error::Error("Invalid 1-based position '{}' at line {}", value, line_num);
  }
}

VariantMap ParseInputFile(const fs::path& input_path, const StrUnorderedMap<size_t>& contig_order) {
  std::ifstream input(input_path);
  if (!input.is_open()) {
    throw error::Error("Cannot open input file: {}", input_path.string());
  }

  std::string line;
  bool header_found = false;
  size_t line_num = 0;
  while (std::getline(input, line)) {
    ++line_num;
    if (line.empty() || line.starts_with('#')) {
      continue;
    }
    header_found = true;
    break;
  }
  if (!header_found) {
    throw error::Error("Cannot find header in input file: {}", input_path.string());
  }

  const vec<std::string> header_fields = string::Split(line, kColumnDelimiter);
  const auto [chrom_idx, pos_idx, ref_idx, alt_idx] = FindColumnIndices(header_fields);
  const size_t min_required_idx = std::max({chrom_idx, pos_idx, ref_idx, alt_idx});

  VariantMap variants;
  while (std::getline(input, line)) {
    ++line_num;
    if (line.empty() || line.starts_with('#')) {
      continue;
    }

    const vec<std::string> fields = string::Split(line, kColumnDelimiter);
    if (fields.size() <= min_required_idx) {
      throw error::Error("Line {} has {} columns, expected at least {}", line_num, fields.size(), min_required_idx + 1);
    }

    const auto& chrom = fields[chrom_idx];
    if (!contig_order.contains(chrom)) {
      Logging::Warn("Skipping variant on unknown contig '{}' at line {}", chrom, line_num);
      continue;
    }

    const size_t contig_idx = contig_order.at(chrom);
    variants[contig_idx][ParseOneBasedPosition(fields[pos_idx], line_num)].push_back(
        {.ref = fields[ref_idx], .alt = fields[alt_idx]});
  }

  return variants;
}

io::VcfHeaderPtr CreateHeader(const io::FastaReader& fasta_reader,
                              StrUnorderedMap<size_t>& contig_order,
                              vec<std::string>& contig_names,
                              const std::string& sample) {
  auto header = io::VcfHeader::Create();
  using enum io::FieldType;

  size_t contig_idx = 0;
  for (const auto& [name, length] : fasta_reader) {
    header->AddContigLine({.id = name, .length = static_cast<u32>(length)});
    contig_order[name] = contig_idx++;
    contig_names.push_back(name);
  }

  header->AddInfoLine(
      {.id = kFieldAc, .description = "Allele count in genotypes", .number = io::kNumberEachAllele, .type = kInteger});
  header->AddInfoLine(
      {.id = kFieldAf, .description = "Allele frequency", .number = io::kNumberEachAllele, .type = kFloat});
  header->AddInfoLine({.id = kFieldAn,
                       .description = "Total number of alleles in called genotypes",
                       .number = io::kNumberOne,
                       .type = kInteger});
  header->AddInfoLine({.id = kFieldDp, .description = "Total read depth", .number = io::kNumberOne, .type = kInteger});
  header->AddFormatLine({.id = kFieldGt, .description = "Genotype", .number = io::kNumberOne, .type = kString});
  header->AddFormatLine({.id = kFieldDp, .description = "Read depth", .number = io::kNumberOne, .type = kInteger});
  header->AddFormatLine({.id = kFieldAd,
                         .description = "Allelic depths for ref and alt alleles",
                         .number = io::kNumberR,
                         .type = kInteger});
  header->AddFormatLine(
      {.id = kFieldGq, .description = "Genotype quality", .number = io::kNumberOne, .type = kInteger});
  header->AddSample(sample);
  header->Sync();
  return header;
}

}  // namespace

void ConvertFeaturesToVcf(const FeaturesToVcfParam& param) {
  Logging::Info("Reading BAM features from '{}'", param.input_file.string());
  Logging::Info("Reading contig lengths from '{}'", param.genome.string());

  const io::FastaReader fasta_reader(param.genome.string());
  StrUnorderedMap<size_t> contig_order;
  vec<std::string> contig_names;
  const auto header = CreateHeader(fasta_reader, contig_order, contig_names, param.sample);
  const auto variants = ParseInputFile(param.input_file, contig_order);

  CreateParentDirectoryIfNotExists(param.output_file);

  const io::VcfIntegerFields empty_int_fields;
  const io::VcfFloatFields empty_float_fields;
  const io::VcfStringFields empty_string_fields;
  const io::TypedVcfFields empty_fields{
      .integer_fields = empty_int_fields, .float_fields = empty_float_fields, .string_fields = empty_string_fields};

  u64 record_count = 0;
  {
    io::VcfWriter writer(param.output_file, header);
    writer.WriteHeader();

    for (const auto& [contig_idx, positions] : variants) {
      const auto& chrom = contig_names[contig_idx];
      for (const auto& [position, alleles] : positions) {
        auto sorted_alleles = alleles;
        std::ranges::sort(sorted_alleles, [](const Allele& lhs, const Allele& rhs) { return lhs.alt < rhs.alt; });

        for (const auto& allele : sorted_alleles) {
          const auto record = writer.CreateRecord(
              chrom, position, ".", {allele.ref, allele.alt}, std::nullopt, "PASS", empty_fields, empty_fields);
          record->SetInfoField(kFieldAc, vec<s32>{kDummyAlleleCount});
          record->SetInfoField(kFieldAf, vec<f32>{kDummyAlleleFrequency});
          record->SetInfoField(kFieldAn, vec<s32>{kDummyTotalAlleles});
          record->SetInfoField(kFieldDp, vec<s32>{kDummyReadDepth});
          record->SetGTField(kDummyGenotype);
          record->SetFormatField(kFieldDp, vec<s32>{kDummyReadDepth});
          record->SetFormatField(kFieldAd, vec<s32>{kDummyRefAlleleDepth, kDummyAltAlleleDepth});
          record->SetFormatField(kFieldGq, vec<s32>{kDummyGenotypeQuality});
          writer.WriteRecord(record);
          ++record_count;
        }
      }
    }
    writer.Flush();
  }

  Logging::Info("Wrote {} records to '{}'", record_count, param.output_file.string());
}

}  // namespace xoos::svc
