#pragma once

#include <string>
#include <tuple>

#include "xoos/io/vcf/vcf-record.h"
#include "xoos/types/int.h"
#include "xoos/types/vec.h"

namespace xoos::svc {

/**
 * @brief Get value from the INFO/FORMAT field.
 * @tparam T Type of the record value to extract (int or float)
 * @tparam U Type of the value to return
 * @param in_header Flag whether the field is present in the VCF header
 * @param values Vector of values from the INFO/FORMAT field
 * @param idx Index of the value in the field
 * @return Field value of the specified type U
 * @note If the field is not present in the header or the index is out of bounds, returns a default-constructed value of
 * type U. If the type T is the same as U, no casting is performed. If T is different from U, the value is cast to type
 * U.
 */
template <typename T, typename U>
constexpr U GetValue(const bool in_header, const vec<T>& values, const u32 idx) {
  if (in_header && std::cmp_less(idx, values.size())) {
    if (std::is_same_v<T, U>) {
      return values.at(idx);
    }
    return static_cast<U>(values.at(idx));
  }
  return U{};
}

/**
 * @brief Get value from the INFO field.
 * @tparam T Type of the record value to extract (int or float)
 * @tparam U Type of the value to return
 * @param in_header Flag whether the field is present in the VCF header
 * @param record VCF record to extract value from
 * @param idx Index of the value in the field
 * @param field Name of the field
 * @return Field value of the specified type U
 */
template <typename T, typename U>
constexpr U GetInfo(const bool in_header, const io::VcfRecordPtr& record, const u32 idx, const std::string& field) {
  return GetValue<T, U>(in_header, record->GetInfoFieldNoCheck<T>(field), idx);
}

/**
 * @brief Get value from the FORMAT field.
 * @tparam T Type of the record value to extract (int or float)
 * @tparam U Type of the value to return
 * @param in_header Flag whether the field is present in the VCF header
 * @param record VCF record to extract value from
 * @param idx Index of the value in the field
 * @param field Name of the field
 * @return Field Value of the specified type U
 */
template <typename T, typename U>
constexpr U GetFormat(const bool in_header, const io::VcfRecordPtr& record, const u32 idx, const std::string& field) {
  return GetValue<T, U>(in_header, record->GetFormatFieldNoCheck<T>(field), idx);
}

/**
 * @brief Get value for the INFO/FORMAT field for REF and ALT alleles.
 * @tparam T Type of the values to extract (int or float)
 * @tparam U Type of the values to return
 * @param in_header Flag whether the field is present in the VCF header
 * @param values Vector of values from the INFO/FORMAT field
 * @param idx Index of the value in the field
 * @return Field values of the specified type U for REF and ALT
 * @note If the field is not present in the header or the index is out of bounds, returns a tuple of default-constructed
 * values of type U. If the type T is the same as U, no casting is performed. If T is different from U, the values
 * are cast to type U.
 * @detail The first value in the tuple is for REF allele, the second value is for the ALT allele.
 */
template <typename T, typename U>
constexpr std::tuple<U, U> GetValueRefAlt(const bool in_header, const vec<T>& values, const u32 idx) {
  if (in_header && std::cmp_less(idx, values.size())) {
    if (std::is_same_v<T, U>) {
      return std::make_tuple(values.at(0), values.at(idx));
    }
    return std::make_tuple(static_cast<U>(values.at(0)), static_cast<U>(values.at(idx)));
  }
  return std::tuple<U, U>{};
}

/**
 * @brief Get values from the INFO field for both REF and ALT.
 * @tparam T Type of the values to extract (int or float)
 * @tparam U Type of the values to return
 * @param in_header Flag whether the field is present in the VCF header
 * @param record VCF record to extract values from
 * @param idx Index of the ALT allele
 * @param field Name of the field
 * @return Field values of the specified type U for REF and ALT
 */
template <typename T, typename U>
constexpr std::tuple<U, U> GetInfoRefAlt(const bool in_header,
                                         const io::VcfRecordPtr& record,
                                         const u32 idx,
                                         const std::string& field) {
  return GetValueRefAlt<T, U>(in_header, record->GetInfoFieldNoCheck<T>(field), idx);
}

/**
 * @brief Get REF and ALT values from the FORMAT field.
 * @tparam T Type of the values to extract (int or float)
 * @tparam U Type of the values to return
 * @param in_header Flag whether the field is present in the VCF header
 * @param record VCF record to extract values from
 * @param idx Index of the ALT allele
 * @param field Name of the field
 * @return Field values of the specified type U for REF and ALT
 */
template <typename T, typename U>
constexpr std::tuple<U, U> GetFormatRefAlt(const bool in_header,
                                           const io::VcfRecordPtr& record,
                                           const u32 idx,
                                           const std::string& field) {
  return GetValueRefAlt<T, U>(in_header, record->GetFormatFieldNoCheck<T>(field), idx);
}

/**
 * @brief Get REF and two ALT values from the FORMAT field for multi-allelic variants.
 * @tparam T Type of the value to extract (int or float)
 * @tparam U Type of the values to return
 * @param in_header Flag whether the field is present in the VCF header
 * @param record VCF record to extract values from
 * @param idx Index of the target ALT allele
 * @param field Name of the field
 * @return Field values of the specified type U for REF and two ALTs
 */
template <typename T, typename U>
constexpr std::tuple<U, U, U> GetFormatThreeAlleles(const bool in_header,
                                                    const io::VcfRecordPtr& record,
                                                    const u32 idx,
                                                    const std::string& field) {
  const auto& values = record->GetFormatFieldNoCheck<T>(field);
  if (in_header && std::cmp_less(idx, values.size())) {
    // the 3rd value is intended for the other ALT allele in multi-allelic variants
    // ALT idx = 1: REF, ALT1, ALT2
    // ALT idx = 2: REF, ALT2, ALT1
    const u32 other_idx = idx > 1 ? 1 : 2;
    if (std::is_same_v<T, U>) {
      return std::make_tuple(values.at(0), values.at(idx), values.at(other_idx));
    }
    return std::make_tuple(
        static_cast<U>(values.at(0)), static_cast<U>(values.at(idx)), static_cast<U>(values.at(other_idx)));
  }
  return std::tuple<U, U, U>{};
}

}  // namespace xoos::svc
