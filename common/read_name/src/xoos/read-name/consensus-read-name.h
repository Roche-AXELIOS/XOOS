#pragma once

#include <optional>
#include <string>
#include <string_view>

#include <xoos/types/int.h>

namespace xoos::read_name {

// A post-inter-molecular-consensus read name encodes the composition of the cluster it was
// collapsed from. It splits into 6 hyphen-delimited fields:
//   {cluster_id}-{fwd_partial}-{rev_partial}-{fwd_full}-{rev_full}-{avg_family_size}
// where {cluster_id} is itself {super_region_id}:{count}. The five trailing fields are all numeric.
constexpr size_t kConsensusReadNameFieldCount = 6;

/**
 * Composition of the read cluster from which a consensus read was generated, parsed from the
 * consensus read name. The counts describe how many partial/full reads on each strand contributed
 * to the cluster, and `avg_family_size` is the mean cluster (family) size.
 */
struct ConsensusReadInfo {
  u32 forward_partial{};
  u32 reverse_partial{};
  u32 forward_full{};
  u32 reverse_full{};
  u32 avg_family_size{};
};

/**
 * @brief Determine whether a read name follows the post-inter-molecular-consensus format.
 *
 * A consensus read name splits into exactly six hyphen-delimited fields, of which the last five
 * (all but the cluster id) are numeric:
 *   {cluster_id}-{fwd_partial}-{rev_partial}-{fwd_full}-{rev_full}-{avg_family_size}
 *
 * @param qname The query name to inspect.
 * @return true if @p qname matches the consensus read name format, false otherwise.
 */
bool IsConsensusReadName(const std::string_view& qname);

/**
 * @brief Parse a post-inter-molecular-consensus read name into its cluster composition fields.
 *
 * @param qname The consensus read name, formatted as
 *   {cluster_id}-{fwd_partial}-{rev_partial}-{fwd_full}-{rev_full}-{avg_family_size}.
 * @return A `ConsensusReadInfo` with the parsed counts.
 * @throws error::Error If @p qname does not follow the consensus read name format.
 */
ConsensusReadInfo ParseConsensusReadName(const std::string_view& qname);

/**
 * @brief A no-throw wrapper around `ParseConsensusReadName`.
 *
 * @param qname The consensus read name to parse.
 * @return The parsed `ConsensusReadInfo`, or `std::nullopt` if @p qname is not a consensus read name.
 */
std::optional<ConsensusReadInfo> ParseConsensusReadNameNoExcept(const std::string_view& qname) noexcept;

/**
 * @brief Build a consensus read name from its cluster id and composition fields.
 *
 * Produces the canonical format
 *   {super_region_id}:{count}-{fwd_partial}-{rev_partial}-{fwd_full}-{rev_full}-{avg_family_size}
 * used by the Read Collapser when emitting inter-molecular consensus reads.
 *
 * @param super_region_id The batch (super region) id component of the cluster id.
 * @param count The cluster index within the super region.
 * @param info The cluster composition counts and average family size.
 * @return The formatted consensus read name.
 */
std::string FormatConsensusReadName(u32 super_region_id, u32 count, const ConsensusReadInfo& info);

}  // namespace xoos::read_name
