#pragma once

#include <xoos/types/int.h>

#include "library-type.h"
#include "workflow.h"

namespace xoos::svc {

/**
 * @brief Concrete configuration profile used to select a complete set of SVC defaults.
 * @details Profile names are JSON keys and need not equal workflow names. In particular, both germline WGS profiles
 * execute Workflow::kGermlineWgs but retain independent duplex and simplex defaults.
 */
enum class ConfigProfile : u8 {
  kGermlineWgsDuplex,
  kGermlineWgsSimplex,
  kTumorOnlyTe,
  kTumorOnlyWgs,
  kTumorNormalWgs,
  kGermlineTagging,
  kGermlineMitochondria,
  kPanelOfNormals,
  kCustom,
};

/**
 * @brief Return the workflow executed by a concrete configuration profile.
 * @param profile Configuration profile to inspect.
 * @return Workflow associated with @p profile.
 * @throws xoos::error::Error if @p profile is not recognized.
 */
Workflow WorkflowForProfile(ConfigProfile profile);

/**
 * @brief Resolve a workflow whose name uniquely identifies one configuration profile.
 * @details Germline WGS is deliberately rejected because its workflow maps to both a duplex and a simplex profile.
 * @param workflow Workflow to resolve.
 * @return The unique configuration profile associated with @p workflow.
 * @throws xoos::error::Error if @p workflow is germline WGS or is not recognized.
 */
ConfigProfile ConfigProfileForWorkflow(Workflow workflow);

/**
 * @brief Select the germline WGS configuration profile for a concrete library type.
 * @details SBX-S selects the simplex profile. SBX-D and SBX-Fast share the duplex profile because they use the same
 * feature schema and training defaults.
 * @param library_type Resolved library type; auto-detect is not a concrete selection.
 * @return Simplex or duplex germline WGS profile.
 * @throws xoos::error::Error if @p library_type is auto-detect or is not recognized.
 */
ConfigProfile GermlineWgsConfigProfile(LibraryType library_type);

}  // namespace xoos::svc
