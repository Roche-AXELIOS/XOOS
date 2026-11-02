#pragma once
#include <xoos/types/int.h>

namespace xoos::svc {

/**
 * @brief Logical SVC execution workflow.
 * @details Workflow identifies algorithm routing, not a concrete configuration profile. The germline WGS workflow can
 * therefore use either the duplex or simplex profile selected from the resolved library type.
 */
enum class Workflow : u8 {
  kGermlineWgs,
  kTumorOnlyTe,
  kTumorOnlyWgs,
  kTumorNormalWgs,
  kGermlineTagging,
  kGermlineMitochondria,
  kPanelOfNormals,
  kCustom
};

/**
 * @brief Determine if the workflow is a germline workflow.
 * @param workflow Workflow enum value
 * @return true if the workflow is germline, false otherwise
 */
bool IsGermlineWorkflow(Workflow workflow);

/**
 * @brief Determine if the feature names for the workflow have sample context (e.g., "tumor_" or "normal_" prefix).
 * @param workflow Workflow enum value
 * @return true if the feature names have sample context, false otherwise
 */
bool FeatureNamesHaveSampleContext(Workflow workflow);

}  // namespace xoos::svc
