#pragma once

namespace xoos::svc {

enum class Workflow {
  kGermline,
  kGermlineMultiSample,
  kTumorOnlyTe,
  kTumorNormalWgs,
  kGermlineTagging,
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
