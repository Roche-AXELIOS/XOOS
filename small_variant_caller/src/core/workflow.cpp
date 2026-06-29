#include "workflow.h"

namespace xoos::svc {

bool IsGermlineWorkflow(const Workflow workflow) {
  using enum Workflow;
  return (workflow == kGermline) || (workflow == kGermlineMultiSample) || (workflow == kGermlineTagging);
}

bool FeatureNamesHaveSampleContext(const Workflow workflow) {
  using enum Workflow;
  return workflow == kTumorNormalWgs || workflow == kGermlineTagging;
}

}  // namespace xoos::svc
