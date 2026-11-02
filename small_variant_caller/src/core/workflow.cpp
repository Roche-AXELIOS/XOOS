#include "workflow.h"

namespace xoos::svc {

bool IsGermlineWorkflow(const Workflow workflow) {
  using enum Workflow;
  return workflow == kGermlineWgs || workflow == kGermlineTagging;
}

bool FeatureNamesHaveSampleContext(const Workflow workflow) {
  using enum Workflow;
  return workflow == kTumorNormalWgs || workflow == kGermlineTagging;
}

}  // namespace xoos::svc
