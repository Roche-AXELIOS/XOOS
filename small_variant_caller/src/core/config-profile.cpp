#include "config-profile.h"

#include <xoos/error/error.h>

namespace xoos::svc {

Workflow WorkflowForProfile(const ConfigProfile profile) {
  using enum ConfigProfile;
  switch (profile) {
    case kGermlineWgsDuplex:
    case kGermlineWgsSimplex:
      return Workflow::kGermlineWgs;
    case kTumorOnlyTe:
      return Workflow::kTumorOnlyTe;
    case kTumorOnlyWgs:
      return Workflow::kTumorOnlyWgs;
    case kTumorNormalWgs:
      return Workflow::kTumorNormalWgs;
    case kGermlineTagging:
      return Workflow::kGermlineTagging;
    case kGermlineMitochondria:
      return Workflow::kGermlineMitochondria;
    case kPanelOfNormals:
      return Workflow::kPanelOfNormals;
    case kCustom:
      return Workflow::kCustom;
    default:
      throw error::Error("Unknown config profile");
  }
}

ConfigProfile ConfigProfileForWorkflow(const Workflow workflow) {
  using enum Workflow;
  switch (workflow) {
    case kGermlineWgs:
      // Unlike other workflows, germline WGS has multiple concrete profiles; workflow alone cannot select one safely.
      throw error::Error("Workflow 'germline-wgs' requires a resolved library type to select its config profile");
    case kTumorOnlyTe:
      return ConfigProfile::kTumorOnlyTe;
    case kTumorOnlyWgs:
      return ConfigProfile::kTumorOnlyWgs;
    case kTumorNormalWgs:
      return ConfigProfile::kTumorNormalWgs;
    case kGermlineTagging:
      return ConfigProfile::kGermlineTagging;
    case kGermlineMitochondria:
      return ConfigProfile::kGermlineMitochondria;
    case kPanelOfNormals:
      return ConfigProfile::kPanelOfNormals;
    case kCustom:
      return ConfigProfile::kCustom;
    default:
      throw error::Error("Unknown workflow");
  }
}

ConfigProfile GermlineWgsConfigProfile(const LibraryType library_type) {
  using enum LibraryType;
  switch (library_type) {
    case kSbxS:
      return ConfigProfile::kGermlineWgsSimplex;
    case kSbxD:
    case kSbxFast:
      // SBX-D and SBX-Fast intentionally share the duplex feature configuration while retaining distinct model tokens.
      return ConfigProfile::kGermlineWgsDuplex;
    case kAutoDetect:
      throw error::Error("Library type must be resolved before selecting a germline-wgs config profile");
    default:
      throw error::Error("Unknown library type");
  }
}

}  // namespace xoos::svc
