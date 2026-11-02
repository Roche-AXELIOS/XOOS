#pragma once
#include <string>
#include <string_view>

#include <xoos/io/metadata-util.h>
#include <xoos/types/vec.h>

namespace xoos::cnc {
const std::string kColumnContig = "Contig";
const std::string kColumnContigAsComment = "#Contig";
const std::string kColumnStart = "Start";
const std::string kColumnEnd = "End";
const std::string kColumnMeanMapq = "MeanMapQ";
const std::string kColumnTotalCoverage = "TotalCoverage";
const std::string kColumnCounts = "Counts";
const std::string kColumnOnTarget = "OnTarget";
const std::string kColumnGCBias = "GCBias";
const std::string kColumnMappability = "Mappability";
const std::string kColumnLogRatio = "LogR";
const std::string kColumnRefAD = "Ref_AD";
const std::string kColumnAltAD = "Alt_AD";
const std::string kColumnBAF = "BAF";

// segment_likelihood.<prior>.tsv column names (Contig/Start/End reuse kColumnContig/kColumnStart/kColumnEnd above)
constexpr std::string_view kSegmentLikelihoodsOutPurity = "Purity";
constexpr std::string_view kSegmentLikelihoodsOutPloidy = "Ploidy";
constexpr std::string_view kSegmentLikelihoodsOutGermlineLogRLikelihood = "GermlineLogRLikelihood";
constexpr std::string_view kSegmentLikelihoodsOutGermlineTotalCopyNumber = "GermlineTotalCopyNumber";
constexpr std::string_view kSegmentLikelihoodsOutSomaticLogRLikelihood = "SomaticLogRLikelihood";
constexpr std::string_view kSegmentLikelihoodsOutSomaticBafLikelihood = "SomaticBafLikelihood";
constexpr std::string_view kSegmentLikelihoodsOutSomaticJointLikelihood = "SomaticJointLikelihood";
constexpr std::string_view kSegmentLikelihoodsOutSomaticTotalCopyNumber = "SomaticTotalCopyNumber";
constexpr std::string_view kSegmentLikelihoodsOutSomaticMinorCopyNumber = "SomaticMinorCopyNumber";
constexpr std::string_view kSegmentLikelihoodsOutClassification = "Classification";
constexpr std::string_view kSegmentLikelihoodsOutContributionLikelihood = "ContributionLikelihood";
constexpr std::string_view kSegmentClassificationGermline = "germline";
constexpr std::string_view kSegmentClassificationSomatic = "somatic";
}  // namespace xoos::cnc
