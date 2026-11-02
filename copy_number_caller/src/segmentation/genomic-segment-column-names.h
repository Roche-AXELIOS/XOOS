#pragma once
#include <string>

#include <xoos/types/vec.h>

namespace xoos::cnc {

static const std::string kGenomicSegColContig = "Contig";
static const std::string kGenomicSegColStart = "Start";
static const std::string kGenomicSegColEnd = "End";
static const std::string kGenomicSegColId = "Id";
static const std::string kGenomicSegColArrStart = "ArrStart";
static const std::string kGenomicSegColArrEnd = "ArrEnd";
static const std::string kGenomicSegColNumLogr = "NumLogRatio";
static const std::string kGenomicSegColNumSnp = "NumSnp";
static const std::string kGenomicSegColAvgMeanMapq = "AvgMeanMapq";
static const std::string kGenomicSegColTotalCopyNumber = "TotalCopyNumber";
static const std::string kGenomicSegColExpectedTotalCopyNumber = "ExpectedTotalCopyNumber";
static const std::string kGenomicSegColMajorCopyNumber = "MajorCopyNumber";
static const std::string kGenomicSegColMinorCopyNumber = "MinorCopyNumber";
static const std::string kGenomicSegColLogrLikelihood = "LogRatioLikelihood";
static const std::string kGenomicSegColBafLikelihood = "BAlleleFrequencyLikelihood";
static const std::string kGenomicSegColJointLikelihood = "JointLikelihood";
static const std::string kGenomicSegColPurity = "Purity";
static const std::string kGenomicSegColPloidy = "Ploidy";
static const std::string kGenomicSegColSex = "Sex";
static const std::string kGenomicSegColFilter = "Filter";
static const std::string kGenomicSegColMBaf = "MirroredBaf";
static const std::string kGenomicSegColInAllosome = "IsInAllosome";
static const std::string kGenomicSegColInPseudoAutosomalRegion = "IsInPseudoAutosomalRegion";
static const std::string kGenomicSegColIsGermlineCnv = "IsGermlineCnv";
static const std::string kGenomicSegColMeanDh = "MeanDh";
static const std::string kGenomicSegColMeanLogr = "MeanLogRatio";

static const vec<std::string> kGenomicSegColsAll{
    kGenomicSegColId,
    kGenomicSegColContig,
    kGenomicSegColStart,
    kGenomicSegColEnd,
    kGenomicSegColArrStart,
    kGenomicSegColArrEnd,
    kGenomicSegColNumLogr,
    kGenomicSegColNumSnp,
    kGenomicSegColAvgMeanMapq,
    kGenomicSegColTotalCopyNumber,
    kGenomicSegColExpectedTotalCopyNumber,
    kGenomicSegColMajorCopyNumber,
    kGenomicSegColMinorCopyNumber,
    kGenomicSegColLogrLikelihood,
    kGenomicSegColBafLikelihood,
    kGenomicSegColJointLikelihood,
    kGenomicSegColPurity,
    kGenomicSegColPloidy,
    kGenomicSegColSex,
    kGenomicSegColFilter,
    kGenomicSegColInAllosome,
    kGenomicSegColInPseudoAutosomalRegion,
    kGenomicSegColMeanDh,
    kGenomicSegColMeanLogr,
    kGenomicSegColIsGermlineCnv,
};

static const vec<std::string> kGenomicSegColsRequiredSeed{kGenomicSegColId,
                                                          kGenomicSegColContig,
                                                          kGenomicSegColStart,
                                                          kGenomicSegColEnd,
                                                          kGenomicSegColInAllosome,
                                                          kGenomicSegColInPseudoAutosomalRegion};

static const vec<std::string> kGenomicSegColsRequiredLogRSegments{kGenomicSegColId,
                                                                  kGenomicSegColContig,
                                                                  kGenomicSegColStart,
                                                                  kGenomicSegColEnd,
                                                                  kGenomicSegColNumLogr,
                                                                  kGenomicSegColInAllosome,
                                                                  kGenomicSegColInPseudoAutosomalRegion,
                                                                  kGenomicSegColMeanLogr};

static const vec<std::string> kGenomicSegColsRequiredBafSegments{kGenomicSegColId,
                                                                 kGenomicSegColContig,
                                                                 kGenomicSegColStart,
                                                                 kGenomicSegColEnd,
                                                                 kGenomicSegColNumLogr,
                                                                 kGenomicSegColNumSnp,
                                                                 kGenomicSegColInAllosome,
                                                                 kGenomicSegColInPseudoAutosomalRegion,
                                                                 kGenomicSegColMeanDh,
                                                                 kGenomicSegColMeanLogr};

static const vec<std::string> kGenomicSegColsRequiredGermline{kGenomicSegColId,
                                                              kGenomicSegColContig,
                                                              kGenomicSegColStart,
                                                              kGenomicSegColEnd,
                                                              kGenomicSegColNumLogr,
                                                              kGenomicSegColTotalCopyNumber,
                                                              kGenomicSegColLogrLikelihood,
                                                              kGenomicSegColAvgMeanMapq,
                                                              kGenomicSegColExpectedTotalCopyNumber,
                                                              kGenomicSegColSex,
                                                              kGenomicSegColFilter,
                                                              kGenomicSegColInAllosome,
                                                              kGenomicSegColInPseudoAutosomalRegion,
                                                              kGenomicSegColMeanLogr};

static const vec<std::string> kGenomicSegColsRequiredSomaticAlleleSpecific{
    kGenomicSegColId,
    kGenomicSegColContig,
    kGenomicSegColStart,
    kGenomicSegColEnd,
    kGenomicSegColNumLogr,
    kGenomicSegColNumSnp,
    kGenomicSegColTotalCopyNumber,
    kGenomicSegColLogrLikelihood,
    kGenomicSegColBafLikelihood,
    kGenomicSegColJointLikelihood,
    kGenomicSegColExpectedTotalCopyNumber,
    kGenomicSegColMajorCopyNumber,
    kGenomicSegColMinorCopyNumber,
    kGenomicSegColPurity,
    kGenomicSegColPloidy,
    kGenomicSegColSex,
    kGenomicSegColMBaf,
    kGenomicSegColInAllosome,
    kGenomicSegColInPseudoAutosomalRegion,
    kGenomicSegColMeanDh,
    kGenomicSegColMeanLogr,
};

static const vec<std::string> kGenomicSegColsRequiredSomaticNoAlleleSpecific{kGenomicSegColId,
                                                                             kGenomicSegColContig,
                                                                             kGenomicSegColStart,
                                                                             kGenomicSegColEnd,
                                                                             kGenomicSegColNumLogr,
                                                                             kGenomicSegColTotalCopyNumber,
                                                                             kGenomicSegColLogrLikelihood,
                                                                             kGenomicSegColExpectedTotalCopyNumber,
                                                                             kGenomicSegColPurity,
                                                                             kGenomicSegColPloidy,
                                                                             kGenomicSegColSex,
                                                                             kGenomicSegColInAllosome,
                                                                             kGenomicSegColInPseudoAutosomalRegion,
                                                                             kGenomicSegColMeanLogr};

// Tumor-only column set mirrors the somatic allele-specific set above with an IsGermlineCnv column that records the
// post-hoc germline-vs-somatic classification. MeanLogRatio remains last so IGV recognizes it as the numeric value.
static const vec<std::string> kGenomicSegColsRequiredSomaticTumorOnlyAlleleSpecific{
    kGenomicSegColId,
    kGenomicSegColContig,
    kGenomicSegColStart,
    kGenomicSegColEnd,
    kGenomicSegColNumLogr,
    kGenomicSegColNumSnp,
    kGenomicSegColTotalCopyNumber,
    kGenomicSegColLogrLikelihood,
    kGenomicSegColBafLikelihood,
    kGenomicSegColJointLikelihood,
    kGenomicSegColExpectedTotalCopyNumber,
    kGenomicSegColMajorCopyNumber,
    kGenomicSegColMinorCopyNumber,
    kGenomicSegColPurity,
    kGenomicSegColPloidy,
    kGenomicSegColSex,
    kGenomicSegColMBaf,
    kGenomicSegColInAllosome,
    kGenomicSegColInPseudoAutosomalRegion,
    kGenomicSegColMeanDh,
    kGenomicSegColIsGermlineCnv,
    kGenomicSegColMeanLogr,
};

static const vec<std::string> kGenomicSegColsOptional{kGenomicSegColArrStart,
                                                      kGenomicSegColArrEnd,
                                                      kGenomicSegColNumLogr,
                                                      kGenomicSegColNumSnp,
                                                      kGenomicSegColAvgMeanMapq,
                                                      kGenomicSegColTotalCopyNumber,
                                                      kGenomicSegColExpectedTotalCopyNumber,
                                                      kGenomicSegColMajorCopyNumber,
                                                      kGenomicSegColMinorCopyNumber,
                                                      kGenomicSegColLogrLikelihood,
                                                      kGenomicSegColBafLikelihood,
                                                      kGenomicSegColJointLikelihood,
                                                      kGenomicSegColPurity,
                                                      kGenomicSegColPloidy,
                                                      kGenomicSegColSex,
                                                      kGenomicSegColFilter,
                                                      kGenomicSegColMeanDh,
                                                      kGenomicSegColMeanLogr,
                                                      kGenomicSegColIsGermlineCnv};

}  // namespace xoos::cnc
