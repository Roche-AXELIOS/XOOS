#pragma once
#include <Eigen/Dense>

#include <xoos/types/vec.h>

#include "coverage/coverage.h"
#include "observations/observations.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {

static const std::string kMeanMapqColumnName = "MeanMapQ";

Observations GetAvgMapqsFromCoverage(const CoverageRecords& coverage);
Eigen::VectorXd GetAvgMeanMapqPerSegment(const vec<GenomicSegment>& segments, const Observations& mapqs);
vec<GenomicSegment>& AssignAvgMeanMapqPerSegment(vec<GenomicSegment>& segments, const Observations& mapqs);

}  // namespace xoos::cnc
