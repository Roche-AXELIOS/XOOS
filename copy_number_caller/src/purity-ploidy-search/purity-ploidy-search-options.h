#pragma once

#include <xoos/types/float.h>
#include <xoos/types/fs.h>

namespace xoos::cnc {

// Maximum LogR value for a segment to be considered
const f64 kSegMaxLogR = 0.9;
const size_t kSegMinNumLogRs = 11;
const size_t kSegMinNumSnps = 11;
const f64 kBafLikelihoodWeight = 1.0;
const f64 kTumorOnlyWgsBafLikelihoodWeight = 1.1;
// When purity is too low for the mBAF-peak WGD test, a genome that is overwhelmingly loss-of-heterozygosity is
// itself strong evidence of whole-genome doubling. If at least this fraction of the assessable autosomal genome
// (by length) is in loss of heterozygosity, the diploid-bias solution is rejected in favour of a uniform prior.
const f64 kMinLohFractionForWgd = 0.35;

struct PurityPloidySearchOptions {
  f64 seg_max_logr = kSegMaxLogR;
  size_t seg_min_num_logrs = kSegMinNumLogRs;
  size_t seg_min_num_snps = kSegMinNumSnps;
  // Minimum fraction of the assessable autosomal genome (by length) that must be loss-of-heterozygosity (minor copy
  // number 0) for the low-purity fallback to reject the diploid-bias solution as whole-genome doubling.
  f64 min_loh_fraction_for_wgd = kMinLohFractionForWgd;
  // Weight applied to the BAF term of the joint likelihood (joint = logr + weight * baf) during the WGD re-run only.
  // The initial diploid-bias pass always uses 1.0. Up-weighting BAF is justified because the WGD QC step that triggers
  // the re-run rejects the initial solution purely on BAF evidence.
  f64 baf_likelihood_weight = kBafLikelihoodWeight;
  // When true, classify segments as germline vs somatic during the purity/ploidy search so that
  // rare germline CNVs do not inflate the purity estimate (tumor-only WGS mode).
  bool tumor_only = false;
  // When true, emit the per-segment, per-grid-cell likelihood diagnostic (one TSV per search pass).
  bool write_segment_likelihood = false;
};
}  // namespace xoos::cnc
