#pragma once

#include <filesystem>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::sex_predict {

namespace fs = std::filesystem;

// Minimum base quality for a pileup base to contribute to a site's allele
// counts. Low-quality bases are the main source of spurious minor alleles.
constexpr u8 kDefaultMinBaseQuality = 20;

// Minimum absolute read depth for a site to be evaluated for heterozygosity.
// Below this, a genuine heterozygous site can be missed by chance alone, which
// would deflate the measured rate rather than merely add noise.
constexpr u32 kMinSiteDepth = 20;

// A site counts as heterozygous when its second-most-common allele holds a
// fraction of total depth within [kMinHetAlleleFraction, kMaxHetAlleleFraction]
// AND both of the top two alleles have at least kMinAlleleSupportingReads
// reads. Requiring two independently supported alleles keeps sequencing error
// -- which rarely reaches 25% with several Q20 bases -- out of the numerator.
constexpr f64 kMinHetAlleleFraction = 0.25;
constexpr f64 kMaxHetAlleleFraction = 0.75;
constexpr u32 kMinAlleleSupportingReads = 3;

// Sites are additionally restricted to a band around the *observed median site
// depth of the same window set*. The upper bound rejects collapsed repeats and
// segmental duplications, which are the dominant source of false heterozygosity
// on a one-copy chromosome. The band must be relative, not absolute: a one-copy
// chrX sits at roughly half the autosome depth, so a fixed band would filter
// the two groups unequally and bias the very ratio being measured.
constexpr f64 kMinRelativeSiteDepth = 0.5;
constexpr f64 kMaxRelativeSiteDepth = 2.0;

// Upper bound of the per-site depth histogram used to derive the median. Sites
// deeper than this are pooled into the final bucket; they are far above any
// plausible WGS median and are excluded by the relative band regardless.
constexpr u32 kDepthHistogramSize = 4096;

/// A 0-based, half-open window on a single contig, the XOOS convention.
struct ScanWindow {
  s32 tid;
  s64 start;
  s64 end;
};

/// Heterozygosity tally over a set of windows.
struct HetCounts {
  /// Sites that passed every depth filter and were tested for heterozygosity.
  u64 evaluated_sites = 0;
  /// Of those, the sites that met the heterozygous-allele criterion.
  u64 het_sites = 0;
  /// Median per-site depth observed across the windows, before the relative
  /// band was applied. Reported for logging and threshold calibration.
  u32 median_site_depth = 0;

  /// Heterozygous fraction, or 0.0 when nothing was evaluated. Callers must
  /// check `evaluated_sites` against a minimum before trusting this.
  f64 Rate() const;
};

/// Count heterozygous sites across `windows` using an exact pileup.
///
/// The depth band is derived from these windows' own median site depth, so each
/// call self-calibrates: pass chrX windows and autosome windows as *separate*
/// calls so that each group is filtered against its own depth distribution.
///
/// Windows are distributed across `num_threads` taskflow workers, each with its
/// own reader (a single htsFile is not safe for concurrent use).
HetCounts CountHeterozygousSites(
    const fs::path& bam_file_path, const vec<ScanWindow>& windows, u8 min_mapq, u8 min_base_quality, u32 num_threads);

}  // namespace xoos::sex_predict
