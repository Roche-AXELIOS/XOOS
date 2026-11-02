#!/usr/bin/env Rscript

# This script generates seed segments

# The libraries to load. This should be as minimal as possible.
loaded_libs <- c()

# Required library packages
required_libs <-
  c(loaded_libs, "dplyr", "GenomicRanges")

# If any packages aren't installed, list them and stop
missing_libs <-
  required_libs[!required_libs %in% rownames(installed.packages())]
if (length(missing_libs) > 0)
  stop(paste(c("Missing required libraries:", missing_libs), collapse = " "))

# Load the necessary packages (only do this in scripts, not library code)
for (lib in loaded_libs)
  suppressPackageStartupMessages(library(lib, character.only = TRUE))


#
# This function generates a set of seed segments to use for segmentation.
#

#' Gets seed segments
#'
#' @param genome String indicating the genome build to generate seed segments 
#'   for
#' @param exclusion_buffer Integer indicating the size of the buffer zone to
#'   exclude upstream and downstream of the centromere gap
get_seed_segs_df <- function(genome, exclusion_buffer = 1e5) {

  if (!genome %in% c("hg19", "hg38")) {
    stop("genome must be either hg19 or hg38")
  }

  if (genome == "hg19") {
    base_seed_segs_df <- rCGH::hg19
  } else if (genome == "hg38") {
    base_seed_segs_df <- rCGH::hg38
  }

  per_chr_seed_segs_df <-
    base_seed_segs_df |>
    tibble::as_tibble() |>
    dplyr::mutate(
      chrom = paste0("chr", .data[["chrom"]])
    ) |>
    dplyr::mutate(
      chrom =
        dplyr::case_when(
          .data[["chrom"]] == "chr23" ~ "chrX",
          .data[["chrom"]] == "chr24" ~ "chrY",
          TRUE ~ .data[["chrom"]]
        )
    )

  # rCGH::hg38 provides correct hg38 chromosome lengths but ships the hg19
  # centromere coordinates, so its centromere positions must not be trusted for
  # hg38. Override them with authoritative hg38 boundaries derived from the UCSC
  # cytoBand acen bands. hg19's rCGH centromere coordinates are correct and are
  # used as-is.
  if (genome == "hg38") {
    per_chr_seed_segs_df <-
      per_chr_seed_segs_df |>
      dplyr::select(-"centromerStart", -"centromerEnd") |>
      dplyr::inner_join(get_hg38_centromere_df(), by = "chrom")
  }

  per_chr_seed_segs_df <-
    per_chr_seed_segs_df |>
    # Exclude buffer zone upstream and downstream of the centromere gap
    dplyr::mutate(
      p_arm_start = 1,
      p_arm_end = .data[["centromerStart"]] - exclusion_buffer,
      q_arm_start = .data[["centromerEnd"]] + exclusion_buffer,
      q_arm_end = .data[["length"]]
    ) 
  num_invalid_p_arms_ends <-
    per_chr_seed_segs_df |>
    dplyr::filter(p_arm_end < 1) |>
    nrow()

  if (num_invalid_p_arms_ends > 0) {
    stop("Invalid p-arm end(s) found. Try reducing the exclusion buffer")
  }

  # Exclude the centromere region
  seed_segs_df <-
    per_chr_seed_segs_df |> 
    dplyr::select("chrom", start = "p_arm_start", end = "p_arm_end") |>
    dplyr::bind_rows(
      per_chr_seed_segs_df |> 
        dplyr::select("chrom", start = "q_arm_start", end = "q_arm_end")
    ) |>
    dplyr::mutate(
      chrom = factor(chrom, levels = paste0("chr", c(1:22, "X", "Y")))
    ) |>
    dplyr::arrange(chrom, start)

  autosomal_seed_segs_df <-
    seed_segs_df |>
    dplyr::filter(!chrom %in% c("chrX", "chrY")) |>
    dplyr::mutate(is_in_par = FALSE)

  # Generate PAR region specific segments
  allosome_seed_segs_gr <-
    seed_segs_df |>
    dplyr::filter(chrom %in% c("chrX", "chrY")) |>
    GenomicRanges::makeGRangesFromDataFrame(
      seqnames.field = "chrom",
      start.field = "start",
      end.field = "end",
    )

  par_gr <- get_par_gr(genome) 
  non_par_gr <- GenomicRanges::setdiff(allosome_seed_segs_gr, par_gr)
  GenomicRanges::mcols(non_par_gr)$is_in_par <- FALSE

  allosome_seed_segs_df <-
    c(par_gr, non_par_gr) |>
    sort() |>
    as.data.frame() |>
    dplyr::rename(chrom = "seqnames") |>
    dplyr::select(chrom, start, end, is_in_par) |>
    dplyr::group_by(chrom) |>
    # Take away the first (before PAR1) and last segments (before PAR2) to 
    # avoid telomeric regions
    dplyr::slice(2:(dplyr::n() - 1))

  seed_segs_df <-
    dplyr::bind_rows(autosomal_seed_segs_df, allosome_seed_segs_df) |>
    dplyr::mutate(
      id = "seed_segment",
      is_in_allosome = chrom %in% c("chrX", "chrY"),
    ) |>
    dplyr::mutate(
      is_in_allosome = ifelse(is_in_allosome, "true", "false"),
      is_in_par = ifelse(is_in_par, "true", "false")
    ) |>
    dplyr::select(
      Id = "id", 
      Contig = "chrom", 
      Start = "start", 
      End = "end",
      IsInAllosome = is_in_allosome,
      IsInPseudoAutosomalRegion = is_in_par
    ) 

  return(seed_segs_df)
}


#' Function to get authoritative hg38 centromere boundaries
#'
#' rCGH::hg38 ships hg19 centromere coordinates, so hg38 centromere positions
#' are taken instead from the UCSC cytoBand `acen` bands
#' (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz).
#' For each chromosome the centromere spans from the minimum `acen` band start
#' to the maximum `acen` band end. These boundaries replace the `centromerStart`
#' and `centromerEnd` columns rCGH would otherwise supply; the `exclusion_buffer`
#' is applied later by `get_seed_segs_df()` when creating the p/q arm segments.
#'
#' @return A tibble with columns chrom, centromerStart, centromerEnd
get_hg38_centromere_df <- function() {
  tibble::tribble(
    ~chrom,  ~centromerStart, ~centromerEnd,
    "chr1",  121700000,       125100000,
    "chr2",  91800000,        96000000,
    "chr3",  87800000,        94000000,
    "chr4",  48200000,        51800000,
    "chr5",  46100000,        51400000,
    "chr6",  58500000,        62600000,
    "chr7",  58100000,        62100000,
    "chr8",  43200000,        47200000,
    "chr9",  42200000,        45500000,
    "chr10", 38000000,        41600000,
    "chr11", 51000000,        55800000,
    "chr12", 33200000,        37800000,
    "chr13", 16500000,        18900000,
    "chr14", 16100000,        18200000,
    "chr15", 17500000,        20500000,
    "chr16", 35300000,        38400000,
    "chr17", 22700000,        27400000,
    "chr18", 15400000,        21500000,
    "chr19", 24200000,        28100000,
    "chr20", 25700000,        30400000,
    "chr21", 10900000,        13000000,
    "chr22", 13700000,        17400000,
    "chrX",  58100000,        63800000,
    "chrY",  10300000,        10600000
  )
}


#' Function to get PAR regions as GRanges
#'
#' @param genome String indicating the genome build to generate PAR regions 
#'   for
get_par_gr <- function(genome) {
  if (genome == "hg38") {
    par_df <-
      tibble::tribble(
        ~name,  ~chrom, ~start,    ~end,      
        "PAR1", "chrX", 10001,     2781479,   
        "PAR1", "chrY", 10001,     2781479,   
        "PAR2", "chrX", 155701383, 156030895, 
        "PAR2", "chrY", 56887903,  57217415,
      ) |>
      dplyr::mutate(is_in_par = TRUE)
  } else if (genome == "hg19") {
    par_df <-
      tibble::tribble(
        ~name,  ~chrom, ~start,    ~end,      
        "PAR1", "chrX", 60001,     2699520,
        "PAR1", "chrY", 10001,     2649520,
        "PAR2", "chrX", 154931044, 155260560,
        "PAR2", "chrY", 59034050,	 59363566
      ) |>
      dplyr::mutate(is_in_par = TRUE)
  } else {
    stop("genome must be either hg19 or hg38")
  }

  par_gr <-
    GenomicRanges::makeGRangesFromDataFrame(
      par_df,
      seqnames.field = "chrom",
      start.field = "start",
      end.field = "end",
      keep.extra.columns = TRUE
    )
  return(par_gr)
}

message("Getting seed segments")
hg19_seed_segs_df <- get_seed_segs_df("hg19")
hg38_seed_segs_df <- get_seed_segs_df("hg38")

message("Writing seed segments")
readr::write_tsv(hg19_seed_segs_df, "hg19_seed.seg")
readr::write_tsv(hg38_seed_segs_df, "hg38_seed.seg")
