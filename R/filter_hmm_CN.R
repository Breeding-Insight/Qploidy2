# Suppress R CMD check notes for dplyr/ggplot2 NSE column references
if (getRversion() >= "2.15.1") utils::globalVariables(c(
  "CN_call", "Sample", "Chr", "n_cn", "dom_cn", "sample_cn", "cell_type",
  "has_mixed", "has_aneuploid", "sample_type"
))


#'
#' Sets \code{CN_call}, \code{post_max}, \code{CN_reliability}, and all
#' \code{post_CN*} columns to \code{NA} for windows that fail any of the
#' supplied quality thresholds.  Masked calls are treated as missing by
#' downstream plotting and analysis functions.  The \code{by_marker} data
#' frame is also updated when present.
#'
#' @param hmm_CN An object of class \code{hmm_CN}.
#' @param min_CN_reliability Numeric \[0,1\]. Minimum joint reliability score
#'   (HMM posterior × BAF support). Default \code{0.15}. Particularly useful
#'   for suppressing high-confidence z-only calls on all-homozygous windows.
#' @param min_post_max Numeric \[0,1\]. Minimum raw HMM posterior probability.
#'   Use when you need a floor on the HMM's own confidence regardless of BAF.
#'   Default \code{NULL} (no filter).
#' @param min_w_baf Numeric \[0,1\]. Minimum BAF weight. Set to a positive
#'   value to exclude windows where BAF evidence is absent or negligible.
#'   Default \code{NULL} (no filter).
#' @param min_n_snps Integer. Minimum number of SNPs per window.  Windows
#'   with very few markers produce noisy z-score and BAF estimates.
#'   Default \code{NULL} (no filter).
#' @param min_n_het Integer. Minimum number of heterozygous markers per window.
#'   A finer control than \code{min_w_baf}: excludes windows where the BAF
#'   weight is non-zero but only a handful of hets drive it.
#'   Default \code{NULL} (no filter).
#' @param min_window_size Numeric. Minimum genomic window span in base pairs.
#'   Very short windows can reflect assembly gaps or edge artefacts.
#'   Default \code{NULL} (no filter).
#' @param max_CN_call Integer. Maximum CN call to retain.  Extreme values
#'   (e.g., > 8 in most polyploid contexts) are often noise rather than
#'   genuine amplifications.  Default \code{NULL} (no filter).
#'
#' @return A modified object of class \code{hmm_CN} with \code{NA} substituted
#'   in \code{by_window} (and \code{by_marker} when present) for all windows
#'   that fail at least one threshold.
#'
#' @examples
#' \dontrun{
#' # Mask any window without BAF support and with low overall reliability
#' filtered <- filter_hmm_CN(hmm_CN,
#'                           min_CN_reliability = 0.3,
#'                           min_w_baf          = 0.05,
#'                           min_n_snps         = 20L)
#' }
#'
#' @export
filter_hmm_CN <- function(hmm_CN,
                          min_CN_reliability = 0.15,
                          min_post_max       = NULL,
                          min_w_baf          = NULL,
                          min_n_snps         = NULL,
                          min_n_het          = NULL,
                          min_window_size    = NULL,
                          max_CN_call        = NULL) {

  stopifnot(inherits(hmm_CN, "hmm_CN"))
  bw <- hmm_CN$by_window

  # ---- build failure mask (TRUE = window fails ≥ 1 criterion) ----
  fail <- rep(FALSE, nrow(bw))

  if (!is.null(min_CN_reliability) && "CN_reliability" %in% names(bw))
    fail <- fail | is.na(bw$CN_reliability) | bw$CN_reliability < min_CN_reliability

  if (!is.null(min_post_max) && "post_max" %in% names(bw))
    fail <- fail | is.na(bw$post_max) | bw$post_max < min_post_max

  if (!is.null(min_w_baf) && "w_baf" %in% names(bw))
    fail <- fail | is.na(bw$w_baf) | bw$w_baf < min_w_baf

  if (!is.null(min_n_snps) && "n_snps" %in% names(bw))
    fail <- fail | is.na(bw$n_snps) | bw$n_snps < as.integer(min_n_snps)

  if (!is.null(min_n_het) && "n_het" %in% names(bw))
    fail <- fail | is.na(bw$n_het) | bw$n_het < as.integer(min_n_het)

  if (!is.null(min_window_size)) {
    win_size <- as.numeric(bw$End) - as.numeric(bw$Start)
    fail <- fail | is.na(win_size) | win_size < min_window_size
  }

  if (!is.null(max_CN_call) && "CN_call" %in% names(bw))
    fail <- fail | (!is.na(bw$CN_call) & bw$CN_call > as.integer(max_CN_call))

  if (!any(fail)) {
    message("No windows failed the thresholds; hmm_CN returned unchanged.")
    return(hmm_CN)
  }

  message(sprintf("%d of %d windows (%.1f%%) masked.",
                  sum(fail), nrow(bw), 100 * mean(fail)))

  # ---- mask by_window ----
  post_cols  <- grep("^post_CN", names(bw), value = TRUE)
  cols_to_na <- intersect(c("CN_call", "post_max", "CN_reliability", post_cols), names(bw))
  bw[fail, cols_to_na] <- NA
  hmm_CN$by_window <- bw

  # ---- propagate to by_marker ----
  if (!is.null(hmm_CN$by_marker) && ".__w__" %in% names(hmm_CN$by_marker)) {
    bm <- hmm_CN$by_marker
    # composite key: Sample + Chr + WindowID avoids collisions across samples
    bw_key  <- paste(bw$Sample,       bw$Chr, bw$WindowID,    sep = "|||")
    bm_key  <- paste(bm$SampleName,   bm$Chr, bm[[".__w__"]], sep = "|||")
    bm_fail <- bm_key %in% bw_key[fail]
    marker_cols <- intersect(c("CN_call", "post_max", "CN_reliability", "w_baf"), names(bm))
    bm[bm_fail, marker_cols] <- NA
    hmm_CN$by_marker <- bm
  }

  hmm_CN
}


#' Count and classify samples by ploidy category
#'
#' Classifies every sample in an \code{hmm_CN} object as euploid, aneuploid,
#' or segmental aneuploidy, prints a summary table, and returns a list of
#' sample-ID vectors for each (sample-level CN, category) combination.
#'
#' Definitions used:
#' \describe{
#'   \item{euploid}{All chromosomes carry the same CN as the sample-wide mode CN.}
#'   \item{aneuploid}{At least one chromosome has a different but uniform CN
#'     (whole-chromosome aneuploidy); no within-chromosome mixed CN.}
#'   \item{segmental}{At least one chromosome has windows with mixed CN
#'     (segmental/partial aneuploidy).}
#' }
#'
#' @param hmm_CN  An object of class \code{hmm_CN}.
#' @param chromosomes  Optional character vector. If supplied, only these
#'   chromosomes are used for classification.
#'
#' @return Invisibly returns a named list.  The top-level names are
#'   \code{"euploid"}, \code{"aneuploid"}, and \code{"segmental"}, each
#'   containing the sample IDs of that type.  A nested element \code{$by_cn}
#'   further splits those vectors by sample-level CN
#'   (e.g. \code{result$by_cn[["2x"]]$euploid}).
#'
#' @importFrom dplyr group_by summarise mutate left_join n_distinct
#' @export
count_types <- function(hmm_CN, chromosomes = NULL) {

  stopifnot(inherits(hmm_CN, "hmm_CN"))
  bw <- hmm_CN$by_window

  if (!is.null(chromosomes))
    bw <- bw[bw$Chr %in% chromosomes, ]
  if (nrow(bw) == 0) stop("No data remaining after chromosome filter.")

  # NA-safe mode: ignores NA values; returns NA_integer_ when all values are NA
  .mode <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) return(NA_integer_)
    tt <- table(x); as.integer(names(tt)[which.max(tt)])
  }

  sample_cn_df <- bw %>%
    group_by(Sample) %>%
    summarise(sample_cn = .mode(CN_call), .groups = "drop") %>%
    filter(!is.na(sample_cn))

  cell_df <- bw %>%
    group_by(Sample, Chr) %>%
    summarise(n_cn = n_distinct(CN_call, na.rm = TRUE), dom_cn = .mode(CN_call), .groups = "drop") %>%
    filter(!is.na(dom_cn)) %>%
    left_join(sample_cn_df, by = "Sample") %>%
    mutate(cell_type = ifelse(n_cn > 1, "mixed",
                       ifelse(dom_cn == sample_cn, "euploid_chr", "aneuploid_chr")))

  sample_class <- cell_df %>%
    group_by(Sample, sample_cn) %>%
    summarise(
      has_mixed     = any(cell_type == "mixed"),
      has_aneuploid = any(cell_type == "aneuploid_chr"),
      .groups = "drop"
    ) %>%
    mutate(
      sample_type = ifelse(has_mixed, "segmental",
                    ifelse(has_aneuploid, "aneuploid", "euploid"))
    )

  # ---- build and return classified object ----
  all_cns <- sort(unique(sample_class$sample_cn))
  result <- list(
    euploid   = sample_class$Sample[sample_class$sample_type == "euploid"],
    aneuploid = sample_class$Sample[sample_class$sample_type == "aneuploid"],
    segmental = sample_class$Sample[sample_class$sample_type == "segmental"],
    by_cn     = setNames(lapply(all_cns, function(cn) {
      sub <- sample_class[sample_class$sample_cn == cn, ]
      list(
        euploid   = sub$Sample[sub$sample_type == "euploid"],
        aneuploid = sub$Sample[sub$sample_type == "aneuploid"],
        segmental = sub$Sample[sub$sample_type == "segmental"]
      )
    }), ifelse(all_cns == 1L, "1x or inbred", paste0(all_cns, "x"))),
    .summary  = sample_class   # stored for print.count_types; not intended for direct use
  )
  class(result) <- "count_types"
  result
}


#' Print method for count_types objects
#'
#' Displays the sample classification summary table produced by
#' \code{\link{count_types}}.
#'
#' @param x  A \code{count_types} object.
#' @param ... Ignored.
#'
#' @method print count_types
#' @export
print.count_types <- function(x, ...) {
  sample_class <- x$.summary
  all_cns   <- sort(unique(sample_class$sample_cn))
  all_types <- c("euploid", "aneuploid", "segmental")
  totals    <- c(euploid = 0L, aneuploid = 0L, segmental = 0L)

  cat("\nSample classification summary\n")
  cat(strrep("=", 52), "\n", sep = "")
  for (cn in all_cns) {
    sub <- sample_class[sample_class$sample_cn == cn, ]
    cn_label <- if (cn == 1L) "1x or inbred" else paste0(cn, "x")
    cat(sprintf("  %s  (n = %d)\n", cn_label, nrow(sub)))
    for (tp in all_types) {
      n <- sum(sub$sample_type == tp)
      if (n > 0) {
        cat(sprintf("    %-22s  %d\n", paste0(cn_label, " - ", tp), n))
        totals[tp] <- totals[tp] + n
      }
    }
  }
  cat(strrep("-", 52), "\n", sep = "")
  cat(sprintf("  Total euploid              %d\n", totals["euploid"]))
  cat(sprintf("  Total aneuploid            %d\n", totals["aneuploid"]))
  cat(sprintf("  Total segmental aneuploidy %d\n", totals["segmental"]))
  cat(sprintf("  Grand total                %d\n", nrow(sample_class)))
  cat(strrep("=", 52), "\n\n", sep = "")
  invisible(x)
}

