if (getRversion() >= "2.15.1") utils::globalVariables(c("CN_mode",
  "V1", "V2", "V4", "V5", "V8", ".__w__",
  "post_CN1", "post_CN2", "post_CN3", "post_CN4", "post_CN5", "post_CN6"))

#' Summarize copy number mode and posterior probability
#'
#' Summarizes the mode of copy number calls and the mean maximum posterior probability per sample,
#' chromosome, or chromosome arm. Handles input as a data.frame or hmm_CN object.
#'
#' @param df Data.frame or hmm_CN object containing windowed CN calls and posterior probabilities.
#' @param level Character. Summarization level: 'sample', 'chromosome', or 'chromosome-arm'.
#' @param centromeres Optional. Named numeric vector or data.frame with centromere positions (required for chromosome-arm).
#' @param cn_col Character. Column name for CN calls (default: 'CN_call').
#' @param post_col Character. Column name for maximum posterior probability (default: 'post_max').
#'
#' @return Data.frame with columns for sample, chromosome (and arm if requested), CN mode, mean max posterior, and window count.
#'
#'
#' @export
summarize_cn_mode <- function(df,
                              level = c("sample", "chromosome", "chromosome-arm"),
                              centromeres = NULL,
                              cn_col = "CN_call",
                              post_col = "post_max") {
  level <- match.arg(level)

  # unwrap if hmm_CN-like object with $by_window
  dat <- if (is.data.frame(df)) df else if (!is.null(df$by_window) && is.data.frame(df$by_window)) df$by_window else
    stop("Input must be a data.frame or an hmm_CN-like object with a data.frame at `$by_window`.")

  # required columns
  req_cols <- c("Sample", "Chr", "Start", "End", cn_col)
  miss <- setdiff(req_cols, names(dat))
  if (length(miss)) stop("Missing required columns in data: ", paste(miss, collapse = ", "))

  if (!is.numeric(dat[[cn_col]]) && !is.integer(dat[[cn_col]])) {
    stop(sprintf("Column '%s' must be numeric/integer.", cn_col))
  }
  if (post_col %in% names(dat) && any(!is.na(dat[[post_col]])) && !is.numeric(dat[[post_col]])) {
    stop(sprintf("Column '%s' must be numeric.", post_col))
  }

  dat <- as.data.frame(dat)

  # chromosome-arm handling -> create chrID.1 (Chr) and chrID.2 (1=p, 2=q)
  if (level == "chromosome-arm") {
    if (is.null(centromeres))
      stop("For level='chromosome-arm', provide `centromeres` (named numeric vector or data.frame with Chr, Centromere).")

    if (is.data.frame(centromeres)) {
      if (!all(c("Chr", "Centromere") %in% names(centromeres)))
        stop("centromeres data.frame must have columns: Chr, Centromere")
      cm <- setNames(centromeres$Centromere, centromeres$Chr)
    } else if (is.numeric(centromeres) && !is.null(names(centromeres))) {
      cm <- centromeres
    } else {
      stop("centromeres must be a named numeric vector or a data.frame with Chr and Centromere.")
    }

    mid <- (dat$Start + dat$End) / 2
    has_cm <- dat$Chr %in% names(cm)

    # Assign arm index: 1 = p (left), 2 = q (right)
    chrID.2 <- rep(NA_integer_, nrow(dat))
    chrID.2[has_cm] <- ifelse(mid[has_cm] < cm[dat$Chr[has_cm]], 1L, 2L)

    if (any(!has_cm)) {
      missing_chr <- unique(dat$Chr[!has_cm])
      warning("No centromere provided for: ", paste(missing_chr, collapse = ", "),
              ". Rows for these chromosomes will be dropped.")
    }

    keep <- has_cm & !is.na(chrID.2)
    dat <- dat[keep, , drop = FALSE]
    chrID.2 <- chrID.2[keep]

    dat$chrID.1 <- dat$Chr
    dat$chrID.2 <- chrID.2
  }

  # grouping variables
  group_vars <- c("Sample")
  if (level %in% c("chromosome", "chromosome-arm")) group_vars <- c(group_vars, "Chr")
  if (level == "chromosome-arm") group_vars <- c(group_vars, "chrID.1", "chrID.2")

  # summarize (base R)
  split_idx <- interaction(dat[group_vars], drop = TRUE, lex.order = TRUE)
  grouped <- split(seq_len(nrow(dat)), split_idx)

  out <- lapply(grouped, function(idx) {
    g <- dat[idx, , drop = FALSE]
    row <- list(
      Sample        = g$Sample[1],
      CN_mode       = mode(g[[cn_col]]),
      mean_max_prob = if (post_col %in% names(g) && any(!is.na(g[[post_col]]))) {
        mean(g[[post_col]], na.rm = TRUE)
      } else {
        NA_real_
      },
      n_windows     = nrow(g)
    )
    if ("Chr" %in% group_vars)     row$Chr      <- g$Chr[1]
    if ("chrID.1" %in% group_vars) row$chrID.1  <- g$chrID.1[1]
    if ("chrID.2" %in% group_vars) row$chrID.2  <- g$chrID.2[1]
    as.data.frame(row, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)

  # nice column ordering
  want <- c("Sample",
            if (level %in% c("chromosome", "chromosome-arm")) "Chr",
            if (level == "chromosome-arm") c("chrID.1", "chrID.2"),
            "CN_mode", "mean_max_prob", "n_windows")
  res <- res[want]

  rownames(res) <- NULL
  res
}

#' Merge CN summary with area-based ploidy estimates
#'
#' Merges summarized HMM CN results with area-based ploidy estimates for each sample and chromosome (or arm).
#' Useful for combining HMM and area-based results for reporting or downstream analysis.
#'
#' @param hmm_summarized Data.frame from summarize_cn_mode.
#' @param qploidy_area_ploidy_estimation Object containing area-based ploidy matrices.
#' @param level Character. Merge level: 'chromosome', 'chromosome-arm', or 'sample'.
#'
#' @return Data.frame with merged columns: sample, chromosome, HMM CN mode, area-based CN, confidence metrics, and window count.
#'
#' @importFrom dplyr rename
#' @importFrom tidyr pivot_longer as_tibble
#'
#' @export
merge_cn_summary_with_estimates <- function(hmm_summarized,
                                            qploidy_area_ploidy_estimation,
                                            level = c("chromosome", "chromosome-arm", "sample")) {
  level <- match.arg(level)

  if (!inherits(qploidy_area_ploidy_estimation, "qploidy_area_ploidy_estimation"))
    stop("`qploidy_area_ploidy_estimation` must be a qploidy_area_ploidy_estimation object.")

  ploidy_mat <- qploidy_area_ploidy_estimation$ploidy
  diff_mat   <- qploidy_area_ploidy_estimation$diff_first_second

  if (is.null(rownames(ploidy_mat)) || is.null(rownames(diff_mat)))
    stop("Matrices in `qploidy_area_ploidy_estimation` must have rownames = sample IDs.")

  mat_long <- function(M, value_name) {
    as_tibble(M, rownames = "Sample") |>
      pivot_longer(cols = -Sample, names_to = "Chr", values_to = value_name)
  }

  if (level == "chromosome") {
    est_ploidy_long <- mat_long(ploidy_mat, "CN_area")
    est_diff_long   <- mat_long(diff_mat,   "area_diff_prob")
    est_long <- merge(est_ploidy_long, est_diff_long, by = c("Sample", "Chr"), all = TRUE, sort = FALSE)

    # hmm_summarized is expected to have CN_mode, mean_max_prob, n_windows from summarize_cn_mode()
    out <- merge(hmm_summarized, est_long, by = c("Sample", "Chr"), all.x = TRUE, sort = FALSE)

  } else if (level == "chromosome-arm") {
    if (!all(c("Sample", "chrID.1", "chrID.2") %in% names(hmm_summarized)))
      stop("For level='chromosome-arm', `hmm_summarized` must include: Sample, chrID.1, chrID.2.")
    hmm_summarized$ChrArm <- paste(hmm_summarized$chrID.1, hmm_summarized$chrID.2, sep = ".")

    est_ploidy_long <- mat_long(ploidy_mat, "CN_area") |>
      rename(ChrArm = Chr)
    est_diff_long   <- mat_long(diff_mat,   "area_diff_prob") |>
      rename(ChrArm = Chr)

    est_long <- merge(est_ploidy_long, est_diff_long, by = c("Sample", "ChrArm"), all = TRUE, sort = FALSE)
    out <- merge(hmm_summarized, est_long, by.x = c("Sample", "ChrArm"), by.y = c("Sample", "ChrArm"),
                 all.x = TRUE, sort = FALSE)

    out$chrID.1 <- NULL
    out$chrID.2 <- NULL
    out$ChrArm  <- NULL

  } else { # level == "sample"
    # Aggregate matrices per sample
    CN_area         <- apply(ploidy_mat, 1, function(x) mode(as.numeric(x)))
    area_diff_prob  <- rowMeans(diff_mat, na.rm = TRUE)

    est_sample <- tibble(
      Sample         = names(CN_area),
      CN_area        = as.numeric(CN_area),
      area_diff_prob = as.numeric(area_diff_prob)
    )

    out <- merge(hmm_summarized, est_sample, by = "Sample", all.x = TRUE, sort = FALSE)
    if (!("Chr" %in% names(out))) out$Chr <- NA_character_
    out <- out[c("Sample", "Chr", setdiff(names(out), c("Sample", "Chr")))]
  }

  # ---- Rename to requested schema ----
  # From summarize_cn_mode(): CN_mode -> CN_HMM; mean_max_prob -> HMM_max_prob; n_windows -> HMM_n_windows
  rename_map <- c(CN_mode = "CN_HMM",
                  mean_max_prob = "HMM_max_prob",
                  n_windows = "HMM_n_windows")
  for (old in names(rename_map)) {
    if (old %in% names(out)) names(out)[names(out) == old] <- rename_map[[old]]
  }

  # Ensure all columns exist and ordered as requested
  want <- c("Sample", "Chr", "CN_HMM", "CN_area", "HMM_max_prob", "area_diff_prob", "HMM_n_windows")
  missing <- setdiff(want, names(out))
  for (m in missing) out[[m]] <- NA
  out <- out[want]

  rownames(out) <- NULL
  out
}

#' Export chromosome-level copy number calls in MAPpoly-compatible format
#'
#' Summarizes the modal copy number per sample and chromosome from an HMM CN object and
#' returns a wide data.frame (samples x chromosomes) suitable for use with MAPpoly.
#'
#' @param hmm_CN_obj An object of class `hmm_CN` as returned by `hmm_estimate_CN_multi`,
#'   or a data.frame with columns `Sample`, `Chr`, `Start`, `End`, `CN_call`, and `post_max`.
#'
#' @return A data.frame with samples as rows and chromosomes as columns, where each cell
#'   contains the modal copy number call for that sample-chromosome combination.
#'
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#'
#' @export
export_mappoly <- function(hmm_CN_obj){
  cn_chr <- summarize_cn_mode(hmm_CN_obj, level = "chromosome")

  aneu <- cn_chr |>
    select(Sample, Chr, CN_mode) |>
    pivot_wider(names_from = Chr, values_from = CN_mode) |>
    column_to_rownames("Sample") |>
    as.data.frame()
  
  return(aneu)
}

#' Import PennCNV output into an HMM result object
#'
#' Parses PennCNV segment output and converts it into the internal `hmm_CN`
#' structure used by Qploidy2, producing both window-level and marker-level
#' summaries compatible with downstream plotting and reporting functions.
#'
#' @param df A data.frame containing PennCNV output rows. It is expected to
#'   include the segment identifier, chromosome, SNP count, CN state, and
#'   confidence columns used by the PennCNV output format.
#' @param geno.pos A data.frame with `MarkerName`, `Chromosome`, and `Position`
#'   columns used to map PennCNV segments back to marker-level genomic positions.
#'
#' @return An object of class `hmm_CN` with populated `by_window` and
#'   `by_marker` entries.
#'
#' @importFrom dplyr mutate row_number left_join rename filter select
#' @importFrom stringr str_extract
#'
#' @export
import_penncnv_results <- function(df, geno.pos) {

  # --- Parse PennCNV columns into by_window format ---
  penncnv_converted <- df %>%
    mutate(
      Sample         = str_extract(V5, "[^/]+$"),
      Chr            = str_extract(V1, "(?<=chr)\\d+"),
      Start          = as.integer(str_extract(V1, "(?<=:)\\d+")),
      End            = as.integer(str_extract(V1, "(?<=-)[0-9]+")),
      n_snps         = as.integer(str_extract(V2, "\\d+")),
      CN_call        = as.integer(str_extract(V4, "(?<=cn=)\\d+")),
      CN_reliability = as.numeric(str_extract(V8, "(?<=conf=)[0-9.]+"))
    ) %>%
    mutate(WindowID = row_number(), .by = c(Sample, Chr)) %>%
    mutate(
      n_het    = NA_integer_,
      z        = NA_real_,
      w_baf    = NA_real_,
      post_max = NA_real_,
      post_CN1 = NA_real_,
      post_CN2 = NA_real_,
      post_CN3 = NA_real_,
      post_CN4 = NA_real_,
      post_CN5 = NA_real_,
      post_CN6 = NA_real_
    ) %>%
    select(Sample, Chr, WindowID, Start, End, n_snps, n_het,
           z, w_baf, CN_call, CN_reliability, post_max,
           post_CN1, post_CN2, post_CN3, post_CN4, post_CN5, post_CN6)

  # --- Expand into by_marker format ---
  penncnv_by_marker <- penncnv_converted %>%
    mutate(Chr = as.integer(Chr)) %>%
    left_join(
      geno.pos %>% rename(Chr = Chromosome) %>% mutate(Chr = as.integer(Chr)),
      by = "Chr",
      relationship = "many-to-many"
    ) %>%
    filter(Position >= Start & Position <= End) %>%
    mutate(
      `.__w__`   = WindowID,
      SampleName = Sample,
      X          = NA_real_,
      Y          = NA_real_,
      R          = 0,
      ratio      = 0.5,
      geno       = NA_integer_,
      baf        = 0.5,
      z          = 0,
      outlier    = FALSE,
      w_baf      = 0,
      post_max   = 0
    ) %>%
    select(Chr, `.__w__`, MarkerName, SampleName,
           X, Y, R, ratio, geno, baf,
           z, Position, outlier, w_baf,
           CN_call, post_max, CN_reliability)

  # --- Assemble into hmm_CN object ---
  params_default <- list(
    cn_grid = NA_integer_,
    exp_ploidy = NA_real_,
    min_snps_per_window = NA_real_,
    distribution = NA_character_,
    mu = NA_real_,
    sigma = NA_real_,
    loglik = NA_real_
  )

  sample_ids <- unique(penncnv_converted$Sample)
  params_samples <- stats::setNames(
    replicate(length(sample_ids), params_default, simplify = FALSE),
    sample_ids
  )

  fake_hmm <- list(
    by_window = penncnv_converted,
    by_marker = penncnv_by_marker,
    params = NULL,
    params_samples = params_samples
  )

  return(structure(fake_hmm, class = "hmm_CN"))
}