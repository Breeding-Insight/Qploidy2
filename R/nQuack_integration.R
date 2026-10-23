#' Convert Qploidy Standardization Data to nQuack Format
#'
#' Converts a Qploidy standardization object (or its \code{$data} slot) into
#' the list-of-matrices format required by the
#' \href{https://github.com/mgaynor1/nQuack}{nQuack} ploidy estimation
#' package. For each sample, heterozygous sites are extracted and one allele
#' count is randomly selected (50/50) to meet nQuack's input expectations.
#'
#' @param qploidy_standardization_object An object of class
#'   \code{qploidy_standardization} or a plain \code{data.frame} containing at
#'   least the columns \code{SampleName}, \code{R} (total read depth),
#'   \code{ratio}, and \code{geno}. When \code{from_standardized_baf = TRUE}
#'   (default) a \code{baf} column (standardized BAF) is also required; when
#'   \code{FALSE}, \code{X} (reference allele count) and \code{Y} (alternative
#'   allele count) are required instead.
#' @param min_depth Integer. Minimum total read depth (\code{R}) required for
#'   a site to be retained. Default is \code{10}.
#' @param max_geno Integer. Maximum genotype dosage considered homozygous
#'   (inclusive). Sites where \code{geno == 0} or \code{geno == max_geno} are
#'   filtered as homozygous. Default is \code{4} (tetraploid).
#' @param seed Integer. Random seed for reproducibility of the random allele selection. Default is \code{123L}.
#' @param from_standardized_baf Logical. If \code{TRUE} (default), allele
#'   counts are derived by inverting the standardization-corrected \code{baf}
#'   column back into reference/alternative counts (\code{alt = round(baf * R)},
#'   \code{ref = R - alt}), the same approach used internally by
#'   \code{.nquack_estimate_ploidy()}. If \code{FALSE}, the raw \code{X}/\code{Y}
#'   allele counts are used directly instead.
#'
#' @return A named \code{list} with one element per sample. Each element is a
#'   two-column integer \code{matrix}:
#'   \itemize{
#'     \item Column 1: total read depth (\code{R}).
#'     \item Column 2: read count of one randomly chosen allele, selected with
#'       equal probability per site (reference/alternative counts inverted
#'       from \code{baf} when \code{from_standardized_baf = TRUE}, or raw
#'       \code{X}/\code{Y} otherwise).
#'   }
#'   Samples with no qualifying sites after filtering return an empty matrix
#'   with zero rows.
#'
#' @details
#' Filtering applied per site (per sample):
#' \enumerate{
#'   \item Sites with \code{ratio == 0} or \code{ratio == 1} (homozygous by
#'     ratio) are removed.
#'   \item Sites with \code{geno == 0} or \code{geno == max_geno}
#'     (homozygous by genotype call) are removed.
#'   \item Sites with total depth \code{R < min_depth} or with a zero count on
#'     the randomly chosen allele are removed.
#' }
#' The random allele selection is vectorised (one \code{sample()} call per
#' sample rather than a per-site loop), making the function efficient for large
#' datasets.
#' 
#' @author adapted from Michelle Gaynor
#' 
#' @seealso \href{https://mlgaynor.com/nQuack/articles/Qploidy2nQuack.html}{nQuack tutorial} for EM model selection using the nQuack package.
#'
#'
#' @export
convert2nQuack <- function(qploidy_standardization_object,
                           min_depth = 10L,
                           max_geno  = NULL,
                           seed      = 123L,
                           from_standardized_baf = TRUE) {

                        
  # Accept qploidy_standardization objects
  if(is.null(max_geno)) {
    max_geno <- as.numeric(qploidy_standardization_object$info["ploidy.standardization"])
  }
 
  if (inherits(qploidy_standardization_object, "qploidy_standardization")) {
    qploidy_standardization_object <- qploidy_standardization_object$data
  }

  dat <- qploidy_standardization_object
  required_cols <- c("SampleName", "R", "ratio", "geno")
  required_cols <- if (from_standardized_baf) c(required_cols, "baf") else c(required_cols, "X", "Y")
  missing_cols  <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    stop("Input is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  samples  <- unique(dat$SampleName)
  templist <- vector("list", length(samples))
  names(templist) <- samples

  for (i in seq_along(samples)) {
    temp <- dat[dat$SampleName == samples[i], ]

    # Remove homozygous sites by ratio and genotype call
    temp <- temp[!is.na(temp$ratio) & temp$ratio != 0 & temp$ratio != 1, ]
    temp <- temp[!is.na(temp$geno)  & temp$geno  != 0 & temp$geno  != max_geno, ]

    if (nrow(temp) == 0L) {
      templist[[i]] <- matrix(integer(0), nrow = 0L, ncol = 2L)
      next
    }

    # Randomly select ref/alt coverage per site (vectorised)
    set.seed(seed)  # For reproducibility
    n       <- nrow(temp)
    use_ref <- sample(c(TRUE, FALSE), size = n, replace = TRUE)

    if (from_standardized_baf) {
      # Invert standardized BAF back into allele counts (mirrors .nquack_estimate_ploidy()),
      # so nQuack receives standardization-corrected allele frequencies instead of raw X/Y.
      alt_count <- as.integer(round(temp$baf * temp$R))
      ref_count <- as.integer(temp$R) - alt_count
      allele    <- ifelse(use_ref, ref_count, alt_count)
    } else {
      allele <- ifelse(use_ref, temp$X, temp$Y)
    }

    xmr <- cbind(R = temp$R, allele = allele)

    # Apply depth and non-zero allele filters
    keep <- !is.na(xmr[, 1L]) & xmr[, 1L] >= min_depth &
            !is.na(xmr[, 2L]) & xmr[, 2L] > 0L
    templist[[i]] <- as.matrix(xmr[keep, , drop = FALSE])
  }

  templist
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Internal helper: per-window BAF log-likelihoods via nQuack EM
#
# When nQuack_EM = TRUE, this replaces compute_baf_likelihoods for each window.
# Instead of a kernel-density template with bandwidth bw, we run the nQuack
# emstep function for each CN state in cn_grid using the same distribution
# family that won the global nQuack pass (nq_dist).
#
# Parameters:
#   win_df       - single-row data.frame from win_df (Chr, WindowID, Start, End)
#   d_all        - full marker data.frame (all chromosomes/windows)
#   baf_col      - BAF column name in d_all
#   cn_grid      - CN states to evaluate (2:6 only supported by nQuack)
#   nq_dist      - winning distribution from global pass, e.g. "beta-uniform"
#   het_range    - c(lo, hi) — sites outside this range are excluded
#   min_het_frac - minimum fraction of het sites required (else ll = -Inf)
#
# Returns a list with:
#   ll_vec   - named numeric vector of log-likelihoods, one per CN in cn_grid
#   prob_vec - softmax probabilities (same length)
#   no_one   - logical, TRUE if only CN=1 evaluated (fallback)
# ---------------------------------------------------------------------------
#' @keywords internal
.nquack_baf_likelihoods_by_window <- function(win_df,
                                              d_all,
                                              baf_col      = "baf",
                                              cn_grid      = 2:6,
                                              nq_dist      = "beta",
                                              het_range    = c(0.05, 0.95),
                                              min_het_frac = 0.01) {
  # Derive base family and uniform flag from the winning distribution string
  nq_dist_lo  <- tolower(nq_dist)
  add_uniform <- grepl("-uniform$", nq_dist_lo)
  base_family <- sub("-uniform$", "", nq_dist_lo)  # "normal", "beta", "beta-binomial"

  # nQuack-supported CN range (2:6); unsupported states get ll = -Inf
  cn_to_idx <- c("2" = 1L, "3" = 2L, "4" = 3L, "5" = 4L, "6" = 5L)
  niter   <- 1000L
  epsilon <- 0.1
  trunc   <- c(0.0, 0.0)

  # Pre-built initial parameter sets (same as nQuack internals, svec = 0.01)
  set <- vector("list", 5L)
  set[[1]] <- list(avec = c(1.00),                            mvec = c(0.50),                            svec = c(0.01))
  set[[2]] <- list(avec = c(0.50, 0.50),                      mvec = c(0.67, 0.33),                      svec = c(0.01, 0.01))
  set[[3]] <- list(avec = c(0.33, 0.33, 0.33),                mvec = c(0.25, 0.50, 0.75),                svec = c(0.01, 0.01, 0.01))
  set[[4]] <- list(avec = c(0.25, 0.25, 0.25, 0.25),          mvec = c(0.20, 0.40, 0.60, 0.80),          svec = c(0.01, 0.01, 0.01, 0.01))
  set[[5]] <- list(avec = c(0.20, 0.20, 0.20, 0.20, 0.20),    mvec = c(0.17, 0.33, 0.50, 0.67, 0.83),    svec = c(0.01, 0.01, 0.01, 0.01, 0.01))

  setU <- vector("list", 5L)
  setU[[1]] <- list(avec = c(0.9, 0.1),                                mvec = c(0.50),                          svec = c(0.01))
  setU[[2]] <- list(avec = c(0.45, 0.45, 0.1),                         mvec = c(0.67, 0.33),                    svec = c(0.01, 0.01))
  setU[[3]] <- list(avec = c(0.3, 0.3, 0.3, 0.1),                      mvec = c(0.25, 0.50, 0.75),              svec = c(0.01, 0.01, 0.01))
  setU[[4]] <- list(avec = c(0.225, 0.225, 0.225, 0.225, 0.1),         mvec = c(0.20, 0.40, 0.60, 0.80),        svec = c(0.01, 0.01, 0.01, 0.01))
  setU[[5]] <- list(avec = c(0.18, 0.18, 0.18, 0.18, 0.18, 0.1),       mvec = c(0.17, 0.33, 0.50, 0.67, 0.83),  svec = c(0.01, 0.01, 0.01, 0.01, 0.01))

  # Subset markers for this window
  win_chr   <- win_df$Chr[1]
  win_start <- win_df$Start[1]
  win_end   <- win_df$End[1]
  d_win <- d_all[!is.na(d_all$Chr) & d_all$Chr == win_chr &
                 !is.na(d_all$Position) &
                 d_all$Position >= win_start & d_all$Position <= win_end, , drop = FALSE]

  # Resolve BAF/ratio column
  bc <- if (baf_col %in% colnames(d_win)) baf_col else if ("ratio" %in% colnames(d_win)) "ratio" else NA_character_
  has_depth <- "R" %in% colnames(d_win)

  # Build xi (allele freq) and xm (matrix) for het sites
  make_empty <- function() {
    # No/insufficient het evidence: the window is uninformative, not impossible.
    # Equal finite log-likelihoods (neutral) let depth (z) alone drive inference;
    # -Inf here would make every state "impossible" and corrupt forward/backward sums.
    ll_vec   <- setNames(rep(0, length(cn_grid)), paste0("CN", cn_grid))
    prob_vec <- setNames(rep(1 / length(cn_grid), length(cn_grid)), paste0("CN", cn_grid))
    list(ll_vec = ll_vec, prob_vec = prob_vec, no_one = FALSE)
  }

  if (is.na(bc) || nrow(d_win) == 0L) return(make_empty())

  baf_raw <- d_win[[bc]]
  valid   <- !is.na(baf_raw) & baf_raw >= het_range[1] & baf_raw <= het_range[2]
  if (sum(valid) < 2L) return(make_empty())
  if (sum(valid) / length(baf_raw) < min_het_frac) return(make_empty())

  xi_all <- baf_raw[valid]

  # For beta-binomial we need R (total depth) and the allele count. Always
  # reconstruct from the standardized/corrected baf_col (xi_all), never the raw
  # Y column, so counts stay consistent with the BAF used for model selection
  # and heterozygote filtering elsewhere.
  xm_all <- NULL
  if (base_family == "beta-binomial") {
    if (!has_depth) {
      warning(".nquack_baf_likelihoods_by_window: 'R' column required for beta-binomial; falling back to -Inf likelihoods.")
      return(make_empty())
    }
    R_vals  <- d_win$R[valid]
    alt_raw <- round(xi_all * R_vals)
    ok_depth <- !is.na(R_vals) & R_vals > 0 & !is.na(alt_raw)
    if (sum(ok_depth) < 2L) return(make_empty())
    xi_all  <- xi_all[ok_depth]
    R_vals  <- as.integer(R_vals[ok_depth])
    alt_cnt <- as.integer(alt_raw[ok_depth])
    xm_all  <- cbind(R_vals, alt_cnt)
  }

  # Compute log-likelihood for each CN state
  ll_vec <- setNames(rep(-Inf, length(cn_grid)), paste0("CN", cn_grid))

  for (cn in cn_grid) {
    cn_char <- as.character(cn)
    idx     <- cn_to_idx[cn_char]
    if (is.na(idx)) next  # outside nQuack range

    parms <- if (add_uniform) setU[[idx]] else set[[idx]]

    ll <- tryCatch({
      if (base_family == "beta") {
        fn  <- if (add_uniform) nQuack::emstepBU else nQuack::emstepB
        res <- fn(parms, xi_all, niter, epsilon, trunc, "free")
        res$loglikelihood
      } else if (base_family == "normal") {
        fn  <- if (add_uniform) nQuack::emstepNU else nQuack::emstepN
        res <- fn(parms, xi_all, niter, epsilon, trunc, "free")
        res$loglikelihood
      } else if (base_family == "beta-binomial") {
        fn  <- if (add_uniform) nQuack::emstepBBU else nQuack::emstepBB
        res <- fn(parms, xm_all, niter, epsilon, trunc, "free")
        res$loglikelihood
      } else {
        -Inf
      }
    }, error = function(e) -Inf)

    ll_vec[paste0("CN", cn)] <- if (is.null(ll) || length(ll) == 0 || !is.finite(ll)) -Inf else ll
  }

  # Softmax probabilities
  ll_finite <- ll_vec[is.finite(ll_vec)]
  if (length(ll_finite) == 0L) {
    prob_vec <- setNames(rep(1 / length(cn_grid), length(cn_grid)), names(ll_vec))
  } else {
    ll_shifted <- ll_vec - max(ll_vec[is.finite(ll_vec)])
    ex         <- exp(pmax(ll_shifted, -700))
    prob_vec   <- ex / sum(ex)
  }

  list(ll_vec = ll_vec, prob_vec = prob_vec, no_one = FALSE)
}

# ---------------------------------------------------------------------------
# Internal helper: run nQuack EM and return estimated ploidy as integer CN
#
# Strategy
# --------
# 1. Resolve BAF frequency column (baf preferred; ratio fallback).
# 2. Apply depth filter (R >= min_depth).
# 3. Apply het-range filter (het_range) to strip homozygous flanking sites —
#    this is the primary guard before passing data to quackBeta/quackBetaBinom.
# 4. Optionally refine with emstepNA(glreturn=TRUE): fit a fixed Normal EM
#    centred on expected dosage peaks; keep sites where rowSums($pir) >=
#    dosage_threshold.  Falls back gracefully if emstepNA fails.
# 5. Build xm = cbind(R, round(xi * R)); nQuack divides col2/col1 internally.
# 6. Derive ploidy name vector from cn_grid (2:6 only); warn for out-of-range.
# 7. Fit quackNormal / quackBeta / quackBetaBinom; select BIC winner via
#    quackit(); map ploidy name to integer CN.
# ---------------------------------------------------------------------------
#' @keywords internal
.nquack_estimate_ploidy <- function(d, sample_id, cn_grid,
                                    baf_col          = "baf",
                                    selected_model   = NULL,
                                    nQuack_cores     = 1L,
                                    het_range        = c(0.05, 0.95),
                                    dosage_threshold = 0.6,
                                    min_depth        = 10L,
                                    fixed_dist       = NULL,
                                    allowed_dists    = NULL,
                                    verbose          = TRUE) {
  if (!requireNamespace("nQuack", quietly = TRUE)) {
    stop("nQuack_EM = TRUE requires the 'nQuack' package. Install with:\n  remotes::install_github(\"mgaynor1/nQuack\")")
  }

  # --- warn if cn_grid is outside nQuack's supported range (2:6) ---
  out_of_range <- cn_grid[cn_grid < 2L | cn_grid > 6L]
  if (length(out_of_range) > 0L) {
    warning(sprintf(
      paste0("nQuack_EM = TRUE: cn_grid value(s) %s are outside nQuack's supported ",
             "range (2:6; haploid, deletion, and ploidy > 6 are not modelled). ",
             "These values will be ignored for ploidy estimation."),
      paste(out_of_range, collapse = ", ")
    ))
  }

  # --- 1. Resolve BAF / frequency column --------------------------------
  if (baf_col %in% names(d) && !all(is.na(d[[baf_col]]))) {
    freq_vec <- d[[baf_col]]
    vmsg("nQuack: using standardized BAF column '%s'", verbose = verbose, level = 1, type = ">>", baf_col)
  } else if ("ratio" %in% names(d)) {
    freq_vec <- d[["ratio"]]
    vmsg("nQuack: BAF column not found, falling back to 'ratio'", verbose = verbose, level = 1, type = ">>")
  } else {
    stop("nQuack_EM = TRUE: data must contain a 'baf' or 'ratio' column.")
  }

  if (!"R" %in% names(d))
    stop("nQuack_EM = TRUE: data must contain column 'R' (total depth).")

  # --- 2. Depth filter --------------------------------------------------
  depth_vec  <- d[["R"]]
  depth_keep <- is.finite(freq_vec) & is.finite(depth_vec) & depth_vec >= min_depth
  xi_depth   <- freq_vec[depth_keep]
  R_depth    <- depth_vec[depth_keep]

  if (length(xi_depth) == 0L)
    stop(sprintf(
      paste0("nQuack_EM = TRUE: no sites with depth >= %d and valid BAF/ratio. ",
             "Check that column 'R' contains raw read counts and '%s' contains valid allele frequencies."),
      min_depth, baf_col
    ))

  vmsg("nQuack: %d sites pass depth >= %d filter", verbose = verbose, level = 1, type = ">>",
       length(xi_depth), min_depth)

  # --- 3. Primary het-range filter --------------------------------------
  # Strip homozygous flanking sites before model fitting.
  # This mirrors het_range used elsewhere in hmm_estimate_CN and is
  # essential to prevent quackBeta from failing at initial parameters.
  hr_keep <- xi_depth >= het_range[1L] & xi_depth <= het_range[2L]
  xi_het  <- xi_depth[hr_keep]
  R_het   <- R_depth[hr_keep]

  vmsg("nQuack: %d / %d sites in het_range [%.2f, %.2f]",
       verbose = verbose, level = 1, type = ">>",
       length(xi_het), length(xi_depth), het_range[1L], het_range[2L])

  if (length(xi_het) < 10L)
    stop(sprintf(
      paste0("nQuack_EM = TRUE: only %d sites remain after depth + het_range filter; ",
             "need >= 10 to fit a mixture model."),
      length(xi_het)
    ))

  # --- 4. Optional emstepNA refinement (secondary het filter) -----------
  # Fit a fixed Normal EM centred on expected dosage peaks; keep sites where
  # rowSums(pir) >= dosage_threshold.  Falls back gracefully if unavailable.
  init_cn <- if (!is.null(selected_model) && !is.null(selected_model$best$best_cn)) {
    as.integer(selected_model$best$best_cn)
  } else {
    cn_grid_valid <- cn_grid[cn_grid >= 2L & cn_grid <= 6L]
    if (length(cn_grid_valid) == 0L) 4L else as.integer(stats::median(cn_grid_valid))
  }
  init_cn <- max(2L, min(6L, init_cn))

  n_het_comp    <- init_cn - 1L
  alpha_init    <- rep(1 / n_het_comp, n_het_comp)
  mean_init     <- seq(1, init_cn - 1) / init_cn
  var_init      <- rep(0.01, n_het_comp)
  parmlist_init <- list(avec = alpha_init, mvec = mean_init, svec = var_init)

  pir_mat <- tryCatch({
    res <- nQuack::emstepNA(
      parmlist = parmlist_init,
      xi       = xi_het,
      niter    = 500L,
      epsilon  = 1e-4,
      trunc    = het_range,
      type     = "fixed",
      glreturn = TRUE
    )
    if (is.list(res) && !is.null(res$pir)) res$pir else NULL
  }, error = function(e) {
    vmsg("nQuack: emstepNA refinement failed (%s); using het_range filter only",
         verbose = verbose, level = 1, type = ">>", conditionMessage(e))
    NULL
  })

  if (!is.null(pir_mat) && is.matrix(pir_mat) && nrow(pir_mat) == length(xi_het)) {
    em_keep <- rowSums(pir_mat, na.rm = TRUE) >= dosage_threshold
    n_em    <- sum(em_keep, na.rm = TRUE)
    vmsg("nQuack: emstepNA kept %d / %d het sites (dosage_threshold = %.2f)",
         verbose = verbose, level = 1, type = ">>",
         n_em, length(xi_het), dosage_threshold)
    if (isTRUE(n_em >= 10L)) {
      xi_het <- xi_het[em_keep]
      R_het  <- R_het[em_keep]
    } else {
      vmsg("nQuack: emstepNA kept < 10 sites; reverting to het_range filter only",
           verbose = verbose, level = 1, type = ">>")
    }
  }

  if (length(xi_het) < 10L)
    stop("nQuack_EM = TRUE: fewer than 10 het sites remain after all filters; cannot fit mixture model.")

  # --- 5. Build xm matrix -----------------------------------------------
  # nQuack expects allele_count to be a *randomly sampled* allele (ref or alt)
  # at each site, not always the alt allele.  Without this random flip the BAF
  # distribution fed to nQuack is asymmetric whenever the sample has more
  # ref-biased het sites than alt-biased ones, which causes nQuack to mistake
  # a diploid (one peak at 0.5) for a higher-ploidy mixture.
  # Solution: for each site independently flip to the ref count with P=0.5,
  # mirroring the approach in the nQuack Qploidy2 vignette (line 98).
  alt_count <- as.integer(round(xi_het * R_het))
  ref_count <- as.integer(R_het) - alt_count
  use_ref   <- sample(c(TRUE, FALSE), size = length(alt_count), replace = TRUE)
  xm_het <- cbind(
    R            = as.integer(R_het),
    allele_count = ifelse(use_ref, ref_count, alt_count)
  )
  vmsg("nQuack: fitting models on %d het sites", verbose = verbose, level = 1, type = ">>", nrow(xm_het))

  # --- 6. Derive ploidy name vector from cn_grid ------------------------
  cn_to_name      <- c("2" = "diploid", "3" = "triploid", "4" = "tetraploid",
                        "5" = "pentaploid", "6" = "hexaploid")
  grid_valid      <- cn_grid[cn_grid >= 2L & cn_grid <= 6L]
  nQuack_mixtures <- unname(cn_to_name[as.character(grid_valid)])
  nQuack_mixtures <- nQuack_mixtures[!is.na(nQuack_mixtures)]
  if (length(nQuack_mixtures) == 0L)
    stop("nQuack_EM = TRUE: cn_grid contains no values in nQuack's supported range (2:6).")
  vmsg("nQuack mixtures: %s", verbose = verbose, level = 1, type = ">>",
       paste(nQuack_mixtures, collapse = ", "))

  # --- 7. Fit quackNormal / quackBeta / quackBetaBinom ------------------
  # If fixed_dist is set (second pass reusing first-pass winner), only run that distribution.
  # fixed_dist is the full out_nq$distribution value e.g. "beta-uniform", "normal", "beta-binomial"
  # Gate each quack* call on whether fixed_dist (if set) matches that family.
  fixed_base    <- if (!is.null(fixed_dist)) sub("-uniform$", "", fixed_dist) else NULL
  # allowed_dists: user-supplied restriction on which families to evaluate
  # e.g. c("normal") → only Normal/Normal+Uniform; NULL → all three families
  allowed_lo <- if (!is.null(allowed_dists)) tolower(trimws(allowed_dists)) else NULL
  run_normal    <- (is.null(fixed_base) || fixed_base == "normal") &&
                   (is.null(allowed_lo)  || "normal" %in% allowed_lo)
  run_beta      <- (is.null(fixed_base) || fixed_base == "beta") &&
                   (is.null(allowed_lo)  || "beta" %in% allowed_lo)
  run_betabinom <- (is.null(fixed_base) || fixed_base == "beta-binomial") &&
                   (is.null(allowed_lo)  || "beta-binomial" %in% allowed_lo)
  if (!run_normal && !run_beta && !run_betabinom)
    stop(sprintf(
      paste0("nQuack_EM = TRUE: nQuack_dists = c(%s) leaves no distribution to evaluate. ",
             "Valid values: 'normal', 'beta', 'beta-binomial'."),
      paste(dQuote(allowed_dists), collapse = ", ")
    ))
  if (!is.null(allowed_dists))
    vmsg("nQuack: restricting distributions to: %s",
         verbose = verbose, level = 1, type = ">>",
         paste(intersect(allowed_lo, c("normal","beta","beta-binomial")), collapse = ", "))
  if (!is.null(fixed_dist))
    vmsg("nQuack: reusing winning distribution from first pass: %s",
         verbose = verbose, level = 1, type = ">>", fixed_dist)
  out_nq <- tryCatch({
    out_parts <- list()
    if (run_normal)
      out_parts$normal   <- nQuack::quackNormal(xm = xm_het, samplename = sample_id,
                                                cores = nQuack_cores, parallel = (nQuack_cores > 1L),
                                                verbose = FALSE)
    if (run_beta)
      out_parts$beta     <- nQuack::quackBeta(xm = xm_het, samplename = sample_id,
                                              cores = nQuack_cores, parallel = (nQuack_cores > 1L),
                                              verbose = FALSE)
    if (run_betabinom)
      out_parts$betabinom <- nQuack::quackBetaBinom(xm = xm_het, samplename = sample_id,
                                                    cores = nQuack_cores, parallel = (nQuack_cores > 1L),
                                                    verbose = FALSE)
    do.call(rbind, out_parts)
  }, error = function(e) {
    stop(sprintf("nQuack_EM = TRUE: nQuack model fitting failed: %s", conditionMessage(e)))
  })

  # --- 8. Select best model by BIC via quackit() ------------------------
  nq_summary  <- nQuack::quackit(model_out         = out_nq,
                                  summary_statistic = "BIC",
                                  mixtures          = nQuack_mixtures)
  winner_name <- nq_summary$winnerBIC[1L]
  vmsg("nQuack BIC winner: %s", verbose = verbose, level = 1, type = ">>", winner_name)

  name_to_cn <- c("diploid" = 2L, "triploid" = 3L, "tetraploid" = 4L,
                   "pentaploid" = 5L, "hexaploid" = 6L)
  nq_cn <- unname(name_to_cn[winner_name])
  if (is.na(nq_cn))
    stop(sprintf("nQuack_EM = TRUE: winner '%s' could not be mapped to a CN value.", winner_name))

  # Identify winning Distribution type (Normal / Beta / Beta-Binomial) for caching.
  # Find in the full model output the row for winner_name with minimum BIC.
  winner_rows <- out_nq[!is.na(out_nq$mixture) & out_nq$mixture == winner_name &
                         !is.na(out_nq$type)    & out_nq$type   != "free", ]
  if (nrow(winner_rows) > 0) {
    best_idx    <- which.min(winner_rows$BIC)
    winner_dist <- winner_rows$distribution[best_idx]
  } else {
    winner_dist <- NULL
  }
  # --- 9. Variance sanity check ----------------------------------------
  # A diploid sample with high BAF variance (e.g. low coverage, noisy mapping)
  # produces one wide peak at 0.5.  With fixed means, nQuack's normal model
  # fits this by assigning a large shared variance; hexaploid wins by BIC simply
  # because 5 components can absorb variance better than 1.  We detect this
  # artefact by refitting the diploid fixed_2 model and inspecting the fitted
  # sigma.  If sigma > nq_sigma_max (default 0.10) the nQuack call is
  # unreliable and we override to diploid.
  nq_sigma_max <- 0.10  # sigma threshold above which nQuack call is distrusted
  if (nq_cn != 2L && run_normal) {
    dip_sigma <- tryCatch({
      dip_res <- nQuack::emstepNA(
        parmlist = list(avec = c(1.0), mvec = c(0.5), svec = c(0.01)),
        xi       = xm_het[, 2L] / xm_het[, 1L],
        niter    = 500L,
        epsilon  = 1e-4,
        trunc    = het_range,
        type     = "fixed_2",
        glreturn = FALSE
      )
      fitted_svec <- dip_res$parm.list$svec
      if (is.numeric(fitted_svec) && length(fitted_svec) > 0L) mean(fitted_svec) else NA_real_
    }, error = function(e) NA_real_)

    if (!is.na(dip_sigma) && is.finite(dip_sigma) && dip_sigma > nq_sigma_max) {
      vmsg(paste0("nQuack: diploid sigma = %.4f > %.2f threshold; ",
                  "BAF variance too high to distinguish ploidy reliably. ",
                  "Overriding BIC winner ('%s') to diploid."),
           verbose = verbose, level = 1, type = ">>",
           dip_sigma, nq_sigma_max, winner_name)
      nq_cn       <- 2L
      winner_name <- "diploid"
      # Re-derive winner_dist from diploid rows
      dip_rows <- out_nq[!is.na(out_nq$mixture) & out_nq$mixture == "diploid" &
                           !is.na(out_nq$type)   & out_nq$type   != "free", ]
      winner_dist <- if (nrow(dip_rows) > 0L) dip_rows$distribution[which.min(dip_rows$BIC)] else winner_dist
    } else {
      vmsg("nQuack: diploid sigma = %.4f (threshold %.2f) — variance OK, keeping BIC winner '%s'",
           verbose = verbose, level = 1, type = ">>",
           if (is.na(dip_sigma)) NaN else dip_sigma, nq_sigma_max, winner_name)
    }
  }

  list(cn = nq_cn, winner = winner_name, dist = winner_dist)
}

# ---------------------------------------------------------------------------
#' Assign BAF dosages using nQuack EM posterior probabilities
#'
#' An alternative to \code{call_BAF_dosages} for the \code{nQuack_EM = TRUE}
#' path.  Runs the appropriate nQuack EM step with \code{glreturn = TRUE} to
#' obtain per-observation posterior component-membership probabilities
#' (\code{pir}), then maps them to integer dosages.
#'
#' @details
#' nQuack mixture components correspond to heterozygous dosage peaks
#' (\code{mvec} positions 1/cn, 2/cn, …, (cn-1)/cn) plus, when a uniform
#' component is present, one final column for the uniform class.  Observations
#' outside \code{het_range} are considered homozygous and assigned dosage 0
#' (BAF < 0.5) or \code{cn} (BAF ≥ 0.5) without involving the mixture model.
#' Observations inside \code{het_range} are assigned the dosage corresponding
#' to \code{argmax(pir[i, 1:ncol_het])} (the uniform column, if present, is
#' excluded from this argmax so that uniform-assigned observations still map to
#' the nearest dosage peak).
#'
#' @param baf_vec   Numeric vector of BAF / allele-frequency values (all
#'   markers in the window, not pre-filtered).
#' @param cn        Integer copy number for the window.
#' @param nq_dist   Distribution family string as returned by
#'   \code{.nquack_estimate_ploidy()} (e.g.\ \code{"beta"},
#'   \code{"beta-uniform"}, \code{"beta-binomial"}).
#' @param het_range Numeric length-2 vector: lower and upper BAF boundaries
#'   for heterozygous markers.  Default \code{c(0.05, 0.95)}.
#' @param R_vec     Integer vector of total read depths, same length as
#'   \code{baf_vec}.  Required when \code{nq_dist} contains
#'   \code{"beta-binomial"}. The alt allele count is always reconstructed as
#'   \code{round(baf_vec * R_vec)}, keeping it consistent with the
#'   standardized/corrected \code{baf_vec} rather than raw read counts.
#' @param max_cn    Optional integer.  When provided, extra
#'   \code{dosage_<k>} columns (filled with \code{NA}) are appended for
#'   dosages \code{cn+1} through \code{max_cn}, matching the column layout
#'   returned by \code{call_BAF_dosages}.
#' @param niter     EM maximum iterations.  Default 1000.
#' @param epsilon   EM convergence tolerance.  Default 0.1.
#'
#' @return A list of class \code{"BAF_dosage_calls"} with element \code{$data},
#'   a \code{data.frame} with the same column layout as \code{call_BAF_dosages}:
#'   \code{BAF}, \code{dosage}, \code{max_prob}, \code{dosage_0}, …,
#'   \code{dosage_<cn>} (plus any padding columns up to \code{max_cn}).
#'
#' @keywords internal
call_BAF_nQuack_dosages <- function(baf_vec,
                                    cn,
                                    nq_dist    = "beta",
                                    het_range  = c(0.05, 0.95),
                                    R_vec      = NULL,
                                    max_cn     = NULL,
                                    niter      = 1000L,
                                    epsilon    = 0.1) {
  if (!requireNamespace("nQuack", quietly = TRUE))
    stop("call_BAF_nQuack_dosages requires the 'nQuack' package.")

  cn <- as.integer(cn)
  stopifnot(is.numeric(baf_vec), cn >= 2L, cn <= 6L)

  nq_dist_lo  <- tolower(nq_dist)
  add_uniform <- grepl("-uniform$", nq_dist_lo)
  base_family <- sub("-uniform$", "", nq_dist_lo)   # "normal", "beta", "beta-binomial"

  trunc <- c(0.0, 0.0)
  cn_to_idx <- c("2" = 1L, "3" = 2L, "4" = 3L, "5" = 4L, "6" = 5L)
  idx <- cn_to_idx[as.character(cn)]
  if (is.na(idx))
    stop(sprintf("call_BAF_nQuack_dosages: cn=%d outside nQuack range (2-6).", cn))

  # Initial parameter sets (identical to nQuack internals)
  set <- vector("list", 5L)
  set[[1]] <- list(avec = c(1.00),                            mvec = c(0.50),                            svec = c(0.01))
  set[[2]] <- list(avec = c(0.50, 0.50),                      mvec = c(0.67, 0.33),                      svec = c(0.01, 0.01))
  set[[3]] <- list(avec = c(0.33, 0.33, 0.33),                mvec = c(0.25, 0.50, 0.75),                svec = c(0.01, 0.01, 0.01))
  set[[4]] <- list(avec = c(0.25, 0.25, 0.25, 0.25),          mvec = c(0.20, 0.40, 0.60, 0.80),          svec = c(0.01, 0.01, 0.01, 0.01))
  set[[5]] <- list(avec = c(0.20, 0.20, 0.20, 0.20, 0.20),    mvec = c(0.17, 0.33, 0.50, 0.67, 0.83),    svec = c(0.01, 0.01, 0.01, 0.01, 0.01))

  setU <- vector("list", 5L)
  setU[[1]] <- list(avec = c(0.9, 0.1),                                mvec = c(0.50),                          svec = c(0.01))
  setU[[2]] <- list(avec = c(0.45, 0.45, 0.1),                         mvec = c(0.67, 0.33),                    svec = c(0.01, 0.01))
  setU[[3]] <- list(avec = c(0.3, 0.3, 0.3, 0.1),                      mvec = c(0.25, 0.50, 0.75),              svec = c(0.01, 0.01, 0.01))
  setU[[4]] <- list(avec = c(0.225, 0.225, 0.225, 0.225, 0.1),         mvec = c(0.20, 0.40, 0.60, 0.80),        svec = c(0.01, 0.01, 0.01, 0.01))
  setU[[5]] <- list(avec = c(0.18, 0.18, 0.18, 0.18, 0.18, 0.1),       mvec = c(0.17, 0.33, 0.50, 0.67, 0.83),  svec = c(0.01, 0.01, 0.01, 0.01, 0.01))

  parms <- if (add_uniform) setU[[idx]] else set[[idx]]

  # The mvec gives the BAF peak positions for het dosages 1..(cn-1).
  # mvec is ordered high→low for diploid (0.67, 0.33) — normalise to ascending order
  # to map column index → dosage integer robustly via closest peak.
  mvec_raw  <- parms$mvec                         # length = cn-1
  mvec_ord  <- order(mvec_raw)                    # ascending sort order
  mvec_sort <- mvec_raw[mvec_ord]                 # sorted mvec (low→high BAF)
  # Expected peaks for dosages 1..(cn-1) in ascending order
  exp_peaks  <- seq_len(cn - 1L) / cn             # e.g. 1/4, 2/4, 3/4 for CN=4
  # Map each sorted mvec column to the nearest dosage integer
  col_to_dosage <- vapply(mvec_sort, function(m) which.min(abs(m - exp_peaks)), integer(1L))

  # Partition observations: outside het_range = homozygous (dosage 0 or cn)
  is_het <- !is.na(baf_vec) & baf_vec >= het_range[1] & baf_vec <= het_range[2]

  # Build output skeleton
  n_obs      <- length(baf_vec)
  dosage_vec <- integer(n_obs)
  max_prob   <- numeric(n_obs)
  prob_mat   <- matrix(0.0, nrow = n_obs, ncol = cn + 1L,
                       dimnames = list(NULL, paste0("dosage_", 0:cn)))

  # --- Homozygous assignments (outside het_range) ---
  hom_idx <- which(!is_het & !is.na(baf_vec))
  if (length(hom_idx) > 0L) {
    d_hom <- ifelse(baf_vec[hom_idx] < 0.5, 0L, cn)
    dosage_vec[hom_idx] <- d_hom
    max_prob[hom_idx]   <- 1.0
    for (k in seq_along(hom_idx)) {
      prob_mat[hom_idx[k], d_hom[k] + 1L] <- 1.0
    }
  }
  # NA baf → dosage NA, leave prob_mat row as 0
  na_idx <- which(is.na(baf_vec))
  dosage_vec <- replace(dosage_vec, na_idx, NA_integer_)

  # --- Heterozygous assignments via nQuack pir ---
  het_idx <- which(is_het)

  if (length(het_idx) >= 2L) {
    xi_het <- baf_vec[het_idx]

    pir <- tryCatch({
      if (base_family == "beta") {
        fn  <- if (add_uniform) nQuack::emstepBU else nQuack::emstepB
        res <- fn(parms, xi_het, niter, epsilon, trunc, "free", glreturn = TRUE)
        res$pir
      } else if (base_family == "normal") {
        fn  <- if (add_uniform) nQuack::emstepNU else nQuack::emstepN
        res <- fn(parms, xi_het, niter, epsilon, trunc, "free", glreturn = TRUE)
        res$pir
      } else if (base_family == "beta-binomial") {
        if (is.null(R_vec))
          stop("call_BAF_nQuack_dosages: R_vec required for beta-binomial.")
        R_het  <- as.integer(R_vec[het_idx])
        ok     <- !is.na(R_het) & R_het > 0L
        if (sum(ok) < 2L) stop("Insufficient depth data for beta-binomial dosage.")
        # Always reconstruct from baf_vec (standardized/corrected), never raw Y,
        # to stay consistent with model selection and heterozygote filtering.
        alt_het <- round(xi_het * R_het)
        xm_het <- cbind(R_het[ok], alt_het[ok])
        fn  <- if (add_uniform) nQuack::emstepBBU else nQuack::emstepBB
        res <- fn(parms, xm_het, niter, epsilon, trunc, "free", glreturn = TRUE)
        # pir rows align only to ok-subset; expand back
        pir_full <- matrix(NA_real_, nrow = length(het_idx), ncol = ncol(res$pir))
        pir_full[ok, ] <- res$pir
        pir_full
      } else {
        NULL
      }
    }, error = function(e) NULL)

    if (!is.null(pir) && nrow(pir) == length(het_idx)) {
      # Number of het-dosage columns (exclude the uniform column if present)
      n_het_cols <- cn - 1L
      pir_het    <- pir[, seq_len(n_het_cols), drop = FALSE]   # [n_het x (cn-1)]

      # argmax over het columns → col index in original (unsorted) mvec order
      # Re-sort back through mvec_ord: column k in pir ↔ mvec_ord[k] in sorted
      # We need the dosage for each col in the ORIGINAL (unsorted) pir layout.
      # pir columns correspond to original parms$mvec order.
      # col_to_dosage was built for sorted mvec; map back via mvec_ord.
      orig_col_to_dosage <- integer(n_het_cols)
      orig_col_to_dosage[mvec_ord] <- col_to_dosage   # invert the sort permutation

      best_col <- max.col(pir_het, ties.method = "first")       # per-row argmax
      d_het    <- orig_col_to_dosage[best_col]                  # dosage integer

      # max_prob: probability for the winning het component (across ALL pir columns
      # including uniform — treat uniform as an additional class)
      mp_het <- pir[cbind(seq_len(nrow(pir)), best_col)]        # prob of winning component

      dosage_vec[het_idx] <- d_het
      max_prob[het_idx]   <- mp_het

      # Fill prob_mat columns for het dosages from pir
      # Each col in pir_het → a dosage integer via orig_col_to_dosage
      for (col_k in seq_len(n_het_cols)) {
        d_col <- orig_col_to_dosage[col_k]
        prob_mat[het_idx, d_col + 1L] <- pir_het[, col_k]
      }
    } else {
      # Fallback: assign each het SNP to nearest dosage peak by BAF distance
      for (k in seq_along(het_idx)) {
        baf_k <- baf_vec[het_idx[k]]
        d_k   <- which.min(abs(baf_k - exp_peaks))  # 1-based → dosage 1..(cn-1)
        dosage_vec[het_idx[k]] <- d_k
        max_prob[het_idx[k]]   <- NA_real_
        prob_mat[het_idx[k], d_k + 1L] <- 1.0
      }
    }
  } else if (length(het_idx) == 1L) {
    # Only one het observation — assign by proximity
    baf_k <- baf_vec[het_idx]
    d_k   <- which.min(abs(baf_k - (seq_len(cn - 1L) / cn)))
    dosage_vec[het_idx] <- d_k
    max_prob[het_idx]   <- NA_real_
    prob_mat[het_idx, d_k + 1L] <- 1.0
  }

  out <- data.frame(BAF = baf_vec, dosage = dosage_vec, max_prob = max_prob,
                    prob_mat, check.names = FALSE, stringsAsFactors = FALSE)

  # Pad with NA columns for dosages cn+1..max_cn (layout compatibility)
  if (!is.null(max_cn) && max_cn > cn) {
    extra <- data.frame(matrix(NA_real_, nrow = n_obs, ncol = max_cn - cn))
    names(extra) <- paste0("dosage_", seq(cn + 1L, max_cn))
    out <- cbind(out, extra)
  }

  structure(list(data = out), class = c("BAF_dosage_calls", "list"))
}
