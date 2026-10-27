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
#'   filtered as homozygous. Default is \code{NULL}, in which case it is
#'   derived from \code{qploidy_standardization_object$info["ploidy.standardization"]}.
#' @param het_range Numeric vector of length 2. When \code{geno} is \code{NA}
#'   for a site, the site is kept as heterozygous only if its \code{ratio}
#'   falls within \code{[het_range[1], het_range[2]]}. Default is
#'   \code{c(0.1, 0.9)}.
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
#'   \item For sites with a valid \code{geno} call: sites where
#'     \code{geno == 0} or \code{geno == max_geno} (homozygous by genotype
#'     call) are removed.
#'   \item For sites where \code{geno} is \code{NA}: sites whose
#'     \code{ratio} falls outside \code{het_range} are removed (treated as
#'     homozygous by allele frequency).
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
                           het_range = c(0.1, 0.9),
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

    # Remove homozygous sites by ratio (applies to all sites)
    temp <- temp[!is.na(temp$ratio) & temp$ratio != 0 & temp$ratio != 1, ]

    # For sites with a valid geno call, remove homozygous genotypes.
    # For sites where geno is NA, fall back to ratio-based het filter.
    has_geno <- !is.na(temp$geno)
    keep_geno <- has_geno & temp$geno != 0 & temp$geno != max_geno
    keep_ratio_fallback <- !has_geno &
      temp$ratio >= het_range[1] & temp$ratio <= het_range[2]
    temp <- temp[keep_geno | keep_ratio_fallback, ]

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
#   min_het_frac - minimum fraction of het sites required (below this, or with
#                  fewer than 2 usable het sites, neutral equal log-likelihoods
#                  are returned instead of fitting any model)
#
# Returns a list with:
#   ll_vec   - named numeric vector of log-likelihoods, one per CN in cn_grid
#   prob_vec - softmax probabilities (same length)
#   no_one   - always FALSE; unused placeholder kept for structural parity
#              with compute_baf_likelihoods()'s equivalent field
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
        res <- fn(parms, xi_all, niter, epsilon, trunc, "fixed-2")
        res$loglikelihood
      } else if (base_family == "normal") {
        fn  <- if (add_uniform) nQuack::emstepNU else nQuack::emstepN
        res <- fn(parms, xi_all, niter, epsilon, trunc, "fixed-2")
        res$loglikelihood
      } else if (base_family == "beta-binomial") {
        fn  <- if (add_uniform) nQuack::emstepBBU else nQuack::emstepBB
        res <- fn(parms, xm_all, niter, epsilon, trunc, "fixed-2")
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
# 4. Build xm = cbind(R, allele_count), randomly flipping each site to its ref
#    or alt count (P=0.5) to keep the BAF distribution symmetric for nQuack;
#    nQuack divides col2/col1 internally.
# 5. Derive ploidy name vector from cn_grid (2:6 only); warn for out-of-range.
# 6. Fit quackNormal / quackBeta / quackBetaBinom; select the global minimum
#    BIC across every evaluated (distribution, type, mixture) combination;
#    map the winning mixture name to an integer CN.
# ---------------------------------------------------------------------------
#' @keywords internal
.nquack_estimate_ploidy <- function(d, sample_id, cn_grid,
                                    baf_col          = "baf",
                                    selected_model   = NULL,
                                    nQuack_cores     = 1L,
                                    het_range        = c(0.05, 0.95),
                                    dosage_threshold = 0,
                                    min_depth        = 10L,
                                    fixed_dist       = NULL,
                                    allowed_dists    = NULL,
                                    verbose          = TRUE) {
  if (!requireNamespace("nQuack", quietly = TRUE)) {
    stop("nQuack_EM = TRUE requires the 'nQuack' package. Install with:\n  remotes::install_github(\"mgaynor1/nQuack\")")
  }
  if (utils::packageVersion("nQuack") < "1.0.5") {
    stop("nQuack_EM = TRUE requires nQuack >= 1.0.5 ('glreturn' argument to emstepNA/emstep* was added in 1.0.5). Update with:\n  remotes::install_github(\"mgaynor1/nQuack@updateGL\")")
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

  if (length(xi_het) < 10L)
    stop("nQuack_EM = TRUE: fewer than 10 het sites remain after all filters; cannot fit mixture model.")

  # --- 4. Build xm matrix -----------------------------------------------
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

  # --- 5. Derive ploidy name vector from cn_grid ------------------------
  cn_to_name      <- c("2" = "diploid", "3" = "triploid", "4" = "tetraploid",
                        "5" = "pentaploid", "6" = "hexaploid")
  grid_valid      <- cn_grid[cn_grid >= 2L & cn_grid <= 6L]
  nQuack_mixtures <- unname(cn_to_name[as.character(grid_valid)])
  nQuack_mixtures <- nQuack_mixtures[!is.na(nQuack_mixtures)]
  if (length(nQuack_mixtures) == 0L)
    stop("nQuack_EM = TRUE: cn_grid contains no values in nQuack's supported range (2:6).")
  vmsg("nQuack mixtures: %s", verbose = verbose, level = 1, type = ">>",
       paste(nQuack_mixtures, collapse = ", "))

  # --- 6. Fit quackNormal / quackBeta / quackBetaBinom ------------------
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
  out_all <- do.call(rbind, out_parts)
  out_all[out_all$type == "fixed_2", ]   # <-- keep only fixed-2 rows for BIC selection
  
  }, error = function(e) {
    stop(sprintf("nQuack_EM = TRUE: nQuack model fitting failed: %s", conditionMessage(e)))
  })

  # --- 7. Select best model by BIC across all evaluated models ----------
  # quackit() only ranks models WITHIN each distribution x type group (its
  # source loops over unique(distribution) and unique(type), picking a local
  # winner per group) - it never returns an overall cross-family winner.
  # Taking winnerBIC[1L] from it just returns the first group's local winner
  # (the Normal/fixed row whenever Normal is evaluated), regardless of every
  # other fit. BIC is explicitly designed to be comparable across models of
  # differing complexity, so the correct cross-family rule is the global
  # minimum BIC across every evaluated (distribution, type, mixture)
  # combination for the requested mixtures - this also guarantees ploidy and
  # distribution are read from the exact same winning row.
  candidates <- out_nq[out_nq$mixture %in% nQuack_mixtures & is.finite(out_nq$BIC), ]
  if (nrow(candidates) == 0L)
    stop("nQuack_EM = TRUE: no finite BIC values available across the evaluated nQuack models.")
  best_idx    <- which.min(candidates$BIC)
  winner_name <- candidates$mixture[best_idx]
  winner_dist <- candidates$distribution[best_idx]
  vmsg("nQuack BIC winner: %s (%s)", verbose = verbose, level = 1, type = ">>", winner_name, winner_dist)

  name_to_cn <- c("diploid" = 2L, "triploid" = 3L, "tetraploid" = 4L,
                   "pentaploid" = 5L, "hexaploid" = 6L)
  nq_cn <- unname(name_to_cn[winner_name])
  if (is.na(nq_cn))
    stop(sprintf("nQuack_EM = TRUE: winner '%s' could not be mapped to a CN value.", winner_name))

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
  if (utils::packageVersion("nQuack") < "1.0.5")
    stop("call_BAF_nQuack_dosages requires nQuack >= 1.0.5 ('glreturn' argument to emstep* was added in 1.0.5). Update with:\n  remotes::install_github(\"mgaynor1/nQuack@updateGL\")")

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
        res <- fn(parms, xi_het, niter, epsilon, trunc, "fixed-2", glreturn = TRUE)
        res$pir
      } else if (base_family == "normal") {
        fn  <- if (add_uniform) nQuack::emstepNU else nQuack::emstepN
        res <- fn(parms, xi_het, niter, epsilon, trunc, "fixed-2", glreturn = TRUE)
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
        res <- fn(parms, xm_het, niter, epsilon, trunc, "fixed-2", glreturn = TRUE)
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

      # max_prob: posterior probability of the winning (non-uniform) dosage component
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


#' Estimate copy-number states using the nQuack EM mixture model
#'
#' Self-contained analogue of \code{\link{hmm_estimate_CN}} that always uses
#' \pkg{nQuack} for both the overall sample ploidy estimate and per-window BAF
#' likelihoods/dosage calling (the "full nQuack" path); it does not call
#' \code{hmm_estimate_CN} internally. \code{cn_grid} must contain at least one
#' value in \code{2:6} (nQuack's supported range) unless \code{z_only = TRUE}.
#'
#' @param qploidy_standarize_result An object of class \code{qploidy_standardization} as returned by \code{standardize()}, or \code{NULL}. When \code{NULL}, raw data must be supplied via \code{data} and \code{geno.pos}.
#' @param sample_id Character scalar. Sample identifier to analyze; rows are filtered as \code{SampleName == sample_id}.
#' @param data Optional. A \code{data.frame} with columns \code{MarkerName}, \code{SampleName}, \code{ratio}, and \code{R}. Required when \code{qploidy_standarize_result} is \code{NULL}.
#' @param geno.pos Optional. A \code{data.frame} with columns \code{MarkerName}, \code{Chromosome}, and \code{Position}. Required when \code{qploidy_standarize_result} is \code{NULL}.
#' @param use_values Character vector of length 2 specifying which columns to use as the BAF-like and depth-like signals. When \code{qploidy_standarize_result} is provided, the first element can be \code{"BAF"} (uses column \code{baf}) or \code{"ratio"} (uses column \code{ratio}), and the second element can be \code{"zscore"} (uses column \code{z}) or \code{"R"} (uses column \code{R}). All four combinations are accepted: \code{c("BAF", "zscore")} (default), \code{c("BAF", "R")}, \code{c("ratio", "zscore")}, \code{c("ratio", "R")}. When \code{qploidy_standarize_result} is \code{NULL}, the only accepted value is \code{c("ratio", "R")}.
#' @param chr Optional. Character or integer vector specifying chromosomes to include. If \code{NULL}, all chromosomes are used.
#' @param segment_zscore Logical. If \code{TRUE}, segment z-scores using changepoint detection to define windows instead of fixed SNP counts. Default \code{TRUE}.
#' @param snps_per_window Integer. Number of SNPs per window when windows are created with fixed-size windowing (\code{segment_zscore = FALSE}). Default \code{50}.
#' @param min_snps_per_window Integer. Minimum SNPs required to keep a window. If \code{NULL}, a dynamic value is chosen based on chromosome size. Windows with fewer SNPs are dropped. Default: \code{NULL} (dynamic), or user-specified value.
#' @param cn_grid Integer vector of copy-number states to consider. Must contain at least one value in \code{2:6} (nQuack's supported range) unless \code{z_only = TRUE}; \code{0} is not allowed (use \code{\link{hmm_estimate_CN}} for CN=0/LOH detection). Unlikely values will be discarded during estimation if their z means are not monotonic with ploidy. Default \code{2:6}.
#' @param max_iter Integer. Maximum EM iterations. Default \code{60}.
#' @param het_quantile Numeric. Quantile used to scale BAF emission weight based on heterozygote count. Default \code{0.8}.
#' @param baf_weight Numeric. Overall weight applied to BAF emission (0–1). Default \code{0.5}.
#' @param z_range Numeric. Padding added to min/max z for initial mean estimation. Default \code{NULL}.
#' @param transition_jump Numeric. Diagonal value for transition matrix (probability to stay in same CN state). Default \code{0.9}.
#' @param initial_prob Numeric. Initial probability for the best CN state in the initial state distribution (pi0). Default \code{0.5}. Sets the prior probability for the expected ploidy (or nQuack's estimated ploidy) at the first window; remaining probability is distributed uniformly across other states.
#' @param z_only Logical. If \code{TRUE}, fit the HMM using the z-emission only (ignores BAF and skips nQuack entirely). Default \code{FALSE}.
#' @param baf_only Logical. If \code{TRUE}, fit the HMM using the BAF emission only (ignores z-score). Mutually exclusive with \code{z_only}. Default \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, print progress messages. Default \code{TRUE}.
#' @param exp_ploidy Numeric. Expected ploidy value. If \code{NA} or \code{NULL}, it is estimated by nQuack's global ploidy fit (or, if that fails, set to the middle value of \code{cn_grid}'s nQuack-supported range). Default \code{NA}.
#' @param rm_outliers Logical. If \code{TRUE}, remove outliers from z-scores before HMM fitting. Default \code{TRUE}.
#' @param outlier_alpha Numeric. Alpha threshold for outlier removal in z-scores. Default \code{0.05}.
#' @param selected_model Optional. A list with a \code{best$best_cn} element giving a previously confirmed/trusted ploidy (e.g., the result of a prior call to this function's internal nQuack ploidy estimation). When supplied, nQuack skips the ploidy search and fits distributions directly for that known ploidy. If \code{NULL} (default), nQuack searches for the best-fitting ploidy.
#' @param correct_scale Logical. If \code{TRUE} (default), the BAF log-likelihood is corrected by the number of markers with valid BAF values in each window, so that windows with different numbers of markers contribute equally to the combined emission. Prevents windows with many markers from dominating the HMM via the BAF term alone.
#' @param min_het_frac Numeric in [0,1]. Threshold for the fraction of BAF values in \code{het_range} considered heterozygous. If the observed heterozygous fraction exceeds this value, CN=1 is excluded during nQuack BAF model selection and per-window likelihood computation, as a meaningful proportion of heterozygous loci makes haploid (CN=1) implausible. Default \code{0.01}.
#' @param het_range Numeric vector of length 2. BAF interval used to define heterozygous loci (default \code{c(0.05, 0.95)}). Used by the \code{min_het_frac} filter and to strip homozygous flanking sites before nQuack fitting.
#' @param dosage_threshold Numeric in [0,1]. Minimum posterior probability required for a dosage call to be counted as a heterozygote when computing the BAF emission weight (\code{w_baf}) per window. Markers whose maximum dosage posterior falls below this threshold are excluded from the heterozygote count. Default \code{0.6}.
#' @param z_no_baf_scale Numeric in [0,1]. Floor weight applied to \code{w_baf} when computing \code{CN_reliability} for windows with no BAF support (\code{w_baf == 0}). A value of 0 means z-only calls have \code{CN_reliability = 0}; 0.25 (default) allows z-only calls to reach at most 25\% of \code{post_max}. Does not affect the HMM emission or CN calling. Applies only when \code{z_only = FALSE}.
#' @param hom_z_sigma_inflate Numeric >= 1. Two-tier correction for reference-bias-driven z elevation in windows with no BAF support (\code{w_baf == 0}). (1) During decoding: sigma is multiplied by this factor for \code{w_baf == 0} windows, widening the z-emission. (2) Post-decoding override: for chromosomes where \emph{every} window has \code{w_baf == 0}, the chromosome mean z is compared to the sample-baseline mu; if the deviation is within \code{hom_z_sigma_inflate * sigma} the entire chromosome is overridden to the sample-mode CN (reference bias); larger deviations keep the decoded CN (genuine change). Default \code{1.5}. Set to \code{1} to disable both corrections.
#' @param nQuack_cores Integer. Number of cores to use for nQuack's internal
#'   parallelism (default: 1).
#' @param nQuack_dists Character vector or \code{NULL}. Restricts which nQuack
#'   mixture-model families are evaluated. Valid values are \code{"normal"},
#'   \code{"beta"}, and \code{"beta-binomial"} (plus their \code{"-uniform"}
#'   variants). \code{NULL} (default) tests all families.
#'
#' @return An object of class \code{hmm_CN}. Unlike \code{hmm_estimate_CN()},
#'   \code{params$bins} is always \code{NULL} and \code{params$bw} is always
#'   \code{NA}, since the nQuack route does not use histogram-binned templates.
#'
#' @export
hmm_estimate_CN_nQuack <- function(
  qploidy_standarize_result = NULL,
  sample_id,
  data = NULL,
  geno.pos = NULL,
  use_values = c("BAF", "zscore"),
  chr = NULL,
  segment_zscore = TRUE,
  snps_per_window = 50,
  min_snps_per_window = NULL,
  cn_grid = 2:6,
  max_iter = 60,
  het_quantile = 0.8,
  baf_weight = 0.5,
  z_range = NULL,
  transition_jump = 0.9,
  initial_prob = 0.5,
  z_only = FALSE,
  baf_only = FALSE,
  verbose = TRUE,
  exp_ploidy = NA,
  rm_outliers = TRUE,
  outlier_alpha = 0.05,
  selected_model = NULL,
  correct_scale = TRUE,
  min_het_frac = 0.01,
  het_range = c(0.05, 0.95),
  dosage_threshold = 0,
  z_no_baf_scale = 0.25,
  hom_z_sigma_inflate = 1.5,
  nQuack_cores = 1L,
  nQuack_dists = NULL
) {
  # ---- Input checks ----
  if (any(cn_grid < 0) | any(cn_grid > 10)) stop("All values in cn_grid must be between 0 and 10.")
  if (0 %in% cn_grid)
    stop("cn_grid cannot include 0 (nullisomy/CN=0): nQuack only models CN 2:6. Use hmm_estimate_CN() for CN=0/LOH detection.")
  if (z_only && baf_only) stop("z_only and baf_only cannot both be TRUE.")
  if (!is.character(sample_id) || length(sample_id) != 1 || nchar(sample_id) == 0)
    stop("sample_id must be a non-empty character scalar.")

  # nQuack only models CN 2:6 for BAF; CN outside that range should not be passed unless z_only
  nq_cn_valid <- cn_grid[cn_grid >= 2L & cn_grid <= 6L]
  if (!z_only && length(nq_cn_valid) == 0)
    stop("cn_grid must contain at least one value in 2:6 for nQuack BAF estimation.")

  vmsg("Preparing inputs and applying initial filters", verbose = verbose, level = 0, type = ">>")

  # ---- Data preparation ----
  valid_uv_std <- list(c("BAF", "zscore"), c("BAF", "R"), c("ratio", "zscore"), c("ratio", "R"))
  valid_uv_raw <- c("ratio", "R")
  if (!is.null(qploidy_standarize_result)) {
    if (!is(qploidy_standarize_result, "qploidy_standardization"))
      stop("Input must be a qploidy_standardization object as returned by standardize().")
    if (!any(sapply(valid_uv_std, identical, use_values)))
      stop("When qploidy_standarize_result is provided, use_values must be one of the four valid combinations.")
    baf_col <- if (use_values[1] == "BAF") "baf" else "ratio"
    z_col   <- if (use_values[2] == "zscore") "z" else "R"
    df <- as.data.frame(qploidy_standarize_result$data)
    if (nrow(df) == 0) stop("Input data is empty.")
    d <- filter(df, SampleName == sample_id)
    if (nrow(d) == 0) stop(sprintf("No data found for sample '%s'.", sample_id))
  } else {
    if (!identical(use_values, valid_uv_raw))
      stop("When qploidy_standarize_result is NULL, use_values must be c(\"ratio\", \"R\").")
    if (is.null(data) || is.null(geno.pos))
      stop("When qploidy_standarize_result is NULL, both 'data' and 'geno.pos' must be provided.")
    req_data <- c("MarkerName", "SampleName", "ratio", "R")
    req_gp   <- c("MarkerName", "Chromosome", "Position")
    miss_d  <- setdiff(req_data, names(data))
    miss_gp <- setdiff(req_gp, names(geno.pos))
    if (length(miss_d))  stop(paste("'data' is missing columns:", paste(miss_d, collapse = ", ")))
    if (length(miss_gp)) stop(paste("'geno.pos' is missing columns:", paste(miss_gp, collapse = ", ")))
    baf_col <- "ratio"
    z_col   <- "R"
    d <- data[data$SampleName == sample_id, , drop = FALSE]
    if (nrow(d) == 0) stop(sprintf("No data found for sample '%s'.", sample_id))
    gp <- geno.pos[, req_gp]
    colnames(gp)[colnames(gp) == "Chromosome"] <- "Chr"
    d <- merge(d, gp, by = "MarkerName", all.x = FALSE)
    if (nrow(d) == 0) stop("No markers remain after merging data with geno.pos.")
  }

  # Subset to requested chromosomes
  if (!is.null(chr)) {
    if (is.numeric(chr)) {
      uniq_chr <- unique(d$Chr)
      uniq_chr <- uniq_chr[!is.na(uniq_chr)]
      chrs <- uniq_chr[chr]
    } else {
      chrs <- chr
    }
    chrs <- chrs[!is.na(chrs)]
    if (length(chrs) == 0) stop("No valid chromosomes found for specified 'chr'.")
    d <- filter(d, !is.na(Chr) & Chr %in% chrs)
    if (nrow(d) == 0) stop("No data found for specified chromosomes after filtering.")
  }

  if ("R" %in% names(d)) d$R[d$R == 0] <- NA

  if (rm_outliers) {
    n_na <- sum(is.na(d[[z_col]]))
    if (z_col != "z") {
      if ("z" %in% names(d)) names(d)[names(d) == "z"] <- ".__z_stash__."
      names(d)[names(d) == z_col] <- "z"
    }
    d <- rm_outlier(d, z = TRUE, alpha = outlier_alpha)
    if (z_col != "z") {
      names(d)[names(d) == "z"] <- z_col
      if (".__z_stash__." %in% names(d)) names(d)[names(d) == ".__z_stash__."] <- "z"
    }
    vmsg("Outliers removed from %s scores: %s", verbose = verbose, level = 1, type = ">>",
         z_col, sum(is.na(d[[z_col]])) - n_na)
  }

  rm <- which(is.na(d[[z_col]]))
  if (length(rm) > 0) d <- d[-rm, , drop = FALSE]
  vmsg("Inputs good to go!", verbose = verbose, level = 1, type = ">>")

  # ---- nQuack global ploidy estimation ----
  vmsg("Running nQuack global ploidy estimation", verbose = verbose, level = 0, type = ">>")

  nq_dist_win <- NULL  # nQuack distribution string for per-window calls

  if (!isTRUE(z_only)) {
    baf_candidates <- c("baf", "BAF", "ratio")
    baf_col_nq  <- baf_candidates[baf_candidates %in% colnames(d)][1]

    nq_result <- tryCatch(
      .nquack_estimate_ploidy(
        d                = d,
        sample_id        = sample_id,
        baf_col          = baf_col_nq,
        cn_grid          = nq_cn_valid,
        selected_model   = selected_model,
        nQuack_cores     = nQuack_cores,
        het_range        = het_range,
        dosage_threshold = dosage_threshold,
        allowed_dists    = nQuack_dists,
        verbose          = verbose
      ),
      error = function(e) {
        vmsg("nQuack ploidy estimation failed: %s", verbose = verbose, level = 1,
             type = ">>", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(nq_result) && !is.null(nq_result$cn) && is.finite(nq_result$cn)) {
      nq_cn       <- as.integer(nq_result$cn)
      nq_dist_win <- nq_result$dist  # e.g. "beta-binomial-uniform"

      if (is.null(exp_ploidy) || (length(exp_ploidy) == 1 && is.na(exp_ploidy))) {
        exp_ploidy <- nq_cn
        vmsg("nQuack estimated exp_ploidy = %d (dist = %s)", verbose = verbose,
             level = 1, type = ">>", nq_cn, nq_dist_win)
      }
      vmsg("nQuack: will use dist = %s for per-window BAF likelihoods",
           verbose = verbose, level = 1, type = ">>", nq_dist_win)
    } else {
      vmsg("nQuack ploidy estimation returned no result; exp_ploidy must be provided",
           verbose = verbose, level = 1, type = ">>")
      if (is.null(exp_ploidy) || (length(exp_ploidy) == 1 && is.na(exp_ploidy)))
        stop("nQuack failed and exp_ploidy is not set. Provide exp_ploidy or check nQuack installation.")
    }
  }

  if (is.null(exp_ploidy) || (length(exp_ploidy) == 1 && is.na(exp_ploidy)))
    stop("exp_ploidy could not be estimated and must be provided explicitly (z_only = TRUE skips nQuack's ploidy estimation).")

  vmsg("exp_ploidy = %s", verbose = verbose, level = 1, type = ">>", exp_ploidy)

  # ---- Build windows ----
  vmsg("Defining windows", verbose = verbose, level = 0, type = ">>")
  d[["Position"]] <- as.numeric(d[["Position"]])
  d <- d[order(d[["Chr"]], d[["Position"]]), ]

  if (is.null(min_snps_per_window)) {
    floor_val <- 5L
    frac      <- 0.20
    mk_by_chrom <- d %>%
      group_by(Chr) %>%
      summarize(n = n())
    x <- min(mk_by_chrom$n)
    m <- max(floor_val, floor(x * frac))
    min_snps_per_window <- min(m, floor(x / 2))
  }

  if (segment_zscore) {
    vmsg("Using z-scores changepoint detection to define windows", verbose = verbose, level = 1, type = ">>")
    if (z_col != "z") {
      if ("z" %in% names(d)) names(d)[names(d) == "z"] <- ".__z_stash__."
      names(d)[names(d) == z_col] <- "z"
    }
    d <- add_changepoint_windows(dat = d, minseglen = min_snps_per_window)
    if (z_col != "z") {
      names(d)[names(d) == "z"] <- z_col
      if (".__z_stash__." %in% names(d)) names(d)[names(d) == ".__z_stash__."] <- "z"
    }
  } else {
    vmsg("Using fixed-size windows", verbose = verbose, level = 1, type = ">>")
    d$.__w__ <- with(d, ave(
      seq_len(nrow(d)), d[["Chr"]],
      FUN = function(k) {
        n <- length(k)
        w <- ceiling(seq_along(k) / snps_per_window)
        if (n > snps_per_window && n %% snps_per_window != 0L) {
          last <- max(w); w[w == last] <- last - 1L
        }
        w
      }
    ))
  }
  win_col <- ".__w__"

  win_df <- d %>%
    group_by(Chr, WindowID = .data[[win_col]]) %>%
    summarize(
      Start  = min(Position),
      End    = max(Position),
      n_snps = n(),
      z_mean = mean(.data[[z_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame()

  if (any(is.nan(win_df$z_mean))) warning("Some windows have no z-score values.")
  if (is.null(win_df) || nrow(win_df) == 0) stop("No windows could be formed.")

  o <- order(win_df$Chr, win_df$Start)
  win_df <- win_df[o, ]
  rownames(win_df) <- NULL
  W <- nrow(win_df)

  # Per-window BAF and depth (R) splits
  d_key   <- paste(d[["Chr"]], d[[win_col]], sep = "\r")
  win_key <- paste(win_df$Chr, win_df$WindowID, sep = "\r")

  baf_split <- split(d[[baf_col]], d_key)
  baf_list  <- baf_split[match(win_key, names(baf_split))]
  baf_list  <- lapply(baf_list, function(x) if (is.null(x)) numeric(0) else x)

  # R_list for beta-binomial dosage calling
  if ("R" %in% names(d)) {
    R_split <- split(d[["R"]], d_key)
    R_list  <- R_split[match(win_key, names(R_split))]
    R_list  <- lapply(R_list, function(x) if (is.null(x)) integer(0) else as.integer(x))
  } else {
    R_list <- vector("list", W)
  }

  for (i in seq_len(W)) {
    if (length(baf_list[[i]]) == 0 || all(is.na(baf_list[[i]]))) {
      warning(sprintf(
        "Window %d (Chr %s, WindowID %s) has no valid BAF values.",
        i, win_df$Chr[i], win_df$WindowID[i]
      ))
    }
  }

  z <- win_df$z_mean
  if (length(z) == 0 || all(is.na(z))) stop("No valid z values found in any window.")
  vmsg("Windows ready", verbose = verbose, level = 1, type = ">>")

  keep <- win_df$n_snps >= min_snps_per_window
  if (!any(keep)) stop("All windows dropped by min_snps_per_window filter.")

  # ---- nQuack per-window BAF likelihoods ----
  if (!z_only) {
    vmsg("Generating nQuack BAF likelihoods by window", verbose = verbose, level = 0, type = ">>")

    baf_results <- lapply(seq_len(nrow(win_df)), function(i) {
      .nquack_baf_likelihoods_by_window(
        win_df   = win_df[i, , drop = FALSE],
        d_all    = d,
        baf_col  = baf_col,
        cn_grid  = nq_cn_valid,
        nq_dist  = nq_dist_win,
        het_range    = het_range,
        min_het_frac = min_het_frac
      )
    })

    ll_baf_matrix <- do.call(rbind, lapply(baf_results, function(r) r$ll_vec))

    # Pad zero-likelihood columns for any cn_grid states outside nQuack's 2:6 range
    cur_cns_nq  <- as.integer(sub("CN", "", colnames(ll_baf_matrix)))
    missing_nq  <- setdiff(cn_grid[cn_grid != 0L], cur_cns_nq)
    if (length(missing_nq) > 0L) {
      zero_cols <- matrix(0, nrow = nrow(ll_baf_matrix), ncol = length(missing_nq),
                          dimnames = list(NULL, paste0("CN", missing_nq)))
      ll_baf_matrix <- cbind(ll_baf_matrix, zero_cols)
    }
    cn_grid_baf_cols <- paste0("CN", intersect(cn_grid, as.integer(sub("CN", "", colnames(ll_baf_matrix)))))
    ll_baf_matrix <- ll_baf_matrix[, cn_grid_baf_cols, drop = FALSE]

    # Best per-window ploidy for dosage calling
    ploidies_temp <- as.integer(sub("CN", "", colnames(ll_baf_matrix)[
      apply(ll_baf_matrix, 1, which.max)]))

    # Final column alignment to full cn_grid
    all_cols <- paste0("CN", cn_grid)
    missing_full <- setdiff(all_cols, colnames(ll_baf_matrix))
    if (length(missing_full) > 0L) {
      z_pad <- matrix(0, nrow = nrow(ll_baf_matrix), ncol = length(missing_full),
                      dimnames = list(NULL, missing_full))
      ll_baf_matrix <- cbind(ll_baf_matrix, z_pad)
    }
    ll_baf_matrix <- ll_baf_matrix[, all_cols, drop = FALSE]
    vmsg("nQuack BAF likelihoods generated", verbose = verbose, level = 1, type = ">>")

    # ---- nQuack per-window dosage calling ----
    vmsg("Calling dosages via nQuack", verbose = verbose, level = 1, type = ">>")
    dosages <- mapply(
      function(bv, rv, cn_win) {
        call_BAF_nQuack_dosages(
          baf_vec   = bv,
          cn        = cn_win,
          nq_dist   = nq_dist_win,
          het_range = het_range,
          R_vec     = if (length(rv) == length(bv)) rv else NULL
        )
      },
      baf_list, R_list, ploidies_temp,
      SIMPLIFY = FALSE
    )

    vmsg("Counting heterozygotes per window", verbose = verbose, level = 1, type = ">>")
    n_het_before_thresh <- vapply(seq_along(dosages), function(i) {
      sum(dosages[[i]]$data$dosage != 0 &
          dosages[[i]]$data$dosage != ploidies_temp[i], na.rm = TRUE)
    }, integer(1))

    n_het_window <- vapply(seq_along(dosages), function(i) {
      sum(dosages[[i]]$data$dosage != 0 &
          dosages[[i]]$data$dosage != ploidies_temp[i] &
          dosages[[i]]$data$max_prob > dosage_threshold, na.rm = TRUE)
    }, integer(1))

    prob_filtered_idx <- which(n_het_before_thresh > 0 & n_het_window == 0)
    if (length(prob_filtered_idx) > 0) {
      prob_filtered_labels <- paste0(win_df$Chr[prob_filtered_idx], ":", win_df$WindowID[prob_filtered_idx])
      vmsg(
        paste0("Warning: %d window(s) had all heterozygous markers discarded ",
               "due to low dosage probability (< dosage_threshold = %.2f).\n",
               "  Affected windows (Chr:WindowID): %s"),
        verbose = verbose, level = 2, type = ">>",
        length(prob_filtered_idx), dosage_threshold,
        paste(prob_filtered_labels, collapse = ", ")
      )
    }

    n_baf <- sapply(baf_list, function(x) sum(!is.na(x)))

    HET_N <- n_het_window
    if (any(HET_N > 0)) {
      ref_q <- suppressWarnings(quantile(HET_N[HET_N > 0], het_quantile, na.rm = TRUE))
      if (!is.finite(ref_q) || ref_q <= 0) ref_q <- max(HET_N, na.rm = TRUE)
      if (!is.finite(ref_q) || ref_q <= 0) ref_q <- 1
      w_baf <- sqrt(pmin(1, HET_N / ref_q))
    } else {
      w_baf <- rep(0, W)
    }
    w_baf[!is.finite(w_baf)] <- 0
    if (round(exp_ploidy, 0) == 1) w_baf <- rep(1, W)
    w_baf <- w_baf * baf_weight
  } else {
    n_het_window <- NA
    n_baf  <- rep(0L, W)
    ll_baf_matrix <- matrix(0, nrow = W, ncol = length(cn_grid))
    w_baf  <- rep(0, W)
    if (round(exp_ploidy, 0) == 1) w_baf <- rep(1, W) * baf_weight
  }
  vmsg("BAF weights defined", verbose = verbose, level = 1, type = ">>")

  # ---- HMM setup ----
  vmsg("Setting initial and transition HMM matrices", verbose = verbose, level = 0, type = ">>")
  K         <- length(cn_grid)
  state_ids <- as.character(cn_grid)
  A <- matrix(1e-3, K, K)
  diag(A) <- transition_jump
  A <- A / rowSums(A)

  pi0 <- if (K > 1) rep((1 - initial_prob) / (K - 1), K) else 1
  best_cn  <- as.numeric(exp_ploidy)
  best_idx <- which(as.numeric(cn_grid) == best_cn)
  if (length(best_idx) == 1) {
    pi0[best_idx] <- initial_prob
  } else {
    pi0 <- rep(1 / K, K)
  }

  vmsg("Defining z-score distribution templates", verbose = verbose, level = 0, type = ">>")
  mu <- define_z_limits(
    z = d[[z_col]], z_window = z, cn_grid = cn_grid,
    exp_ploidy = exp_ploidy, z_range = z_range, verbose = verbose
  )

  sig <- sd(z, na.rm = TRUE)
  if (!is.finite(sig) || sig <= 1e-6) sig <- 0.1
  W <- length(z)

  # Sort windows by descending z-deviation for EM stability
  chr_z_dev  <- abs(tapply(z, win_df$Chr, mean, na.rm = TRUE) - mean(z, na.rm = TRUE))
  chr_em_ord <- names(sort(chr_z_dev, decreasing = TRUE))
  em_perm    <- order(match(win_df$Chr, chr_em_ord), win_df$Start)
  z_em       <- z[em_perm]
  ll_baf_em  <- ll_baf_matrix[em_perm, , drop = FALSE]
  n_baf_em   <- n_baf[em_perm]
  w_baf_em   <- w_baf[em_perm]

  vmsg("Initial z-score mean by state: %s", verbose = verbose, level = 2, type = ">>",
       paste0(round(mu, 3), collapse = ", "))
  vmsg("Initial z-score SD: %s", verbose = verbose, level = 2, type = ">>", sig)

  # ---- EM loop ----
  vmsg("Starting EM loop", verbose = verbose, level = 0, type = ">>")
  rm_res <- em_hmm_cn(
    cn_grid, mu, K, state_ids, sig, z_em,
    z_only, ll_baf_em, n_baf_em, w_baf_em,
    correct_scale, A, pi0, W, max_iter, verbose,
    update_pi0 = FALSE, baf_only = baf_only
  )
  list2env(rm_res, envir = environment())

  # Monotonicity pruning
  idx_mono <- 0L
  mu_tol   <- 1e-6
  while (any(abs(mu - sort(mu)) > mu_tol) && idx_mono < 10L) {
    idx_mono <- idx_mono + 1L
    exp_idx  <- which(as.numeric(names(mu)) == as.numeric(exp_ploidy))
    mu_exp   <- mu[exp_idx]
    lower_ploidies <- which(as.numeric(names(mu)) < as.numeric(exp_ploidy))
    keep_lower <- lower_ploidies[order(-as.numeric(names(mu)[lower_ploidies]))]
    if (length(keep_lower) > 0) {
      last_mu    <- mu_exp; valid_lower <- c()
      for (ii in keep_lower) {
        if (mu[ii] < last_mu) { valid_lower <- c(valid_lower, ii); last_mu <- mu[ii] }
      }
      lower_idx <- valid_lower
    } else { lower_idx <- integer(0) }
    higher_ploidies <- which(as.numeric(names(mu)) > as.numeric(exp_ploidy))
    keep_higher <- higher_ploidies[order(as.numeric(names(mu)[higher_ploidies]))]
    if (length(keep_higher) > 0) {
      last_mu     <- mu_exp; valid_higher <- c()
      for (ii in keep_higher) {
        if (mu[ii] > last_mu) { valid_higher <- c(valid_higher, ii); last_mu <- mu[ii] }
      }
      higher_idx <- valid_higher
    } else { higher_idx <- integer(0) }
    valid_idx <- sort(c(lower_idx, exp_idx, higher_idx))
    cn_grid   <- cn_grid[valid_idx]
    mu <- define_z_limits(
      z = d[[z_col]], z_window = z, cn_grid = cn_grid,
      exp_ploidy = exp_ploidy, z_range_out = FALSE, verbose = verbose
    )
    K             <- length(cn_grid)
    ll_baf_matrix <- ll_baf_matrix[, valid_idx, drop = FALSE]
    ll_baf_em     <- ll_baf_matrix[em_perm, , drop = FALSE]
    pi0       <- pi0[valid_idx]
    state_ids <- state_ids[valid_idx]
    A         <- A[valid_idx, valid_idx]
    rm_res <- em_hmm_cn(
      cn_grid, mu, K, state_ids, sig, z_em,
      z_only, as.matrix(ll_baf_em), n_baf_em, w_baf_em,
      correct_scale, as.matrix(A), pi0, W, max_iter, verbose,
      update_pi0 = FALSE, baf_only = baf_only
    )
    list2env(rm_res, envir = environment())
  }

  # Collapsed-state pruning
  collapse_frac <- 0.10
  if (K > 1) {
    mu_sorted  <- sort(mu)
    gaps       <- diff(mu_sorted)
    keep_names <- names(mu_sorted)
    drop_mask  <- rep(FALSE, K)
    exp_name   <- as.character(round(exp_ploidy))
    for (gi in seq_along(gaps)) {
      if (!drop_mask[gi] && abs(gaps[gi]) < collapse_frac * sig) {
        candidate <- keep_names[gi + 1]
        if (candidate == exp_name) next
        drop_mask[gi + 1] <- TRUE
        vmsg("Collapsed CN state dropped: CN%s", verbose = verbose, level = 2, type = ">>", candidate)
      }
    }
    if (any(drop_mask)) {
      keep_state_names <- keep_names[!drop_mask]
      keep_state_idx   <- which(names(mu) %in% keep_state_names)
      mu            <- mu[keep_state_idx]
      cn_grid       <- cn_grid[keep_state_idx]
      K             <- length(cn_grid)
      state_ids     <- state_ids[keep_state_idx]
      pi0           <- pi0[keep_state_idx] / sum(pi0[keep_state_idx])
      A             <- A[keep_state_idx, keep_state_idx]
      ll_baf_matrix <- ll_baf_matrix[, keep_state_idx, drop = FALSE]
    }
  }
  vmsg("EM loop complete", verbose = verbose, level = 1, type = ">>")

  # ---- Decoding ----
  w_baf_dec <- w_baf
  if (!z_only && 0 %in% cn_grid && "0" %in% names(mu) && "1" %in% names(mu)) {
    cn0_mid <- mu["0"] + 0.5 * (mu["1"] - mu["0"])
    w_baf_dec[z <= cn0_mid] <- 0
  }
  sig_dec <- rep(sig, W)
  if (!z_only && !baf_only && hom_z_sigma_inflate > 1)
    sig_dec <- ifelse(w_baf_dec == 0, sig * hom_z_sigma_inflate, sig)
  no_baf <- !z_only & (n_baf == 0)
  if (!z_only && any(no_baf)) w_baf_dec[no_baf] <- 0

  ll_em <- matrix(NA_real_, W, K, dimnames = list(NULL, state_ids))
  for (k in seq_len(K)) {
    c_k <- cn_grid[k]
    llz <- dnorm(z, mean = mu[as.character(c_k)], sd = sig_dec, log = TRUE)
    llz[is.nan(llz)] <- 0
    if (z_only) {
      ll_em[, k] <- llz
    } else if (baf_only) {
      llb <- ll_baf_matrix[, k]
      if (correct_scale) llb <- llb / n_baf
      llb[no_baf] <- 0
      ll_em[, k] <- llb
    } else {
      llb <- ll_baf_matrix[, k]
      if (correct_scale) llb <- llb / n_baf
      llb[no_baf] <- 0
      ll_em[, k] <- llz + w_baf_dec * llb
    }
  }

  llz_mat <- matrix(NA_real_, W, K)
  llb_mat <- matrix(NA_real_, W, K)
  if (!z_only && K > 1) {
    for (k in seq_len(K)) {
      c_k <- cn_grid[k]
      llz_mat[, k] <- dnorm(z, mean = mu[as.character(c_k)], sd = sig_dec, log = TRUE)
      llb_k <- ll_baf_matrix[, k]
      if (correct_scale) llb_k <- llb_k / n_baf
      llb_k[no_baf] <- 0
      llb_mat[, k] <- llb_k
    }
  }

  # Two-signal consistency gating
  if (!z_only && !baf_only && K > 1) {
    baf_gate <- 0.1
    gate_log <- character(0)
    for (chr_i in unique(win_df$Chr)) {
      idx         <- which(win_df$Chr == chr_i)
      mean_w_baf  <- mean(w_baf_dec[idx], na.rm = TRUE)
      if (mean_w_baf < baf_gate) next
      chr_z_sum   <- colSums(llz_mat[idx, , drop = FALSE], na.rm = TRUE)
      chr_baf_sum <- colSums(llb_mat[idx, , drop = FALSE], na.rm = TRUE)
      best_z_chr  <- which.max(chr_z_sum)
      best_baf_chr <- which.max(chr_baf_sum)
      if (best_z_chr == best_baf_chr) {
        ll_em[idx, best_z_chr] <- ll_em[idx, best_z_chr] + log(3) * mean_w_baf
        gate_log <- c(gate_log, sprintf("%s:agree(CN%d)", chr_i, cn_grid[best_z_chr]))
      } else {
        ll_em[idx, ] <- llz_mat[idx, , drop = FALSE]
        gate_log <- c(gate_log, sprintf("%s:disagree(z=CN%d,baf=CN%d)",
          chr_i, cn_grid[best_z_chr], cn_grid[best_baf_chr]))
      }
    }
    message(sprintf("[GATE] W=%d K=%d baf_gate=%.2f | %s",
      W, K, baf_gate, if (length(gate_log)) paste(gate_log, collapse = " ") else "none"))
  }

  vmsg("Decoding CN path", verbose = verbose, level = 0, type = ">>")
  A_dec <- matrix(1e-3, K, K)
  diag(A_dec) <- transition_jump
  A_dec <- A_dec / rowSums(A_dec)

  gamma_dec <- matrix(0, W, K, dimnames = list(NULL, colnames(ll_em)))
  for (chr_i in unique(win_df$Chr)) {
    idx <- which(win_df$Chr == chr_i)
    logpi0_chr <- NULL
    if (!z_only) {
      mean_w_baf_chr <- mean(w_baf_dec[idx], na.rm = TRUE)
      if (is.nan(mean_w_baf_chr)) mean_w_baf_chr <- 0
      chr_z_sum <- colSums(llz_mat[idx, , drop = FALSE], na.rm = TRUE)
      chr_z_sum <- chr_z_sum - max(chr_z_sum)
      lse_z     <- logsumexp(chr_z_sum)
      logpi0_z  <- if (is.finite(lse_z)) chr_z_sum - lse_z else log(pi0 / sum(pi0))
      chr_baf_sum <- colSums(ll_baf_matrix[idx, , drop = FALSE], na.rm = TRUE)
      chr_baf_sum <- chr_baf_sum - max(chr_baf_sum)
      lse_baf     <- logsumexp(chr_baf_sum)
      logpi0_baf  <- if (is.finite(lse_baf)) chr_baf_sum - lse_baf else log(pi0 / sum(pi0))
      logpi0_chr  <- log(
        (1 - mean_w_baf_chr) * exp(logpi0_z) + mean_w_baf_chr * exp(logpi0_baf)
      )
      max_lp     <- max(logpi0_chr)
      logpi0_chr <- if (is.finite(max_lp)) logpi0_chr - max_lp else log(pi0 / sum(pi0))
    } else {
      logpi0_chr <- log(pi0 / sum(pi0))
    }
    gamma_dec[idx, ] <- fb_smooth(ll_em[idx, , drop = FALSE], log(A_dec), logpi0_chr)
  }
  if (!is.matrix(gamma_dec) || ncol(gamma_dec) == 0L)
    stop("Internal error: gamma_dec is not a valid matrix after fb_smooth.")

  safe_which_max <- function(x) { idx <- which.max(x); if (length(idx) == 0L) 1L else idx }
  cn_call  <- cn_grid[vapply(seq_len(nrow(gamma_dec)), function(i) safe_which_max(gamma_dec[i, ]), integer(1L))]
  post_max <- apply(gamma_dec, 1, max)
  cn0_idx  <- cn_call == 0
  cn_call[cn0_idx]  <- NA_integer_
  post_max[cn0_idx] <- NA_real_
  z[cn0_idx]        <- NA_real_
  if (!z_only) w_baf[cn0_idx] <- NA_real_

  # Post-processing: override weak z-only calls back to sample mode
  if (!z_only && !baf_only && hom_z_sigma_inflate > 1) {
    baf_gate_post   <- 0.1
    het_idx         <- w_baf > 0
    mode_cn         <- as.integer(if (any(het_idx)) mode(cn_call[het_idx]) else round(exp_ploidy))
    mode_col        <- which(cn_grid == mode_cn)
    baseline_mu_val <- if (length(mode_col) == 1) mu[mode_col] else mu[which.min(abs(cn_grid - mode_cn))]
    threshold       <- sig * hom_z_sigma_inflate
    for (chr_i in unique(win_df$Chr)) {
      idx             <- which(win_df$Chr == chr_i)
      mean_w_baf_chr  <- mean(w_baf[idx], na.rm = TRUE)
      chr_z_mean      <- mean(z[idx], na.rm = TRUE)
      chr_cn_mode     <- as.integer(mode(cn_call[idx]))
      if (!is.finite(chr_z_mean)) next
      if (chr_cn_mode == mode_cn) next
      if (abs(chr_z_mean - baseline_mu_val) >= threshold) next
      cn_call[idx] <- mode_cn
      if (length(mode_col) == 1) post_max[idx] <- gamma_dec[idx, mode_col]
      message(sprintf("[POST-GATE] %s: CN%d->CN%d (mean_w_baf=%.3f, z_dev=%.3f < %.3f*sig)",
        chr_i, chr_cn_mode, mode_cn, mean_w_baf_chr,
        abs(chr_z_mean - baseline_mu_val), hom_z_sigma_inflate))
    }
  }

  # CN_reliability
  w_baf_use  <- if (z_only) rep(0, W) else w_baf
  w_baf_norm <- pmax(z_no_baf_scale, pmin(1, w_baf_use / baf_weight))

  baf_entropy_bins <- 10L
  baf_breaks <- seq(0, 1, length.out = baf_entropy_bins + 1L)
  baf_info <- vapply(seq_len(W), function(i) {
    if (isTRUE(is.na(w_baf_use[i])) || isTRUE(w_baf_use[i] == 0)) return(z_no_baf_scale)
    bv <- baf_list[[i]]
    bv <- pmax(0, pmin(1, bv[is.finite(bv)]))
    if (length(bv) == 0L) return(z_no_baf_scale)
    counts <- tabulate(findInterval(bv, baf_breaks, rightmost.closed = TRUE),
                       nbins = baf_entropy_bins)
    s <- sum(counts)
    if (s == 0L) return(z_no_baf_scale)
    p <- counts[counts > 0L] / s
    H <- -sum(p * log(p))
    H_max <- log(baf_entropy_bins)
    max(0, min(1, 1 - H / H_max))
  }, numeric(1L))

  if (!z_only && !baf_only && K > 1 && exists("llz_mat") && exists("llb_mat")) {
    softmax_row <- function(x) { x <- x - max(x); e <- exp(x); e / sum(e) }
    zbaf_agree <- vapply(seq_len(W), function(i) {
      if (isTRUE(is.na(w_baf_use[i])) || isTRUE(w_baf_use[i] == 0)) return(1)
      pz <- softmax_row(llz_mat[i, ])
      pb <- softmax_row(llb_mat[i, ])
      if (!all(is.finite(pz)) || !all(is.finite(pb))) return(1)
      1 - 0.5 * sum(abs(pz - pb))
    }, numeric(1L))
  } else {
    zbaf_agree <- rep(1, W)
  }

  CN_reliability <- post_max * baf_info * zbaf_agree * w_baf_norm
  CN_reliability[!is.finite(CN_reliability)] <- NA_real_
  vmsg("CN path decoded", verbose = verbose, level = 1, type = ">>")

  # ---- Assemble output ----
  post_df <- as.data.frame(gamma_dec)
  names(post_df) <- paste0("post_CN", cn_grid)
  result <- cbind(
    data.frame(
      Sample         = sample_id,
      Chr            = win_df$Chr,
      WindowID       = win_df$WindowID,
      Start          = win_df$Start,
      End            = win_df$End,
      n_snps         = win_df$n_snps,
      n_het          = n_het_window,
      z              = z,
      w_baf          = if (z_only) 0 else w_baf,
      CN_call        = cn_call,
      CN_reliability = CN_reliability,
      post_max       = post_max,
      stringsAsFactors = FALSE
    ),
    post_df
  )

  params <- list(
    cn_grid             = cn_grid,
    distribution        = if (!is.null(nq_dist_win)) nq_dist_win else "nQuack",
    mu                  = if (baf_only) NA else mu,
    sigma               = if (baf_only) NA else sig,
    A                   = A,
    pi0                 = pi0,
    bins                = NULL,
    bw                  = NA_real_,
    loglik              = ll_hist[length(ll_hist)],
    z_range             = z_range,
    het_quantile        = het_quantile,
    baf_weight          = baf_weight,
    transition_jump     = transition_jump,
    z_only              = z_only,
    baf_only            = baf_only,
    exp_ploidy          = exp_ploidy,
    rm_outliers         = rm_outliers,
    outlier_alpha       = outlier_alpha,
    segment_zscore      = segment_zscore,
    snps_per_window     = snps_per_window,
    min_snps_per_window = min_snps_per_window,
    nQuack_dists        = nQuack_dists
  )

  if (!is.null(result) && !is.null(d)) {
    map_df <- result[, c("Chr", "WindowID", "w_baf", "CN_call", "post_max", "CN_reliability")]
    names(map_df)[names(map_df) == "WindowID"] <- ".__w__"
    d <- merge(d, map_df, by = c("Chr", ".__w__"), all.x = TRUE, sort = FALSE)
    d <- d[order(match(seq_len(nrow(d)), as.integer(rownames(d)))), ]
    rownames(d) <- NULL
  }

  vmsg("Done!", verbose = verbose, level = 1, type = ">>")
  return(structure(list(by_window = result, by_marker = d, params = params), class = "hmm_CN"))
}