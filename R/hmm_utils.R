##' @keywords internal
# Suppress global variable warnings for non-standard evaluation in dplyr/ggplot2
globalVariables(c(
  "Sample", "Chr", "Start", "End", "CN_call", "prob_call", "w_baf", "Mid",
  "n_vlines", "n_regions", "region_id"
))

.check_max_cn_dev <- function(max_cn_dev) {
  if (!is.null(max_cn_dev) &&
      (!is.numeric(max_cn_dev) || length(max_cn_dev) != 1 || is.na(max_cn_dev) ||
       max_cn_dev < 0 || max_cn_dev != round(max_cn_dev)))
    stop("`max_cn_dev` must be NULL or a single non-negative integer.")
  invisible(NULL)
}

.check_sex_depth_args <- function(sex_min_depth, sex_min_markers) {
  if (!is.numeric(sex_min_depth) || length(sex_min_depth) != 1 || is.na(sex_min_depth) ||
      sex_min_depth < 0)
    stop("`sex_min_depth` must be a single non-negative number.")
  if (!is.numeric(sex_min_markers) || length(sex_min_markers) != 1 || is.na(sex_min_markers) ||
      sex_min_markers < 1 || sex_min_markers != round(sex_min_markers))
    stop("`sex_min_markers` must be a single positive integer.")
  invisible(NULL)
}

# Presence of the heterogametic sex chromosome (Y for XY, W for ZW) from depth.
# Y/W-specific markers have no reference clusters in single-sex standardizations,
# so their CN cannot be called; presence is scored from raw depth (R) instead.
# `d` must be the sample's data before any marker filtering. Returns NULL when
# not applicable (no XY/ZW system, no Y/W group in sex_chr, or no R column).
.sex_het_depth <- function(d, sex_system, sex_chr, sex_min_depth = 6, sex_min_markers = 2) {
  het <- switch(sex_system, XY = "Y", ZW = "W", NULL)
  if (is.null(het) || is.null(sex_chr[[het]]) || !"R" %in% names(d)) return(NULL)
  on_het <- !is.na(d$Chr) & as.character(d$Chr) %in% sex_chr[[het]]
  r <- d$R[on_het]
  r[is.na(r)] <- 0  # markers with missing depth count as uncovered
  depth <- if (any(on_het)) tapply(r, d$MarkerName[on_het], max) else numeric(0)
  n_markers <- length(depth)
  n_covered <- sum(depth > sex_min_depth)
  list(
    chr         = het,
    n_covered   = n_covered,
    n_markers   = n_markers,
    score       = sprintf("%d/%d", n_covered, n_markers),
    present     = if (n_markers == 0) NA else n_covered >= sex_min_markers,
    min_depth   = sex_min_depth,
    min_markers = sex_min_markers
  )
}

# Keep only CN states within max_cn_dev of exp_ploidy (CN 0 kept if keep_zero)
.restrict_cn_grid <- function(cn_grid, exp_ploidy, max_cn_dev, keep_zero = FALSE, verbose = TRUE) {
  if (is.null(max_cn_dev)) return(cn_grid)
  keep <- abs(cn_grid - as.numeric(exp_ploidy)) <= max_cn_dev | (keep_zero & cn_grid == 0)
  new_grid <- cn_grid[keep]
  if (length(new_grid) == 0)
    stop(sprintf("No cn_grid values within max_cn_dev = %d of exp_ploidy = %s.", as.integer(max_cn_dev), exp_ploidy))
  if (length(new_grid) < length(cn_grid))
    vmsg("cn_grid restricted by max_cn_dev = %d around exp_ploidy = %s: %s",
         verbose = verbose, level = 1, type = ">>",
         as.integer(max_cn_dev), exp_ploidy, paste(new_grid, collapse = ", "))
  new_grid
}

# Resolve min_snps_per_window into one threshold per chromosome in `d`.
# Accepts NULL (dynamic), a single value, or a vector named by chromosome
# (unnamed chromosomes fall back to the dynamic default).
.resolve_min_snps <- function(min_snps_per_window, d, chr_subset = NULL, verbose = TRUE) {
  chrs  <- unique(as.character(d$Chr))
  n_chr <- table(as.character(d$Chr))[chrs]

  x <- min(n_chr)
  dyn <- as.integer(min(max(5L, floor(x * 0.20)), floor(x / 2)))

  if (is.null(min_snps_per_window)) {
    vmsg("min_snps_per_window not provided, using: %d (based on the smallest chromosome, %d markers)",
         verbose = verbose, level = 1, type = ">>", dyn, as.integer(x))
    return(stats::setNames(rep(dyn, length(chrs)), chrs))
  }

  if (!is.numeric(min_snps_per_window) || any(is.na(min_snps_per_window)) ||
      any(min_snps_per_window < 1))
    stop("`min_snps_per_window` must be NULL, a positive number, or a positive numeric vector named by chromosome.")

  nms <- names(min_snps_per_window)
  if (is.null(nms)) {
    if (length(min_snps_per_window) != 1)
      stop("An unnamed `min_snps_per_window` must be a single value; name the elements by chromosome to set it per chromosome.")
    return(stats::setNames(rep(as.integer(min_snps_per_window), length(chrs)), chrs))
  }
  if (any(nms == "") || anyDuplicated(nms))
    stop("All elements of `min_snps_per_window` must have unique chromosome names.")

  unknown <- setdiff(nms, chrs)
  if (length(unknown) > 0 && is.null(chr_subset))
    warning("`min_snps_per_window` names not found in the data: ",
            paste(unknown, collapse = ", "), call. = FALSE)

  out <- stats::setNames(rep(dyn, length(chrs)), chrs)
  given <- intersect(nms, chrs)
  out[given] <- as.integer(min_snps_per_window[given])
  default_chrs <- setdiff(chrs, given)
  if (length(default_chrs) > 0)
    vmsg("min_snps_per_window not provided for %s; using: %d (based on the smallest chromosome, %d markers)",
         verbose = verbose, level = 1, type = ">>",
         paste(default_chrs, collapse = ", "), dyn, as.integer(x))
  out
}

# Compact form for params/printing: a single value when all chromosomes share it
.compact_min_snps <- function(v) {
  if (length(unique(v)) == 1) unname(v[1]) else v
}

.format_min_snps <- function(v, short = FALSE, max_exceptions = 2L) {
  if (is.null(v)) return("NA")
  if (length(v) == 1) return(as.character(v))
  if (!short) return(paste(names(v), v, sep = "=", collapse = ", "))
  # Most common value, listing only the few chromosomes that differ
  common <- as.integer(names(which.max(table(v))))
  exc <- v[v != common]
  if (length(exc) <= max_exceptions)
    return(paste0(common, " (", paste(names(exc), exc, sep = "=", collapse = ", "), ")"))
  paste0(min(v), "-", max(v), " (by chr)")
}

#' Numerically stable log-sum-exp
#'
#' Computes the logarithm of the sum of exponentials, \eqn{\log\left(\sum_i e^{x_i}\right)}, in a numerically stable way by subtracting the maximum element before exponentiating. This avoids overflow/underflow for large or small values in \code{x}.
#'
#' @param x Numeric vector. Values to exponentiate and sum.
#'
#' @return Numeric scalar. The stable computation of \code{log(sum(exp(x)))}.
#'
#' @details
#' This function uses the identity \eqn{\log\sum_i e^{x_i} = m + \log\sum_i e^{x_i - m}}, where \eqn{m=\max_i x_i}, to ensure numerical stability.
#'
#'
#' @keywords internal
logsumexp <- function(x) {
  m <- max(x); m + log(sum(exp(x - m)))
}


#' Viterbi decoding for a first-order HMM in log-space
#'
#' Computes the most likely state sequence (Viterbi path) for a Hidden Markov Model given per-position emission log-likelihoods and log transition/prior probabilities. Operates in log-space for numerical stability.
#'
#' @param ll_em Numeric matrix (W x K). Emission log-likelihoods, where \code{W} is the number of positions/windows and \code{K} is the number of states. Column order must match \code{logA}.
#' @param logA Numeric matrix (K x K). Log transition probabilities, where \code{logA[i, j]} is \eqn{\log p(s_t=j | s_{t-1}=i)}. Each row should represent a valid log-probability distribution.
#' @param logpi0 Numeric vector of length \code{K}. Initial state log-probabilities at position 1.
#'
#' @return Integer vector of length \code{W}. 1-based indices of the most likely states at each position (Viterbi path).
#'
#' @details
#' Stores backpointers to reconstruct the optimal path after dynamic programming. All calculations are performed in log-space for stability.
#'
#' @examples
#' set.seed(1)
#' W <- 5; K <- 3
#' ll_em <- matrix(rnorm(W*K), W, K)
#' A <- matrix(1/K, K, K); diag(A) <- 0.8; A <- A / rowSums(A)
#' logA <- log(A)
#' logpi0 <- log(rep(1/K, K))
#' vpath <- viterbi(ll_em, logA, logpi0)
#' vpath
#'
#' @export
viterbi <- function(ll_em, logA, logpi0) {
  W <- nrow(ll_em); K <- ncol(ll_em)
  delta <- matrix(-Inf, W, K)
  psi <- matrix(NA_integer_, W, K)
  delta[1, ] <- logpi0 + ll_em[1, ]
  for (i in 2:W) {
    for (k in 1:K) {
      tmp <- delta[i-1, ] + logA[,k]
      psi[i, k] <- which.max(tmp)
      delta[i, k] <- ll_em[i, k] + max(tmp)
    }
  }
  path <- integer(W)
  path[W] <- which.max(delta[W, ])
  for (i in (W-1):1) path[i] <- psi[i+1, path[i+1]]
  path
}

#' Smoothed posteriors via forward-backward
#'
#' Runs one forward-backward pass on the converged emission matrix.  If
#' \code{logpi0} is supplied the forward pass is initialised from that
#' distribution; otherwise a uniform distribution is used.  Passing a
#' BAF-derived \code{logpi0} prevents a noisy first-window z-score from
#' anchoring the forward pass in the wrong CN state for the whole chromosome.
#'
#' @param ll_em  Numeric matrix (W x K). Final emission log-likelihoods from EM.
#' @param logA   Numeric matrix (K x K). Converged log transition probabilities.
#' @param logpi0 Numeric vector of length K, or NULL for uniform. Log initial
#'   state probabilities for the forward pass.
#'
#' @return Numeric matrix (W x K). Row-normalised posterior probabilities.
#'
#' @keywords internal
#' @noRd
fb_smooth <- function(ll_em, logA, logpi0 = NULL) {
  W <- nrow(ll_em); K <- ncol(ll_em)
  if (is.null(logpi0)) logpi0 <- rep(-log(K), K)

  if (W == 1L) {
    log_e1 <- logpi0 + ll_em[1, ]
    log_gamma <- matrix(log_e1 - logsumexp(log_e1), nrow = 1)
    return(exp(log_gamma))
  }

  log_alpha <- matrix(-Inf, W, K)
  log_alpha[1, ] <- logpi0 + ll_em[1, ]
  for (i in 2:W) {
    for (k in 1:K) {
      log_alpha[i, k] <- ll_em[i, k] + logsumexp(log_alpha[i - 1, ] + logA[, k])
    }
  }

  log_beta <- matrix(0, W, K)
  for (i in (W - 1):1) {
    for (k in 1:K) {
      log_beta[i, k] <- logsumexp(logA[k, ] + ll_em[i + 1, ] + log_beta[i + 1, ])
    }
  }

  log_gamma <- log_alpha + log_beta
  row_lse <- apply(log_gamma, 1, logsumexp)
  # If a row's logsumexp is non-finite (all -Inf or NaN), the posteriors for that
  # window are undefined.  Fall back to uniform so the matrix stays numeric.
  degenerate <- !is.finite(row_lse)
  if (any(degenerate)) row_lse[degenerate] <- 0  # exp(0)=1 -> uniform after subtraction
  log_gamma <- sweep(log_gamma, 1, row_lse, "-")
  # Force degenerate rows to uniform (all -Inf -> -log(K) after sweep is still -Inf;
  # handle by writing log(1/K) directly).
  if (any(degenerate)) log_gamma[degenerate, ] <- -log(K)

  exp(log_gamma)
}

#' Bidirectional Viterbi decoding (max-product forward-backward)
#'
#' Combines a standard forward max-pass (delta) with a backward max-pass (beta)
#' so that every window, including the first, is scored using evidence from the
#' full observation chain.  At each position the state is chosen as
#' \code{argmax_k delta[t,k] + beta[t,k]}.  The forward pass uses a uniform
#' initial distribution so window 1 is treated symmetrically with all others
#' (the backward pass already starts uniformly at window W).  logpi0 is accepted
#' but intentionally not used in the forward initialisation.
#'
#' @param ll_em Numeric matrix (W x K). Emission log-likelihoods.
#' @param logA  Numeric matrix (K x K). Log transition probabilities (row = from, col = to).
#' @param logpi0 Numeric vector of length K. Accepted for interface compatibility; not used.
#'
#' @return Integer vector of length W. 1-based state indices.
#'
#' @keywords internal
#' @noRd
viterbi_bidi <- function(ll_em, logA, logpi0) {
  W <- nrow(ll_em); K <- ncol(ll_em)

  if (W == 1L) return(which.max(ll_em[1, ]))

  # Forward max pass -- uniform init so window 1 is not biased by pi0
  delta <- matrix(-Inf, W, K)
  delta[1, ] <- ll_em[1, ]
  for (i in 2:W) {
    for (k in 1:K) {
      delta[i, k] <- ll_em[i, k] + max(delta[i - 1, ] + logA[, k])
    }
  }

  # Backward max pass: beta[t,k] = max log-prob of obs t+1..W given state k at t
  beta <- matrix(0, W, K)          # beta[W, ] = log(1) = 0
  for (i in (W - 1):1) {
    for (k in 1:K) {
      beta[i, k] <- max(logA[k, ] + ll_em[i + 1, ] + beta[i + 1, ])
    }
  }

  # Decode: each position independently picks the globally best-supported state
  apply(delta + beta, 1, which.max)
}

#' Internal worker for parallel HMM CN estimation
#'
#' Runs hmm_estimate_CN for a single sample within a parallel loop, forwarding arguments.
#' Returns the result or NULL on error (with a warning).
#'
#' @param sid Character. Sample identifier.
#' @param obj Standardized input object for copy-number estimation (usually of class 'qploidy_standardization').
#' @param dots Named list of additional arguments to pass to hmm_estimate_CN.
#'
#' @return List as returned by hmm_estimate_CN, or NULL if an error occurs.
#'
#' @keywords internal
#' @noRd
worker <- function(sid, obj, dots, data = NULL, geno.pos = NULL, use_values = c("BAF", "zscore")) {
  # Route to hmm_estimate_CN_nQuack when nQuack-specific params are present in dots
  nq_params <- c("nQuack_dists", "nQuack_type", "nQuack_uniform")
  fn <- if (any(nq_params %in% names(dots))) hmm_estimate_CN_nQuack else hmm_estimate_CN
  collected_warnings <- character(0)
  result <- withCallingHandlers(
    tryCatch({
      if (!is.null(obj)) {
        do.call(fn,
                c(list(qploidy_standarize_result = obj, sample_id = sid, use_values = use_values), dots))
      } else {
        do.call(fn,
                c(list(qploidy_standarize_result = NULL, sample_id = sid, data = data, geno.pos = geno.pos, use_values = use_values), dots))
      }
    }, error = function(e) {
      collected_warnings <<- c(
        collected_warnings,
        sprintf("Sample '%s' failed: %s", sid, conditionMessage(e))
      )
      NULL
    }),
    warning = function(w) {
      collected_warnings <<- c(
        collected_warnings,
        sprintf("[Sample '%s'] %s", sid, conditionMessage(w))
      )
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warnings = collected_warnings)
}


##' Initialize monotonic z-score means for HMM ploidy states
##'
##' Computes a monotonic ramp of z-score means (mu) for each ploidy state, centered on the expected ploidy, with padding to avoid boundary collapse.
##' Used for initializing HMM emission parameters in ploidy estimation.
##'
##' @param z Numeric vector. Window-level z-scores.
##' @param cn_grid Integer vector. Copy-number states to consider.
##' @param exp_ploidy Numeric. Expected ploidy value (center of ramp).
##' @param z_range Numeric. Padding added to min/max z for ramp initialization. If NULL, estimated from z interquartile range.
##' @param verbose Logical. If TRUE, prints estimated z_range. Default FALSE.
##' @param z_range_out Logical. If TRUE, expand range (min - z_range, max + z_range). If FALSE, reduce range (min + z_range, max - z_range).
##'
##' @return Named numeric vector of z-score means (mu) for each ploidy state in cn_grid.
##'
##' @details
##' The ramp is constructed as: mu_c = z_mean + step * (c - exp_ploidy), where step is chosen so the lowest and highest ploidy states fit within the padded z range.
##' Ensures monotonic initialization for EM fitting.
##' @keywords internal
##' @export
define_z_limits <- function(z, z_window, cn_grid, exp_ploidy, z_range = NULL, verbose = FALSE, z_range_out = TRUE) {
  state_ids <- as.character(cn_grid)

  if (is.null(z_range) || (length(z_range) == 1 && is.na(z_range))) {
    z_range <- (1/length(cn_grid)) * (as.numeric(quantile(z, probs = 0.75)) - as.numeric(quantile(z, probs = 0.25)))
    vmsg("Estimated z_range from data: %f", verbose = verbose, level = 2, type = ">>", z_range)
  }
  z_mean <- mean(z_window, na.rm = TRUE)
  if (z_range_out) {
    z_lo <- min(z_window, na.rm = TRUE) - z_range
    z_hi <- max(z_window, na.rm = TRUE) + z_range
  } else {
    z_lo <- min(z_window, na.rm = TRUE) + z_range
    z_hi <- max(z_window, na.rm = TRUE) - z_range
  }
  cmin <- min(cn_grid)
  cmax <- max(cn_grid)

  step_lo <- if (exp_ploidy > cmin) max((z_mean - z_lo) / (exp_ploidy - cmin), 1e-6) else 1e-6
  step_hi <- if (exp_ploidy < cmax) max((z_hi  - z_mean) / (cmax - exp_ploidy), 1e-6) else 1e-6

  # Use directional steps so that asymmetric z-deviations (e.g. loss stronger than
  # gain in a homeologous exchange) do not misplace mu on the smaller-deviation side.
  mu_vec <- z_mean + ifelse(
    as.numeric(cn_grid) <= exp_ploidy,
    step_lo * (as.numeric(cn_grid) - exp_ploidy),
    step_hi * (as.numeric(cn_grid) - exp_ploidy)
  )
  mu     <- setNames(mu_vec, as.character(cn_grid))
  mu <- mu[state_ids]
  return(mu)
}


#' Update a multi-sample hmm_CN_multi object with a new single-sample hmm_CN result
#'
#' Replaces or adds the results for a given sample in a multi-sample HMM CN object.
#' If the sample already exists, its data and parameters are replaced; otherwise, the new sample is appended.
#'
#' @param hmm_CN_multi An object of class 'hmm_CN_multi' (list with by_window, by_marker, params_samples).
#' @param hmm_CN An object of class 'hmm_CN' (single-sample result with by_window, by_marker, params).
#' @param rm_sample Logical. If TRUE, removes the sample from hmm_CN_multi without adding the new hmm_CN results (useful for cleanup). Default FALSE.
#'
#' @return An updated hmm_CN_multi object with the new or replaced sample's results.
#' @details
#' This function is useful for incrementally building or updating a multi-sample HMM CN result object
#' as new samples are processed. It ensures that only one entry per sample is present in each component.
#'
#' @export
update_hmm_CN_multi <- function(hmm_CN_multi, hmm_CN, rm_sample = FALSE){

  if(!inherits(hmm_CN_multi, "hmm_CN") || !all(c("by_window", "by_marker", "params_samples") %in% names(hmm_CN_multi)))
    stop("hmm_CN multi must be of class 'hmm_CN' with components: by_window, by_marker, params_samples")

  if(!inherits(hmm_CN, "hmm_CN") || !all(c("by_window", "by_marker", "params") %in% names(hmm_CN)))
    stop("hmm_CN must be of class 'hmm_CN' with components: by_window, by_marker, params")

  hmm_CN_multi_new <- hmm_CN_multi

  sample <- unique(hmm_CN$by_window$Sample)
  sample1 <- unique(hmm_CN$by_marker$SampleName)
  if(any(sample != sample1)) stop("Sample in hmm by window and by marker differ")

  if(any(hmm_CN_multi$by_window$Sample %in% sample)){
    idx <- which(hmm_CN_multi$by_window$Sample %in% sample)
    hmm_CN_multi_new$by_window <- hmm_CN_multi$by_window[-idx,]
  }
  if(any(hmm_CN_multi$by_marker$SampleName %in% sample)){
    idx <- which(hmm_CN_multi$by_marker$SampleName %in% sample)
    hmm_CN_multi_new$by_marker <- hmm_CN_multi$by_marker[-idx,]
  }
  if(any(names(hmm_CN_multi_new$params_samples) %in% sample)){
    idx <- which(names(hmm_CN_multi_new$params_samples) %in% sample)
    hmm_CN_multi_new$params_samples <- hmm_CN_multi_new$params_samples[-idx]
  }

  if(rm_sample) return(hmm_CN_multi_new)

  if(length(sample) == 1) {
    params <- list(hmm_CN$params)
    names(params) <- sample
  } else params <- hmm_CN$params

  by_window <- bind_rows(list(hmm_CN_multi_new$by_window, hmm_CN$by_window))
  by_marker <- bind_rows(list(hmm_CN_multi_new$by_marker, hmm_CN$by_marker))

  idx <- which(colnames(by_window) == "post_max")
  idx1 <- order(colnames(by_window)[(idx + 1):ncol(by_window)])
  by_window <- by_window[,c(1:idx, idx+idx1)]

  hmm_CN_multi_new$params_samples <- c(hmm_CN_multi_new$params_samples, params)
  hmm_CN_multi_new$by_window <- by_window
  hmm_CN_multi_new$by_marker <- by_marker

  return(hmm_CN_multi_new)
}
