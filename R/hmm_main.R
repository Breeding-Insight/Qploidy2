if (getRversion() >= "2.15.1") utils::globalVariables(c("ll_em", "ll_hist"))

#' Polyploid copy-number HMM (z + BAF) per sample and window
#'
#' Aggregates SNP-level data into windows, builds BAF histograms, and fits a
#' Hidden Markov Model that combines a Gaussian emission for depth/z and a
#' polyploid BAF “comb” emission to infer copy number (CN) per window.
#' Returns a tidy data frame with CN calls and posterior probabilities that is
#' ready for plotting (e.g., with \code{\link{plot_cn_track}}) or merging across
#' chromosomes/batches.
#'
#' @param qploidy_standarize_result An object of class \code{qploidy_standardization} as returned by \code{standardize()}, or \code{NULL}. When \code{NULL}, raw data must be supplied via \code{data} and \code{geno.pos}.
#' @param data Optional. A \code{data.frame} with columns \code{MarkerName}, \code{SampleName}, \code{X}, \code{Y}, \code{R}, and \code{ratio}. Required when \code{qploidy_standarize_result} is \code{NULL}.
#' @param geno.pos Optional. A \code{data.frame} with columns \code{MarkerName}, \code{Chromosome}, and \code{Position}. Required when \code{qploidy_standarize_result} is \code{NULL}.
#' @param use_values Character vector of length 2 specifying which columns to use as the BAF-like and depth-like signals. When \code{qploidy_standarize_result} is provided, the first element can be \code{"BAF"} (uses column \code{baf}) or \code{"ratio"} (uses column \code{ratio}), and the second element can be \code{"zscore"} (uses column \code{z}) or \code{"R"} (uses column \code{R}). All four combinations are accepted: \code{c("BAF", "zscore")} (default), \code{c("BAF", "R")}, \code{c("ratio", "zscore")}, \code{c("ratio", "R")}. When \code{qploidy_standarize_result} is \code{NULL}, the only accepted value is \code{c("ratio", "R")}.
#' @param sample_id Character scalar. Sample identifier to analyze; rows are filtered as \code{SampleName == sample_id}.
#' @param chr Optional. Character or integer vector specifying chromosomes to include. If \code{NULL}, all chromosomes are used.
#' @param segment_zscore Logical. If \code{TRUE}, segment z-scores using changepoint detection to define windows instead of fixed SNP counts. Default \code{TRUE}.
#' @param snps_per_window Integer. Number of SNPs per window when windows are created. Default \code{50}.
#' @param reflect Logical. If TRUE (default), apply reflection for continuous
#'   kernels to keep mass within [0,1] (useful for gaussian).
#' @param add_uniform Logical. If TRUE, add a uniform noise component to the
#'   mixture before renormalization. Default FALSE.
#' @param min_snps_per_window Integer. Minimum SNPs required to keep a window. If NULL, a dynamic value is chosen based on chromosome size (see code for details). Windows with fewer SNPs are dropped. Default: NULL (dynamic), or user-specified value.
#' @param cn_grid Integer vector of copy-number states to consider (e.g., \code{1:4}). They are limited to 0-10. Include \code{0} to allow nullisomy calls; see also \code{cn0_background_frac}. Unlikely values will be discarded during estimation if their z means are not monotonic with ploidy.
#' @param cn0_background_frac Numeric in [0, 1]. Array probes emit residual background signal even at CN=0, which elevates the z-score above the linear extrapolation. This parameter shifts the initial \code{mu[0]} upward by \code{cn0_background_frac * (mu[1] - mu[0])}, placing it at the expected background-elevated z rather than the theoretical zero-signal floor. Ignored when \code{0} is not in \code{cn_grid}. Default \code{0} (no background correction; suitable for sequencing data).
#' @param M Integer. Number of BAF histogram bins on [0,1]. Default \code{100}.
#' @param max_iter Integer. Maximum EM iterations. Default \code{60}.
#' @param het_quantile Numeric. Quantile used to scale BAF emission weight based on heterozygote count. Default \code{0.8}.
#' @param baf_weight Numeric. Overall weight applied to BAF emission (0–1). Default \code{0.5}.
#' @param z_range Numeric. Padding added to min/max z for initial mean estimation. Default \code{0.2}.
#' @param transition_jump Numeric. Diagonal value for transition matrix (probability to stay in same CN state). Default \code{0.995}.
#' @param initial_prob Numeric. Initial probability for the best CN state in the initial state distribution (pi0). Default \code{0.5}. Sets the prior probability for the expected ploidy (or best CN from BAF model) at the first window; remaining probability is distributed uniformly across other states. If the best CN is not found, pi0 is uniform across all states.
#' @param z_only Logical. If \code{TRUE}, fit the HMM using the z-emission only (ignores BAF). Default \code{FALSE}.
#' @param baf_only Logical. If \code{TRUE}, fit the HMM using the BAF emission only (ignores z-score). Mutually exclusive with \code{z_only}. Default \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, print progress messages. Default \code{TRUE}.
#' @param exp_ploidy Numeric. Expected ploidy value. If \code{NA} or \code{NULL}, it is set to the best CN from the BAF model. Default \code{NA}.
#' @param rm_outliers Logical. If \code{TRUE}, remove outliers from z-scores before HMM fitting. Default \code{TRUE}.
#' @param outlier_alpha Numeric. Alpha threshold for outlier removal in z-scores. Default 0.05.
#' @param selected_model Optional. An object of class \code{selected_BAF_model}. If provided, it is used instead of running \code{select_best_baf_model}.
#' @param dists Character vector. Distribution families to test in BAF model selection (e.g., c("gaussian", "beta", ...)).
#' @param bw_grid Numeric vector. Bandwidth (SD or concentration) values to try for kernel smoothing in BAF model selection. Default: c(0.02, 0.03, 0.04).
#' @param add_uniform_grid Logical vector. Whether to include a uniform noise component in the BAF template. Default: c(FALSE, TRUE).
#' @param uniform_weight_grid Numeric vector. Mixture weights for the uniform component (when enabled). Default: c(0.01, 0.03, 0.05, 0.10, 0.15).
#' @param param_count Optional named integer vector. Number of free parameters per distribution (for BIC penalty in BAF model selection). If NULL, defaults to 0 for all.
#' @param count_grid_as_params Logical. If TRUE (default), adds +1 to BIC penalty for each hyperparameter tuned by grid search (bw, and uniform_weight if used).
#' @param correct_scale Logical. If TRUE (default), the BAF log-likelihood is corrected by the number of markers with valid BAF values in each window, so that windows with different numbers of markers contribute equally to the combined emission. Prevents windows with many markers from dominating the HMM via the BAF term alone.
#' @param min_het_frac Numeric in [0,1]. Threshold for the fraction of BAF values in \code{het_range} considered heterozygous. If the observed heterozygous fraction exceeds this value, CN=1 is excluded from \code{cn_grid} during BAF model selection and per-window likelihood computation, as a meaningful proportion of heterozygous loci makes haploid (CN=1) implausible. Default \code{0.01}.
#' @param het_range Numeric vector of length 2. BAF interval used to define heterozygous loci (default \code{c(0.05, 0.95)}). Used by the \code{min_het_frac} filter and by \code{plot_heterozygosity}.
#' @param dosage_threshold Numeric in [0,1]. Minimum posterior probability required for a dosage call to be counted as a heterozygote when computing the BAF emission weight (\code{w_baf}) per window. Markers whose maximum dosage posterior falls below this threshold are excluded from the heterozygote count. Default \code{0.6}.
#' @param z_no_baf_scale Numeric in [0,1]. Floor weight applied to \code{w_baf} when computing \code{CN_reliability} for windows with no BAF support (\code{w_baf == 0}). A value of 0 means z-only calls have \code{CN_reliability = 0}; 0.25 (default) allows z-only calls to reach at most 25\% of \code{post_max}. Does not affect the HMM emission or CN calling. Applies only when \code{z_only = FALSE}.
#' @param hom_z_sigma_inflate Numeric >= 1. Two-tier correction for reference-bias-driven z elevation in windows with no BAF support (\code{w_baf == 0}). (1) During decoding: sigma is multiplied by this factor for \code{w_baf == 0} windows, widening the z-emission. (2) Post-decoding override: for chromosomes where \emph{every} window has \code{w_baf == 0}, the chromosome mean z is compared to the sample-baseline mu; if the deviation is within \code{hom_z_sigma_inflate * sigma} the entire chromosome is overridden to the sample-mode CN (reference bias); larger deviations keep the decoded CN (genuine change). Default \code{1.5}. Set to \code{1} to disable both corrections.
#' @param nQuack_EM Logical. If \code{TRUE}, use the \pkg{nQuack} package to estimate \code{exp_ploidy} from
#'   whole-sample read-coverage data (total depth and one randomly chosen allele count per heterozygous site).
#'   The three mixture-model distributions shipped with nQuack (\code{quackNormal}, \code{quackBeta},
#'   \code{quackBetaBinom}) are each fitted and the winner by BIC is used to determine the ploidal level, which
#'   becomes \code{exp_ploidy} for the HMM run. Per-window BAF likelihoods are still computed with the existing
#'   template method. Requires nQuack to be installed (\code{remotes::install_github("mgaynor1/nQuack")}) and
#'   the input data to contain columns \code{R} (total depth), \code{X} (ref allele count), \code{Y} (alt
#'   allele count), \code{ratio}, and \code{geno}. Only takes effect when \code{selected_model} is \code{NULL}
#'   and \code{exp_ploidy} is \code{NA}/\code{NULL}. Default \code{FALSE}.
#'@param nQuack_dists Character vector or \code{NULL}. Restricts which nQuack mixture-model
#'  families are evaluated when \code{nQuack_EM = TRUE}. Valid values are \code{"normal"},
#'  \code{"beta"}, and \code{"beta-binomial"}. Both the base family and its
#'  \code{+uniform} variant are automatically tested for each entry. If \code{NULL} (default),
#'  all three families are evaluated and the best by BIC is selected. Example:
#'  \code{nQuack_dists = "normal"} tests only Normal and Normal+Uniform.
#' @param nQuack_exp_ploidy Logical. If \code{TRUE}, use the \pkg{nQuack} package to
#'   estimate \code{exp_ploidy} only (same as \code{nQuack_EM}), but then use the
#'   template BAF model (\code{select_best_baf_model}) for per-window BAF likelihoods
#'   and dosage assignment. This is a middle-ground option: nQuack provides a more
#'   robust ploidy estimate while the template method handles the window-level inference.
#'   Requires the same input columns as \code{nQuack_EM}. Ignored when \code{nQuack_EM
#'   = TRUE}. Default \code{FALSE}.
#' @param rerun_overall_ploidy Logical. If \code{TRUE}, after the initial HMM run the function identifies markers whose window-level CN call disagrees with the sample-wide mode CN, removes them, and reruns \code{hmm_estimate_CN} (with \code{rerun_overall_ploidy = FALSE}) to obtain a cleaner overall ploidy estimate. Useful when a few aneuploid segments would otherwise bias the genome-wide BAF model. Default \code{FALSE}.
#' @param recycled_obj_rerun_overall_ploidy Named list for internal use during the recursive call triggered by
#'   \code{rerun_overall_ploidy = TRUE}. Contains two elements:
#'   \describe{
#'     \item{\code{rm_mks}}{Character vector of marker names to exclude from the BAF model selection step
#'       (markers whose CN call disagreed with the mode in the first pass).}
#'     \item{\code{recycled_d}}{The pre-processed marker-level data frame (\code{d}) from the first pass,
#'       recycled to skip redundant input validation and outlier removal in the second pass.}
#'   }
#'   Users should never need to set this directly. Default \code{NULL}.
#'
#' @return An object of class \code{hmm_CN}, a list with two elements:
#'   \describe{
#'     \item{result}{A \code{data.frame} with one row per window and columns:
#'       \itemize{
#'         \item \code{Sample}, \code{Chr}, \code{WindowID}, \code{Start}, \code{End}
#'         \item \code{n_snps}, \code{n_het} – counts per window
#'         \item \code{z} – mean z per window
#'         \item \code{w_baf} – weight applied to the BAF emission (0–1)
#'         \item \code{CN_call} – Viterbi copy-number call (factor/integer)
#'         \item \code{post_max} – posterior probability of the called state
#'         \item \code{post_CN<k>} – posterior probability for each state in \code{cn_grid}
#'       }
#'     }
#'     \item{params}{A list of fitted HMM parameters:
#'       \itemize{
#'         \item \code{cn_grid} – copy-number states
#'         \item \code{mu} – estimated z means per state
#'         \item \code{sigma} – estimated z standard deviation
#'         \item \code{A} – transition matrix
#'         \item \code{pi0} – initial state probabilities
#'         \item \code{bins} – number of BAF bins (M)
#'         \item \code{bw} – BAF template bandwidth
#'         \item \code{loglik} – final log-likelihood
#'       }
#'     }
#'   }
#' The \code{result} data frame is suitable for downstream plotting or merging. The \code{params} attribute contains all HMM model parameters for reproducibility and inspection.
#'
#' @details
#' \strong{Windowing.} If \code{chr} is specified, only those chromosomes are analyzed. Windows are formed within each chromosome by grouping every \code{snps_per_window} consecutive SNPs. Windows with fewer than \code{min_snps_per_window} SNPs are removed.
#'
#' \strong{Emissions.} For state \code{c} (copy number), the z-emission is \eqn{\mathcal{N}(\mu_c, \sigma^2)}, where \eqn{\sigma} is a single standard deviation shared across all CN states (it is not scaled per state, even though depth variance is in reality expected to increase with copy number). The BAF emission is a multinomial log-likelihood comparing the observed BAF histogram to a state-specific template (a “comb” with peaks at \eqn{d/c}). Windows with few heterozygotes are down-weighted via \code{w_baf} (derived from \code{n_het}). The combined per-window, per-state log-likelihood is a tempered log-product, \code{llz + w_baf * llb}: the z emission always contributes fully, assuming conditional independence of z and BAF given the CN state, while the BAF emission's evidence is scaled down (tempered) by \code{w_baf} when fewer heterozygotes support it.
#'
#' \strong{Numerical safety.} The implementation guards against non-finite emissions, enforces stochastic \code{A} rows, and lower-bounds \code{sigma}.
#'
#' @section Required helpers:
#' This function calls \code{\link{generate_baf_template}}, \code{\link{baf_log_likelihood}}, \code{\link{logsumexp}}, and \code{\link{em_hmm_cn}} which must be available in your package/namespace.
#'
#' @examples
#' \dontrun{
#' # Suppose qploidy_standarize_result is from standardize(),
#' # with columns: SampleName, Chr, Position, baf, z
#' res <- hmm_estimate_CN(
#'   qploidy_standarize_result = qploidy_standarize_result,
#'   sample_id = "ASample",
#'   chr = c("Chr1", "Chr2"),
#'   snps_per_window = 500,
#'   cn_grid = 2:6
#' )
#' head(res$by_window)
#' res$params
#' }
#'
#' @seealso \code{\link{generate_baf_template}}, \code{\link{baf_log_likelihood}},
#'   \code{\link{viterbi}}, \code{\link{plot_cn_track}}
#' @importFrom stats dnorm quantile sd
#' @importFrom graphics hist
#' @importFrom utils tail
#' @importFrom dplyr filter
#'
#' @export
hmm_estimate_CN <- function(
  qploidy_standarize_result = NULL,
  sample_id,
  data = NULL,
  geno.pos = NULL,
  use_values = c("BAF", "zscore"),
  chr = NULL,
  reflect = TRUE,
  add_uniform = FALSE,
  segment_zscore = TRUE,
  snps_per_window = 50,
  min_snps_per_window = NULL,
  cn_grid = 1:4,
  M = 100,
  max_iter = 60,
  het_quantile = 0.8, # increase this value to reduce the weight of baf when few hets
  baf_weight = 0.5,
  z_range = NULL,
  transition_jump = 0.9, # decrease this value if you think there changes in CN is likely
  initial_prob = 0.5, # Initial probability for the best CN state in the initial state distribution (pi0). Default 0.5. Sets the prior probability for the expected ploidy (or best CN from BAF model) at the first window; remaining probability is distributed uniformly across other states. If the best CN is not found, pi0 is uniform across all states.
  z_only = FALSE,
  baf_only = FALSE,
  verbose = TRUE,
  exp_ploidy = NA,
  rm_outliers = TRUE,
  outlier_alpha = 0.05,
  selected_model = NULL,
  dists = c("gaussian", "beta_binomial", "negative_binomial"),
  bw_grid = c(0.01, 0.02, 0.03, 0.04),
  add_uniform_grid = FALSE,
  uniform_weight_grid = c(0.01, 0.03, 0.05, 0.10, 0.15),
  param_count = NULL,
  count_grid_as_params = TRUE,
  correct_scale = TRUE,
  min_het_frac = 0.01,
  het_range = c(0.05, 0.95),
  dosage_threshold = 0.6,
  z_no_baf_scale = 0.25,
  hom_z_sigma_inflate = 1.5,
  cn0_background_frac = 0,
  nQuack_EM = FALSE,
  nQuack_cores = 1L,
  nQuack_dists = NULL,
  nQuack_exp_ploidy = FALSE,
  rerun_overall_ploidy = TRUE,
  recycled_obj_rerun_overall_ploidy = NULL
) {
  # --- input checks ---
  # If second cycle of rerun_overall_ploidy, check that recycled_obj_rerun_overall_ploidy is provided
  if (!rerun_overall_ploidy && !is.null(recycled_obj_rerun_overall_ploidy)) {
    d <- recycled_obj_rerun_overall_ploidy$recycled_d
    baf_col <- recycled_obj_rerun_overall_ploidy$recycled_baf_col
    z_col <- recycled_obj_rerun_overall_ploidy$recycled_z_col
    cn_grid <- recycled_obj_rerun_overall_ploidy$recycled_cn_grid
  } else {
    vmsg("Preparing inputs and applying initial filters", verbose = verbose, level = 0, type = ">>")

    if(any(cn_grid < 0) | any(cn_grid > 10)) stop("All values in cn_grid must be between 0 and 10.")
    if (z_only && baf_only) stop("z_only and baf_only cannot both be TRUE.")

    if (!is.character(sample_id) || length(sample_id) != 1 || nchar(sample_id) == 0) {
      stop("sample_id must be a non-empty character scalar.")
    }

    # Validate and resolve use_values
    valid_uv_std <- list(c("BAF", "zscore"), c("BAF", "R"), c("ratio", "zscore"), c("ratio", "R"))
    valid_uv_raw <- c("ratio", "R")
    if (!is.null(qploidy_standarize_result)) {
      if (!is(qploidy_standarize_result, "qploidy_standardization")) {
        stop("Input must be a qploidy_standardization object as returned by standardize().")
      }
      if (!any(sapply(valid_uv_std, identical, use_values))) {
        stop("When qploidy_standarize_result is provided, use_values must be one of: c(\"BAF\", \"zscore\"), c(\"BAF\", \"R\"), c(\"ratio\", \"zscore\"), c(\"ratio\", \"R\").")
      }
      baf_col <- if (use_values[1] == "BAF") "baf" else "ratio"
      z_col <- if (use_values[2] == "zscore") "z" else "R"
      df <- as.data.frame(qploidy_standarize_result$data)
      if (nrow(df) == 0) stop("Input data is empty. No SNPs found for any sample.")
      d <- filter(df, SampleName == sample_id)
      if (nrow(d) == 0) stop(sprintf("No data found for sample '%s'.", sample_id))
    } else {
      if (!identical(use_values, valid_uv_raw)) {
        stop("When qploidy_standarize_result is NULL, use_values must be c(\"ratio\", \"R\").")
      }
      if (is.null(data) || is.null(geno.pos)) {
        stop("When qploidy_standarize_result is NULL, both 'data' and 'geno.pos' must be provided.")
      }
      req_data <- c("MarkerName", "SampleName", "ratio", "R")
      req_gp <- c("MarkerName", "Chromosome", "Position")
      miss_d <- setdiff(req_data, names(data))
      miss_gp <- setdiff(req_gp, names(geno.pos))
      if (length(miss_d)) stop(paste("'data' is missing columns:", paste(miss_d, collapse = ", ")))
      if (length(miss_gp)) stop(paste("'geno.pos' is missing columns:", paste(miss_gp, collapse = ", ")))
      baf_col <- "ratio"
      z_col <- "R"
      d <- data[data$SampleName == sample_id, , drop = FALSE]
      if (nrow(d) == 0) stop(sprintf("No data found for sample '%s'.", sample_id))
      gp <- geno.pos[, req_gp]
      colnames(gp)[colnames(gp) == "Chromosome"] <- "Chr"
      d <- merge(d, gp, by = "MarkerName", all.x = FALSE)
      if (nrow(d) == 0) stop("No markers remain after merging data with geno.pos.")
    }

    # subset to chromosomes
    if (!is.null(chr)) {
      if (is.numeric(chr)) {
        # index into the observed (non-NA) chromosome labels, not raw positions,
        # so an out-of-range/NA entry in unique(d$Chr) can't leak into chrs
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

    # Replace R == 0 with NA (zero depth is uninformative and distorts z/R distributions)
    if ("R" %in% names(d)) d$R[d$R == 0] <- NA

    # Remove outliers from depth/z scores
    if (rm_outliers) {
      n_na <- sum(is.na(d[[z_col]]))
      if (z_col != "z") {
        # Stash any existing "z" column to avoid creating a duplicate when renaming z_col -> "z"
        if ("z" %in% names(d)) names(d)[names(d) == "z"] <- ".__z_stash__."
        names(d)[names(d) == z_col] <- "z"
      }
      d <- rm_outlier(d, z = TRUE, alpha = outlier_alpha)
      if (z_col != "z") {
        names(d)[names(d) == "z"] <- z_col
        if (".__z_stash__." %in% names(d)) names(d)[names(d) == ".__z_stash__."] <- "z"
      }
      vmsg("Outliers removed from %s scores: %s", verbose = verbose, level = 1, type = ">>", z_col, sum(is.na(d[[z_col]])) - n_na)
    }

    # Remove markers with missing BAF/ratio and z/R
    rm <- which(is.na(d[[z_col]]))
    if (length(rm) > 0) d <- d[-rm, , drop = FALSE]

    vmsg("Inputs good to go!", verbose = verbose, level = 1, type = ">>")

    if (rerun_overall_ploidy) {
      keep_d <- d
      keep_baf_col <- baf_col
      keep_z_col <- z_col
      keep_cn_grid <- cn_grid

    } else {
      keep_d <- NULL
    }
  }
  # --- set expected ploidy ---
  # Calculate expected ploidy using sample-level BAF distribution
  vmsg("Setting internal parameters", verbose = verbose, level = 0, type = ">>")

  if (!is.null(selected_model)) {
    if (!inherits(selected_model, "selected_BAF_model")) {
      stop("selected_model argument must be of class 'selected_BAF_model'.")
    }
    # Use user-provided selected_model
    vmsg("Using user-provided selected_model object", verbose = verbose, level = 1, type = ">>")
  } else {
    # This step was added to be activated if rerun_overall_ploidy is TRUE
    if (!is.null(recycled_obj_rerun_overall_ploidy$rm_mks)) {
      bafs_overall <- d[[baf_col]][-which(d$MarkerName %in% recycled_obj_rerun_overall_ploidy$rm_mks)]
    } else {
      bafs_overall <- d[[baf_col]]
    }

    if (isTRUE(nQuack_EM) && !z_only) {
      # nQuack_EM = TRUE (full nQuack path): skip template BAF model fitting.
      # A synthetic selected_model will be built from nQuack results below.
      # nQuack_exp_ploidy is irrelevant here: nQuack_EM takes full precedence.
      selected_model <- NULL
    } else {
      selected_model <- select_best_baf_model(
        baf_vec = bafs_overall,
        sample = sample_id,
        cn_grid = cn_grid[cn_grid >= 1L],  # BAF is undefined at CN=0
        dists = dists,
        M = M,
        reflect = reflect,
        bw_grid = bw_grid,
        add_uniform_grid = add_uniform_grid,
        uniform_weight_grid = uniform_weight_grid,
        param_count = param_count,
        count_grid_as_params = count_grid_as_params,
        min_het_frac = min_het_frac,
        het_range = het_range
      )
    }
  }

  # Initialise nQuack outputs to NULL so downstream exists() guards are safe
  # even when the user supplies exp_ploidy directly (skipping the nQuack block).
  if (!exists("nq_dist"))           nq_dist           <- NULL
  if (!exists("nq_ploidy"))         nq_ploidy         <- NULL
  if (!exists("nq_winner"))         nq_winner         <- NULL
  if (!exists("nq_uniform_weight")) nq_uniform_weight <- NULL

  # Use exp_ploidy argument if provided, otherwise use selected_model$best$best_cn
  if (is.null(exp_ploidy) || (length(exp_ploidy) == 1 && is.na(exp_ploidy))) {
    # --- nQuack EM alternative: use nQuack to determine exp_ploidy ---
    if ((isTRUE(nQuack_EM) || isTRUE(nQuack_exp_ploidy)) && !z_only) {
      # On second pass (rerun_overall_ploidy), reuse the winning distribution from the first
      # pass so only that one distribution is re-fit on the cleaner data (faster, consistent).
      nq_cached_dist <- recycled_obj_rerun_overall_ploidy$recycled_nq_dist
      if (!is.null(nq_cached_dist)) {
        vmsg("nQuack: re-estimating exp_ploidy with winning distribution (%s) only",
             verbose = verbose, level = 1, type = ">>", nq_cached_dist)
      } else {
        vmsg("nQuack: estimating exp_ploidy via nQuack EM mixture model", verbose = verbose, level = 1, type = ">>")
      }
      nq_result <- .nquack_estimate_ploidy(
        d                = d,
        sample_id        = sample_id,
        cn_grid          = cn_grid,
        baf_col          = baf_col,
        selected_model   = selected_model,
        nQuack_cores     = nQuack_cores,
        het_range        = het_range,
        dosage_threshold = dosage_threshold,
        fixed_dist       = nq_cached_dist,
        allowed_dists    = nQuack_dists,
        verbose          = verbose
      )
      nq_ploidy <- nq_result$cn
      nq_winner <- nq_result$winner
      nq_dist   <- nq_result$dist
      if (!is.null(nq_ploidy) && is.numeric(nq_ploidy) && is.finite(nq_ploidy)) {
        exp_ploidy <- nq_ploidy
        vmsg("nQuack estimated exp_ploidy: %s", verbose = verbose, level = 1, type = ">>", exp_ploidy)
      } else {
        stop(sprintf(
          paste0("nQuack: nQuack could not determine exp_ploidy for sample '%s'. ",
                 "nQuack returned: cn=%s, dist=%s, winner=%s. ",
                 "Check that the BAF data has sufficient heterozygous markers. ",
                 "Set nQuack_EM = FALSE and nQuack_exp_ploidy = FALSE to use the template BAF model instead."),
          sample_id,
          if (is.null(nq_ploidy)) "NULL" else as.character(nq_ploidy),
          if (is.null(nq_dist))   "NULL" else as.character(nq_dist),
          if (is.null(nq_winner)) "NULL" else as.character(nq_winner)
        ))
      }

      # On the nQuack_exp_ploidy path the template selected_model was already built
      # above by select_best_baf_model; only override it on the full nQuack_EM path.
      if (isTRUE(nQuack_EM)) {
      # Build a synthetic selected_BAF_model-like object from nQuack results so
      # that all downstream code that reads selected_model$best$* continues to work.
      nq_add_uniform     <- grepl("-uniform$", tolower(nq_dist))
      nq_dist_base       <- sub("-uniform$", "", tolower(nq_dist))
      nq_dist_qploidy    <- switch(nq_dist_base,
        "normal"       = "gaussian",
        "beta"         = "beta",
        "beta-binomial"= "beta_binomial",
        nq_dist_base          # fallback: keep as-is
      )
      # Default uniform_weight if uniform component is present (nQuack does not report it)
      nq_uniform_weight  <- if (nq_add_uniform) 0.05 else 0
      # bw is NA on the nQuack path; dosage assignment uses call_BAF_nQuack_dosages instead.
      selected_model <- structure(
        list(
          best = list(
            dist           = nq_dist_qploidy,
            bw             = NA_real_,
            add_uniform    = nq_add_uniform,
            uniform_weight = nq_uniform_weight,
            best_cn        = exp_ploidy
          )
        ),
        class = "selected_BAF_model"
      )
      vmsg("Synthetic selected_BAF_model built from nQuack results (dist=%s, add_uniform=%s, best_cn=%s)",
           verbose = verbose, level = 1, type = ">>",
           nq_dist_qploidy, nq_add_uniform, exp_ploidy)
      } # end if (isTRUE(nQuack_EM))
    } else {
      exp_ploidy <- selected_model$best$best_cn
      vmsg("exp_ploidy not provided, using best CN from BAF model: %s", verbose = verbose, level = 1, type = ">>", exp_ploidy)
    }
  } else {
    vmsg("exp_ploidy provided by user: %s", verbose = verbose, level = 1, type = ">>", exp_ploidy)
    # When nQuack_EM = TRUE but exp_ploidy was supplied by the user, nQuack is not
    # run and selected_model is still NULL (select_best_baf_model was skipped).
    # Build a minimal synthetic selected_BAF_model so downstream code that reads
    # selected_model$best$* does not crash.  bw = NA because dosage assignment will
    # use call_BAF_nQuack_dosages (nq_dist will be set on first nQuack window call)
    # or, if nQuack_EM = FALSE, this branch is unreachable without a model.
    if (isTRUE(nQuack_EM) && is.null(selected_model)) {
      selected_model <- structure(
        list(best = list(
          dist           = "gaussian",
          bw             = NA_real_,
          add_uniform    = FALSE,
          uniform_weight = 0,
          best_cn        = exp_ploidy
        )),
        class = "selected_BAF_model"
      )
      vmsg("nQuack_EM = TRUE with user-supplied exp_ploidy: synthetic selected_BAF_model built",
           verbose = verbose, level = 1, type = ">>")
    }
  }

  vmsg("Parameters ready", verbose = verbose, level = 1, type = ">>")
  # --- build windows ---
  vmsg("Defining windows", verbose = verbose, level = 0, type = ">>")
  d[["Position"]] <- as.numeric(d[["Position"]])
  d <- d[order(d[["Chr"]], d[["Position"]]), ]

  if (is.null(min_snps_per_window)) {
    # If min_snp_per_window is not defined, the minimum will be set based on the smaller chromosome number
    floor <- 5 # hard default
    frac <- 0.20 # hard default
    mk_by_chrom <- d %>%
      group_by(Chr) %>%
      summarize(n = n())
    x <- min(mk_by_chrom$n)
    m <- max(floor, floor(x * frac))
    min_snps_per_window <- min(m, floor(x / 2))
  }

  # Segmented z-score
  if (segment_zscore) {
    vmsg("Using z-scores changepoint detection to define windows", verbose = verbose, level = 1, type = ">>")
    if (z_col != "z") {
      # Stash any existing "z" column to avoid creating a duplicate when renaming z_col -> "z"
      if ("z" %in% names(d)) names(d)[names(d) == "z"] <- ".__z_stash__."
      names(d)[names(d) == z_col] <- "z"
    }
    d <- add_changepoint_windows(dat = d, minseglen = min_snps_per_window)
    if (z_col != "z") {
      names(d)[names(d) == "z"] <- z_col
      if (".__z_stash__." %in% names(d)) names(d)[names(d) == ".__z_stash__."] <- "z"
    }
  } else {
    # simple fixed-size windows
    vmsg("Using user-defined fixed-size intervals to define windows", verbose = verbose, level = 1, type = ">>")
    # equal-SNP windows within chromosome
    d$.__w__ <- with(
      d,
      ave(
        seq_len(nrow(d)),
        d[["Chr"]],
        FUN = function(k) {
          n <- length(k)

          # initial windowing (same idea as before)
          w <- ceiling(seq_along(k) / snps_per_window)

          # if there is a remainder AND we have more than one full window,
          # merge the last (small) window into the previous one
          if (n > snps_per_window && n %% snps_per_window != 0L) {
            last <- max(w)
            w[w == last] <- last - 1L
          }
          w
        }
      )
    )
  }
  win_col <- ".__w__"

  # summarize per window
  vmsg("Summarizing data by window", verbose = verbose, level = 1, type = ">>")

  # Grouped summarise (hash-based grouping) instead of by()+do.call(rbind), which is quadratic in the number of groups.
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
  if (any(is.nan(win_df$z_mean))) {
    warning("Some windows have no z-score values.")
  }
  if (is.null(win_df) || nrow(win_df) == 0) stop("No windows could be formed. Check input data and window parameters.")

  # order windows first by chr then start
  o <- order(win_df$Chr, win_df$Start)
  win_df <- win_df[o, ]
  rownames(win_df) <- NULL
  W <- nrow(win_df) # W = number of windows

  # extract BAF vectors per window (list) and z vector
  # Split once by (Chr, WindowID) instead of rescanning all of 'd' per window (O(N) instead of O(W x N)).
  d_key     <- paste(d[["Chr"]], d[[win_col]], sep = "\r")
  win_key   <- paste(win_df$Chr, win_df$WindowID, sep = "\r")
  baf_split <- split(d[[baf_col]], d_key)
  baf_list  <- baf_split[match(win_key, names(baf_split))]
  baf_list  <- lapply(baf_list, function(x) if (is.null(x)) numeric(0) else x)
  for (i in seq_len(W)) {
    if (length(baf_list[[i]]) == 0 || all(is.na(baf_list[[i]]))) {
      warning(sprintf(
        "Window %d (Chr %s, WindowID %s) has no valid BAF values. Only z-score likelihood will considered for those",
        i, win_df$Chr[i], win_df$WindowID[i]
      ))
    }
  }
  z <- win_df$z_mean
  if (length(z) == 0 || all(is.na(z))) stop("No valid z values found in any window.")

  vmsg("Windows ready", verbose = verbose, level = 1, type = ">>")

  keep <- win_df$n_snps >= min_snps_per_window
  if (!any(keep)) stop("All windows dropped by min_snps_per_window filter. Try lowering min_snps_per_window or check input data.")
  # n_het is now filled later using vectorized dosages
  # Handle single-window case: assign CN by BAF likelihood only, skip HMM/EM
  if (sum(keep) == 1) {
    vmsg("Only one window remains after filtering. Assigning CN by BAF likelihood only", verbose = verbose, level = 1, type = ">>")
    # BAF is undefined at CN=0, so it must be excluded from the grid used here
    cn_grid_baf <- cn_grid[cn_grid >= 1L]
    win_idx <- which(keep)
    use_nq_single <- isTRUE(nQuack_EM) &&
      !is.null(nq_dist) && !is.na(nq_dist)
    if (use_nq_single) {
      ll_baf_matrix <- .nquack_baf_likelihoods_by_window(
        win_df    = win_df[win_idx[1], , drop = FALSE],
        d_all     = d,
        baf_col   = baf_col,
        cn_grid   = cn_grid_baf,
        nq_dist   = nq_dist,
        het_range = het_range,
        min_het_frac = min_het_frac
      )
    } else {
      ll_baf_matrix <- compute_baf_likelihoods(baf_list[keep][[1]],
        cn_grid_baf,
        M = M,
        bw = selected_model$best$bw,
        plot = FALSE,
        dist = selected_model$best$dist,
        reflect = reflect,
        add_uniform = selected_model$best$add_uniform,
        uniform_weight = selected_model$best$uniform_weight,
        min_het_frac = min_het_frac,
        het_range = het_range
      )
    }

    cn_call <- cn_grid_baf[which.max(ll_baf_matrix$ll_vec)]
    post_max <- rep(1, 1)
    post_df <- as.data.frame(matrix(0, nrow = 1, ncol = length(cn_grid)))
    names(post_df) <- paste0("post_CN", cn_grid)
    post_df[1, paste0("post_CN", cn_grid_baf[which.max(ll_baf_matrix$prob_vec)])] <- 1
    # n_het is now calculated using dosages
    if (use_nq_single) {
      # nQuack path: use pir-based dosage assignment
      win_chr_s   <- win_df$Chr[win_idx[1]]
      win_start_s <- win_df$Start[win_idx[1]]
      win_end_s   <- win_df$End[win_idx[1]]
      d_win_s <- d[!is.na(d$Chr) & d$Chr == win_chr_s &
                   !is.na(d$Position) &
                   d$Position >= win_start_s & d$Position <= win_end_s, , drop = FALSE]
      bc_s <- if (baf_col %in% colnames(d_win_s)) baf_col else if ("ratio" %in% colnames(d_win_s)) "ratio" else NA_character_
      dosages <- list(call_BAF_nQuack_dosages(
        baf_vec   = if (!is.na(bc_s)) d_win_s[[bc_s]] else baf_list[keep][[1]],
        cn        = max(2L, cn_call),  # nQuack requires cn >= 2; CN=1 windows use cn=2
        nq_dist   = nq_dist,
        het_range = het_range,
        R_vec     = if ("R" %in% colnames(d_win_s)) d_win_s$R else NULL,
        Y_vec     = if ("Y" %in% colnames(d_win_s)) d_win_s$Y else NULL
      ))
    } else {
      dosages <- mapply(function(x, y) {
        call_BAF_dosages(x,
          cn = y,
          bw = selected_model$best$bw,
          plot = FALSE,
          dist = selected_model$best$dist,
          add_uniform = selected_model$best$add_uniform,
          uniform_weight = selected_model$best$uniform_weight
        )
      }, baf_list[keep], cn_call, SIMPLIFY = FALSE)
    }
    n_het <- sum(dosages[[1]]$data$dosage != 0 & dosages[[1]]$data$dosage != cn_call[1], na.rm = TRUE)
    result <- cbind(
      data.frame(
        Sample = sample_id,
        Chr = win_df$Chr[keep],
        WindowID = win_df$WindowID[keep],
        Start = win_df$Start[keep],
        End = win_df$End[keep],
        n_snps = win_df$n_snps[keep],
        n_het = n_het,
        z = win_df$z_mean[keep],
        w_baf = 1,
        CN_call = cn_call,
        CN_reliability = post_max,
        post_max = post_max,
        stringsAsFactors = FALSE
      ),
      post_df
    )
    params <- list(
      cn_grid = cn_grid,
      distribution = selected_model$best$dist,
      mu = NA,
      sigma = NA,
      A = NA,
      pi0 = NA,
      bins = M,
      bw = if (isTRUE(nQuack_EM)) NA_real_ else selected_model$best$bw,
      loglik = NA,
      z_range = z_range,
      het_quantile = het_quantile,
      baf_weight = baf_weight,
      transition_jump = transition_jump,
      z_only = z_only,
      exp_ploidy = exp_ploidy,
      rm_outliers = rm_outliers,
      outlier_alpha = outlier_alpha,
      segment_zscore = segment_zscore,
      snps_per_window = snps_per_window,
      min_snps_per_window = min_snps_per_window,
      add_uniform = add_uniform,
      uniform_weight = selected_model$best$uniform_weight
    )
    window_map <- match(d$.__w__, result$WindowID)
    d$w_baf          <- result$w_baf[window_map]
    d$CN_call        <- result$CN_call[window_map]
    d$post_max       <- result$post_max[window_map]
    d$CN_reliability <- result$CN_reliability[window_map]
    vmsg("Done!", verbose = verbose, level = 1, type = ">>")
    return(structure(list(by_window = result, by_marker = d, params = params), class = "hmm_CN"))
  }

  # When z_only = FALSE and CN=0 is in the grid, call hmm_estimate_CN with z_only=TRUE
  # to identify nullisomic windows, then strip their BAF before the main run.
  if (!z_only && 0 %in% cn_grid) {
    vmsg("Running z_only=TRUE pre-pass to detect CN=0 windows", verbose = verbose, level = 0, type = ">>")
    pre_z <- hmm_estimate_CN(
      qploidy_standarize_result = qploidy_standarize_result,
      data            = data,
      geno.pos        = geno.pos,
      use_values      = use_values,
      sample_id       = sample_id,
      chr             = chr,
      segment_zscore  = segment_zscore,
      snps_per_window = snps_per_window,
      reflect         = reflect,
      add_uniform     = add_uniform,
      min_snps_per_window = min_snps_per_window,
      cn_grid         = cn_grid,
      cn0_background_frac = cn0_background_frac,
      M               = M,
      max_iter        = max_iter,
      het_quantile    = het_quantile,
      baf_weight      = baf_weight,
      z_range         = z_range,
      transition_jump = transition_jump,
      initial_prob    = initial_prob,
      z_only          = TRUE,
      baf_only        = FALSE,
      verbose         = FALSE,
      exp_ploidy      = exp_ploidy,
      rm_outliers     = rm_outliers,
      outlier_alpha   = outlier_alpha,
      selected_model  = selected_model,
      correct_scale   = correct_scale,
      min_het_frac    = min_het_frac,
      het_range       = het_range,
      dosage_threshold = dosage_threshold,
      z_no_baf_scale  = z_no_baf_scale,
      hom_z_sigma_inflate = hom_z_sigma_inflate,
      nQuack_EM       = FALSE,   # z_only pre-pass: never use nQuack
      rerun_overall_ploidy = FALSE
    )
    cn0_key <- with(pre_z$by_window[is.na(pre_z$by_window$CN_call), ],
                    paste(Chr, WindowID, sep = "__"))
    win_key <- paste(win_df$Chr, win_df$WindowID, sep = "__")
    cn0_win_idx <- which(win_key %in% cn0_key)
    if (length(cn0_win_idx) > 0) {
      vmsg("CN=0 pre-pass: dropping %d window(s) from analysis", verbose = verbose, level = 1, type = ">>", length(cn0_win_idx))
      cn0_win_df  <- win_df[cn0_win_idx, , drop = FALSE]
      baf_list    <- baf_list[-cn0_win_idx]
      z           <- z[-cn0_win_idx]
      win_df      <- win_df[-cn0_win_idx, , drop = FALSE]
    } else {
      cn0_win_df <- NULL
    }
    cn_grid <- cn_grid[cn_grid != 0]  # Remove CN=0 from the grid for the main run
  } else {
    cn0_win_df <- NULL
  }

  # Generate BAF likelihoods and probabilities per window
  # Uses parameters from selected_model (template) or nQuack EM per-window
  if (!z_only) {
    vmsg("Generating BAF likelihoods by window", verbose = verbose, level = 0, type = ">>")
    grid1 <- unique(c(1L, cn_grid[cn_grid >= 1L])) # BAF undefined at CN=0; always test CN=1 for LOH

    # Determine which window-level BAF likelihood method to use
    use_nquack_windows <- isTRUE(nQuack_EM) &&
      !is.null(nq_dist) && !is.na(nq_dist)

    if (use_nquack_windows) {
      vmsg("nQuack_EM: computing per-window BAF likelihoods via nQuack EM (%s)",
           verbose = verbose, level = 1, type = ">>", nq_dist)
      baf_results <- lapply(seq_along(baf_list), function(i) {
        .nquack_baf_likelihoods_by_window(
          win_df    = win_df[i, , drop = FALSE],
          d_all     = d,
          baf_col   = baf_col,
          cn_grid   = grid1,
          nq_dist   = nq_dist,
          het_range = het_range,
          min_het_frac = min_het_frac
        )
      })
    } else {
      baf_results <- lapply(baf_list, function(baf_vec) {
        compute_baf_likelihoods(baf_vec,
          grid1,
          M = M,
          bw = selected_model$best$bw,
          plot = FALSE,
          dist = selected_model$best$dist,
          reflect = reflect,
          add_uniform = selected_model$best$add_uniform,
          uniform_weight = selected_model$best$uniform_weight,
          het_range = het_range,
          min_het_frac = min_het_frac
        )
      })
    }

    ll_baf_matrix <- do.call(rbind, lapply(baf_results, function(res) res$ll_vec))
    ploidies_temp <- apply(ll_baf_matrix, 1, which.max)
    ploidies_temp <- grid1[ploidies_temp]

    colnames(ll_baf_matrix) <- paste0("CN", grid1)
    if (!any(cn_grid == 1)) {
      ll_baf_matrix <- ll_baf_matrix[, -which(colnames(ll_baf_matrix) == "CN1"), drop = FALSE]
    }
    # Pad zero-likelihood for cn_grid states excluded from BAF (e.g. CN=0) and align columns to cn_grid order
    cur_cns <- as.integer(sub("CN", "", colnames(ll_baf_matrix)))
    missing_cns <- setdiff(cn_grid, cur_cns)
    if (length(missing_cns) > 0) {
      zero_cols <- matrix(0, nrow = nrow(ll_baf_matrix), ncol = length(missing_cns),
                          dimnames = list(NULL, paste0("CN", missing_cns)))
      ll_baf_matrix <- cbind(ll_baf_matrix, zero_cols)[, paste0("CN", cn_grid), drop = FALSE]
    }
    vmsg("BAF likelihoods generated", verbose = verbose, level = 1, type = ">>")

    # baf_only: replace global likelihoods with per-chromosome fits so each
    # chromosome uses a bandwidth optimised for its own BAF pattern.
    # When nQuack_EM = TRUE the per-window likelihoods above already come from
    # nQuack; re-fitting with the template here would use bw = NA and fail.
    if (baf_only && !isTRUE(nQuack_EM)) {
      vmsg("baf_only: refitting per-chromosome BAF models", verbose = verbose, level = 1, type = ">>")
      for (chr_i in unique(win_df$Chr)) {
        chr_idx  <- which(win_df$Chr == chr_i)
        chr_baf  <- Filter(Negate(is.na), unlist(baf_list[chr_idx]))
        if (length(chr_baf) < 10) next
        chr_mod  <- select_best_baf_model(
          baf_vec = chr_baf, sample = chr_i,
          cn_grid = cn_grid[cn_grid >= 1L],
          dists = dists, M = M, reflect = reflect,
          bw_grid = bw_grid, add_uniform_grid = add_uniform_grid,
          uniform_weight_grid = uniform_weight_grid,
          param_count = param_count, count_grid_as_params = count_grid_as_params,
          min_het_frac = min_het_frac, het_range = het_range
        )
        chr_res <- lapply(baf_list[chr_idx], function(bv) {
          compute_baf_likelihoods(bv, grid1, M = M, bw = chr_mod$best$bw,
            plot = FALSE, dist = chr_mod$best$dist, reflect = reflect,
            add_uniform = chr_mod$best$add_uniform,
            uniform_weight = chr_mod$best$uniform_weight,
            het_range = het_range, min_het_frac = min_het_frac)
        })
        chr_ll <- do.call(rbind, lapply(chr_res, function(r) r$ll_vec))
        colnames(chr_ll) <- paste0("CN", grid1)
        if (!any(cn_grid == 1L))
          chr_ll <- chr_ll[, colnames(chr_ll) != "CN1", drop = FALSE]
        for (cn_nm in setdiff(colnames(ll_baf_matrix), colnames(chr_ll)))
          chr_ll <- cbind(chr_ll, matrix(0, nrow(chr_ll), 1, dimnames = list(NULL, cn_nm)))
        ll_baf_matrix[chr_idx, ] <- chr_ll[, colnames(ll_baf_matrix), drop = FALSE]
      }
    }

    vmsg("Generating BAF weights by window", verbose = verbose, level = 0, type = ">>")

    vmsg("Calling dosages", verbose = verbose, level = 1, type = ">>")
    use_nquack_dosages <- isTRUE(nQuack_EM) &&
      !is.null(nq_dist) && !is.na(nq_dist)
    if (use_nquack_dosages) {
      # nQuack path: use pir-based dosage assignment per window
      # Build per-window data frames (BAF + depth) from d for call_BAF_nQuack_dosages
      bc_main <- if (baf_col %in% colnames(d)) baf_col else if ("ratio" %in% colnames(d)) "ratio" else NA_character_
      has_R   <- "R" %in% colnames(d)
      has_Y   <- "Y" %in% colnames(d)
      dosages <- lapply(seq_along(baf_list), function(i) {
        chr_i   <- win_df$Chr[i]
        start_i <- win_df$Start[i]
        end_i   <- win_df$End[i]
        d_wi <- d[!is.na(d$Chr) & d$Chr == chr_i &
                  !is.na(d$Position) &
                  d$Position >= start_i & d$Position <= end_i, , drop = FALSE]
        baf_i <- if (!is.na(bc_main) && nrow(d_wi) > 0) d_wi[[bc_main]] else baf_list[[i]]
        call_BAF_nQuack_dosages(
          baf_vec   = baf_i,
          cn        = max(2L, ploidies_temp[i]),  # nQuack requires cn >= 2; CN=1 windows use cn=2
          nq_dist   = nq_dist,
          het_range = het_range,
          R_vec     = if (has_R && nrow(d_wi) > 0) d_wi$R else NULL,
          Y_vec     = if (has_Y && nrow(d_wi) > 0) d_wi$Y else NULL
        )
      })
    } else {
      dosages <- mapply(
        function(x, y) {
          call_BAF_dosages(x,
            cn = y,
            bw = selected_model$best$bw,
            plot = FALSE,
            dist = selected_model$best$dist,
            add_uniform = selected_model$best$add_uniform,
            uniform_weight = selected_model$best$uniform_weight
          )
        },
        baf_list, ploidies_temp,
        SIMPLIFY = FALSE
      )
    }

    # For each window, count heterozygotes: dosage != 0 & dosage != ploidies_temp & dosage_prob > dosage_threshold
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

    # Warn about windows where all heterozygotes were discarded by the dosage_threshold filter
    prob_filtered_idx <- which(n_het_before_thresh > 0 & n_het_window == 0)
    if (length(prob_filtered_idx) > 0) {
      prob_filtered_labels <- paste0(
        win_df$Chr[prob_filtered_idx], ":", win_df$WindowID[prob_filtered_idx]
      )
      vmsg(
        paste0(
          "Warning: %d window(s) had all heterozygous markers discarded due to low dosage ",
          "probability (< dosage_threshold = %.2f), resulting in BAF weight = 0 for those windows. ",
          "Only z-score will contribute to CN calling there.\n  Affected windows (Chr:WindowID): %s"
        ),
        verbose = verbose, level = 2, type = ">>",
        length(prob_filtered_idx),
        dosage_threshold,
        paste(prob_filtered_labels, collapse = ", ")
      )
    }

    # Count number of markers with BAF values in each window
    n_baf <- sapply(baf_list, function(x) sum(!is.na(x)))

    # BAF weights from heterozygote counts
    # Lower the number of heterozygotes in the window, lower the weight of the BAF emission
    # If there are no heterozygous, only z-score will be considered
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
    if (round(exp_ploidy, 0) == 1) w_baf <- rep(1, W) # if expected ploidy is 1, we don't expect heterozygotes, so BAF likelihood should be taken fully into account regardless of heterozygote count
    w_baf <- w_baf * baf_weight
  } else {
    n_het_window <- NA
    n_baf <- rep(0L, W)
    ll_baf_matrix <- matrix(0, nrow = W, ncol = length(cn_grid))
    w_baf <- rep(0, W)
    if (round(exp_ploidy, 0) == 1) {
      w_baf <- rep(1, W) # if expected ploidy is 1, we don't expect heterozygotes, so BAF likelihood should be taken fully into account regardless of heterozygote count
      w_baf <- w_baf * baf_weight
    }
  }
  vmsg("BAF weights defined", verbose = verbose, level = 1, type = ">>")

  # --- HMM for all chromosomes together ---
  # fit a single HMM across all windows (all chromosomes)
  vmsg("Setting initial and transitions HMM matrices", verbose = verbose, level = 0, type = ">>")

  # Setup for all windows
  K <- length(cn_grid)
  state_ids <- as.character(cn_grid)

  # Transition matrix and initial state distribution
  # A is K×K. Before normalizing, every off-diagonal entry is 0.001 (rare jumps), and the diagonal is 0.995 (strong tendency to stay in the same CN state).
  A <- matrix(1e-3, K, K)
  diag(A) <- transition_jump
  # A <- A / rowSums(A) makes each row sum to 1 (a valid Markov matrix). After this, values are almost unchanged (just scaled so each row sums exactly to 1).
  A <- A / rowSums(A)
  # To consider: If is possible to change this matrix to favor specific small jumps, like if changing from 2->4 is more likely than 2->6
  # pi0 is the initial state distribution at the first window; here it’s uniform (no prior preference for any CN at the start).
  # Set pi0 to favor exp_ploidy (if provided) or best CN from BAF model
  pi0 <- if (K > 1) rep((1 - initial_prob) / (K - 1), K) else 1
  if (!is.null(exp_ploidy) && !is.na(exp_ploidy)) {
    best_cn <- as.numeric(exp_ploidy)
  } else {
    best_cn <- as.numeric(selected_model$best$best_cn)
  }
  best_idx <- which(as.numeric(cn_grid) == best_cn)
  if (length(best_idx) == 1) {
    pi0[best_idx] <- initial_prob
  } else {
    pi0 <- rep(1 / K, K) # fallback to uniform if best CN not found
  }

  vmsg("Initial and transition HMM matrices ready", verbose = verbose, level = 1, type = ">>")


  vmsg("Defining z-score distribution templates", verbose = verbose, level = 0, type = ">>")
  mu <- define_z_limits(
    z = d[[z_col]], z_window = z, cn_grid = cn_grid,
    exp_ploidy = exp_ploidy, z_range = z_range, verbose = verbose
  )
  # Array background raises CN=0 z above the zero-signal floor; shift mu[0] accordingly.
  if (0 %in% cn_grid && cn0_background_frac > 0 && "0" %in% names(mu) && "1" %in% names(mu))
    mu["0"] <- mu["0"] + cn0_background_frac * (mu["1"] - mu["0"])

  # sig is the (shared) standard deviation of the z emission across states.
  # It starts at the sample SD of z, with a safety floor of 0.1 to avoid zero/near-zero variance that would blow up log-likelihoods.
  # Assumption: a single sigma is used for every CN state (not state-specific), even though depth
  # variance is in reality expected to scale with copy number; em_hmm_cn() re-floors sig to 1e-3 each M-step.
  sig <- sd(z, na.rm = TRUE)
  if (!is.finite(sig) || sig <= 1e-6) sig <- 0.1
  W <- length(z)
  # Sort windows so chromosomes with the largest z-deviation from the genome mean
  # come first in the EM chain.  This ensures aneuploid chromosomes are processed
  # early in the forward pass, improving mu calibration regardless of chr naming.
  chr_z_dev <- abs(tapply(z, win_df$Chr, mean, na.rm = TRUE) - mean(z, na.rm = TRUE))
  chr_em_ord <- names(sort(chr_z_dev, decreasing = TRUE))
  em_perm     <- order(match(win_df$Chr, chr_em_ord), win_df$Start)
  z_em        <- z[em_perm]
  ll_baf_em   <- ll_baf_matrix[em_perm, , drop = FALSE]
  n_baf_em    <- n_baf[em_perm]
  w_baf_em    <- w_baf[em_perm]
  vmsg("Initial z-score mean by state: %s", verbose = verbose, level = 2, type = ">>", paste0(round(mu, 3), collapse = ", "))
  vmsg("Initial z-score SD: %s", verbose = verbose, level = 2, type = ">>", sig)

  # --- EM loop ---
  vmsg("Starting EM loop", verbose = verbose, level = 0, type = ">>")

  rm_res <- em_hmm_cn(
    cn_grid, mu, K, state_ids, sig, z_em,
    z_only, ll_baf_em, n_baf_em, w_baf_em,
    correct_scale, A, pi0, W, max_iter, verbose,
    update_pi0 = FALSE, baf_only = baf_only
  )

  list2env(rm_res, envir = environment())
  vmsg("Updated z-score mean by state: %s", verbose = verbose, level = 2, type = ">>", paste(sprintf("CN%d: %.3f", cn_grid, mu), collapse = "; "))
  vmsg("Updated z-score SD: %s", verbose = verbose, level = 2, type = ">>", round(sig, 3))

  # If z mean is not from the lowest to the highest follow lower ploidy to higher ploidy
  # It means that user tested unlikely ploidies, in this case, modify cn_grid and run again
  idx <- 0
  mu_tol <- 1e-6  # tolerance to avoid floating-point noise falsely triggering non-monotonicity
  while (any(abs(mu - sort(mu)) > mu_tol) && idx < 10) {
    # Avoid infinite loop
    idx <- idx + 1
    # Identify valid ploidies - for lower ploidies than the expected should have lower mu, higher ploidies than expect should have higher mu
    # If not, the ploidy is not valid and should be removed from the grid, and the estimation should be rerun
    exp_idx <- which(as.numeric(names(mu)) == as.numeric(exp_ploidy))
    mu_exp <- mu[exp_idx]
    # For lower ploidies, keep only those with strictly decreasing mu
    lower_ploidies <- which(as.numeric(names(mu)) < as.numeric(exp_ploidy))
    keep_lower <- lower_ploidies[order(-as.numeric(names(mu)[lower_ploidies]))] # descending order
    if (length(keep_lower) > 0) {
      last_mu <- mu_exp
      valid_lower <- c()
      for (ii in keep_lower) {
        if (mu[ii] < last_mu) {
          valid_lower <- c(valid_lower, ii)
          last_mu <- mu[ii]
        }
      }
      lower_idx <- valid_lower
    } else {
      lower_idx <- integer(0)
    }
    # For higher ploidies, keep only those with strictly increasing mu
    higher_ploidies <- which(as.numeric(names(mu)) > as.numeric(exp_ploidy))
    keep_higher <- higher_ploidies[order(as.numeric(names(mu)[higher_ploidies]))]
    if (length(keep_higher) > 0) {
      last_mu <- mu_exp
      valid_higher <- c()
      for (ii in keep_higher) {
        if (mu[ii] > last_mu) {
          valid_higher <- c(valid_higher, ii)
          last_mu <- mu[ii]
        }
      }
      higher_idx <- valid_higher
    } else {
      higher_idx <- integer(0)
    }
    valid_idx <- c(lower_idx, exp_idx, higher_idx)
    valid_idx <- sort(valid_idx)
    cn_grid <- cn_grid[valid_idx]
    mu <- define_z_limits(
      z = d[[z_col]], z_window = z, cn_grid = cn_grid,
      exp_ploidy = exp_ploidy, z_range_out = FALSE, verbose = verbose
    ) # redefine mu with the new cn_grid, but without z_range to avoid changing the limits too much and keep the same order of the ploidies, which is already checked in the previous steps
    if (0 %in% cn_grid && cn0_background_frac > 0 && "0" %in% names(mu) && "1" %in% names(mu))
      mu["0"] <- mu["0"] + cn0_background_frac * (mu["1"] - mu["0"])
    K <- length(cn_grid)
    ll_baf_matrix <- ll_baf_matrix[, valid_idx]
    ll_baf_em     <- ll_baf_matrix[em_perm, , drop = FALSE]
    pi0 <- pi0[valid_idx]
    state_ids <- state_ids[valid_idx]
    A <- A[valid_idx, valid_idx]
    vmsg("Some ploidies were removed due to non-monotonic z means", verbose = verbose, level = 2, type = ">>")
    vmsg("Rerunning EM with updated cn_grid", verbose = verbose, level = 1, type = ">>")
    rm_res <- em_hmm_cn(
      cn_grid, mu, K, state_ids, sig, z_em,
      z_only, as.matrix(ll_baf_em), n_baf_em, w_baf_em,
      correct_scale, as.matrix(A), pi0, W, max_iter, verbose,
      update_pi0 = FALSE, baf_only = baf_only
    )
    list2env(rm_res, envir = environment())
    vmsg("Updated z-score mean by state: %s", verbose = verbose, level = 2, type = ">>", paste(sprintf("CN%d: %.3f", cn_grid, mu), collapse = "; "))
    vmsg("Updated z-score SD: %s", verbose = verbose, level = 2, type = ">>", sprintf("%.3f", sig))
  }

  # Post-EM collapsed-state pruning (single pass, no re-run).
  # After EM convergence, adjacent CN means that differ by less than
  # collapse_frac * sigma are functionally indistinguishable to the HMM decoder
  # and will cause arbitrary state flipping.  Remove the higher CN state in each
  # such pair (it is less likely to be the true signal in a predominantly diploid
  # sample).  This only fires on collapsed states; a normal diploid sample with a
  # genuine z_range will never trigger it.
  # Collapsed-state threshold: only drop states whose means are virtually identical
  # (< 0.10 * sigma apart). This catches the degenerate case where EM collapses
  # two means to the same value (gap ~ 0.02sigma) without touching states that are
  # genuinely close due to low z_range (gap ~ 0.3sigma is real signal).
  # The exp_ploidy state is never dropped regardless of gap.
  collapse_frac <- 0.10  # min gap as fraction of sigma; only truly degenerate pairs
  if (K > 1) {
    mu_sorted   <- sort(mu)  # ascending: CN1 < CN2 < ...
    gaps        <- diff(mu_sorted)  # gap[i] = mu[i+1] - mu[i]
    keep_names  <- names(mu_sorted)
    drop_mask   <- rep(FALSE, K)
    exp_name    <- as.character(round(exp_ploidy))
    for (gi in seq_along(gaps)) {
      if (!drop_mask[gi] && abs(gaps[gi]) < collapse_frac * sig) {
        # Never drop the exp_ploidy state
        candidate <- keep_names[gi + 1]
        if (candidate == exp_name) next
        drop_mask[gi + 1] <- TRUE
        vmsg("Collapsed CN state dropped: CN%s (|mu_gap|=%.4f < %.4f*sigma=%.4f)",
          verbose = verbose, level = 2, type = ">>",
          candidate, abs(gaps[gi]), collapse_frac, sig)
      }
    }
    if (any(drop_mask)) {
      keep_state_names <- keep_names[!drop_mask]
      keep_state_idx   <- which(names(mu) %in% keep_state_names)
      mu         <- mu[keep_state_idx]
      cn_grid    <- cn_grid[keep_state_idx]
      K          <- length(cn_grid)
      state_ids  <- state_ids[keep_state_idx]
      pi0        <- pi0[keep_state_idx]
      pi0        <- pi0 / sum(pi0)
      A          <- A[keep_state_idx, keep_state_idx]
      ll_baf_matrix <- ll_baf_matrix[, keep_state_idx, drop = FALSE]
      vmsg("Updated cn_grid after collapse pruning: %s", verbose = verbose, level = 2, type = ">>",
        paste(sprintf("CN%d: %.3f", cn_grid, mu), collapse = "; "))
    }
  }

  vmsg("EM loop complete", verbose = verbose, level = 1, type = ">>")

  # Recompute ll_em in ORIGINAL window order from converged mu/sig.
  # The EM ran on sorted windows; the decoder needs original order for correct chr labels.
  # Per-window sigma inflation for all-hom windows (w_baf == 0) reduces z-emission
  # discriminability in the forward-backward pass.
  # When CN=0 is in the grid and z_only = FALSE, BAF noise on nullisomic windows
  # can inflate w_baf, weakening the CN=0 z-emission and boosting other states.
  # Force w_baf = 0 for windows whose z falls below the midpoint of mu[0] and mu[1].
  w_baf_dec <- w_baf
  if (!z_only && 0 %in% cn_grid && "0" %in% names(mu) && "1" %in% names(mu)) {
    cn0_mid <- mu["0"] + 0.5 * (mu["1"] - mu["0"])
    w_baf_dec[z <= cn0_mid] <- 0
  }
  sig_dec <- rep(sig, W)
  if (!z_only && !baf_only && hom_z_sigma_inflate > 1)
    sig_dec <- ifelse(w_baf_dec == 0, sig * hom_z_sigma_inflate, sig)
  # Windows with no valid BAF values carry no BAF information and must be excluded
  # from BAF scoring rather than treated as a log-likelihood of 0 (i.e. probability 1).
  no_baf <- !z_only & (n_baf == 0)
  if (!z_only && any(no_baf)) w_baf_dec[no_baf] <- 0
  ll_em <- matrix(NA_real_, W, K, dimnames = list(NULL, state_ids))
  for (k in seq_len(K)) {
    c_k  <- cn_grid[k]
    llz  <- dnorm(z, mean = mu[as.character(c_k)], sd = sig_dec, log = TRUE)
    llz[is.nan(llz)] <- 0
    if (z_only) {
      ll_em[, k] <- llz
    } else if (baf_only) {
      llb <- ll_baf_matrix[, k]
      if (correct_scale) llb <- llb / n_baf
      llb[no_baf] <- 0  # excluded: neutral emission, equal across all states
      ll_em[, k] <- llb
    } else {
      llb <- ll_baf_matrix[, k]
      if (correct_scale) llb <- llb / n_baf
      llb[no_baf] <- 0  # excluded: BAF term dropped, w_baf_dec forced to 0 above
      # Tempered log-product: z always contributes fully, BAF contributes a fraction of its evidence.
      ll_em[, k] <- llz + w_baf_dec * llb
    }
  }

  # Per-window z and BAF log-likelihood matrices (used for gating and pi0 blending below)
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

  # Two-signal consistency gating (chromosome level)
  # Compare the chromosome-level z consensus vs BAF consensus for each chromosome.
  # Per-window z-scores are too noisy to gate individually (especially with wide CN grids),
  # so we aggregate across all windows on each chromosome before deciding.
  #
  # - If the chromosome-level mean BAF weight is low (< 0.1): BAF is sparse/absent,
  #   z already dominates; no change.
  # - If BAF weight is meaningful AND chromosome z and BAF agree on the best CN state:
  #   reinforce with a log(3) * mean_w_baf bonus on the agreed state for all windows.
  # - If BAF weight is meaningful AND they disagree: BAF is contradicting the depth
  #   signal at the chromosome level. Suppress BAF for all windows on that chromosome
  #   (replace ll_em with z-only emission).
  if (!z_only && !baf_only && K > 1) {
    baf_gate <- 0.1  # minimum mean BAF weight for a chromosome to be gated
    gate_log <- character(0)

    for (chr_i in unique(win_df$Chr)) {
      idx         <- which(win_df$Chr == chr_i)
      mean_w_baf  <- mean(w_baf_dec[idx], na.rm = TRUE)

      if (mean_w_baf < baf_gate) next  # BAF sparse: z already dominates, no change

      # Chromosome-level consensus: sum log-likelihoods across windows
      chr_z_sum   <- colSums(llz_mat[idx, , drop = FALSE], na.rm = TRUE)
      chr_baf_sum <- colSums(llb_mat[idx, , drop = FALSE], na.rm = TRUE)
      best_z_chr  <- which.max(chr_z_sum)
      best_baf_chr <- which.max(chr_baf_sum)

      if (best_z_chr == best_baf_chr) {
        # Signals agree at chromosome level: reinforce the agreed state
        ll_em[idx, best_z_chr] <- ll_em[idx, best_z_chr] + log(3) * mean_w_baf
        gate_log <- c(gate_log, sprintf("%s:agree(CN%d)", chr_i, cn_grid[best_z_chr]))
      } else {
        # Signals disagree at chromosome level: suppress BAF, use z-only emission
        ll_em[idx, ] <- llz_mat[idx, , drop = FALSE]
        gate_log <- c(gate_log, sprintf("%s:disagree(z=CN%d,baf=CN%d)",
          chr_i, cn_grid[best_z_chr], cn_grid[best_baf_chr]))
      }
    }
    message(sprintf("[GATE] W=%d K=%d baf_gate=%.2f | %s",
      W, K, baf_gate, if (length(gate_log)) paste(gate_log, collapse = " ") else "none"))
  }

  vmsg("Decoding CN path", verbose = verbose, level = 0, type = ">>")
  # EM-learned A can have degenerate uniform rows for states never visited:
  # a missing-state can "jump for free" to the best future window, giving it a
  # ~log(K) backward bonus that can flip window 1's call.  Using the prior A
  # (strong diagonal, rare off-diagonal) costs ~log(K/1e-3) ≈ 8 for any
  # off-diagonal transition, making a degenerate state unable to overcome the
  # correct state's backward advantage regardless of emission noise.
  A_dec <- matrix(1e-3, K, K)
  diag(A_dec) <- transition_jump
  A_dec <- A_dec / rowSums(A_dec)
  # Decode per chromosome. EM parameters (mu, sig) are global so z-scores are
  # comparable across chromosomes. pi0 per chromosome is derived from the sum of
  # BAF log-likelihoods across all windows of that chromosome: this anchors the
  # forward pass to the ploidy the BAF data supports and prevents a high z-score
  # in the first window from locking the whole chromosome into the wrong CN state.
  gamma_dec <- matrix(0, W, K, dimnames = list(NULL, colnames(ll_em)))
  for (chr_i in unique(win_df$Chr)) {
    idx <- which(win_df$Chr == chr_i)
    logpi0_chr <- NULL  # default: uniform
    if (!z_only) {
      # pi0 is blended between z-score consensus and BAF consensus based on mean
      # BAF weight for this chromosome. When BAF weight is low (sparse heterozygous
      # markers), the BAF log-likelihoods are unreliable as a prior, and the z-score
      # sum across windows is a better anchor for the starting state distribution.
      mean_w_baf_chr <- mean(w_baf_dec[idx], na.rm = TRUE)
      # Guard: mean(na.rm=TRUE) of an all-NA vector returns NaN, not 0.
      # Treat NaN mean as 0 (no BAF data) so the blend below doesn't propagate NaN.
      if (is.nan(mean_w_baf_chr)) mean_w_baf_chr <- 0

      # z-score based pi0: sum of z log-likelihoods across chromosome windows
      chr_z_sum <- colSums(llz_mat[idx, , drop = FALSE], na.rm = TRUE)
      chr_z_sum <- chr_z_sum - max(chr_z_sum)
      lse_z     <- logsumexp(chr_z_sum)
      logpi0_z  <- if (is.finite(lse_z)) chr_z_sum - lse_z else log(pi0 / sum(pi0))

      # BAF-based pi0: sum of raw BAF log-likelihoods (original behaviour)
      chr_baf_sum <- colSums(ll_baf_matrix[idx, , drop = FALSE], na.rm = TRUE)
      chr_baf_sum <- chr_baf_sum - max(chr_baf_sum)
      lse_baf     <- logsumexp(chr_baf_sum)
      logpi0_baf  <- if (is.finite(lse_baf)) chr_baf_sum - lse_baf else log(pi0 / sum(pi0))

      # Blend: weight BAF pi0 by mean_w_baf_chr, z pi0 by (1 - mean_w_baf_chr)
      # so chromosomes with little BAF data are anchored by depth, not BAF noise.
      logpi0_chr <- log(
        (1 - mean_w_baf_chr) * exp(logpi0_z) + mean_w_baf_chr * exp(logpi0_baf)
      )
      # Numeric stability: subtract max.  If all entries are -Inf or NaN (can
      # happen when both logpi0_z and logpi0_baf are degenerate), fall back to
      # the uniform/exp_ploidy prior rather than propagating NaN into fb_smooth.
      max_lp <- max(logpi0_chr)
      logpi0_chr <- if (is.finite(max_lp)) logpi0_chr - max_lp else log(pi0 / sum(pi0))
    } else {
      # z_only: fix per-chromosome pi0 at the exp_ploidy prior; deviates only if
      # within-chromosome evidence is strong enough to overcome it.
      logpi0_chr <- log(pi0 / sum(pi0))
    }
    gamma_dec[idx, ] <- fb_smooth(ll_em[idx, , drop = FALSE], log(A_dec), logpi0_chr)
  }
  if (!is.matrix(gamma_dec) || ncol(gamma_dec) == 0L)
    stop(sprintf("Internal error: gamma_dec is not a valid matrix after fb_smooth (dim=%s, K=%d). Please report this.", paste(dim(gamma_dec), collapse="x"), K))
  # which.max on a row of all-NaN returns integer(0), and apply() then returns a
  # list (mixed lengths).  Use a safe wrapper that falls back to state 1 (exp_ploidy)
  # when all posterior values are non-finite, so the result is always an integer vector.
  safe_which_max <- function(x) { idx <- which.max(x); if (length(idx) == 0L) 1L else idx }
  cn_call  <- cn_grid[vapply(seq_len(nrow(gamma_dec)), function(i) safe_which_max(gamma_dec[i, ]), integer(1L))]
  post_max <- apply(gamma_dec, 1, max)
  # CN=0 windows carry only background array signal; treat as missing like sequencing gaps.
  cn0_idx <- cn_call == 0
  cn_call[cn0_idx]  <- NA_integer_
  post_max[cn0_idx] <- NA_real_
  z[cn0_idx]        <- NA_real_
  if (!z_only) w_baf[cn0_idx] <- NA_real_

  # Post-processing: override mis-called chromosomes when BAF is absent/sparse
  # and the z evidence is too weak to support a call different from the expected ploidy.
  #
  # Two tiers:
  #   Tier 1 (all-hom): every window has w_baf == 0; apply original threshold.
  #   Tier 2 (sparse BAF): mean_w_baf < baf_gate AND called CN != mode_CN AND
  #           |chr_z_mean - mu[mode_CN]| < hom_z_sigma_inflate * sig.
  #           These are chromosomes where a slight z fluctuation drove a CN call
  #           that has no BAF support; override back to the sample-mode CN.
  if (!z_only && !baf_only && hom_z_sigma_inflate > 1) {
    baf_gate_post   <- 0.1  # same gate as in gating block above
    het_idx         <- w_baf > 0
    mode_cn         <- as.integer(if (any(het_idx)) mode(cn_call[het_idx]) else round(exp_ploidy))
    mode_col        <- which(cn_grid == mode_cn)
    baseline_mu_val <- if (length(mode_col) == 1) mu[mode_col] else mu[which.min(abs(cn_grid - mode_cn))]
    threshold       <- sig * hom_z_sigma_inflate

    for (chr_i in unique(win_df$Chr)) {
      idx            <- which(win_df$Chr == chr_i)
      mean_w_baf_chr <- mean(w_baf[idx], na.rm = TRUE)
      chr_z_mean     <- mean(z[idx], na.rm = TRUE)
      chr_cn_mode    <- as.integer(mode(cn_call[idx]))

      # Only act when chromosome CN differs from sample mode and z is within threshold
      if (!is.finite(chr_z_mean)) next
      if (chr_cn_mode == mode_cn) next
      if (abs(chr_z_mean - baseline_mu_val) >= threshold) next  # genuine z deviation

      # Override regardless of BAF weight: the z deviation is within noise of the
      # sample mode, so the call is not supported by depth signal. When BAF is also
      # present and was not strong enough to prevent this call, it too is within noise.
      cn_call[idx] <- mode_cn
      if (length(mode_col) == 1)
        post_max[idx] <- gamma_dec[idx, mode_col]
      message(sprintf("[POST-GATE] %s: CN%d→CN%d (mean_w_baf=%.3f, z_dev=%.3f < %.3f*sig)",
        chr_i, chr_cn_mode, mode_cn, mean_w_baf_chr,
        abs(chr_z_mean - baseline_mu_val), hom_z_sigma_inflate))
    }
  }

  # CN_reliability: joint confidence combining four components:
  #   post_max   - HMM posterior certainty for the called state
  #   baf_info   - BAF informativeness: 1 - normalised Shannon entropy of a
  #                10-bin BAF histogram on [0,1] per window.  Values are clipped
  #                to [0,1] before binning so ratio columns are handled safely.
  #                Uniform/dispersed BAF -> 0; peaks at heterozygous positions -> 1.
  #                Windows with no BAF data use z_no_baf_scale as the floor.
  #   zbaf_agree - per-window z vs BAF agreement: 1 - total variation distance
  #                between the softmax of llz_mat and softmax of llb_mat for that
  #                window (both available from the gating block). 1 = identical
  #                distributions; 0 = completely opposite. Windows with w_baf==0
  #                or when gating is skipped (z_only/baf_only/K==1) default to 1.
  #   w_baf_norm - existing BAF-weight scaling (unchanged)
  w_baf_use <- if (z_only) rep(0, W) else w_baf
  w_baf_norm <- pmax(z_no_baf_scale, pmin(1, w_baf_use / baf_weight))

  # --- baf_info: normalised BAF entropy per window ---
  baf_entropy_bins <- 10L
  baf_breaks <- seq(0, 1, length.out = baf_entropy_bins + 1L)
  baf_info <- vapply(seq_len(W), function(i) {
    if (isTRUE(is.na(w_baf_use[i])) || isTRUE(w_baf_use[i] == 0)) return(z_no_baf_scale)
    bv <- baf_list[[i]]
    bv <- pmax(0, pmin(1, bv[is.finite(bv)]))  # clip to [0,1]; handles ratio columns
    if (length(bv) == 0L) return(z_no_baf_scale)
    counts <- tabulate(findInterval(bv, baf_breaks, rightmost.closed = TRUE),
                       nbins = baf_entropy_bins)
    s <- sum(counts)
    if (s == 0L) return(z_no_baf_scale)
    p <- counts[counts > 0L] / s
    H <- -sum(p * log(p))           # Shannon entropy (nats)
    H_max <- log(baf_entropy_bins)  # entropy of uniform distribution
    max(0, min(1, 1 - H / H_max))  # 0 = uniform/dispersed, 1 = perfectly peaked
  }, numeric(1L))

  # --- zbaf_agree: per-window z vs BAF posterior agreement ---
  # Uses llz_mat and llb_mat computed in the gating block above.
  # Softmax each row independently, then take 1 - total variation distance.
  if (!z_only && !baf_only && K > 1 && exists("llz_mat") && exists("llb_mat")) {
    softmax_row <- function(x) { x <- x - max(x); e <- exp(x); e / sum(e) }
    zbaf_agree <- vapply(seq_len(W), function(i) {
      if (isTRUE(is.na(w_baf_use[i])) || isTRUE(w_baf_use[i] == 0)) return(1)
      pz <- softmax_row(llz_mat[i, ])
      pb <- softmax_row(llb_mat[i, ])
      if (!all(is.finite(pz)) || !all(is.finite(pb))) return(1)
      1 - 0.5 * sum(abs(pz - pb))  # 1 - total variation distance; in [0,1]
    }, numeric(1L))
  } else {
    zbaf_agree <- rep(1, W)
  }

  CN_reliability <- post_max * baf_info * zbaf_agree * w_baf_norm
  CN_reliability[!is.finite(CN_reliability)] <- NA_real_

  vmsg("CN path decoded", verbose = verbose, level = 1, type = ">>")

  # Prepare output
  post_df <- as.data.frame(gamma_dec)
  names(post_df) <- paste0("post_CN", cn_grid)
  result <- cbind(
    data.frame(
      Sample = sample_id,
      Chr = win_df$Chr,
      WindowID = win_df$WindowID,
      Start = win_df$Start,
      End = win_df$End,
      n_snps = win_df$n_snps,
      n_het = n_het_window,
      z = z,
      w_baf = if (z_only) 0 else w_baf,
      CN_call = cn_call,
      CN_reliability = CN_reliability,
      post_max = post_max,
      stringsAsFactors = FALSE
    ),
    post_df
  )
  # Re-insert CN=0 windows removed by the pre-pass with NA values.
  if (!is.null(cn0_win_df) && nrow(cn0_win_df) > 0) {
    cn0_post <- as.data.frame(matrix(NA_real_, nrow = nrow(cn0_win_df), ncol = ncol(post_df)))
    names(cn0_post) <- names(post_df)
    cn0_rows <- cbind(
      data.frame(
        Sample = sample_id, Chr = cn0_win_df$Chr, WindowID = cn0_win_df$WindowID,
        Start = cn0_win_df$Start, End = cn0_win_df$End, n_snps = cn0_win_df$n_snps,
        n_het = NA_integer_, z = NA_real_, w_baf = NA_real_,
        CN_call = NA_integer_, CN_reliability = NA_real_, post_max = NA_real_,
        stringsAsFactors = FALSE
      ),
      cn0_post
    )
    result <- result[order(result$Chr, result$Start), ]
    result <- rbind(result, cn0_rows)
    result <- result[order(result$Chr, result$Start), ]
  }
  params <- list(
    cn_grid = cn_grid,
    distribution = if ((isTRUE(nQuack_EM) || isTRUE(nQuack_exp_ploidy)) && !is.null(nq_dist) && !is.na(nq_dist)) {
                     nq_pretty <- c(
                       "normal"                = "Normal",
                       "normal-uniform"        = "Normal + Uniform",
                       "beta"                  = "Beta",
                       "beta-uniform"          = "Beta + Uniform",
                       "beta-binomial"         = "Beta-Binomial",
                       "beta-binomial-uniform" = "Beta-Binomial + Uniform"
                     )
                     pretty <- unname(nq_pretty[tolower(nq_dist)])
                     if (is.na(pretty)) paste0("nQuack - ", nq_dist) else paste0("nQuack - ", pretty)
                   } else selected_model$best$dist,
    mu = if (baf_only) NA else mu,
    sigma = if (baf_only) NA else sig,
    A = A,
    pi0 = pi0,
    bins = M,
    bw = if (isTRUE(nQuack_EM)) NA_real_ else selected_model$best$bw,
    loglik = ll_hist[length(ll_hist)],
    z_range = z_range,
    het_quantile = het_quantile,
    baf_weight = baf_weight,
    transition_jump = transition_jump,
    z_only = z_only,
    baf_only = baf_only,
    exp_ploidy = exp_ploidy,
    rm_outliers = rm_outliers,
    outlier_alpha = outlier_alpha,
    segment_zscore = segment_zscore,
    snps_per_window = snps_per_window,
    min_snps_per_window = min_snps_per_window,
    add_uniform = add_uniform,
    uniform_weight = if (isTRUE(nQuack_EM) && !is.null(nq_uniform_weight)) nq_uniform_weight else selected_model$best$uniform_weight,
    cn0_background_frac = cn0_background_frac
  )

  if (!is.null(result) && !is.null(d)) {
    map_df <- result[, c("Chr", "WindowID", "w_baf", "CN_call", "post_max", "CN_reliability")]
    names(map_df)[names(map_df) == "WindowID"] <- ".__w__"
    d <- merge(d, map_df, by = c("Chr", ".__w__"), all.x = TRUE, sort = FALSE)
    d <- d[order(match(seq_len(nrow(d)), as.integer(rownames(d)))), ]
    rownames(d) <- NULL
  }

  rm_markers <- d$MarkerName[which(d$CN_call != mode(d$CN_call))]

  if (rerun_overall_ploidy) {
    if(length(rm_markers) > 0){
    vmsg("Rerunning overall ploidy estimation with markers that disagree with the initial CN call removed", verbose = verbose, level = 0, type = ">>")

    result <- hmm_estimate_CN(
      qploidy_standarize_result = qploidy_standarize_result,
      sample_id = sample_id,
      data = data,
      geno.pos = geno.pos,
      use_values = use_values,
      chr = chr,
      reflect = reflect,
      add_uniform = add_uniform,
      segment_zscore = segment_zscore,
      snps_per_window = snps_per_window,
      min_snps_per_window = min_snps_per_window,
      cn_grid = cn_grid,
      M = M,
      max_iter = max_iter,
      het_quantile = het_quantile, # increase this value to reduce the weight of baf when few hets
      baf_weight = baf_weight,
      z_range = z_range,
      transition_jump = transition_jump, # decrease this value if you think there changes in CN is likely
      initial_prob = initial_prob, # Initial probability for the best CN state in the initial state distribution (pi0). Default 0.15. Sets the prior probability for the expected ploidy (or best CN from BAF model) at the first window; remaining probability is distributed uniformly across other states. If the best CN is not found, pi0 is uniform across all states.
      z_only = z_only,
      baf_only = baf_only,
      verbose = verbose,
      exp_ploidy = NULL,
      rm_outliers = rm_outliers,
      outlier_alpha = outlier_alpha,
      selected_model = NULL,
      dists = dists,
      bw_grid = bw_grid,
      add_uniform_grid = add_uniform_grid,
      uniform_weight_grid = uniform_weight_grid,
      param_count = param_count,
      count_grid_as_params = count_grid_as_params,
      correct_scale = correct_scale,
      min_het_frac = min_het_frac,
      het_range = het_range,
      dosage_threshold = dosage_threshold,
      z_no_baf_scale = z_no_baf_scale,
      hom_z_sigma_inflate = hom_z_sigma_inflate,
      cn0_background_frac = cn0_background_frac,
      nQuack_EM       = nQuack_EM,
      nQuack_cores    = nQuack_cores,
      nQuack_dists      = nQuack_dists,
      nQuack_exp_ploidy = nQuack_exp_ploidy,
      rerun_overall_ploidy = FALSE,
      recycled_obj_rerun_overall_ploidy = list(rm_mks = rm_markers,
                                               recycled_d = keep_d,
                                               recycled_baf_col = keep_baf_col,
                                               recycled_z_col = keep_z_col,
                                               recycled_cn_grid = keep_cn_grid,
                                               recycled_nq_exp_ploidy = if (isTRUE(nQuack_EM) && !is.null(nq_ploidy)) nq_ploidy else NULL,
                                               recycled_nq_winner     = if (isTRUE(nQuack_EM) && !is.null(nq_winner)) nq_winner else NULL,
                                               recycled_nq_dist       = if (isTRUE(nQuack_EM) && !is.null(nq_dist))   nq_dist   else NULL)
    )

    # The rerun has its own pre-pass, but apply the outer call's cn0_win_df as
    # a hard override to guarantee CN=0 chromosomes are NA in the final output.
    if (!is.null(cn0_win_df) && nrow(cn0_win_df) > 0) {
      cn0_chr_win <- paste(cn0_win_df$Chr, cn0_win_df$WindowID, sep = "__")
      bw_key <- paste(result$by_window$Chr, result$by_window$WindowID, sep = "__")
      cn0_bw  <- bw_key %in% cn0_chr_win
      result$by_window$CN_call[cn0_bw]       <- NA_integer_
      result$by_window$z[cn0_bw]             <- NA_real_
      result$by_window$w_baf[cn0_bw]         <- NA_real_
      result$by_window$post_max[cn0_bw]      <- NA_real_
      result$by_window$CN_reliability[cn0_bw] <- NA_real_
      bm_key <- paste(result$by_marker$Chr, result$by_marker$`.__w__`, sep = "__")
      cn0_bm  <- bm_key %in% cn0_chr_win
      result$by_marker$CN_call[cn0_bm]       <- NA_integer_
      result$by_marker$w_baf[cn0_bm]         <- NA_real_
      result$by_marker$post_max[cn0_bm]      <- NA_real_
      result$by_marker$CN_reliability[cn0_bm] <- NA_real_
    }
    return(result)
    } else {
      vmsg(paste("Rerun option chosen but no CN different than mode", mode(d$CN_call),"was found. Rerun not needed."), verbose = verbose, level = 0, type = ">>")
    }
  }

  vmsg("Done!", verbose = verbose, level = 1, type = ">>")

  return(structure(list(by_window = result, by_marker = d, params = params), class = "hmm_CN"))
}

#' Run hmm_estimate_CN in parallel for multiple samples (using parLapply)
#'
#' Runs copy-number HMM for each sample in a user-defined vector (or all samples if "all" is specified).
#' Returns a combined data.frame with results for all samples.
#'
#' @param qploidy_standarize_result An object of class qploidy_standardization as returned by standardize().
#' @param data Optional data.frame with raw input data (MarkerName, SampleName, ratio, R) if qploidy_standarize_result is not provided.
#' @param use_values Character vector of column names to use for CN estimation. If qploidy_standarize_result is provided, these should match the columns in qploidy_standarize_result$data (e.g., c("BAF", "zscore")). If qploidy_standarize_result is NULL, use_values must be c("ratio", "R") and these columns must be present in the provided 'data'.
#' @param geno.pos Optional data.frame with marker genomic positions (MarkerName, Chromosome, Position) if qploidy_standarize_result is not provided.
#' @param sample_ids Character vector of sample IDs to analyze, or "all" for all samples in the data.
#' @param n_cores Number of cores to use (default: 2).
#' @param parallel_type Character. Parallel backend to use: \code{"FORK"}, \code{"PSOCK"}, or \code{"auto"} (default). \code{"auto"} selects \code{"FORK"} on Unix/macOS (faster; workers inherit the parent environment automatically) and \code{"PSOCK"} on Windows (the only option available there). Use \code{"PSOCK"} explicitly if you need cross-platform reproducibility or are debugging worker crashes.
#' @param ... Additional arguments passed to hmm_estimate_CN (e.g., chr, snps_per_window, etc).
#' @return A data.frame with results for all samples, as returned by hmm_estimate_CN$by_window, combined.
#'
#' @importFrom parallel makeCluster parLapply stopCluster clusterExport
#' @importFrom dplyr bind_rows
#'
#' @export
hmm_estimate_CN_multi <- function(qploidy_standarize_result = NULL,
                                  data = NULL,
                                  geno.pos = NULL,
                                  use_values = c("BAF", "zscore"),
                                  sample_ids = "all",
                                  n_cores = 2,
                                  parallel_type = "auto",
                                  ...) {
  # sanity check: require either a standardization object or raw data + geno.pos
  if (!is.null(qploidy_standarize_result)) {
    if (!inherits(qploidy_standarize_result, "qploidy_standardization")) {
      stop("Input must be a qploidy_standardization object as returned by standardize().")
    }
    all_samples <- unique(as.data.frame(qploidy_standarize_result$data)$SampleName)
  } else {
    if (is.null(data) || is.null(geno.pos)) {
      stop("When qploidy_standarize_result is NULL, both 'data' and 'geno.pos' must be provided.")
    }
    req_data <- c("MarkerName", "SampleName", "ratio", "R")
    req_gp <- c("MarkerName", "Chromosome", "Position")
    miss_d <- setdiff(req_data, names(data))
    miss_gp <- setdiff(req_gp, names(geno.pos))
    if (length(miss_d)) stop(paste("'data' is missing columns:", paste(miss_d, collapse = ", ")))
    if (length(miss_gp)) stop(paste("'geno.pos' is missing columns:", paste(miss_gp, collapse = ", ")))
    if (!identical(use_values, c("ratio", "R"))) {
      stop("When qploidy_standarize_result is NULL, use_values must be c(\"ratio\", \"R\").")
    }
    all_samples <- unique(data$SampleName)
  }

  # resolve parallel backend
  allowed_types <- c("auto", "FORK", "PSOCK")
  if (!parallel_type %in% allowed_types) {
    stop(sprintf("'parallel_type' must be one of: %s.", paste(allowed_types, collapse = ", ")))
  }
  if (parallel_type == "auto") {
    parallel_type <- if (.Platform$OS.type == "windows") "PSOCK" else "FORK"
  }

  if (identical(sample_ids, "all")) {
    sample_ids <- all_samples
  } else {
    sample_ids <- intersect(sample_ids, all_samples)
    if (length(sample_ids) == 0) stop("No valid sample IDs found in input data.")
  }

  # capture dots ONCE
  dots <- list(...)

  cl <- makeCluster(n_cores, type = parallel_type)
  on.exit(stopCluster(cl), add = TRUE)

  if (parallel_type == "PSOCK") {
    # PSOCK workers are blank R processes: must reload packages and export symbols
    clusterEvalQ(cl, {
      options(warn = 1)
      library(methods)
      library(Qploidy2)
      NULL
    })
    clusterExport(cl,
      varlist = c(
        "worker", "hmm_estimate_CN", "generate_baf_template",
        "baf_log_likelihood", "logsumexp", "viterbi"
      ),
      envir = environment()
    )
  }
  # FORK workers inherit everything from the parent process — no export needed

  results_list <- parLapply(
    cl, sample_ids, worker,
    qploidy_standarize_result, dots, data, geno.pos, use_values
  )

  # Re-issue warnings collected inside workers on the main session
  for (worker_res in results_list) {
    for (w in worker_res$warnings) warning(w, call. = FALSE)
  }
  results_list <- lapply(results_list, `[[`, "result")

  parameters <- lapply(results_list, function(x) x$params)
  names(parameters) <- sample_ids
  by_window <- lapply(results_list, function(x) x$by_window)
  by_marker <- lapply(results_list, function(x) x$by_marker)

  by_window <- Filter(Negate(is.null), by_window)
  if (length(by_window) == 0) stop("No results returned for any sample.")

  by_window <- bind_rows(by_window)
  by_marker <- bind_rows(by_marker)

  idx <- which(colnames(by_window) == "post_max")
  idx1 <- order(colnames(by_window)[(idx + 1):ncol(by_window)])
  by_window <- by_window[, c(1:idx, idx + idx1)]

  rownames(by_window) <- NULL
  return(structure(list(by_window = by_window, by_marker = by_marker, params_samples = parameters), class = "hmm_CN"))
}


#' Print method for hmm_CN objects
#'
#' Prints a summary of the HMM copy-number estimation results, including key parameters.
#' If the object contains a params_samples list (multi-sample), prints the number of samples.
#' Otherwise, prints details from the params list.
#'
#' @param x An object of class 'hmm_CN'.
#' @param ... Additional arguments (ignored).
#'
#' @method print hmm_CN
#'
#' @export
print.hmm_CN <- function(x, ...) {
  if (!inherits(x, "hmm_CN")) {
    stop("Object is not of class 'hmm_CN'.")
  }
  if (!is.null(x$params_samples)) {
    cat("hmm_CN multi-sample result\n")
    print(count_types(x))
    invisible(x)
    return()
  }
  params <- x$params
  cat("hmm_CN result\n")
  cat("  Copy-number grid:", paste(params$cn_grid, collapse = ", "), "\n")
  cat("  Expected ploidy:", params$exp_ploidy, "\n")
  cat("  Minimum SNPs per window:", params$min_snps_per_window, "\n")
  cat("  Initial state probabilities (pi0):", paste(round(params$pi0, 3), collapse = ", "), "\n")
  if (!isTRUE(params$baf_only)) {
    cat("  Estimated z means per CN:", paste(round(params$mu, 3), collapse = ", "), "\n")
    cat("  Estimated z mean:", mean(x$by_window$z, na.rm = TRUE), "\n")
    cat("  Estimated z sigma:", round(params$sigma, 3), "\n")
  }
  cat("  BAF Emission distribution:", if (!is.null(params$distribution)) params$distribution else "(not specified)", "\n")
  cat("  Final log-likelihood:", params$loglik, "\n")
  # Print range of CN_call values
  if (!is.null(x$by_window) && !is.null(x$by_window$CN_call)) {
    cn_vals <- x$by_window$CN_call
    cat("  Range of CN_call (window-level):", paste(range(as.numeric(cn_vals), na.rm = TRUE), collapse = " - "), "\n")
  }
  invisible(x)
}

#' Write hmm_CN object to three CSV files
#'
#' Writes the three main components of an hmm_CN object (by_window, by_marker, params) to CSV files.
#' The user must provide a file prefix; files will be named <prefix>_by_window.csv, <prefix>_by_marker.csv, <prefix>_params.csv.
#' @param hmm_CN An object of class 'hmm_CN'.
#' @param prefix File prefix for output files (character scalar).
#'
#' @importFrom data.table fwrite
#' @export
write_hmm_CN <- function(hmm_CN, prefix) {
  stopifnot(inherits(hmm_CN, "hmm_CN"))
  stopifnot(is.character(prefix) && length(prefix) == 1)
  fwrite(hmm_CN$by_window, paste0(prefix, "_by_window.csv.gz"), row.names = FALSE, compress = "gzip")
  fwrite(hmm_CN$by_marker, paste0(prefix, "_by_marker.csv.gz"), row.names = FALSE, compress = "gzip")
  # Save params or params_samples as RDS, depending on which is present
  if (!is.null(hmm_CN$params_samples)) {
    saveRDS(list(params_samples = hmm_CN$params_samples), file = paste0(prefix, "_params.rds"))
  } else if (!is.null(hmm_CN$params)) {
    saveRDS(hmm_CN$params, file = paste0(prefix, "_params.rds"))
  } else {
    stop("hmm_CN object must have either 'params' or 'params_samples'.")
  }
}

#' Read hmm_CN object from three files
#'
#' Reads the three main components of an hmm_CN object (by_window, by_marker, params or params_samples) from files.
#' The user must provide the three file paths explicitly.
#' @param by_marker_file Path to the by_marker CSV file.
#' @param by_window_file Path to the by_window CSV file.
#' @param params_file Path to the params RDS file.
#' @importFrom data.table fread
#'
#' @return An object of class 'hmm_CN'.
#' @export
read_hmm_CN <- function(by_marker_file, by_window_file , params_file) {
  by_window <- fread(by_window_file, data.table = FALSE)
  by_marker <- fread(by_marker_file, data.table = FALSE)

  # Input checks for by_window columns
  required_cols_window <- c(
    "Sample", # Identifier for the sample
    "Chr", # Chromosome name
    "WindowID", # Unique identifier for the genomic window
    "Start", # Start position of the window in base pairs
    "End", # End position of the window in base pairs
    "n_snps", # Number of SNPs within the window
    "n_het", # Number of heterozygous SNPs within the window
    "z", # Z-score for the window
    "w_baf", # B-allele frequency HMM weight for the window
    "CN_call", # Copy number call for the window
    "post_max" # Maximum posterior probability for the copy number call
  )
  missing_cols_window <- setdiff(required_cols_window, colnames(by_window))
  if (length(missing_cols_window) > 0) {
    stop(sprintf(
      "by_window file is missing required column(s): %s",
      paste(missing_cols_window, collapse = ", ")
    ))
  }
  # Check for at least one post_CN column
  post_cn_cols <- grep("^post_CN[0-9]+$", colnames(by_window), value = TRUE)
  if (length(post_cn_cols) == 0) {
    stop("by_window file must contain at least one 'post_CN[number]' column (e.g., post_CN1, post_CN2, ...)")
  }

  # Input checks for by_marker columns
  required_cols_marker <- c(
    "Chr", # Chromosome name
    ".__w__", # Window
    "MarkerName", # Name of the marker
    "SampleName", # Identifier for the sample
    "X", # X-coordinate for the marker
    "Y", # Y-coordinate for the marker
    "R", # Intensity value for the marker
    "ratio", # Ratio value for the marker
    "geno", # Genotype call for the marker
    "baf", # B-allele frequency for the marker
    "z", # Z-score for the marker
    "Position", # Genomic position of the marker
    "outlier", # Outlier status for the marker
    "w_baf", # Weighted B-allele frequency for the marker
    "CN_call", # Copy number call for the marker
    "post_max" # Maximum posterior probability for the copy number call
  )
  missing_cols_marker <- setdiff(required_cols_marker, colnames(by_marker))
  if (length(missing_cols_marker) > 0) {
    stop(sprintf(
      "by_marker file is missing required column(s): %s",
      paste(missing_cols_marker, collapse = ", ")
    ))
  }

  # Or make it flexible to handle both local files and URLs
  if (grepl("^http", params_file)) {
    params_obj <- readRDS(gzcon(url(params_file)))
  } else {
    params_obj <- readRDS(params_file)
  }
  # Multi-sample: params_obj is a list with params_samples, or just the params_samples list itself
  required_params <- c(
    "cn_grid", "distribution", "mu", "sigma", "A", "pi0", "bins", "bw", "loglik", "z_range",
    "het_quantile", "baf_weight", "transition_jump", "z_only", "exp_ploidy", "rm_outliers",
    "outlier_alpha", "segment_zscore", "snps_per_window", "min_snps_per_window", "add_uniform", "uniform_weight"
  )

  check_params_list <- function(param_list) {
    missing <- setdiff(required_params, names(param_list))
    length(missing) > 0
  }

  if (is.list(params_obj) && (!is.null(params_obj$params_samples) || (is.null(names(params_obj)) && all(sapply(params_obj, function(x) is.list(x) && all(required_params %in% names(x))))))) {
    params_samples <- if (!is.null(params_obj$params_samples)) params_obj$params_samples else params_obj
    # Check all internal lists
    check <- sapply(params_samples, check_params_list)
    if (any(check)) {
      warning(paste("Sample/s", names(check)[which(check)], " don't have all parameters information"))
    }
    structure(list(by_window = by_window, by_marker = by_marker, params_samples = params_samples), class = "hmm_CN")
  } else if (is.list(params_obj) && !is.null(names(params_obj)) && all(required_params %in% names(params_obj))) {
    check_params_list(params_obj)
    structure(list(by_window = by_window, by_marker = by_marker, params = params_obj), class = "hmm_CN")
  } else {
    stop("params_file does not contain a recognizable params or params_samples object with all required items.")
  }
}