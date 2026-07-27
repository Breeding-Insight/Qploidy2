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
#'   \code{X} (reference allele count), \code{Y} (alternative allele count),
#'   \code{ratio}, and \code{geno}.
#' @param min_depth Integer. Minimum total read depth (\code{R}) required for
#'   a site to be retained. Default is \code{10}.
#' @param max_geno Integer. Maximum genotype dosage considered homozygous
#'   (inclusive). Sites where \code{geno == 0} or \code{geno == max_geno} are
#'   filtered as homozygous. Default is \code{4} (tetraploid).
#' @param seed Integer. Random seed for reproducibility of the random allele selection. Default is \code{123L}.
#'
#' @return A named \code{list} with one element per sample. Each element is a
#'   two-column integer \code{matrix}:
#'   \itemize{
#'     \item Column 1: total read depth (\code{R}).
#'     \item Column 2: read count of one randomly chosen allele (X or Y,
#'       selected with equal probability per site).
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
                           seed      = 123L) {


  # Accept qploidy_standardization objects
  if(is.null(max_geno)) {
    max_geno <- as.numeric(qploidy_standardization_object$info["ploidy.standardization"])
  }

  if (inherits(qploidy_standardization_object, "qploidy_standardization")) {
    qploidy_standardization_object <- qploidy_standardization_object$data
  }

  dat <- qploidy_standardization_object
  required_cols <- c("SampleName", "R", "X", "Y", "ratio", "geno")
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

    # Randomly select X or Y coverage per site (vectorised)
    set.seed(seed)  # For reproducibility
    n       <- nrow(temp)
    choices <- sample(c(1L, 2L), size = n, replace = TRUE)
    allele  <- ifelse(choices == 1L, temp$X, temp$Y)

    xmr <- cbind(R = temp$R, allele = allele)

    # Apply depth and non-zero allele filters
    keep <- !is.na(xmr[, 1L]) & xmr[, 1L] >= min_depth &
      !is.na(xmr[, 2L]) & xmr[, 2L] > 0L
    templist[[i]] <- as.matrix(xmr[keep, , drop = FALSE])
  }

  templist
}


##' Select the best nQuack EM model for a Qploidy standardization object
##'
##' Runs the three nQuack EM model families (Normal, Beta, Beta-Binomial) on
##' the reference samples (those with a non-\code{NA} genotype call in the
##' standardization object) and identifies which distribution/type combination
##' best recovers the known ploidy. Intended to be used upstream of
##' \code{\link[nQuack]{quackNormal}} and friends to choose the most
##' appropriate model before running genome-wide ploidy estimation.
##'
##' @param qploidy_standardization An object of class
##'   \code{qploidy_standardization} as returned by \code{standardize()}.
##'   Must contain a \code{$data} slot with a \code{geno} column and an
##'   \code{$info} slot with a \code{ploidy.standardization} entry.
##' @param n_cores Integer. Number of parallel cores passed to the nQuack EM
##'   functions. Default \code{1}.
##'
##' @return An object of class \code{nQuack_model_selected}, a named list with:
##'   \describe{
##'     \item{\code{all_tested_models}}{A \code{data.frame} with BIC scores for
##'       every distribution/type/ploidy combination tested on each reference
##'       sample.}
##'     \item{\code{test_summary}}{A \code{data.frame} summarising accuracy
##'       (fraction of reference samples with the correct winning ploidy) per
##'       distribution and type.}
##'     \item{\code{best_model}}{A single-row \code{data.frame} giving the
##'       distribution and type with the highest accuracy.}
##'   }
##'
##' @details
##' Only samples with a non-\code{NA} \code{geno} value in the standardization
##' object are used as reference. Ploidy support is limited to 2–6 (diploid
##' through hexaploid) because nQuack models are parameterised for those levels.
##'
##' @seealso \code{\link{convert2nQuack}},
##'   \href{https://mlgaynor.com/nQuack/}{nQuack package}.
##'
##' @importFrom dplyr mutate group_by summarize
##'
##' @export
select_best_nQuack_model <- function(qploidy_standardization,
                                     n_cores = 1) {

    # --- 0. Input checks & derived values ------------------------------------
    if (!inherits(qploidy_standardization, "qploidy_standardization"))
        stop("`qploidy_standardization` must be a qploidy_standardization object.")

    ploidy <- as.numeric(qploidy_standardization$info["ploidy.standardization"])
    if (is.na(ploidy) || ploidy < 2 || ploidy > 6)
        stop("nQuack implementation supports ploidies 2 to 6 only. ",
             "Got: ", ploidy)

    # --- 1. Reference samples (non-NA geno) ----------------------------------
    samples_genos <- unique(
        qploidy_standardization$data$SampleName[
            !is.na(qploidy_standardization$data$geno)
        ]
    )

    if (length(samples_genos) == 0)
        stop("No reference samples found. Ensure the `genos` samples passed to ",
             "`standardize()` are represented in `qploidy_standardization$data`.")

    # --- 2. Convert to nQuack format -----------------------------------------
    nQuack_list <- convert2nQuack(qploidy_standardization,
                                  seed     = 123L,
                                  max_geno = ploidy)

    # --- 3. Run EM models on each reference sample ---------------------------
    em_resultsNormal   <- vector("list", length(samples_genos))
    em_resultsBeta     <- vector("list", length(samples_genos))
    em_resultsBetaBinom <- vector("list", length(samples_genos))

    for (i in seq_along(samples_genos)) {
        xm <- nQuack_list[[samples_genos[i]]]   # index by name, not position

        em_resultsNormal[[i]]    <- quackNormal(xm         = xm,
                                                samplename = samples_genos[i],
                                                cores      = n_cores,
                                                parallel   = n_cores > 1L)

        em_resultsBeta[[i]]      <- quackBeta(xm         = xm,
                                              samplename = samples_genos[i],
                                              cores      = n_cores,
                                              parallel   = n_cores > 1L)

        em_resultsBetaBinom[[i]] <- quackBetaBinom(xm         = xm,
                                                   samplename = samples_genos[i],
                                                   cores      = n_cores,
                                                   parallel   = n_cores > 1L)
    }

    em_results_df <- do.call(rbind, c(em_resultsNormal,
                                      em_resultsBeta,
                                      em_resultsBetaBinom))

    # --- 4. Summarise accuracy per distribution/type -------------------------
    ploidy_names <- c("diploid", "triploid", "tetraploid", "pentaploid", "hexaploid")
    ploidy_name  <- ploidy_names[ploidy - 1L]

    summary_list <- vector("list", length(samples_genos))
    for (i in seq_along(samples_genos)) {
        one_sample       <- em_results_df[em_results_df$sample == samples_genos[i], ]
        summary_list[[i]] <- quackit(model_out = one_sample)
    }
    summary_all <- do.call(rbind, summary_list)

    alloutputcombo <- summary_all %>%
        mutate(accuracy = ifelse(winnerBIC == ploidy_name, 1, 0))

    sumcheck <- alloutputcombo %>%
        group_by(Distribution, Type) %>%
        summarize(total = n(), correct = sum(accuracy), .groups = "drop")

    best_model <- as.data.frame(sumcheck[order(sumcheck$correct, decreasing = TRUE)[1L], ])
    best_model$Uniform <- grepl("-uniform", best_model$Distribution)
    best_model$Distribution <- gsub("-uniform","", best_model$Distribution)

    structure(
        list(all_tested_models = em_results_df,
             test_summary      = sumcheck,
             best_model        = best_model),
        class = "nQuack_model_selected"
    )
}
