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
#' @seealso \code{\link[https://mlgaynor.com/nQuack/articles/Qploidy2nQuack.html]{nQuack tutorial}} for ploidy estimation using the nQuack package.
#'
#'
#' @export
convert2nQuack <- function(qploidy_standardization_object,
                           min_depth = 10L,
                           max_geno  = NULL,
                           seed      = 123L) {
  # Accept qploidy_standardization objects
  if (inherits(qploidy_standardization_object, "qploidy_standardization")) {
    qploidy_standardization_object <- qploidy_standardization_object$data
  }

  if(is.null(max_geno)) {
    max_geno <- as.numeric(qploidy_standardization_object$info$ploidy.standardization)
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
