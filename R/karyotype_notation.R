utils::globalVariables(c(
    "seg_id", "Ploidy", "Chr_CN", "ISCN_segment"
))

##' Build ISCN-style karyotype notation from HMM copy-number calls
##'
##' Summarises per-marker copy-number calls into a compact karyotype table.
##' For each sample and chromosome, the function derives the sample-wide ploidy
##' (mode of \code{CN_call} across all markers), the chromosome-level copy
##' number (mode of \code{CN_call} within the chromosome), and identifies
##' contiguous segments where the copy number deviates, reporting them in an
##' ISCN-inspired genomic notation (\code{Chr:g.start_end×CN}).
##'
##' @param df A \code{data.frame} with at least the columns \code{SampleName},
##'   \code{Chr}, \code{Position} (base-pair coordinates), and \code{CN_call}
##'   (integer copy-number per marker). Can also be an \code{hmm_CN} object,
##'   in which case the \code{by_marker} slot is used automatically.
##' @param sample_name Character vector of sample names to include, or
##'   \code{NULL} (default) to return notation for all samples in \code{df}.
##'
##' @return A \code{data.frame} with one row per sample–chromosome combination
##'   and the following columns:
##'   \describe{
##'     \item{\code{ChrID}}{Chromosome identifier.}
##'     \item{\code{Ploidy}}{Sample-wide ploidy (modal \code{CN_call}).}
##'     \item{\code{ChrCN}}{Chromosome-level copy number (modal
##'       \code{CN_call} within the chromosome).}
##'     \item{\code{Segment}}{Semicolon-separated ISCN-style notation for
##'       segments that deviate from the chromosome-level copy number, e.g.
##'       \code{"Chr01:g.1000000_5000000×2"}. \code{NA} when no deviating
##'       segments are found.}
##'   }
##'
##' @importFrom dplyr arrange group_by mutate summarise left_join rename case_when
##' @importFrom stats na.omit
##'
##' @examples
##' \dontrun{
##' # All samples
##' karyotype_notation(hmm_CN)
##'
##' # Single sample
##' karyotype_notation(hmm_CN, sample_name = "S1")
##' }
##'
##' @export
karyotype_notation <- function(df, sample_name = NULL) {

    if (inherits(df, "hmm_CN")) {
        df <- df$by_marker
        if (!is.null(sample_name))
            df <- df[df$SampleName %in% sample_name, ]
    } else if (is.list(df) && !is.data.frame(df) &&
               all(vapply(df, inherits, logical(1L), "hmm_CN"))) {
        # Filter the list by name before combining (avoids loading unused samples)
        if (!is.null(sample_name) && !is.null(names(df)))
            df <- df[names(df) %in% sample_name]
        df <- dplyr::bind_rows(lapply(df, `[[`, "by_marker"))
        if (!is.null(sample_name))
            df <- df[df$SampleName %in% sample_name, ]
    }

    if (!is.null(sample_name) && is.data.frame(df))
        df <- df[df$SampleName %in% sample_name, ]

  # --- Sample-level ploidy (mode of CN_call across all markers) ---
  sample_ploidy <- df %>%
    group_by(SampleName) %>%
    summarise(Ploidy = mode(na.omit(CN_call)), .groups = "drop")

  # --- Chromosome-level ploidy (mode of CN_call per chromosome) ---
  chr_ploidy <- df %>%
    group_by(SampleName, Chr) %>%
    summarise(Chr_CN = mode(na.omit(CN_call)), .groups = "drop")

  # --- Segment-level: identify runs of same CN_call per chr/sample ---
  segments <- df %>%
    arrange(SampleName, Chr, Position) %>%
    group_by(SampleName, Chr) %>%
    mutate(
      CN_call_filled = dplyr::coalesce(CN_call, dplyr::lag(CN_call, default = dplyr::first(CN_call))),
      seg_id = cumsum(CN_call_filled != dplyr::lag(CN_call_filled, default = dplyr::first(CN_call_filled))) + 1L,
      CN_call_filled = NULL
    ) %>%
    group_by(SampleName, Chr, seg_id) %>%
    summarise(
      start    = min(Position),
      end      = max(Position),
      seg_CN   = mode(na.omit(CN_call)),
      .groups  = "drop"
    )

  # --- Join chromosome-level CN into segments for ISCN notation ---
  segments <- segments %>%
    left_join(chr_ploidy, by = c("SampleName", "Chr")) %>%
    left_join(sample_ploidy, by = "SampleName")

  # --- Build ISCN-style notation ---
  # Only flag segments that deviate from chromosome-level ploidy
  segments <- segments %>%
    mutate(
      ISCN_segment = case_when(
        seg_CN == Ploidy ~ NA_character_,
        seg_CN == Chr_CN ~ paste0(Chr, ":g.", start, "_", end, "\u00d7", seg_CN),
        TRUE             ~ paste0(Chr, ":g.", start, "_", end, "\u00d7", seg_CN)
      )
    )

  # --- Collapse segments per Sample + Chr into one string ---
  result <- segments %>%
    group_by(SampleName, Ploidy, Chr, Chr_CN) %>%
    summarise(
      Segment = {
        segs <- na.omit(ISCN_segment)
        if (length(segs) == 0) NA_character_ else paste(segs, collapse = "; ")
      },
      .groups = "drop"
    ) %>%
    rename(ChrID = Chr, ChrCN = Chr_CN)

  return(result)
}
