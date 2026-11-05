utils::globalVariables(c(
    "seg_id", "Ploidy", "Chr_CN", "ISCN_segment", "CN_call_filled",
    "ChrType", "ExpectedCN", "exp_hom", "exp_het"
))

##' Build ISCN-style karyotype notation from HMM copy-number calls
##'
##' Summarises per-marker copy-number calls into a compact karyotype table.
##' For each sample and chromosome, the function derives the sample-wide ploidy
##' (mode of \code{CN_call} across all markers), the chromosome-level copy
##' number (mode of \code{CN_call} within the chromosome), and identifies
##' contiguous segments where the copy number deviates, reporting them in an
##' ISCN-inspired genomic notation (\code{Chr:g.start_endxCN}).
##'
##' @param df A \code{data.frame} with at least the columns \code{SampleName},
##'   \code{Chr}, \code{Position} (base-pair coordinates), and \code{CN_call}
##'   (integer copy-number per marker). Can also be an \code{hmm_CN} object,
##'   in which case the \code{by_marker} slot is used automatically.
##' @param sample_name Character vector of sample names to include, or
##'   \code{NULL} (default) to return notation for all samples in \code{df}.
##' @param sex_system Character or \code{NULL}. Sex-determination system of the species.
##'   When \code{NULL} (default), the value stored in the \code{params} of an
##'   \code{hmm_CN} input (set through \code{sex_system} in
##'   \code{\link{hmm_estimate_CN}}) is used, falling back to \code{"none"}.
##'   Supplying a value overrides the stored one. Options:
##'   \describe{
##'     \item{\code{"none"}}{All chromosomes are treated as autosomes.}
##'     \item{\code{"XY"}}{Females XX, males XY. A sample carrying at least one
##'       Y is reported as male.}
##'     \item{\code{"ZW"}}{Males ZZ, females ZW (e.g. birds, Lepidoptera). A
##'       sample carrying at least one W is reported as female.}
##'     \item{\code{"X0"}}{Females XX, males X0. Sex is inferred from the X copy
##'       number: equal to the ploidy is female, half the ploidy is male.}
##'     \item{\code{"haplodiploid"}}{Males haploid, females diploid (e.g. bees,
##'       ants, wasps). Sex is inferred from the sample ploidy (1 = male,
##'       2 = female); no \code{sex_chr} is needed. Diploid males cannot be
##'       distinguished from females.}
##'   }
##' @param sex_chr Named list mapping sex-chromosome groups to chromosome IDs
##'   in \code{df$Chr}. When \code{NULL}, the value stored in the \code{hmm_CN}
##'   \code{params} is used if its \code{sex_system} matches. Allowed names: \code{X}, \code{Y}, \code{PAR} for
##'   \code{"XY"}; \code{Z}, \code{W}, \code{PAR} for \code{"ZW"}; \code{X} for
##'   \code{"X0"}. \code{X}/\code{Z} and \code{Y}/\code{W} hold markers
##'   specific to that chromosome; \code{PAR} holds pseudoautosomal markers
##'   shared by both. Each element may contain several IDs. For
##'   \code{"XY"}/\code{"ZW"}, the copy number of a missing group is derived
##'   from the PAR (e.g. Y = PAR - X); when all three groups are given, the
##'   relation X + Y = PAR is checked and reported in \code{SexChrConsistent}.
##'   Detecting the absence of Y/W (CN = 0) requires \code{0} in the
##'   \code{cn_grid} used by the HMM.
##'
##' @return A \code{data.frame} with one row per sample-chromosome combination
##'   and the following columns:
##'   \describe{
##'     \item{\code{ChrID}}{Chromosome identifier.}
##'     \item{\code{Ploidy}}{Sample-wide ploidy (modal \code{CN_call}; sex
##'       chromosomes are excluded when \code{sex_system} is set).}
##'     \item{\code{ChrCN}}{Chromosome-level copy number (modal
##'       \code{CN_call} within the chromosome).}
##'     \item{\code{Segment}}{Semicolon-separated ISCN-style notation for
##'       segments that deviate from the expected copy number, e.g.
##'       \code{"Chr01:g.1000000_5000000x2"}. \code{NA} when no deviating
##'       segments are found.}
##'   }
##'   When \code{sex_system != "none"}, the following columns are added:
##'   \describe{
##'     \item{\code{Sex}}{\code{"male"}, \code{"female"}, or \code{NA} when it
##'       cannot be inferred.}
##'     \item{\code{SexChr}}{Sex-chromosome complement, e.g. \code{"XY"},
##'       \code{"XXXY"}, \code{"ZW"}, \code{"X0"}. \code{NA} if unknown.}
##'     \item{\code{SexChrConsistent}}{For \code{"XY"}/\code{"ZW"} with all
##'       three groups defined: whether X + Y (or Z + W) equals the PAR copy
##'       number. \code{NA} otherwise.}
##'     \item{\code{ChrType}}{\code{"autosome"}, \code{"X"}, \code{"Y"},
##'       \code{"Z"}, \code{"W"}, or \code{"PAR"}.}
##'     \item{\code{ExpectedCN}}{Copy number expected for the chromosome given
##'       the ploidy and the inferred sex-chromosome complement. Segments are
##'       flagged against this value, so e.g. a single X in an XY male is not
##'       reported as a deviation, while an extra X in an XXY diploid is.}
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
##'
##' # XY system with X-specific, Y-specific, and pseudoautosomal contigs
##' karyotype_notation(hmm_CN, sex_system = "XY",
##'                    sex_chr = list(X = "ChrX", Y = "ChrY", PAR = "ChrX_PAR"))
##' }
##'
##' @export
karyotype_notation <- function(df, sample_name = NULL,
                               sex_system = NULL,
                               sex_chr = NULL) {

    sex <- .resolve_sex_args(df, sex_system, sex_chr)
    sex_system <- sex$sex_system
    sex_chr    <- sex$sex_chr
    sex_roles  <- .sex_chr_roles(sex_system, sex_chr)

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

  hom <- switch(sex_system, XY = "X", ZW = "Z", X0 = "X", NA_character_)
  het <- switch(sex_system, XY = "Y", ZW = "W", NA_character_)

  df$ChrType <- "autosome"
  if (!is.null(sex_roles)) {
    idx <- match(df$Chr, sex_roles$Chr)
    df$ChrType[!is.na(idx)] <- sex_roles$ChrType[idx[!is.na(idx)]]
    if (all(is.na(idx)))
      warning("None of the chromosome IDs in `sex_chr` were found in `df$Chr`.")
  }

  # --- Sample-level ploidy (mode of CN_call across autosomal/PAR markers) ---
  sample_ploidy <- df %>%
    group_by(SampleName) %>%
    summarise(Ploidy = {
      x <- CN_call[!ChrType %in% c(hom, het)]
      if (all(is.na(x))) x <- CN_call
      mode(na.omit(x))
    }, .groups = "drop")

  # --- Sex inference and expected copy numbers of sex chromosomes ---
  sex_tbl <- df %>%
    group_by(SampleName) %>%
    summarise(
      cn_hom = .mode_or_na(CN_call[ChrType %in% hom]),
      cn_het = .mode_or_na(CN_call[ChrType %in% het]),
      cn_par = .mode_or_na(CN_call[ChrType == "PAR"]),
      .groups = "drop"
    ) %>%
    left_join(sample_ploidy, by = "SampleName")
  sex_tbl <- .infer_sex(sex_tbl, sex_system, hom, het)

  # --- Chromosome-level ploidy (mode of CN_call per chromosome) ---
  chr_ploidy <- df %>%
    group_by(SampleName, Chr, ChrType) %>%
    summarise(Chr_CN = mode(na.omit(CN_call)), .groups = "drop") %>%
    left_join(sample_ploidy, by = "SampleName") %>%
    left_join(sex_tbl[, c("SampleName", "exp_hom", "exp_het")], by = "SampleName") %>%
    mutate(ExpectedCN = case_when(
      ChrType %in% hom ~ exp_hom,
      ChrType %in% het ~ exp_het,
      TRUE             ~ as.integer(Ploidy)
    )) %>%
    select(-"Ploidy", -"exp_hom", -"exp_het")

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
  # Only flag segments that deviate from the expected copy number; for sex
  # chromosomes of unknown complement fall back to the chromosome-level CN
  segments <- segments %>%
    mutate(
      ISCN_segment = case_when(
        seg_CN == dplyr::coalesce(ExpectedCN, as.integer(Chr_CN)) ~ NA_character_,
        TRUE ~ paste0(Chr, ":g.", start, "_", end, "\u00d7", seg_CN)
      )
    )

  # --- Collapse segments per Sample + Chr into one string ---
  result <- segments %>%
    group_by(SampleName, Ploidy, Chr, Chr_CN, ChrType, ExpectedCN) %>%
    summarise(
      Segment = {
        segs <- na.omit(ISCN_segment)
        if (length(segs) == 0) NA_character_ else paste(segs, collapse = "; ")
      },
      .groups = "drop"
    ) %>%
    rename(ChrID = Chr, ChrCN = Chr_CN)

  if (sex_system == "none") {
    result <- result[, c("SampleName", "Ploidy", "ChrID", "ChrCN", "Segment")]
  } else {
    result <- result %>%
      left_join(sex_tbl[, c("SampleName", "Sex", "SexChr", "SexChrConsistent")],
                by = "SampleName") %>%
      dplyr::relocate("Segment", .after = dplyr::last_col())
  }

  return(result)
}

# Resolve sex_system/sex_chr: explicit arguments win, then hmm_CN params, then "none"
.resolve_sex_args <- function(x, sex_system = NULL, sex_chr = NULL) {
  choices <- c("none", "XY", "ZW", "X0", "haplodiploid")

  plist <- NULL
  if (inherits(x, "hmm_CN")) {
    plist <- if (!is.null(x$params_samples)) x$params_samples else list(x$params)
  } else if (is.list(x) && !is.data.frame(x) &&
             length(x) > 0 && all(vapply(x, inherits, logical(1L), "hmm_CN"))) {
    plist <- lapply(x, `[[`, "params")
  }

  stored_sys <- NULL
  stored_chr <- NULL
  if (!is.null(plist)) {
    sys <- unique(unlist(lapply(plist, `[[`, "sex_system")))
    if (length(sys) > 1)
      warning("Samples in the hmm_CN object have different `sex_system` values (",
              paste(sys, collapse = ", "), "); using \"", sys[1], "\".", call. = FALSE)
    if (length(sys) > 0) {
      stored_sys <- sys[1]
      hit <- Filter(function(p) identical(p$sex_system, stored_sys), plist)
      stored_chr <- hit[[1]]$sex_chr
    }
  }

  if (is.null(sex_system)) sex_system <- if (is.null(stored_sys)) "none" else stored_sys
  sex_system <- match.arg(sex_system, choices)
  if (is.null(sex_chr) && identical(sex_system, stored_sys)) sex_chr <- stored_chr

  list(sex_system = sex_system, sex_chr = sex_chr)
}

# Validate `sex_chr` against `sex_system`; returns data.frame(Chr, ChrType) or NULL
.sex_chr_roles <- function(sex_system, sex_chr) {
  if (sex_system %in% c("none", "haplodiploid")) {
    if (!is.null(sex_chr))
      warning("`sex_chr` is ignored when sex_system = \"", sex_system, "\".")
    return(NULL)
  }

  allowed <- switch(sex_system,
                    XY = c("X", "Y", "PAR"),
                    ZW = c("Z", "W", "PAR"),
                    X0 = "X")

  if (!is.list(sex_chr) || length(sex_chr) == 0 || is.null(names(sex_chr)) ||
      any(names(sex_chr) == ""))
    stop("`sex_chr` must be a named list with elements among: ",
         paste(allowed, collapse = ", "), ".")
  bad <- setdiff(names(sex_chr), allowed)
  if (length(bad) > 0)
    stop("Invalid `sex_chr` name(s) for sex_system = \"", sex_system, "\": ",
         paste(bad, collapse = ", "), ". Allowed: ", paste(allowed, collapse = ", "), ".")
  if (!all(vapply(sex_chr, function(x) is.character(x) && length(x) > 0, logical(1))))
    stop("Each element of `sex_chr` must be a non-empty character vector of chromosome IDs.")

  ids <- unlist(sex_chr, use.names = FALSE)
  if (anyDuplicated(ids))
    stop("Chromosome IDs cannot belong to more than one `sex_chr` group: ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "), ".")

  g <- names(sex_chr)
  if (sex_system == "X0" && !"X" %in% g)
    stop("sex_system = \"X0\" requires `sex_chr$X`.")
  if (sex_system %in% c("XY", "ZW")) {
    hom <- allowed[1]; het <- allowed[2]
    # The heterogametic count must be observable directly or via PAR - hom
    if (!(het %in% g || all(c(hom, "PAR") %in% g)))
      stop("sex_system = \"", sex_system, "\" requires `sex_chr$", het,
           "`, or both `sex_chr$", hom, "` and `sex_chr$PAR`.")
  }

  data.frame(Chr = ids,
             ChrType = rep(g, lengths(sex_chr)),
             stringsAsFactors = FALSE)
}

.mode_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_integer_ else as.integer(mode(x))
}

# Adds Sex, SexChr, SexChrConsistent, exp_hom, exp_het to the per-sample table
.infer_sex <- function(tbl, sex_system, hom, het) {
  P <- as.integer(tbl$Ploidy)
  n <- nrow(tbl)
  tbl$Sex <- NA_character_
  tbl$SexChr <- NA_character_
  tbl$SexChrConsistent <- NA
  tbl$exp_hom <- NA_integer_
  tbl$exp_het <- NA_integer_

  if (sex_system %in% c("XY", "ZW")) {
    n_hom <- dplyr::coalesce(tbl$cn_hom, tbl$cn_par - tbl$cn_het)
    n_het <- dplyr::coalesce(tbl$cn_het, tbl$cn_par - tbl$cn_hom)
    valid <- !is.na(n_hom) & !is.na(n_het) & n_hom >= 0 & n_het >= 0 & (n_hom + n_het) > 0
    has_het <- n_het >= 1
    tbl$Sex <- ifelse(!valid, NA_character_,
                      ifelse(has_het == (sex_system == "XY"), "male", "female"))
    tbl$SexChr <- ifelse(valid, paste0(strrep(hom, pmax(n_hom, 0L)),
                                       strrep(het, pmax(n_het, 0L))), NA_character_)
    all3 <- !is.na(tbl$cn_hom) & !is.na(tbl$cn_het) & !is.na(tbl$cn_par)
    tbl$SexChrConsistent <- ifelse(all3, tbl$cn_hom + tbl$cn_het == tbl$cn_par, NA)
    # A euploid sex complement fills the ploidy: hom + het == Ploidy
    tbl$exp_hom <- ifelse(valid, pmax(P - n_het, 0L), NA_integer_)
    tbl$exp_het <- ifelse(valid, pmax(P - n_hom, 0L), NA_integer_)
  } else if (sex_system == "X0") {
    n_hom <- tbl$cn_hom
    female <- !is.na(n_hom) & n_hom == P
    male   <- !is.na(n_hom) & 2L * n_hom == P
    tbl$Sex <- ifelse(female, "female", ifelse(male, "male", NA_character_))
    tbl$SexChr <- ifelse(is.na(n_hom), NA_character_,
                         paste0(strrep("X", pmax(n_hom, 0L)), ifelse(male, "0", "")))
    tbl$exp_hom <- ifelse(female, P, ifelse(male, P %/% 2L, NA_integer_))
  } else if (sex_system == "haplodiploid") {
    tbl$Sex <- ifelse(P %in% 1L, "male", ifelse(P %in% 2L, "female", NA_character_))
  }

  tbl$exp_hom <- as.integer(tbl$exp_hom)
  tbl$exp_het <- as.integer(tbl$exp_het)
  tbl
}
