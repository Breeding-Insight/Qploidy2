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
##'   Detecting the absence of Y/W (CN = 0) from CN calls requires \code{0} in
##'   the \code{cn_grid} used by the HMM; the depth-based score (see
##'   \code{sex_het_depth}) does not.
##' @param sex_het_depth \code{data.frame} or \code{NULL}. Presence of Y (XY) or
##'   W (ZW) scored from depth, with columns \code{SampleName}, \code{present}
##'   (logical), and optionally \code{score} (e.g. \code{"5/9"}). When
##'   \code{NULL} (default) and \code{df} is an \code{hmm_CN} object, it is
##'   taken from \code{params$sex_het_depth}, computed by
##'   \code{\link{hmm_estimate_CN}} from the depth of Y/W-specific markers
##'   (\code{sex_min_depth}, \code{sex_min_markers}). When available for a
##'   sample, it decides whether Y/W is present (and so the sex) instead of the
##'   Y/W CN calls, which are usually missing because Y/W-specific markers are
##'   removed by single-sex standardization. The number of Y/W copies is then
##'   taken from CN: called Y/W CN, else PAR - X, else Ploidy - X, else 1;
##'   unknown when neither X nor Y/W has a CN call.
##'   Samples without a score fall back to the CN-based inference. Only used
##'   when \code{sex_system} is \code{"XY"}/\code{"ZW"} and \code{sex_chr}
##'   defines \code{Y}/\code{W}; otherwise it is ignored (with a warning if
##'   supplied explicitly).
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
##'     \item{\code{Sex}}{\code{"male"}, \code{"female"}, \code{"?"} when
##'       \code{SexChr} is \code{"?"} (conflicting sex-chromosome evidence, to
##'       be investigated), or \code{NA} when it cannot be inferred.}
##'     \item{\code{SexChr}}{Sex-chromosome complement, e.g. \code{"XY"},
##'       \code{"XXXY"}, \code{"XXYY"}, \code{"ZW"}, \code{"X0"}. Each
##'       \code{"0"} marks a copy missing relative to the ploidy (e.g.
##'       \code{"XXX0"}: tetraploid with PAR = X = 3 and no Y). \code{"?"} when
##'       Y/W is absent by depth but PAR > X/Z, so the extra PAR copies are not
##'       explained, or when Y/W is present by depth but neither X/Z nor Y/W
##'       has a CN call, so PAR cannot be split between them (e.g. PAR 4 could
##'       be XXXY, XXYY, or XYYY), or whenever the depth evidence contradicts
##'       PAR vs X/Z (\code{SexChrConsistent} is \code{FALSE}): Y/W present
##'       with PAR == X/Z (likely contamination), Y/W absent with PAR > X/Z, or
##'       PAR < X/Z. These should be investigated; deviating X/Z copies are
##'       still flagged. When X/Z is segmental (several CNs) and its modal CN
##'       contradicts the depth evidence while another X/Z segment agrees
##'       with it, that segment sets the complement and the rest is reported
##'       in \code{Segment} as a segmental gain or loss (e.g. Y absent, PAR 2,
##'       X 2 on one arm and 1 on the larger other: \code{"XX"}, a female with
##'       a partial X deletion). \code{NA} if unknown.}
##'     \item{\code{SexChrConsistent}}{For \code{"XY"}/\code{"ZW"}: whether the
##'       sex-chromosome evidence agrees. With a Y/W depth score
##'       (\code{sex_het_depth}), Y/W present must go with PAR > X/Z
##'       (\code{TRUE}), and Y/W absent with PAR == X/Z. Y/W present with
##'       PAR == X (likely contamination), PAR > X without Y/W, or PAR < X
##'       (impossible) give \code{FALSE}, and \code{SexChr} is then
##'       \code{"?"}. Without a score, whether X + Y (or
##'       Z + W) equals the PAR copy number when all three groups are called.
##'       \code{NA} when the needed groups are missing.}
##'     \item{\code{SexHetCN}}{Estimated number of Y (XY) or W (ZW) copies. Y/W
##'       has no CN calls of its own, so with a depth score this is 0 when
##'       Y/W is absent and PAR - X/Z when present (e.g. ploidy 4, PAR 4, X 2
##'       gives 2), falling back to ploidy - X/Z and then 1. \code{NA} when Y/W
##'       is present but X/Z has no CN call.}
##'     \item{\code{SexHetScore}}{Y/W depth score: Y/W-specific markers with
##'       depth above \code{sex_min_depth} over all Y/W-specific markers (e.g.
##'       \code{"5/9"}). \code{NA} when not available.}
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
                               sex_chr = NULL,
                               sex_het_depth = NULL) {

    sex <- .resolve_sex_args(df, sex_system, sex_chr)
    sex_system <- sex$sex_system
    sex_chr    <- sex$sex_chr
    sex_roles  <- .sex_chr_roles(sex_system, sex_chr)
    # Depth-based Y/W presence only applies when the Y (XY) or W (ZW) group is
    # defined in sex_chr; otherwise sex is inferred from CN calls alone
    if (!.het_depth_applies(sex_system, sex_chr)) {
      if (!is.null(sex_het_depth))
        warning("`sex_het_depth` is ignored: it requires sex_system \"XY\" or \"ZW\" ",
                "with `sex_chr$Y` (or `sex_chr$W`) defined.", call. = FALSE)
      sex_het_depth <- NULL
    } else if (is.null(sex_het_depth)) {
      sex_het_depth <- .het_depth_table(df)
    }
    sex_het_depth <- .check_het_depth_table(sex_het_depth)

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

  df$ChrType <- rep("autosome", nrow(df))
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
      # All X/Z CNs, to resolve segmental X/Z against the Y/W depth evidence
      hom_cns = list(sort(unique(CN_call[ChrType %in% hom & !is.na(CN_call)]))),
      .groups = "drop"
    ) %>%
    left_join(sample_ploidy, by = "SampleName") %>%
    left_join(sex_het_depth, by = "SampleName")
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
      left_join(sex_tbl[, c("SampleName", "Sex", "SexChr", "SexChrConsistent",
                            "SexHetCN", "SexHetScore")],
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

# Depth-based Y/W presence applies only to XY/ZW with the Y/W group in sex_chr
.het_depth_applies <- function(sex_system, sex_chr) {
  het_grp <- switch(sex_system, XY = "Y", ZW = "W", NULL)
  !is.null(het_grp) && het_grp %in% names(sex_chr)
}

# Per-sample Y/W depth scores stored by hmm_estimate_CN() in params$sex_het_depth.
# Returns data.frame(SampleName, present, score) or NULL.
.het_depth_table <- function(x) {
  plist <- NULL
  if (inherits(x, "hmm_CN")) {
    if (!is.null(x$params_samples)) {
      plist <- x$params_samples
    } else if (!is.null(x$params)) {
      s <- unique(c(x$by_window$Sample, x$by_marker$SampleName))
      if (length(s) == 1) plist <- stats::setNames(list(x$params), s)
    }
  } else if (is.list(x) && !is.data.frame(x) && length(x) > 0 &&
             all(vapply(x, inherits, logical(1L), "hmm_CN"))) {
    plist <- lapply(x, `[[`, "params")
    names(plist) <- vapply(x, function(o)
      as.character(unique(c(o$by_window$Sample, o$by_marker$SampleName))[1]), character(1))
  }
  if (is.null(plist)) return(NULL)
  rows <- lapply(names(plist), function(s) {
    h <- plist[[s]]$sex_het_depth
    if (is.null(h) || is.null(h$present)) return(NULL)
    data.frame(SampleName = s, present = as.logical(h$present),
               score = as.character(h$score), stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) NULL else do.call(rbind, rows)
}

# Validate/normalize a sex_het_depth table to data.frame(SampleName, het_present, SexHetScore)
.check_het_depth_table <- function(x) {
  empty <- data.frame(SampleName = character(0), het_present = logical(0),
                      SexHetScore = character(0), stringsAsFactors = FALSE)
  if (is.null(x)) return(empty)
  if (!is.data.frame(x) || !all(c("SampleName", "present") %in% names(x)))
    stop("`sex_het_depth` must be a data.frame with columns `SampleName` and `present`.")
  out <- data.frame(SampleName = as.character(x$SampleName),
                    het_present = as.logical(x$present),
                    SexHetScore = if ("score" %in% names(x)) as.character(x$score) else NA_character_,
                    stringsAsFactors = FALSE)
  out[!duplicated(out$SampleName), , drop = FALSE]
}

# Base X/Z copy number given Y/W presence by depth and all X/Z CNs of the sample.
# Keeps the modal X/Z unless it contradicts the depth rule (present: PAR > X;
# absent: PAR == X) and another X/Z segment satisfies it.
.base_hom_cn <- function(cn_hom, cn_par, present, cns) {
  if (is.na(present) || is.na(cn_hom) || is.na(cn_par) || length(cns) < 2)
    return(cn_hom)
  if (present) {
    if (cn_hom < cn_par) return(cn_hom)
    below <- cns[cns < cn_par]
    if (length(below) > 0) max(below) else cn_hom
  } else {
    if (cn_hom == cn_par || !cn_par %in% cns) cn_hom else cn_par
  }
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
    dep <- if ("het_present" %in% names(tbl)) as.logical(tbl$het_present) else rep(NA, n)
    # Segmental X/Z: when the modal X/Z contradicts the depth evidence but another
    # X/Z segment agrees with it, that segment gives the base X/Z copy number and
    # the rest is a segmental gain/loss (flagged in Segment). E.g. Y absent,
    # PAR 2, X 2 on one arm and 1 on the (larger) other: XX female with a deletion
    if ("hom_cns" %in% names(tbl))
      tbl$cn_hom <- as.integer(unname(mapply(.base_hom_cn, tbl$cn_hom, tbl$cn_par,
                                             dep, tbl$hom_cns)))

    # CN-based counts (used when no depth score is available)
    n_hom <- dplyr::coalesce(tbl$cn_hom, tbl$cn_par - tbl$cn_het)
    n_het <- dplyr::coalesce(tbl$cn_het, tbl$cn_par - tbl$cn_hom)

    # Depth-based presence of Y/W overrides the CN-based one; copy counts still from CN
    if (any(!is.na(dep))) {
      pos <- function(x) ifelse(!is.na(x) & x >= 1L, as.integer(x), NA_integer_)
      het_if_present <- dplyr::coalesce(pos(tbl$cn_het), pos(tbl$cn_par - tbl$cn_hom),
                                        pos(P - tbl$cn_hom), 1L)
      # Neither X/Z nor Y/W called: PAR only gives X + Y, so the split is unknown
      het_if_present[is.na(tbl$cn_hom) & is.na(tbl$cn_het)] <- NA_integer_
      n_het_dep <- ifelse(dep, het_if_present, 0L)
      n_hom_dep <- ifelse(dep,
                          dplyr::coalesce(tbl$cn_hom, tbl$cn_par - n_het_dep, P - n_het_dep),
                          dplyr::coalesce(tbl$cn_hom, tbl$cn_par, P))
      n_het <- ifelse(is.na(dep), n_het, n_het_dep)
      n_hom <- ifelse(is.na(dep), n_hom, n_hom_dep)
    }

    valid <- !is.na(n_hom) & !is.na(n_het) & n_hom >= 0 & n_het >= 0 & (n_hom + n_het) > 0
    has_het <- ifelse(is.na(dep), ifelse(valid, n_het >= 1, NA), dep)
    tbl$Sex <- ifelse(is.na(has_het), NA_character_,
                      ifelse(has_het == (sex_system == "XY"), "male", "female"))
    # Estimated Y/W copies (0 when a depth score says Y/W is absent)
    tbl$SexHetCN <- ifelse(valid, as.integer(n_het), NA_integer_)
    # Sex-chromosome complement; each "0" marks a copy missing relative to the
    # ploidy (e.g. XXX0 in a tetraploid with PAR = X = 3 and no Y)
    n_missing <- ifelse(is.na(P), 0L, pmax(P - n_hom - n_het, 0L))
    tbl$SexChr <- ifelse(valid, paste0(strrep(hom, pmax(n_hom, 0L)),
                                       strrep(het, pmax(n_het, 0L)),
                                       strrep("0", n_missing)), NA_character_)
    # Consistency of the sex-chromosome evidence:
    #  - with a depth score, Y/W presence must agree with the PAR vs X/Z copy
    #    number: present <=> PAR > X (PAR - X >= 1 copies of Y/W), absent <=>
    #    PAR == X. Anything else is inconsistent: Y/W present with PAR == X
    #    (likely contamination), PAR > X without Y/W, or PAR < X (impossible);
    #  - otherwise (no score) X + Y must equal PAR when all three are called.
    par_minus_hom <- tbl$cn_par - tbl$cn_hom
    all3 <- !is.na(tbl$cn_hom) & !is.na(tbl$cn_het) & !is.na(tbl$cn_par)
    cn_rule <- ifelse(all3, tbl$cn_hom + tbl$cn_het == tbl$cn_par, NA)
    depth_rule <- ifelse(is.na(dep) | is.na(par_minus_hom), NA,
                         ifelse(dep, par_minus_hom >= 1, par_minus_hom == 0))
    tbl$SexChrConsistent <- dplyr::coalesce(as.logical(depth_rule), as.logical(cn_rule))
    # Complement cannot be resolved ("?", to be investigated) when:
    #  - the depth evidence contradicts PAR vs X/Z (depth_rule FALSE above);
    #    Y/W copies and expected CNs are still estimated (ploidy - X, else 1
    #    when present), so the X/Z excess stays flagged;
    #  - Y/W is present but its number of copies is unknown (no X/Z or Y/W call)
    unresolved <- (!is.na(depth_rule) & !depth_rule) |
                  (!is.na(dep) & dep & is.na(n_het))
    tbl$SexChr[unresolved] <- "?"
    # An unresolved complement leaves the sex unresolved too
    tbl$Sex[unresolved] <- "?"
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
  if (!"SexHetScore" %in% names(tbl)) tbl$SexHetScore <- NA_character_
  if (!"SexHetCN" %in% names(tbl)) tbl$SexHetCN <- NA_integer_
  tbl
}
