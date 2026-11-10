library(testthat)
library(ggplot2)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

# Minimal marker data frame: 2 samples, 3 chromosomes, uniform CN
make_marker_df <- function() {
  data.frame(
    SampleName = rep(c("S1", "S2"), each = 30),
    Chr        = rep(rep(c("Chr01", "Chr02", "Chr03"), each = 10), 2),
    Position   = rep(seq(1e6, 10e6, length.out = 10), 6),
    CN_call    = rep(4L, 60),
    dosage     = sample(0:4, 60, replace = TRUE),
    w_baf      = runif(60, 0, 1),
    post_max_CN     = runif(60, 0.5, 1),
    post_max_dosage = runif(60, 0.5, 1),
    stringsAsFactors = FALSE
  )
}

# Data with a CN change on Chr01 for S1 (CN = 4 then CN = 2)
make_marker_df_with_segment <- function() {
  df <- make_marker_df()
  # First 5 markers of Chr01 for S1 get CN = 2
  idx <- df$SampleName == "S1" & df$Chr == "Chr01"
  df$CN_call[idx][1:5] <- 2L
  df
}

# ---------------------------------------------------------------------------
# karyotype_notation
# ---------------------------------------------------------------------------

test_that("karyotype_notation returns a data.frame with correct columns", {
  df  <- make_marker_df()
  out <- karyotype_notation(df)

  expect_s3_class(out, "data.frame")
  expect_named(out, c("SampleName", "Ploidy", "ChrID", "ChrCN", "Segment"),
               ignore.order = TRUE)
})

test_that("karyotype_notation returns one row per sample-chromosome", {
  df  <- make_marker_df()
  out <- karyotype_notation(df)

  # 2 samples × 3 chromosomes = 6 rows
  expect_equal(nrow(out), 6L)
})

test_that("karyotype_notation infers sample-wide ploidy correctly", {
  df  <- make_marker_df()   # all CN = 4
  out <- karyotype_notation(df)

  expect_true(all(out$Ploidy == 4L))
})

test_that("karyotype_notation sets Segment to NA when CN is uniform", {
  df  <- make_marker_df()   # no CN changes
  out <- karyotype_notation(df)

  expect_true(all(is.na(out$Segment)))
})

test_that("karyotype_notation reports deviating segment in ISCN notation", {
  df  <- make_marker_df_with_segment()
  out <- karyotype_notation(df)

  seg_row <- out[out$SampleName == "S1" & out$ChrID == "Chr01", ]
  expect_false(is.na(seg_row$Segment))
  # Should contain the chr name, genomic coordinates, and the copy number
  expect_equal(seg_row$Segment, "Chr01:g.1e+06_5e+06×2")
})

test_that("karyotype_notation Segment is NA for unaffected chromosomes", {
  df  <- make_marker_df_with_segment()
  out <- karyotype_notation(df)

  # Chr02 and Chr03 for S1 are unaffected
  unaffected <- out[out$SampleName == "S1" & out$ChrID %in% c("Chr02", "Chr03"), ]
  expect_true(all(is.na(unaffected$Segment)))
})

test_that("karyotype_notation handles multiple samples independently", {
  df  <- make_marker_df_with_segment()
  out <- karyotype_notation(df)

  # S2 has no CN changes
  s2_segs <- out[out$SampleName == "S2", ]
  expect_true(all(is.na(s2_segs$Segment)))
})

# ---------------------------------------------------------------------------
# plot_karyotype
# ---------------------------------------------------------------------------

test_that("plot_karyotype returns a ggplot object", {
  df <- make_marker_df()
  p  <- plot_karyotype(df, sample_name = "S1", color_by = "dosage")
  expect_s3_class(p, "gg")
})

test_that("plot_karyotype works for each color_by option", {
  df <- make_marker_df()
  for (cb in c("dosage", "w_baf", "post_max_CN", "post_max_dosage")) {
    p <- plot_karyotype(df, sample_name = "S1", color_by = cb)
    expect_s3_class(p, "gg")
  }
})

test_that("plot_karyotype works with color_by = NULL (black tick marks)", {
  df <- make_marker_df()
  p  <- plot_karyotype(df, sample_name = "S1", color_by = NULL)
  expect_s3_class(p, "gg")
})

test_that("plot_karyotype works with color_by = NA", {
  df <- make_marker_df()
  p  <- plot_karyotype(df, sample_name = "S1", color_by = NA)
  expect_s3_class(p, "gg")
})

test_that("plot_karyotype respects nrow argument", {
  df <- make_marker_df()
  p  <- plot_karyotype(df, sample_name = "S1", color_by = "dosage", nrow = 2)
  expect_s3_class(p, "gg")
})

test_that("plot_karyotype works with notation argument", {
  df       <- make_marker_df_with_segment()
  notation <- karyotype_notation(df)
  p        <- plot_karyotype(df, sample_name = "S1", color_by = "dosage",
                              notation = FALSE)
  expect_s3_class(p, "gg")
})

test_that("plot_karyotype shows segment ranges only with show_segments = TRUE", {
  df  <- make_marker_df_with_segment()
  seg <- karyotype_notation(df, sample_name = "S1")
  chr <- seg$ChrID[!is.na(seg$Segment)][1]
  strip <- function(p) {
    lab <- ggplot2::ggplot_build(p)$layout$facet$params$labeller
    lab(data.frame(Chr = chr))[[1]][[1]]
  }
  off <- strip(plot_karyotype(df, sample_name = "S1"))
  on  <- strip(plot_karyotype(df, sample_name = "S1", show_segments = TRUE))
  expect_false(grepl("g\\.", off))
  expect_match(off, "CN: ")
  expect_match(on, "g\\.")
  expect_error(plot_karyotype(df, "S1", show_segments = NA), "show_segments")
})

# --- Input validation -------------------------------------------------------

test_that("plot_karyotype errors when df is not a data.frame", {
  expect_error(plot_karyotype(list(), "S1"), "`df` must be a data.frame")
})

test_that("plot_karyotype errors on missing required columns", {
  df <- make_marker_df()
  df$CN_call <- NULL
  expect_error(plot_karyotype(df, "S1", color_by = "dosages"))
})

test_that("plot_karyotype errors when sample_name is not found", {
  df <- make_marker_df()
  expect_error(plot_karyotype(df, "UNKNOWN"), "not found")
})

test_that("plot_karyotype errors when sample_name is not a scalar string", {
  df <- make_marker_df()
  expect_error(plot_karyotype(df, c("S1", "S2")), "single character string")
})

test_that("plot_karyotype errors when nrow is invalid", {
  df <- make_marker_df()
  expect_error(plot_karyotype(df, "S1", nrow = 0),   "positive integer")
  expect_error(plot_karyotype(df, "S1", nrow = 1.5), "positive integer")
})

# ---------------------------------------------------------------------------
# Sex chromosomes
# ---------------------------------------------------------------------------

# Marker-level CN calls as the pipeline really produces them: 2 autosomes plus
# X-specific (hom) and PAR contigs. Y/W-specific markers are dropped by
# standardization (a single-sex reference has no clusters for them), so the
# het chromosome has NO rows and its presence is given separately, from raw
# depth, through het_depth(). Only supply Y/W in `cn` to emulate the legacy case
# of a CN-called Y/W; omit X/Z too to emulate a dropped X-specific chromosome.
# cn is a named list per sample: list(X = , PAR = ) [, Y = ]
make_sex_df <- function(cn, hom = "X", het = "Y", ploidy = 2L) {
  do.call(rbind, lapply(names(cn), function(s) {
    v <- cn[[s]]
    calls <- c(Chr01 = ploidy, Chr02 = ploidy)
    if (!is.null(v[[hom]])) calls[paste0("Chr", hom)] <- v[[hom]]
    if (!is.null(v[[het]])) calls[paste0("Chr", het)] <- v[[het]]
    calls["ChrPAR"] <- v$PAR
    data.frame(
      SampleName = s,
      Chr        = rep(names(calls), each = 10),
      Position   = rep(seq(1e6, 10e6, length.out = 10), length(calls)),
      CN_call    = rep(as.integer(calls), each = 10),
      stringsAsFactors = FALSE
    )
  }))
}
xy_chr <- list(X = "ChrX", Y = "ChrY", PAR = "ChrPAR")

# Depth-based Y/W presence as stored by hmm_estimate_CN() (params$sex_het_depth):
# `present` is a named logical (sample -> Y/W present); the score mimics 9
# Y/W-specific markers of which 6 are covered when present, none otherwise.
het_depth <- function(present) {
  data.frame(SampleName = names(present), present = unname(present),
             score = ifelse(present, "6/9", "0/9"), stringsAsFactors = FALSE)
}
het_params <- function(present, score = ifelse(present, "6/9", "0/9"), ...) {
  list(sex_het_depth = list(present = present, score = score), ...)
}

test_that("karyotype_notation XY infers sex and complement from Y depth", {
  # No ChrY rows: Y-specific markers never get a CN call
  df  <- make_sex_df(list(F = list(X = 2L, PAR = 2L),
                          M = list(X = 1L, PAR = 2L)))
  expect_false("ChrY" %in% df$Chr)
  out <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr,
                            sex_het_depth = het_depth(c(F = FALSE, M = TRUE)))

  expect_true(all(c("Sex", "SexChr", "SexChrConsistent", "SexHetScore", "ChrType",
                    "ExpectedCN") %in% names(out)))
  sex <- unique(out[, c("SampleName", "Sex", "SexChr", "SexChrConsistent", "SexHetScore")])
  expect_equal(sex$Sex[sex$SampleName == "F"], "female")
  expect_equal(sex$Sex[sex$SampleName == "M"], "male")
  expect_equal(sex$SexChr, c("XX", "XY"))
  expect_equal(sex$SexHetScore, c("0/9", "6/9"))
  # Y present <=> PAR > X; Y absent <=> PAR == X
  expect_true(all(sex$SexChrConsistent))
  # Expected sex-chromosome dosages are not reported as deviations
  expect_true(all(is.na(out$Segment)))
  expect_true(all(out$Ploidy == 2L))
  expect_false("ChrY" %in% out$ChrID)
})

test_that("karyotype_notation XY flags X aneuploidy; Y copies come from PAR - X", {
  df  <- make_sex_df(list(K = list(X = 2L, PAR = 3L)))
  out <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr,
                            sex_het_depth = het_depth(c(K = TRUE)))

  expect_equal(unique(out$SexChr), "XXY")
  expect_equal(unique(out$Sex), "male")
  # A single Y is expected, so the extra X is flagged
  expect_false(is.na(out$Segment[out$ChrID == "ChrX"]))
  # Y has no CN call, so it has no row and Y aneuploidy cannot be reported
  expect_false("ChrY" %in% out$ChrID)
})

test_that("depth-based Y presence is only used when sex_chr defines Y", {
  # Depth says no Y; with Y defined this is "?" (PAR > X unexplained)
  df  <- make_sex_df(list(K = list(X = 2L, PAR = 3L)))
  no_y_dep <- het_depth(c(K = FALSE))
  with_y <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr,
                               sex_het_depth = no_y_dep)
  expect_equal(unique(with_y$SexChr), "?")
  expect_equal(unique(with_y$Sex), "?")

  # Without sex_chr$Y the score is ignored: Y = PAR - X from CN alone
  x_par <- list(X = "ChrX", PAR = "ChrPAR")
  expect_warning(
    no_y <- karyotype_notation(df, sex_system = "XY", sex_chr = x_par,
                               sex_het_depth = no_y_dep),
    "sex_het_depth` is ignored")
  expect_equal(unique(no_y$SexChr), "XXY")
  expect_equal(unique(no_y$Sex), "male")
  expect_equal(unique(no_y$SexHetCN), 1L)
  expect_true(all(is.na(no_y$SexHetScore)))
  expect_true(all(is.na(no_y$SexChrConsistent)))

  # Same for a score stored in params when sex_chr is overridden without Y (silent)
  hmm <- structure(list(by_window = NULL, by_marker = df,
                        params = c(list(sex_system = "XY", sex_chr = xy_chr),
                                   het_params(FALSE))),
                   class = "hmm_CN")
  expect_equal(unique(karyotype_notation(hmm)$SexChr), "?")
  # print() works on a minimal object and shows the stored depth score
  expect_output(print(hmm), "Y/W depth score: 0/9 \\(absent\\)")
  expect_output(print(hmm), "pi0\\): NA")
  expect_silent(stored <- karyotype_notation(hmm, sex_chr = x_par))
  expect_equal(unique(stored$SexChr), "XXY")
  expect_true(all(is.na(stored$SexHetScore)))
  # ...and the wrappers that forward the stored score stay silent too
  bw <- data.frame(Sample = df$SampleName, Chr = df$Chr, WindowID = 1L,
                   Start = df$Position, End = df$Position, CN_call = df$CN_call)
  hmm$by_window <- bw
  expect_silent(ct <- count_types(hmm, sex_chr = x_par))
  expect_equal(ct$sex$SexChr, "XXY")
  expect_silent(plot_karyotype(hmm, sample_name = "K", sex_chr = x_par))
  expect_silent(plot_karyotype(hmm, sample_name = "K", sex_system = "none", notation = TRUE))

  # Other sex systems ignore it as well
  expect_warning(karyotype_notation(df[df$Chr != "ChrPAR", ], sex_system = "X0",
                                    sex_chr = list(X = "ChrX"), sex_het_depth = no_y_dep),
                 "ignored")
})

test_that("karyotype_notation XY derives X and Y copies from PAR and ploidy", {
  # X-specific markers kept: Y = PAR - X
  df <- make_sex_df(list(M = list(X = 3L, PAR = 4L)), ploidy = 4L)
  dep <- het_depth(c(M = TRUE))
  with_x <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr, sex_het_depth = dep)
  expect_equal(unique(with_x$SexChr), "XXXY")
  expect_equal(unique(with_x$Sex), "male")
  expect_true(all(with_x$SexChrConsistent))

  # X-specific markers dropped too: PAR 4 = X + Y could be XXXY, XXYY or XYYY
  df_nox <- make_sex_df(list(M = list(PAR = 4L)), ploidy = 4L)
  no_x <- karyotype_notation(df_nox, sex_system = "XY",
                             sex_chr = list(Y = "ChrY", PAR = "ChrPAR"),
                             sex_het_depth = dep)
  expect_equal(unique(no_x$Sex), "?")
  expect_equal(unique(no_x$SexChr), "?")
  expect_true(all(is.na(no_x$SexHetCN)))
  expect_true(all(is.na(no_x$SexChrConsistent)))   # needs X to compare with PAR
  # same in a diploid: XY vs YY cannot be told apart without X
  no_x2 <- karyotype_notation(make_sex_df(list(M = list(PAR = 2L))), sex_system = "XY",
                              sex_chr = xy_chr, sex_het_depth = dep)
  expect_equal(unique(no_x2$SexChr), "?")

  # Absent Y: all PAR copies are X
  no_y <- karyotype_notation(df_nox, sex_system = "XY",
                             sex_chr = list(Y = "ChrY", PAR = "ChrPAR"),
                             sex_het_depth = het_depth(c(M = FALSE)))
  expect_equal(unique(no_y$SexChr), "XXXX")
  expect_equal(unique(no_y$Sex), "female")
})

test_that("SexChrConsistent compares Y depth presence with PAR vs X copy number", {
  consistent <- function(x, par, present, ploidy = 2L) {
    df <- make_sex_df(list(S = list(X = x, PAR = par)), ploidy = ploidy)
    out <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr,
                              sex_het_depth = het_depth(c(S = present)))
    list(flag = unique(out$SexChrConsistent), out = out)
  }

  # Y present and PAR > X: consistent
  expect_true(consistent(1L, 2L, TRUE)$flag)
  expect_true(consistent(3L, 4L, TRUE, ploidy = 4L)$flag)
  # Y absent and PAR == X: consistent
  expect_true(consistent(2L, 2L, FALSE)$flag)

  # Y present but PAR == X == ploidy: sex-chromosome aneuploidy (XXY), flagged
  # Y present but PAR == X: no PAR copy left for Y, so the complement is "?"
  # (to be investigated); X is still expected at ploidy - 1 and flagged
  for (p in c(2L, 4L)) {
    an <- consistent(p, p, TRUE, ploidy = p)
    expect_false(an$flag)
    expect_equal(unique(an$out$Sex), "?")
    expect_equal(unique(an$out$SexChr), "?")
    expect_equal(unique(an$out$SexHetCN), 1L)   # ploidy - X = 0, so at least 1
    expect_false(is.na(an$out$Segment[an$out$ChrID == "ChrX"]))
    expect_equal(an$out$ExpectedCN[an$out$ChrID == "ChrX"], p - 1L)
  }
  # PAR == X below the ploidy with Y present is unresolved too
  expect_equal(unique(consistent(3L, 3L, TRUE, ploidy = 4L)$out$SexChr), "?")

  # Y absent but PAR exceeds X: copies unexplained by X, and PAR < X is impossible
  expect_false(consistent(1L, 2L, FALSE)$flag)
  expect_false(consistent(3L, 2L, TRUE)$flag)
  expect_false(consistent(3L, 2L, FALSE)$flag)

  # Aneuploid classification follows the same rule
  mk <- function(s, x, par) data.frame(
    Sample = s, Chr = rep(c("Chr01", "Chr02", "ChrX", "ChrPAR"), each = 2),
    WindowID = rep(1:2, 4), Start = rep(c(1, 5e6), 4), End = rep(c(4e6, 9e6), 4),
    CN_call = rep(c(2L, 2L, x, par), each = 2), stringsAsFactors = FALSE)
  hmm <- structure(list(by_window = rbind(mk("XXY", 2L, 2L), mk("XY", 1L, 2L)),
                        params_samples = list(XXY = het_params(TRUE), XY = het_params(TRUE))),
                   class = "hmm_CN")
  ct <- count_types(hmm, sex_system = "XY", sex_chr = xy_chr)
  # Y present with PAR == X is "?": flagged only in $sex, not as aneuploid
  expect_setequal(ct$euploid, c("XXY", "XY"))
  expect_length(ct$aneuploid, 0)
  expect_equal(ct$sex$SexChr[ct$sex$Sample == "XXY"], "?")
  expect_equal(ct$sex$Sex[ct$sex$Sample == "XXY"], "?")

  # ...but mixed CN within a chromosome of a "?" sample is still segmental
  seg <- mk("XXY", 2L, 2L)
  seg$CN_call[seg$Chr == "ChrX"][1] <- 1L
  hmm$by_window <- rbind(seg, mk("XY", 1L, 2L))
  ct2 <- count_types(hmm, sex_system = "XY", sex_chr = xy_chr)
  expect_equal(ct2$segmental, "XXY")
  # and autosomal aneuploidy still counts
  an <- mk("XXY", 2L, 2L)
  an$CN_call[an$Chr == "Chr01"] <- 3L
  hmm$by_window <- rbind(an, mk("XY", 1L, 2L))
  expect_equal(count_types(hmm, sex_system = "XY", sex_chr = xy_chr)$aneuploid, "XXY")
})

test_that("karyotype_notation estimates Y copies from depth presence, PAR and X", {
  cases <- data.frame(
    # One row per example of the decision table (x = NA: X not called)
    present = c(rep(TRUE, 8), rep(FALSE, 8)),
    ploidy  = c(4L, 4L, 4L, 2L,   4L, 2L,   4L,   4L,     4L, 4L, 2L,   4L, 4L, 2L,   2L,   4L),
    par     = c(4L, 4L, 3L, 2L,   4L, 2L,   2L,   4L,     4L, 3L, 1L,   4L, 4L, 2L,   2L,   4L),
    x       = c(3L, 2L, 2L, 1L,   4L, 2L,   3L,   NA,     4L, 3L, 1L,   2L, 3L, 1L,   3L,   NA),
    sexchr  = c("XXXY", "XXYY", "XXY0", "XY",   # present, PAR > X
                "?", "?",                       # present, PAR = X (likely contamination)
                "?",                            # present, PAR < X (impossible)
                "?",                            # present, X not called
                "XXXX", "XXX0", "X0",           # absent,  PAR = X
                "?", "?", "?",                  # absent,  PAR > X
                "?",                            # absent,  PAR < X (impossible)
                "XXXX"),                        # absent,  X not called
    y_cn    = c(1L, 2L, 1L, 1L,   1L, 1L,   1L,   NA,     0L, 0L, 0L,   0L, 0L, 0L,   0L,   0L),
    consistent = c(TRUE, TRUE, TRUE, TRUE,  FALSE, FALSE,  FALSE,  NA,
                   TRUE, TRUE, TRUE,  FALSE, FALSE, FALSE,  FALSE,  NA),
    stringsAsFactors = FALSE)

  for (i in seq_len(nrow(cases))) {
    cs <- cases[i, ]
    x <- if (is.na(cs$x)) NULL else cs$x
    df <- make_sex_df(list(S = list(X = x, PAR = cs$par)), ploidy = cs$ploidy)
    out <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr,
                              sex_het_depth = het_depth(c(S = cs$present)))
    info <- sprintf("ploidy %d, PAR %d, X %s, Y present %s", cs$ploidy, cs$par, cs$x, cs$present)
    expect_equal(unique(out$SexChr), cs$sexchr, info = info)
    expect_equal(unique(out$SexHetCN), cs$y_cn, info = info)
    expect_equal(unique(out$SexChrConsistent), cs$consistent, info = info)
    # Sex is "?" whenever the complement is "?"
    expect_equal(unique(out$Sex),
                 if (cs$sexchr == "?") "?" else if (cs$present) "male" else "female",
                 info = info)
  }

  # Missing copies are flagged: XXX0 in a tetraploid (X and PAR both expected 4)
  df <- make_sex_df(list(S = list(X = 3L, PAR = 3L)), ploidy = 4L)
  out <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr,
                            sex_het_depth = het_depth(c(S = FALSE)))
  expect_equal(out$ExpectedCN[out$ChrID == "ChrX"], 4L)
  expect_false(is.na(out$Segment[out$ChrID == "ChrX"]))
  expect_false(is.na(out$Segment[out$ChrID == "ChrPAR"]))

  # count_types reports the estimate too
  bw <- data.frame(Sample = df$SampleName, Chr = df$Chr, WindowID = 1L,
                   Start = df$Position, End = df$Position, CN_call = df$CN_call)
  hmm <- structure(list(by_window = bw, params_samples = list(S = het_params(FALSE))),
                   class = "hmm_CN")
  ct <- count_types(hmm, sex_system = "XY", sex_chr = xy_chr)
  expect_equal(ct$sex$SexHetCN, 0L)
  expect_equal(ct$sex$SexChr, "XXX0")
})

test_that("segmental X is resolved with the Y depth evidence", {
  # X has CN `first` on its first 4 markers and the modal `rest` on the other 6
  seg_x <- function(first, rest, par, present) {
    df <- make_sex_df(list(S = list(X = rest, PAR = par)))
    df$CN_call[df$Chr == "ChrX"][1:4] <- first
    karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr,
                       sex_het_depth = het_depth(c(S = present)))
  }
  # Y absent, PAR 2, X 2 on one arm and 1 on the larger rest: XX female with
  # a partial X deletion (modal X = 1 alone would give "?")
  del <- seg_x(2L, 1L, 2L, FALSE)
  expect_equal(unique(del$Sex), "female")
  expect_equal(unique(del$SexChr), "XX")
  expect_true(unique(del$SexChrConsistent))
  expect_equal(unique(del$SexHetCN), 0L)
  expect_equal(del$ExpectedCN[del$ChrID == "ChrX"], 2L)
  expect_match(del$Segment[del$ChrID == "ChrX"], "×1$")   # the deleted region

  # Y present, PAR 2, X 1 on one arm and 2 on the larger rest: XY male with a
  # partial X duplication (modal X = 2 alone would give "?")
  dup <- seg_x(1L, 2L, 2L, TRUE)
  expect_equal(unique(dup$Sex), "male")
  expect_equal(unique(dup$SexChr), "XY")
  expect_true(unique(dup$SexChrConsistent))
  expect_equal(dup$ExpectedCN[dup$ChrID == "ChrX"], 1L)
  expect_match(dup$Segment[dup$ChrID == "ChrX"], "×2$")   # the duplicated region

  # No X segment agrees with the depth evidence: still "?"
  expect_equal(unique(seg_x(3L, 1L, 2L, FALSE)$SexChr), "?")
  # Without a depth score the modal X keeps the CN-based inference (Y = PAR - X)
  df <- make_sex_df(list(S = list(X = 1L, PAR = 2L)))
  df$CN_call[df$Chr == "ChrX"][1:4] <- 2L
  expect_equal(unique(karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr)$SexChr), "XY")

  # count_types: resolved sex, sample classed as segmental
  bw <- data.frame(Sample = "S", Chr = df$Chr, WindowID = seq_len(nrow(df)),
                   Start = df$Position, End = df$Position + 1, CN_call = df$CN_call)
  hmm <- structure(list(by_window = bw, params_samples = list(S = het_params(FALSE))),
                   class = "hmm_CN")
  ct <- count_types(hmm, sex_system = "XY", sex_chr = xy_chr)
  expect_equal(ct$sex$SexChr, "XX")
  expect_equal(ct$segmental, "S")
})

test_that("karyotype_notation XY cannot infer sex without a Y depth score", {
  # Neither X nor Y specific markers survive: PAR alone says nothing about sex
  df <- make_sex_df(list(F = list(PAR = 2L), M = list(PAR = 2L)))
  out <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr)
  expect_true(all(is.na(out$Sex)))
  expect_true(all(is.na(out$SexHetScore)))
})

test_that("karyotype_notation falls back to CN-called Y when there is no depth score", {
  # Legacy / non-standardized input where Y-specific markers were kept and called
  df <- make_sex_df(list(F = list(X = 2L, Y = 0L, PAR = 2L),
                         M = list(X = 1L, Y = 1L, PAR = 2L)))
  out <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr)
  sex <- unique(out[, c("SampleName", "Sex", "SexChr", "SexChrConsistent")])
  expect_equal(sex$Sex, c("female", "male"))
  expect_equal(sex$SexChr, c("XX", "XY"))
  expect_true(all(sex$SexChrConsistent))
  expect_true(all(is.na(out$Segment)))

  bad <- make_sex_df(list(M = list(X = 1L, Y = 1L, PAR = 3L)))
  expect_false(unique(karyotype_notation(bad, sex_system = "XY", sex_chr = xy_chr)$SexChrConsistent))
})

test_that("karyotype_notation ZW: W present is female (W from depth)", {
  # No ChrW rows: W-specific markers are dropped by standardization
  df  <- make_sex_df(list(F = list(Z = 1L, PAR = 2L),
                          M = list(Z = 2L, PAR = 2L)),
                     hom = "Z", het = "W")
  expect_false("ChrW" %in% df$Chr)
  out <- karyotype_notation(df, sex_system = "ZW",
                            sex_chr = list(Z = "ChrZ", W = "ChrW", PAR = "ChrPAR"),
                            sex_het_depth = het_depth(c(F = TRUE, M = FALSE)))
  sex <- unique(out[, c("SampleName", "Sex", "SexChr")])
  expect_equal(sex$Sex[sex$SampleName == "F"], "female")
  expect_equal(sex$SexChr[sex$SampleName == "F"], "ZW")
  expect_equal(sex$Sex[sex$SampleName == "M"], "male")
  expect_equal(sex$SexChr[sex$SampleName == "M"], "ZZ")
  expect_true(all(is.na(out$Segment)))
})

test_that("karyotype_notation X0 infers sex from X dosage", {
  df <- make_sex_df(list(F = list(X = 2L, PAR = 2L),
                         M = list(X = 1L, PAR = 2L)))
  df <- df[df$Chr %in% c("Chr01", "Chr02", "ChrX"), ]
  out <- karyotype_notation(df, sex_system = "X0", sex_chr = list(X = "ChrX"))
  sex <- unique(out[, c("SampleName", "Sex", "SexChr")])
  expect_equal(sex$SexChr[sex$SampleName == "M"], "X0")
  expect_equal(sex$Sex[sex$SampleName == "M"], "male")
  expect_equal(sex$Sex[sex$SampleName == "F"], "female")
  expect_true(all(is.na(out$Segment)))
})

test_that("karyotype_notation haplodiploid infers sex from ploidy", {
  df <- make_marker_df()
  df$CN_call[df$SampleName == "S1"] <- 1L
  df$CN_call[df$SampleName == "S2"] <- 2L
  out <- karyotype_notation(df, sex_system = "haplodiploid")
  expect_equal(unique(out$Sex[out$SampleName == "S1"]), "male")
  expect_equal(unique(out$Sex[out$SampleName == "S2"]), "female")
})

test_that("karyotype_notation validates sex_chr", {
  df <- make_sex_df(list(M = list(X = 1L, PAR = 2L)))
  expect_error(karyotype_notation(df, sex_system = "XY", sex_chr = list(Z = "ChrX")),
               "Invalid")
  expect_error(karyotype_notation(df, sex_system = "XY", sex_chr = list(X = "ChrX")),
               "requires")
  expect_error(karyotype_notation(df, sex_system = "XY",
                                  sex_chr = list(X = "ChrX", Y = "ChrX")),
               "more than one")
  expect_error(karyotype_notation(df, sex_system = "XY"), "named list")
})

test_that("count_types does not count XY males as aneuploid", {
  # Two windows per chromosome; F = XX female, M = XY male, K = XXY diploid
  mk <- function(s, cn) data.frame(
    Sample = s, Chr = rep(names(cn), each = 2), WindowID = rep(1:2, length(cn)),
    Start = rep(c(1, 5e6), length(cn)), End = rep(c(4e6, 9e6), length(cn)),
    CN_call = rep(unlist(cn), each = 2), stringsAsFactors = FALSE)
  # No ChrY windows (Y-specific markers are dropped); Y presence comes from depth
  bw <- rbind(
    mk("F", list(Chr01 = 2L, Chr02 = 2L, ChrX = 2L, ChrPAR = 2L)),
    mk("M", list(Chr01 = 2L, Chr02 = 2L, ChrX = 1L, ChrPAR = 2L)),
    mk("K", list(Chr01 = 2L, Chr02 = 2L, ChrX = 2L, ChrPAR = 3L)))
  hmm <- structure(list(by_window = bw,
                        params_samples = list(F = het_params(FALSE), M = het_params(TRUE),
                                              K = het_params(TRUE))),
                   class = "hmm_CN")

  plain <- count_types(hmm)
  expect_true("M" %in% plain$aneuploid)

  sexed <- count_types(hmm, sex_system = "XY", sex_chr = xy_chr)
  expect_setequal(sexed$euploid, c("F", "M"))
  expect_equal(sexed$aneuploid, "K")
  expect_equal(sexed$sex$Sex[match(c("F", "M", "K"), sexed$sex$Sample)],
               c("female", "male", "male"))
  expect_equal(sexed$sex$SexHetScore[match(c("F", "M"), sexed$sex$Sample)], c("0/9", "6/9"))
  expect_output(print(sexed), "Sex: female 1 \\| male 2")
})

test_that("sex settings stored in hmm_CN params are used by default", {
  df <- make_sex_df(list(F = list(X = 2L, PAR = 2L),
                         M = list(X = 1L, PAR = 2L)))
  bw <- data.frame(Sample = df$SampleName, Chr = df$Chr, WindowID = 1L,
                   Start = df$Position, End = df$Position, CN_call = df$CN_call)
  p  <- list(sex_system = "XY", sex_chr = xy_chr)
  hmm <- structure(list(by_window = bw, by_marker = df,
                        params_samples = list(F = c(p, het_params(FALSE)),
                                              M = c(p, het_params(TRUE)))),
                   class = "hmm_CN")

  kn <- karyotype_notation(hmm)
  expect_true("Sex" %in% names(kn))
  expect_equal(unique(kn$Sex[kn$SampleName == "M"]), "male")

  ct <- count_types(hmm)
  expect_setequal(ct$euploid, c("F", "M"))

  # Explicit argument overrides the stored setting
  expect_false("Sex" %in% names(karyotype_notation(hmm, sex_system = "none")))
  expect_true("M" %in% count_types(hmm, sex_system = "none")$aneuploid)

  # No stored setting keeps the previous behaviour
  hmm$params_samples <- list(F = list(), M = list())
  expect_false("Sex" %in% names(karyotype_notation(hmm)))
})

test_that(".sex_het_depth scores Y/W presence from raw depth", {
  d <- data.frame(MarkerName = c(paste0("y", 1:5), "a1"),
                  Chr = c(rep("ChrY", 5), "Chr01"),
                  R = c(10, 7, 6, NA, 0, 50))
  h <- .sex_het_depth(d, "XY", xy_chr, sex_min_depth = 6, sex_min_markers = 2)
  expect_equal(h$n_covered, 2)            # depth strictly above 6; NA and 0 uncovered
  expect_equal(h$n_markers, 5)            # autosomal marker not counted
  expect_equal(h$score, "2/5")
  expect_true(h$present)
  expect_false(.sex_het_depth(d, "XY", xy_chr, 6, sex_min_markers = 3)$present)
  expect_true(.sex_het_depth(d, "XY", xy_chr, sex_min_depth = 5, 3)$present)

  # Not applicable / no evidence
  expect_null(.sex_het_depth(d, "none", NULL))
  expect_null(.sex_het_depth(d, "XY", list(X = "ChrX", PAR = "ChrPAR")))
  expect_null(.sex_het_depth(d[, c("MarkerName", "Chr")], "XY", xy_chr))
  none <- .sex_het_depth(d[d$Chr != "ChrY", ], "XY", xy_chr)
  expect_equal(none$score, "0/0")
  expect_true(is.na(none$present))
  # ZW uses W
  zw <- .sex_het_depth(transform(d, Chr = sub("ChrY", "ChrW", Chr)), "ZW",
                       list(Z = "ChrZ", W = "ChrW", PAR = "ChrPAR"))
  expect_equal(zw$score, "2/5")

  expect_error(.check_sex_depth_args(-1, 2), "sex_min_depth")
  expect_error(.check_sex_depth_args(6, 0), "sex_min_markers")
  expect_error(.check_sex_depth_args(6, 1.5), "sex_min_markers")
})

test_that("depth-based Y presence overrides CN calls and handles sample subsets", {
  # Tetraploids; X-specific markers kept, Y-specific ones dropped
  df <- make_sex_df(list(F = list(X = 4L, PAR = 4L),
                         M = list(X = 3L, PAR = 4L)), ploidy = 4L)
  dep <- data.frame(SampleName = c("F", "M"), present = c(FALSE, TRUE), score = c("0/9", "5/9"))
  out <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr, sex_het_depth = dep)
  sex <- unique(out[, c("SampleName", "Sex", "SexChr", "SexHetScore")])
  expect_equal(sex$Sex, c("female", "male"))
  expect_equal(sex$SexChr, c("XXXX", "XXXY"))
  expect_equal(sex$SexHetScore, c("0/9", "5/9"))
  expect_true(all(is.na(out$Segment)))

  # More than one Y copy: Y = PAR - X (XXYY); no score column value without a score
  df2 <- make_sex_df(list(M = list(X = 2L, PAR = 4L)), ploidy = 4L)
  out2 <- karyotype_notation(df2, sex_system = "XY", sex_chr = xy_chr,
                             sex_het_depth = data.frame(SampleName = "M", present = TRUE))
  expect_equal(unique(out2$SexChr), "XXYY")
  expect_true(all(is.na(out2$SexHetScore)))

  # X monosomy with no Y: PAR - X alone would imply one Y (male); depth says Y is absent
  df5 <- make_sex_df(list(T = list(X = 1L, PAR = 2L)))
  expect_equal(unique(karyotype_notation(df5, sex_system = "XY", sex_chr = xy_chr)$Sex), "male")
  out5 <- karyotype_notation(df5, sex_system = "XY", sex_chr = xy_chr,
                             sex_het_depth = het_depth(c(T = FALSE)))
  expect_equal(unique(out5$Sex), "?")
  expect_equal(unique(out5$SexChr), "?")      # the 2nd PAR copy is unexplained
  expect_equal(unique(out5$SexHetCN), 0L)
  expect_false(is.na(out5$Segment[out5$ChrID == "ChrX"]))   # X deficit is flagged

  # Depth wins over a (legacy) contradicting Y CN call
  df3 <- make_sex_df(list(S = list(X = 1L, Y = 1L, PAR = 2L)))
  out3 <- karyotype_notation(df3, sex_system = "XY", sex_chr = xy_chr,
                             sex_het_depth = het_depth(c(S = FALSE)))
  # no Y by depth, so PAR > X is unexplained: complement and sex are "?"
  expect_equal(unique(out3$SexHetCN), 0L)
  expect_equal(unique(out3$Sex), "?")
  # with Y absent and PAR == X, depth makes it female despite the Y call
  df3b <- make_sex_df(list(S = list(X = 2L, Y = 1L, PAR = 2L)))
  out3b <- karyotype_notation(df3b, sex_system = "XY", sex_chr = xy_chr,
                              sex_het_depth = het_depth(c(S = FALSE)))
  expect_equal(unique(out3b$Sex), "female")

  # Samples without a score keep the CN-based inference
  out4 <- karyotype_notation(df3, sex_system = "XY", sex_chr = xy_chr,
                             sex_het_depth = het_depth(c(other = FALSE)))
  expect_equal(unique(out4$Sex), "male")

  expect_error(karyotype_notation(df3, sex_system = "XY", sex_chr = xy_chr,
                                  sex_het_depth = data.frame(Sample = "S")), "SampleName")
})

test_that("depth-based Y presence stored in params reaches all sex-aware functions", {
  df <- make_sex_df(list(F = list(X = 2L, PAR = 2L),
                         M = list(X = 1L, PAR = 2L)))
  bw <- data.frame(Sample = df$SampleName, Chr = df$Chr, WindowID = 1L,
                   Start = df$Position, End = df$Position, CN_call = df$CN_call)
  het <- function(p, s) c(list(sex_system = "XY", sex_chr = xy_chr), het_params(p, s))
  hmm <- structure(list(by_window = bw, by_marker = df,
                        params_samples = list(F = het(FALSE, "1/9"), M = het(TRUE, "7/9"))),
                   class = "hmm_CN")

  kn <- karyotype_notation(hmm)
  expect_equal(unique(kn$Sex[kn$SampleName == "M"]), "male")
  expect_equal(unique(kn$SexHetScore[kn$SampleName == "M"]), "7/9")
  expect_false("ChrY" %in% kn$ChrID)

  ct <- count_types(hmm)
  expect_equal(ct$sex$Sex[ct$sex$Sample == "M"], "male")
  expect_equal(ct$sex$SexHetScore[ct$sex$Sample == "F"], "1/9")
  expect_setequal(ct$euploid, c("F", "M"))

  # single-sample hmm_CN (params, not params_samples)
  one <- structure(list(by_window = bw[bw$Sample == "M", ], by_marker = df[df$SampleName == "M", ],
                        params = het(TRUE, "7/9")), class = "hmm_CN")
  expect_match(plot_karyotype(one, sample_name = "M")$labels$subtitle, "Sex: male \\(XY\\)")
})

test_that("compare_cn_track_summary adds a Sex column from stored settings", {
  df <- make_sex_df(list(F = list(X = 2L, PAR = 2L),
                         M = list(X = 1L, PAR = 2L)))
  bw <- data.frame(Sample = df$SampleName, Chr = df$Chr, WindowID = 1L,
                   Start = df$Position, End = df$Position + 1, CN_call = df$CN_call)
  p  <- list(sex_system = "XY", sex_chr = xy_chr)
  hmm <- structure(list(by_window = bw,
                        params_samples = list(F = c(p, het_params(FALSE)),
                                              M = c(p, het_params(TRUE)))),
                   class = "hmm_CN")

  has_sex_layer <- function(g) any(vapply(g$layers, function(l)
    is.data.frame(l$data) && "lbl" %in% names(l$data) && any(grepl("^M \\(XY\\)$", l$data$lbl)),
    logical(1)))
  block_labels <- function(g) unlist(lapply(g$layers, function(l)
    if (is.data.frame(l$data) && "mid_y" %in% names(l$data)) l$data$lbl))

  g <- compare_cn_track_summary(hmm)
  expect_true(has_sex_layer(g))
  expect_false(any(grepl("aneuploid", block_labels(g))))

  g_off <- compare_cn_track_summary(hmm, sex_system = NULL)
  expect_false(has_sex_layer(g_off))
  expect_true(any(grepl("aneuploid", block_labels(g_off))))
})

test_that("count_types, plot_karyotype and compare_cn_track_summary report Y depth checks", {
  # Tetraploids: M = XXYY (Y present, PAR > X), C = Y present but PAR == X (?)
  df <- make_sex_df(list(M = list(X = 2L, PAR = 4L),
                         C = list(X = 4L, PAR = 4L)), ploidy = 4L)
  bw <- data.frame(Sample = df$SampleName, Chr = df$Chr, WindowID = 1L,
                   Start = df$Position, End = df$Position + 1, CN_call = df$CN_call)
  het <- function(s) c(list(sex_system = "XY", sex_chr = xy_chr), het_params(TRUE, s))
  hmm <- structure(list(by_window = bw, by_marker = df,
                        params_samples = list(M = het("5/9"), C = het("7/9"))),
                   class = "hmm_CN")

  # count_types: SexChrConsistent reported; the "?" sample is only flagged in $sex
  ct <- count_types(hmm)
  sx <- ct$sex[match(c("M", "C"), ct$sex$Sample), ]
  expect_equal(sx$SexChr, c("XXYY", "?"))
  expect_equal(sx$SexChrConsistent, c(TRUE, FALSE))
  expect_equal(sx$SexHetCN, c(2L, 1L))
  expect_setequal(ct$euploid, c("M", "C"))
  expect_length(ct$aneuploid, 0)
  expect_equal(sx$Sex, c("male", "?"))
  expect_output(print(ct), "Sex: female 0 \\| male 1 \\| \\? 1 \\| unknown 0")
  expect_output(print(ct), "Sex chromosomes to check: 1")

  # plot_karyotype: depth score in the subtitle, check flag only for C
  sub_m <- plot_karyotype(hmm, sample_name = "M")$labels$subtitle
  sub_c <- plot_karyotype(hmm, sample_name = "C")$labels$subtitle
  expect_match(sub_m, "Sex: male \\(XXYY\\), Y markers covered 5/9")
  expect_false(grepl("check", sub_m))
  expect_match(sub_c, "Sex: \\? \\(\\?\\), Y markers covered 7/9 \\[check sex chromosomes\\]")

  # compare_cn_track_summary: C's label in red, explained in the caption
  g <- compare_cn_track_summary(hmm)
  lyr <- Filter(function(l) is.data.frame(l$data) && "check" %in% names(l$data), g$layers)[[1]]
  expect_equal(lyr$data$check[match(c("M (XXYY)", "? (?)"), lyr$data$lbl)], c(FALSE, TRUE))
  expect_equal(lyr$aes_params$colour[lyr$data$lbl == "? (?)"], "#cb181d")
  expect_match(g$labels$caption, "Red Sex labels")
  # no caption when nothing needs checking
  ok <- hmm; ok$by_window <- bw[bw$Sample == "M", ]; ok$params_samples <- ok$params_samples["M"]
  expect_null(compare_cn_track_summary(ok)$labels$caption)
})

test_that("count_types and compare_cn_track_summary follow karyotype_notation for sex", {
  # PAR: one 8-marker window at CN 2 and two 1-marker windows at CN 3. Per
  # marker PAR = 2 (XY); counting windows would give PAR = 3 (XYY)
  win <- data.frame(Chr = c("Chr01", "Chr02", "ChrX", "ChrPAR", "ChrPAR", "ChrPAR"),
                    CN = c(2L, 2L, 1L, 2L, 3L, 3L), n = c(10, 10, 10, 8, 1, 1))
  bm <- do.call(rbind, lapply(seq_len(nrow(win)), function(i) data.frame(
    SampleName = "S", Chr = win$Chr[i], Position = i * 1e7 + seq_len(win$n[i]) * 1e3,
    CN_call = win$CN[i], stringsAsFactors = FALSE)))
  bw <- do.call(rbind, lapply(seq_len(nrow(win)), function(i) {
    p <- bm$Position[bm$Chr == win$Chr[i] & bm$Position %/% 1e7 == i]
    data.frame(Sample = "S", Chr = win$Chr[i], WindowID = i, Start = min(p), End = max(p),
               n_snps = win$n[i], CN_call = win$CN[i], stringsAsFactors = FALSE)
  }))
  p <- list(sex_system = "XY", sex_chr = xy_chr)
  hmm <- structure(list(by_window = bw, by_marker = bm,
                        params_samples = list(S = c(p, het_params(TRUE, "5/9")))),
                   class = "hmm_CN")

  kn <- karyotype_notation(hmm)
  expect_equal(unique(kn$SexChr), "XY")
  ct <- count_types(hmm)
  expect_equal(ct$sex$SexChr, unique(kn$SexChr))
  expect_equal(ct$sex$Sex, unique(kn$Sex))
  g <- compare_cn_track_summary(hmm)
  lyr <- Filter(function(l) is.data.frame(l$data) && "check" %in% names(l$data), g$layers)[[1]]
  expect_equal(lyr$data$lbl, "M (XY)")

  # Without by_marker, windows are weighted by n_snps
  no_bm <- hmm; no_bm$by_marker <- NULL
  expect_equal(count_types(no_bm)$sex$SexChr, "XY")

  # Samples missing from by_marker fall back to their windows (no error)
  two <- hmm
  two$by_window <- rbind(bw, transform(bw, Sample = "T"))
  two$params_samples$T <- c(p, het_params(TRUE, "5/9"))
  ct2 <- count_types(two)
  expect_equal(ct2$sex$SexChr[match(c("S", "T"), ct2$sex$Sample)], c("XY", "XY"))
  only_t <- compare_cn_track_summary(two, samples_to_plot = "T")
  lyr <- Filter(function(l) is.data.frame(l$data) && "check" %in% names(l$data),
                only_t$layers)[[1]]
  expect_equal(lyr$data$lbl, "M (XY)")
})

test_that("plot_karyotype supports sex chromosomes", {
  df <- make_sex_df(list(M = list(X = 1L, PAR = 2L)))
  # Y presence comes from the depth score stored in the hmm_CN params
  hmm <- structure(list(by_window = NULL, by_marker = df,
                        params = c(list(sex_system = "XY", sex_chr = xy_chr), het_params(TRUE))),
                   class = "hmm_CN")
  p  <- plot_karyotype(hmm, sample_name = "M")
  expect_s3_class(p, "gg")
  expect_match(p$labels$subtitle, "Sex: male \\(XY\\)")
  # No ghost second copy for the male X
  ghost_x <- vapply(p$layers, function(l) {
    d <- l$data
    is.data.frame(d) && all(c("grey_type", "Chr") %in% names(d)) &&
      any(d$grey_type == "Expected (ghost)" & d$Chr == "ChrX")
  }, logical(1))
  expect_false(any(ghost_x))
})

