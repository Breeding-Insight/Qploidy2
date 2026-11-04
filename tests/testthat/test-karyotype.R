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

# Diploid samples with 2 autosomes plus X, Y and PAR contigs.
# cn is a named list per sample: list(X = , Y = , PAR = )
make_sex_df <- function(cn, hom = "X", het = "Y") {
  chrs <- c("Chr01", "Chr02", paste0("Chr", hom), paste0("Chr", het), "ChrPAR")
  do.call(rbind, lapply(names(cn), function(s) {
    v <- cn[[s]]
    data.frame(
      SampleName = s,
      Chr        = rep(chrs, each = 10),
      Position   = rep(seq(1e6, 10e6, length.out = 10), length(chrs)),
      CN_call    = rep(c(2L, 2L, v[[hom]], v[[het]], v$PAR), each = 10),
      stringsAsFactors = FALSE
    )
  }))
}
xy_chr <- list(X = "ChrX", Y = "ChrY", PAR = "ChrPAR")

test_that("karyotype_notation XY infers sex and complement", {
  df  <- make_sex_df(list(F = list(X = 2L, Y = 0L, PAR = 2L),
                          M = list(X = 1L, Y = 1L, PAR = 2L)))
  out <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr)

  expect_true(all(c("Sex", "SexChr", "SexChrConsistent", "ChrType", "ExpectedCN") %in% names(out)))
  sex <- unique(out[, c("SampleName", "Sex", "SexChr", "SexChrConsistent")])
  expect_equal(sex$Sex[sex$SampleName == "F"], "female")
  expect_equal(sex$Sex[sex$SampleName == "M"], "male")
  expect_equal(sex$SexChr[sex$SampleName == "M"], "XY")
  expect_true(all(sex$SexChrConsistent))
  # Expected sex-chromosome dosages are not reported as deviations
  expect_true(all(is.na(out$Segment)))
  expect_true(all(out$Ploidy == 2L))
})

test_that("karyotype_notation XY flags sex chromosome aneuploidy", {
  df  <- make_sex_df(list(K = list(X = 2L, Y = 1L, PAR = 3L)))
  out <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr)

  expect_equal(unique(out$SexChr), "XXY")
  expect_equal(unique(out$Sex), "male")
  expect_false(is.na(out$Segment[out$ChrID == "ChrX"]))
  expect_false(is.na(out$Segment[out$ChrID == "ChrY"]))
})

test_that("karyotype_notation XY derives missing groups from PAR", {
  df <- make_sex_df(list(M = list(X = 3L, Y = 1L, PAR = 4L)))

  no_y <- karyotype_notation(df[df$Chr != "ChrY", ], sex_system = "XY",
                             sex_chr = list(X = "ChrX", PAR = "ChrPAR"))
  expect_equal(unique(no_y$SexChr), "XXXY")
  expect_equal(unique(no_y$Sex), "male")
  expect_true(all(is.na(no_y$SexChrConsistent)))

  no_x <- karyotype_notation(df[df$Chr != "ChrX", ], sex_system = "XY",
                             sex_chr = list(Y = "ChrY", PAR = "ChrPAR"))
  expect_equal(unique(no_x$SexChr), "XXXY")
})

test_that("karyotype_notation XY reports PAR inconsistency", {
  df  <- make_sex_df(list(M = list(X = 1L, Y = 1L, PAR = 3L)))
  out <- karyotype_notation(df, sex_system = "XY", sex_chr = xy_chr)
  expect_false(unique(out$SexChrConsistent))
})

test_that("karyotype_notation ZW: W present is female", {
  df  <- make_sex_df(list(F = list(Z = 1L, W = 1L, PAR = 2L),
                          M = list(Z = 2L, W = 0L, PAR = 2L)),
                     hom = "Z", het = "W")
  out <- karyotype_notation(df, sex_system = "ZW",
                            sex_chr = list(Z = "ChrZ", W = "ChrW", PAR = "ChrPAR"))
  sex <- unique(out[, c("SampleName", "Sex", "SexChr")])
  expect_equal(sex$Sex[sex$SampleName == "F"], "female")
  expect_equal(sex$SexChr[sex$SampleName == "F"], "ZW")
  expect_equal(sex$Sex[sex$SampleName == "M"], "male")
})

test_that("karyotype_notation X0 infers sex from X dosage", {
  df <- make_sex_df(list(F = list(X = 2L, Y = 0L, PAR = 2L),
                         M = list(X = 1L, Y = 0L, PAR = 2L)))
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
  df <- make_sex_df(list(M = list(X = 1L, Y = 1L, PAR = 2L)))
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
  bw <- rbind(
    mk("F", list(Chr01 = 2L, Chr02 = 2L, ChrX = 2L, ChrY = 0L, ChrPAR = 2L)),
    mk("M", list(Chr01 = 2L, Chr02 = 2L, ChrX = 1L, ChrY = 1L, ChrPAR = 2L)),
    mk("K", list(Chr01 = 2L, Chr02 = 2L, ChrX = 2L, ChrY = 1L, ChrPAR = 3L)))
  hmm <- structure(list(by_window = bw), class = "hmm_CN")

  plain <- count_types(hmm)
  expect_true("M" %in% plain$aneuploid)

  sexed <- count_types(hmm, sex_system = "XY", sex_chr = xy_chr)
  expect_setequal(sexed$euploid, c("F", "M"))
  expect_equal(sexed$aneuploid, "K")
  expect_equal(sexed$sex$Sex[match(c("F", "M", "K"), sexed$sex$Sample)],
               c("female", "male", "male"))
  expect_output(print(sexed), "Sex: female 1 \\| male 2")
})

test_that("sex settings stored in hmm_CN params are used by default", {
  df <- make_sex_df(list(F = list(X = 2L, Y = 0L, PAR = 2L),
                         M = list(X = 1L, Y = 1L, PAR = 2L)))
  bw <- data.frame(Sample = df$SampleName, Chr = df$Chr, WindowID = 1L,
                   Start = df$Position, End = df$Position, CN_call = df$CN_call)
  p  <- list(sex_system = "XY", sex_chr = xy_chr)
  hmm <- structure(list(by_window = bw, by_marker = df,
                        params_samples = list(F = p, M = p)), class = "hmm_CN")

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

test_that("plot_karyotype supports sex chromosomes", {
  df <- make_sex_df(list(M = list(X = 1L, Y = 1L, PAR = 2L)))
  p  <- plot_karyotype(df, sample_name = "M", sex_system = "XY", sex_chr = xy_chr)
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

