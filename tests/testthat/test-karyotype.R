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

