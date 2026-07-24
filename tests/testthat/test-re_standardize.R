library(testthat)
library(Qploidy2)

# ── Shared fixtures ────────────────────────────────────────────────────────────
# Build once and reuse across all tests in this file.

local_fixtures <- local({
  env <- new.env(parent = emptyenv())

  setup <- function() {
    if (!is.null(env$multi_res)) return(env)   # already built

    set.seed(42)
    vcf_file <- tempfile(fileext = ".vcf")
    simulate_vcf(
      seed = 42, file_path = vcf_file,
      n_tetraploid = 10, n_diploid = 2, n_triploid = 2,
      n_markers = 60
    )

    env$data     <- qploidy_read_vcf(vcf_file)
    env$genos    <- qploidy_read_vcf(vcf_file, geno = TRUE)
    env$genos    <- env$genos[grep("Tetraploid", env$genos$SampleName), ]
    env$geno_pos <- qploidy_read_vcf(vcf_file, geno.pos = TRUE)

    env$std <- standardize(
      data                  = env$data,
      genos                 = env$genos,
      geno.pos              = env$geno_pos,
      ploidy.standardization = 4,
      threshold.n.clusters  = 3,
      n.cores               = 1,
      verbose               = FALSE
    )

    env$multi_res <- hmm_estimate_CN_multi(
      qploidy_standarize_result = env$std,
      sample_ids = unique(env$std$data$SampleName),
      n_cores    = 1,
      chr        = 1,
      snps_per_window     = 15,
      min_snps_per_window = 5,
      cn_grid             = c(2, 3, 4),
      M                   = 21,
      exp_ploidy          = 4,
      verbose             = FALSE
    )

    env
  }

  setup
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("re_standardize errors when hmm_CN_multi is not hmm_CN class", {
  fixtures <- local_fixtures()
  expect_error(
    re_standardize(
      data           = fixtures$data,
      geno.pos       = fixtures$geno_pos,
      hmm_CN_multi   = list(x = 1),           # wrong class
      ploidy.standardization = 4,
      verbose        = FALSE
    ),
    "class 'hmm_CN'"
  )
})

test_that("re_standardize errors when hmm_CN_multi lacks params_samples (single-sample object)", {
  fixtures <- local_fixtures()

  # hmm_estimate_CN returns a single-sample hmm_CN without params_samples
  single <- hmm_estimate_CN(
    qploidy_standarize_result = fixtures$std,
    sample_id    = unique(fixtures$std$data$SampleName)[1],
    chr          = 1,
    snps_per_window     = 15,
    min_snps_per_window = 5,
    cn_grid      = c(2, 3, 4),
    M            = 21,
    verbose      = FALSE
  )

  expect_error(
    re_standardize(
      data           = fixtures$data,
      geno.pos       = fixtures$geno_pos,
      hmm_CN_multi   = single,
      ploidy.standardization = 4,
      verbose        = FALSE
    ),
    "hmm_estimate_CN_multi"
  )
})

# ── Return-value structure ─────────────────────────────────────────────────────

test_that("re_standardize returns a valid qploidy_standardization object", {
  fixtures <- local_fixtures()

  result <- re_standardize(
    data                   = fixtures$data,
    geno.pos               = fixtures$geno_pos,
    hmm_CN_multi           = fixtures$multi_res,
    ploidy.standardization = 4,
    threshold.n.clusters   = 3,
    n.cores                = 1,
    threshold.geno.prob    = 0.5,
    threshold.missing.geno = 0.9,
    verbose                = FALSE
  )

  expect_s3_class(result, "qploidy_standardization")
  expect_named(result, c("info", "filters", "data"))
  expect_true(is.data.frame(result$data) || inherits(result$data, "tbl_df"))
  expect_true(nrow(result$data) > 0)

  multi <- hmm_estimate_CN_multi(qploidy_standarize_result = result,
                                 data = fixtures$data,
                                 geno.pos = fixtures$geno_pos)
  pl <- summarize_cn_mode(multi, level = "sample")
  expect_equal(as.vector(table(pl$CN_mode)), c(2,2,10))

})

test_that("re_standardize output contains expected columns", {
  fixtures <- local_fixtures()

  result <- re_standardize(
    data                   = fixtures$data,
    geno.pos               = fixtures$geno_pos,
    hmm_CN_multi           = fixtures$multi_res,
    ploidy.standardization = 4,
    n.cores                = 1,
    verbose                = FALSE
  )

  expected_cols <- c("MarkerName", "SampleName", "baf", "z", "Chr", "Position")
  expect_true(all(expected_cols %in% names(result$data)))
})

test_that("re_standardize filters slot has required names", {
  fixtures <- local_fixtures()

  result <- re_standardize(
    data                   = fixtures$data,
    geno.pos               = fixtures$geno_pos,
    hmm_CN_multi           = fixtures$multi_res,
    ploidy.standardization = 4,
    n.cores                = 1,
    verbose                = FALSE
  )

  required_filter_names <- c(
    "n.markers.start", "miss.rm", "clusters.rm",
    "no.geno.info.rm", "n.markers.end"
  )
  expect_true(all(required_filter_names %in% names(result$filters)))
  expect_true(result$filters["n.markers.start"] >= result$filters["n.markers.end"])
})

test_that("re_standardize info slot reflects supplied ploidy", {
  fixtures <- local_fixtures()

  result <- re_standardize(
    data                   = fixtures$data,
    geno.pos               = fixtures$geno_pos,
    hmm_CN_multi           = fixtures$multi_res,
    ploidy.standardization = 4,
    n.cores                = 1,
    verbose                = FALSE
  )

  expect_equal(as.numeric(result$info["ploidy.standardization"]), 4)
})

test_that("re_standardize uses mode of CN_call when ploidy.standardization is NULL", {
  fixtures <- local_fixtures()

  result <- suppressMessages(re_standardize(
    data                   = fixtures$data,
    geno.pos               = fixtures$geno_pos,
    hmm_CN_multi           = fixtures$multi_res,
    ploidy.standardization = NULL,   # should be inferred
    n.cores                = 1,
    verbose                = TRUE    # prints the inferred value
  ))

  expect_s3_class(result, "qploidy_standardization")
  # inferred ploidy should be a positive integer stored in info
  inferred <- as.numeric(result$info["ploidy.standardization"])
  expect_true(is.finite(inferred) && inferred >= 1)
})

test_that("re_standardize writes output file when out_filename is provided", {
  fixtures <- local_fixtures()
  tmp <- tempfile(fileext = ".tsv.gz")

  re_standardize(
    data                   = fixtures$data,
    geno.pos               = fixtures$geno_pos,
    hmm_CN_multi           = fixtures$multi_res,
    ploidy.standardization = 4,
    n.cores                = 1,
    out_filename           = tmp,
    verbose                = FALSE
  )

  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 0)
})

test_that("re_standardize print method works without error", {
  fixtures <- local_fixtures()

  result <- re_standardize(
    data                   = fixtures$data,
    geno.pos               = fixtures$geno_pos,
    hmm_CN_multi           = fixtures$multi_res,
    ploidy.standardization = 4,
    n.cores                = 1,
    verbose                = FALSE
  )

  printed <- expect_no_error(print(result))
  expect_true(is.data.frame(printed) || is.null(printed))
})

test_that("re_standardize with selected_model = NULL errors from call_hmm_dosages", {
  fixtures <- local_fixtures()

  # call_hmm_dosages requires either a selected_BAF_model or explicit bw/dist;
  # passing NULL without those args must error informatively.
  expect_error(
    re_standardize(
      data                   = fixtures$data,
      geno.pos               = fixtures$geno_pos,
      hmm_CN_multi           = fixtures$multi_res,
      selected_model         = NULL,
      ploidy.standardization = 4,
      n.cores                = 1,
      verbose                = FALSE
    )
  )
})
