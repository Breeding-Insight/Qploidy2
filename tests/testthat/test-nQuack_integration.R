# Test suite for nQuack integration functions
library(testthat)
library(Qploidy2)

# ---------------------------------------------------------------------------
# Shared simulation setup
# ---------------------------------------------------------------------------
# NOTE: genos is filtered to Tetraploid samples only (matching test-hmm_cn.R).
# This means standardize() assigns geno calls only for tetraploid-consistent
# sites. Diploid and triploid samples end up with geno == 0 or geno == max_geno
# (homozygous by genotype call) at virtually all sites, so convert2nQuack()
# correctly returns empty matrices for those samples — they have no usable
# heterozygous sites after the homozygous filters are applied.
setup_simu <- function(seed = 123) {
  set.seed(seed)
  vcf_file <- tempfile(fileext = ".vcf")
  simulate_vcf(
    seed = seed, file_path = vcf_file,
    n_tetraploid = 5, n_diploid = 2, n_triploid = 2,
    n_markers = 50
  )
  data      <- qploidy_read_vcf(vcf_file)
  genos     <- qploidy_read_vcf(vcf_file, geno = TRUE)
  genos     <- genos[grep("Tetraploid", genos$SampleName), ]
  genos.pos <- qploidy_read_vcf(vcf_file, geno.pos = TRUE)

  stand_file <- tempfile(fileext = ".tsv.gz")
  simu_data_standardized <- standardize(
    data = data,
    genos = genos,
    geno.pos = genos.pos,
    ploidy.standardization = 4,
    threshold.n.clusters = 3,
    n.cores = 1,
    out_filename = stand_file,
    verbose = FALSE
  )
  list(std = simu_data_standardized, vcf = vcf_file)
}

# ---------------------------------------------------------------------------
# convert2nQuack tests
# ---------------------------------------------------------------------------
test_that("convert2nQuack returns a named list with one matrix per sample", {
  env <- setup_simu()
  std <- env$std

  result <- convert2nQuack(std)

  expect_type(result, "list")
  # One entry per sample in the standardized data
  expected_samples <- unique(std$data$SampleName)
  expect_equal(sort(names(result)), sort(expected_samples))

  # Every element must be a two-column numeric matrix
  for (nm in names(result)) {
    mat <- result[[nm]]
    expect_true(is.matrix(mat), info = paste("Element", nm, "should be a matrix"))
    expect_equal(ncol(mat), 2L, info = paste("Element", nm, "should have 2 columns"))
    expect_true(is.numeric(mat), info = paste("Element", nm, "should be numeric"))
  }
})

test_that("convert2nQuack: tetraploids have data; diploids/triploids use ratio fallback", {
  env <- setup_simu()
  std <- env$std

  result <- convert2nQuack(std)

  # Tetraploids were used as reference genotypes → they have het geno calls →
  # convert2nQuack should return non-empty matrices for them.
  tetra_names <- grep("Tetraploid", names(result), value = TRUE)
  for (nm in tetra_names) {
    expect_gt(nrow(result[[nm]]), 0L,
              label = paste(nm, "should have rows after filtering"))
  }

  # Diploids and triploids have geno == NA (not used as reference genotypes).
  # They fall back to the ratio-based het filter [0.1, 0.9].
  # With het_range = c(0.1, 0.9) (default) sites with ratio in that interval
  # are kept, so these samples may have non-empty matrices.
  non_tetra_names <- grep("Diploid|Triploid", names(result), value = TRUE)
  for (nm in non_tetra_names) {
    mat <- result[[nm]]
    expect_true(is.matrix(mat), info = paste(nm, "should still be a matrix"))
    expect_equal(ncol(mat), 2L, info = paste(nm, "should have 2 columns"))
  }
})

test_that("convert2nQuack: het_range controls fallback filter for geno=NA samples", {
  env <- setup_simu()
  std <- env$std

  # Wide range: more sites pass for geno=NA samples
  result_wide   <- convert2nQuack(std, het_range = c(0.1, 0.9))
  # Narrow range: fewer sites pass
  result_narrow <- convert2nQuack(std, het_range = c(0.4, 0.6))
  # Impossible range: no sites pass → all geno=NA samples are empty
  result_empty  <- convert2nQuack(std, het_range = c(0.6, 0.4))

  non_tetra_names <- grep("Diploid|Triploid", names(result_wide), value = TRUE)
  for (nm in non_tetra_names) {
    expect_true(nrow(result_wide[[nm]])   >= nrow(result_narrow[[nm]]),
                label = paste(nm, ": wide >= narrow row count"))
    expect_equal(nrow(result_empty[[nm]]), 0L,
                 label = paste(nm, ": impossible range yields empty matrix"))
  }
})

test_that("convert2nQuack column names are R and allele", {
  env    <- setup_simu()
  result <- convert2nQuack(env$std)

  for (nm in names(result)) {
    expect_equal(colnames(result[[nm]]), c("R", "allele"),
                 info = paste("Column names check for sample", nm))
  }
})

test_that("convert2nQuack: from_standardized_baf = FALSE uses raw X/Y columns", {
  env <- setup_simu()
  std <- env$std

  result_baf <- convert2nQuack(std, from_standardized_baf = TRUE)
  result_raw <- convert2nQuack(std, from_standardized_baf = FALSE)

  expect_type(result_raw, "list")
  # Same samples returned in both cases
  expect_equal(sort(names(result_baf)), sort(names(result_raw)))

  for (nm in names(result_raw)) {
    mat <- result_raw[[nm]]
    expect_true(is.matrix(mat))
    expect_equal(ncol(mat), 2L)
  }
})

test_that("convert2nQuack respects min_depth filter", {
  env <- setup_simu()
  std <- env$std

  result_low  <- convert2nQuack(std, min_depth = 1L)
  result_high <- convert2nQuack(std, min_depth = 1000L)

  # Higher depth threshold → fewer or equal rows per sample
  for (nm in names(result_low)) {
    expect_true(nrow(result_low[[nm]]) >= nrow(result_high[[nm]]),
                label = paste(nm, ": low min_depth >= high min_depth row count"))
  }
})

test_that("convert2nQuack accepts a plain data.frame instead of qploidy_standardization", {
  env <- setup_simu()
  df  <- env$std$data  # plain data.frame

  result <- convert2nQuack(df, max_geno = 4L)
  expect_type(result, "list")
  expect_equal(sort(names(result)), sort(unique(df$SampleName)))
})

test_that("convert2nQuack errors on missing required columns", {
  env        <- setup_simu()
  df_missing <- env$std$data[, setdiff(names(env$std$data), "R")]

  expect_error(convert2nQuack(df_missing, max_geno = 4L),
               regexp = "missing required columns")
})

test_that("convert2nQuack is reproducible with the same seed", {
  env <- setup_simu()
  r1  <- convert2nQuack(env$std, seed = 42L)
  r2  <- convert2nQuack(env$std, seed = 42L)

  for (nm in names(r1)) {
    expect_equal(r1[[nm]], r2[[nm]],
                 info = paste("Reproducibility check for sample", nm))
  }
})

# ---------------------------------------------------------------------------
# hmm_estimate_CN_nQuack
# ---------------------------------------------------------------------------
test_that("hmm_estimate_CN_nQuack returns valid hmm_CN object", {
  skip_if_not_installed("nQuack")

  env    <- setup_simu()
  std    <- env$std
  # Use a tetraploid sample — it has het sites that nQuack can model
  sample <- grep("Tetraploid", unique(std$data$SampleName), value = TRUE)[1]

  res <- hmm_estimate_CN_nQuack(
    qploidy_standarize_result = std,
    sample_id    = sample,
    chr          = 1,
    cn_grid      = c(2, 3, 4, 5, 6),
    verbose      = TRUE
  )

  expect_s3_class(res, "hmm_CN")
  expect_true("by_window" %in% names(res))
  expect_true("by_marker" %in% names(res))
  expect_true("params"    %in% names(res))
  expect_true(is.data.frame(res$by_window))
  expect_true(nrow(res$by_window) > 0)
  expect_true(all(c("Sample", "Chr", "WindowID", "CN_call") %in% names(res$by_window)))
  # Distribution should be recorded in params as a non-empty string
  expect_true(!is.null(res$params$distribution) && nchar(res$params$distribution) > 0)
  # bw is NA — nQuack builds synthetic selected_model, no kernel-density template
  expect_true(is.na(res$params$bw))
  expect_equal(as.numeric(summarize_cn_mode(res)[2]), 4)
})

test_that("hmm_estimate_CN_nQuack with restricted nQuack_dists", {
  skip_if_not_installed("nQuack")

  env    <- setup_simu()
  std    <- env$std
  sample <- grep("Tetraploid", unique(std$data$SampleName), value = TRUE)[1]

  res <- hmm_estimate_CN_nQuack(
    qploidy_standarize_result = std,
    sample_id    = sample,
    chr          = 1,
    cn_grid      = c(2, 3, 4),
    nQuack_dists = "normal",
    verbose      = TRUE
  )

  expect_s3_class(res, "hmm_CN")
  # Winning distribution should be Normal or Normal+Uniform
  expect_true(grepl("normal", res$params$distribution, ignore.case = TRUE))
  expect_equal(as.numeric(summarize_cn_mode(res)[2]), 4)

  sample <- grep("Diploid", unique(std$data$SampleName), value = TRUE)[1]

  res <- hmm_estimate_CN_nQuack(
    qploidy_standarize_result = std,
    sample_id    = sample,
    chr          = 1,
    cn_grid      = c(2, 3, 4),
    nQuack_dists = "normal", nQuack_uniform = 1,
    verbose      = TRUE, dosage_threshold = 0.6
  )

  expect_s3_class(res, "hmm_CN")
  # Winning distribution should be Normal or Normal+Uniform
  expect_true(grepl("normal", res$params$distribution, ignore.case = TRUE))
  expect_equal(as.numeric(summarize_cn_mode(res)[2]), 2)

  sample <- grep("Triploid", unique(std$data$SampleName), value = TRUE)[1]

  res <- hmm_estimate_CN_nQuack(
    qploidy_standarize_result = std,
    sample_id    = sample,
    chr          = 1,
    nQuack_uniform = 1,
    cn_grid      = c(2, 3, 4),
    nQuack_dists = "normal",
    verbose      = TRUE
  )

  expect_s3_class(res, "hmm_CN")
  # Winning distribution should be Normal or Normal+Uniform
  expect_true(grepl("normal", res$params$distribution, ignore.case = TRUE))
  expect_equal(as.numeric(summarize_cn_mode(res)[2]), 3)
})

# ---------------------------------------------------------------------------
# hmm_estimate_CN_multi with nQuack options
# ---------------------------------------------------------------------------
test_that("hmm_estimate_CN_multi with nQuack runs without error", {
  skip_if_not_installed("nQuack")

  env <- setup_simu()
  std <- env$std

  multi_res <- hmm_estimate_CN_multi(
    qploidy_standarize_result = std,
    n_cores      = 1,
    chr          = 1,
    cn_grid      = c(2, 3, 4),
    nQuack_dists = "normal",
    nQuack_uniform = 1
  )

  p <- plot_cn_track(multi_res, sample_id = "Diploid1")
  p$arranged

  expect_s3_class(multi_res, "hmm_CN")
  expect_true(is.data.frame(multi_res$by_window))
  expect_true(nrow(multi_res$by_window) > 0)
  expect_true(all(c("Sample", "Chr", "WindowID", "CN_call") %in% names(multi_res$by_window)))
  expect_equal(summarize_cn_mode(multi_res)[,2], c(2,2,4,4,4,4,4,3,3))

})

# ---------------------------------------------------------------------------
# hmm_estimate_CN_nQuack: argument variants, pinned to the simulated ploidy
# ---------------------------------------------------------------------------
# Simulated truth for setup_simu(): DiploidN = 2, TriploidN = 3, TetraploidN = 4.
true_cn <- c(Diploid1 = 2L, Triploid1 = 3L, Tetraploid1 = 4L)

run_nq <- function(std, sample, ...) {
  set.seed(42)
  suppressWarnings(suppressMessages(
    hmm_estimate_CN_nQuack(qploidy_standarize_result = std, sample_id = sample,
                           verbose = FALSE, ...)
  ))
}

test_that("hmm_estimate_CN_nQuack input validation", {
  skip_if_not_installed("nQuack")
  std <- setup_simu()$std

  expect_error(run_nq(std, "Tetraploid1", chr = 1, cn_grid = 2:8),
               "subset of 0:6")
  expect_error(run_nq(std, "Tetraploid1", chr = 1, cn_grid = 2:4,
                      nQuack_dists = c("normal", "beta")),
               "single distribution family")
  expect_error(run_nq(std, "Tetraploid1", chr = 1, cn_grid = 2:4, nQuack_dists = NULL),
               "single distribution family")
  # Invalid type/uniform make the global fit fail; without exp_ploidy this must stop
  expect_error(run_nq(std, "Tetraploid1", chr = 1, cn_grid = 2:4, nQuack_type = "free"),
               "exp_ploidy")
  expect_error(run_nq(std, "Tetraploid1", chr = 1, cn_grid = 2:4, nQuack_uniform = 2),
               "exp_ploidy")
})

test_that("hmm_estimate_CN_nQuack (normal, fixed, uniform = 1) recovers simulated ploidy", {
  skip_if_not_installed("nQuack")
  std <- setup_simu()$std

  for (s in names(true_cn)) {
    res <- run_nq(std, s, chr = 1, cn_grid = 2:4,
                  nQuack_dists = "normal", nQuack_type = "fixed", nQuack_uniform = 1)
    expect_true(all(res$by_window$CN_call == true_cn[[s]]), info = s)
    expect_equal(res$params$exp_ploidy, true_cn[[s]], info = s)
    expect_equal(res$params$distribution, "normal-uniform", info = s)
    expect_equal(res$params$nQuack_type, "fixed", info = s)
    expect_equal(res$params$nQuack_uniform, 1, info = s)
  }
})

test_that("hmm_estimate_CN_nQuack baf_only = TRUE recovers simulated ploidy and ignores z", {
  skip_if_not_installed("nQuack")
  std <- setup_simu()$std

  for (s in names(true_cn)) {
    res <- run_nq(std, s, chr = 1, cn_grid = 2:4, baf_only = TRUE)
    expect_equal(res$params$exp_ploidy, true_cn[[s]], info = s)
    # Diploid BAF peaks (0, 0.5, 1) are a subset of CN4's, so short z-derived windows can call CN4
    if (s != "Diploid1")
      expect_true(all(res$by_window$CN_call == true_cn[[s]]), info = s)
    expect_true(all(is.na(res$params$mu)), info = s)
    # z-based collapsed-state pruning must not remove BAF-supported states
    expect_equal(res$params$cn_grid, 2:4, info = s)
  }
})

test_that("hmm_estimate_CN_nQuack use_values = c('ratio', 'zscore') uses ratio for the global fit", {
  skip_if_not_installed("nQuack")
  std <- setup_simu()$std

  for (s in c("Diploid1", "Triploid1")) {
    res <- run_nq(std, s, chr = 1, cn_grid = 2:4, use_values = c("ratio", "zscore"))
    expect_true(all(res$by_window$CN_call == true_cn[[s]]), info = s)
    expect_equal(res$params$exp_ploidy, true_cn[[s]], info = s)
  }
})

test_that("hmm_estimate_CN_nQuack nQuack_type = 'fixed_2' recovers simulated ploidy", {
  skip_if_not_installed("nQuack")
  std <- setup_simu()$std

  for (s in c("Triploid1", "Tetraploid1")) {
    res <- run_nq(std, s, chr = 1, cn_grid = 2:4, nQuack_type = "fixed_2")
    expect_true(all(res$by_window$CN_call == true_cn[[s]]), info = s)
    expect_equal(res$params$nQuack_type, "fixed_2", info = s)
  }
})

test_that("hmm_estimate_CN_nQuack gives the same calls for a reordered cn_grid", {
  skip_if_not_installed("nQuack")
  std <- setup_simu()$std

  for (s in c("Diploid1", "Triploid1")) {
    res_sorted    <- run_nq(std, s, chr = 1, cn_grid = 2:4)
    res_reordered <- run_nq(std, s, chr = 1, cn_grid = c(4, 2, 3))
    expect_true(all(res_reordered$by_window$CN_call == true_cn[[s]]), info = s)
    expect_equal(res_reordered$by_window$CN_call, res_sorted$by_window$CN_call, info = s)
  }
})

test_that("hmm_estimate_CN_nQuack falls back to the requested distribution when the global fit fails", {
  skip_if_not_installed("nQuack")
  std <- setup_simu()$std

  # het_range this narrow leaves < 10 het sites, so the global nQuack fit fails
  res <- run_nq(std, "Tetraploid1", chr = 1, cn_grid = 2:4,
                het_range = c(0.49, 0.491), exp_ploidy = 4)
  expect_true(all(res$by_window$CN_call == 4L))
  expect_equal(res$params$distribution, "normal-uniform")

  expect_error(run_nq(std, "Tetraploid1", chr = 1, cn_grid = 2:4,
                      het_range = c(0.49, 0.491)),
               "exp_ploidy")
})

test_that("hmm_estimate_CN_nQuack keeps windows below min_snps_per_window", {
  skip_if_not_installed("nQuack")
  std <- setup_simu()$std

  for (s in c("Triploid1", "Tetraploid1")) {
    res <- run_nq(std, s, cn_grid = 2:4, segment_zscore = FALSE,
                  snps_per_window = 9, min_snps_per_window = 10)
    expect_true(any(res$by_window$n_snps < 10), info = s)
    expect_equal(sum(res$by_window$n_snps), nrow(res$by_marker), info = s)
    expect_equal(res$params$exp_ploidy, true_cn[[s]], info = s)
  }
})

test_that("hmm_estimate_CN_nQuack supports CN = 1 in cn_grid", {
  skip_if_not_installed("nQuack")
  std <- setup_simu()$std

  # Heterozygous samples: CN1 is ruled out by min_het_frac and calls are unchanged
  for (s in c("Diploid1", "Tetraploid1")) {
    res <- run_nq(std, s, chr = 1, cn_grid = 1:4)
    expect_equal(res$params$exp_ploidy, true_cn[[s]], info = s)
    expect_false(any(res$by_window$CN_call == 1L, na.rm = TRUE), info = s)
  }

  # Haploid-like sample (no heterozygous BAF): exp_ploidy falls back to 1
  hap <- std
  idx <- hap$data$SampleName == "Diploid1"
  hap$data$baf[idx] <- round(hap$data$baf[idx])
  res <- run_nq(hap, "Diploid1", chr = 1, cn_grid = 1:4)
  expect_equal(res$params$exp_ploidy, 1L)
  expect_true(all(!is.na(res$by_window$CN_call)))
})

test_that("hmm_estimate_CN_nQuack supports CN = 0 via a z-only pre-pass", {
  skip_if_not_installed("nQuack")
  std <- setup_simu()$std
  chrs <- sort(unique(as.character(std$data$Chr)))

  # Collapse depth on the second chromosome to mimic nullisomy
  idx <- std$data$SampleName == "Tetraploid1" & std$data$Chr == chrs[2]
  std$data$z[idx] <- std$data$z[idx] - 20

  res <- run_nq(std, "Tetraploid1", cn_grid = 0:4)
  expect_true(all(is.na(res$by_window$CN_call[res$by_window$Chr == chrs[2]])))
  expect_true(all(res$by_window$CN_call[res$by_window$Chr == chrs[1]] == 4L))
  expect_true(chrs[2] %in% res$by_window$Chr)
})

test_that(".nquack_estimate_ploidy reports the correct ploidy for CN 2-6 (hexa/penta labels)", {
  skip_if_not_installed("nQuack")

  # High-depth BAF drawn from dosages 1..(cn-1) of a CN-copy genome
  sim_d <- function(cn, seed, n = 1500) {
    set.seed(seed)
    dos <- sample(seq_len(cn - 1), n, replace = TRUE)
    R   <- sample(800:1000, n, replace = TRUE)
    data.frame(baf = stats::rbinom(n, R, dos / cn) / R, R = R)
  }

  for (cn in 2:6) {
    for (u in c(0, 1)) {
      set.seed(1)
      est <- .nquack_estimate_ploidy(sim_d(cn, cn), "s", cn_grid = 2:6,
                                     uniform = u, verbose = FALSE)
      expect_equal(est$cn, cn, info = sprintf("cn = %d, uniform = %d", cn, u))
    }
  }

  # Subsets / reordered grids must not swap labels
  for (g in list(c(5, 6), c(6, 5, 4))) {
    set.seed(1)
    expect_equal(.nquack_estimate_ploidy(sim_d(5, 5), "s", cn_grid = g, uniform = 0,
                                         verbose = FALSE)$cn, 5L)
    set.seed(1)
    expect_equal(.nquack_estimate_ploidy(sim_d(6, 6), "s", cn_grid = g, uniform = 0,
                                         verbose = FALSE)$cn, 6L)
  }
})
