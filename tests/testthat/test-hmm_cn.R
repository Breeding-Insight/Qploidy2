# Test suite for hmm_estimate_CN and plot_cn_track
library(testthat)
library(Qploidy2)

test_that("hmm_estimate_CN, plot_cn_track, and other HMM functions work as expected", {
  set.seed(123)
  vcf_file <- tempfile(fileext = ".vcf")
  simulate_vcf(
    seed = 123, file_path = vcf_file,
    n_tetraploid = 5, n_diploid = 2, n_triploid = 2,
    n_markers = 50
  )
  data <- qploidy_read_vcf(vcf_file)
  genos <- qploidy_read_vcf(vcf_file, geno = TRUE)
  genos <- genos[grep("Tetraploid", genos$SampleName), ]
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
  sample <- unique(simu_data_standardized$data$SampleName)[1]

  # hmm_estimate_CN test
  res <- hmm_estimate_CN(
    qploidy_standarize_result = simu_data_standardized,
    sample_id = sample,
    chr = 1,
    segment_zscore = TRUE,
    snps_per_window = 10,
    min_snps_per_window = 5,
    cn_grid = c(2, 3, 4),
    M = 21,
    exp_ploidy = 4
  )
  expect_type(res, "list")
  expect_true("by_window" %in% names(res))
  expect_true("params" %in% names(res))
  expect_s3_class(res, "hmm_CN")
  expect_true(is.data.frame(res$by_window))
  expect_true(nrow(res$by_window) > 0)
  expect_true(all(c("Sample", "Chr", "WindowID", "CN_call") %in% names(res$by_window)))
  expect_equal(round(sum(res$by_window$post_CN4),2), 6)
  expect_equal(round(sum(res$by_marker$post_max),2), 50)

  res <- hmm_estimate_CN(
    qploidy_standarize_result = simu_data_standardized,
    sample_id = sample,
    chr = 1,
    segment_zscore = TRUE,
    snps_per_window = 10,
    min_snps_per_window = 5,
    cn_grid = c(2, 3, 4),
    M = 21,
    exp_ploidy = NA
  )

  expect_type(res, "list")
  expect_true("by_window" %in% names(res))
  expect_true("params" %in% names(res))
  expect_s3_class(res, "hmm_CN")
  expect_true(is.data.frame(res$by_window))
  expect_true(nrow(res$by_window) > 0)
  expect_true(all(c("Sample", "Chr", "WindowID", "CN_call") %in% names(res$by_window)))
  expect_equal(round(sum(res$by_window$post_CN4),2), 6)
  expect_equal(round(sum(res$by_marker$post_max),2), 50)

  res <- hmm_estimate_CN(
    qploidy_standarize_result = simu_data_standardized,
    sample_id = sample,
    chr = 1,
    segment_zscore = FALSE,
    snps_per_window = 10,
    min_snps_per_window = 5,
    cn_grid = c(2, 3, 4),
    M = 21,
    exp_ploidy = NA,
    rm_outliers = FALSE,
  )

  expect_type(res, "list")
  expect_true("by_marker" %in% names(res))
  expect_true("by_window" %in% names(res))
  expect_true("params" %in% names(res))
  expect_s3_class(res, "hmm_CN")
  expect_true(is.data.frame(res$by_window))
  expect_true(nrow(res$by_window) > 0)
  expect_true(all(c("Sample", "Chr", "WindowID", "CN_call") %in% names(res$by_window)))
  expect_equal(round(sum(res$by_window$post_CN4),2), 5)
  expect_equal(round(sum(res$by_marker$post_max),2), 50)

  samples <- unique(simu_data_standardized$data$SampleName)
  res <- hmm_estimate_CN(
    qploidy_standarize_result = simu_data_standardized,
    sample_id = samples[1],
    chr = 1,
    cn_grid = c(1:4)
  )

  expect_type(res, "list")
  expect_true("by_marker" %in% names(res))
  expect_true("by_window" %in% names(res))
  expect_true("params" %in% names(res))
  expect_s3_class(res, "hmm_CN")
  expect_true(is.data.frame(res$by_window))
  expect_true(nrow(res$by_window) > 0)
  expect_true(all(c("Sample", "Chr", "WindowID", "CN_call") %in% names(res$by_window)))
  expect_equal(round(sum(res$by_window$post_CN4),2), 4)
  expect_equal(round(sum(res$by_marker$post_max),2), 50.000)

  # plot_cn_track test (should return a gg object)
  p <- plot_cn_track(hmm_CN = res,
                     qploidy_standarize_result = simu_data_standardized,
                     sample_id = samples[1],
                     show_window_lines = TRUE)
  expect_type(p, "list")
  expect_true("gg" %in% class(p$arranged) || "ggarrange" %in% class(p$arranged))

  # Test hmm_estimate_CN_multi
  multi_res <- hmm_estimate_CN_multi(
    qploidy_standarize_result = simu_data_standardized,
    n_cores = 1,
    chr = 1,
    bw = c(0.02,0.03,0.04, 0.06),
    cn_grid = c(2, 3, 4),
  )

  expect_type(multi_res, "list")
  expect_true(length(multi_res) >= 1)
  expect_true(inherits(multi_res, "hmm_CN"))
  expect_equal(mean(multi_res$by_window$CN_call), 3.14, tolerance = 1e-3)
  expect_equal(mean(multi_res$by_marker$CN_call), 3.333, tolerance = 1e-3)

  # Test summarize_cn_mode
  summ <- summarize_cn_mode(df = multi_res, level = "sample")
  expect_true(is.data.frame(summ))
  expect_true("CN_mode" %in% names(summ))
  expect_equal(summ$CN_mode, c(2,2,4,4,4,4,4,3,3))

  multi_res <- hmm_estimate_CN_multi(
    qploidy_standarize_result = simu_data_standardized,
    sample_ids = "all",
    n_cores = 1,
    chr = 1,
    bw = c(0.02,0.03,0.04, 0.06),
    cn_grid = c(1:4)
  )

  expect_type(multi_res, "list")
  expect_true(length(multi_res) >= 1)
  expect_true(inherits(multi_res, "hmm_CN"))
  expect_equal(mean(multi_res$by_window$CN_call), 3.14, tolerance = 1e-3)
  expect_equal(mean(multi_res$by_marker$CN_call), 3.333, tolerance = 1e-3)

  oneCN <- hmm_estimate_CN(
    qploidy_standarize_result = simu_data_standardized,
    sample_id = "Triploid1",
    chr = 1,
    segment_zscore = TRUE,
    snps_per_window =10,
    min_snps_per_window = 5,
    cn_grid = c(2, 3, 4),
    exp_ploidy = 3
  )

  updated_CN <- update_hmm_CN_multi(multi_res, oneCN)

  temp_dir <- tempdir()
  formatted_date <- format(Sys.Date(), "%m_%d_%Y")

  write_hmm_CN(updated_CN, prefix = paste0(temp_dir, "/", formatted_date))

  # Test summarize_cn_mode
  summ <- summarize_cn_mode(df = multi_res, level = "sample")
  expect_true(is.data.frame(summ))
  expect_true("CN_mode" %in% names(summ))
  expect_equal(summ$CN_mode, c(2,2,4,4,4,4,4,3,3))

  summ <- summarize_cn_mode(multi_res, level = "chromosome")
  expect_true(is.data.frame(summ))
  expect_true("CN_mode" %in% names(summ))
  expect_equal(summ$CN_mode, c(2,2,4,4,4,4,4,3,3))

  p <- compare_cn_track(multi_res, samples_to_plot = samples)

  expect_true("gg" %in% class(p))

})

test_that("hmm_estimate_CN additional argument coverage", {
  set.seed(123)
  vcf_file <- tempfile(fileext = ".vcf")
  simulate_vcf(
    seed = 123, file_path = vcf_file,
    n_tetraploid = 5, n_diploid = 2, n_triploid = 2,
    n_markers = 50
  )
  data    <- qploidy_read_vcf(vcf_file)
  genos   <- qploidy_read_vcf(vcf_file, geno = TRUE)
  genos   <- genos[grep("Tetraploid", genos$SampleName), ]
  genos.pos <- qploidy_read_vcf(vcf_file, geno.pos = TRUE)

  simu_std <- standardize(
    data = data, genos = genos, geno.pos = genos.pos,
    ploidy.standardization = 4, threshold.n.clusters = 3,
    n.cores = 1, verbose = FALSE
  )
  sample <- unique(simu_std$data$SampleName)[1]

  # --- z_only = TRUE: BAF emission is disabled; w_baf must be 0 everywhere ---
  res_zo <- hmm_estimate_CN(
    qploidy_standarize_result = simu_std,
    sample_id = sample, chr = 1,
    cn_grid = c(2, 3, 4), M = 21, exp_ploidy = 4,
    z_only = TRUE, verbose = FALSE
  )
  expect_s3_class(res_zo, "hmm_CN")
  expect_true(is.data.frame(res_zo$by_window))
  expect_true(nrow(res_zo$by_window) > 0)
  expect_true(all(res_zo$by_window$w_baf == 0))
  expect_true(res_zo$params$z_only)

  # --- rerun_overall_ploidy = TRUE: second pass with cleaner BAF model ---
  res_rr <- hmm_estimate_CN(
    qploidy_standarize_result = simu_std,
    sample_id = sample, chr = 1,
    cn_grid = c(2, 3, 4), M = 21, exp_ploidy = 4,
    rerun_overall_ploidy = TRUE, verbose = FALSE
  )
  expect_s3_class(res_rr, "hmm_CN")
  expect_true(is.data.frame(res_rr$by_window))
  expect_true(nrow(res_rr$by_window) > 0)
  expect_true(all(c("Sample", "Chr", "WindowID", "CN_call") %in% names(res_rr$by_window)))

  # --- selected_model: pre-computed BAF model is passed in directly ---
  baf_vec <- simu_std$data$baf[simu_std$data$SampleName == sample]
  sel_mod <- select_best_baf_model(baf_vec = baf_vec, cn_grid = c(2, 3, 4), M = 21)
  res_sm <- hmm_estimate_CN(
    qploidy_standarize_result = simu_std,
    sample_id = sample, chr = 1,
    cn_grid = c(2, 3, 4), M = 21, exp_ploidy = 4,
    selected_model = sel_mod, verbose = FALSE
  )
  expect_s3_class(res_sm, "hmm_CN")
  expect_equal(res_sm$params$bw,           sel_mod$best$bw)
  expect_equal(res_sm$params$distribution, sel_mod$best$dist)

  # --- use_values = c("BAF", "R"): R used as depth metric instead of zscore ---
  res_br <- hmm_estimate_CN(
    qploidy_standarize_result = simu_std,
    sample_id = sample, chr = 1,
    use_values = c("BAF", "R"),
    cn_grid = c(2, 3, 4), M = 21, exp_ploidy = 4,
    verbose = FALSE
  )
  expect_s3_class(res_br, "hmm_CN")
  expect_true(is.data.frame(res_br$by_window))
  expect_true(nrow(res_br$by_window) > 0)

  # --- transition_jump and initial_prob are stored in params ---
  res_tp <- hmm_estimate_CN(
    qploidy_standarize_result = simu_std,
    sample_id = sample, chr = 1,
    cn_grid = c(2, 3, 4), M = 21, exp_ploidy = 4,
    transition_jump = 0.99, initial_prob = 0.5,
    verbose = FALSE
  )
  expect_s3_class(res_tp, "hmm_CN")
  expect_equal(res_tp$params$transition_jump, 0.99)

  # --- input validation: non-existent sample_id must error ---
  expect_error(
    hmm_estimate_CN(
      qploidy_standarize_result = simu_std,
      sample_id = "DoesNotExist",
      chr = 1, cn_grid = c(2, 3, 4)
    )
  )

  # --- input validation: selected_model wrong class must error ---
  expect_error(
    hmm_estimate_CN(
      qploidy_standarize_result = simu_std,
      sample_id = sample, chr = 1,
      cn_grid = c(2, 3, 4),
      selected_model = list(not_a_model = TRUE)
    )
  )
})

# ---------------------------------------------------------------------------
# hmm_estimate_CN: argument variants, pinned to the simulated ploidy
# ---------------------------------------------------------------------------
test_that("hmm_estimate_CN baf_only, reordered cn_grid and window filter recover simulated ploidy", {
  set.seed(123)
  vcf_file <- tempfile(fileext = ".vcf")
  simulate_vcf(seed = 123, file_path = vcf_file,
               n_tetraploid = 5, n_diploid = 2, n_triploid = 2, n_markers = 50)
  data      <- qploidy_read_vcf(vcf_file)
  genos     <- qploidy_read_vcf(vcf_file, geno = TRUE)
  genos     <- genos[grep("Tetraploid", genos$SampleName), ]
  genos.pos <- qploidy_read_vcf(vcf_file, geno.pos = TRUE)
  simu_std  <- standardize(data = data, genos = genos, geno.pos = genos.pos,
                           ploidy.standardization = 4, threshold.n.clusters = 3,
                           n.cores = 1, out_filename = tempfile(fileext = ".tsv.gz"),
                           verbose = FALSE)

  # Simulated truth. Diploids are not pinned here: the grid route standardized
  # on tetraploid genotypes currently calls them CN = 4.
  true_cn <- c(Triploid1 = 3L, Tetraploid1 = 4L)

  run_grid <- function(sample, ...) {
    set.seed(42)
    suppressWarnings(suppressMessages(
      hmm_estimate_CN(qploidy_standarize_result = simu_std, sample_id = sample,
                      verbose = FALSE, ...)
    ))
  }

  for (s in names(true_cn)) {
    # baf_only: BAF-only emission and prior, z-based pruning skipped
    res <- run_grid(s, chr = 1, cn_grid = 2:4, baf_only = TRUE)
    expect_true(all(res$by_window$CN_call == true_cn[[s]]), info = paste("baf_only", s))
    expect_true(all(is.na(res$params$mu)), info = paste("baf_only", s))

    # Reordered cn_grid gives the same calls as the sorted one
    res_sorted    <- run_grid(s, chr = 1, cn_grid = 2:4)
    res_reordered <- run_grid(s, chr = 1, cn_grid = c(4, 2, 3))
    expect_true(all(res_sorted$by_window$CN_call == true_cn[[s]]), info = paste("sorted", s))
    expect_equal(res_reordered$by_window$CN_call, res_sorted$by_window$CN_call,
                 info = paste("reordered", s))

    # Windows smaller than min_snps_per_window are kept
    res <- run_grid(s, cn_grid = 2:4, segment_zscore = FALSE,
                    snps_per_window = 9, min_snps_per_window = 10)
    expect_true(any(res$by_window$n_snps < 10), info = paste("min_snps", s))
    expect_equal(sum(res$by_window$n_snps), nrow(res$by_marker), info = paste("min_snps", s))
  }

  # Per-chromosome min_snps_per_window is stored and no chromosome is lost
  chrs <- sort(unique(as.character(simu_std$data$Chr)))
  res_chr <- run_grid("Tetraploid1", cn_grid = 2:4,
                      min_snps_per_window = stats::setNames(c(40, 5), chrs[1:2]))
  expect_setequal(unique(res_chr$by_window$Chr), chrs[1:2])
  expect_equal(res_chr$params$min_snps_per_window, stats::setNames(c(40L, 5L), chrs[1:2]))

  # max_cn_dev keeps only CN states near exp_ploidy (CN 0 always kept)
  res_dev <- run_grid("Tetraploid1", chr = 1, cn_grid = 0:6, max_cn_dev = 1)
  expect_equal(res_dev$params$max_cn_dev, 1)
  expect_true(all(res_dev$params$cn_grid %in% c(0, 3:5)))
  expect_true(all(res_dev$by_window$CN_call %in% c(3:5, NA)))
})

test_that(".resolve_min_snps handles NULL, scalar and named inputs", {
  d <- data.frame(Chr = rep(c("A", "B", "X"), c(100, 60, 20)))

  expect_equal(.resolve_min_snps(NULL, d, verbose = FALSE), c(A = 5L, B = 5L, X = 5L))
  expect_equal(.resolve_min_snps(30, d, verbose = FALSE), c(A = 30L, B = 30L, X = 30L))
  expect_equal(.resolve_min_snps(c(A = 30, X = 3), d, verbose = FALSE),
               c(A = 30L, B = 5L, X = 3L))
  expect_warning(.resolve_min_snps(c(Y = 3), d, verbose = FALSE), "not found")
  expect_silent(.resolve_min_snps(c(Y = 3), d, chr_subset = "A", verbose = FALSE))
  expect_error(.resolve_min_snps(c(10, 20), d, verbose = FALSE), "single value")
  expect_error(.resolve_min_snps(0, d, verbose = FALSE), "positive")
})

test_that("hmm_estimate_CN scores Y presence from raw depth before filtering", {
  vcf_file <- tempfile(fileext = ".vcf")
  simulate_vcf(seed = 123, file_path = vcf_file,
               n_tetraploid = 5, n_diploid = 2, n_triploid = 2, n_markers = 50)
  data  <- qploidy_read_vcf(vcf_file)
  genos <- qploidy_read_vcf(vcf_file, geno = TRUE)
  genos <- genos[grep("Tetraploid", genos$SampleName), ]
  gp    <- qploidy_read_vcf(vcf_file, geno.pos = TRUE)
  st <- suppressWarnings(standardize(data = data, genos = genos, geno.pos = gp,
                                     ploidy.standardization = 4, threshold.n.clusters = 3,
                                     n.cores = 1, out_filename = NULL, verbose = FALSE))
  chrs <- sort(unique(as.character(st$data$Chr)))
  # Treat the 2nd chromosome as Y-specific; "female" Tetraploid2 has no reads on it
  st$data$R[st$data$SampleName == "Tetraploid2" & st$data$Chr == chrs[2]] <- 0
  sex_chr <- list(Y = chrs[2], PAR = chrs[1])
  n_y <- length(unique(st$data$MarkerName[st$data$Chr == chrs[2]]))

  run <- function(s, ...) suppressWarnings(hmm_estimate_CN(
    st, sample_id = s, chr = chrs[1], snps_per_window = 10, min_snps_per_window = 5,
    cn_grid = 2:4, exp_ploidy = 4, verbose = FALSE,
    sex_system = "XY", sex_chr = sex_chr, ...))

  # chr = chrs[1] drops the Y chromosome, but the score still uses its markers
  m <- run("Tetraploid1")
  expect_false(chrs[2] %in% m$by_marker$Chr)
  expect_equal(m$params$sex_het_depth$n_markers, n_y)
  expect_true(m$params$sex_het_depth$present)
  expect_equal(m$params$sex_het_depth$min_depth, 6)
  expect_equal(m$params$sex_het_depth$min_markers, 2)

  f <- run("Tetraploid2")
  expect_equal(f$params$sex_het_depth$score, sprintf("0/%d", n_y))
  expect_false(f$params$sex_het_depth$present)

  # No X in sex_chr: with Y present PAR cannot be split, so SexChr and Sex are "?"
  expect_equal(unique(karyotype_notation(m)$Sex), "?")
  expect_equal(unique(karyotype_notation(m)$SexChr), "?")
  expect_equal(unique(karyotype_notation(f)$Sex), "female")

  # thresholds are respected and validated
  strict <- run("Tetraploid1", sex_min_markers = n_y + 1)
  expect_false(strict$params$sex_het_depth$present)
  expect_error(run("Tetraploid1", sex_min_depth = -1), "sex_min_depth")

  # no XY system -> no score
  expect_null(suppressWarnings(hmm_estimate_CN(
    st, sample_id = "Tetraploid1", chr = chrs[1], snps_per_window = 10,
    min_snps_per_window = 5, cn_grid = 2:4, exp_ploidy = 4, verbose = FALSE))$params$sex_het_depth)

  # passed through hmm_estimate_CN_multi(...)
  multi <- suppressWarnings(suppressMessages(hmm_estimate_CN_multi(
    st, sample_ids = c("Tetraploid1", "Tetraploid2"), n_cores = 1, chr = chrs[1],
    snps_per_window = 10, min_snps_per_window = 5, cn_grid = 2:4, exp_ploidy = 4,
    verbose = FALSE, sex_system = "XY", sex_chr = sex_chr, sex_min_depth = 6)))
  ct <- count_types(multi)
  expect_equal(ct$sex$Sex[match(c("Tetraploid1", "Tetraploid2"), ct$sex$Sample)],
               c("?", "female"))
})

test_that("prep_baf_ll neutralizes no-BAF windows and penalizes failed fits", {
  # rows: informative, CN3 fit failed, no BAF, every state failed
  m <- rbind(c(-5, -8, -12), c(-5, -Inf, -12), c(0, 0, 0), c(-Inf, -Inf, -Inf))
  n <- c(10, 10, 0, 10)

  for (cs in c(TRUE, FALSE)) {
    out <- sapply(1:3, function(k) prep_baf_ll(m, k, n, correct_scale = cs))
    expect_true(all(is.finite(out)), info = paste("correct_scale", cs))
    # failed state is strictly worse than every successful fit in its window
    expect_lt(out[2, 2], min(out[2, -2]))
    expect_equal(out[2, 2], -1e6)
    # no-BAF and all-failed windows are neutral (equal across states)
    expect_equal(out[3, ], c(0, 0, 0))
    expect_equal(out[4, ], c(0, 0, 0))
  }
  # informative windows are scaled by n_baf only when requested
  expect_equal(prep_baf_ll(m, 1, n, TRUE)[1], -0.5)
  expect_equal(prep_baf_ll(m, 1, n, FALSE)[1], -5)
  # NaN (e.g. failed fit divided by n_baf) is treated as a failure too
  m_nan <- m; m_nan[1, 3] <- NaN
  expect_equal(prep_baf_ll(m_nan, 3, n, TRUE)[1], -1e6)
})

test_that("em_hmm_cn does not let a failed CN-specific BAF fit win", {
  cn_grid <- 2:4
  K <- length(cn_grid)
  state_ids <- as.character(cn_grid)
  W <- 6
  # CN3 fit failed (-Inf) in windows 1-4, CN4 is the best successful fit;
  # window 5 has no BAF, window 6 has every state failed
  ll_baf <- cbind(c(-5, -5, -5, -5, 0, -Inf),
                  c(-Inf, -Inf, -Inf, -Inf, 0, -Inf),
                  c(-3, -3, -3, -3, 0, -Inf))
  n_baf <- c(10, 10, 10, 10, 0, 10)
  A <- matrix(0.05, K, K); diag(A) <- 0.9; A <- A / rowSums(A)
  pi0 <- rep(1 / K, K)
  mu <- stats::setNames(c(-1, 0, 1), state_ids)
  z <- rep(0.9, W)

  run_em <- function(baf_only, w_baf) {
    em_hmm_cn(cn_grid = cn_grid, mu = mu, K = K, state_ids = state_ids,
              sig = 0.5, z = z, z_only = FALSE, ll_baf_matrix = ll_baf,
              n_baf = n_baf, w_baf = w_baf, correct_scale = TRUE,
              A = A, pi0 = pi0, W = W, max_iter = 20, verbose = FALSE,
              baf_only = baf_only)
  }

  for (bo in c(TRUE, FALSE)) {
    res <- run_em(baf_only = bo, w_baf = rep(1, W))
    expect_true(all(is.finite(res$ll_em)), info = paste("baf_only", bo))
    expect_equal(rowSums(res$gamma), rep(1, W), tolerance = 1e-8)
    # the failed state never gets posterior mass where BAF was observed
    expect_true(all(res$gamma[1:4, 2] < 1e-6), info = paste("baf_only", bo))  # col 2 = CN3
    expect_equal(unname(apply(res$gamma[1:4, ], 1, which.max)), rep(3L, 4))
    # neutralized windows (no BAF / all states failed): equal emission across states
    if (bo) expect_equal(unname(res$ll_em[5:6, ]), matrix(0, 2, K))
  }

  # combined emission with zero BAF weight: no NaN from 0 * -Inf
  res0 <- run_em(baf_only = FALSE, w_baf = c(0, 0, 1, 1, 0, 0))
  expect_true(all(is.finite(res0$ll_em)))
  expect_true(all(is.finite(res0$gamma)))
})
