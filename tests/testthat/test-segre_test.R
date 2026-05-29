library(testthat)
library(Qploidy2)

# ── segreg_poly_cn ─────────────────────────────────────────────────────────────

test_that("segreg_poly_cn returns a data.frame with expected columns", {
  res <- segreg_poly_cn(0, 0, pop_ploidy = 4)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("progeny_ploidy", "cn", "prob"))
})

test_that("segreg_poly_cn probabilities sum to 1 for normal diploid parents", {
  res <- segreg_poly_cn(0, 0, pop_ploidy = 4)
  expect_equal(sum(res$prob), 1, tolerance = 1e-9)
})

test_that("segreg_poly_cn normal 4x parents produce only cn = 0 progeny", {
  res <- segreg_poly_cn(0, 0, pop_ploidy = 4)
  expect_equal(nrow(res), 1)
  expect_equal(res$cn, 0)
  expect_equal(res$prob, 1)
})

test_that("segreg_poly_cn probabilities sum to 1 for informative parents", {
  res <- segreg_poly_cn(-1, -1, pop_ploidy = 4)
  expect_equal(sum(res$prob), 1, tolerance = 1e-9)
  expect_true(all(res$prob > 0))
  expect_true(nrow(res) > 1)
})

test_that("segreg_poly_cn per-homolog vector input matches scalar input", {
  scalar <- segreg_poly_cn(-1, -1, pop_ploidy = 4)
  vector <- segreg_poly_cn(c(-1, 0, 0, 0), c(-1, 0, 0, 0))
  # same cn classes and probs (rows may differ in progeny_ploidy for odd ploidy)
  scalar_agg <- aggregate(prob ~ cn, data = scalar, FUN = sum)
  vector_agg <- aggregate(prob ~ cn, data = vector, FUN = sum)
  scalar_agg <- scalar_agg[order(scalar_agg$cn), ]
  vector_agg <- vector_agg[order(vector_agg$cn), ]
  expect_equal(scalar_agg$cn,   vector_agg$cn)
  expect_equal(scalar_agg$prob, vector_agg$prob, tolerance = 1e-9)
})

test_that("segreg_poly_cn interploid cross probabilities sum to 1", {
  res <- segreg_poly_cn(-1, -1, ploidy_P = 4, ploidy_Q = 2)
  expect_equal(sum(res$prob), 1, tolerance = 1e-9)
})

test_that("segreg_poly_cn errors when deletion exceeds ploidy", {
  expect_error(segreg_poly_cn(-3, 0, pop_ploidy = 2), "cannot be less than -ploidy")
  # duplications above ploidy are allowed (e.g. cn=3 for ploidy=2)
  expect_no_error(segreg_poly_cn(3, 0, pop_ploidy = 2))
})

# ── simulate_cn_segregation ────────────────────────────────────────────────────

test_that("simulate_cn_segregation returns expected list structure", {
  sim <- simulate_cn_segregation(n_parents_P = 2, n_parents_Q = 2,
                                 n_markers = 50, n_progeny = 20, seed = 1)
  expect_type(sim, "list")
  expect_named(sim, c("dosages", "true_P", "true_Q",
                      "parent_names_P", "parent_names_Q"))
})

test_that("simulate_cn_segregation dosages has correct columns", {
  sim <- simulate_cn_segregation(n_markers = 50, n_progeny = 20, seed = 1)
  expect_s3_class(sim$dosages, "data.frame")
  expect_true(all(c("MarkerName", "SampleName", "CN_call") %in%
                    colnames(sim$dosages)))
})

test_that("simulate_cn_segregation produces correct number of rows", {
  np <- 2; nq <- 2; nm <- 30; nprog <- 15
  sim <- simulate_cn_segregation(n_parents_P = np, n_parents_Q = nq,
                                 n_markers = nm, n_progeny = nprog, seed = 2)
  expected_rows <- (np + nq + nprog) * nm
  expect_equal(nrow(sim$dosages), expected_rows)
})

test_that("simulate_cn_segregation true parent names are within parent vectors", {
  sim <- simulate_cn_segregation(n_parents_P = 3, n_parents_Q = 3,
                                 true_P = 2, true_Q = 3, seed = 5)
  expect_equal(sim$true_P, "Parent_P2")
  expect_equal(sim$true_Q, "Parent_Q3")
  expect_true(sim$true_P %in% sim$parent_names_P)
  expect_true(sim$true_Q %in% sim$parent_names_Q)
})

test_that("simulate_cn_segregation CN_call values are non-negative integers", {
  sim <- simulate_cn_segregation(n_markers = 50, n_progeny = 30, seed = 3)
  cn <- sim$dosages$CN_call
  expect_true(is.numeric(cn))
  expect_true(all(cn >= 0))
})

test_that("simulate_cn_segregation is reproducible with same seed", {
  sim1 <- simulate_cn_segregation(n_markers = 50, n_progeny = 20, seed = 99)
  sim2 <- simulate_cn_segregation(n_markers = 50, n_progeny = 20, seed = 99)
  expect_identical(sim1$dosages$CN_call, sim2$dosages$CN_call)
})

# ── test_cn_segregation ────────────────────────────────────────────────────────

test_that("test_cn_segregation returns a list with summary and cn_detail", {
  sim <- simulate_cn_segregation(n_parents_P = 2, n_parents_Q = 2,
                                 n_markers = 100, n_progeny = 60, seed = 7)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  expect_s3_class(res, "tested_cn_segregation")
  expect_type(res, "list")
  expect_named(res, c("summary", "cn_detail"))
})

test_that("test_cn_segregation summary has expected columns and one row per marker", {
  sim <- simulate_cn_segregation(n_parents_P = 2, n_parents_Q = 2,
                                 n_markers = 100, n_progeny = 60, seed = 7)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  expect_true(all(c("parent_P", "parent_Q", "marker",
                    "cn_dev_P", "cn_dev_Q", "p.value") %in%
                    colnames(res$summary)))
  # one unique marker per combination in summary
  n_dupes <- sum(duplicated(paste(res$summary$parent_P,
                                  res$summary$parent_Q,
                                  res$summary$marker)))
  expect_equal(n_dupes, 0L)
})

test_that("test_cn_segregation cn_detail has expected columns", {
  sim <- simulate_cn_segregation(n_parents_P = 2, n_parents_Q = 2,
                                 n_markers = 100, n_progeny = 60, seed = 7)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  expect_true(all(c("parent_P", "parent_Q", "marker",
                    "cn", "prob", "count") %in% colnames(res$cn_detail)))
})

test_that("test_cn_segregation p-values are in [0, 1]", {
  sim <- simulate_cn_segregation(n_parents_P = 2, n_parents_Q = 2,
                                 n_markers = 100, n_progeny = 60, seed = 7)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  pv <- res$summary$p.value[!is.na(res$summary$p.value)]
  expect_true(all(pv >= 0 & pv <= 1))
})

test_that("test_cn_segregation result has tested_cn_segregation class", {
  sim <- simulate_cn_segregation(n_parents_P = 2, n_parents_Q = 2,
                                 n_markers = 50, n_progeny = 40, seed = 8)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  expect_s3_class(res, "tested_cn_segregation")
  expect_type(res, "list")
})

test_that("test_cn_segregation covers all P x Q combinations", {
  np <- 3; nq <- 2
  sim <- simulate_cn_segregation(n_parents_P = np, n_parents_Q = nq,
                                 n_markers = 80, n_progeny = 50, seed = 9)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  combos_found <- unique(paste(res$summary$parent_P, res$summary$parent_Q))
  expect_equal(length(combos_found), np * nq)
})

test_that("test_cn_segregation cn_detail has multiple cn rows per marker when informative", {
  sim <- simulate_cn_segregation(n_parents_P = 1, n_parents_Q = 1,
                                 n_markers = 50, n_progeny = 80,
                                 prop_informative = 1.0, seed = 21)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  rows_per_marker <- table(res$cn_detail$marker[
    res$cn_detail$parent_P == sim$parent_names_P[1] &
    res$cn_detail$parent_Q == sim$parent_names_Q[1]])
  expect_true(any(rows_per_marker > 1))
})

test_that("test_cn_segregation summary and cn_detail share the same markers", {
  sim <- simulate_cn_segregation(n_parents_P = 2, n_parents_Q = 2,
                                 n_markers = 60, n_progeny = 50, seed = 22)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  summary_keys <- sort(unique(paste(res$summary$parent_P,
                                    res$summary$parent_Q,
                                    res$summary$marker)))
  detail_keys  <- sort(unique(paste(res$cn_detail$parent_P,
                                    res$cn_detail$parent_Q,
                                    res$cn_detail$marker)))
  expect_equal(summary_keys, detail_keys)
})

test_that("true parent combination has higher median p-value than decoys", {
  sim <- simulate_cn_segregation(n_parents_P = 3, n_parents_Q = 3,
                                 n_markers = 300, n_progeny = 100,
                                 prop_informative = 0.5, seed = 42)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)

  s <- res$summary
  s$combination <- paste(s$parent_P, s$parent_Q)
  true_combo    <- paste(sim$true_P, sim$true_Q)

  medians <- tapply(s$p.value, s$combination, median, na.rm = TRUE)

  # True combination should have the highest median p-value
  expect_equal(names(which.max(medians)), true_combo)
})

test_that("test_cn_segregation errors on invalid error_rate", {
  sim <- simulate_cn_segregation(n_markers = 30, n_progeny = 20, seed = 1)
  expect_error(
    test_cn_segregation(sim$dosages, sim$parent_names_P,
                        sim$parent_names_Q, ploidy = 4, error_rate = 0),
    "error_rate"
  )
  expect_error(
    test_cn_segregation(sim$dosages, sim$parent_names_P,
                        sim$parent_names_Q, ploidy = 4, error_rate = 1.5),
    "error_rate"
  )
})

# ── plot_tested_cn_segregation ─────────────────────────────────────────────────

test_that("plot_tested_cn_segregation returns a ggplot object (facet = TRUE)", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_cn_segregation(n_parents_P = 2, n_parents_Q = 2,
                                 n_markers = 80, n_progeny = 50, seed = 11)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  p <- plot_tested_cn_segregation(res, facet = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_tested_cn_segregation returns a ggplot object (facet = FALSE)", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_cn_segregation(n_parents_P = 2, n_parents_Q = 2,
                                 n_markers = 80, n_progeny = 50, seed = 11)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  p <- plot_tested_cn_segregation(res, facet = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_tested_cn_segregation works with significance_line = NULL", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_cn_segregation(n_markers = 60, n_progeny = 40, seed = 12)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  expect_s3_class(plot_tested_cn_segregation(res, significance_line = NULL),
                  "ggplot")
})

test_that("plot_tested_cn_segregation works with bonferroni = FALSE", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_cn_segregation(n_markers = 60, n_progeny = 40, seed = 13)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  expect_s3_class(plot_tested_cn_segregation(res, bonferroni = FALSE),
                  "ggplot")
})

test_that("plot_tested_cn_segregation with only one combination of parents", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_cn_segregation(n_parents_P = 1, n_parents_Q = 1,
                                 n_markers = 100, n_progeny = 40, seed = 13)
  res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
                             sim$parent_names_Q, ploidy = 4)
  p <- plot_tested_cn_segregation(res, bonferroni =  TRUE)
  expect_s3_class(p, "ggplot")
})
