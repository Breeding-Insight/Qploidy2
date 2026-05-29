#' Compute Expected Copy-Number Segregation in Polyploid Crosses
#'
#' Calculates the expected distribution of progeny copy-number (CN) states
#' resulting from a cross between two polyploid parents, supporting even and
#' odd ploidies and interploid (e.g. 4x × 2x) crosses.
#'
#' Each parent's CN state can be supplied either as a per-homolog integer vector
#' (e.g. \code{c(-1, 0, 0, 0)} for a tetraploid carrying one deletion) or as a
#' single integer representing the total CN deviation spread one-per-homolog
#' (e.g. \code{-2} expands to \code{c(-1, -1, 0, 0)} for a tetraploid). The
#' scalar form requires the ploidy to be known via \code{pop_ploidy},
#' \code{ploidy_P}, or \code{ploidy_Q}.
#'
#' Gametes are formed by randomly sampling \code{floor(ploidy/2)} homologs (even
#' ploidy) or both \code{floor(ploidy/2)} and \code{ceiling(ploidy/2)} homologs
#' with equal probability (odd ploidy). Progeny CN is the sum of the two
#' parental gamete contributions.
#'
#' @param cn_P CN state of parent P. Either an integer vector of per-homolog
#'   values (length = ploidy of P) or a single integer representing the total CN
#'   deviation (requires \code{pop_ploidy} or \code{ploidy_P}).
#' @param cn_Q CN state of parent Q. Same format as \code{cn_P}.
#' @param pop_ploidy Integer. Shared ploidy for both parents. Used to expand
#'   scalar \code{cn_P}/\code{cn_Q} inputs when both parents have the same
#'   ploidy. Overridden by \code{ploidy_P}/\code{ploidy_Q} when those are set.
#' @param ploidy_P Integer. Ploidy of parent P. Overrides \code{pop_ploidy} for
#'   parent P. Required when passing a scalar \code{cn_P} in an interploid
#'   cross.
#' @param ploidy_Q Integer. Ploidy of parent Q. Overrides \code{pop_ploidy} for
#'   parent Q. Required when passing a scalar \code{cn_Q} in an interploid
#'   cross.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{progeny_ploidy}{Integer. Ploidy of the progeny class (sum of gamete
#'       sizes from each parent).}
#'     \item{cn}{Integer. Copy-number deviation of the progeny class.}
#'     \item{prob}{Numeric. Expected probability of that progeny class.}
#'   }
#'   Rows are sorted by \code{progeny_ploidy} then \code{cn}.
#'
#' @examples
#' # Tetraploid x tetraploid, both parents carry one deletion
#' segreg_poly_cn(-1, -1, pop_ploidy = 4)
#'
#' # Equivalent using per-homolog vectors
#' segreg_poly_cn(c(-1, 0, 0, 0), c(-1, 0, 0, 0))
#'
#' # Interploid cross: tetraploid x diploid
#' segreg_poly_cn(-1, -1, ploidy_P = 4, ploidy_Q = 2)
#'
#' # Triploid x tetraploid (odd ploidy supported)
#' segreg_poly_cn(-1, -1, ploidy_P = 3, ploidy_Q = 4)
#'
#' @export
segreg_poly_cn <- function(cn_P, cn_Q, pop_ploidy = NULL, ploidy_P = NULL, ploidy_Q = NULL) {
  # cn_P and cn_Q can be:
  #   - integer vectors of per-homolog CN (length = ploidy), e.g. c(-1,0,0,0)
  #   - single integers representing total CN spread across homologs (ploidy required),
  #     e.g. segreg_poly_cn(-1, -1, pop_ploidy = 4) expands to c(-1,0,0,0) each
  # pop_ploidy: shared ploidy for both parents (shorthand when parents have same ploidy)
  # ploidy_P / ploidy_Q: individual parent ploidies (override pop_ploidy per parent)

  expand_cn <- function(cn, ploidy) {
    if (length(cn) == 1) {
      # Deletions cannot exceed ploidy (can't lose more copies than you have)
      if (cn < -ploidy) stop("cn cannot be less than -ploidy (", -ploidy, ")")
      if (cn == 0) return(integer(ploidy))
      if (cn < 0) {
        # e.g. -2 -> c(-1,-1,0,0) for ploidy=4
        v <- integer(ploidy)
        v[seq_len(-cn)] <- -1L
        return(v)
      }
      # Duplications: distribute cn evenly across homologs.
      # e.g. cn=5, ploidy=4 -> c(2,1,1,1); cn=8, ploidy=4 -> c(2,2,2,2)
      base    <- cn %/% ploidy
      remainder <- cn %% ploidy
      v <- rep(base, ploidy)
      if (remainder > 0) v[seq_len(remainder)] <- v[seq_len(remainder)] + 1L
      return(v)
    }
    cn
  }

  # Resolve per-parent ploidy: explicit ploidy_P/Q > pop_ploidy > vector length
  resolve_ploidy <- function(cn, ploidy_arg) {
    if (!is.null(ploidy_arg)) ploidy_arg else if (!is.null(pop_ploidy)) pop_ploidy else length(cn)
  }

  pl_P <- resolve_ploidy(cn_P, ploidy_P)
  pl_Q <- resolve_ploidy(cn_Q, ploidy_Q)

  cn_P <- expand_cn(cn_P, pl_P)
  cn_Q <- expand_cn(cn_Q, pl_Q)

  if (length(cn_P) != pl_P) stop("length of cn_P does not match ploidy_P")
  if (length(cn_Q) != pl_Q) stop("length of cn_Q does not match ploidy_Q")

  # For even ploidy: gametes of size ploidy/2.
  # For odd ploidy: gametes of size floor(ploidy/2) or ceil(ploidy/2),
  #   all configurations equally probable across both sizes.
  # Returns a data.frame with columns: cn, gamete_size, prob
  gamete_dist <- function(cn_parent) {
    p <- length(cn_parent)
    k_vals <- unique(c(floor(p / 2), ceiling(p / 2)))
    combos <- unlist(lapply(k_vals, function(k) {
      combn(p, k, simplify = FALSE)
    }), recursive = FALSE)
    cn_vals   <- sapply(combos, function(idx) sum(cn_parent[idx]))
    size_vals <- sapply(combos, length)
    df <- aggregate(prob ~ cn + gamete_size,
                    data = data.frame(cn = cn_vals, gamete_size = size_vals,
                                      prob = 1 / length(combos)),
                    FUN = sum)
    return(df)
  }

  gP <- gamete_dist(cn_P)
  gQ <- gamete_dist(cn_Q)

  # --- Convolve gamete distributions ---
  # All pairwise combinations of P and Q gametes
  combos <- merge(gP, gQ, by = NULL, suffixes = c("_P", "_Q"))
  combos$progeny_ploidy <- combos$gamete_size_P + combos$gamete_size_Q
  combos$cn             <- combos$cn_P + combos$cn_Q
  combos$prob           <- combos$prob_P * combos$prob_Q

  progeny <- aggregate(prob ~ progeny_ploidy + cn, data = combos, FUN = sum)
  progeny <- progeny[order(progeny$progeny_ploidy, progeny$cn), ]
  rownames(progeny) <- NULL

  return(progeny)
}


#' Test Copy-Number Segregation Against Expected Frequencies
#'
#' For every combination of parent P and parent Q, computes the expected CN
#' segregation frequencies using \code{\link{segreg_poly_cn}}, then tests
#' whether the observed progeny CN counts match those expectations via a
#' chi-square goodness-of-fit test.
#'
#' @param dosages A \code{data.frame} containing at minimum columns for sample
#'   name, marker name, and CN call (see \code{sample_col}, \code{marker_col},
#'   \code{cn_col}).
#' @param parent_names_P Character vector of sample names belonging to parent
#'   group P (e.g. all clones of one parent).
#' @param parent_names_Q Character vector of sample names belonging to parent
#'   group Q (e.g. all clones of the other parent).
#' @param ploidy Integer. Expected ploidy of the population. CN deviations are
#'   computed as \code{CN_call - ploidy}. Used as the fallback when
#'   \code{chrom_col} is \code{NULL} or a chromosome mode cannot be determined.
#' @param chrom_col Name of the column in \code{dosages} containing chromosome
#'   identifiers (e.g. \code{"Chr"}). When provided, the modal CN across
#'   progeny is computed per chromosome and used as the effective ploidy for
#'   each marker's CN-deviation and for \code{\link{segreg_poly_cn}}. This
#'   means a marker with \code{CN_call == mode(chrom)} will have
#'   \code{cn_dev = 0} and \code{ploidy_P}/\code{ploidy_Q} = mode(chrom), while
#'   a marker with a different CN uses \code{cn_dev = CN_call - mode(chrom)}
#'   with \code{ploidy_P}/\code{ploidy_Q} = mode(chrom). If \code{NULL}
#'   (default), the global \code{ploidy} is used uniformly.
#' @param sample_col Name of the column containing sample identifiers. Default
#'   \code{"SampleName"}.
#' @param marker_col Name of the column containing marker identifiers. Default
#'   \code{"MarkerName"}.
#' @param cn_col Name of the column containing integer CN calls. Default
#'   \code{"CN_call"}.
#' @param progeny_names Character vector of sample names to treat as progeny.
#'   If \code{NULL} (default), all samples not listed in \code{parent_names_P}
#'   or \code{parent_names_Q} are used as progeny.
#' @param error_rate Numeric in (0, 1). A small probability substituted for any
#'   zero-probability expected class that has at least one observed count.
#'   This allows the chi-square test to run even when the theoretical model
#'   assigns zero probability to an observed CN state (e.g. due to genotyping
#'   errors or model mis-specification). After substitution all probabilities
#'   are rescaled to sum to 1. Default \code{0.005}.
#'
#' @return A named list of class \code{tested_cn_segregation} with two
#'   elements:
#'   \describe{
#'     \item{summary}{A \code{data.frame} with one row per (parent combination,
#'       marker), with columns \code{parent_P}, \code{parent_Q}, \code{marker},
#'       \code{cn_dev_P}, \code{cn_dev_Q}, \code{ploidy_P}, \code{ploidy_Q}, \code{p.value}.}
#'     \item{cn_detail}{A \code{data.frame} with one row per (parent
#'       combination, marker, CN state), with columns \code{parent_P},
#'       \code{parent_Q}, \code{marker}, \code{cn}, \code{prob}, \code{count}.}
#'   }
#'
#' @examples
#' \dontrun{
#' test_cn_segregation(
#'   dosages        = my_dosages,
#'   parent_names_P = c("ParentA_rep1", "ParentA_rep2"),
#'   parent_names_Q = c("ParentB_rep1"),
#'   ploidy         = 4
#' )
#' }
#'
#' @export
test_cn_segregation <- function(dosages,
                                parent_names_P,
                                parent_names_Q,
                                ploidy,
                                sample_col    = "SampleName",
                                marker_col    = "MarkerName",
                                cn_col        = "CN_call",
                                chrom_col     = NULL,
                                progeny_names = NULL,
                                error_rate    = 0.0001) {

  if (!is.numeric(error_rate) || length(error_rate) != 1 ||
      error_rate <= 0 || error_rate >= 1)
    stop("`error_rate` must be a single numeric value in (0, 1).")

  all_parent_names <- c(parent_names_P, parent_names_Q)

  # ── 1. Split progenies and parents ──────────────────────────────────────────
  if (is.null(progeny_names)) {
    progenies_data <- dosages[!dosages[[sample_col]] %in% all_parent_names, ]
  } else {
    if (!is.character(progeny_names))
      stop("`progeny_names` must be a character vector or NULL.")
    progenies_data <- dosages[ dosages[[sample_col]] %in% progeny_names, ]
  }
  parents_data <- dosages[dosages[[sample_col]] %in% all_parent_names, ]

  # ── 2. Build per-marker effective ploidy (chromosome mode or global ploidy) ──
  if (!is.null(chrom_col)) {
    if (!chrom_col %in% colnames(dosages))
      stop("'chrom_col' column '", chrom_col, "' not found in dosages.")
    # Compute modal CN per chromosome from progeny data
    chr_mode <- tapply(
      progenies_data[[cn_col]],
      progenies_data[[chrom_col]],
      mode
    )
    # Build marker -> chromosome lookup from all data (parents + progenies)
    marker_chr <- unique(dosages[, c(marker_col, chrom_col)])
    colnames(marker_chr) <- c("marker", "chr")
    # Map marker to effective ploidy via chromosome mode
    marker_eff_ploidy <- chr_mode[marker_chr$chr]
    names(marker_eff_ploidy) <- marker_chr$marker
    # Fall back to global ploidy for markers whose chromosome has no mode
    marker_eff_ploidy[is.na(marker_eff_ploidy)] <- ploidy
  } else {
    # All markers use the global ploidy
    all_markers <- unique(dosages[[marker_col]])
    marker_eff_ploidy <- setNames(rep(ploidy, length(all_markers)), all_markers)
  }

  # ── 3. Observed counts per marker × CN deviation ────────────────────────────
  progenies_eff_ploidy <- marker_eff_ploidy[progenies_data[[marker_col]]]
  obs_counts <- aggregate(
    list(count = rep(1L, nrow(progenies_data))),
    by  = list(
      marker = progenies_data[[marker_col]],
      cn_dev = progenies_data[[cn_col]] - progenies_eff_ploidy
    ),
    FUN = length
  )

  # ── 4. Parent CN deviations per marker ──────────────────────────────────────
  parents_eff_ploidy <- marker_eff_ploidy[parents_data[[marker_col]]]
  parents_cn <- data.frame(
    sample     = parents_data[[sample_col]],
    marker     = parents_data[[marker_col]],
    cn_dev     = parents_data[[cn_col]] - parents_eff_ploidy,
    eff_ploidy = parents_eff_ploidy,
    stringsAsFactors = FALSE
  )

  # ── 5. Iterate over all P × Q parent combinations ───────────────────────────
  combinations <- expand.grid(
    parent_P = parent_names_P,
    parent_Q = parent_names_Q,
    stringsAsFactors = FALSE
  )

  results <- vector("list", nrow(combinations))

  for (ci in seq_len(nrow(combinations))) {
    pP <- combinations$parent_P[ci]
    pQ <- combinations$parent_Q[ci]

    cn_P_markers <- parents_cn[parents_cn$sample == pP, c("marker", "cn_dev", "eff_ploidy")]
    cn_Q_markers <- parents_cn[parents_cn$sample == pQ, c("marker", "cn_dev", "eff_ploidy")]

    # Markers with non-missing CN calls for both parents
    shared <- merge(cn_P_markers, cn_Q_markers,
                    by = "marker", suffixes = c("_P", "_Q"))
    # eff_ploidy_P / eff_ploidy_Q are kept separate — parents may differ
    shared <- shared[!is.na(shared$cn_dev_P) & !is.na(shared$cn_dev_Q), ]

    if (nrow(shared) == 0) next

    # ── A. Compute expected distributions for unique (cn_dev_P, cn_dev_Q, eff_ploidy_P, eff_ploidy_Q) tuples
    #       (many markers share the same parental state — compute once, reuse)
    unique_pairs <- unique(shared[, c("cn_dev_P", "cn_dev_Q", "eff_ploidy_P", "eff_ploidy_Q")])

    exp_list <- lapply(seq_len(nrow(unique_pairs)), function(i) {
      exp_df <- tryCatch(
        segreg_poly_cn(unique_pairs$cn_dev_P[i], unique_pairs$cn_dev_Q[i],
                       ploidy_P = unique_pairs$eff_ploidy_P[i],
                       ploidy_Q = unique_pairs$eff_ploidy_Q[i]),
        error = function(e) NULL
      )
      if (is.null(exp_df) || nrow(exp_df) == 0) return(NULL)
      data.frame(cn_dev_P     = unique_pairs$cn_dev_P[i],
                 cn_dev_Q     = unique_pairs$cn_dev_Q[i],
                 eff_ploidy_P = unique_pairs$eff_ploidy_P[i],
                 eff_ploidy_Q = unique_pairs$eff_ploidy_Q[i],
                 cn           = exp_df$cn,
                 prob         = exp_df$prob)
    })
    exp_all <- do.call(rbind, exp_list)
    if (is.null(exp_all) || nrow(exp_all) == 0) next

    # ── B. Attach expected prob to every (marker, cn) row — single join
    shared_exp <- merge(shared[, c("marker", "cn_dev_P", "cn_dev_Q", "eff_ploidy_P", "eff_ploidy_Q")],
                        exp_all,
                        by = c("cn_dev_P", "cn_dev_Q", "eff_ploidy_P", "eff_ploidy_Q"))
    # shared_exp columns: marker, cn_dev_P, cn_dev_Q, cn, prob

    # ── C. Full join with observed counts across all markers at once
    obs_sub <- obs_counts[obs_counts$marker %in% shared$marker, ]
    full_data <- merge(shared_exp[, c("marker", "cn", "prob")],
                       obs_sub,
                       by.x = c("marker", "cn"),
                       by.y = c("marker", "cn_dev"),
                       all = TRUE)
    full_data$count[is.na(full_data$count)] <- 0L
    full_data$prob[is.na(full_data$prob)]   <- 0

    # ── D. Apply error_rate and rescale per marker
    full_data$prob[full_data$prob == 0] <- error_rate
    full_data$prob <- ave(full_data$prob, full_data$marker,
                          FUN = function(p) p / sum(p))

    # Drop markers where all observed counts are zero
    mk_totals <- tapply(full_data$count, full_data$marker, sum)
    full_data  <- full_data[full_data$marker %in% names(mk_totals[mk_totals > 0]), ]

    if (nrow(full_data) == 0) next

    # ── E. Run chi-square test per marker
    mk_split  <- split(full_data, full_data$marker)
    p_values  <- vapply(mk_split, function(d) {
      # Single-cell case: 0 degrees of freedom, chi-square statistic = 0 → p = 1
      # (all observations fall in the one expected class — perfect fit)
      if (nrow(d) == 1L) return(1)
      test <- tryCatch(
        chisq.test(d$count, p = d$prob, rescale.p = TRUE),
        error = function(e) NULL
      )
      if (is.null(test)) NA_real_ else test$p.value
    }, numeric(1))

    # ── F. Build two tidy outputs ────────────────────────────────────────────
    marker_cn_devs <- unique(shared[, c("marker", "cn_dev_P", "cn_dev_Q", "eff_ploidy_P", "eff_ploidy_Q")])

    # summary: one row per marker — parental CN deviations + p.value
    summary_combo <- merge(
      marker_cn_devs,
      data.frame(marker  = names(p_values),
                 p.value = p_values,
                 stringsAsFactors = FALSE,
                 row.names = NULL),
      by = "marker"
    )
    summary_combo$parent_P <- pP
    summary_combo$parent_Q <- pQ
    colnames(summary_combo)[colnames(summary_combo) == "eff_ploidy_P"] <- "ploidy_P"
    colnames(summary_combo)[colnames(summary_combo) == "eff_ploidy_Q"] <- "ploidy_Q"
    summary_combo <- summary_combo[, c("parent_P", "parent_Q", "marker",
                                       "cn_dev_P", "cn_dev_Q", "ploidy_P", "ploidy_Q", "p.value")]
    summary_combo <- summary_combo[order(summary_combo$marker), ]

    # cn_detail: one row per marker × cn — expected prob + observed count
    detail_combo <- full_data[, c("marker", "cn", "prob", "count")]
    detail_combo$parent_P <- pP
    detail_combo$parent_Q <- pQ
    detail_combo <- detail_combo[, c("parent_P", "parent_Q", "marker",
                                     "cn", "prob", "count")]
    detail_combo <- detail_combo[order(detail_combo$marker, detail_combo$cn), ]

    results[[ci]] <- list(summary = summary_combo, cn_detail = detail_combo)
  }

  summary_out <- do.call(rbind, lapply(results, `[[`, "summary"))
  detail_out  <- do.call(rbind, lapply(results, `[[`, "cn_detail"))
  rownames(summary_out) <- NULL
  rownames(detail_out)  <- NULL

  return(structure(
    list(summary = summary_out, cn_detail = detail_out),
    class = "tested_cn_segregation"
  ))
}


#' Plot P-value Histograms from CN Segregation Tests
#'
#' Plots a histogram of p-values returned by \code{\link{test_cn_segregation}},
#' with one facet (or fill colour) per parent combination.
#'
#' @param x A \code{tested_cn_segregation} object (the list returned by
#'   \code{\link{test_cn_segregation}}). The \code{$summary} element is used
#'   for plotting.
#' @param bins Integer. Number of histogram bins. Default \code{30}.
#' @param alpha Numeric in (0, 1]. Bar transparency. Default \code{0.7}.
#' @param facet Logical. If \code{TRUE} (default), each parent combination gets
#'   its own facet panel. If \code{FALSE}, all combinations are overlaid in a
#'   single panel with colours.
#' @param significance_line Numeric or \code{NULL}. Nominal significance
#'   threshold before Bonferroni correction. A vertical dashed line is drawn
#'   at \code{significance_line / n_tests}, where \code{n_tests} is the number
#'   of markers tested per parent combination (or across all combinations when
#'   \code{facet = FALSE}). Set to \code{NULL} to suppress the line.
#'   Default \code{0.05}.
#' @param color_by Character. What to use for colours/fills. One of:
#'   \describe{
#'     \item{\code{"combination"} (default)}{Colour by parent name combination
#'       (e.g. \code{"ParentA × ParentB"}).}
#'     \item{\code{"cn_pair"}}{Colour by the canonical parental CN-deviation
#'       pair (e.g. \code{"-1 × 0"}). The pair is sorted so that
#'       \code{"-1 × 0"} and \code{"0 × -1"} receive the same colour.}
#'   }
#' @param bonferroni Logical. If \code{TRUE} (default), the significance line
#'   is drawn at \code{significance_line / n_tests}. If \code{FALSE}, the raw
#'   \code{significance_line} value is used.
#'
#' @return A \code{ggplot} object.
#'
#' @import ggplot2
#'
#' @export
plot_tested_cn_segregation <- function(x,
                                       bins = 30,
                                       alpha = 0.7,
                                       facet = TRUE,
                                       color_by = c("cn_pair", "combination"),
                                       significance_line = 0.05,
                                       bonferroni = TRUE) {

  color_by <- match.arg(color_by)

  # Use the summary element (one row per marker) for plotting
  x <- x$summary

  # Build colour grouping variable
  if (color_by == "combination") {
    x$colour_group <- paste(x$parent_P, "\u00d7", x$parent_Q)
    legend_title <- "Parent combination"
  } else {
    # Canonical CN pair: sort the two deviations so -1×0 == 0×-1
    x$colour_group <- apply(
      cbind(x$cn_dev_P, x$cn_dev_Q), 1,
      function(v) paste(sort(v), collapse = " \u00d7 ")
    )
    legend_title <- "Parental CN pair"
  }

  # For facets keep combination label regardless of colour_by
  x$combination <- paste(x$parent_P, "\u00d7", x$parent_Q)

  # Compute Bonferroni-corrected threshold per facet group (or globally)
  get_threshold <- function(df, alpha, bonferroni) {
    if (is.null(alpha)) return(NULL)
    if (!bonferroni) return(alpha)
    n <- nrow(df)
    if (n == 0) return(alpha)
    alpha / n
  }

  if (facet) {
    # One threshold per combination panel
    thresholds <- lapply(split(x, x$combination), get_threshold,
                         alpha = significance_line, bonferroni = bonferroni)
    vline_df <- data.frame(
      combination = names(thresholds),
      xintercept  = unlist(thresholds),
      stringsAsFactors = FALSE
    )
  } else {
    # Single global threshold across all tests
    global_thresh <- get_threshold(x, significance_line, bonferroni)
  }

  # Discrete viridis palette scales to any number of levels
  colour_scale_fill  <- scale_fill_viridis_d(option = "turbo",
                                              name = legend_title)
  colour_scale_color <- scale_color_viridis_d(option = "turbo",
                                               name = legend_title)

  if (facet) {
    p <- ggplot(x, aes(x = p.value, fill = colour_group)) +
      geom_histogram(bins = bins, alpha = alpha,
                              colour = "white", linewidth = 0.2) +
      colour_scale_fill +
      facet_wrap(~combination, scales = "free_y") +
      theme(
        strip.text      = element_text(size = 9),
        legend.position = if (color_by == "cn_pair") "right" else "none"
      )
  } else {
    p <- ggplot(x, aes(x = p.value, colour = colour_group)) +
      geom_freqpoly(bins = bins, linewidth = 0.8, alpha = alpha) +
      colour_scale_color +
      ggplot2::theme(
        legend.position = "right"
      )
  }

  p <- p +
    labs(
      x        = "p-value",
      y        = "Number of markers",
      title    = "CN segregation test: p-value distribution",
      subtitle = if (!is.null(significance_line) && bonferroni)
        paste0("Red dashed line: Bonferroni-corrected \u03b1 = ",
               significance_line, " / n markers per combination")
        else if (!is.null(significance_line))
        paste0("Red dashed line: \u03b1 = ", significance_line)
        else NULL
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 9)
    )

  if (!is.null(significance_line)) {
    if (facet) {
      p <- p + geom_vline(
        data     = vline_df,
        aes(xintercept = xintercept),
        linetype  = "dashed",
        colour    = "red",
        linewidth = 0.6
      )
    } else {
      p <- p + geom_vline(
        xintercept = global_thresh,
        linetype   = "dashed",
        colour     = "red",
        linewidth  = 0.6
      )
    }
  }

  return(p)
}


#' Simulate a Dosages Object for CN Segregation Testing
#'
#' Generates a synthetic \code{data.frame} in the format expected by
#' \code{\link{test_cn_segregation}}. Exactly one P × Q parent pair is the
#' "true" cross; progeny CN calls for that combination are drawn from the
#' expected segregation distribution. All other parent combinations produce
#' progeny with CN calls drawn from a uniform distribution, ensuring they are
#' unlikely to pass the chi-square test.
#'
#' @param n_parents_P Integer. Number of parent P candidates. Default \code{3}.
#' @param n_parents_Q Integer. Number of parent Q candidates. Default \code{3}.
#' @param n_markers Integer. Number of markers to simulate. Default \code{500}.
#' @param n_progeny Integer. Number of progeny samples. Default \code{100}.
#' @param ploidy Integer. Population ploidy. Default \code{4}.
#' @param true_P Integer. Index (within 1:\code{n_parents_P}) of the true
#'   parent P. Default \code{1}.
#' @param true_Q Integer. Index (within 1:\code{n_parents_Q}) of the true
#'   parent Q. Default \code{1}.
#' @param prop_informative Numeric in (0, 1]. Proportion of markers that carry
#'   a CN deviation in at least one parent (i.e. are informative for
#'   segregation). The rest are fixed at CN = \code{ploidy} in all parents.
#'   Default \code{0.4}.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#'   Default \code{42}.
#'
#' @return A \code{list} with:
#'   \describe{
#'     \item{dosages}{A \code{data.frame} with columns \code{MarkerName},
#'       \code{SampleName}, \code{Chr}, \code{Position}, \code{X}, \code{Y},
#'       \code{baf}, \code{z}, \code{CN_call}, \code{post_max_CN},
#'       \code{dosage}, \code{post_max_dosage}, matching the Qploidy format.}
#'     \item{true_P}{Name of the true parent P.}
#'     \item{true_Q}{Name of the true parent Q.}
#'     \item{parent_names_P}{Character vector of all parent P names.}
#'     \item{parent_names_Q}{Character vector of all parent Q names.}
#'   }
#'
#' @examples
#' sim <- simulate_cn_segregation(n_parents_P = 3, n_parents_Q = 2,
#'                                n_markers = 200, n_progeny = 80)
#' res <- test_cn_segregation(sim$dosages, sim$parent_names_P,
#'                            sim$parent_names_Q, ploidy = 4)
#' plot_tested_cn_segregation(res)
#'
#' @export
simulate_cn_segregation <- function(n_parents_P    = 3,
                                    n_parents_Q    = 3,
                                    n_markers      = 500,
                                    n_progeny      = 100,
                                    ploidy         = 4,
                                    true_P         = 1,
                                    true_Q         = 1,
                                    prop_informative = 0.4,
                                    seed           = 42) {

  if (!is.null(seed)) set.seed(seed)

  parent_names_P <- paste0("Parent_P", seq_len(n_parents_P))
  parent_names_Q <- paste0("Parent_Q", seq_len(n_parents_Q))
  true_P_name    <- parent_names_P[true_P]
  true_Q_name    <- parent_names_Q[true_Q]
  progeny_names  <- paste0("Progeny_", seq_len(n_progeny))
  marker_names   <- paste0("MK_", seq_len(n_markers))

  # Possible CN deviations: -ploidy/2 to +ploidy/2, excluding impossible
  max_dev  <- ploidy %/% 2
  dev_pool <- seq(-max_dev, max_dev)

  # ── Assign CN deviations to parents per marker ──────────────────────────────
  # True parents get informative deviations; decoys are fixed at 0
  n_inform <- round(n_markers * prop_informative)
  inform_idx <- sample(n_markers, n_inform)

  # True parent CN deviations per marker (scalar, single integer)
  true_P_cn <- integer(n_markers)  # 0 = normal CN
  true_Q_cn <- integer(n_markers)
  true_P_cn[inform_idx] <- sample(dev_pool[dev_pool != 0], n_inform,
                                  replace = TRUE)
  true_Q_cn[inform_idx] <- sample(dev_pool[dev_pool != 0], n_inform,
                                  replace = TRUE)

  # Decoy parents are fixed at CN deviation = 0 (always "normal")
  make_parent_rows <- function(pname, cn_devs) {
    data.frame(
      MarkerName      = marker_names,
      SampleName      = pname,
      Chr             = 1L,
      Position        = seq_len(n_markers) * 1000L,
      X               = NA_real_,
      Y               = NA_real_,
      baf             = NA_real_,
      z               = NA_real_,
      CN_call         = cn_devs + ploidy,
      post_max_CN     = NA_real_,
      dosage          = NA_real_,
      post_max_dosage = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  parent_rows <- vector("list", n_parents_P + n_parents_Q)
  k <- 1L
  for (pn in parent_names_P) {
    devs <- if (pn == true_P_name) true_P_cn else integer(n_markers)
    parent_rows[[k]] <- make_parent_rows(pn, devs)
    k <- k + 1L
  }
  for (pn in parent_names_Q) {
    devs <- if (pn == true_Q_name) true_Q_cn else integer(n_markers)
    parent_rows[[k]] <- make_parent_rows(pn, devs)
    k <- k + 1L
  }

  # ── Simulate progeny CN calls ─────────────────────────────────────────────
  # For each marker, draw progeny CN from the expected segregation distribution
  progeny_cn_mat <- matrix(ploidy, nrow = n_markers, ncol = n_progeny,
                           dimnames = list(marker_names, progeny_names))

  for (mi in seq_len(n_markers)) {
    cp <- true_P_cn[mi]
    cq <- true_Q_cn[mi]

    exp_df <- tryCatch(
      segreg_poly_cn(cp, cq, pop_ploidy = ploidy),
      error = function(e) NULL
    )

    if (is.null(exp_df) || nrow(exp_df) == 0) next

    # Aggregate over progeny_ploidy: sum probs for the same cn deviation
    # (relevant when gametes of different sizes produce the same cn)
    exp_agg <- aggregate(prob ~ cn, data = exp_df, FUN = sum)

    # CN deviations → absolute CN; clamp to valid range [0, 2*ploidy]
    exp_agg$cn_abs <- exp_agg$cn + ploidy
    exp_agg <- exp_agg[exp_agg$cn_abs >= 0 & exp_agg$cn_abs <= 2 * ploidy, ]

    if (nrow(exp_agg) == 0 || sum(exp_agg$prob) == 0) next

    cn_states <- exp_agg$cn_abs
    probs     <- exp_agg$prob / sum(exp_agg$prob)

    # Use sample.int on indices to avoid the R scalar-sample trap:
    # sample(x, ...) when length(x)==1 samples from 1:x, not from c(x)
    idx <- sample.int(length(cn_states), n_progeny, replace = TRUE, prob = probs)
    progeny_cn_mat[mi, ] <- cn_states[idx]
  }

  # Convert matrix to long data.frame
  progeny_rows <- data.frame(
    MarkerName      = rep(marker_names, times = n_progeny),
    SampleName      = rep(progeny_names, each  = n_markers),
    Chr             = 1L,
    Position        = rep(seq_len(n_markers) * 1000L, times = n_progeny),
    X               = NA_real_,
    Y               = NA_real_,
    baf             = NA_real_,
    z               = NA_real_,
    CN_call         = as.integer(progeny_cn_mat),
    post_max_CN     = NA_real_,
    dosage          = NA_real_,
    post_max_dosage = NA_real_,
    stringsAsFactors = FALSE
  )

  dosages <- do.call(rbind, c(parent_rows, list(progeny_rows)))
  rownames(dosages) <- NULL

  return(list(
    dosages        = dosages,
    true_P         = true_P_name,
    true_Q         = true_Q_name,
    parent_names_P = parent_names_P,
    parent_names_Q = parent_names_Q
  ))
}
