#' Split Markers by Sub-genome Specificity
#'
#' Classifies markers according to which sub-genome(s) they carry signal in,
#' based on signal intensity (R) observed in a set of diploid reference samples
#' representing each sub-genome.  Works for any number of sub-genomes (A, B, C,
#' …) and returns every combination (A, B, AB, AC, ABC, …) that is present.
#'
#' For each marker the median R is computed within every sub-genome group.  A
#' sub-genome is considered "active" for that marker when its median R is at
#' least \code{threshold} times the maximum median R observed across all
#' sub-genomes.  Markers are then labelled by the concatenation of all active
#' sub-genome identifiers (e.g. \code{"AB"}, \code{"C"}).
#'
#' @param data A \code{data.frame} containing at least columns for sample ID,
#'   marker name, and signal intensity R.
#' @param subgenome_samples A named list of character vectors.  Each element
#'   gives the sample IDs of the diploid reference individuals for one
#'   sub-genome, and the element name is the sub-genome identifier (e.g.
#'   \code{list(A = c("s1","s2"), B = c("s3","s4"))}).
#' @param R_col Name of the column in \code{data} that holds signal intensity
#'   values.  Default \code{"R"}.
#' @param marker_col Name of the column in \code{data} that holds marker
#'   identifiers.  Default \code{"MarkerName"}.
#' @param sample_col Name of the column in \code{data} that holds sample
#'   identifiers.  Default \code{"SampleName"}.
#' @param threshold Numeric in (0, 1].  A sub-genome is considered active for a
#'   marker when \code{median_R / max_median_R >= threshold}.  Default
#'   \code{0.8}.
#' @param plot Logical.  If \code{TRUE}, produces a boxplot of R values per
#'   diploid sub-genome group faceted by marker-type combination.  Default
#'   \code{TRUE}.
#' @param plot_filename Character path for saving the plot (e.g.
#'   \code{"markers.png"}).  Any format recognised by \code{ggplot2::ggsave}
#'   is accepted.  When \code{NULL} (default) the plot is printed to the
#'   current graphics device.
#'
#' @return A named list of character vectors.  Each element contains the marker
#'   names belonging to a particular sub-genome combination (e.g. \code{"A"},
#'   \code{"B"}, \code{"AB"}).  List items are ordered first by the number of
#'   sub-genomes in the combination, then alphabetically.  Markers with no
#'   active sub-genome are returned under the name \code{""} if any exist.
#'
#' @export
split_mk_type <- function(data,
                          subgenome_samples,
                          R_col         = "R",
                          marker_col    = "MarkerName",
                          sample_col    = "SampleName",
                          threshold     = 0.8,
                          plot          = TRUE,
                          plot_filename = NULL) {

  subgenome_names <- names(subgenome_samples)
  if (is.null(subgenome_names) || any(nchar(subgenome_names) == 0L)) {
    stop("'subgenome_samples' must be a named list with non-empty names.")
  }
  if (!threshold > 0 || !threshold <= 1) {
    stop("'threshold' must be in (0, 1].")
  }
  for (col in c(R_col, marker_col, sample_col)) {
    if (!col %in% colnames(data)) stop("Column '", col, "' not found in 'data'.")
  }

  # Subset to diploid reference samples only
  all_ref <- unlist(subgenome_samples, use.names = FALSE)
  df <- data[data[[sample_col]] %in% all_ref, c(sample_col, marker_col, R_col),
             drop = FALSE]

  # Assign sub-genome label to each row
  df$subgenome <- NA_character_
  for (sg in subgenome_names) {
    df$subgenome[df[[sample_col]] %in% subgenome_samples[[sg]]] <- sg
  }

  # Median R per marker x sub-genome
  med_R <- aggregate(
    df[[R_col]],
    by  = list(marker = df[[marker_col]], subgenome = df$subgenome),
    FUN = function(x) median(x, na.rm = TRUE)
  )
  names(med_R)[3] <- "med_R"

  # Pivot to wide: rows = markers, columns = sub-genomes
  med_wide <- reshape(med_R,
                      idvar     = "marker",
                      timevar   = "subgenome",
                      direction = "wide")
  colnames(med_wide) <- sub("^med_R\\.", "", colnames(med_wide))

  # Ensure a column exists for every sub-genome (fill absent ones with NA)
  for (sg in subgenome_names) {
    if (!sg %in% colnames(med_wide)) med_wide[[sg]] <- NA_real_
  }

  med_mat <- as.matrix(med_wide[, subgenome_names, drop = FALSE])
  rownames(med_mat) <- med_wide$marker

  # Per-marker maximum; guard against all-NA rows
  row_max <- apply(med_mat, 1L, function(x) {
    m <- suppressWarnings(max(x, na.rm = TRUE))
    if (is.infinite(m)) NA_real_ else m
  })

  # Active sub-genome: ratio to maximum >= threshold (NA → inactive)
  ratio_mat  <- med_mat / row_max
  active_mat <- !is.na(ratio_mat) & ratio_mat >= threshold
  colnames(active_mat) <- subgenome_names

  # Build combination label for every marker
  combo_labels <- apply(active_mat, 1L, function(row) {
    paste(subgenome_names[row], collapse = "")
  })

  result <- split(med_wide$marker, combo_labels)

  # Sort: single sub-genomes first, then pairs, then triples, etc.; alpha within
  result <- result[order(nchar(names(result)), names(result))]

  if (plot) {
    # Attach combination label to each observation for faceting
    marker_type <- data.frame(
      marker    = med_wide$marker,
      mk_type   = combo_labels,
      stringsAsFactors = FALSE
    )
    plot_df <- merge(df, marker_type,
                     by.x = marker_col, by.y = "marker",
                     all.x = FALSE)
    plot_df$mk_type <- factor(
      plot_df$mk_type,
      levels = names(result)
    )
    plot_df$subgenome <- factor(plot_df$subgenome, levels = subgenome_names)

    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = subgenome, y = .data[[R_col]])) +
      ggplot2::geom_boxplot(fill        = "steelblue",
                            color       = "black",
                            outlier.shape = 21) +
      ggplot2::facet_wrap(~ mk_type) +
      ggplot2::labs(
        title = "R values by diploid sub-genome group and marker type",
        x     = "Sub-genome (diploid group)",
        y     = "R"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "lightgray"),
        strip.text       = ggplot2::element_text(face = "bold"),
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
      )

    if (!is.null(plot_filename)) {
      ggplot2::ggsave(plot_filename, plot = p)
    } else {
      print(p)
    }
  }

  return(result)
}