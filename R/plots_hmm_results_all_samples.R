# Suppress R CMD check notes for non-standard evaluation (dplyr/ggplot2)
if (getRversion() >= "2.15.1") utils::globalVariables(c(
  # compare_cn_track
  "CN_call", "CN_reliability", "Sample", "Chr", "Start", "End", "post_max",
  "w_baf", "n_het", "n_snps", "total_bp", "w",
  # compare_cn_track_summary
  "w_len", "tot", "n_cn", "dom_cn", "sample_cn", "cell_type", "cn_label",
  "has_mixed", "has_aneuploid", "sample_type", "y_pos", "prev_cn", "prev_type",
  "line_y", "mid_y", "lbl", "x_pos"
))

#' Compare CNV tracks across samples
#'
#' Plots CNV windows as horizontal segments for multiple samples, faceted by chromosome.
#' Segments are colored by copy number (CN_call) using a high-contrast palette:
#' baseline CN (most frequent, weighted by window length) is black; losses are blue;
#' gains are red. Transparency reflects posterior confidence (post_max).
#'
#' @param hmm_CN An object of class hmm_CN (output from hmm_estimate_CN), containing a result data frame with columns: Sample, Chr, Start, End, CN_call, post_max, etc.
#' @param samples_to_plot Character vector of sample IDs to include. If NULL, the first sample in the data is plotted.
#' @param chromosomes Optional character vector of chromosomes to include (matching result$Chr). If NULL, plots all chromosomes present after sample filtering.
#' @param facet_ncol Number of columns for facet_wrap. If NULL, will be determined automatically unless facet_nrow is set.
##' @param facet_nrow Number of rows for facet_wrap. If set, facet_ncol will be determined automatically unless both are set.
##' @param gray_CN Integer or NULL. If provided, this copy-number value is colored gray (used as the baseline color). All CNs below it are colored blue and above red. If NULL (default), the baseline is auto-detected as the most frequent CN weighted by window length.
#' @param add_het Logical. If TRUE, a heterozygosity column is added to the right side of the plot, with one square per sample colored from blue (low) to red (high). Requires \code{hmm_dosage_calls}. Default is TRUE.
#' @param hmm_dosage_calls An object of class \code{hmm_dosage_calls} (inherits \code{data.frame}) with columns \code{SampleName}, \code{Chr}, \code{dosage}, and optionally others. Used to compute per-sample heterozygosity as the fraction of markers with dosage not equal to 0 or 1 (among non-NA dosage values). Required when \code{add_het = TRUE}.
#' @param interactive Logical. If TRUE, returns an interactive \code{plotly} figure instead of a static ggplot. Hovering over segments shows: Sample, Chr, Start, End, CN, post_max, w_baf, and Heterozygosity (if \code{add_het = TRUE} and \code{hmm_dosage_calls} is provided). Requires the \pkg{plotly} package. Default is FALSE.
#'
#' @return A ggplot/ggarrange object (static), or a plotly figure when \code{interactive = TRUE}.
#'
#' @importFrom dplyr filter group_by summarise arrange desc mutate left_join
#' @importFrom ggplot2 ggplot geom_segment geom_tile geom_point facet_wrap
#' @importFrom ggplot2 scale_x_continuous scale_color_manual scale_alpha scale_fill_gradient
#' @importFrom ggplot2 labs theme_bw theme aes
#' @importFrom ggplot2 element_blank element_rect element_text
#' @importFrom ggpubr ggarrange get_legend
#'
#' @export
compare_cn_track <- function(hmm_CN,
                             samples_to_plot = NULL,
                             chromosomes = NULL,
                             facet_ncol = NULL,
                             facet_nrow = NULL,
                             gray_CN = NULL,
                             add_het = TRUE,
                             hmm_dosage_calls = NULL,
                             interactive = FALSE) {

  stopifnot(inherits(hmm_CN, "hmm_CN"))
  cnv_df <- hmm_CN$by_window

  req_cols <- c("Sample", "Chr", "Start", "End", "CN_call", "post_max")
  missing_cols <- setdiff(req_cols, colnames(cnv_df))
  if (length(missing_cols) > 0) {
    stop("hmm_CN$by_window is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  # CN_reliability combines HMM confidence and BAF support; fall back to post_max for old objects
  if (!"CN_reliability" %in% names(cnv_df)) cnv_df$CN_reliability <- cnv_df$post_max

  # If samples_to_plot is NULL, use the first sample found
  if (is.null(samples_to_plot) || length(samples_to_plot) == 0) {
    samples_to_plot <- unique(cnv_df$Sample)[1]
  }

  # ---- local helper: dependency-free, dark, high-contrast palette with black baseline ----
  make_cn_palette_dark <- function(cn_levels, baseline) {

    cn_num <- sort(unique(as.integer(as.character(cn_levels))))
    cn_chr <- as.character(cn_num)

    below <- cn_num[cn_num < baseline]
    above <- cn_num[cn_num > baseline]

    vals <- setNames(rep(NA_character_, length(cn_num)), cn_chr)

    # baseline: gray
    if (as.character(baseline) %in% cn_chr) {
      vals[as.character(baseline)] <- "#888888"  # gray
    }

    # predefined dark ramps (long enough for typical CN ranges; saturates beyond)
    blues <- c(
      "#08306b", "#2171b5", "#6baed6"  # dark, medium, light blue
    )
    reds <- c(
      "#67000d", "#cb181d", "#fc9272"  # dark, medium, light red
    )
    # If more levels, interpolate with greater color jumps
    if (length(below) > length(blues)) {
      blues <- colorRampPalette(blues, space = "Lab")(length(below))
    } else {
      blues <- blues[seq_len(length(below))]
    }
    if (length(above) > length(reds)) {
      reds <- colorRampPalette(reds, space = "Lab")(length(above))
    } else {
      reds <- reds[seq_len(length(above))]
    }
    if (length(below) > 0) {
      vals[as.character(below)] <- blues
    }
    if (length(above) > 0) {
      vals[as.character(above)] <- reds
    }
    vals
  }

  # ---- filter data ----
  plot_df <- cnv_df |>
    filter(.data$Sample %in% samples_to_plot)

  if (!is.null(chromosomes)) {
    plot_df <- plot_df |>
      filter(.data$Chr %in% chromosomes)
  }

  if (nrow(plot_df) == 0) {
    stop("No rows left after filtering. Check samples_to_plot and chromosomes.")
  }

  # order samples top-to-bottom in the user-provided order (first listed = top)
  samples_present <- intersect(samples_to_plot, unique(plot_df$Sample))
  plot_df$Sample <- factor(plot_df$Sample, levels = rev(samples_present))

  # CN_call as integer -> ordered factor
  plot_df$CN_call <- as.integer(as.character(plot_df$CN_call))

  # baseline CN = user-supplied gray_CN, or most frequent CN weighted by window length
  if (!is.null(gray_CN)) {
    baseline_cn <- as.integer(gray_CN)
  } else {
    plot_df$w <- pmax(0, plot_df$End - plot_df$Start)
    cn_totals <- plot_df |>
      group_by(.data$CN_call) |>
      summarise(total_bp = sum(.data$w, na.rm = TRUE), .groups = "drop") |>
      arrange(desc(.data$total_bp))
    baseline_cn <- cn_totals$CN_call[1]
  }

  cn_levels <- sort(unique(plot_df$CN_call))
  plot_df$CN_call <- factor(plot_df$CN_call, levels = cn_levels)

  # palette
  cn_pal <- make_cn_palette_dark(levels(plot_df$CN_call), baseline_cn)

  # Order chromosomes by extracting the numeric part and sorting accordingly
  chr_levels <- unique(plot_df$Chr)
  chr_nums <- suppressWarnings(as.numeric(gsub("[^0-9]+", "", chr_levels)))
  chr_order <- chr_levels[order(chr_nums, chr_levels)]
  plot_df$Chr <- factor(plot_df$Chr, levels = chr_order)

  # ---- build plot ----
  n_facets <- length(unique(plot_df$Chr))
  ncol_final <- facet_ncol
  nrow_final <- facet_nrow
  if (!is.null(facet_nrow) && is.null(facet_ncol)) {
    ncol_final <- ceiling(n_facets / facet_nrow)
  }
  if (!is.null(facet_ncol) && is.null(facet_nrow)) {
    nrow_final <- NULL
  }
  p_main <- ggplot(plot_df) +
    geom_segment(
      aes(
        x = .data$Start, xend = .data$End,
        y = .data$Sample, yend = .data$Sample,
        color = .data$CN_call, alpha = .data$CN_reliability
      ),
      linewidth = 4,
      lineend = "butt"
    ) +
    facet_wrap(~ Chr, scales = "free_x", ncol = ncol_final, nrow = nrow_final) +
    scale_x_continuous(
      labels = function(x) format(x / 1e6, trim = TRUE),
      name = "Position (Mb)"
    ) +
    scale_color_manual(
      values = cn_pal,
      drop = FALSE,
      name = "CN"
    ) +
    scale_alpha(range = c(0.35, 1), guide = "none") +
    labs(y = NULL) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = NA),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  # ---- compute het_df (shared by both interactive and static paths) ----
  het_df <- NULL
  if (add_het && !is.null(hmm_dosage_calls)) {
    if (!inherits(hmm_dosage_calls, "hmm_dosage_calls")) {
      stop("hmm_dosage_calls must be an object of class 'hmm_dosage_calls'")
    }
    het_raw <- hmm_dosage_calls |>
      filter(.data$SampleName %in% samples_present)
    if (!is.null(chromosomes)) {
      het_raw <- het_raw |> filter(.data$Chr %in% chromosomes)
    }
    het_df <- het_raw |>
      group_by(.data$SampleName) |>
      summarise(
        het = sum(!is.na(.data$dosage) & .data$dosage != 0 & .data$dosage != 1,
                  na.rm = TRUE) / max(sum(!is.na(.data$dosage)), 1L),
        .groups = "drop"
      )
    het_df$SampleName <- factor(het_df$SampleName, levels = levels(plot_df$Sample))
    het_df$label <- "Het."
  }

  # ---- interactive (plotly) output ----
  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for interactive = TRUE. Install with: install.packages('plotly')")
    }

    # Build per-segment tooltip lookup for het
    het_lookup <- if (!is.null(het_df)) {
      setNames(round(het_df$het, 3), as.character(het_df$SampleName))
    } else NULL

    plot_df$Mid <- (as.numeric(plot_df$Start) + as.numeric(plot_df$End)) / 2
    plot_df$tooltip_text <- paste0(
      "Sample: ", plot_df$Sample,
      "<br>Chr: ", plot_df$Chr,
      "<br>Start: ", format(as.numeric(plot_df$Start), big.mark = ",", scientific = FALSE),
      "<br>End: ", format(as.numeric(plot_df$End), big.mark = ",", scientific = FALSE),
      "<br>CN: ", plot_df$CN_call,
      "<br>post_max: ", round(plot_df$post_max, 3),
      "<br>CN_reliability: ", round(plot_df$CN_reliability, 3),
      "<br>w_baf: ", if ("w_baf" %in% names(plot_df)) round(plot_df$w_baf, 3) else "NA",
      "<br>w_het: ", round(plot_df$n_het/plot_df$n_snps, 3)
    )

    # map CN colors to hex per row
    plot_df$seg_color <- cn_pal[as.character(plot_df$CN_call)]

    chrs <- levels(plot_df$Chr)[levels(plot_df$Chr) %in% unique(as.character(plot_df$Chr))]
    n_chr <- length(chrs)

    # one subplot column per chromosome.
    # Draw one add_segments trace per CN level per chromosome so every CN level
    # gets a proper legend entry (shown only on the first chr where that CN appears).
    # Use only observed CN levels to avoid NA in showlegend.
    cn_levels_all <- as.character(sort(unique(plot_df$CN_call)))

    # For each CN level, find the first chromosome where it has data
    first_chr_for_cn <- vapply(cn_levels_all, function(cn) {
      found <- Filter(function(ch) {
        any(as.character(plot_df$CN_call[plot_df$Chr == ch]) == cn)
      }, chrs)
      if (length(found) == 0L) chrs[1] else found[1]
    }, character(1))

    chr_figs <- lapply(seq_along(chrs), function(i) {
      ch <- chrs[i]
      d <- plot_df[plot_df$Chr == ch, ]
      fig <- plotly::plot_ly()

      # One trace per CN level — guarantees a legend swatch for every level
      for (cn in cn_levels_all) {
        dsub <- d[as.character(d$CN_call) == as.character(cn), ]
        if (nrow(dsub) == 0) next
        fig <- fig |>
          plotly::add_segments(
            x    = as.numeric(dsub$Start),
            xend = as.numeric(dsub$End),
            y    = as.character(dsub$Sample),
            yend = as.character(dsub$Sample),
            line = list(color = cn_pal[[as.character(cn)]], width = 10),
            name = as.character(cn),
            legendgroup = as.character(cn),
            showlegend = (ch == first_chr_for_cn[[as.character(cn)]]),
            opacity = 1,
            hoverinfo = "none"
          )
      }

      # Invisible tooltip anchors at segment midpoints
      sample_order <- rev(levels(plot_df$Sample))
      fig <- fig |>
        plotly::add_markers(
          x    = d$Mid,
          y    = as.character(d$Sample),
          text = d$tooltip_text,
          hoverinfo = "text",
          marker = list(color = "rgba(0,0,0,0)", size = 6, line = list(width = 0)),
          showlegend = FALSE
        ) |>
        plotly::layout(
          xaxis = list(title = ch, tickformat = ".2s", tickangle = 270),
          yaxis = list(
            title = "",
            categoryarray = sample_order,
            categoryorder = "array",
            showticklabels = (i == 1)
          )
        )
      fig
    })

    fig_main <- plotly::subplot(chr_figs,
                                nrows = 1,
                                shareY = TRUE,
                                margin = 0.002)

    if (!is.null(het_df)) {
      het_df$tooltip_text <- paste0(
        "Sample: ", het_df$SampleName,
        "<br>Heterozygosity: ", round(het_df$het, 3)
      )

      sample_order <- rev(levels(plot_df$Sample))
      fig_het <- plotly::plot_ly() |>
        plotly::add_markers(
          x    = het_df$label,
          y    = as.character(het_df$SampleName),
          text = het_df$tooltip_text,
          hoverinfo = "text",
          marker = list(
            symbol    = "square",
            size      = 20,
            color     = het_df$het,
            colorscale = list(list(0, "#2166ac"), list(1, "#d6604d")),
            cmin      = 0,
            cmax      = 1,
            showscale = TRUE,
            colorbar  = list(title = "Het.", len = 0.4, y = 0.5),
            line      = list(width = 0)
          ),
          showlegend = FALSE
        ) |>
        plotly::layout(
          xaxis = list(title = ""),
          yaxis = list(
            title = "",
            categoryarray = sample_order,
            categoryorder = "array",
            showticklabels = FALSE
          )
        )
      return(plotly::subplot(fig_main, fig_het,
                             nrows = 1, widths = c(0.92, 0.08),
                             shareY = FALSE, margin = 0.0001))
    }
    return(fig_main)
  }

  # ---- static (ggplot) output ----
  if (!is.null(het_df)) {
    het_df$facet_dummy <- " "  # blank strip to match p_main chromosome strip height

    p_het <- ggplot(het_df, aes(x = .data$label, y = .data$SampleName, fill = .data$het)) +
      geom_tile(color = "white") +
      scale_fill_gradient(
        low = "#2166ac", high = "#d6604d",
        name = "Heterozygosity",
        limits = c(0, 1)
      ) +
      facet_wrap(~ facet_dummy) +
      theme_bw() +
      labs(x = "", y = NULL) +
      theme(
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        panel.grid       = element_blank(),
        legend.position  = "top",
        axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.background = element_rect(fill = "grey95"),
        strip.text       = element_text(colour = "grey95")  # invisible text, keeps height
      )

    p_main <- p_main + theme(legend.position = "top")

    leg_main <- get_legend(p_main)
    leg_het  <- get_legend(p_het)

    plots_row <- ggarrange(
      p_main + theme(legend.position = "none"),
      p_het  + theme(legend.position = "none"),
      ncol = 2, widths = c(15, 1), align = "hv"
    )
    legends_row <- ggarrange(leg_main, leg_het, ncol = 2, widths = c(7, 4))
    return(ggarrange(legends_row, plots_row, nrow = 2, heights = c(2, 20)))
  }

  return(p_main)
}


#' Karyotype-style copy-number overview across all samples
#'
#' One tile per sample × chromosome: samples in rows ordered by sample-level
#' CN (ascending), chromosomes in columns.  Within each CN block samples are
#' further sub-ordered as euploid → aneuploid → segmental.  Chromosomes whose
#' windows carry mixed CN (segmental aneuploidy) are shown in purple; others
#' use the same blue/gray/red palette as \code{\link{compare_cn_track}}.
#' Strong black lines separate CN blocks; medium dashed lines separate
#' sub-blocks.  Sub-block labels appear to the right of the matrix.
#'
#' @param hmm_CN     An object of class \code{hmm_CN}.
#' @param chromosomes  Optional character vector of chromosomes to include.
#' @param samples_to_plot  Optional character vector of sample IDs to include.
#' @param gray_CN    Integer or \code{NULL}. Baseline CN shown in gray; auto-
#'   detected as the most window-length-weighted CN when \code{NULL}.
#' @param filter_type  Optional character vector. Keep only samples whose
#'   ploidy category matches. Accepted values: \code{"euploid"},
#'   \code{"aneuploid"}, \code{"segmental"}. \code{NULL} (default) keeps all.
#' @param filter_cn  Optional integer vector. Keep only samples whose
#'   sample-level CN is in this vector (e.g. \code{c(2L, 4L)}). \code{NULL}
#'   (default) keeps all.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom dplyr group_by summarise mutate left_join arrange filter lag desc
#'   n_distinct
#' @importFrom ggplot2 ggplot aes geom_tile geom_hline geom_text scale_fill_manual
#'   scale_x_continuous scale_y_continuous coord_cartesian expansion labs
#'   theme_bw theme element_blank element_text element_line margin guide_legend
#' @export
compare_cn_track_summary <- function(hmm_CN,
                             chromosomes     = NULL,
                             samples_to_plot = NULL,
                             gray_CN         = NULL,
                             filter_type     = NULL,
                             filter_cn       = NULL) {

  stopifnot(inherits(hmm_CN, "hmm_CN"))
  cnv_df <- hmm_CN$by_window

  if (!is.null(samples_to_plot))
    cnv_df <- cnv_df[cnv_df$Sample %in% samples_to_plot, ]
  if (!is.null(chromosomes))
    cnv_df <- cnv_df[cnv_df$Chr %in% chromosomes, ]
  if (nrow(cnv_df) == 0) stop("No data remaining after filtering.")

  # NA-safe mode: ignores NA values; returns NA_integer_ when all values are NA
  .mode <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) return(NA_integer_)
    tt <- table(x); as.integer(names(tt)[which.max(tt)])
  }

  # ---- chromosome order ----
  chr_lvls <- unique(cnv_df$Chr)
  chr_lvls <- chr_lvls[order(
    suppressWarnings(as.numeric(gsub("[^0-9]+", "", chr_lvls))), chr_lvls)]
  n_chr <- length(chr_lvls)

  # ---- baseline CN ----
  cnv_df$w_len <- pmax(0, as.numeric(cnv_df$End) - as.numeric(cnv_df$Start))
  if (!is.null(gray_CN)) {
    baseline_cn <- as.integer(gray_CN)
  } else {
    bdf <- cnv_df %>%
      group_by(CN_call) %>%
      summarise(tot = sum(w_len), .groups = "drop") %>%
      arrange(desc(tot))
    baseline_cn <- bdf$CN_call[1]
  }

  # ---- CN colour palette ----
  all_cns <- sort(unique(as.integer(cnv_df$CN_call)))
  below   <- all_cns[all_cns < baseline_cn]
  above   <- all_cns[all_cns > baseline_cn]
  blues   <- c("#08306b", "#2171b5", "#6baed6")
  reds    <- c("#67000d", "#cb181d", "#fc9272")
  if (length(below) > 3) blues <- colorRampPalette(blues, space = "Lab")(length(below))
  else if (length(below) > 0) blues <- blues[seq_len(length(below))]
  if (length(above) > 3) reds  <- colorRampPalette(reds,  space = "Lab")(length(above))
  else if (length(above) > 0) reds  <- reds[seq_len(length(above))]
  cn_pal <- c(
    if (length(below) > 0) setNames(blues, as.character(below)),
    setNames("#888888", as.character(baseline_cn)),
    if (length(above) > 0) setNames(reds,  as.character(above))
  )
  mixed_col <- "#9B59B6"

  # ---- per-sample and per-cell classification ----
  sample_cn_df <- cnv_df %>%
    group_by(Sample) %>%
    summarise(sample_cn = .mode(CN_call), .groups = "drop") %>%
    filter(!is.na(sample_cn))

  cell_df <- cnv_df %>%
    group_by(Sample, Chr) %>%
    summarise(
      n_cn   = n_distinct(CN_call, na.rm = TRUE),
      dom_cn = .mode(CN_call),
      .groups = "drop"
    ) %>%
    filter(!is.na(dom_cn)) %>%
    left_join(sample_cn_df, by = "Sample") %>%
    mutate(
      cell_type = ifelse(n_cn > 1, "mixed",
                  ifelse(dom_cn == sample_cn, "euploid_chr", "aneuploid_chr")),
      cn_label  = ifelse(cell_type == "mixed", "Mixed", paste0("CN=", dom_cn))
    )

  # ---- per-sample ploidy category ----
  sample_class <- cell_df %>%
    group_by(Sample, sample_cn) %>%
    summarise(
      has_mixed     = any(cell_type == "mixed"),
      has_aneuploid = any(cell_type == "aneuploid_chr"),
      .groups = "drop"
    ) %>%
    mutate(
      sample_type = ifelse(has_mixed, "segmental",
                    ifelse(has_aneuploid, "aneuploid", "euploid")),
      sample_type = factor(sample_type, levels = c("euploid", "aneuploid", "segmental"))
    )

  # ---- apply type / CN filters before ordering ----
  if (!is.null(filter_type)) {
    valid_types <- c("euploid", "aneuploid", "segmental")
    bad <- setdiff(filter_type, valid_types)
    if (length(bad)) stop("filter_type must be one or more of: ", paste(valid_types, collapse = ", "))
    sample_class <- sample_class[as.character(sample_class$sample_type) %in% filter_type, ]
  }
  if (!is.null(filter_cn)) {
    sample_class <- sample_class[sample_class$sample_cn %in% as.integer(filter_cn), ]
  }
  if (nrow(sample_class) == 0)
    stop("No samples remain after filtering. Check filter_type and filter_cn.")
  cell_df <- cell_df[cell_df$Sample %in% sample_class$Sample, ]

  # ---- row order: ascending CN, then euploid < aneuploid < segmental ----
  sample_ordered <- sample_class %>% arrange(sample_cn, sample_type)
  sample_order   <- sample_ordered$Sample
  n_sam          <- length(sample_order)
  # first sample (lowest CN, euploid) → top → highest y value
  y_map <- setNames(rev(seq_len(n_sam)), sample_order)

  # ---- plot data ----
  cn_label_levels <- c(paste0("CN=", all_cns), "Mixed")
  plot_data <- cell_df %>%
    mutate(
      y_pos    = y_map[as.character(Sample)],
      x_pos    = match(Chr, chr_lvls),
      cn_label = factor(cn_label, levels = cn_label_levels)
    )

  fill_vals <- c(setNames(unname(cn_pal), paste0("CN=", all_cns)), Mixed = mixed_col)
  fill_vals <- fill_vals[names(fill_vals) %in% levels(plot_data$cn_label)]

  # ---- line positions ----
  row_meta <- sample_ordered %>%
    mutate(y_pos = y_map[as.character(Sample)]) %>%
    arrange(desc(y_pos)) %>%                     # top → bottom
    mutate(prev_cn = lag(sample_cn), prev_type = lag(sample_type))

  cn_lines <- row_meta %>%
    filter(!is.na(prev_cn) & sample_cn != prev_cn) %>%
    mutate(line_y = y_pos + 0.5)

  sub_lines <- row_meta %>%
    filter(!is.na(prev_type) & sample_type != prev_type & sample_cn == prev_cn) %>%
    mutate(line_y = y_pos + 0.5)

  # ---- right-side sub-block labels ----
  label_df <- row_meta %>%
    group_by(sample_cn, sample_type) %>%
    summarise(mid_y = mean(y_pos), .groups = "drop") %>%
    mutate(
      lbl = paste0(
        ifelse(sample_cn == 1L, "1x or inbred", paste0(sample_cn, "x")),
        " - ",
        ifelse(sample_type == "euploid", "euploid",
        ifelse(sample_type == "aneuploid", "aneuploid",
               "segmental aneuploidy"))
      )
    )

  # ---- assemble plot ----
  ggplot() +
    geom_tile(
      data = plot_data,
      aes(x = x_pos, y = y_pos, fill = cn_label),
      color = "white", linewidth = 0.3
    ) +
    scale_fill_manual(values = fill_vals, name = "CN call",
                      guide = guide_legend(nrow = 1)) +
    geom_hline(                                  # strong: between CN blocks
      data = cn_lines, aes(yintercept = line_y),
      color = "black", linewidth = 1.4, inherit.aes = FALSE
    ) +
    geom_hline(                                  # medium: between sub-blocks
      data = sub_lines, aes(yintercept = line_y),
      color = "black", linewidth = 0.55, linetype = "dashed", inherit.aes = FALSE
    ) +
    geom_text(
      data = label_df, aes(x = n_chr + 0.65, y = mid_y, label = lbl),
      hjust = 0, size = 2.8, inherit.aes = FALSE
    ) +
    scale_x_continuous(
      breaks = seq_len(n_chr), labels = chr_lvls,
      expand = expansion(add = c(0, 0.5)), position = "top"
    ) +
    scale_y_continuous(
      breaks = unname(y_map), labels = names(y_map),
      expand = expansion(add = c(0.5, 0.5))
    ) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_blank(),
      axis.ticks      = element_line(colour = "grey70"),
      axis.text.x     = element_text(angle = 45, hjust = 0, size = 8),
      axis.text.y     = element_text(size = 7),
      legend.position = "top", 
      plot.margin     = margin(5, 140, 5, 5)
    )
}
