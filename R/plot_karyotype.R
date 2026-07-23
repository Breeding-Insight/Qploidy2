utils::globalVariables(c(
    "segment_id", "Position_Mb", "color_value", "copy_num",
    "x_pos", "width", "y_start", "y_end", "grey_type",
    "x", "x_start", "x_end", "CN_call"
))

##' Plot a karyotype for a single sample
##'
##' Draws a chromosome-copy karyotype for one sample. Each chromosome is shown
##' as a vertical facet; within each facet, one column is drawn per chromosome
##' copy as estimated by the HMM copy-number call (\code{CN_call}). Consecutive
##' regions with the same \code{CN_call} are grouped into segments; if the CN
##' changes along a chromosome, each segment shows only as many columns as its
##' copy number. Marker positions are overlaid as horizontal tick marks, coloured
##' by the chosen \code{color_by} variable. Optionally, karyotype notation from
##' \code{\link{karyotype_notation}} can be supplied to annotate facet strip
##' labels with chromosome-level copy number and segment deviations, and to
##' display the sample ploidy in the subtitle.
##'
##' @param df A \code{data.frame} with at least the columns: \code{SampleName},
##'   \code{Chr}, \code{Position} (base-pair), \code{CN_call} (integer copy
##'   number per marker), and whichever column is requested via \code{color_by}.
##'   Typically the \code{by_marker} slot from an \code{hmm_CN} object after
##'   joining with BAF/dosage information.
##' @param sample_name Character scalar. The value of \code{SampleName} to
##'   display. An error is thrown if no rows match.
##' @param color_by Character scalar or \code{NULL}. Column used to colour the
##'   marker tick marks. Pass \code{NULL} or \code{NA} to draw all tick marks in
##'   black. Otherwise one of:
##'   \describe{
##'     \item{\code{"dosage"}}{Allele dosage (0 to ploidy).}
##'     \item{\code{"post_max_CN"}}{Posterior probability of the most likely copy-number state.}
##'     \item{\code{"post_max_dosage"}}{Posterior probability of the most likely dosage associated with the most likely
##'       copy-number state.}
##'     \item{\code{"w_baf"}}{B-allele frequency weight applied on the HMM according to the number of heterozygous markers in the window.}
##'   }
##' @param nrow Integer. Number of rows passed to \code{\link[ggplot2]{facet_wrap}}
##'   when arranging chromosome panels. Default \code{1} places all chromosomes
##'   in a single row.
##' @param notation Logical. If \code{TRUE} (default), \code{\link{karyotype_notation}}
##'   is run internally on \code{df} for \code{sample_name} and used to:
##'   \itemize{
##'     \item Annotate each chromosome strip with the chromosome-level copy
##'       number (\code{CN: X}) and any deviating segments
##'       (\code{g.start_endxXX}), one segment per line.
##'     \item Add \code{Ploidy: X} to the plot subtitle.
##'   }
##'   Set to \code{FALSE} to omit the annotation.
##'
##' @importFrom ggplot2 ggplot aes geom_rect geom_segment geom_point scale_color_viridis_c scale_color_viridis_d scale_color_identity scale_fill_manual scale_y_reverse scale_x_continuous facet_wrap labeller label_value labs theme_minimal theme element_blank element_text unit margin guide_legend
##' @importFrom scales squish
##' @importFrom dplyr filter arrange mutate group_by summarise ungroup left_join distinct
##' @importFrom tidyr unnest
##'
##' @return A \code{\link[ggplot2]{ggplot}} object. The y-axis is reversed so
##'   chromosome position increases downward (top = chromosome start). The x-axis
##'   shows copy indices (1, 2, ..., CN). Panels are faceted by chromosome in
##'   genomic order.
##'
##' @examples
##' \dontrun{
##' # Default: notation annotated automatically
##' plot_karyotype(hmm_CN, sample_name = "S1", color_by = "dosage")
##'
##' # No coloring, no notation
##' plot_karyotype(hmm_CN, sample_name = "S1", color_by = NULL, notation = FALSE)
##'
##' # Two rows of chromosome panels
##' plot_karyotype(hmm_CN, sample_name = "S1", color_by = "w_baf", nrow = 2)
##' }
##'
##' @export
plot_karyotype <- function(df, sample_name,
                           color_by = NULL,
                           nrow = 1,
                           notation = TRUE) {
    if (is.null(color_by) || (length(color_by) == 1 && is.na(color_by))) {
        color_by <- NULL
    } else {
        color_by <- match.arg(color_by,
                              choices = c("dosage", "post_max_CN", "post_max_dosage", "w_baf"))
    }

    # --- 0. Input checks -----------------------------------------------------
    if (inherits(df, "hmm_CN")) {
        df <- df$by_marker
    } else if (is.list(df) && !is.data.frame(df) &&
               all(vapply(df, inherits, logical(1L), "hmm_CN"))) {
        df <- df[[sample_name]]$by_marker
    }

    # Normalise column names: hmm_CN$by_marker uses "post_max" for the
    # posterior probability of the CN call; alias it so color_by = "post_max_CN"
    # works whether the input came from hmm_CN or call_hmm_dosages().
    if ("post_max" %in% names(df) && !"post_max_CN" %in% names(df))
        df[["post_max_CN"]] <- df[["post_max"]]

    if (!is.data.frame(df))
        stop("`df` must be a data.frame or an hmm_CN object.")

    required_cols <- c("SampleName", "Chr", "Position", "CN_call")
    if (!is.null(color_by)) required_cols <- c(required_cols, color_by)
    missing_cols  <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
        dosage_cols <- intersect(missing_cols, c("dosage", "post_max_dosage"))
        if (length(dosage_cols) > 0)
            stop("Column(s) ", paste(sQuote(dosage_cols), collapse = ", "),
                 " not found in `df`. ",
                 "Run `call_hmm_dosages()` on your hmm_CN object first and pass ",
                 "the result as `df`.")
        stop("The following required columns are missing from `df`: ",
             paste(missing_cols, collapse = ", "))
    }

    if (!is.character(sample_name) || length(sample_name) != 1)
        stop("`sample_name` must be a single character string.")

    if (!is.numeric(nrow) || length(nrow) != 1 || nrow < 1 || nrow != as.integer(nrow))
        stop("`nrow` must be a single positive integer.")

    if (!is.logical(notation) || length(notation) != 1)
        stop("`notation` must be TRUE or FALSE.")


    df_sample <- df %>%
        filter(SampleName == sample_name) %>%
        arrange(Chr, Position) %>%
        mutate(
            Position_Mb = Position / 1e6,
            color_value = if(!is.null(color_by)) .data[[color_by]] else NA
        )

    if (nrow(df_sample) == 0) stop(paste("Sample", sample_name, "not found."))

    # --- 2. Chromosome ordering ----------------------------------------------
    chr_levels <- unique(df_sample$Chr) %>%
        .[order(as.numeric(gsub("[^0-9]", "", .)))]
    df_sample$Chr <- factor(df_sample$Chr, levels = chr_levels)

    # --- 3. Assign segment IDs (consecutive runs of same CN_call per Chr) -----
    df_sample <- df_sample %>%
        group_by(Chr) %>%
        mutate(
            segment_id = cumsum(CN_call != lag(CN_call, default = first(CN_call)))
        ) %>%
        ungroup()

    # --- 4. Create chromosome segments (copies) per CN segment ---------------
    chr_segments <- df_sample %>%
        group_by(Chr, segment_id) %>%
        summarise(
            y_start = min(Position_Mb),
            y_end = max(Position_Mb),
            CN_call = first(CN_call),
            color_value = first(color_value),
            .groups = "drop"
        ) %>%
        rowwise() %>%
        mutate(copy_num = list(seq_len(CN_call))) %>%
        ungroup() %>%
        unnest(copy_num) %>%
        mutate(
            x_pos      = copy_num,
            width      = 0.8,
            grey_type  = "Gap (no markers)"
        )

    # --- 5. Create marker positions for each copy ---------------------------
    # Each marker appears on all copies of its CN segment
    marker_segments <- df_sample %>%
        left_join(
            chr_segments %>%
                select(Chr, segment_id, copy_num) %>%
                distinct(),
            by = c("Chr", "segment_id"),
            relationship = "many-to-many"
        ) %>%
        mutate(
            x_start = copy_num - 0.4,
            x_end   = copy_num + 0.4
        )

    # Detect whether any NA color values exist (drives NA legend entry)
    has_na_color <- !is.null(color_by) && any(is.na(marker_segments$color_value))
    na_legend_df <- if (has_na_color)
        data.frame(x = NA_real_, y = NA_real_, grey_type = "Color NA value")
    else
        NULL

    # --- 6. Color scale setup -----------------------------------------------
    if (!is.null(color_by)) {
      color_label <- switch(color_by,
                            post_max_CN      = "P(CN call)",
                            post_max_dosage  = "P(Dosage call)",
                            w_baf            = "BAF Weight",
                            dosage           = "Dosage"
      )

      is_discrete <- color_by %in% c("dosage")

      if (is_discrete) {
        # Convert to ordered factor so the discrete colour scale is used
        dosage_levels <- sort(unique(na.omit(marker_segments$color_value)))
        marker_segments <- marker_segments %>%
            mutate(color_value = factor(color_value, levels = dosage_levels))
        color_scale <- scale_color_viridis_d(
          na.value = "grey50",
          name     = color_label,
          option   = "plasma"
        )
      } else {
        mid_val <- median(marker_segments$color_value, na.rm = TRUE)
        color_scale <- scale_color_viridis_c(
          option = if(color_label == "BAF Weight") "plasma" else "D",
          direction = -1,
          limits = c(0, 1),
          oob = squish,
          name = color_label
        )
      }
    } else {
        color_label <- NULL
        color_scale <- scale_color_identity(guide = "none")
    }

    # --- 7. Notation: per-chromosome strip labels and ploidy text -----------
    if (isTRUE(notation)) {
        notation_sample <- karyotype_notation(df, sample_name = sample_name)

        if (nrow(notation_sample) == 0)
            warning("karyotype_notation() returned no rows for sample '", sample_name,
                    "'. Strip labels will not be annotated.")

        chr_labels <- setNames(
            vapply(levels(df_sample$Chr), function(chr) {
                row <- notation_sample[notation_sample$ChrID == chr, ]
                if (nrow(row) == 0) return(as.character(chr))
                cn_line <- paste0("CN: ", row$ChrCN)
                if (!is.na(row$Segment)) {
                    # Split multiple segments, strip "ChrXX:" prefix from each
                    segs       <- trimws(strsplit(row$Segment, ";", fixed = TRUE)[[1]])
                    segs_clean <- sub(paste0(chr, ":"), "", segs, fixed = TRUE)
                    paste(c(chr, cn_line, segs_clean), collapse = "\n")
                } else {
                    paste(chr, cn_line, sep = "\n")
                }
            }, character(1)),
            levels(df_sample$Chr)
        )
        facet_labeller <- labeller(Chr = chr_labels)
        ploidy_text    <- if (nrow(notation_sample) > 0)
            paste("Ploidy:", unique(notation_sample$Ploidy)[1])
        else
            NULL
    } else {
        facet_labeller <- label_value
        ploidy_text    <- NULL
    }

    # --- 6. Plot -------------------------------------------------------------
    p <- ggplot() +
        # Chromosome background rectangles (grey80 = gap / no markers)
        geom_rect(
            data = chr_segments,
            aes(
                xmin = x_pos - width / 2, xmax = x_pos + width / 2,
                ymin = y_start, ymax = y_end,
                fill = grey_type
            ),
            color = "grey30", linewidth = 0.3
        ) +
        # Invisible point that adds the "Color NA value" swatch to the fill legend
        {if(has_na_color) { geom_point(
            data  = na_legend_df,
            aes(x = x, y = y, fill = grey_type),
            shape = 22, size = 4, color = NA,
            show.legend = TRUE
        ) }} +
        scale_fill_manual(
            values = c("Gap (no markers)" = "grey80", "Color NA value" = "grey50"),
            name   = NULL,
            guide  = guide_legend(
                order        = 2,
                override.aes = list(size = 4, color = NA)
            )
        ) +
        # Marker tick marks colored by color_value (black when color_by is NULL)
        geom_segment(
            data = marker_segments,
            aes(
                x = x_start, xend = x_end,
                y = Position_Mb, yend = Position_Mb,
                color = if (!is.null(color_by)) color_value else "black"
            ),
            linewidth = 0.4
        ) +
        color_scale +
        # Reverse y-axis so chromosome starts at top
        scale_y_reverse(
            name = "Position (Mb)",
            breaks = scales::pretty_breaks(n = 8)
        ) +
        scale_x_continuous(
            name = "Copy Number",
            breaks = function(x) seq(1, max(x, na.rm = TRUE)),
            expand = c(0.1, 0.1)
        ) +
        # Facet by chromosome
        facet_wrap(~Chr, nrow = nrow, labeller = facet_labeller) +
        labs(
            title    = paste("Karyotype \u2014", sample_name),
            subtitle = paste(
                c(if (!is.null(color_label)) paste("Colored by:", color_label),
                  ploidy_text),
                collapse = "  |  "
            ),
            caption  = "Higher copy numbers are represented as parallel columns. Note that duplicated copies may in reality be arranged in tandem within the chromosome."
        ) +
        theme_minimal(base_size = 11) +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 7),
            strip.text = element_text(face = "bold", size = 9),
            plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 10, color = "grey40"),
            legend.position = "right",
            panel.spacing = unit(0.8, "lines"),
            plot.caption  = element_text(size = 7, color = "grey50", hjust = 0,
                                         margin = margin(t = 6))
        )

    return(p)
}
