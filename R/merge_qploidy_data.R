#' Merge multiple Qploidy standardization and HMM objects into one
#'
#' Combines a list of \code{qploidy_standardization} objects and (optionally) a list
#' of \code{hmm_CN} objects into single merged objects. Duplicate sample names across
#' batches are detected and renamed with a numeric suffix (\code{_1}, \code{_2}, …)
#' with a warning.
#'
#' @param qploidy_list A list of objects of class \code{qploidy_standardization} (as
#'   returned by \code{standardize()} or \code{read_qploidy_standardization()}).
#'   Must contain at least one element.
#' @param hmm_list Optional. A list of objects of class \code{hmm_CN} (as returned by
#'   \code{hmm_estimate_CN()} or \code{hmm_estimate_CN_multi()}). When provided, its
#'   length must match \code{qploidy_list}. If \code{NULL} (default), only the merged
#'   \code{qploidy_standardization} is returned.
#'
#' @return A named list with:
#'   \describe{
#'     \item{\code{qploidy}}{A merged \code{qploidy_standardization} object. The
#'       \code{data} slot is the row-bind of all input data frames (with renamed sample
#'       names where needed). The \code{info} slot is taken from the first object (a
#'       warning is issued if any info fields differ across objects). The \code{filters}
#'       slot contains the element-wise sum of all numeric filter vectors.}
#'     \item{\code{hmm}}{A merged \code{hmm_CN} object (class \code{hmm_CN}), or
#'       \code{NULL} if \code{hmm_list} was not provided. \code{by_window} and
#'       \code{by_marker} are row-bound; post-CN columns that are absent in some
#'       objects are filled with \code{NA}. All per-sample parameter lists are
#'       collected into a single \code{params_samples} named list.}
#'   }
#'
#' @details
#' \strong{Duplicate sample detection.}
#' Sample names are gathered from \code{qploidy_standardization$data$SampleName} (and
#' from \code{hmm_CN$by_window$Sample} / \code{hmm_CN$by_marker$SampleName} if
#' \code{hmm_list} is provided). Any sample name that appears in more than one batch
#' is renamed in every occurrence: the copy from batch \emph{i} receives the suffix
#' \code{_i} (e.g.\ \code{MySample} → \code{MySample_1}, \code{MySample_2}). A
#' warning listing all affected sample names is issued.
#'
#' \strong{Filter aggregation.}
#' Only named elements present in \emph{all} filter vectors are summed; elements
#' missing from some objects are dropped with a warning.
#'
#' @examples
#' \dontrun{
#' merged <- merge_qploidy_datas(
#'   qploidy_list = list(std1, std2),
#'   hmm_list     = list(hmm1, hmm2)
#' )
#' merged$qploidy
#' merged$hmm
#' }
#'
#' @importFrom dplyr bind_rows
#' @export
merge_qploidy_datas <- function(qploidy_list, hmm_list = NULL) {

  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.list(qploidy_list) || length(qploidy_list) == 0)
    stop("'qploidy_list' must be a non-empty list.")
  for (i in seq_along(qploidy_list)) {
    if (!inherits(qploidy_list[[i]], "qploidy_standardization"))
      stop(sprintf("Element %d of 'qploidy_list' is not a 'qploidy_standardization' object.", i))
    if (!"SampleName" %in% names(as.data.frame(qploidy_list[[i]]$data)))
      stop(sprintf("Element %d of 'qploidy_list' has no 'SampleName' column in $data.", i))
  }

  if (!is.null(hmm_list)) {
    if (!is.list(hmm_list) || length(hmm_list) == 0)
      stop("'hmm_list' must be a non-empty list or NULL.")
    if (length(hmm_list) != length(qploidy_list))
      stop("'hmm_list' and 'qploidy_list' must have the same length.")
    for (i in seq_along(hmm_list)) {
      if (!inherits(hmm_list[[i]], "hmm_CN"))
        stop(sprintf("Element %d of 'hmm_list' is not an 'hmm_CN' object.", i))
    }
  }

  n <- length(qploidy_list)

  # ── Collect sample names per batch ────────────────────────────────────────
  q_samples <- lapply(qploidy_list, function(obj)
    unique(as.character(as.data.frame(obj$data)$SampleName)))

  # All sample names across batches
  all_names <- unlist(q_samples)
  dup_names <- unique(all_names[duplicated(all_names)])

  # Build per-batch rename maps (only for duplicated names)
  rename_maps <- vector("list", n)
  for (i in seq_len(n)) {
    m <- setNames(q_samples[[i]], q_samples[[i]])   # identity map
    for (s in intersect(q_samples[[i]], dup_names)) {
      m[s] <- paste0(s, "_", i)
    }
    rename_maps[[i]] <- m
  }

  if (length(dup_names) > 0) {
    warning(
      "Duplicate sample name(s) detected across batches; appending batch index suffix: ",
      paste(dup_names, collapse = ", "),
      call. = FALSE
    )
  }

  # ── Merge qploidy_standardization ─────────────────────────────────────────

  # info: use first; warn if fields differ
  merged_info <- qploidy_list[[1]]$info
  for (i in seq_len(n)[-1]) {
    other <- qploidy_list[[i]]$info
    common <- intersect(names(merged_info), names(other))
    differ <- common[merged_info[common] != other[common]]
    if (length(differ) > 0)
      warning(
        sprintf(
          "info field(s) differ between batch 1 and batch %d: %s. Using values from batch 1.",
          i, paste(differ, collapse = ", ")
        ),
        call. = FALSE
      )
  }

  # filters: element-wise sum over common names
  all_filter_names <- Reduce(intersect, lapply(qploidy_list, function(obj) names(obj$filters)))
  dropped_filter <- setdiff(
    unique(unlist(lapply(qploidy_list, function(obj) names(obj$filters)))),
    all_filter_names
  )
  if (length(dropped_filter) > 0)
    warning(
      "Filter counter(s) not present in all batches and were dropped from the merge: ",
      paste(dropped_filter, collapse = ", "),
      call. = FALSE
    )
  merged_filters <- Reduce("+",
    lapply(qploidy_list, function(obj) obj$filters[all_filter_names])
  )

  # data: rename samples, then rbind
  merged_data <- bind_rows(lapply(seq_len(n), function(i) {
    d <- as.data.frame(qploidy_list[[i]]$data)
    d$SampleName <- rename_maps[[i]][d$SampleName]
    d
  }))

  merged_qploidy <- structure(
    list(info = merged_info, filters = merged_filters, data = merged_data),
    class = "qploidy_standardization"
  )

  # ── Merge hmm_CN (optional) ───────────────────────────────────────────────
  merged_hmm <- NULL
  if (!is.null(hmm_list)) {

    # Helper: extract params_samples from either variant of hmm_CN
    get_params_samples <- function(obj, batch_idx, rmap) {
      if (!is.null(obj$params_samples)) {
        ps <- obj$params_samples
      } else {
        # single-sample object: wrap params in a one-element named list
        sample_name <- if (!is.null(obj$by_window) && "Sample" %in% names(obj$by_window))
          as.character(obj$by_window$Sample[1])
        else
          paste0("sample_", batch_idx)
        ps <- setNames(list(obj$params), sample_name)
      }
      # Rename duplicate sample keys
      names(ps) <- vapply(names(ps), function(nm) {
        if (nm %in% names(rmap)) rmap[[nm]] else nm
      }, character(1))
      ps
    }

    # Rename Sample column in by_window
    rename_by_window <- function(df, rmap) {
      if ("Sample" %in% names(df))
        df$Sample <- ifelse(df$Sample %in% names(rmap), rmap[df$Sample], df$Sample)
      df
    }

    # Rename SampleName column in by_marker
    rename_by_marker <- function(df, rmap) {
      if ("SampleName" %in% names(df))
        df$SampleName <- ifelse(df$SampleName %in% names(rmap), rmap[df$SampleName], df$SampleName)
      df
    }

    by_window_list   <- lapply(seq_len(n), function(i)
      rename_by_window(as.data.frame(hmm_list[[i]]$by_window), rename_maps[[i]]))
    by_marker_list   <- lapply(seq_len(n), function(i)
      rename_by_marker(as.data.frame(hmm_list[[i]]$by_marker), rename_maps[[i]]))
    params_list_all  <- lapply(seq_len(n), function(i)
      get_params_samples(hmm_list[[i]], i, rename_maps[[i]]))

    merged_by_window   <- bind_rows(by_window_list)
    merged_by_marker   <- bind_rows(by_marker_list)
    merged_params      <- do.call(c, params_list_all)

    rownames(merged_by_window) <- NULL
    rownames(merged_by_marker) <- NULL

    # Re-order post_CN columns to end (mirrors hmm_estimate_CN_multi)
    idx <- which(colnames(merged_by_window) == "post_max")
    if (length(idx) == 1 && idx < ncol(merged_by_window)) {
      idx1 <- order(colnames(merged_by_window)[(idx + 1):ncol(merged_by_window)])
      merged_by_window <- merged_by_window[, c(seq_len(idx), idx + idx1)]
    }

    merged_hmm <- structure(
      list(by_window = merged_by_window,
           by_marker = merged_by_marker,
           params_samples = merged_params),
      class = "hmm_CN"
    )
  }

  list(qploidy = merged_qploidy, hmm = merged_hmm)
}
