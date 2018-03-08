select_rng <- function(.data, .drop_ranges, dots) {

  if (.drop_ranges) {
    var_names <- ranges_vars(.data)
    selected_vars <- tidyselect::vars_select(var_names, UQS(dots))
  } else {
    var_names <- names(mcols(.data))
    if (length(var_names) == 0L) {
      stop("No metadata columns to select", call. = FALSE)
    }

    if (is(.data, "Ranges")) {
      exclude <- c("start", "end", "width")
    } else {
      exclude <- c("start", "end", "width", "strand", "seqnames")
    }

    selected_vars <- tidyselect::vars_select(var_names, UQS(dots),
                                             .exclude = exclude)
  }

  selected_vars <- syms(selected_vars)
  rng_os <- overscope_ranges(.data)
  on.exit(overscope_clean(rng_os))
  rng_list <- lapply(selected_vars, overscope_eval_next, overscope = rng_os)

  if (.drop_ranges) {
    return(do.call("DataFrame", rng_list))
  } else {
    if (length(rng_list) > 0) {
      mcols(.data) <- do.call("DataFrame", rng_list)
    } else {
      mcols(.data) <- NULL
    }

    return(.data)
  }
}


#' Select metadata columns of the Ranges object by name or position
#'
#' @param .data a `Ranges` object
#' @param ... One or more metadata column names.
#' @param .drop_ranges If TRUE select will always return a tibble. In this
#' case, you may select columns that form the core part of the Ranges object.
#' @details Note that by default select only acts on the metadata columns (and
#' will therefore return a Ranges object) if a core component of a Ranges is dropped or selected
#' without the other required components (this includes the seqnames, strand, start, end,
#' width names), then select will throw an error unless .drop_ranges is set to TRUE.
#' @return a Ranges object or a tibble
#' @seealso [dplyr::select()]
#' @importFrom dplyr select
#' @importFrom tidyselect vars_select
#' @rdname ranges-select
#' @method select GenomicRanges
#' @examples
#' df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
#' strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10), counts = rpois(10, 2))
#' rng <- as_granges(df)
#' select(rng, -gc)
#' select(rng, gc)
#' select(rng, counts, gc)
#' select(rng, 2:1)
#' select(rng, seqnames, strand, .drop_ranges = TRUE)
#' @export
select.GenomicRanges <- function(.data, ..., .drop_ranges = FALSE) {

  dots <- quos(...)
  # no selection? - return .data
  if (length(dots) == 0) {
    return(.data)
  }

  # if the quo is named return an error
  if (length(names(dots)) == 0) {
    stop("select does not support renaming variables", call. = FALSE)
  }

  select_rng(.data, .drop_ranges, dots)
}

#' @rdname ranges-select
#' @method select Ranges
#' @export
select.Ranges <- function(.data, ..., .drop_ranges = FALSE) {
  dots <- quos(...)

  # no selection? - return .data
  if (length(dots) == 0) {
    return(.data)
  }

  # if the quo is named return an error
  if (length(names(dots)) == 0) {
    stop("select does not support renaming variables", call. = FALSE)
  }

  select_rng(.data, .drop_ranges,  dots)
}
