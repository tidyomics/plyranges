# filter.R
filter_rng <- function(.data, expr) {
  expr <- UQ(expr)
  gr_env <- as.env(.data, parent.frame())
  overscope <- new_overscope(gr_env, parent.env(gr_env))

  r <- overscope_eval_next(overscope, expr)
  on.exit(overscope_clean(overscope))

  if (!is.logical(r)) {
    if (!(is(r, "Rle") & is.logical(runValue(r)))) {
      stop("expr must be evaluate to a logical vector or logical-Rle")
    }
  }
  # drop missing values in r
  .data[r & !is.na(r),]
}


#' Subset a \code{Ranges} object
#'
#' @param .data A \code{Ranges} object
#' @param expr a valid logical expression to act on .data
#
#' @description Unlike \pkg{dplyr}'s filter, the filter for a \code{Ranges}
#' is only valid for a single logical expression. For any Ranges objects
#' filter can act on all core components of the class including start, end,
#' width (for IRanges) or seqnames and strand (for GRanges) in addition to
#' metadaat columns. If the Ranges object is grouped filter will act on each
#' parition of the data.
#'
#' @return a Ranges object
#'
#' @importFrom dplyr filter
#' @importFrom IRanges as.env
#' @importFrom S4Vectors runValue endoapply
#' @importFrom rlang enquo UQ new_overscope overscope_eval_next overscope_clean eval_bare
#' @seealso \code{\link[dplyr]{filter}}
#' @method filter GRanges
#' @name filter-ranges
#' @rdname filter-ranges
#' @export
filter.GRanges <- function(.data, expr) {
  expr <- enquo(expr)
  filter_rng(.data, expr)
}

#' @rdname filter-ranges
#' @method filter IRanges
#' @export
filter.IRanges <- function(.data, expr) {

  expr <- enquo(expr)
  filter_rng(.data, expr)

}

#' @rdname filter-ranges
#' @method filter GRangesGrouped
#' @export
filter.GRangesGrouped <- function(.data, expr) {
    expr <- enquo(expr)
    split_ranges <- split_groups(.data)
    unlist(endoapply(split_ranges, filter_rng, expr), use.names = FALSE)

}

#' @rdname filter-ranges
#' @method filter IRangesGrouped
#' @export
filter.IRangesGrouped <- function(.data, expr) {
  expr <- enquo(expr)
  split_ranges <- split_groups(.data)
  unlist(endoapply(split_ranges, filter_rng, expr), use.names = FALSE)
}
