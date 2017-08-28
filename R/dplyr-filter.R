# filter.R
#' Subset a \code{GRanges} object
#'
#' @param .data A \code{GRanges} object
#' @param expr a valid logical expression to act on .data
#'
#' @return a GRanges object
#'
#' @description Unlike \pkg{dplyr}'s filter, the filter for \code{GRanges}
#' is only valid for a single logical expression.
#' @importFrom dplyr filter
#' @importFrom IRanges as.env runValue
#' @seealso \link[dplyr]{filter}
filter.GRanges <- function(.data, expr) {

  expr <- enexpr(expr)
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


filter.GroupedGRanges <- function(.data, expr) {

    expr <- enexpr(expr)
    groups <- groups(.data)
    gr_env <- as.env(.data, parent.frame())
    value <- lapply(groups, eval_bare, env = gr_env)
    value <- as(value, "RleList")
    split_ranges <- GRangesList(splitAsList(GRanges(.data), value, drop = TRUE))
    split_ranges
    by_group <- endoapply(split_ranges, filter.GRanges, !!expr)
    new("GroupedGRanges", unlist(by_group), groups = groups)

}
