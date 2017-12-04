# filter.R
filter_rng <- function(.data, dots) {
  dots <- UQS(dots)
  found_n <- is_n(dots)
  dots <- check_n(dots)
  overscope <- overscope_ranges(.data, bind_data = found_n)

  r <- lapply(dots, overscope_eval_next, overscope = overscope)
  on.exit(overscope_clean(overscope))
  r <- Reduce(`&`, r)

  if (!is.logical(r)) {
    if (!is(r, "Rle")) {
      stop("Argument to filter condition must evaluate to a logical vector or logical-Rle")
    }
  }
  if (length(r) == 1) {
    if (r) {
      return(.data)
    } else {
      new(class(.data))
    }
  }
  # drop missing values in r
  .data[r & !is.na(r),]
}


#' Subset a \code{Ranges} object
#'
#' @param .data A \code{Ranges} object
#' @param ... a valid logical expression to act on .data, if multiple expression
#' are added, this is equivalent to performing an AND operation.
#
#' @description  For any Ranges objects
#' filter can act on all core components of the class including start, end,
#' width (for IRanges) or seqnames and strand (for GRanges) in addition to
#' metadata columns. If the Ranges object is grouped filter will act on each
#' parition of the data.
#'
#' @return a Ranges object
#'
#' @importFrom dplyr filter
#' @importFrom IRanges as.env
#' @importFrom S4Vectors runValue endoapply
#' @importFrom rlang enquo UQ new_overscope overscope_eval_next overscope_clean eval_bare
#' @seealso \code{\link[dplyr]{filter}}
#' @method filter GenomicRanges
#' @name filter-ranges
#' @rdname filter-ranges
#'
#' @examples
#' set.seed(100)
#' df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
#' strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10))
#' rng <- Ranges(df)
#' rng %>% filter(strand == "+")
#' rng %>% filter(gc > mean(gc))
#' rng %>% group_by(strand) %>% filter(gc > mean(gc))
#'
#' @export
filter.GenomicRanges <- function(.data, ...) {
  dots <- quos(...)
  filter_rng(.data, dots)
}

#' @rdname filter-ranges
#' @method filter Ranges
#' @export
filter.Ranges <- function(.data, ...) {
  dots <- quos(...)
  filter_rng(.data, dots)
}

#' @rdname filter-ranges
#' @method filter GRangesGrouped
#' @export
filter.GRangesGrouped <- function(.data, ...) {
    dots <- quos(...)
    split_ranges <- split_groups(.data)
    unlist(endoapply(split_ranges, filter_rng, dots), use.names = FALSE)

}

#' @rdname filter-ranges
#' @method filter IRangesGrouped
#' @export
filter.IRangesGrouped <- function(.data, ...) {
  dots <- quos(...)
  split_ranges <- split_groups(.data)
  unlist(endoapply(split_ranges, filter_rng, dots), use.names = FALSE)
}
