# filter.R
filter_rng <- function(.data, dots) {
  overscope <- overscope_ranges(.data)

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
#' @param ...  valid logical predictates to subset .data by. These
#' are determined by variables in \code{.data}. If more than
#' one condition is supplied, the conditions are combined with \code{&}. Only
#' rows where the condition evaluates to \code{TRUE} are kept.
#
#' @details  For any Ranges objects
#' \code{filter} can act on all core components of the class including start, end,
#' width (for IRanges) or seqnames and strand (for GRanges) in addition to
#' metadata columns. If the Ranges object is grouped, \code{filter} will act
#' seperately on each  parition of the data.
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
#' df <- data.frame(start = 1:10,
#'                  width = 5,
#'                  seqnames = "seq1",
#'                  strand = sample(c("+", "-", "*"), 10, replace = TRUE),
#'                  gc = runif(10))
#'
#' rng <- as_granges(df)
#'
#' filter(rng, strand == "+")
#' filter(rng, gc > 0.5)
#'
#' # multiple criteria
#' filter(rng, strand == "+" | start > 5)
#' filter(rng, strand == "+" & start > 5)
#'
#' # multiple conditions are the same as and
#' filter(rng, strand == "+", start > 5)
#'
#' # grouping acts on each subset of the data
#' rng %>%
#'   group_by(strand) %>%
#'   filter(gc > 0.5)
#'
#' @export
filter.GenomicRanges <- function(.data, ...) {
  dots <- quos(...)
  filter_rng(.data, dots)
}

#' @method filter Ranges
#' @export
filter.Ranges <- function(.data, ...) {
  dots <- quos(...)
  filter_rng(.data, dots)
}

#' @method filter GRangesGrouped
#' @export
filter.GRangesGrouped <- function(.data, ...) {
    dots <- quos(...)
    split_ranges <- split_groups(.data)
    new("GRangesGrouped",
        unlist(endoapply(split_ranges, filter_rng, dots), use.names = FALSE),
        groups = groups(.data)
    )

}

#' @method filter IRangesGrouped
#' @export
filter.IRangesGrouped <- function(.data, ...) {
  dots <- quos(...)
  split_ranges <- split_groups(.data)
  new("IRangesGrouped",
      unlist(endoapply(split_ranges, filter_rng, dots), use.names = FALSE),
      groups = groups(.data)
  )
}
