# filter.R
filter_rng <- function(.data, dots) {
  overscope <- overscope_ranges(.data)
  r <- overscope_eval_update(overscope, dots, bind_envir = FALSE)
  r <- Reduce(`&`, r)
  if (!(is.logical(r) || is(r, "LogicalList"))) {
    if (!(is(r, "Rle") || is(r, "RleList"))) {
      stop("Argument to filter condition must evaluate to a logical vector or logical-Rle")
    }
  }
  # drop missing values in r
  r & !is.na(r)
}

filter_grp <- function(.data, ...) {
    dots <- set_dots_unnamed(...)
    ii <- filter_rng(.data, dots)
    inx <- S4Vectors::split(
      seq_along(.data@delegate),
      .data@group_indices
    )
    inx_update <- inx[ii]
    rng <- .data@delegate[sort(unlist(inx_update))]
    group_by(rng, !!!groups(.data))
}

#' Subset a `Ranges` object
#'
#' @param .data A `Ranges` object
#' @param ...  valid logical predictates to subset .data by. These
#' are determined by variables in `.data`. If more than
#' one condition is supplied, the conditions are combined with `&`. Only
#' rows where the condition evaluates to `TRUE` are kept.
#' @param .preserve when FALSE (the default) grouping structure is recalculated, TRUE is currently not implemented.
#'
#'
#' @details  For any Ranges objects
#' `filter` can act on all core components of the class including start, end,
#' width (for IRanges) or seqnames and strand (for GRanges) in addition to
#' metadata columns. If the Ranges object is grouped, `filter` will act
#' seperately on each  parition of the data.
#'
#' @return a Ranges object
#'
#' @importFrom dplyr filter
#' @importFrom IRanges as.env
#' @importFrom S4Vectors runValue endoapply
#' @seealso [dplyr::filter()]
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
#' @method filter Ranges
#' @export
filter.Ranges <- function(.data, ..., .preserve = FALSE) {
  dots <- set_dots_unnamed(...)
  .data[filter_rng(.data, dots)]
}

#' @method filter DelegatingGenomicRanges
#' @export
filter.DelegatingGenomicRanges <- function(.data, ..., .preserve = FALSE) {
  dots <- set_dots_unnamed(...)
  delegate <- .data@delegate
  .data@delegate <- delegate[filter_rng(delegate, dots)]
  return(.data)
}

#' @method filter DelegatingIntegerRanges
#' @export
filter.DelegatingIntegerRanges <- filter.DelegatingGenomicRanges

#' @method filter GroupedGenomicRanges
#' @export
filter.GroupedGenomicRanges <- function(.data, ..., .preserve = FALSE) {
  filter_grp(.data, ...)
}

#' @method filter GroupedIntegerRanges
#' @export
filter.GroupedIntegerRanges <- filter.GroupedGenomicRanges