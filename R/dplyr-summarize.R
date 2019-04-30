summarize_rng <- function(.data, dots) {
  overscope <- overscope_ranges(.data)
  results <- DataFrame(overscope_eval_update(overscope, dots, FALSE))
  rownames(results) <- NULL
  results
}

#' Aggregate a Ranges object
#'
#' @param .data a Ranges object
#' @param ... Name-value pairs of summary functions.
#'
#' @return a [S4Vectors::DataFrame()]
#' @seealso [dplyr::summarise()]
#' @importFrom S4Vectors rbind cbind
#' @importFrom dplyr summarise summarize
#' @rdname ranges-summarise
#' @examples
#' df <- data.frame(start = 1:10, width = 5,  seqnames = "seq1",
#' strand = sample(c("+", "-", "*"), 10, replace = TRUE), gc = runif(10))
#' rng <- as_granges(df)
#' rng %>% summarise(gc = mean(gc))
#' rng %>% group_by(strand) %>% summarise(gc = mean(gc))
#' @method summarise Ranges
#' @rdname ranges-summarise
#' @export
summarise.Ranges <- function(.data, ...) {
  dots <- set_dots_named(...)
  summarize_rng(.data, dots)
}

#' @method summarise DelegatingGenomicRanges
#' @export
summarise.DelegatingGenomicRanges <- function(.data, ...) {
  dots <- set_dots_named(...)
  delegate <- .data@delegate
  summarize_rng(delegate, dots)
}

#' @method summarise DelegatingGenomicRanges
#' @export
summarise.DelegatingIntegerRanges <- function(.data, ...) {
  dots <- set_dots_named(...)
  delegate <- .data@delegate
  summarize_rng(delegate, dots)
}

#' @importFrom rlang !!! enquos
#' @importFrom dplyr bind_cols bind_rows
#' @method summarise GroupedGenomicRanges
#' @export
summarise.GroupedGenomicRanges <- function(.data, ...) {
  dots <- set_dots_named(...)
  cbind(mcols(.data@inx, use.names = FALSE), summarize_rng(.data, dots))
}

#' @method summarise GroupedIntegerRanges
#' @export
summarise.GroupedIntegerRanges <- summarise.GroupedGenomicRanges