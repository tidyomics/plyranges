#' Ranges friendly coverage method
#'
#' @param x a \code{Ranges} object
#' @param shift shift how much should each range in x be shifted by? (default = 0L)
#' @param width width how long should the returned coverage score be?
#' This must be either a  positive integer, NA or NULL (default = NULL)
#' @param weight weight how much weight should be assigned to each range? Either
#' an integer or numeric vector or a column in x. (default = 1L)
#'
#' @return An expanded Ranges object with a score column corresponding to
#' the coverage value over that interval. Note that compute_coverage
#' drops metadata associated with the orginal ranges.
#' @seealso \link[IRanges]{coverage}
#' @importFrom IRanges coverage ranges runValue
#' @export
#' @rdname set_coverage
set_coverage <- function(x, shift, width, weight) {
  UseMethod("set_coverage")
}

#' @rdname set_coverage
#' @export
set_coverage.GenomicRanges <- function(x, shift = 0L, width = NULL, weight = 1L) {
  GRanges(coverage(x, shift, width, weight, method = "auto"))
}

set_coverage.Ranges <- function(x, shift = 0L, width = NULL, weight = 1L) {
  cvg <- coverage(x, shift, width, weight, method)
  rng <- ranges(cvg)
  mcols(rng)[["score"]] <- runValue(cvg)
  rng
}
