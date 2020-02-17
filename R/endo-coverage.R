#' Compute coverage over a Ranges object
#'
#' @param x a `Ranges` object
#' @param shift shift how much should each range in x be shifted by? (default = 0L)
#' @param width width how long should the returned coverage score be?
#' This must be either a  positive integer or NULL (default = NULL)
#' @param weight weight how much weight should be assigned to each range? Either
#' an integer or numeric vector or a column in x. (default = 1L)
#' @param ... other optional parameters to pass to coverage
#'
#' @return An expanded Ranges object with a score column corresponding to
#' the coverage value over that interval. Note that compute_coverage
#' drops metadata associated with the orginal ranges.
#' @seealso [IRanges::coverage()], [GenomicRanges::coverage()]
#' @examples
#' rng <- as_iranges(data.frame(start = 1:10, width = 5))
#' compute_coverage(rng)
#' compute_coverage(rng, shift = 14L)
#' compute_coverage(rng, width = 10L)
#' @importFrom IRanges coverage ranges
#' @importFrom S4Vectors runValue
#' @export
compute_coverage <- function(x, shift, width, weight, ...) {
  UseMethod("compute_coverage")
}

#' @export
compute_coverage.default <- function(x, shift = 0L, width = NULL, weight = 1L, ...) {
  as_ranges(coverage(x, shift, width, weight, ...))
}


setMethod("coverage", "DelegatingGenomicRanges",
          function(x, shift = 0L, width = NULL, weight = 1L, ...) {
            coverage(load_delegate(x), shift, width, weight, ...)
          })
