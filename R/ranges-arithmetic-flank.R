# ranges-arithmetic-flank.R

#' Generate flanking regions
#' @description Find flanking regions  to the left or right or
#' upstream or downstream of a Ranges object.

#' @param x a Ranges object.
#' @param width the width of the flanking region relative to the ranges in
#' `x`. Either an integer vector of length 1 or an integer vector
#' the same length as x. The width can be negative in which case the
#' flanking region is reversed.
#'
#' @details The function
#' `flank_left` will create the flanking region to the left of starting
#' coordinates in `x`, while `flank_right` will create the flanking
#' region to the right of the starting coordinates in `x`. The function
#' `flank_upstream` will `flank_left` if the strand of rows in `x` is
#' not negative and will `flank_right` if the strand of rows in `x` is
#' negative. The function `flank_downstream` will `flank_right` if the strand of rows in `x` is
#' not negative and will `flank_leftt` if the strand of rows in `x` is
#' negative.
#'
#' By default `flank_left` and `flank_right` will
#' ignore strandedness of any ranges, while `flank_upstream` and
#' `flank_downstream` will take into account the strand of `x`.
#'
#' @return A Ranges object of same length as `x`.
#'
#' @seealso [IRanges::flank()][GenomicRanges::flank()]
#' @importFrom IRanges flank
#' @examples
#' gr <- as_granges(data.frame(start = 10:15,
#'                             width = 5,
#'                             seqnames = "seq1",
#'                             strand = c("+", "+", "-", "-", "+", "*")))
#' flank_left(gr, width = 5L)
#' flank_right(gr, width = 5L)
#' flank_upstream(gr, width = 5L)
#' flank_downstream(gr, width = 5L)
#' @rdname flank-ranges
#' @export
flank_left <- function(x, width = 0L) {
  stopifnot(is.numeric(width))
  flank(x, width, start = TRUE, ignore.strand = TRUE)
}

#' @rdname flank-ranges
#' @export
flank_right <- function(x, width = 0L) {
  stopifnot(is.numeric(width))
  flank(x, width, start = FALSE, ignore.strand = TRUE)
}

#' @rdname flank-ranges
#' @export
flank_upstream <- function(x, width = 0L) {
  stopifnot(is.numeric(width))
  stopifnot(is(x, "GenomicRanges"))
  flank(x, width, start = TRUE, ignore.strand = FALSE)
}

#' @rdname flank-ranges
#' @export
flank_downstream <- function(x, width = 0L) {
  stopifnot(is.numeric(width))
  stopifnot(is(x, "GenomicRanges"))
  flank(x, width, start = FALSE, ignore.strand = FALSE)
}
