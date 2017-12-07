# ranges-arithmetic-flank.R

#' Shift all coordinates in a genomic interval left or right, upstream or downstream
#'
#' @param x a Ranges object .
#' @param width the width of the flanking range relative to the genomic interval
#' Ranges object by. Either an integer vector of length 1 or an integer vector
#' the same length as x. The width can be negative in which case the
#' flanking region is reversed.
#'
#' @details By default \code{flank_left} and \code{flank_right} will
#' ignore strandedness of any ranges, while flank upstream and
#' @seealso \code{\link[IRanges]{flank}}
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
