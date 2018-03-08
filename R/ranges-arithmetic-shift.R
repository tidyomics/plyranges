#' Shift all coordinates in a genomic interval left or right, upstream or downstream
#'
#' @param x a Ranges object .
#' @param shift the amount to move the genomic interval in the Ranges object by.
#' Either a non-negative integer vector of length 1 or an integer vector
#' the same length as x.
#'
#' @details Shifting left or right will ignore any strand information
#' in the Ranges object, while shifting upstream/downstream will shift coordinates
#' on the positive strand left/right and the negative strand right/left. By
#' default, unstranded features are treated as positive.
#'
#' @return a Ranges object with start and end coordinates shifted.
#'
#' @seealso [IRanges::shift()]
#' @importFrom IRanges shift
#' @rdname shift-ranges
#' @examples
#' ir <- as_iranges(data.frame(start = 10:15, width = 5))
#' shift_left(ir, 5L)
#' shift_right(ir, 5L)
#' gr <- as_granges(data.frame(start = 10:15,
#'                             width = 5,
#'                             seqnames = "seq1",
#'                             strand = c("+", "+", "-", "-", "+", "*")))
#' shift_upstream(gr, 5L)
#' shift_downstream(gr, 5L)
#' @export
shift_left <- function(x, shift = 0L) {
  stopifnot(all(shift > 0) && is.numeric(shift))
  shift_l <- -1L * shift
  shift(x, shift_l)
}

#' @rdname shift-ranges
#' @export
shift_right <- function(x, shift = 0L) {
  stopifnot(all(shift > 0) && is.numeric(shift))
  shift(x, shift)
}

#' @rdname shift-ranges
#' @export
shift_upstream <- function(x, shift = 0L) {
  stopifnot(is(x, "GenomicRanges"))
  x[strand(x) == "-"] <- shift_right(x[strand(x) == "-"], shift)
  x[strand(x) == "+"] <- shift_left(x[strand(x) == "+"], shift)
  x
}

#' @rdname shift-ranges
#' @export
shift_downstream <- function(x, shift = 0L) {
  stopifnot(is(x, "GenomicRanges"))
  x[strand(x) == "-"] <- shift_left(x[strand(x) == "-"], shift)
  x[strand(x) == "+"] <- shift_right(x[strand(x) == "+"], shift)
  x
}
