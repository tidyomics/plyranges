#' Shift all coordinates in a genomic interval left or right, upstream or downstream
#'
#' @param x a Ranges object .
#' @param shift the amount to move the genomic interval in the Ranges object by.
#' Either a non-negative integer vector of length 1 or an integer vector the same length as x.
#' @seealso \code{\link[IRanges]{shift}}
#' @rdname shift-ranges
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
  x[strand(x) == "+"] <- shift_left(x[strand(x) == "+"], shift)
  x[strand(x) == "-"] <- shift_right(x[strand(x) == "-"], shift)
  x
}

#' @rdname shift-ranges
#' @export
shift_downstream <- function(x, shift = 0L) {
  x[strand(x) == "+"] <- shift_right(x[strand(x) == "+"], shift)
  x[strand(x) == "-"] <- shift_left(x[strand(x) == "-"], shift)
  x
}
