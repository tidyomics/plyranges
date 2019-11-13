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
#' default, unstranded features are treated as positive. When
#' using [shift_upstream()] or [shift_downstream()] when the `shift` argument is
#' indexed by the strandedness of the input ranges.
#'  
#'
#' @return a Ranges object with start and end coordinates shifted.
#'
#' @seealso [IRanges::shift()]
#' @importFrom IRanges shift
#' @importFrom BiocGenerics "%in%" which
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
  .norm_args_shift(x, shift)
  shift_l <- -1L * shift
  shift(x, shift_l)
}

#' @rdname shift-ranges
#' @export
shift_right <- function(x, shift = 0L) {
  .norm_args_shift(x, shift)
  shift(x, shift)
}

#' @rdname shift-ranges
#' @export
shift_upstream <- function(x, shift = 0L) {
  stopifnot(is(x, "GenomicRanges"))
  .norm_args_shift(x, shift)
  neg <- strand(x) == "-"
  pos <- strand(x) %in% c("+", "*")
  if (length(x) == length(shift)) {
    shift_neg <- shift[which(neg)]
    shift_pos <- shift[which(pos)]
    x[neg] <- shift_right(x[neg], shift_neg)
    x[pos] <- shift_left(x[pos], shift_pos)
  } else {
    x[neg] <- shift_right(x[neg], shift)
    x[pos] <- shift_left(x[pos], shift)
  }
  x
}

#' @rdname shift-ranges
#' @export
shift_downstream <- function(x, shift = 0L) {
  stopifnot(is(x, "GenomicRanges"))
  .norm_args_shift(x, shift)
  neg <- strand(x) == "-"
  pos <- strand(x) %in% c("+", "*")
  if (length(x) == length(shift)) {
    shift_neg <- shift[which(neg)]
    shift_pos <- shift[which(pos)]
    x[neg] <- shift_left(x[neg], shift_neg)
    x[pos] <- shift_right(x[pos], shift_pos)
  } else {
    x[neg] <- shift_left(x[neg], shift)
    x[pos] <- shift_right(x[pos], shift)
  }
  x
}

.norm_args_shift <- function(x, shift) {
  stopifnot(all(shift >= 0) && is.numeric(shift))
  if (length(shift) != 1 && length(shift) != length(x)) {
    stop("shift must be of length equal to 1 or x")
  }
}