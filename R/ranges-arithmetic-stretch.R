#' Stretch a genomic interval
#'
#' @param x a Ranges object, to fix by either the start, end or center
#' of an interval use `anchor_start(x)`, `anchor_end(x)`,
#' `anchor_center(x)`. To fix by strand use `anchor_3p(x)` or
#' `anchor_5p(x)`.
#' @param extend the amount to alter the width of a Ranges object by. Either an
#' integer vector of length 1 or an integer vector the same length as x.
#' @description By default, `stretch(x)` will anchor by the center of a Ranges
#' object. This means that half of the value of `extend` will be added to
#' the end of the range and the remaining half subtracted from the start of 
#' the Range. The other anchors will leave the start/end fixed and stretch
#' the end/start respectively. 
#' @export
#' @return a Ranges object with modified start or end (or both) coordinates
#' @seealso `anchor()`, `mutate()`
#' @examples
#' rng <- as_iranges(data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0)))
#' rng2 <- stretch(anchor_center(rng), 10)
#' stretch(anchor_start(rng2), 10)
#' stretch(anchor_end(rng2), 10)
#' grng <- as_granges(data.frame(seqnames = "chr1",
#'                          strand = c("+", "-", "-", "+", "+", "-", "+"),
#'                          start=c(2:-1, 13:15),
#'                          width=c(0:3, 2:0)))
#' stretch(anchor_3p(grng), 10)
#' stretch(anchor_5p(grng), 10)
stretch <- function(x, extend) { UseMethod("stretch") }

#' @export
stretch.Ranges <- function(x, extend = 0L) {
  stopifnot(length(extend) == 1L | length(extend) == length(x))
  stretch_center(x, extend)
}

#' @export
stretch.AnchoredIntegerRanges <- function(x, extend = 0L) {
  anchor <- anchor(x)
  rng <- x@delegate
  stopifnot(length(extend) == 1L | length(extend) == length(rng))
  switch(
    anchor,
    start = mutate(rng, end = end + extend),
    end = mutate(rng, start = start - extend),
    center = stretch_center(rng, extend)
  )
}


#' @export
stretch.AnchoredGenomicRanges <- function(x, extend = 0L) {
  anchor <- anchor(x)
  rng <- x@delegate
  stopifnot(length(extend) == 1L | length(extend) == length(rng))
  
  switch(
    anchor,
    start = mutate(rng, end = end + extend),
    end = mutate(rng, start = start - extend),
    center = stretch_center(rng, extend),
    "3p" = stretch_by_strand(rng, extend, "3p"),
    "5p" = stretch_by_strand(rng, extend, "5p")
  )
}


#' @importFrom IRanges mid
stretch_center <- function(x, extend) {
  x + extend / 2
}

stretch_by_strand <- function(x, extend, anchor) {
  # anchor by 3p start is fixed for negative strand
  # and end is fixed for positive strand
  if (anchor == "3p") {
    x_pos <- x[strand(x) == "+" | strand(x) == "*"]
    start(x_pos) <- start(x_pos) - extend
    x_neg <- x[strand(x) == "-"]
    end(x_neg) <- end(x_neg) + extend
  } else {
    x_pos <- x[strand(x) == "+" | strand(x) == "*"]
    end(x_pos) <- end(x_pos) + extend
    x_neg <- x[strand(x) == "-"]
    start(x_neg) <- start(x_neg) - extend
  }
  x[strand(x) == "+" | strand(x) == "*"] <- x_pos
  x[strand(x) == "-"] <- x_neg
  return(x)
}