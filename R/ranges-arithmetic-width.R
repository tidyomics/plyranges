#' Set the width of a genomic interval
#'
#' @param x a Ranges object, to fix by either the start, end or center
#' of an interval use `anchor_start(x)`, `anchor_end(x)`,
#' `anchor_center(x)`. To fix by start for negative strand and
#' end by positive strand use `anchor_3p(x)`. To fix by start for positive
#' strand and end by negative strand use `anchor_end(x)`.
#' @param width the amount to alter the width of a Ranges object by. Either an
#' integer vector of length 1 or an integer vector the same length as x.
#'
#' @details By default unstranded features are treated as
#' positive strand features.
#'
#' @return a Ranges object with modified width
#' @importFrom IRanges resize
#' @seealso [IRanges::resize()]
#' @export
#' @examples
#' rng <- as_iranges(data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0)))
#' set_width(rng, width = 10)
#' set_width(anchor_start(rng), 10)
#' set_width(anchor_end(rng), 10)
#' set_width(anchor_center(rng), 10)
#' grng <- as_granges(data.frame(seqnames = "chr1",
#'                          strand = c("+", "-", "-", "+", "+", "-", "+"),
#'                          start=c(2:-1, 13:15),
#'                          width=c(0:3, 2:0)))
#' set_width(anchor_3p(grng), 10)
#' set_width(anchor_5p(grng), 10)
set_width <- function(x, width) UseMethod("set_width")

#' @export
set_width.Ranges <- function(x, width = 0L) {
  width(x) <- width
  x
}

#' @export
set_width.AnchoredIntegerRanges <- function(x, width = 0L) {
  anchor <- anchor(x)
  rng <- x@delegate
  switch(anchor,
         start = resize(rng, width, fix = "start"),
         end = resize(rng, width, fix = "end"),
         center = resize(rng, width, fix = "center")
         )
}


#' @export
set_width.AnchoredGenomicRanges <- function(x, width = 0L) {
  anchor <- anchor(x)
  rng <- x@delegate
  switch(anchor,
              start = resize(rng, width, fix = "start", ignore.strand = TRUE),
              end = resize(rng, width, fix = "end", ignore.strand = TRUE),
              center = resize(rng, width, fix = "center", ignore.strand = TRUE),
              "3p" = resize_by_strand(rng, width, "3p"),
              "5p" = resize_by_strand(rng, width, "5p")
  )
}

# anchor by strand version of resize
resize_by_strand <- function(x, width, anchor) {
  # anchor_3p == fix by start for negative strand
  if (anchor == "3p") {
    return(resize(x, width, fix = "end", ignore.strand = FALSE))
  } else {
    return(resize(x, width, fix = "start", ignore.strand = FALSE))
  }
}
