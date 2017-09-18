#' Set the width of a genomic interval
#'
#' @param x a Ranges object, to fix by either the start, end or center
#' of an interval use \code{anchor_start(x)}, \code{anchor_end(x)},
#' \code{anchor_center(x)}. To leave strand fixed use either \code{anchor_3p(x)}
#' or \code{anchor_end(x)}.
#' @param width the amount to alter the width of a Ranges object by. Either an
#' integer vector of length 1 or an integer vector the same length as x.
#' @importFrom IRanges resize
#' @seealso \code{\link[IRanges]{resize}}
#' @export
#' @examples
#' rng <- Ranges(data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0)))
#' set_width(rng, width = 10)
#' set_width(anchor_start(rng), 10)
#' set_width(anchor_end(rng), 10)
#' set_width(anchor_center(rng), 10)
#' grng <- Ranges(data.frame(seqnames = "chr1",
#'                          strand = c("+", "-", "-", "+", "+", "-", "+"),
#'                          start=c(2:-1, 13:15),
#'                          width=c(0:3, 2:0)))
#' set_width(anchor_3p(grng), 10)
#' set_width(anchor_5p(grng), 10)
set_width <- function(x, width) UseMethod("set_width")

set_width.Ranges <- function(x, width = 0L) {

  anchor <- anchors(x)
  if (length(anchor) == 0) {
    width(x) <- width
    return(x)
  }

  if (any(anchor %in% c("3p", "5p"))) {
    stop("Unable to anchor by strand for IRanges", call. = FALSE)
  }

  for (i in seq_along(anchor)) {
    x <- switch(anchor[[i]],
                start = resize(x, width, fix = "start"),
                end = resize(x, width, fix = "end"),
                center = resize(x, width, fix = "center"))
  }

  metadata(x)[["anchor"]] <- NULL

  return(x)
}


set_width.GenomicRanges <- function(x, width = 0L) {
  anchor <- anchors(x)
  if (length(anchor) == 0) {
    width(x) <- width
    return(x)
  }

  for (i in seq_along(anchor)) {
    x <- switch(anchor[[i]],
                start = resize(x, width, fix = "start", ignore.strand = TRUE),
                end = resize(x, width, fix = "end", ignore.strand = TRUE),
                center = resize(x, width, fix = "center", ignore.strand = TRUE),
                "3p" = resize_by_strand(x, width, "+"),
                "5p" = resize_by_strand(x, width, "-"))
  }

  metadata(x)[["anchor"]] <- NULL

  return(x)

}

# anchor by strand version of resize
# filter by strand then resize
resize_by_strand <- function(x, width, strand) {
  x_sub <- x[strand(x) == strand]
  x_sub <- resize(x_sub, width = width)
  x[strand(x) == strand] <- x_sub
  return(x)
}
