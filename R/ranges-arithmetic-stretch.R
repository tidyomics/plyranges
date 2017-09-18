#' Stretch a genomic interval
#'
#' @param x a Ranges object, to fix by either the start, end or center
#' of an interval use \code{anchor_start(x)}, \code{anchor_end(x)},
#' \code{anchor_center(x)}. To fix by strand use \code{anchor_3p(x)} or
#' \code{anchor_5p(x)}.
#' @param extend the amount to alter the width of a Ranges object by. Either an
#' integer vector of length 1 or an integer vector the same length as x.
#' @description Without anchoring, this function will extend the interval
#' in either direction by the integer vector in extend.
#' @export
#' @examples
#' rng <- Ranges(data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0)))
#' rng2 <- stretch(anchor_center(rng), 10)
#' stretch(anchor_start(rng2), 10)
#' stretch(anchor_end(rng2), 10)
#' grng <- Ranges(data.frame(seqnames = "chr1",
#'                          strand = c("+", "-", "-", "+", "+", "-", "+"),
#'                          start=c(2:-1, 13:15),
#'                          width=c(0:3, 2:0)))
#' strech(anchor_3p(grng), 10)
#' strech(anchor_5p(grng), 10)
stretch <- function(x, extend) UseMethod("stretch")

stretch.Ranges <- function(x, extend = 0L) {

  anchor <- metadata(x)[["anchor"]]
  if (length(anchor) == 0) {
    start(x) <- start(x) - extend
    end(x) <- end(x) + extend
    return(x)
  }

  if (any(anchor %in% c("3p", "5p"))) {
    stop("Unable to anchor by strand for IRanges", call. = FALSE)
  }


  for (i in seq_along(anchor)) {
    x <- switch(anchor[[i]],
                start = mutate(x, end = end + extend),
                end = mutate(x, start = start + extend),
                center = stretch_center(x, extend))
  }
  metadata(x)[["anchor"]] <- NULL
  return(x)
}


stretch.GenomicRanges <- function(x, extend = 0L) {
  anchor <- metadata(x)[["anchor"]]
  if (length(anchor) == 0) {
    start(x) <- start(x) - extend
    end(x) <- end(x) + extend
    return(x)
  }

  for (i in seq_along(anchor)) {
    x <- switch(anchor[[i]],
                start = mutate(x, end = end + extend),
                end = mutate(x, start = start + extend),
                center = stretch_center(x, extend),
                "3p" = stretch_by_strand(x, extend, "+"),
                "5p" = stretch_by_strand(x, extend, "-"))
  }
  metadata(x)[["anchor"]] <- NULL
  return(x)
}


stretch_center <- function(x, extend) {
  m <- (end(x) + start(x)) / 2
  ns <-  floor(m - extend)
  ne <- ceiling(m + extend)
  start(x) <- ns
  end(x) <- ne
  x
}

stretch_by_strand <- function(x, extend, strand) {
  x_sub <- x[strand(x) == strand]
  start(x_sub) <- start(x_sub) - extend
  end(x_sub) <- end(x_sub) + extend
  x[strand(x) == strand] <- x_sub
  x
}


