# arithemtic-ranges.R
is_ranges <- function(x) {
  is(class(x), "IRanges") || is(class(x), "GRanges")
}

anchor_start <- function(x) {
  stopifnot(!is_ranges(x))
  metadata(x) <- list(anchor = "start")
  x
}

anchor_end <- function(x) {
  stopifnot(!is_ranges(x))
  metadata(x) <- list(anchor = "end")
  x
}

anchor_center <- function(x) {
  stopifnot(!is_ranges(x))
  metadata(x) <- list(anchor = "center")
  x
}

anchor_3p <- function(x) {
  stopifnot(!is_ranges(x))
  metadata(x) <- list(anchor = "3p")
  x
}

anchor_5p <- function(x) {
  stopifnot(!is_ranges(x))
  metadata(x) <- list(anchor = "5p")
  x
}

#' Set the width of a genomic interval
#'
#' @param x a Ranges object, to fix by either the start, end or center
#' of an interval use \code{anchor_start(x)}, \code{anchor_end(x)},
#' \code{anchor_center(x)}.
#' @param width the amount to alter the width of a Ranges object by. Either an
#' integer vector of length 1 or an integer vector the same length as x.
#' @importFrom rlang quo enquo eval_tidy
#' @importFrom IRanges resize
#' @seealso \code{\link[IRanges]{resize}}
#' @examples
#' rng <- Ranges(data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0)))
#' set_width(rng, width = 10L)
#' set_width(anchor_start(rng), 10L)
#' set_width(anchor_end(rng), 10L)
#' set_width(anchor_center(rng), 10L)
set_width <- function(x, width = 0L) {

  anchor <- metadata(x)[["anchor"]]
  if (is.null(anchor)) {
    width(x) <- width
    return(x)
  }

  if (inherits(x,"IRanges")) {
    return(
      switch(anchor,
                start = resize(x, width = width, fix = "start"),
                end = resize(x, width = width, fix = "end"),
                center = resize(x, width = width, fix = "center"),
                "3p" = stop("Unable to anchor by strand for IRanges", call. = FALSE),
                "5p" = stop("Unable to anchor by strand for IRanges", call. = FALSE)
             )
      )
  }

  if (inherits(x, "GRanges")) {
    return(
      switch(anchor,
           start = resize(x, width = width, fix = "start", ignore.strand = TRUE),
           end = resize(x, width = width, fix = "end", ignore.strand = TRUE),
           center = resize(x, width = width, fix = "center", ignore.strand = TRUE),
           "3p" = {
             x_5p <- x[strand(x) == "+"]
             x_5p <- resize(x_5p, width = width)
             x[strand(x) == "+"] <- x_5p
             x
           },
           "5p" = {
             x_3p <- x[strand(x) == "-"]
             x_3p <- resize(x_3p, width = width)
             x[strand(x) == "-"] <- x_3p
             x

           }
           )
    )
  }

}

#' Stretch a genomic interval
#'
#' @param x a Ranges object, to fix by either the start, end or center
#' of an interval use \code{anchor_start(x)}, \code{anchor_end(x)},
#' \code{anchor_center(x)}.
#' @param extend the amount to alter the width of a Ranges object by. Either an
#' integer vector of length 1 or an integer vector the same length as x.
#' @description Without anchoring, this function will extend the interval
#' in either direction by the integer vector in extend.
#' @importFrom rlang quo enquo eval_tidy
#' @examples
#' rng <- Ranges(data.frame(start=c(2:-1, 13:15), width=c(0:3, 2:0)))
#' rng2 <- stretch(anchor_center(rng), 10)
#' stretch(anchor_start(rng2), 10)
#' stretch(anchor_end(rng2), 10)
stretch <- function(x, extend = 0L) {

  anchor <- metadata(x)[["anchor"]]
  if (is.null(anchor)) {
    start(x) <- start(x) - extend
    end(x) <- end(x) + extend
    return(x)
  }


  if (inherits(x,"IRanges")) {
    return(
      switch(anchor,
                start = mutate(x, end = end + extend),
                end = mutate(x, start = start + extend),
                center = stretch_center(x, extend),
                "3p" = stop("Unable to anchor by strand for IRanges", call. = FALSE),
                "5p" = stop("Unable to anchor by strand for IRanges", call. = FALSE)
             )
      )
  }

  if (inherits(x,"GRanges")) {
    return(
      switch(anchor,
                start = mutate(x, end = end + extend),
                end = mutate(x, start = start + extend),
                center = stretch_center(x, extend),
                "3p" = {
                  x_5p <- x[strand(x) == "+"]
                  start(x_5p) <- start(x_5p) - extend
                  end(x_5p) <- end(x_5p) + extend
                  x[strand(x) == "+"] <- x_5p
                  x
                },
                "5p" = {
                  x_3p <- x[strand(x) == "-"]
                  start(x_3p) <- start(x_3p) - extend
                  end(x_3p) <- end(x_3p) + extend
                  x[strand(x) == "-"] <- x_3p
                  x

                }
             )
    )
  }


}

stretch_center <- function(x, extend) {
  m <- (end(x) + start(x)) / 2
  ns <-  floor(m - extend)
  ne <- ceiling(m + extend)
  start(x) <- ns
  end(x) <- ne
  return(x)
}


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
  x_3p <- x[strand(x) == "-"]
  x_3p <- shift(x_3p, shift)
  x[strand(x) == "-"] <- x_3p
  x
}

#' @rdname shift-ranges
#' @export
shift_downstream <- function(x, shift = 0L) {
  x_5p <- x[strand(x) == "+"]
  x_5p <- shift(x_5p, shift)
  x[strand(x) == "+"] <- x_5p
  x
}

