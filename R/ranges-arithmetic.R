# arithemtic-ranges.R
anchor_start <- function(x) {
  cl <- match.call()
  cl$x <- enquo(x)
  list(cl = cl, x = x)

}

anchor_end <- function(x) {
  cl <- match.call()
  cl$x <- enquo(x)
  list(cl = cl, x = x)
}

anchor_center <- function(x) {
  cl <- match.call()
  cl$x <- enquo(x)
  list(cl = cl, x = x)
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

  expr <- enquo(x)
  qexpr <- eval_tidy(UQ(expr))
  if (!inherits(qexpr, "Ranges")) {
    fun_name <- gsub("\\(.*", "",qexpr$cl[[1]])
    x <- qexpr$x
  } else {
    x <- qexpr
    width(x) <- width
    return(x)
  }

  resize_expr <- switch(fun_name,
                        anchor_start = quo(resize(!!x, width = !!width, fix = "start")),
                        anchor_end = quo(resize(!!x, width = !!width, fix = "end")),
                        anchor_center = quo(resize(!!x, width = !!width, fix = "center"))
                        )
  eval_tidy(resize_expr)

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
  expr <- enquo(x)
  qexpr <- eval_tidy(UQ(expr))

  if (!inherits(qexpr, "Ranges")) {
    fun_name <- gsub("\\(.*", "",qexpr$cl[[1]])
    x <- qexpr$x
  } else {
    x <- qexpr
    start(x) <- start(x) - extend
    end(x) <- end(x) + extend
    return(x)
  }

  resize_expr <- switch(fun_name,
                        anchor_start = quo(mutate(!!x, end = end + !!extend)),
                        anchor_end = quo(mutate(!!x, start = start + !!extend)),
                        anchor_center = quo(stretch_center(!!x, !!extend))
  )

  eval_tidy(resize_expr)

}

stretch_center <- function(x, extend) {
  m <- (end(x) + start(x)) / 2
  ns <-  floor(m - extend)
  ne <- ceiling(m + extend)
  start(x) <- ns
  end(x) <- ne
  x
}


#' Shift all coordinates in a genomic interval either left or right
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

