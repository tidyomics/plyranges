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

#' Extend the width of a genomic interval
#'
#' @param x a Ranges object, to fix by either the start, end or center
#' of an interval use \code{anchor_start(x)}, \code{anchor_end(x)},
#' \code{anchor_center(x)}.
#' @param width the amount to stretch a Ranges object by. Either an
#' integer vector of length 1 or the same length as x.
#' @importFrom rlang quo enquo eval_tidy
#' @importFrom IRanges resize
#' @seealso \code{\link[IRanges]{resize}}
stretch <- function(x, width = 0L) {

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
                        anchor_center = quo(resize(!!x, width = !!width, fix = "center")),
                        )
  eval_tidy(resize_expr)

}


