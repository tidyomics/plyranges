# Element wise set operations

#' Row-wise set operations on Ranges objects
#' @param x,y Ranges objects
#' @details Each of these functions acts on the rows of Ranges object.
#' \code{%union%} will return the entire range between two ranges assuming there
#' are no gaps, if there are gaps use
#' \code{span} instead. \code{%intersect%} will update the Ranges object
#' with a hit column indicating whether or not the two ranges intersect.
#' \code{%setdiff} will return the ranges for each row in x that are not in
#' the corresponding row of y. \code{between} will return the gaps between
#' two ranges.
#' @return A Ranges object
#' @importFrom IRanges punion pintersect pgap psetdiff
#' @export
#' @rdname element-setops
`%union%` <- function(x, y) {
  punion(x,y)
}

#' @export
#' @rdname element-setops
`%intersect%` <- function(x,y) {
  pintersect(x,y)
}
#' @export
#' @rdname element-setops
`%setdiff%` <- function(x, y) {
  psetdiff(x,y)
}

#' @export
#' @rdname element-setops
between <- function(x, y) {
  pgap(x,y)
}

#' @export
#' @rdname element-setops
span <- function(x, y) {
  punion(x, y, fill.gap = TRUE)
}
