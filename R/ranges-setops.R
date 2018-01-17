# Element wise set operations

#' Row-wise set operations on Ranges objects
#'
#' @param x,y Ranges objects
#'
#'
#' @details Each of these functions acts on the rows between pairs of
#' Ranges object.
#' The function \code{\link{\%union\%}}
#' will return the entire range between two ranges objects assuming there
#' are no gaps, if you would like to force gaps use \code{\link{span}} instead.
#' The function \code{\link{\%intersect\%}} will create a new ranges object
#' with a hit column indicating whether or not the two ranges intersect.
#' The function\code{\link{\%setdiff\%}} will return the ranges for each
#' row in x that are not in the corresponding row of y.
#' The function \code{\link{between}} will return the gaps between
#' two ranges.
#'
#' @return A Ranges object
#' @importFrom IRanges punion pintersect pgap psetdiff
#' @seealso \link[IRanges]{punion} \link[IRanges]{pintersect}
#' \link[IRanges]{pgap} \link[IRanges]{psetdiff}
#' @examples
#' x <- as_iranges(data.frame(start = 1:10, width = 5))
#' # stretch x by 3 on the right
#' y <- stretch(anchor_start(x), 3)
#' # take the rowwise union
#' x %union% y
#' # take the rowwise intersection
#' x %intersect% y
#' # asymetric difference
#' y %setdiff% x
#' x %setdiff% y
#' # if there are gaps between the rows of each range use span
#' y <- as_iranges(data.frame(start = c(20:15, 2:5),
#' width = c(10:15,1:4)))
#' # fill in the gaps and take the rowwise union
#' span(x,y)
#' # find the gaps
#' between(x,y)
#' @export
#' @rdname element-setops
`%union%` <- function(x, y) {
  punion(x,y, fill.gap = FALSE)
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
