# ranges-join-precede.R

#' Find preceding Ranges
#'
#' @param x,y Ranges objects, which ranges in x precede those in y.
#'
#' @details By default \code{join_precede} will return the ranges
#' in x that come before the ranges in y and ignore any strand information.
#' The function \code{join_precede_right} will find all ranges in y
#' that are on the right-hand side of the ranges in x ignoring any strand
#' information. Finally, \code{join_precede_downstream} will find all ranges in y
#' that are that are downstream of the ranges in x. On the positive strand this
#' will result in ranges in y that are right of those in x and on the negative
#' strand it will result in ranges in y that are left of those in x.
#'
#' @return A Ranges object with a metadata column called precede that
#' contains the corresponding Ranges from y that are preceded by the ranges in x.
#'
#' @rdname precede-ranges
#' @importFrom IRanges precede
#' @export
join_precede <- function(x,y) { UseMethod("join_precede") }

#' @export
join_precede.Ranges <- function(x,y) {
  hits <- precede(x,y)
  nearest_rng(x,y,hits, type = "precedes")
}

#' @export
join_precede.GenomicRanges <- function(x,y) {
  hits <- precede(x,y, ignore.strand = TRUE)
  nearest_rng(x,y, hits, type = "precedes")
}


#' @rdname precede-ranges
#' @importFrom IRanges precede
#' @export
join_precede_right <- function(x,y) { UseMethod("join_precede_right") }

#' @export
join_precede_right.Ranges <- function(x,y) {
  hits <- precede(x,y, select = "all")
  nearest_rng_all(x,y, hits, type = "precedes")
}

#' @export
join_precede_right.GenomicRanges <- function(x,y) {
  hits <- precede(x,y, select = "all", ignore.strand = TRUE)
  nearest_rng_all(x,y, hits, type = "precedes")
}


#' @rdname precede-ranges
#' @importFrom IRanges precede
#' @export
join_precede_downstream <- function(x,y) {UseMethod("join_precede_downstream")}

#' @export
join_precede_downstream.GenomicRanges <- function(x,y) {
  hits <- precede(x,y, select = "all", ignore.strand = FALSE)
  nearest_rng_all(x, y, hits, type = "precedes")
}
