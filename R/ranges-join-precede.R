# ranges-join-precede.R

#' Find preceding Ranges
#'
#' @param x,y Ranges objects, which ranges in x precede those in y.
#'
#' @details By default \code{join_precede} will find abritrary ranges
#' in y that are preceded by ranges in x and ignore any strand information.
#'
#' @return A Ranges object with a metadata column called precede that
#' contains the corresponding Ranges that is precedes the ranges x.
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

# does precede left/right make sense? precede will always be right...
# likewise precede will produce ranges downstream

#' Find all preceding Ranges
#'
#' @param x,y Ranges objects, which ranges y precede those on x.
#'
#' @details By definition,  \code{join_precede_right} will find all ranges in y
#' that are on the right-hand side of the ranges in x ignoring any strand
#' information.
#'
#' @return A Ranges object with a metadata column called precedes that
#' contains the corresponding Ranges that is precedes the ranges x.
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

#' Find all preceding Ranges downstream
#'
#' @param x,y Ranges objects, which ranges y precede those on x.
#'
#' @details By definition,  \code{join_precede_right} will find all ranges in y
#' that are that are downstream of the ranges in x. On the positive strand this
#' will result in ranges in y that are right of those in x and on the negative
#' strand it will result in ranges in y that left of those in x.
#'
#' @return A Ranges object with a metadata column called precedes that
#' contains the corresponding Ranges that is precedes the ranges x.
#' @rdname precede-ranges
#' @importFrom IRanges precede
#' @export
join_precede_downstream <- function(x,y) {UseMethod("join_precede_downstream")}

#' @export
join_precede_downstream.GenomicRanges <- function(x,y) {
  hits <- precede(x,y, select = "all", ignore.strand = FALSE)
  nearest_rng_all(x, y, hits, type = "precedes")
}
